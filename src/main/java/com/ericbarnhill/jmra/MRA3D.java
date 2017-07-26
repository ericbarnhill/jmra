package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

class MRA3D extends MRA<double[][][], boolean[][][], double[]> {

    private int w;
    private int h;
    private int d;
    private int area;
    private int volume;
    private int wPad;
    private int hPad;
    private int dPad;
    private int areaPad;
    private int volumePad;
    private int stride;
    private double[][][] paddedData;
    private boolean[][][] paddedMask;
    private MRA1D mra1d;

    public MRA3D(double[][][] originalData, boolean[][][] maskData, ArrayList<ArrayList<double[]>> filterBank, int decompositionLevels) {
        super(originalData, maskData, filterBank, decompositionLevels);
        this.w = originalData.length;
        this.h = originalData[0].length;
        this.d = originalData[0][0].length;
        this.area = w*h;
        this.volume = w*h*d;
        this.wPad = (int)nextPwr2(w);
        this.hPad = (int)nextPwr2(h);
        this.dPad = (int)nextPwr2(d);
        this.areaPad = wPad * hPad;
        this.volumePad = wPad * hPad * dPad;
        this.paddedData = JVCLUtils.zeroPadBoundaries(originalData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.paddedMask = JVCLUtils.zeroPadBoundaries(maskData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.stride = 8;
        ArrayList<Integer> LA = new ArrayList<Integer>(2);
        ArrayList<Integer> LS = new ArrayList<Integer>(2);
        this.L = new ArrayList<ArrayList<Integer>>();
        L.add(LA);
        L.add(LS);
        mra1d = new MRA1D();
        setFilterLengths();
    }

    public MRA3D(double[][][] originalData, ArrayList<ArrayList<double[]>> filterBank, int decompositionLevels) {
        this(originalData, ArrayMath.fillWithTrue(originalData.length,originalData[0].length, originalData[0][0].length), filterBank, decompositionLevels);
    }

    private void setFilterLengths() {
        this.filterBank = filterBank;
        for (int i = 0; i < filterBank.size(); i++) {
            for (double[] d : filterBank.get(i)) {
                L.get(i).add(d.length);
            }
        }
    }

    public void dwt() {
        scalingData = ArrayMath.deepCopy(paddedData);
        for (int decompositionLevel = 0; decompositionLevel < decompositionLevels; decompositionLevel++) {
            decompose(scalingData, 2);
            scalingData = waveletData.get(decompositionLevel*stride);
        }
    }

    public void idwt() {
        for (int decompositionLevel = decompositionLevels-1; decompositionLevel >= 0; decompositionLevel--) {
            recompose(decompositionLevel, 2);
        }
        filteredData = JVCLUtils.stripBorderPadding(waveletData.get(0), (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
    }

    void decompose(double[][][] data, int dimensionLevel) {
        double[][][] lo = AFB(data, h0);
        double[][][] hi = AFB(data, h1);
        lo = ArrayMath.shiftDim(lo);
        hi = ArrayMath.shiftDim(hi);
        if (dimensionLevel == 0) {
            waveletData.add(lo);
            waveletData.add(hi);
        } else {
            decompose(lo, dimensionLevel-1);
            decompose(hi, dimensionLevel-1);
        }
    }

    double[][][] AFB(double[][][] data, double[] filter) {
        final int fi = data.length;
        final int fj = data[0].length;
        final int fk = data[0][0].length;
        final int fk2 = fk / 2;
        double[][][] filtData = new double[fi][fj][fk];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                filtData[i][j] = mra1d.AFB(data[i][j], filter);
            }
        }
        return filtData;
    }


    void recompose(int decompositionLevel, int dimensionLevel) {
        int localStride = (int)Math.pow(2, 2 - dimensionLevel);
        int localPair = localStride / 2;
        System.out.println(localStride);
        int localIndex = stride*decompositionLevel;
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[][][] lo = waveletData.get(ind);
            double[][][] hi = waveletData.get(ind + localPair);
            switch (localStride) {
                case 4: 
                    lo = ArrayMath.shiftDim(lo, 2);
                    hi = ArrayMath.shiftDim(hi, 2);
                    break;
                case 2:
                    lo = ArrayMath.shiftDim(lo, 1);
                    hi = ArrayMath.shiftDim(hi, 1);
                    break;
            }
            double[][][] y = SFB(lo, hi, g0, g1);
            switch (localStride) {
                case 4: 
                    y = ArrayMath.shiftDim(y, 1);
                    break;
                case 2:
                    y = ArrayMath.shiftDim(y, 2);
                    break;
            }
            waveletData.set(ind, y);
        }
        if (dimensionLevel > 0) {
            recompose(decompositionLevel, dimensionLevel-1);
        } else if (decompositionLevel > 0) {
            waveletData.set(stride*(decompositionLevel-1), waveletData.get(stride*decompositionLevel));
        }
    }
    
    double[][][] SFB(double[][][] lo, double[][][] hi, double[] g0, double[] g1) {
        final int fi = lo.length;
        final int fj = lo[0].length;
        final int fk = lo[0][0].length*2;
        double[][][] y = new double[fi][fj][fk];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                y[i][j] = mra1d.SFB(lo[i][j], hi[i][j], g0, g1);
            }
        }
        return y;
    }

    public void threshold(Threshold.ThreshMeth threshMeth, Threshold.NoiseEstMeth noiseEstMeth) {
        // loop through each subband, pass method
        for (int i = 0; i < waveletData.size(); i++) {
            // avoid scaling datas
            if (i % stride != 0) {
                int level = (int)Math.floor(i / stride);
                int decimFac = (int)Math.pow(2, level+1);
                boolean[][][] maskDownsampled = ArrayMath.decimate(paddedMask, decimFac);
                double[] waveletVec = ArrayMath.vectorize(waveletData.get(i));
                boolean[] maskVec = ArrayMath.vectorize(maskDownsampled);
                waveletData.set(i, 
                    ArrayMath.devectorize(
                        Threshold.threshold(
                            waveletVec, maskVec, threshMeth, noiseEstMeth)
                        ,wPad/decimFac, hPad/decimFac)
                    );
            }
        }
    }
    // for debugging and testing
    public void data2File(double[][][] data, String path) {
        int w = data.length;
        int h = data[0].length;
        int d = data[0][0].length;
        ImageStack is = new ImageStack(w,h);
        for (int k = 0; k < d; k++) {
            FloatProcessor fp = new FloatProcessor(w,h);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                            fp.putPixelValue(i,j,data[i][j][k]);
                }
            }
            is.addSlice(fp);
        }
        ImagePlus ip = new ImagePlus("", is);
        FileSaver fs = new FileSaver(ip);
        fs.saveAsTiff(path);
    }



        


}
