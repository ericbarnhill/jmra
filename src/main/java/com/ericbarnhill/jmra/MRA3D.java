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

     int w;
     int h;
     int d;
     int area;
     int volume;
     int wPad;
     int hPad;
     int dPad;
     int areaPad;
     int volumePad;
     int stride;
     double[][][] paddedData;
     boolean[][][] paddedMask;
     MRA1D mra1d;

    public MRA3D(double[][][] originalData, boolean[][][] maskData, ArrayList<ArrayList<double[]>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convolutionType) {
        super(originalData, maskData, filterBank, decompLvls, convolutionType);
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
        this.paddedData = ArrayMath.zeroPadBoundaries(originalData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.dimLvls = 3;
        this.stride = 8;
        mra1d = new MRA1D(convolutionType);
        initializeWaveletData();
    }

    public MRA3D(double[][][] originalData, ArrayList<ArrayList<double[]>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convolutionType) {
        this(originalData, ArrayMath.fillWithTrue(originalData.length,originalData[0].length, originalData[0][0].length), filterBank, decompLvls, convolutionType);
    }

    void initializeWaveletData() {
        for (int n = 0; n < stride * decompLvls; n++) {
            waveletData.add(new double[0][][]);
        }
    }

    @Override
    void decompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl);
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl; // starting point
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[][][] x = new double[0][][];
            // figure out where the scaling image is coming from
            if (dimLvl == 0) {
                if (decompLvl == 0) {
                    x = ArrayMath.deepCopy(paddedData);
                } else {
                    x = waveletData.get(localIndex - stride);
                }
            } else {
                x = waveletData.get(ind);
            }
            // decompose into lo and hi
            switch(localStride) { // shift dim for y processing. done as a switch block so higher dim code can all use the identical form
                case 4: 
                    x = ArrayMath.shiftDim(x, 2);
                    break;
                case 2:
                    x = ArrayMath.shiftDim(x, 1);
                    break;
            }
            double[][][] lo = AFB(x, afl, decompLvl);
            double[][][] hi = AFB(x, afh, decompLvl);
            switch(localStride) {
                case 4: 
                    lo = ArrayMath.shiftDim(lo, 1);
                    hi = ArrayMath.shiftDim(hi, 1);
                    break;
                case 2:
                    lo = ArrayMath.shiftDim(lo, 2);
                    hi = ArrayMath.shiftDim(hi, 2);
                    break;
            }
            waveletData.set(ind, lo);
            waveletData.set(ind + localPair, hi);
        }
        if (dimLvl < dimLvls - 1) {
            decompose(decompLvl, dimLvl+1);
       }
    }

    @Override
    void recompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl);
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl;
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[][][] lo = waveletData.get(ind);
            double[][][] hi = waveletData.get(ind + localPair);
            switch (localStride) {
                case 4: 
                    lo = ArrayMath.shiftDim(lo, 1);
                    hi = ArrayMath.shiftDim(hi, 1);
                    break;
                case 2:
                    lo = ArrayMath.shiftDim(lo, 2);
                    hi = ArrayMath.shiftDim(hi, 2);
                    break;
            }
            double[][][] y = SFB(lo, hi, sfl, sfh, decompLvl);
            switch (localStride) {
                case 4: 
                    y = ArrayMath.shiftDim(y, 2);
                    break;
                case 2:
                    y = ArrayMath.shiftDim(y, 1);
                    break;
            }
            waveletData.set(ind, y);
        }
        if (dimLvl > 0) {
            recompose(decompLvl, dimLvl-1);
        } else if (decompLvl > 0) {
            waveletData.set(stride*(decompLvl-1), waveletData.get(stride*decompLvl));
        }
    }
    
    @Override
    double[][][] AFB(double[][][] data, double[] filter, int decompLvl) {
        final int fi = data.length;
        final int fj = data[0].length;
        final int fk = data[0][0].length;
        final int fk2 = fk / 2;
        double[][][] filtData = new double[fi][fj][fk];
        for (int i = 0; i < fi; i++) { 
            for (int j = 0; j < fj; j++) {
                filtData[i][j] = mra1d.AFB(data[i][j], filter, decompLvl); 
            } 
        }
        return filtData;
    } 

    @Override
    double[][][] SFB(double[][][] lo, double[][][] hi, double[] sfl, double[] sfh, int decompLvl) {
        final int fi = lo.length;
        final int fj = lo[0].length;
        final int fk = lo[0][0].length*2;
        double[][][] y = new double[fi][fj][fk];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                y[i][j] = mra1d.SFB(lo[i][j], hi[i][j], sfl, sfh, decompLvl);
            }
        }
        return y;
    }

    @Override
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
