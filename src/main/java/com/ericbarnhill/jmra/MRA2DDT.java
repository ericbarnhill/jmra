package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

class MRA2DDT extends MRA2D {

    int w;
    int h;
    int area;
    int wPad;
    int hPad;
    int areaPad;
    int stride;
    MRA1D mra1d;
    double[][] paddedData;
    boolean[][] paddedMask;
    boolean undecimated;
    double[] fafl;
    double[] fafh;
    double[] fsfl;
    double[] fsfh;

    public MRA2DDT(double[][] originalData, boolean[][] maskData, ArrayList<ArrayList<ArrayList<double[]>>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convolutionType) {
        super(originalData, maskData, filterBank.get(0), decompLvls, convolutionType);
        ArrayList<ArrayList<double[]>> firstFilterBank = filterBank.get(1);
        ArrayList<double[]> faf = firstFilterBank.get(0);
        ArrayList<double[]> fsf = firstFilterBank.get(1);
        double[] fafl = faf.get(0);
        double[] fafh = faf.get(1);
        double[] fsfl = fsf.get(0);
        double[] fsfh = fsf.get(1);
    }

    public MRA2DDT(double[][] originalData, ArrayList<ArrayList<ArrayList<double[]>>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convolutionType) {
        this(originalData, ArrayMath.fillWithTrue(originalData.length,originalData[0].length), filterBank, decompLvls, convolutionType);
    }

    @Override
    void decompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl); // 4 for dimLvl 0, 2 for dimLvl 1
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl; // starting point
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[][] x = new double[0][];
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
                case 2:
                    x = ArrayMath.shiftDim(x);
                    break;
            }
            double[][] lo = new double[0][];
            double[][] hi = new double[0][];
            if (decompLvl == 0) {
                lo = AFB(x, fafl, decompLvl);
                hi = AFB(x, fafh, decompLvl);
            } else {
                lo = AFB(x, afl, decompLvl);
                hi = AFB(x, afh, decompLvl);
            }    
            switch(localStride) {
                case 2:
                    lo = ArrayMath.shiftDim(lo);
                    hi = ArrayMath.shiftDim(hi);
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
            double[][] lo = waveletData.get(ind);
            double[][] hi = waveletData.get(ind + localPair);
            switch (localStride) {
                case 2:
                    lo = ArrayMath.shiftDim(lo);
                    hi = ArrayMath.shiftDim(hi);
                    break;
            }
            double[][] y = new double[0][]; 
            if (decompLvl == 0) {
                y = SFB(lo, hi, fsfl, fsfh, decompLvl);
            } else {
                y = SFB(lo, hi, sfl, sfh, decompLvl);
            }
            switch (localStride) {
                case 2:
                    y = ArrayMath.shiftDim(y);
                    break;
            }
            waveletData.set(ind, y);
        }
        if (dimLvl > 0) {
            recompose(decompLvl, dimLvl-1);
        }
        if (decompLvl > 0) {
            waveletData.set(stride*(decompLvl-1), waveletData.get(stride*decompLvl));
        }
        if (decompLvl == 0 && dimLvl == 0) {
            filteredData = ArrayMath.stripBorderPadding(waveletData.get(0), (wPad-w)/2, (hPad-h)/2);
        }
    }

    @Override
  double[][] AFB(double[][] data, double[] filter, int decompLvl) { 
      final int fi = data.length;
      final int fj = data[0].length;
      double[][] filtData = new double[fi][];
      for (int i = 0; i < fi; i++) {
              filtData[i] = mra1d.AFB(data[i], filter, decompLvl);
      }
      return filtData;
    }


    @Override
    double[][] SFB(double[][] lo, double[][] hi, double[] sfl, double[] sfh, int decompLvl) {
        final int fi = lo.length;
        double[][] y = new double[fi][];
        for (int i = 0; i < fi; i++) {
            y[i] = mra1d.SFB(lo[i], hi[i], sfl, sfh, decompLvl);
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
                boolean[][] maskDownsampled = ArrayMath.decimate(paddedMask, decimFac);
                double[] waveletVec = ArrayMath.vectorize(waveletData.get(i));
                boolean[] maskVec = ArrayMath.vectorize(maskDownsampled);
                waveletData.set(i, 
                    ArrayMath.devectorize(
                        Threshold.threshold(
                            waveletVec, maskVec, threshMeth, noiseEstMeth)
                        ,w/decimFac)
                    );
            }
        }
    }
    // for debugging and testing
    public void data2File(double[][] data, String path) {
        int w = data.length;
        int h = data[0].length;
        FloatProcessor fp = new FloatProcessor(w,h);
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                fp.putPixelValue(i,j,data[i][j]);
            }
        }
        ImagePlus ip = new ImagePlus("", fp);
        FileSaver fs = new FileSaver(ip);
        fs.saveAsTiff(path);
    }

}
