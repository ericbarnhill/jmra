package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

class MRA2D extends MRA<double[][], boolean[][], double[]> {

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

    public MRA2D() {
        super();
    }

    public MRA2D(double[][] origData, boolean[][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        this.af = filterBank.af;
        this.sf = filterBank.sf;
        this.w = origData.length;
        this.h = origData[0].length;
        this.area = w*h;
        this.wPad = (int)nextPwr2(w);
        this.hPad = (int)nextPwr2(h);
        this.areaPad = wPad * hPad;
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, (wPad-w)/2, (hPad-h)/2);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, (wPad-w)/2, (hPad-h)/2);
        this.stride = 4;
        this.dimLvls = 2;
        this.waveletData = new ArrayList<double[][]>(20);
        mra1d = new MRA1D(convType);
        initializeWaveletData();
    }

    public MRA2D(double[][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length), filterBank, decompLvls, convType);
    }

    void initializeWaveletData() {
        for (int n = 0; n < stride * decompLvls; n++) {
            waveletData.add(new double[0][]);
        }
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
            double[][] lo = AFB(x, fb.af.lo, decompLvl);
            double[][] hi = AFB(x, fb.af.hi, decompLvl);
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
            double[][] y = SFB(lo, hi, fb.sf.lo, fb.sf.hi, decompLvl);
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

    public void accept(Threshold threshold) {
        threshold.visit(this);
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
