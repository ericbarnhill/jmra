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

public class MRA2D extends MRA<double[][], boolean[][], double[]> {

    public int w;
    public int h;
    public int area;
    public int wPad;
    public int hPad;
    public int areaPad;
    public MRA1D mra1d;
    public double[][] paddedData;
    public boolean[][] paddedMask;

    public MRA2D() {
        super();
    }

    public MRA2D(double[][] origData, boolean[][] maskData, FilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, fb, decompLvls, convType);
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

    public MRA2D(double[][] origData, boolean[][] maskData, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, decompLvls, convType);
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

    public MRA2D(double[][] origData, FilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length), fb, decompLvls, convType);
    }

    void initializeWaveletData() {
        for (int n = 0; n < stride * decompLvls; n++) {
            waveletData.add(new double[0][]);
        }
    }

    @Override
    public ArrayList<double[][]> getDecomposition(int localIndex, int ind, int decompLvl, int dimLvl, int localStride) {
        double[][] x = new double[0][];
        if (dimLvl == 0) {
            if (decompLvl == 0) {
                x = ArrayMath.deepCopy(paddedData);
            } else {
                x = waveletData.get(localIndex - stride);
            }
        } else {
            x = waveletData.get(ind);
        }
        switch(localStride) {
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
        ArrayList<double[][]> loAndHi = new ArrayList<double[][]>();
        loAndHi.add(lo);
        loAndHi.add(hi);
        return loAndHi;
    }

    @Override
    public double[][] getRecomposition(int localPair, int ind, int decompLvl, int dimLvl, int localStride) {
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
        return y;
    }

    @Override
    public double[][] AFB(double[][] data, double[] filter, int decompLvl) { 
      final int fi = data.length;
      final int fj = data[0].length;
      double[][] filtData = new double[fi][];
      for (int i = 0; i < fi; i++) {
              filtData[i] = mra1d.AFB(data[i], filter, decompLvl);
      }
      return filtData;
    }


    @Override
    public double[][] SFB(double[][] lo, double[][] hi, double[] sfl, double[] sfh, int decompLvl) {
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
