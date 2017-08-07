package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

public class MRA3D extends MRA<double[][][], boolean[][][], double[]> {

     public int w;
     public int h;
     public int d;
     public int area;
     public int volume;
     public int wPad;
     public int hPad;
     public int dPad;
     public int areaPad;
     public int volumePad;
     public double[][][] paddedData;
     public boolean[][][] paddedMask;
     public MRA1D mra1d;

     public MRA3D() {
         super();
     }

    public MRA3D(double[][][] origData, boolean[][][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        this.w = origData.length;
        this.h = origData[0].length;
        this.d = origData[0][0].length;
        this.area = w*h;
        this.volume = w*h*d;
        this.wPad = (int)nextPwr2(w);
        this.hPad = (int)nextPwr2(h);
        this.dPad = (int)nextPwr2(d);
        this.areaPad = wPad * hPad;
        this.volumePad = wPad * hPad * dPad;
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.dimLvls = 3;
        this.stride = 8;
        mra1d = new MRA1D(convType);
        initializeWaveletData();
    }

    public MRA3D(double[][][] origData, boolean[][][] maskData, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, decompLvls, convType);
        this.w = origData.length;
        this.h = origData[0].length;
        this.d = origData[0][0].length;
        this.area = w*h;
        this.volume = w*h*d;
        this.wPad = (int)nextPwr2(w);
        this.hPad = (int)nextPwr2(h);
        this.dPad = (int)nextPwr2(d);
        this.areaPad = wPad * hPad;
        this.volumePad = wPad * hPad * dPad;
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.dimLvls = 3;
        this.stride = 8;
        mra1d = new MRA1D(convType);
        initializeWaveletData();
    }

    public MRA3D(double[][][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length, origData[0][0].length), filterBank, decompLvls, convType);
    }

    void initializeWaveletData() {
        for (int n = 0; n < stride * decompLvls; n++) {
            waveletData.add(new double[0][][]);
        }
    }

    @Override
    public ArrayList<double[][][]> getDecomposition(int localIndex, int ind, int decompLvl, int dimLvl, int localStride) {
        double[][][] x = new double[0][][];
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
            case 4: 
                x = ArrayMath.shiftDim(x, 1);
                break;
            case 2:
                x = ArrayMath.shiftDim(x, 2);
                break;
        }
        double[][][] lo = AFB(x, fb.af.lo, decompLvl);
        double[][][] hi = AFB(x, fb.af.hi, decompLvl);
        switch(localStride) {
            case 4: 
                lo = ArrayMath.shiftDim(lo, 2);
                hi = ArrayMath.shiftDim(hi, 2);
                break;
            case 2:
                lo = ArrayMath.shiftDim(lo, 1);
                hi = ArrayMath.shiftDim(hi, 1);
                break;
        }
        ArrayList<double[][][]> loAndHi = new ArrayList<double[][][]>();
        loAndHi.add(lo);
        loAndHi.add(hi);
        return loAndHi;
    }

    @Override
    public double[][][] getRecomposition(int localPair, int ind, int decompLvl, int dimLvl, int localStride) {
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
        double[][][] y = SFB(lo, hi, fb.sf.lo, fb.sf.hi, decompLvl);
        switch (localStride) {
            case 4: 
                y = ArrayMath.shiftDim(y, 2);
                break;
            case 2:
                y = ArrayMath.shiftDim(y, 1);
                break;
        }
        return y;
    }
    
    @Override
    public double[][][] AFB(double[][][] data, double[] filter, int decompLvl) {
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
    public double[][][] SFB(double[][][] lo, double[][][] hi, double[] sfl, double[] sfh, int decompLvl) {
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


    public void accept(Threshold threshold) {
        threshold.visit(this);
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
