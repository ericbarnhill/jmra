package com.ericbarnhill.jmra.dualTree;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

class MRA3DDT extends MRA3D {

    int w;
    int h;
    int area;
    int wPad;
    int hPad;
    int areaPad;
    int stride;
    MRA1D mra1d;
    double[][][] paddedData;
    boolean[][][] paddedMask;
    boolean undecimated;
    DTFilterBank fb;

    public MRA3DDT(double[][][] origData, boolean[][][] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, decompLvls, convType);
        this.fb = fb;
    }

    public MRA3DDT(double[][][] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length, origData[0][0].length), fb, decompLvls, convType);
    }

    @Override
    public void decompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl); // 4 for dimLvl 0, 2 for dimLvl 1
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
            double[][][] lo = new double[0][][];
            double[][][] hi = new double[0][][];
            switch(localStride) { // shift dim for y processing. done as a switch block so higher dim code can all use the identical form
                case 8: 
                    if (decompLvl == 0) {
                        lo = AFB(x, fb.faf.get(0).lo, decompLvl);
                        hi = AFB(x, fb.faf.get(0).hi, decompLvl);
                    } else {
                        lo = AFB(x, fb.af.get(0).lo, decompLvl);
                        hi = AFB(x, fb.af.get(0).hi, decompLvl);
                    }    
                    break;
                case 4: 
                    x = ArrayMath.shiftDim(x, 2);
                    if (decompLvl == 0) {
                        lo = AFB(x, fb.faf.get(1).lo, decompLvl);
                        hi = AFB(x, fb.faf.get(1).hi, decompLvl);
                    } else {
                        lo = AFB(x, fb.af.get(1).lo, decompLvl);
                        hi = AFB(x, fb.af.get(1).hi, decompLvl);
                    }    
                    lo = ArrayMath.shiftDim(lo, 1);
                    hi = ArrayMath.shiftDim(hi, 1);
                    break;
                case 2:
                    x = ArrayMath.shiftDim(x, 1);
                    if (decompLvl == 0) {
                        lo = AFB(x, fb.faf.get(2).lo, decompLvl);
                        hi = AFB(x, fb.faf.get(2).hi, decompLvl);
                    } else {
                        lo = AFB(x, fb.af.get(2).lo, decompLvl);
                        hi = AFB(x, fb.af.get(2).hi, decompLvl);
                    }    
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
    public void recompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl);
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl;
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[][][] lo = waveletData.get(ind);
            double[][][] hi = waveletData.get(ind + localPair);
            double[][][] y = new double[0][][]; 
            switch (localStride) {
                case 8:
                    if (decompLvl == 0) {
                        y = SFB(lo, hi, fb.fsf.get(0).lo, fb.fsf.get(0).hi, decompLvl);
                    } else {
                        y = SFB(lo, hi, fb.sf.get(0).lo, fb.sf.get(0).hi, decompLvl);
                    }
                    break;
                case 4: 
                    lo = ArrayMath.shiftDim(lo, 1);
                    hi = ArrayMath.shiftDim(hi, 1);
                    if (decompLvl == 0) {
                        y = SFB(lo, hi, fb.fsf.get(1).lo, fb.fsf.get(1).hi, decompLvl);
                    } else {
                        y = SFB(lo, hi, fb.sf.get(1).lo, fb.sf.get(1).hi, decompLvl);
                    }
                    y = ArrayMath.shiftDim(y, 2);
                    break;
                case 2:
                    lo = ArrayMath.shiftDim(lo, 2);
                    hi = ArrayMath.shiftDim(hi, 2);
                    if (decompLvl == 0) {
                        y = SFB(lo, hi, fb.fsf.get(2).lo, fb.fsf.get(2).hi, decompLvl);
                    } else {
                        y = SFB(lo, hi, fb.sf.get(2).lo, fb.sf.get(2).hi, decompLvl);
                    }
                    y = ArrayMath.shiftDim(y, 1);
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

}
