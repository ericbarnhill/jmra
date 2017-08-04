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
    FilterPair faf;
    FilterPair fsf;

    public MRA2DDT(double[][] origData, boolean[][] maskData, DTFilterSet dtfs, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, new FilterBank(dtfs.f1, dtfs.f2), decompLvls, convType);
        faf = dtfs.ff1;
        fsf = dtfs.ff2;
    }

    public MRA2DDT(double[][] origData, DTFilterSet dtfs, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length), dtfs, decompLvls, convType);
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
            double[][] lo = new double[0][];
            double[][] hi = new double[0][];
            switch(localStride) { // shift dim for y processing. done as a switch block so higher dim code can all use the identical form
                case 2:
                    x = ArrayMath.shiftDim(x);
                    if (decompLvl == 0) {
                        lo = AFB(x, faf2.lo, decompLvl);
                        hi = AFB(x, faf2.hi, decompLvl);
                    } else {
                        lo = AFB(x, af2.lo, decompLvl);
                        hi = AFB(x, af2.hi, decompLvl);
                    }    
                    lo = ArrayMath.shiftDim(lo);
                    hi = ArrayMath.shiftDim(hi);
                    break;
                case 4:
                    if (decompLvl == 0) {
                        lo = AFB(x, faf1.lo, decompLvl);
                        hi = AFB(x, faf1.hi, decompLvl);
                    } else {
                        lo = AFB(x, af1.lo, decompLvl);
                        hi = AFB(x, af1.hi, decompLvl);
                    }    
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
            double[][] y = new double[0][]; 
            switch (localStride) {
                case 2:
                    lo = ArrayMath.shiftDim(lo);
                    hi = ArrayMath.shiftDim(hi);
                    if (decompLvl == 0) {
                        y = SFB(lo, hi, fsf2.lo, fsf2.hi, decompLvl);
                    } else {
                        y = SFB(lo, hi, sf2.lo, sf2.hi, decompLvl);
                    }
                    y = ArrayMath.shiftDim(y);
                    break;
                case 4:
                    if (decompLvl == 0) {
                        y = SFB(lo, hi, fsf1.lo, fsf1.hi, decompLvl);
                    } else {
                        y = SFB(lo, hi, sf1.lo, sf1.hi, decompLvl);
                    }
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
