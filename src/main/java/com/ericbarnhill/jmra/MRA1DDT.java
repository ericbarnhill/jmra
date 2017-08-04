package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;

public class MRA1DDT extends MRA1D {

    private int w;
    private int wPad;
    double[] paddedData;
    boolean[] paddedMask;
    DTFilterBank dtfb;

    public MRA1DDT() { super(); }

    public MRA1DDT(ConvolverFactory.ConvolutionType convType) { super(convType); }

    public MRA1DDT(double[] origData, boolean[] maskData, DTFilterBank dtfb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, dtfb.get(0), decompLvls, convType);
        this.w = origData.length;
        this.wPad = (int)nextPwr2(w);
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, wPad);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, wPad);
        this.stride = 2;
    }

    public MRA1DDT(double[] origData, DTFilterBank dtfb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length), dtfb.get(0), decompLvls, convType);
    }

    @Override
    void decompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, 2 - dimLvl); // 4 for dimLvl 0, 2 for dimLvl 1
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl; // starting point
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[] x = new double[0];
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
            double[] lo = new double[0];
            double[] hi = new double[0];
            if (decompLvl == 0) {
                lo = AFB(x, faf.lo, decompLvl);
                hi = AFB(x, faf.hi, decompLvl);
            } else {
                lo = AFB(x, af.lo, decompLvl);
                hi = AFB(x, af.hi, decompLvl);
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
        int localStride = (int)Math.pow(2, 2 - dimLvl);
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl;
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[] lo = waveletData.get(ind);
            double[] hi = waveletData.get(ind + localPair);
            double[] y  = new double[0];
            if (decompLvl == 0) {
                y = SFB(lo, hi, fsf.lo, fsf.hi, decompLvl);
            } else {
                y = SFB(lo, hi, sf.lo, sf.hi, decompLvl);
            }
            waveletData.set(ind, y);
        }
        if (dimLvl > 0) {
            recompose(decompLvl, dimLvl-1);
        } else if (decompLvl > 0) {
            System.out.println("Moving from "+stride*decompLvl+" to "+stride*(decompLvl-1));
            waveletData.set(stride*(decompLvl-1), waveletData.get(stride*decompLvl));
        }
    }

    @Override
     void accept(Threshold threshold) {
         threshold.visit(this);
    }

    @Override
     double[] AFB(double[] y, double[] filter, int decompLvl) {
        final int N = y.length/2;
        final int L = filter.length;
        y = Shifter.circShift(y, -L/2);
        y = upFirDn.upFirDn(y,filter, 1, 2);
        for (int i = 0; i < L/2; i++) {
            y[i] += y[N + i];
        }
        y = ArrayMath.crop(y,N);
        return y;
    }

    @Override
     double[] SFB(double[] lo, double[] hi, double[] sfl, double[] sfh, int decompLvl) {
        final int N = 2*lo.length;
        final int L0 = sfl.length;
        lo = upFirDn.upFirDn(lo, sfl, 2, 1);
        hi = upFirDn.upFirDn(hi, sfh, 2, 1);
        double[] y = ArrayMath.add(lo, hi);
        for (int i = 0; i < L0 - 2; i++) {
            y[i] += y[N + i];
        }
        y = ArrayMath.crop(y,N);
        y = Shifter.circShift(y, 1-L0/2);
        return y;
     }


}
