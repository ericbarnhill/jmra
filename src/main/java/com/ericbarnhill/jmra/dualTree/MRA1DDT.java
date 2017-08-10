package com.ericbarnhill.jmra.dualTree;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;

public class MRA1DDT extends MRA1D {

    private int w;
    private int wPad;
    double[] paddedData;
    boolean[] paddedMask;
    DTFilterBank fb;

    public MRA1DDT() { super(); }

    public MRA1DDT(ConvolverFactory.ConvolutionType convType) { super(convType); }

    public MRA1DDT(double[] origData, boolean[] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, decompLvls, convType);
        this.w = origData.length;
        this.wPad = (int)nextPwr2(w);
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, wPad);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, wPad);
        this.stride = 2;
        this.fb = fb;
    }

    public MRA1DDT(double[] origData, int decompLvls, DTFilterBank fb, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length), fb, decompLvls, convType);
    }

    @Override
    public ArrayList<double[]> getDecomposition(int localIndex, int ind, int decompLvl, int dimLvl, int localStride) {
        double[] x = new double[0];
        if (dimLvl == 0) {
            if (decompLvl == 0) {
                x = ArrayMath.deepCopy(paddedData);
            } else {
                x = waveletData.get(localIndex - stride);
            }
        } else {
            x = waveletData.get(ind);
        }
        double[] lo = new double[0];
        double[] hi = new double[0];
        if (decompLvl == 0) {
            lo = AFB(x, fb.faf.get(0).lo, decompLvl);
            hi = AFB(x, fb.faf.get(0).hi, decompLvl);
        } else {
            lo = AFB(x, fb.af.get(0).lo, decompLvl);
            hi = AFB(x, fb.af.get(0).hi, decompLvl);
        }    
        ArrayList<double[]> loAndHi = new ArrayList<double[]>();
        loAndHi.add(lo);
        loAndHi.add(hi);
        return loAndHi;
    }

    @Override
    public double[] getRecomposition(int localPair, int ind, int decompLvl, int dimLvl, int localStride) {
        double[] lo = waveletData.get(ind);
        double[] hi = waveletData.get(ind + localPair);
        double[] y  = new double[0];
        if (decompLvl == 0) {
            y = SFB(lo, hi, fb.fsf.get(0).lo, fb.fsf.get(0).hi, decompLvl);
        } else {
            y = SFB(lo, hi, fb.sf.get(0).lo, fb.sf.get(0).hi, decompLvl);
        }
        return y;
    }

    @Override
     public void accept(Threshold threshold) {
         threshold.visit(this);
    }

    @Override
     public double[] AFB(double[] y, double[] filter, int decompLvl) {
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
     public double[] SFB(double[] lo, double[] hi, double[] sfl, double[] sfh, int decompLvl) {
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
