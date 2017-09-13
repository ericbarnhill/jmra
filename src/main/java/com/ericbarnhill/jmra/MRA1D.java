package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;

public class MRA1D extends MRA<double[], boolean[], double[]> {

    FilterPair af;
    FilterPair sf;
    private int w;
    private int wPad;
    double[] paddedData;
    boolean[] paddedMask;

    public MRA1D() { super(); }

    public MRA1D(ConvolverFactory.ConvolutionType convType) { super(convType); }

    public MRA1D(double[] origData, boolean[] maskData, FilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, fb, decompLvls, convType);
        this.af = fb.af;
        this.sf = fb.sf;
        this.w = origData.length;
        this.wPad = (int)nextPwr2(w);
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, wPad);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, wPad);
        this.stride = 2;
    }

    public MRA1D(double[] origData, boolean[] maskData, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, decompLvls, convType);
        this.af = fb.af;
        this.sf = fb.sf;
        this.w = origData.length;
        this.wPad = (int)nextPwr2(w);
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, wPad);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, wPad);
        this.stride = 2;
    }

    public MRA1D(double[] origData, FilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length), fb, decompLvls, convType);
    }

    @Override
    public ArrayList<double[]> getDecomposition(int localIndex, int ind, int decompLvl, int dimLvl, int localStride) {
        double[] x = new double[0];
        if (dimLvl == 0) {
            if (decompLvl == 0) {
                x = ArrayMath.deepCopy(paddedData);
            } else {
                x = getData(localIndex - stride);
            }
        } else {
            x = getData(ind);
        }
        double[] lo = analysis(x, fb.af.lo, decompLvl);
        double[] hi = analysis(x, fb.af.hi, decompLvl);
        ArrayList<double[]> loAndHi = new ArrayList<double[]>();
        loAndHi.add(lo);
        loAndHi.add(hi);
        return loAndHi;
    }

    @Override
    public double[] getRecomposition(int localPair, int ind, int decompLevel, int decompLvl, int dimLvl) {
        double[] lo = getData(ind);
        double[] hi = getData(ind + localPair);
        double[] y = synthesis(lo, hi, fb.sf.lo, fb.sf.hi, decompLvl);
        return y;
    }

    @Override
     public void accept(Threshold threshold) {
         threshold.visit(this);
    }

    @Override
     public double[] analysis(double[] y, double[] filter, int decompLvl) {
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
     public double[] synthesis(double[] lo, double[] hi, double[] sfl, double[] sfh, int decompLvl) {
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

    @Override
    public double[] getData(int index) {
        double[] data = waveletData.get(index);
        return data;
    }

    @Override
    public void setData(int index, double[] data) {
        waveletData.set(index, data);
    }

    @Override
    public double[] getFilteredData() {
        return ArrayMath.stripBorderPadding(waveletData.get(0), (wPad - w)/2);
    }


}
