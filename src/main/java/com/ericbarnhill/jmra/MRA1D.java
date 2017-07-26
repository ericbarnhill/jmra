package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import java.util.ArrayList;
import java.util.Arrays;

class MRA1D extends MRA<double[], boolean[], double[]> {

    private int w;
    private int wPad;
    private int stride;

    public MRA1D() { super(); }

    public MRA1D(double[] originalData, boolean[] maskData, ArrayList<ArrayList<double[]>> filterBank, int decompositionLevels) {
        super(originalData, maskData, filterBank, decompositionLevels);
        this.w = originalData.length;
        this.wPad = (int)nextPwr2(w);
        this.paddedData = JVCLUtils.zeroPadBoundaries(originalData, wPad);
        this.stride = 2;
        ArrayList<Integer> LA = new ArrayList<Integer>(2);
        ArrayList<Integer> LS = new ArrayList<Integer>(2);
        this.L = new ArrayList<ArrayList<Integer>>();
        L.add(LA);
        L.add(LS);
        setFilterLengths();
    }

    public MRA1D(double[] originalData, ArrayList<ArrayList<double[]>> filterBank, int decompositionLevels) {
        this(originalData, ArrayMath.fillWithTrue(originalData.length), filterBank, decompositionLevels);
    }

    private void setFilterLengths() {
        this.filterBank = filterBank;
        for (int i = 0; i < filterBank.size(); i++) {
            for (double[] d : filterBank.get(i)) {
                L.get(i).add(d.length);
            }
        }
    }

    public void dwt() {
        scalingData = ArrayMath.deepCopy(originalData);
        for (int decompositionLevel = 0; decompositionLevel < decompositionLevels; decompositionLevel++) {
            decompose(scalingData, 1);
            scalingData = waveletData.get(decompositionLevel*stride);
        }
    }

    public void idwt() {
        for (int decompositionLevel = decompositionLevels-1; decompositionLevel >= 0; decompositionLevel--) {
            recompose(decompositionLevel, 1);
        }
        filteredData = ArrayMath.deepCopy(waveletData.get(0));
    }

    void decompose(double[] data, int dimensionLevel) {
        final int fi = data.length/2;
        double[] lo = new double[fi];
        double[] hi = new double[fi];
        lo = AFB(data, h0);
        hi = AFB(data, h1);
        if (dimensionLevel == 0) {
            waveletData.add(lo);
            waveletData.add(hi);
        } else {
            throw new RuntimeException("MRA1D: Dimension level error");
        }
    }

    void recompose(int decompositionLevel, int dimensionLevel) {
        int localStride = (int)Math.pow(2, 2 - dimensionLevel);
        int localPair = localStride / 2;
        int localIndex = stride*decompositionLevel;
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[] lo = waveletData.get(ind);
            double[] hi = waveletData.get(ind + localPair);
            final int fi = lo.length;
            double[] y = SFB(lo, hi, g0, g1);
            waveletData.set(ind, y);
        }
        if (dimensionLevel > 0) {
            recompose(decompositionLevel, dimensionLevel-1);
        } else if (decompositionLevel > 0) {
            System.out.println("Moving from "+stride*decompositionLevel+" to "+stride*(decompositionLevel-1));
            waveletData.set(stride*(decompositionLevel-1), waveletData.get(stride*decompositionLevel));
        }
    }

    public void threshold(Threshold.ThreshMeth threshMeth, Threshold.NoiseEstMeth noiseEstMeth) {
        // loop through each subband, pass method
        for (int i = 0; i < waveletData.size(); i++) {
            // avoid scaling datas
            if (i % stride != 0) {
                int level = (int)Math.floor(i / stride);
                boolean[] maskDownsampled = ArrayMath.decimate(maskData, (int)Math.pow(2,level));
                waveletData.set(i, Threshold.threshold(waveletData.get(i), maskDownsampled, threshMeth, noiseEstMeth));
            }
        }
    }

     double[] AFB(double[] y, double[] filter) {
        final int N = y.length/2;
        final int L = filter.length;
        y = ArrayMath.deepCopy(y);
        y = Shifter.circShift(y, -L/2);
        y = UpFirDn.upFirDn(y,filter, 1, 2);
        for (int i = 0; i < L/2; i++) {
            y[i] += y[N + i];
        }
        y = ArrayMath.crop(y,N);
        return y;
    }

     double[] SFB(double[] lo, double[] hi, double[] g0, double[] g1) {
        final int N = 2*lo.length;
        final int L0 = g0.length;
        lo = UpFirDn.upFirDn(lo, g0, 2, 1);
        hi = UpFirDn.upFirDn(hi, g1, 2, 1);
        double[] y = ArrayMath.add(lo, hi);
        for (int i = 0; i < L0 - 2; i++) {
            y[i] += y[N + i];
        }
        y = ArrayMath.crop(y,N);
        y = Shifter.circShift(y, 1-L0/2);
        return y;
    }
    

        


}
