package com.ericbarnhill.jmra; 
import java.util.ArrayList;
import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.UpFirDn;

public abstract class MRA<N, B, V> {
    // N is ND array of numeric type
    // B is ND boolean array
    // V is 1D vector array of numeric type
     N originalData;
     N scalingData;
     B maskData;
     ArrayList<N> waveletData;
     N filteredData;
     int decompositionLevels; // also used for stride
     int dimensionLevels;
     ArrayList<ArrayList<V>> filterBank; // order h0, h1, g0, g1
     ArrayList<V> analysisFilters;
     V h0;
     V h1;
     V g0;
     V g1;
     ArrayList<V> synthesisFilters;
     ArrayList<ArrayList<Integer>> L;

     public MRA() {}
    public MRA(N originalData, B maskData, ArrayList<ArrayList<V>> filterBank, int decompositionLevels) {
        this.originalData = originalData;
        this.scalingData = originalData; 
        this.maskData = maskData;
        this.filterBank = filterBank;
        this.analysisFilters = filterBank.get(0);
        this.h0 = analysisFilters.get(0);
        this.h1 = analysisFilters.get(1);
        this.synthesisFilters = filterBank.get(1);
        this.g0 = synthesisFilters.get(0);
        this.g1 = synthesisFilters.get(1);
        this.decompositionLevels = decompositionLevels;
        waveletData = new ArrayList<N>();
    } 

    public abstract void dwt();
    public abstract void idwt();
    public abstract void threshold(Threshold.ThreshMeth threshMeth, Threshold.NoiseEstMeth noiseEstMeth);
    abstract void decompose(N data, int recursionLevel);
    public ArrayList<N> getDecomposition() {
        return waveletData;
    }
    public N getOriginalData() {
        return originalData;
    }
    public N getScalingData() {
        return scalingData;
    }
    public N getFilteredData() {
        return filteredData;
    }
    long nextPwr2(int n) {
        double logn = Math.log(n) / Math.log(2);
        return (long)Math.pow(2,(int)Math.ceil(logn));
    }
    abstract N AFB(N y, V filter);
    abstract N SFB(N lo, N hi, V g0, V g1);

}
    
