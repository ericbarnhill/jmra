package com.ericbarnhill.jmra; 
import java.util.ArrayList;
import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

public abstract class MRA<N, B, V> {
    // N is ND array of numeric type
    // B is ND boolean array
    // V is 1D vector array of numeric type
     N originalData;
     N scalingData;
     B maskData;
     ArrayList<N> waveletData;
     N filteredData;
     int decompLvls; // also used for stride
     int dimLvls;
     int stride;
     ArrayList<ArrayList<V>> filterBank; // order af, afh, sfl, sfh
     ArrayList<V> analysisFilters;
     V afl;
     V afh;
     V sfl;
     V sfh;
     ArrayList<V> synthesisFilters;
     ConvolverFactory.ConvolutionType convolutionType;
     UpFirDn upFirDn;

     public MRA() {
     }

     public MRA(ConvolverFactory.ConvolutionType convolutionType) {
        this.convolutionType =  convolutionType;
        upFirDn = new UpFirDn(convolutionType); 
     }

    public MRA(N originalData, B maskData, ArrayList<ArrayList<V>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convolutionType) {
        this.originalData = originalData;
        this.scalingData = originalData; 
        this.maskData = maskData;
        this.filterBank = filterBank;
        this.analysisFilters = filterBank.get(0);
        this.afl = analysisFilters.get(0);
        this.afh = analysisFilters.get(1);
        this.synthesisFilters = filterBank.get(1);
        this.sfl = synthesisFilters.get(0);
        this.sfh = synthesisFilters.get(1);
        this.decompLvls = decompLvls;
        waveletData = new ArrayList<N>();
        this.convolutionType =  convolutionType;
        upFirDn = new UpFirDn(convolutionType); 
    } 

    void dwt() {
        for (int decompLvl = 0; decompLvl < decompLvls; decompLvl++) {
            decompose(decompLvl, 0);
        }
    }

    void idwt() {
        for (int decompLvl = decompLvls-1; decompLvl >= 0; decompLvl--) {
            recompose(decompLvl, dimLvls-1);
        }
    }

    abstract void decompose(int decompLvl, int dimLvl);
    abstract void recompose(int decompLvl, int dimLvl);

    abstract N AFB(N y, V filter, int decompLvl);
    abstract N SFB(N lo, N hi, V sfl, V sfh, int decompLvl); 

    public abstract void threshold(Threshold.ThreshMeth threshMeth, Threshold.NoiseEstMeth noiseEstMeth);

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
}
    
