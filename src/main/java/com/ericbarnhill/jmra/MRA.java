package com.ericbarnhill.jmra; 
import java.util.ArrayList;
import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

public abstract class MRA<N, B, V> {
    // N is ND array of numeric type
    // B is ND boolean array
    // V is 1D vector array of numeric type
     N origData;
     B maskData;
     public ArrayList<N> waveletData; // sometimes manipulated and put back
     public int decompLvls; // also used for stride
     public int dimLvls;
     public int stride;
     public FilterBank fb;
     public ConvolverFactory.ConvolutionType convType;
     public UpFirDn upFirDn;

     public MRA() {
     }

     public MRA(ConvolverFactory.ConvolutionType convType) {
        this.convType =  convType;
        upFirDn = new UpFirDn(convType); 
     }

     // for extensions that require a filter bank with different specifications
    public MRA(N origData, B maskData, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this.origData = origData;
        this.maskData = maskData;
        this.fb = fb;
        this.decompLvls = decompLvls;
        waveletData = new ArrayList<N>();
        this.convType =  convType;
        upFirDn = new UpFirDn(convType); 
    } 

    public MRA(N origData, B maskData, FilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this.origData = origData;
        this.maskData = maskData;
        this.fb = fb;
        this.decompLvls = decompLvls;
        waveletData = new ArrayList<N>();
        this.convType =  convType;
        upFirDn = new UpFirDn(convType); 
    } 

    final public void dwt() {
        for (int decompLvl = 0; decompLvl < decompLvls; decompLvl++) {
            decompose(decompLvl, 0);
        }
    }

    final public void idwt() {
        for (int decompLvl = decompLvls-1; decompLvl >= 0; decompLvl--) {
            recompose(decompLvl, dimLvls-1);
        }
    }

    abstract public void decompose(int decompLvl, int dimLvl);

    abstract public void recompose(int decompLvl, int dimLvl);

    abstract public N AFB(N y, V filter, int decompLvl);

    abstract public N SFB(N lo, N hi, V sfl, V sfh, int decompLvl); 

    abstract public void accept(Threshold threshold);

    public N getOriginalData() {
        return origData;
    }

    public N getFilteredData() {
        return waveletData.get(0);
    }

    public ArrayList<N> getDecomposition() {
        return waveletData;
    }

    public long nextPwr2(int n) {
        double logn = Math.log(n) / Math.log(2);
        return (long)Math.pow(2,(int)Math.ceil(logn));
    }
}
    
