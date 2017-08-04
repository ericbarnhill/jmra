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
     int decompLvls; // also used for stride
     int dimLvls;
     int stride;
     FilterBank fb;
     ConvolverFactory.ConvolutionType convType;
     UpFirDn upFirDn;

     public MRA() {
     }

     public MRA(ConvolverFactory.ConvolutionType convType) {
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

    final void dwt() {
        for (int decompLvl = 0; decompLvl < decompLvls; decompLvl++) {
            decompose(decompLvl, 0);
        }
    }

    final void idwt() {
        for (int decompLvl = decompLvls-1; decompLvl >= 0; decompLvl--) {
            recompose(decompLvl, dimLvls-1);
        }
    }

    abstract void decompose(int decompLvl, int dimLvl);

    abstract void recompose(int decompLvl, int dimLvl);

    abstract N AFB(N y, V filter, int decompLvl);

    abstract N SFB(N lo, N hi, V sfl, V sfh, int decompLvl); 

    abstract void accept(Threshold threshold);

    public N getOriginalData() {
        return origData;
    }

    public N getFilteredData() {
        return waveletData.get(0);
    }

    long nextPwr2(int n) {
        double logn = Math.log(n) / Math.log(2);
        return (long)Math.pow(2,(int)Math.ceil(logn));
    }
}
    
