/* 
 * Copyright (C) 2018 Eric Barnhill
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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

/**
 * MRA is an abstract base class for multi-resolution analysis.
 * JMRA uses a Template design pattern, with concrete classes
 * for number of dimensions  and type of transform (decimated 
 * or undecimated, etc.) Operations such as thresholding and 
 * serializing are performed with Visitor patterns.
 */ 
public abstract class MRA<N, B, V> {
    // N is ND array of numeric type
    // B is ND boolean array
    // V is 1D vector array of numeric type
     public N origData;
     public B maskData;
     public ArrayList<N> waveletData;
     public int decompLvls; // also used for stride
     public int dimLvls;
     public int stride;
     public FilterBank fb;
     public ConvolverFactory.ConvolutionType convType;
     public UpFirDn upFirDn;


     /** Empty constructor. */
     public MRA() {
     }

     /** Constructor containing only convolution operator. */
     public MRA(ConvolverFactory.ConvolutionType convType) {
        this.convType =  convType;
        upFirDn = new UpFirDn(convType); 
     }

     /** for extensions that require a filter bank with different specifications */
    public MRA(N origData, B maskData, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this.origData = origData;
        this.maskData = maskData;
        this.fb = fb;
        this.decompLvls = decompLvls;
        waveletData = new ArrayList<N>();
        this.convType =  convType;
        upFirDn = new UpFirDn(convType); 
    } 

     /** full constructor, containing image, mask, filter bank, number of levels 
      * of the decomposition, and convolution type (CPU, GPU, or FFT)
      */
    public MRA(N origData, B maskData, FilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this.origData = origData;
        this.maskData = maskData;
        this.fb = fb;
        this.decompLvls = decompLvls;
        waveletData = new ArrayList<N>();
        this.convType =  convType;
        upFirDn = new UpFirDn(convType); 
    } 

    /** perform discrete wavelet transform on image data */
    final public void dwt() {
        for (int decompLvl = 0; decompLvl < decompLvls; decompLvl++) {
            decompose(decompLvl, 0);
        }
    }

    /** perform inverse wavelet transform on image data */
    final public void idwt() {
        for (int decompLvl = decompLvls-1; decompLvl >= 0; decompLvl--) {
            recompose(decompLvl, dimLvls-1);
        }
    }

    /** implementation of recursive decomposition */
    public void decompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl);
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl; // starting point
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            ArrayList<N> loAndHi = getDecomposition(localIndex, ind, decompLvl, dimLvl, localStride);
            N lo = loAndHi.get(0);
            N hi = loAndHi.get(1);
            setData(ind, lo);
            setData(ind + localPair, hi);
        }
        if (dimLvl < dimLvls - 1) {
            decompose(decompLvl, dimLvl+1);
       }
    }

    /** specifics of the decomposition are left abstract */
    abstract public ArrayList<N> getDecomposition(int localIndex, int ind, int decompLvl, int dimLvl, int localStride);

    /** implementation of recursive recomposition */
    public void recompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl);
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl;
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            N y = getRecomposition(localPair, ind, decompLvl, dimLvl, localStride);
            setData(ind, y);
        }
        if (dimLvl > 0) {
            recompose(decompLvl, dimLvl-1);
        }
        if (decompLvl > 0) {
            waveletData.set(stride*(decompLvl-1), waveletData.get(stride*decompLvl));
        }
    }

    /** specifics of recomposition are left abstract */
    abstract public N getRecomposition(int localPair, int ind, int decompLvl, int dimLvl, int localStride);

    /** abstract method for the specific wavelet analysis used in the 
     * decomposition
     */
    abstract public N analysis(N y, V filter, int decompLvl);

    /** abstract method for the specific wavelet synthesis used in the
     * recomposition
     */
    abstract public N synthesis(N lo, N hi, V sfl, V sfh, int decompLvl); 

    /** returns the image data of this MRA */
    abstract public N getData(int index);

    /** sets the image data of this MRA */
    abstract public void setData(int index, N data);

    abstract public void accept(Threshold threshold);

    abstract public N getFilteredData();

    public ArrayList<N> getDecomposition() {
        return waveletData;
    }

    public long nextPwr2(int n) {
        double logn = Math.log(n) / Math.log(2);
        return (long)Math.pow(2,(int)Math.ceil(logn));
    }
}
    
