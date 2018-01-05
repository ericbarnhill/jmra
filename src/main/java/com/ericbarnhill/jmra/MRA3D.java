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

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;

/** 3D multi-resolution analysis */
public class MRA3D extends MRA<double[][][], boolean[][][], double[]> {

     public int w;
     public int h;
     public int d;
     public int area;
     public int volume;
     public int wPad;
     public int hPad;
     public int dPad;
     public int areaPad;
     public int volumePad;
     public double[][][] paddedData;
     public boolean[][][] paddedMask;
     public MRA1D mra1d;

     public MRA3D() {
         super();
     }

    public MRA3D(double[][][] origData, boolean[][][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        this.w = origData.length;
        this.h = origData[0].length;
        this.d = origData[0][0].length;
        this.area = w*h;
        this.volume = w*h*d;
        this.wPad = (int)nextPwr2(w);
        this.hPad = (int)nextPwr2(h);
        this.dPad = (int)nextPwr2(d);
        this.areaPad = wPad * hPad;
        this.volumePad = wPad * hPad * dPad;
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.dimLvls = 3;
        this.stride = 8;
        mra1d = new MRA1D(convType);
        initializeWaveletData();
    }

    public MRA3D(double[][][] origData, boolean[][][] maskData, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, decompLvls, convType);
        this.w = origData.length;
        this.h = origData[0].length;
        this.d = origData[0][0].length;
        this.area = w*h;
        this.volume = w*h*d;
        this.wPad = (int)nextPwr2(w);
        this.hPad = (int)nextPwr2(h);
        this.dPad = (int)nextPwr2(d);
        this.areaPad = wPad * hPad;
        this.volumePad = wPad * hPad * dPad;
        this.paddedData = ArrayMath.zeroPadBoundaries(origData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.paddedMask = ArrayMath.zeroPadBoundaries(maskData, (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
        this.dimLvls = 3;
        this.stride = 8;
        mra1d = new MRA1D(convType);
        initializeWaveletData();
    }

    public MRA3D(double[][][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length, origData[0][0].length), filterBank, decompLvls, convType);
    }

    void initializeWaveletData() {
        for (int n = 0; n < stride * decompLvls; n++) {
            waveletData.add(new double[0][][]);
        }
    }

    @Override
    public ArrayList<double[][][]> getDecomposition(int localIndex, int ind, int decompLvl, int dimLvl, int localStride) {
        double[][][] x = new double[0][][];
        if (dimLvl == 0) {
            if (decompLvl == 0) {
                x = ArrayMath.deepCopy(paddedData);
            } else {
                x = getData(localIndex - stride);
            }
        } else {
            x = getData(ind);
        }
        switch(localStride) {
            case 4: 
                x = ArrayMath.shiftDim(x, 1);
                break;
            case 2:
                x = ArrayMath.shiftDim(x, 2);
                break;
        }
        double[][][] lo = analysis(x, fb.af.lo, decompLvl);
        double[][][] hi = analysis(x, fb.af.hi, decompLvl);
        switch(localStride) {
            case 4: 
                lo = ArrayMath.shiftDim(lo, 2);
                hi = ArrayMath.shiftDim(hi, 2);
                break;
            case 2:
                lo = ArrayMath.shiftDim(lo, 1);
                hi = ArrayMath.shiftDim(hi, 1);
                break;
        }
        ArrayList<double[][][]> loAndHi = new ArrayList<double[][][]>();
        loAndHi.add(lo);
        loAndHi.add(hi);
        return loAndHi;
    }

    @Override
    public double[][][] getRecomposition(int localPair, int ind, int decompLvl, int dimLvl, int localStride) {
        double[][][] lo = getData(ind);
        double[][][] hi = getData(ind + localPair);
        switch (localStride) {
            case 4: 
                lo = ArrayMath.shiftDim(lo, 1);
                hi = ArrayMath.shiftDim(hi, 1);
                break;
            case 2:
                lo = ArrayMath.shiftDim(lo, 2);
                hi = ArrayMath.shiftDim(hi, 2);
                break;
        }
        double[][][] y = synthesis(lo, hi, fb.sf.lo, fb.sf.hi, decompLvl);
        switch (localStride) {
            case 4: 
                y = ArrayMath.shiftDim(y, 2);
                break;
            case 2:
                y = ArrayMath.shiftDim(y, 1);
                break;
        }
        return y;
    }
    
    @Override
    public double[][][] analysis(double[][][] data, double[] filter, int decompLvl) {
        final int fi = data.length;
        final int fj = data[0].length;
        final int fk = data[0][0].length;
        final int fk2 = fk / 2;
        double[][][] filtData = new double[fi][fj][];
        for (int i = 0; i < fi; i++) { 
            for (int j = 0; j < fj; j++) {
                filtData[i][j] = mra1d.analysis(data[i][j], filter, decompLvl); 
            } 
        }
        return filtData;
    } 

    @Override
    public double[][][] synthesis(double[][][] lo, double[][][] hi, double[] sfl, double[] sfh, int decompLvl) {
        final int fi = lo.length;
        final int fj = lo[0].length;
        final int fk = lo[0][0].length*2;
        double[][][] y = new double[fi][fj][];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                y[i][j] = mra1d.synthesis(lo[i][j], hi[i][j], sfl, sfh, decompLvl);
            }
        }
        return y;
    }

    public double[][][] getData(int index) {
        return waveletData.get(index);
    }

    public void setData(int index, double[][][] data) {
        waveletData.set(index, data);
    }

    public void accept(Threshold threshold) {
        threshold.visit(this);
    }

    @Override
    public double[][][] getFilteredData() {
        return ArrayMath.stripBorderPadding(getData(0), (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
    }

}
