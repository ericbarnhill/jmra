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

package com.ericbarnhill.jmra.dualTree;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

/** 3D MRA within a 3D dual-tree analysis */
public class MRA3DDT extends MRA3D {

    DTFilterBank fb;

    public MRA3DDT(double[][][] origData, boolean[][][] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, decompLvls, convType);
        this.fb = fb;
    }

    public MRA3DDT(double[][][] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length, origData[0][0].length), fb, decompLvls, convType);
    }

    @Override
    public ArrayList<double[][][]> getDecomposition(int localIndex, int ind, int decompLvl, int dimLvl, int localStride) {
        double[][][] x = new double[0][][];
        if (dimLvl == 0) {
            if (decompLvl == 0) {
                x = ArrayMath.deepCopy(paddedData);
            } else {
                x = waveletData.get(localIndex - stride);
            }
        } else {
            x = waveletData.get(ind);
        }
        double[][][] lo = new double[0][][];
        double[][][] hi = new double[0][][];
        switch(localStride) {
            case 8: 
                if (decompLvl == 0) {
                    lo = analysis(x, fb.faf.get(0).lo, decompLvl);
                    hi = analysis(x, fb.faf.get(0).hi, decompLvl);
                } else {
                    lo = analysis(x, fb.af.get(0).lo, decompLvl);
                    hi = analysis(x, fb.af.get(0).hi, decompLvl);
                }    
                break;
            case 4: 
                x = ArrayMath.shiftDim(x, 1);
                if (decompLvl == 0) {
                    lo = analysis(x, fb.faf.get(1).lo, decompLvl);
                    hi = analysis(x, fb.faf.get(1).hi, decompLvl);
                } else {
                    lo = analysis(x, fb.af.get(1).lo, decompLvl);
                    hi = analysis(x, fb.af.get(1).hi, decompLvl);
                }    
                lo = ArrayMath.shiftDim(lo, 2);
                hi = ArrayMath.shiftDim(hi, 2);
                break;
            case 2:
                x = ArrayMath.shiftDim(x, 2);
                if (decompLvl == 0) {
                    lo = analysis(x, fb.faf.get(2).lo, decompLvl);
                    hi = analysis(x, fb.faf.get(2).hi, decompLvl);
                } else {
                    lo = analysis(x, fb.af.get(2).lo, decompLvl);
                    hi = analysis(x, fb.af.get(2).hi, decompLvl);
                }    
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
        double[][][] lo = waveletData.get(ind);
        double[][][] hi = waveletData.get(ind + localPair);
        double[][][] y = new double[0][][]; 
        switch (localStride) {
            case 8:
                if (decompLvl == 0) {
                    y = synthesis(lo, hi, fb.fsf.get(0).lo, fb.fsf.get(0).hi, decompLvl);
                } else {
                    y = synthesis(lo, hi, fb.sf.get(0).lo, fb.sf.get(0).hi, decompLvl);
                }
                break;
            case 4: 
                lo = ArrayMath.shiftDim(lo, 1);
                hi = ArrayMath.shiftDim(hi, 1);
                if (decompLvl == 0) {
                    y = synthesis(lo, hi, fb.fsf.get(1).lo, fb.fsf.get(1).hi, decompLvl);
                } else {
                    y = synthesis(lo, hi, fb.sf.get(1).lo, fb.sf.get(1).hi, decompLvl);
                }
                y = ArrayMath.shiftDim(y, 2);
                break;
            case 2:
                lo = ArrayMath.shiftDim(lo, 2);
                hi = ArrayMath.shiftDim(hi, 2);
                if (decompLvl == 0) {
                    y = synthesis(lo, hi, fb.fsf.get(2).lo, fb.fsf.get(2).hi, decompLvl);
                } else {
                    y = synthesis(lo, hi, fb.sf.get(2).lo, fb.sf.get(2).hi, decompLvl);
                }
                y = ArrayMath.shiftDim(y, 1);
                break;
        }
        return y;
    }

}
