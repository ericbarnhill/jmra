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
import java.util.ArrayList;
import java.util.Arrays;

/**
 * 1D undecimated multi-resolution analysis.
 */
public class MRA1DU extends MRA1D { 

     public MRA1DU(ConvolverFactory.ConvolutionType convType) {
         super(convType);
     }

     @Override
     public double[] analysis(double[] y, double[] filter, int decompLvl) {
        double[] filterATrous = ArrayMath.divide(aTrous(filter, decompLvl), Math.sqrt(2));
        y = upFirDn.upFirDn(y,filterATrous, 1, 1);
        return y;
     }

     @Override
     public double[] synthesis(double[] lo, double[] hi, double[] sfl, double[] sfh, int decompLvl) {
        double[] sflATrous = ArrayMath.divide(aTrous(sfl, decompLvl), Math.sqrt(2));
        double[] sfhATrous = ArrayMath.divide(aTrous(sfh, decompLvl), Math.sqrt(2));
        int M = (int)Math.pow(2, decompLvl);
        int N = sfl.length + sfh.length;
        lo = upFirDn.upFirDn(lo, sflATrous, 1, 1);
        hi = upFirDn.upFirDn(hi, sfhATrous, 1, 1);
        // in undec, crop the lo to fit the hi arrays copy of range
        lo = Arrays.copyOfRange(lo, 0, hi.length);
        double[] y = ArrayMath.add(lo, hi);
        try {
            y = Arrays.copyOfRange(y, M*(N/2-1), y.length - M*(N/2-1)); 
        } catch (Exception e) {
            throw new RuntimeException("Error in synthesis -- perhaps too many " +
                    "decomposition levels for the dimensions?");
        }
        
        return y;
    }

    private double[] aTrous(double[] filter, int decompLvl) {
        final int fi = filter.length;
        final int scaleFactor = (int)Math.pow(2, decompLvl);
        double[] aTrousFilter = new double[(fi-1)*scaleFactor + 1];
        for (int i = 0; i < fi; i++) {
            aTrousFilter[i*scaleFactor] = filter[i];
        }
        return aTrousFilter;
    }

    @Override
    public double[] getFilteredData() {
        return getData(0);
    }
}
