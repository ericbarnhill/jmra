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
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import org.apache.commons.math4.stat.descriptive.moment.Mean;

/** Undecimated 3D multi-resolution analysis. */
public class MRA3DU extends MRA3D {

     MRA1DU mra1du;

     public MRA3DU(ConvolverFactory.ConvolutionType convType) {
         super();
         mra1du = new MRA1DU(convType);
     }

    public MRA3DU(double[][][] origData, boolean[][][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        mra1du = new MRA1DU(convType);
        paddedData = origData;
    }

    public MRA3DU(double[][][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length, origData[0][0].length), filterBank, decompLvls, convType);
    }

    @Override
    public double[][][] analysis(double[][][] data, double[] filter, int decompositionLevel) {
        final int fi = data.length;
        final int fj = data[0].length;
        double[][][] filtData = new double[fi][fj][];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                filtData[i][j] = mra1du.analysis(data[i][j], filter, decompositionLevel);
                //filtData[i][j] = mra1d.analysis(data[i][j], filter, decompositionLevel);
            }
        }
        return filtData;
    }
    
    @Override
    public double[][][] synthesis(double[][][] lo, double[][][] hi, double[] sfl, double[] sfh, int decompositionLevel) {
        final int fi = lo.length;
        final int fj = lo[0].length;
        double[][][] y = new double[fi][fj][];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                y[i][j] = mra1du.synthesis(lo[i][j], hi[i][j], sfl, sfh, decompositionLevel);
                //y[i][j] = mra1d.synthesis(lo[i][j], hi[i][j], sfl, sfh, decompositionLevel);
            }
        }
        return y;
    }

    @Override
    // need to override so data is not cropped
    public double[][][] getFilteredData() {
        return waveletData.get(0);
    }

}
