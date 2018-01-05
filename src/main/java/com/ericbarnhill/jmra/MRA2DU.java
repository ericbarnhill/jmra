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
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import com.ericbarnhill.jmra.filters.*;

/**
 * Undecimated 2D multi-resolution analysis.
 */
public class MRA2DU extends MRA2D {

    MRA1DU mra1du;

    public MRA2DU(ConvolverFactory.ConvolutionType convType) {
        super(convType);
        mra1du = new MRA1DU(convType);
    }

    public MRA2DU(double[][] origData, boolean[][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        mra1du = new MRA1DU(convType);
    }

    public MRA2DU(double[][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length), filterBank, decompLvls, convType);
    }

    @Override
    public double[][] analysis(double[][] data, double[] filter, int decompositionLevel) {
        final int J = decompLvls - decompositionLevel;
        final int fi = data.length;
        final int fj = data[0].length;
        final int N = filter.length;
        double[][] filtData = new double[fi][];
        for (int i = 0; i < fi; i++) {
            //filtData[i] = mra1du.analysis(data[i], filter, J);
            filtData[i] = mra1d.analysis(data[i], filter, J);
        }
        return filtData;
    }

    @Override
    public double[][] synthesis(double[][] lo, double[][] hi, double[] sfl, double[] sfh, int decompositionLevel) {
        final int J = decompLvls - decompositionLevel;
        final int N = sfl.length + sfh.length;
        final int fi = lo.length;
        final int fj = lo[0].length;
        double[][] y = new double[fi][];
        for (int i = 0; i < fi; i++) {
            //y[i] = mra1du.synthesis(lo[i], hi[i], sfl, sfh, J);
            y[i] = mra1d.synthesis(lo[i], hi[i], sfl, sfh, J);
        }
        return y;
    }

    public void accept(Threshold threshold) {
        threshold.visit(this);
    }
}
