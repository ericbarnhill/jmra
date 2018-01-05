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

/** Undecimated 1D MRA within a dual-tree analysis. */
public class MRA1DDTU extends MRA1DDT {

    MRA1DU mra1du;

    public MRA1DDTU(double[] origData, boolean[] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, fb, decompLvls, convType);
        mra1du = new MRA1DU(convType);
    }

    public MRA1DDTU(double[] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length), fb, decompLvls, convType);
    }

    @Override
    public double[] analysis(double[] x, double[] filter, int decompLvl) {
        return mra1du.analysis(x, filter, decompLvl);
    }

    @Override
    public double[] synthesis(double[] lo, double[] hi, double[] sfl, double[] sfh, int decompLvl) {
        return mra1du.synthesis(lo, hi, sfl, sfh, decompLvl);
    }

}
