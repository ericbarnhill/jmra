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

import java.util.ArrayList; 
import com.ericbarnhill.arrayMath.*;
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;

/** 3D complex dual-tree analysis with serial data storage. */
public class DualTree3DCplxSerial extends DualTree3DCplx {

    public DualTree3DCplxSerial(double[][][] origData, boolean[][][] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        super(origData, maskData, fb, decompLvls, convType, undecimated);
        this.stride = 8;
    }

    public DualTree3DCplxSerial(double[][][] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        this(origData, ArrayMath.fillWithTrue(origData.length, origData[0].length, origData[0][0].length), fb, decompLvls, convType, undecimated);
    }

    @Override
    public void setTrees() {
        int[][] bankIndices = { {0, 0, 0}, {1, 0, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 0}, {1, 1, 0}, {0, 1, 1}, {1, 1, 1} };
        for (int[] indices : bankIndices) {
            ArrayList<FilterPair> faf = new ArrayList<FilterPair>();
            ArrayList<FilterPair> fsf = new ArrayList<FilterPair>();
            ArrayList<FilterPair> af = new ArrayList<FilterPair>();
            ArrayList<FilterPair> sf = new ArrayList<FilterPair>();
            for (int j = 0; j < 3; j++) {
                faf.add(fb.faf.get(indices[j]));
                fsf.add(fb.fsf.get(indices[j]));
                af.add(fb.af.get(indices[j]));
                sf.add(fb.sf.get(indices[j]));
            }
            banks.add(new DTFilterBank(faf, fsf, af, sf));
        }
        for (DTFilterBank bank : banks) {
            if (undecimated) {
                trees.add(new MRA3DDTUSerial(origData, maskData, bank, decompLvls, convType));
            } else {
                trees.add(new MRA3DDTSerial(origData, maskData, bank, decompLvls, convType));
            }
        }
    }

    @Override
    public void accept(Threshold threshold) {
        threshold.visitSerial(this);
    }
}
