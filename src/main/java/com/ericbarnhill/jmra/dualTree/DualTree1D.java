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
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.arrayMath.*;

/** 1D dual-tree analysis. */ 
public class DualTree1D extends DualTree<double[], boolean[], double[]> {

    public DualTree1D(double[] origData, boolean[] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        super(origData, maskData, fb, decompLvls, convType, undecimated);
    }

    public DualTree1D(double[] origData, boolean[] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, maskData, fb, decompLvls, convType, false);
    }

    public void setTrees() {
        // tabled until I address the Complex DualTree first
        /*
        DTFilterBank tree1Bank = new DTFilterBank(fb.faf.lo, fb.faf.lo, fb.fsf.lo, fb.fsf.lo, fb.af1, fb.af1, fb.sf1, fb.sf1);
        DTFilterBank tree2Bank = new DTFilterBank(fb.faf.hi, fb.faf.hi, fb.fsf.hi, fb.fsf.hi, fb.af2, fb.af2, fb.sf2, fb.sf2);
        trees.add(new MRA1DDT(origData, maskData, tree1Bank, decompLvls, convType));
        trees.add(new MRA1DDT(origData, maskData, tree2Bank, decompLvls, convType));
        */
    }

    public double[] getFilteredData() {
        // to be addressed after complex dualtrees are up and running
        /*
        return ArrayMath.divide(ArrayMath.add(tree1Bank.getFilteredData(), tree2Bank.getFilteredData()), Math.sqrt(2));
        */
        return new double[0];
    }
    
    public void addSubtract(boolean fwd) {}

    public void dwt() {
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).dwt();
        }
    }

    public void idwt() {
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).idwt();
        }
    }
}
