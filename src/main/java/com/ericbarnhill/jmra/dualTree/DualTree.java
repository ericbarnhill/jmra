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
import com.ericbarnhill.jmra.filters.*;

/** 
 * Abstract class for implementation of Dual-Tree MRA. 
 */ 
public abstract class DualTree<N, B, V> {

    public ArrayList<MRA<N, B, V>> trees;
    public ArrayList<DTFilterBank> banks;

    N origData;
    B maskData;
    DTFilterBank fb;
    int decompLvls;
    public ConvolverFactory.ConvolutionType convType;
    public boolean undecimated;

    /** 
     * Base constructor. Note that the DualTree takes a DTFilterBank 
     * instead of the FilterBank used in the MRA class. 
     */
    public DualTree(N origData, B maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        this.origData = origData;
        this.maskData = maskData;
        this.fb = fb;
        this.decompLvls = decompLvls;
        this.convType = convType;
        this.undecimated = undecimated;
        trees = new ArrayList<MRA<N, B, V>>();
        banks = new ArrayList<DTFilterBank>();
    }

    /**
     * Build trees out of filter bank filters.
     */
    public abstract void setTrees();

    /** Discrete wavelet transform. */
    public abstract void dwt();

    /** Inverse discrete wavelet transform. */
    public abstract void idwt();

    /** Normalizes the dual-tree transforms.  */
    public abstract void addSubtract(boolean fwd);

    /** Returns filtered data from the DualTree. */
    public abstract N getFilteredData();

    /** Returns original data from the DualTree. */
    public N getOrigData() {
        return origData;
    }

    /** Returns wavelet decompositions from the trees */
    public ArrayList<ArrayList<N>> getDecomposition() {
        ArrayList<ArrayList<N>> decomposition = new ArrayList<ArrayList<N>>();
        for (int i = 0; i < trees.size(); i++) {
            decomposition.add(trees.get(i).getDecomposition());
        }
        return decomposition;
    }

}
