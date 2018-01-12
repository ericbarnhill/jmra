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

import com.ericbarnhill.arrayMath.*;
import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;

/** 2D complex dual-tree analysis. */
public class DualTree2DCplx extends DualTree<double[][], boolean[][], double[]> {


    private final int stride;

    public DualTree2DCplx(double[][] origData, boolean[][] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        super(origData, maskData, fb, decompLvls, convType, undecimated);
        this.stride = 4;
    }

    public DualTree2DCplx(double[][] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        this(origData, ArrayMath.fillWithTrue(origData.length, origData[0].length) , fb, decompLvls, convType, undecimated);
    }

    public DualTree2DCplx(double[][] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length, origData[0].length) , fb, decompLvls, convType, false);
    }

    @Override
    public void setTrees() {
        int[][] bankIndices = { {0, 0}, {1, 0}, {0, 0}, {1, 0}, {0, 0}, {1, 0}, {0, 0}, {1, 0} };
        for (int[] indices : bankIndices) {
            ArrayList<FilterPair> faf = new ArrayList<FilterPair>();
            ArrayList<FilterPair> fsf = new ArrayList<FilterPair>();
            ArrayList<FilterPair> af = new ArrayList<FilterPair>();
            ArrayList<FilterPair> sf = new ArrayList<FilterPair>();
            for (int j = 0; j < 2; j++) {
                faf.add(fb.faf.get(indices[j]));
                fsf.add(fb.fsf.get(indices[j]));
                af.add(fb.af.get(indices[j]));
                sf.add(fb.sf.get(indices[j]));
            }
            banks.add(new DTFilterBank(faf, fsf, af, sf));
        }
        double[][] normData = ArrayMath.divide(origData, 2); 
        for (DTFilterBank bank : banks) {
            if (undecimated) {
                trees.add(new MRA2DDTU(normData, maskData, bank, decompLvls, convType));
            } else {
                trees.add(new MRA2DDT(normData, maskData, bank, decompLvls, convType));
            }
        }
    }

    public double[][] getFilteredData() {
        double[][] bankSum = new double[origData.length][origData[0].length];
        for (MRA<double[][], boolean[][], double[]> tree : trees) {
            bankSum = ArrayMath.add(bankSum, tree.getFilteredData());
        }
        //bankSum = ArrayMath.divide(bankSum, 2);
        bankSum = ArrayMath.divide(bankSum, 8);
        return bankSum;
    }

    public void dwt() {
        setTrees();
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).dwt();
        }
        addSubtract(true); // operation is the same in 2D, boolean added to fit interface
    }

    public void idwt() {
        addSubtract(false); // operation is the same in 2D, boolean added to fit interface
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).idwt();
        }
    }

    public void accept(Threshold threshold) {
        threshold.visit(this);
    }

    public void addSubtract(boolean fwd) {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images
                ArrayList<double[][]> pm1 = pm(trees.get(0).getData(i), 
                        trees.get(3).getData(i));
                trees.get(0).setData(i, pm1.get(0));
                trees.get(3).setData(i, pm1.get(1));
                ArrayList<double[][]> pm2 = pm(trees.get(1).getData(i), 
                        trees.get(2).getData(i));
                trees.get(1).setData(i, pm2.get(0));
                trees.get(2).setData(i, pm2.get(1));
            }
        }
    }
    
	public static ArrayList<double[][]> pm (double[][] u, double[][] v) {
        double[][] p = ArrayMath.divide(
                    ArrayMath.add(
                            ArrayMath.deepCopy(u),
                            ArrayMath.deepCopy(v)
                    ),
                    Math.sqrt(2)
                );
		double[][] m = ArrayMath.divide(
					ArrayMath.subtract(
						ArrayMath.deepCopy(u),
						ArrayMath.deepCopy(v)
					),
				    Math.sqrt(2)
                );
        //double[][] p = ArrayMath.fillWithRandom(u.length, u[0].length);
        //double[][] m = ArrayMath.fillWithRandom(u.length, u[0].length);
		ArrayList<double[][]> pm = new ArrayList<double[][]>();
		pm.add(p);
		pm.add(m);
		return pm;
	}
}
