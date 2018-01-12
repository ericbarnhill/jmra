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

public class DualTree3DCplx extends DualTree<double[][][], boolean[][][], double[]> {

    int stride;

    public DualTree3DCplx(double[][][] origData, boolean[][][] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        super(origData, maskData, fb, decompLvls, convType, undecimated);
        this.stride = 8;
        origData = ArrayMath.divide(origData, Math.sqrt(8));
    }

    public DualTree3DCplx(double[][][] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        this(origData, ArrayMath.fillWithTrue(origData.length, origData[0].length, origData[0][0].length), fb, decompLvls, convType, undecimated);
    }

    public DualTree3DCplx(double[][][] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length, origData[0].length, origData[0][0].length), fb, decompLvls, convType, true);
    }

    public DualTree3DCplx(double[][][] origData) {
        this(origData, ArrayMath.fillWithTrue(origData.length, origData[0].length, origData[0][0].length), Wavelets.getFarrasKingsbury(), 3, ConvolverFactory.ConvolutionType.FDCPU, true);
    }

    @Override
    public void setTrees() {
        int[][] bankIndices = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1} };
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
                trees.add(new MRA3DDTU(origData, maskData, bank, decompLvls, convType));
            } else {
                trees.add(new MRA3DDT(origData, maskData, bank, decompLvls, convType));
            }
        }
    }

    public double[][][] getFilteredData() {
        double[][][] bankSum = new double[origData.length][origData[0].length][origData[0][0].length];
        for (MRA<double[][][], boolean[][][], double[]> tree : trees) {
            bankSum = ArrayMath.add(bankSum, tree.getFilteredData());
        }
        bankSum = ArrayMath.divide(bankSum, Math.sqrt(8));
        //bankSum = ArrayMath.divide(bankSum, 16);
        return bankSum;
    }

    public void dwt() { 
        setTrees();
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).dwt();
        }
        addSubtract(true);
    }

    public void idwt() {
        addSubtract(false);
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).idwt();
        }
    }

    public void accept(Threshold threshold) {
        threshold.visit(this);
    }

    public void addSubtract(boolean fwd) {
        if (fwd) {
            addSubtractFwd();
        } else {
            addSubtractInv();
        }
    }

    void addSubtractFwd() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images
                ArrayList<double[][][]> pm1 = pm4(
                        trees.get(0).getData(i), 
                        trees.get(5).getData(i),
                        trees.get(3).getData(i),
                        trees.get(6).getData(i));

                trees.get(0).setData(i, pm1.get(0));
                trees.get(5).setData(i, pm1.get(1));
                trees.get(3).setData(i, pm1.get(2));
                trees.get(6).setData(i, pm1.get(3));

                ArrayList<double[][][]> pm2 = pm4(
                        trees.get(7).getData(i), 
                        trees.get(2).getData(i),
                        trees.get(4).getData(i),
                        trees.get(1).getData(i));

                trees.get(7).setData(i, pm1.get(0));
                trees.get(2).setData(i, pm1.get(1));
                trees.get(4).setData(i, pm1.get(2));
                trees.get(1).setData(i, pm1.get(3));

            }
        }
    }

    void addSubtractInv() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images

                ArrayList<double[][][]> pm1 = pm4inv(
                        trees.get(0).getData(i), 
                        trees.get(5).getData(i),
                        trees.get(3).getData(i),
                        trees.get(6).getData(i));

                trees.get(0).setData(i, pm1.get(0));
                trees.get(5).setData(i, pm1.get(1));
                trees.get(3).setData(i, pm1.get(2));
                trees.get(6).setData(i, pm1.get(3));

                ArrayList<double[][][]> pm2 = pm4inv(
                        trees.get(7).getData(i), 
                        trees.get(2).getData(i),
                        trees.get(4).getData(i),
                        trees.get(1).getData(i));

                trees.get(7).setData(i, pm1.get(0));
                trees.get(2).setData(i, pm1.get(1));
                trees.get(4).setData(i, pm1.get(2));
                trees.get(1).setData(i, pm1.get(3));

            }
        }
    }

	static ArrayList<double[][][]> pm4 (double[][][] a, double[][][] b,
			double[][][] c, double[][][] d) {
		 double[][][] p = 	
				ArrayMath.divide(
					ArrayMath.subtract(
						ArrayMath.subtract(
							ArrayMath.subtract(
								ArrayMath.deepCopy(a),
								ArrayMath.deepCopy(b)
							), ArrayMath.deepCopy(c)
						), ArrayMath.deepCopy(d)
					), 2
			);
         double[][][] q = 
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.add(
							ArrayMath.subtract(
								ArrayMath.deepCopy(a),
								ArrayMath.deepCopy(b)
							), ArrayMath.deepCopy(c)
						), ArrayMath.deepCopy(d)
					), 2
			);
         double[][][] r = 
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.subtract(
							ArrayMath.add(
								ArrayMath.deepCopy(a),
								ArrayMath.deepCopy(b)
							), ArrayMath.deepCopy(c)
						), ArrayMath.deepCopy(d)
					), 2
			);
         double[][][] s = 
				ArrayMath.divide(
					ArrayMath.subtract(
						ArrayMath.add(
							ArrayMath.add(
								ArrayMath.deepCopy(a),
								ArrayMath.deepCopy(b)
							), ArrayMath.deepCopy(c)
						), ArrayMath.deepCopy(d)
					), 2
			);
		ArrayList<double[][][]> pm = new ArrayList<double[][][]>();
		pm.add(p);
		pm.add(q);
		pm.add(r);
		pm.add(s);
		return pm;
	}

	static ArrayList<double[][][]> pm4inv (double[][][] a, double[][][] b,
			double[][][] c, double[][][] d) {
        double[][][] p = 
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.add(
							ArrayMath.add(
								ArrayMath.deepCopy(a),
								ArrayMath.deepCopy(b)
							), ArrayMath.deepCopy(c)
						), ArrayMath.deepCopy(d)
					), 2
			);
        double[][][] q = 
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.add(
							ArrayMath.subtract(
								ArrayMath.multiply(ArrayMath.deepCopy(a), -1),
								ArrayMath.deepCopy(b)
							), ArrayMath.deepCopy(c)
						), ArrayMath.deepCopy(d)
					), 2
			);
        double[][][] r = 
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.subtract(
							ArrayMath.add(
								ArrayMath.multiply(ArrayMath.deepCopy(a), -1),
								ArrayMath.deepCopy(b)
							), ArrayMath.deepCopy(c)
						), ArrayMath.deepCopy(d)
					), 2
			);
        double[][][] s = 
				ArrayMath.divide(
					ArrayMath.subtract(
						ArrayMath.add(
							ArrayMath.add(
								ArrayMath.multiply(ArrayMath.deepCopy(a), -1),
								ArrayMath.deepCopy(b)
							), ArrayMath.deepCopy(c)
						), ArrayMath.deepCopy(d)
					), 2
			);
		ArrayList<double[][][]> pm = new ArrayList<double[][][]>();
		pm.add(p);
		pm.add(q);
		pm.add(r);
		pm.add(s);
		return pm;
	}
}
