package com.ericbarnhill.jmra.dualTree;

import java.util.ArrayList; 
import com.ericbarnhill.arrayMath.*;
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;

public class DualTree3DCplx extends DualTree<double[][][], boolean[][][], double[]> {


    private final int stride;

    public DualTree3DCplx(double[][][] origData, boolean[][][] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, fb, decompLvls, convType);
        this.stride = 8;
    }

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
            trees.add(new MRA3DDT(origData, maskData, bank, decompLvls, convType));
        }
    }

    public double[][][] getFilteredData() {
        double[][][] bankSum = new double[origData.length][origData[0].length][origData[0][0].length];
        for (MRA<double[][][], boolean[][][], double[]> tree : trees) {
            bankSum = ArrayMath.add(bankSum, tree.getFilteredData());
        }
        bankSum = ArrayMath.divide(bankSum, Math.sqrt(8));
        return bankSum;
    }

    public void dwt() {
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

    public void addSubtract(boolean fwd) {
        if (fwd) {
            addSubtractFwd();
        } else {
            addSubtractInv();
        }
    }

    private void addSubtractFwd() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images

                ArrayList<double[][][]> pm1 = pm4(
                        trees.get(0).waveletData.get(i), 
                        trees.get(5).waveletData.get(i),
                        trees.get(3).waveletData.get(i),
                        trees.get(6).waveletData.get(i));

                trees.get(0).waveletData.set(i, pm1.get(0));
                trees.get(5).waveletData.set(i, pm1.get(5));
                trees.get(3).waveletData.set(i, pm1.get(3));
                trees.get(6).waveletData.set(i, pm1.get(6));

                ArrayList<double[][][]> pm2 = pm4(
                        trees.get(8).waveletData.get(i), 
                        trees.get(2).waveletData.get(i),
                        trees.get(4).waveletData.get(i),
                        trees.get(1).waveletData.get(i));

                trees.get(8).waveletData.set(i, pm1.get(8));
                trees.get(2).waveletData.set(i, pm1.get(2));
                trees.get(4).waveletData.set(i, pm1.get(4));
                trees.get(1).waveletData.set(i, pm1.get(1));

            }
        }
    }

    private void addSubtractInv() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images

                ArrayList<double[][][]> pm1 = pm4inv(
                        trees.get(0).waveletData.get(i), 
                        trees.get(5).waveletData.get(i),
                        trees.get(3).waveletData.get(i),
                        trees.get(6).waveletData.get(i));

                trees.get(0).waveletData.set(i, pm1.get(0));
                trees.get(5).waveletData.set(i, pm1.get(5));
                trees.get(3).waveletData.set(i, pm1.get(3));
                trees.get(6).waveletData.set(i, pm1.get(6));

                ArrayList<double[][][]> pm2 = pm4inv(
                        trees.get(8).waveletData.get(i), 
                        trees.get(2).waveletData.get(i),
                        trees.get(4).waveletData.get(i),
                        trees.get(1).waveletData.get(i));

                trees.get(8).waveletData.set(i, pm1.get(8));
                trees.get(2).waveletData.set(i, pm1.get(2));
                trees.get(4).waveletData.set(i, pm1.get(4));
                trees.get(1).waveletData.set(i, pm1.get(1));

            }
        }
    }

	private static ArrayList<double[][][]> pm4 (double[][][] a, double[][][] b,
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

	private static ArrayList<double[][][]> pm4inv (double[][][] a, double[][][] b,
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
