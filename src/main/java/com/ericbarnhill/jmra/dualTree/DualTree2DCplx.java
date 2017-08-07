package com.ericbarnhill.jmra.dualTree;

import com.ericbarnhill.arrayMath.*;
import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;

public class DualTree2DCplx extends DualTree<double[][], boolean[][], double[]> {


    private final int stride;

    public DualTree2DCplx(double[][] origData, boolean[][] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, fb, decompLvls, convType);
        this.stride = 4;
    }

    public DualTree2DCplx(double[][] origData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length, origData[0].length) , fb, decompLvls, convType);
    }

    public void setTrees() {
        int[][] bankIndices = { {0, 0}, {1, 0}, {0, 0}, {1, 0}, {0, 1}, {1, 1}, {0, 1}, {1, 1} };
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
        for (DTFilterBank bank : banks) {
            trees.add(new MRA2DDT(origData, maskData, bank, decompLvls, convType));
        }
    }

    public double[][] getFilteredData() {
        double[][] bankSum = new double[origData.length][origData[0].length];
        for (MRA<double[][], boolean[][], double[]> tree : trees) {
            bankSum = ArrayMath.add(bankSum, tree.getFilteredData());
        }
        bankSum = ArrayMath.divide(bankSum, 2);
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
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images
                ArrayList<double[][]> pm1 = pm(trees.get(0).waveletData.get(i), 
                        trees.get(3).waveletData.get(i));
                trees.get(0).waveletData.set(i, pm1.get(0));
                trees.get(3).waveletData.set(i, pm1.get(1));
                ArrayList<double[][]> pm2 = pm(trees.get(1).waveletData.get(i), 
                        trees.get(2).waveletData.get(i));
                trees.get(1).waveletData.set(i, pm2.get(0));
                trees.get(2).waveletData.set(i, pm2.get(1));
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
		ArrayList<double[][]> pm = new ArrayList<double[][]>();
		pm.add(p);
		pm.add(m);
		return pm;
	}
}
