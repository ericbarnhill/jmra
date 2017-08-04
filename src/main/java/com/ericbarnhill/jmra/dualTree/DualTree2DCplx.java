package com.ericbarnhill.jmra.dualTree;

import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;

public class DualTree2DCplx extends DualTree<double[][], boolean[][], double[]> {


    private final int stride;

    public DualTree2DCplx(double[][] origData, boolean[][] maskData, DTFilterSet2D fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, fb, decompLvls, convType);
        this.stride = 4;
    }

    public void setTrees() {
        int[][] bankIndices = { {0, 0}, {1, 0}, {0, 0}, {1, 0}, {0, 1}, {1, 1}, {0, 1}, {1, 1} };
        for (int[] indices : bankIndices) {
            ArrayList<FilterPair> faf = new ArrayList<FilterPair>();
            ArrayList<FilterPair> fsf = new ArrayList<FilterPair>();
            ArrayList<FilterPair> af = new ArrayList<FilterPair>();
            ArrayList<FilterPair> sf = new ArrayList<FilterPair>();
            for (int j = 0; j < 2; j++) {
                faf.add(fb.faf(indices[j]));
                fsf.add(fb.fsf(indices[j]));
                af.add(fb.af(indices[j]));
                sf.add(fb.sf(indices[j]));
            }
            banks.add(new DTFilterBank(faf, fsf, af, sf));
        }
        for (DTFilterBank bank : banks) {
            trees.add(new MRA2DDT(origData, maskData, bank, decompLvls, convType));
        }
    }

    public double[][] getFilteredData() {
        double[][] bankSum = ArrayMath.fillWithZeros(origData.length, origData[0].length);
        for (MRA2DDT tree : trees) {
            bankSum = ArrayMath.add(bankSum, tree.getFilteredData);
        }
        bankSum = ArrayMath.divide(bankSum, 2);
        return bankSum;
    }

    public void dwt() {
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).dwt();
        }
        addSubtract();
    }

    public void idwt() {
        addSubtract();
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).idwt();
        }
    }

    private void addSubtract() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images
                ArrayList<ArrayList<double[][]>> pm1 = pm(trees.get(0).waveletData.get(i), 
                        trees.get(3).waveletData.get(i));
                trees.get(0).waveletData.set(i, pm1.get(0));
                trees.get(3).waveletData.set(i, pm1.get(3));
                ArrayList<ArrayList<double[][]>> pm2 = pm(trees.get(1).waveletData.get(i), 
                        trees.get(2).waveletData.get(i));
                trees.get(1).waveletData.set(i, pm2.get(1));
                trees.get(2).waveletData.set(i, pm2.get(2));
            }
        }
    }
    
	public static ArrayList<ArrayList<double[][]>> pm (ArrayList<double[][]> u, ArrayList<double[][]> v) {
		ArrayList<double[][]> p = new ArrayList<double[][]>();
		ArrayList<double[][]> m = new ArrayList<double[][]>();
		// no operation on low pass
		p.add(ArrayMath.deepCopy(u.get(0)));
		m.add(ArrayMath.deepCopy(v.get(0)));
		// sum / diff on high pass
		for (int i = 1; i < u.size(); i++) {
			p.add(
				ArrayMath.divide(
					ArrayMath.add(
							ArrayMath.deepCopy(u.get(i)),
							ArrayMath.deepCopy(v.get(i))
					),
					Math.sqrt(2)
				)
			);
			m.add(
				ArrayMath.divide(
					ArrayMath.subtract(
						ArrayMath.deepCopy(u.get(i)),
						ArrayMath.deepCopy(v.get(i))
					),
					Math.sqrt(2)
				)
			);
		}
		ArrayList<ArrayList<double[][]>> pm = new ArrayList<ArrayList<double[][]>>();
		pm.add(p);
		pm.add(m);
		return pm;
	}
}
