package com.ericbarnhill.jmra.dualTree;

import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;

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
                faf.add(fb.faf(indices[j]));
                fsf.add(fb.fsf(indices[j]));
                af.add(fb.af(indices[j]));
                sf.add(fb.sf(indices[j]));
            }
            banks.add(new DTFilterBank(faf, fsf, af, sf));
        }
        for (DTFilterBank bank : banks) {
            trees.add(new MRA3DDT(origData, maskData, bank, decompLvls, convType));
        }
    }

    public double[][] getFilteredData() {
        double[][] bankSum = ArrayMath.fillWithZeros(origData.length, origData[0].length);
        for (MRA3DDT tree : trees) {
            bankSum = ArrayMath.add(bankSum, tree.getFilteredData);
        }
        bankSum = ArrayMath.divide(bankSum, Math.sqrt(8));
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

    private void addSubtractFwd() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images

                ArrayList<ArrayList<double[][]>> pm1 = pm4(
                        trees.get(0).waveletData.get(i), 
                        trees.get(5).waveletData.get(i),
                        trees.get(3).waveletData.get(i),
                        trees.get(6).waveletData.get(i));

                trees.get(0).waveletData.set(pm1.get(0));
                trees.get(5).waveletData.set(pm1.get(5));
                trees.get(3).waveletData.set(pm1.get(3));
                trees.get(6).waveletData.set(pm1.get(6));

                ArrayList<ArrayList<double[][]>> pm2 = pm4(
                        trees.get(8).waveletData.get(i), 
                        trees.get(2).waveletData.get(i),
                        trees.get(4).waveletData.get(i),
                        trees.get(1).waveletData.get(i));

                trees.get(8).waveletData.set(pm1.get(8));
                trees.get(2).waveletData.set(pm1.get(2));
                trees.get(4).waveletData.set(pm1.get(4));
                trees.get(1).waveletData.set(pm1.get(1));

            }
        }
    }

    private void addSubtractInv() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images

                ArrayList<ArrayList<double[][]>> pm1 = pm4inv(
                        trees.get(0).waveletData.get(i), 
                        trees.get(5).waveletData.get(i),
                        trees.get(3).waveletData.get(i),
                        trees.get(6).waveletData.get(i));

                trees.get(0).waveletData.set(pm1.get(0));
                trees.get(5).waveletData.set(pm1.get(5));
                trees.get(3).waveletData.set(pm1.get(3));
                trees.get(6).waveletData.set(pm1.get(6));

                ArrayList<ArrayList<double[][]>> pm2 = pm4inv(
                        trees.get(8).waveletData.get(i), 
                        trees.get(2).waveletData.get(i),
                        trees.get(4).waveletData.get(i),
                        trees.get(1).waveletData.get(i));

                trees.get(8).waveletData.set(pm1.get(8));
                trees.get(2).waveletData.set(pm1.get(2));
                trees.get(4).waveletData.set(pm1.get(4));
                trees.get(1).waveletData.set(pm1.get(1));

            }
        }
    }

	private static ArrayList<ArrayList<double[][][]>> pm4 (ArrayList<double[][][]> a, ArrayList<double[][][]> b,
			ArrayList<double[][][]> c, ArrayList<double[][][]> d) {
		ArrayList<double[][][]> p = new ArrayList<double[][][]>();
		ArrayList<double[][][]> q = new ArrayList<double[][][]>();
		ArrayList<double[][][]> r = new ArrayList<double[][][]>();
		ArrayList<double[][][]> s = new ArrayList<double[][][]>();
		// no operation on lo pass
		p.add(JVCLUtils.deepCopy(a.get(0)));
		q.add(JVCLUtils.deepCopy(b.get(0)));
		r.add(JVCLUtils.deepCopy(c.get(0)));
		s.add(JVCLUtils.deepCopy(d.get(0)));
		// pm on high pass
		for (int i = 1; i < a.size(); i++) {
			p.add(
				ArrayMath.divide(
					ArrayMath.subtract(
						ArrayMath.subtract(
							ArrayMath.subtract(
								JVCLUtils.deepCopy(a.get(i)),
								JVCLUtils.deepCopy(b.get(i))
							), JVCLUtils.deepCopy(c.get(i))
						), JVCLUtils.deepCopy(d.get(i))
					), 2
				)
			);
			q.add(
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.add(
							ArrayMath.subtract(
								JVCLUtils.deepCopy(a.get(i)),
								JVCLUtils.deepCopy(b.get(i))
							), JVCLUtils.deepCopy(c.get(i))
						), JVCLUtils.deepCopy(d.get(i))
					), 2
				)
			);
			r.add(
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.subtract(
							ArrayMath.add(
								JVCLUtils.deepCopy(a.get(i)),
								JVCLUtils.deepCopy(b.get(i))
							), JVCLUtils.deepCopy(c.get(i))
						), JVCLUtils.deepCopy(d.get(i))
					), 2
				)
			);
			s.add(
				ArrayMath.divide(
					ArrayMath.subtract(
						ArrayMath.add(
							ArrayMath.add(
								JVCLUtils.deepCopy(a.get(i)),
								JVCLUtils.deepCopy(b.get(i))
							), JVCLUtils.deepCopy(c.get(i))
						), JVCLUtils.deepCopy(d.get(i))
					), 2
				)
			);
		}
		ArrayList<ArrayList<double[][][]>> pm = new ArrayList<ArrayList<double[][][]>>();
		pm.add(p);
		pm.add(q);
		pm.add(r);
		pm.add(s);
		return pm;
	}

	private static ArrayList<ArrayList<double[][][]>> pm4inv (ArrayList<double[][][]> a, ArrayList<double[][][]> b,
			ArrayList<double[][][]> c, ArrayList<double[][][]> d) {
		ArrayList<double[][][]> p = new ArrayList<double[][][]>();
		ArrayList<double[][][]> q = new ArrayList<double[][][]>();
		ArrayList<double[][][]> r = new ArrayList<double[][][]>();
		ArrayList<double[][][]> s = new ArrayList<double[][][]>();
		// no operation on lo pass
		p.add(JVCLUtils.deepCopy(a.get(0)));
		q.add(JVCLUtils.deepCopy(b.get(0)));
		r.add(JVCLUtils.deepCopy(c.get(0)));
		s.add(JVCLUtils.deepCopy(d.get(0)));
		// pm on high pass
		for (int i = 1; i < a.size(); i++) {
			p.add(
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.add(
							ArrayMath.add(
								JVCLUtils.deepCopy(a.get(i)),
								JVCLUtils.deepCopy(b.get(i))
							), JVCLUtils.deepCopy(c.get(i))
						), JVCLUtils.deepCopy(d.get(i))
					), 2
				)
			);
			q.add(
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.add(
							ArrayMath.subtract(
								ArrayMath.multiply(JVCLUtils.deepCopy(a.get(i)), -1),
								JVCLUtils.deepCopy(b.get(i))
							), JVCLUtils.deepCopy(c.get(i))
						), JVCLUtils.deepCopy(d.get(i))
					), 2
				)
			);
			r.add(
				ArrayMath.divide(
					ArrayMath.add(
						ArrayMath.subtract(
							ArrayMath.add(
								ArrayMath.multiply(JVCLUtils.deepCopy(a.get(i)), -1),
								JVCLUtils.deepCopy(b.get(i))
							), JVCLUtils.deepCopy(c.get(i))
						), JVCLUtils.deepCopy(d.get(i))
					), 2
				)
			);
			s.add(
				ArrayMath.divide(
					ArrayMath.subtract(
						ArrayMath.add(
							ArrayMath.add(
								ArrayMath.multiply(JVCLUtils.deepCopy(a.get(i)), -1),
								JVCLUtils.deepCopy(b.get(i))
							), JVCLUtils.deepCopy(c.get(i))
						), JVCLUtils.deepCopy(d.get(i))
					), 2
				)
			);
		}
		ArrayList<ArrayList<double[][][]>> pm = new ArrayList<ArrayList<double[][][]>>();
		pm.add(p);
		pm.add(q);
		pm.add(r);
		pm.add(s);
		return pm;
	}
}
