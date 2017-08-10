package com.ericbarnhill.jmra.dualTree;

import java.util.ArrayList; 
import com.ericbarnhill.arrayMath.*;
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;

public class DualTree3DCplxSerial extends DualTree3DCplx {

    public ArrayList<MRA3DDTSerial> trees;

    public DualTree3DCplxSerial(double[][][] origData, boolean[][][] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType, boolean undecimated) {
        super(origData, maskData, fb, decompLvls, convType, undecimated);
        this.stride = 8;
        this.trees = new ArrayList<MRA3DDTSerial>();
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
    public void dwt() {
        setTrees();
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).dwt();
        }
        addSubtract(true);
    }

    @Override
    public void idwt() {
        addSubtract(false);
        for (int i = 0; i < trees.size(); i++) {
            trees.get(i).idwt();
        }
    }

    @Override
    public double[][][] getFilteredData() {
        double[][][] bankSum = new double[origData.length][origData[0].length][origData[0][0].length];
        for (MRA3DDTSerial tree : trees) {
            bankSum = ArrayMath.add(bankSum, tree.getFilteredData());
        }
        bankSum = ArrayMath.divide(bankSum, Math.sqrt(8));
        return bankSum;
    }

    @Override
    void addSubtractFwd() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images

                ArrayList<double[][][]> pm1 = pm4(
                        Serializer.loadData(trees.get(0).waveletTempFiles.get(i)), 
                        Serializer.loadData(trees.get(5).waveletTempFiles.get(i)),
                        Serializer.loadData(trees.get(3).waveletTempFiles.get(i)),
                        Serializer.loadData(trees.get(6).waveletTempFiles.get(i)));

                Serializer.writeData(pm1.get(0), trees.get(0).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(1), trees.get(5).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(2), trees.get(3).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(3), trees.get(6).waveletTempFiles.get(i));

                ArrayList<double[][][]> pm2 = pm4(
                        Serializer.loadData(trees.get(7).waveletTempFiles.get(i)), 
                        Serializer.loadData(trees.get(2).waveletTempFiles.get(i)),
                        Serializer.loadData(trees.get(4).waveletTempFiles.get(i)),
                        Serializer.loadData(trees.get(1).waveletTempFiles.get(i)));

                Serializer.writeData(pm1.get(0), trees.get(7).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(1), trees.get(2).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(2), trees.get(4).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(3), trees.get(1).waveletTempFiles.get(i));

            }
        }
    }

    @Override
    void addSubtractInv() {
        for (int i = 0; i < stride*decompLvls; i++) {
            if (i % stride != 0) { // skip low pass images

                ArrayList<double[][][]> pm1 = pm4inv(
                        Serializer.loadData(trees.get(0).waveletTempFiles.get(i)), 
                        Serializer.loadData(trees.get(5).waveletTempFiles.get(i)),
                        Serializer.loadData(trees.get(3).waveletTempFiles.get(i)),
                        Serializer.loadData(trees.get(6).waveletTempFiles.get(i)));

                Serializer.writeData(pm1.get(0), trees.get(0).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(1), trees.get(5).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(2), trees.get(3).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(3), trees.get(6).waveletTempFiles.get(i));

                ArrayList<double[][][]> pm2 = pm4inv(
                        Serializer.loadData(trees.get(7).waveletTempFiles.get(i)), 
                        Serializer.loadData(trees.get(2).waveletTempFiles.get(i)),
                        Serializer.loadData(trees.get(4).waveletTempFiles.get(i)),
                        Serializer.loadData(trees.get(1).waveletTempFiles.get(i)));

                Serializer.writeData(pm1.get(0), trees.get(7).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(1), trees.get(2).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(2), trees.get(4).waveletTempFiles.get(i));
                Serializer.writeData(pm1.get(3), trees.get(1).waveletTempFiles.get(i));

            }
        }
    }

    @Override
    public void accept(Threshold threshold) {
        threshold.visitSerial(this);
    }
}
