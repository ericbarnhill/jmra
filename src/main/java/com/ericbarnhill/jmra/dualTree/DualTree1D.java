package com.ericbarnhill.jmra.dualTree;

import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;

public class DualTree1D extends DualTree<double[], boolean[], double[]> {

    public DualTree1D(double[] origData, boolean[] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, fb, decompLvls, convType);
    }

    public void setTrees() {
        DTFilterBank tree1Bank = new DTFilterBank(fb.faf1, fb.faf1, fb.fsf1, fb.fsf1, fb.af1, fb.af1, fb.sf1, fb.sf1);
        DTFilterBank tree2Bank = new DTFilterBank(fb.faf2, fb.faf2, fb.fsf2, fb.fsf2, fb.af2, fb.af2, fb.sf2, fb.sf2);
        trees.add(new MRA1DDT(origData, maskData, tree1Bank, decompLvls, convType));
        trees.add(new MRA1DDT(origData, maskData, tree2Bank, decompLvls, convType));
    }

    public double[] getFilteredData() {
        return ArrayMath.divide(ArrayMath.add(tree1Bank.getFilteredData(), tree2Bank.getFilteredData()), Math.sqrt(2));
    }
    
    public void addSubtract() {}

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
