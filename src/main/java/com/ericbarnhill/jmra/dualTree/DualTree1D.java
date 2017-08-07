package com.ericbarnhill.jmra.dualTree;

import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.arrayMath.*;

public class DualTree1D extends DualTree<double[], boolean[], double[]> {

    public DualTree1D(double[] origData, boolean[] maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, fb, decompLvls, convType);
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
