package com.ericbarnhill.jmra.dualTree;

import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;

abstract class DualTree<N, B, V> {

    public ArrayList<MRA<N, B, V>> trees;
    public ArrayList<DTFilterBank> banks;

    N origData;
    B maskData;
    DTFilterBank fb;
    int decompLvls;
    public ConvolverFactory.ConvolutionType convType;
    public boolean undecimated;

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

    public abstract void setTrees();

    public abstract void dwt();

    public abstract void idwt();

    public abstract void addSubtract(boolean fwd);

    public abstract N getFilteredData();

    public N getOrigData() {
        return origData;
    }

}
