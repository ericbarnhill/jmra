package com.ericbarnhill.jmra.dualTree;

import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;

public class DualTree<N, B, V> {

    public ArrayList<MRA<N, B, V>> trees;
    public ArrayList<DTFilterBank> banks;

    N origData;
    B maskData;
    DTFilterBank fb;
    int decompLvls;
    ConvolverFactory.ConvolutionType convType;

    public DualTree(N origData, B maskData, DTFilterBank fb, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this.origData = origData;
        this.maskData = maskData;
        this.fb = fb;
        this.decompLvls = decompLvls;
        this.convType = convType;
        setTrees();
    }

    public abstract void setTrees();

    public abstract void dwt();

    public abstract void idwt();

    public abstract void addSubtract();

    public N getFilteredData();

}
