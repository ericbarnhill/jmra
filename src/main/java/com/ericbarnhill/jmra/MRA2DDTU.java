package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

class MRA2DDTU extends MRA2DDT {

    MRA2DU mra2du;

    public MRA2DDTU(double[][] origData, boolean[][] maskData, DTFilterSet dtfs, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, dtfs, decompLvls, convType);
        mra2du = new MRA2DU();
    }

    public MRA2DDTU(double[][] origData, DTFilterSet dtfs, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length), dtfs, decompLvls, convType);
    }

    @Override
    double[][] AFB(double[][] x, double[] filter, int decompLvl) {
        return mra2du.AFB(x, filter, decompLvl);
    }

    @Override
    double[][] SFB(double[][] lo, double[][] hi, double[] sfl, double[] sfh, int decompLvl) {
        return mra2du.SFB(lo, hi, sfl, sfh, decompLvl);
    }

}
