package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

class MRA1DDTU extends MRA1DDT {

    MRA1DU mra1du;

    public MRA1DDTU(double[] origData, boolean[] maskData, ArrayList<ArrayList<ArrayList<double[]>>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        mra1du = new MRA1DU();
    }

    public MRA1DDTU(double[] origData, ArrayList<ArrayList<ArrayList<double[]>>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length), filterBank, decompLvls, convType);
    }

    @Override
    double[] AFB(double[] x, double[] filter, int decompLvl) {
        return mra1du.AFB(x, filter, decompLvl);
    }

    @Override
    double[] SFB(double[] lo, double[] hi, double[] sfl, double[] sfh, int decompLvl) {
        return mra1du.SFB(lo, hi, sfl, sfh, decompLvl);
    }

}
