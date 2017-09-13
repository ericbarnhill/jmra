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
import com.ericbarnhill.jmra.filters.*;

public class MRA2DU extends MRA2D {

    MRA1DU mra1du;

    public MRA2DU(ConvolverFactory.ConvolutionType convType) {
        super(convType);
        mra1du = new MRA1DU(convType);
    }

    public MRA2DU(double[][] origData, boolean[][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        mra1du = new MRA1DU(convType);
    }

    public MRA2DU(double[][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length), filterBank, decompLvls, convType);
    }

    @Override
    public double[][] analysis(double[][] data, double[] filter, int decompositionLevel) {
        final int J = decompLvls - decompositionLevel;
        final int fi = data.length;
        final int fj = data[0].length;
        final int N = filter.length;
        double[][] filtData = new double[fi][];
        for (int i = 0; i < fi; i++) {
            //filtData[i] = mra1du.analysis(data[i], filter, J);
            filtData[i] = mra1d.analysis(data[i], filter, J);
        }
        return filtData;
    }

    @Override
    public double[][] synthesis(double[][] lo, double[][] hi, double[] sfl, double[] sfh, int decompositionLevel) {
        final int J = decompLvls - decompositionLevel;
        final int N = sfl.length + sfh.length;
        final int fi = lo.length;
        final int fj = lo[0].length;
        double[][] y = new double[fi][];
        for (int i = 0; i < fi; i++) {
            //y[i] = mra1du.synthesis(lo[i], hi[i], sfl, sfh, J);
            y[i] = mra1d.synthesis(lo[i], hi[i], sfl, sfh, J);
        }
        return y;
    }

    public void accept(Threshold threshold) {
        threshold.visit(this);
    }
}
