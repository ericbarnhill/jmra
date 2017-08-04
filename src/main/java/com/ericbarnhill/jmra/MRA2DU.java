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

class MRA2DU extends MRA2D {

    MRA1DU mra1du;

    public MRA2DU() {
        super();
    }

    public MRA2DU(double[][] origData, boolean[][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        mra1du = new MRA1DU(convType);
    }

    public MRA2DU(double[][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length), filterBank, decompLvls, convType);
    }

    @Override
    double[][] AFB(double[][] data, double[] filter, int decompositionLevel) {
        final int J = decompLvls - decompositionLevel;
        final int fi = data.length;
        final int fj = data[0].length;
        final int N = filter.length;
        double[][] filtData = new double[fi][];
        for (int i = 0; i < fi; i++) {
            filtData[i] = mra1du.AFB(data[i], filter, J);
        }
        return filtData;
    }

    @Override
    double[][] SFB(double[][] lo, double[][] hi, double[] sfl, double[] sfh, int decompositionLevel) {
        final int J = decompLvls - decompositionLevel;
        final int N = sfl.length + sfh.length;
        final int fi = lo.length;
        final int fj = lo[0].length;
        double[][] y = new double[fi][];
        for (int i = 0; i < fi; i++) {
            y[i] = mra1du.SFB(lo[i], hi[i], sfl, sfh, J);
        }
        return y;
    }

    public void accept(Threshold threshold) {
        threshold.visit(this);
    }

    // for debugging and testing
    public void data2File(double[][] data, String path) {
        int w = data.length;
        int h = data[0].length;
        FloatProcessor fp = new FloatProcessor(w,h);
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                fp.putPixelValue(i,j,data[i][j]);
            }
        }
        ImagePlus ip = new ImagePlus("", fp);
        FileSaver fs = new FileSaver(ip);
        fs.saveAsTiff(path);
    }

}
