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

class MRA2DU extends MRA2D {

    MRA1DU mra1du;

    public MRA2DU(double[][] originalData, boolean[][] maskData, ArrayList<ArrayList<double[]>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convolutionType) {
        super(originalData, maskData, filterBank, decompLvls, convolutionType);
        mra1du = new MRA1DU(convolutionType);
    }

    public MRA2DU(double[][] originalData, ArrayList<ArrayList<double[]>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convolutionType) {
        this(originalData, ArrayMath.fillWithTrue(originalData.length,originalData[0].length), filterBank, decompLvls, convolutionType);
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

    @Override
    public void threshold(Threshold.ThreshMeth threshMeth, Threshold.NoiseEstMeth noiseEstMeth) {
        // loop through each subband, pass method
        for (int i = 0; i < waveletData.size(); i++) {
            // avoid scaling datas
            if (i % stride != 0) {
                int level = (int)Math.floor(i / stride);
                double[] waveletVec = ArrayMath.vectorize(waveletData.get(i));
                boolean[] maskVec = ArrayMath.vectorize(maskData);
                waveletData.set(i, 
                    ArrayMath.devectorize(
                        Threshold.threshold(
                            waveletVec, maskVec, threshMeth, noiseEstMeth)
                        ,w)
                    );
            }
        }
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
