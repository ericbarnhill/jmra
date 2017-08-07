package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

public class MRA3DU extends MRA3D {

     MRA1DU mra1du;

     public MRA3DU() {
         super();
     }

    public MRA3DU(double[][][] origData, boolean[][][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        mra1du = new MRA1DU(convType);
    }

    public MRA3DU(double[][][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length, origData[0][0].length), filterBank, decompLvls, convType);
    }

    @Override
    public double[][][] AFB(double[][][] data, double[] filter, int decompositionLevel) {
        final int fi = data.length;
        final int fj = data[0].length;
        final int fk = data[0][0].length;
        double[][][] filtData = new double[fi][fj][fk];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                filtData[i][j] = mra1du.AFB(data[i][j], filter, decompositionLevel);
            }
        }
        return filtData;
    }
    
    @Override
    public double[][][] SFB(double[][][] lo, double[][][] hi, double[] sfl, double[] sfh, int decompositionLevel) {
        final int fi = lo.length;
        final int fj = lo[0].length;
        double[][][] y = new double[fi][fj][];
        //DEBUGGING
        ArrayMath.displaySize(lo);
        ArrayMath.displaySize(hi);
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                y[i][j] = mra1du.SFB(lo[i][j], hi[i][j], sfl, sfh, decompositionLevel);
            }
        }
        return y;
    }

    public void accept(Threshold threshold) {
        threshold.visit(this);
    }

    // for debugging and testing
    public void data2File(double[][][] data, String path) {
        int w = data.length;
        int h = data[0].length;
        int d = data[0][0].length;
        ImageStack is = new ImageStack(w,h);
        for (int k = 0; k < d; k++) {
            FloatProcessor fp = new FloatProcessor(w,h);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                            fp.putPixelValue(i,j,data[i][j][k]);
                }
            }
            is.addSlice(fp);
        }
        ImagePlus ip = new ImagePlus("", is);
        FileSaver fs = new FileSaver(ip);
        fs.saveAsTiff(path);
    }



        


}
