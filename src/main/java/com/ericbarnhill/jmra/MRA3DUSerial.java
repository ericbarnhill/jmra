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
import org.apache.commons.math4.stat.descriptive.moment.Mean;

public class MRA3DUSerial extends MRA3DSerial {

     MRA1DU mra1du;

     public MRA3DUSerial(ConvolverFactory.ConvolutionType convType) {
         super();
         mra1du = new MRA1DU(convType);
     }

    public MRA3DUSerial(double[][][] origData, boolean[][][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        mra1du = new MRA1DU(convType);
        paddedData = origData;
    }

    public MRA3DUSerial(double[][][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        this(origData, ArrayMath.fillWithTrue(origData.length,origData[0].length, origData[0][0].length), filterBank, decompLvls, convType);
    }

    @Override
    public double[][][] AFB(double[][][] data, double[] filter, int decompositionLevel) {
        final int fi = data.length;
        final int fj = data[0].length;
        double[][][] filtData = new double[fi][fj][];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                filtData[i][j] = mra1du.AFB(data[i][j], filter, decompositionLevel);
                //filtData[i][j] = mra1d.AFB(data[i][j], filter, decompositionLevel);
            }
        }
        return filtData;
    }
    
    @Override
    public double[][][] SFB(double[][][] lo, double[][][] hi, double[] sfl, double[] sfh, int decompositionLevel) {
        final int fi = lo.length;
        final int fj = lo[0].length;
        double[][][] y = new double[fi][fj][];
        for (int i = 0; i < fi; i++) {
            for (int j = 0; j < fj; j++) { 
                y[i][j] = mra1du.SFB(lo[i][j], hi[i][j], sfl, sfh, decompositionLevel);
                //y[i][j] = mra1d.SFB(lo[i][j], hi[i][j], sfl, sfh, decompositionLevel);
            }
        }
        return y;
    }

    @Override
    public double[][][] getFilteredData() {
        return Serializer.loadData(waveletTempFiles.get(0));
    }

}
