package com.ericbarnhill.jmra;

import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Min;
import org.apache.commons.math4.stat.descriptive.rank.Median;

public class MRANoThresh{


    static double[][] image2Array(String path) {
        try {
            ImagePlus ip = new Opener().openImage(path);
            ImageProcessor ipr = ip.getProcessor();
            int w = ip.getWidth();
            int h = ip.getHeight();
            double[][] array = new double[w][h];
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    array[i][j] = (double)ipr.getPixelValue(i,j);
                }
            }
            return array;
        } catch (Exception e) {
            return null;
        }
    }

    static void array2Image(double[][] array, String path) {
        int w = array.length;
        int h = array[0].length;
        FloatProcessor fp = new FloatProcessor(w,h);
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                fp.putPixelValue(i,j,array[i][j]);
            }
        }
        ImagePlus ip = new ImagePlus("", fp);
        FileSaver fs = new FileSaver(ip);
        fs.saveAsTiff(path);
    }

    @Test
    public void MRA2DTest() {
        // PREP IMAGE
        String root = "/home/ericbarnhill/Documents/code/jmra/test-images/"; 
        String testFile = root + "lena.tif";
        double[][] image = image2Array(testFile);
        // PREP MRA
        FilterBank fb = Wavelets.getFarras();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA2D mra = new MRA2D(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        List<double[][]> decomp = mra.getDecomposition();
        System.out.println("Displaying decomp sizes: "); 
        for (int n = 0; n < decomp.size(); n++) {
            ArrayMath.displaySize(decomp.get(n));
            String path = root + Integer.toString(n)+ "_before.tif";
            array2Image(decomp.get(n), path);
        }
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);

        mra.idwt();
        /*
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_after.tif";
            array2Image(decomp.get(n), path);
        }
        */
        image = mra.getFilteredData();
        String resultFile = root + "lena_idwt.tif";
        array2Image(image, resultFile);
    }
}

        

