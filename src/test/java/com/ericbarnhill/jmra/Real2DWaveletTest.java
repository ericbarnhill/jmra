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
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Min;
import org.apache.commons.math4.stat.descriptive.rank.Median;

public class Real2DWaveletTest {


    public static double[][] image2Array(String path) {
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

    public static void array2Image(double[][] array, String path) {
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
    public void testWaveletRecon() {
        String file = "/home/ericbarnhill/Documents/code/lena.tif";
        String file2 = "/home/ericbarnhill/Documents/code/lenafilt.tif";
        double[][] array = image2Array(file);
        double[] h0 = Wavelets.afLoNoTree;
        double[] h1 = Wavelets.afHiNoTree;
        double[] g0 = Wavelets.sfLoNoTree;
        double[] g1 = Wavelets.sfHiNoTree;
        ArrayList<double[]> analysisFilters = new ArrayList<double[]>();
        analysisFilters.add(h0);
        analysisFilters.add(h1);
        ArrayList<double[]> synthesisFilters = new ArrayList<double[]>();
        synthesisFilters.add(g0);
        synthesisFilters.add(g1);
        ArrayList<ArrayList<double[]>> filterBank = new ArrayList<ArrayList<double[]>>();
        filterBank.add(analysisFilters);
        filterBank.add(synthesisFilters);

        double[][] image = image2Array(file);
        double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        //image = ArrayMath.add(ArrayMath.multiply(noise, 10), image);
        array2Image(image, "/home/ericbarnhill/Documents/code/lenanoise.tif");
        MRA2D mra = new MRA2D(image, filterBank, 3);
        mra.dwt();
        List<double[][]> decomp = mra.getDecomposition();
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_before.tif";
            array2Image(decomp.get(n), path);
        }
        mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);

        mra.idwt();
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_after.tif";
            array2Image(decomp.get(n), path);
        }
        image = mra.getFilteredData();
        array2Image(image, file2);
    }
}

        

