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
import com.ericbarnhill.niftijio.*;
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Min;
import org.apache.commons.math4.stat.descriptive.rank.Median;

public class TestMRAUNoThresh{


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
    public void MRA2DUTest() {
        // PREP IMAGE
        String root = "/home/ericbarnhill/Documents/code/jmra/test-images/"; 
        String testFile = root + "lena.tif";
        double[][] image = image2Array(testFile);
        // PREP MRA
        FilterBank fb = Wavelets.getFarras();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA2DU mra = new MRA2DU(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        List<double[][]> decomp = mra.getDecomposition();
        /*
        System.out.println("Displaying decomp sizes: "); 
        for (int n = 0; n < decomp.size(); n++) {
            ArrayMath.displaySize(decomp.get(n));
            String path = root + Integer.toString(n)+ "_before.tif";
            array2Image(decomp.get(n), path);
        }
        */
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);

        mra.idwt();
        /*
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_after.tif";
            array2Image(decomp.get(n), path);
        }
        */
        double[][] filtImage = mra.getFilteredData();
        String resultFile = root + "lena_idwt_u.tif";
        array2Image(filtImage, resultFile);
    }

    @Test
    public void MRA3DUTest() {
        System.out.println("MRA 3DU Test");
        // PREP IMAGE
        String filepath = "/home/ericbarnhill/Documents/MATLAB/ericbarnhill/projects/2017-07-06-florian-new-protocol/scratch/fieldmaps/3.nii";
        String root = "/home/ericbarnhill/Documents/code/jmra/test-images/"; 
        String outputpath = "/home/ericbarnhill/Documents/code/jmra/test-images/test_3d_u.nii";
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(filepath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][][][] niftiArray = nv.data.toArray();
        double[][][] image = ArrayMath.convert4dto3d(niftiArray);
        // PREP MRA
        FilterBank fb = Wavelets.getFarras();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA3DU mra = new MRA3DU(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        ArrayList<double[][][]> decomp = mra.getDecomposition();
        System.out.println("3d decomp size "+decomp.size());
        for (int n = 0; n < decomp.size(); n++) {
            String path = root + Integer.toString(n)+ "_before_3d.tif";
            mra.data2File(decomp.get(n), path);
        }
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        mra.idwt();
        double[][][] filteredData = mra.getFilteredData();
        mra.data2File(filteredData, root+"filtdata_3d_u.tif");
        /*
        nv.data = new FourDimensionalArray(ArrayMath.convertTo4d(filteredData));
        nv.header.dim[4] = 1;
        try {
            nv.write(outputpath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        */
    }
}

        

