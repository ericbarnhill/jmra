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
import com.ericbarnhill.jmra.dualTree.*;
import com.ericbarnhill.niftijio.*;
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Min;
import org.apache.commons.math4.stat.descriptive.rank.Median;

public class TestDualTreeNoThresh{


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
    public void MRA2DTTest() {
        System.out.println("MRA2DDT Test");
        // PREP IMAGE
        String root = "/home/ericbarnhill/Documents/code/jmra/test-images/"; 
        String testFile = root + "lena.tif";
        double[][] image = image2Array(testFile);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA2DDT mra = new MRA2DDT(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        List<double[][]> decomp = mra.getDecomposition();
        
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);

        mra.idwt();
        /*
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_after.tif";
            array2Image(decomp.get(n), path);
        }
        */
        image = mra.getFilteredData();
        String resultFile = root + "lena_dt_idwt.tif";
        array2Image(image, resultFile);
    }

    @Test
    public void DualTree2DTest() {
        System.out.println("DualTree 2D Test");
        // PREP IMAGE
        String root = "/home/ericbarnhill/Documents/code/jmra/test-images/"; 
        String testFile = root + "lena.tif";
        double[][] image = image2Array(testFile);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        DualTree2DCplx dt = new DualTree2DCplx(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        dt.dwt();
        /*
        System.out.println("Displaying trees and decomp sizes: "); 
        for (int i = 0; i < dt.trees.size(); i++) {
            MRA<double[][], boolean[][], double[]> tree = dt.trees.get(i);
            ArrayList<double[][]> decomp = tree.getDecomposition();
            for (int j = 0; j < decomp.size(); j++) {
                ArrayMath.displaySize(decomp.get(j));
                String path = root + Integer.toString(j)+ "_before_dualtree.tif";
                array2Image(decomp.get(j), path);
            }
        }
        */
       
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);

        dt.idwt();
        /*
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_after.tif";
            array2Image(decomp.get(n), path);
        }
        */
        image = dt.getFilteredData();
        String resultFile = root + "lena_dualtree_idwt.tif";
        array2Image(image, resultFile);
    }

    @Test
    public void MRA3DDTTest() {
        // PREP IMAGE
        String root = "/home/ericbarnhill/Documents/code/jmra/test-images/"; 
        String filepath = "/home/ericbarnhill/Documents/MATLAB/ericbarnhill/projects/2017-07-06-florian-new-protocol/scratch/fieldmaps/3.nii";
        String outputpath = "/home/ericbarnhill/Documents/code/jmra/test-images/test_3DDT.nii";
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(filepath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][][][] niftiArray = nv.data.toArray();
        double[][][] image = ArrayMath.convert4dto3d(niftiArray);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        DTFilterBank fb3d = new DTFilterBank(fb.faf.get(0), fb.faf.get(0), fb.faf.get(0), fb.fsf.get(0), fb.fsf.get(0), fb.fsf.get(0), fb.af.get(0), fb.af.get(0), fb.af.get(0), fb.sf.get(0), fb.sf.get(0), fb.sf.get(0));
        MRA3DDT mra = new MRA3DDT(image, fb3d, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        ArrayList<double[][][]> decomp = mra.getDecomposition();
        mra.idwt();
        double[][][] filteredData = mra.getFilteredData();
        mra.data2File(filteredData, root+"filtdata_3d_dt.tif");
    }

    @Test
    public void DualTree3DTest() {
        String root = "/home/ericbarnhill/Documents/code/jmra/test-images/"; 
        String filepath = "/home/ericbarnhill/Documents/MATLAB/ericbarnhill/projects/2017-07-06-florian-new-protocol/scratch/fieldmaps/3.nii";
        String outputpath = "/home/ericbarnhill/Documents/code/jmra/test-images/test_3DDT.nii";
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(filepath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][][][] niftiArray = nv.data.toArray();
        double[][][] image = ArrayMath.convert4dto3d(niftiArray);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        DualTree3DCplx dt = new DualTree3DCplx(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        dt.dwt();
        /*
        System.out.println("Displaying trees and decomp sizes: "); 
        for (int i = 0; i < dt.trees.size(); i++) {
            MRA<double[][], boolean[][], double[]> tree = dt.trees.get(i);
            ArrayList<double[][]> decomp = tree.getDecomposition();
            for (int j = 0; j < decomp.size(); j++) {
                ArrayMath.displaySize(decomp.get(j));
                String path = root + Integer.toString(j)+ "_before_dualtree.tif";
                array2Image(decomp.get(j), path);
            }
        }
        */
       
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);

        dt.idwt();
        /*
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_after.tif";
            array2Image(decomp.get(n), path);
        }
        */
        image = dt.getFilteredData();
        String resultFile = root + "dualtree3d_filt.tif";
        new MRA3D().data2File(image, resultFile);
    }

}

        

