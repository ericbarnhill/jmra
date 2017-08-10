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
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.dualTree.*;
import com.ericbarnhill.niftijio.*;
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Min;
import org.apache.commons.math4.stat.descriptive.rank.Median;

public class TestDualTreeThresh{


    @Test
    public void MRA2DTTest() {
        System.out.println("MRA2DDT Test");
        double[][] image = FilePaths.image2Array(FilePaths.image2D);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA2DDT mra = new MRA2DDT(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        List<double[][]> decomp = mra.getDecomposition();
        
        mra.idwt();
        image = mra.getFilteredData();
        String resultFile = FilePaths.root + "lena_mra2ddt_idwt.tif";
        FilePaths.array2Image(image, resultFile);
    }

    @Test
    public void MRA2DTUTest() {
        System.out.println("MRA2DDTU Test");
        double[][] image = FilePaths.image2Array(FilePaths.image2D);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA2DDTU mra = new MRA2DDTU(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        List<double[][]> decomp = mra.getDecomposition();
        mra.idwt();
        image = mra.getFilteredData();
        String resultFile = FilePaths.root + "lena_dtu_idwt.tif";
        FilePaths.array2Image(image, resultFile);
    }

    @Test
    public void DualTree2DTest() {
        System.out.println("Thresh DT2D");
        double[][] image = FilePaths.image2Array(FilePaths.image2D);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        double[][] noise = ArrayMath.multiply(ArrayMath.fillWithRandom(image.length, image[0].length),100);
        double[][] noisyImg = ArrayMath.add(image, noise);
        DualTree2DCplx dt = new DualTree2DCplx(noisyImg, fb, 3, ConvolverFactory.ConvolutionType.FDCPU, true);
        dt.dwt();
        Threshold t = new Threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        dt.accept(t);
        dt.idwt();
        image = dt.getFilteredData();
        String resultFile = FilePaths.root + "lena_mra2ddtu_idwt_thresh.tif";
        FilePaths.array2Image(image, resultFile);
    }

    @Ignore
    public void MRA3DDTTest() {
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(FilePaths.nifti3D);
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
        mra.data2File(filteredData, FilePaths.root+"filtdata_3d_u.tif");
    }

    @Test
    public void DualTree3DTest() {
        long t1 = System.currentTimeMillis();
        System.out.println("Dualtree 3D non serial with thresh");
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(FilePaths.nifti3D);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][][][] niftiArray = nv.data.toArray();
        double[][][] image = ArrayMath.convert4dto3d(niftiArray);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        DualTree3DCplx dt = new DualTree3DCplx(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU, true);
        dt.dwt();
        Visualizer.dumpDecomposition(dt);
        Threshold threshold = new Threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        threshold.visit(dt);
        dt.idwt();
        image = dt.getFilteredData();
        String resultFile = FilePaths.root + "dualtree3d_nonserial_thresh.tif";
        new MRA3D().data2File(image, resultFile);
        long t2 = System.currentTimeMillis();
        System.out.println("Time unserialized " + ((t2-t1)/1000) + " sec");

    }


}

        

