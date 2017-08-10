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

    @Test
    public void MRA2DTTest() {
        System.out.println("MRA2DDT Test");
        // PREP IMAGE
        double[][] image = FilePaths.image2Array(FilePaths.image2D);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA2DDT mra = new MRA2DDT(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        List<double[][]> decomp = mra.getDecomposition();
        
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);

        mra.idwt();
        image = mra.getFilteredData();
        String resultFile = FilePaths.root + "lena_dt_idwt.tif";
        FilePaths.array2Image(image, resultFile);
    }

    @Test
    public void DualTree2DTest() {
        System.out.println("No Thresh DualTree 2D Test");
        // PREP IMAGE
        double[][] image = FilePaths.image2Array(FilePaths.image2D);
        // PREP MRA
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        DualTree2DCplx dt = new DualTree2DCplx(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU, true);
        dt.dwt();
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        dt.idwt();
        image = dt.getFilteredData();
        String resultFile = FilePaths.root + "lena_dualtree_idwt.tif";
        FilePaths.array2Image(image, resultFile);
    }

    @Test
    public void MRA3DDTTest() {
        System.out.println("No Thresh MRE3DDT Test");
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
        mra.data2File(filteredData, FilePaths.root+"filtdata_3d_dt.tif");
    }
    @Test
    public void DualTree3DSerialTest() {
        long t1 = System.currentTimeMillis();
        System.out.println("No Thresh DT3D Serial");
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
        DualTree3DCplxSerial dt = new DualTree3DCplxSerial(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU, false);
        dt.dwt();
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        dt.idwt();
        image = dt.getFilteredData();
        String resultFile = FilePaths.root + "dualtree3d_serial_filt.tif";
        new MRA3D().data2File(image, resultFile);
        long t2 = System.currentTimeMillis();
        System.out.println("Time serialized " + ((t2-t1)/1000) + " sec");
    }

}

        

