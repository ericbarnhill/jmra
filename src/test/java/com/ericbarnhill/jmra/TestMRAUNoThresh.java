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

    @Test
    public void MRA2DUTest() {
        // PREP IMAGE
        double[][] image = FilePaths.image2Array(FilePaths.image2D);
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
        String resultFile = FilePaths.root + "lena_idwt_u.tif";
        FilePaths.array2Image(filtImage, resultFile);
    }

    @Ignore
    public void MRA3DUTest() {
        System.out.println("MRA 3DU Test");
        // PREP IMAGE
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(FilePaths.nifti3D);
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
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        mra.idwt();
        double[][][] filteredData = mra.getFilteredData();
        FilePaths.data2File(filteredData, FilePaths.root+"filtdata_3d_u.tif");
    }
}

        

