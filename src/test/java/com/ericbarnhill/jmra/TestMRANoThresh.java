package com.ericbarnhill.jmra;

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

public class TestMRANoThresh{

    @Test
    public void MRA2DTest() {
        System.out.println("MRA 2D Test");
        // PREP IMAGE
        double[][] image = FilePaths.image2Array(FilePaths.image2D);
        // PREP MRA
        FilterBank fb = Wavelets.getFarras();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA2D mra = new MRA2D(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        List<double[][]> decomp = mra.getDecomposition();
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);

        mra.idwt();
        image = mra.getFilteredData();
        String resultFile = FilePaths.root + "lena_idwt.tif";
        FilePaths.array2Image(image, resultFile);
    }

    @Test
    public void MRA3DTest() {
        System.out.println("MRA 3D Test");
        // PREP IMAGE
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(FilePaths.nifti3D);
        } catch (Exception e) {}
        double[][][][] niftiArray = nv.data.toArray();
        double[][][] image = ArrayMath.convert4dto3d(niftiArray);
        // PREP MRA
        FilterBank fb = Wavelets.getFarras();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA3D mra = new MRA3D(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        ArrayList<double[][][]> decomp = mra.getDecomposition();
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        mra.idwt();
        double[][][] filteredData = mra.getFilteredData();
        mra.data2File(filteredData, FilePaths.root+"filtdata_3d.tif");
    }
}

        

