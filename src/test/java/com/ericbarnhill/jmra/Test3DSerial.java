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

public class Test3DSerial{
    @Test
    public void MRA3DSerialTest() {
        System.out.println("MRA 3D Serial Test");
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(FileOps.nifti3D);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][][][] niftiArray = nv.data.toArray();
        double[][][] image = ArrayMath.convert4dto3d(niftiArray);
        // PREP MRA
        FilterBank fb = Wavelets.getFarras();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA3DSerial mra = new MRA3DSerial(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        ArrayList<double[][][]> decomp = mra.getDecomposition();
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        mra.idwt();
        double[][][] filteredData = mra.getFilteredData();
        FileOps.data2Nifti(filteredData, FileOps.imgDir+"filtdata_3d_serial.tif");
    }

    @Test
    public void MRA3DUSerialTest() {
        System.out.println("MRA 3DU Serial Test");
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(FileOps.nifti3D);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][][][] niftiArray = nv.data.toArray();
        double[][][] image = ArrayMath.convert4dto3d(niftiArray);
        // PREP MRA
        FilterBank fb = Wavelets.getFarras();
        //double[][] noise = ArrayMath.fillWithRandom(image.length, image[0].length);
        MRA3DUSerial mra = new MRA3DUSerial(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        ArrayList<double[][][]> decomp = mra.getDecomposition();
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        mra.idwt();
        double[][][] filteredData = mra.getFilteredData();
        FileOps.data2Nifti(filteredData, FileOps.imgDir+"filtdata_3d_serial_u.tif");
    }
}
        

