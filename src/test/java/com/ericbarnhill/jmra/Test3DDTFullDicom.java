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

public class Test3DDTFullDicom{


    @Ignore
    public void DualTree3DFullDicomTest() {
        long t1 = System.currentTimeMillis();
        System.out.println("DT3D Full Dicom");
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
        DualTree3DCplxSerial dt = new DualTree3DCplxSerial(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU, true);
        dt.dwt();
        Threshold threshold = new Threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        //threshold.visitSerial(dt);

        dt.idwt();
        /*
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_after.tif";
            array2Image(decomp.get(n), path);
        }
        */
        image = dt.getFilteredData();
        String resultFile = FilePaths.root + "dualtree3d_fulldicom.tif";
        new MRA3D().data2File(image, resultFile);
        long t2 = System.currentTimeMillis();
        System.out.println("Time for fullsize dicom processing:  " + ((t2-t1)/1000) + " sec");

    }
}

        

