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
import com.ericbarnhill.niftijio.*;
import com.ericbarnhill.jvcl.*;
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Min;
import org.apache.commons.math4.stat.descriptive.rank.Median;

public class TestUndec3D {


    @Test
    public void testWaveletReconUndec() {
        String filepath = "/home/ericbarnhill/Documents/MATLAB/ericbarnhill/projects/2017-07-06-florian-new-protocol/scratch/fieldmaps/3.nii";
        String outputpath = "/home/ericbarnhill/Documents/code/test.nii";
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(filepath);
        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][][][] niftiArray = nv.data.toArray();
        int w = niftiArray.length;
        int h = niftiArray[0].length;
        int d = niftiArray[0][0].length;
        double[][][] image = new double[w][h][d];
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                for (int k = 0; k < d; k++) {
                    image[i][j][k] = niftiArray[i][j][k][0];
                }
            }
        }
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

        MRA3DU mra = new MRA3DU(image, filterBank, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        List<double[][][]> decomp = mra.getDecomposition();
        System.out.println("Displaying decomp sizes: "); 
        for (int n = 0; n < decomp.size(); n++) {
            ArrayMath.displaySize(decomp.get(n));
        }
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_3D_before.tif";
            mra.data2File(decomp.get(n), path);
        }
        //mra.threshold(Threshold.ThreshMeth.SOFT, Threshold.NoiseEstMeth.VISU_SHRINK);
        mra.idwt();
        for (int n = 0; n < decomp.size(); n++) {
            String path = "/home/ericbarnhill/Documents/code/" + Integer.toString(n)+ "_3D_after.tif";
            mra.data2File(decomp.get(n), path);
        }
        nv.data = new FourDimensionalArray(ArrayMath.convertTo4d(image));
        nv.header.dim[4] = 1;
        
        

        try {
            nv.write(outputpath);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}

        

