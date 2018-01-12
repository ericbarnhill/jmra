/* 
 * Copyright (C) 2018 Eric Barnhill
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

package com.ericbarnhill.jmra;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.filters.*;
import com.ericbarnhill.jmra.dualTree.*;
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Min;
import org.apache.commons.math4.stat.descriptive.moment.Variance;
import org.apache.commons.math4.stat.descriptive.rank.Median;

/**
 * Tests whether hard, soft and nng thresholding methods
 * are in accord with the Matlab gold standard implementation.
 */
public class TestThresholding {

static final double EPS = 0.1;

    @Test
    public void SoftThresholdingTest2D() {
        final double[][] image = FileOps.image2Array(FileOps.imgSrcDir + "/boat_noisy.tif");
        FilterBank fb = Wavelets.getFarras();
        final double SIGMA = 14.3549;
        Threshold threshold = new Threshold(Threshold.ThreshMeth.SOFT, SIGMA);
        final double[][] goldStdImage  = FileOps.image2Array(FileOps.imgSrcDir + "/boat_soft_dwt_matlab.tif");
        MRA2D mra = new MRA2D(image, fb, 4, ConvolverFactory.ConvolutionType.FDCPU);
        MRA2D mrau = new MRA2DU(image, fb, 4, ConvolverFactory.ConvolutionType.FDCPU);
        ArrayList<MRA2D> mras = new ArrayList<MRA2D>();
        mras.add(mra);
        mras.add(mrau);
        System.out.println("-- Test MRA 2D Soft Thresholding ");
        int fileInd = 0;
        for (MRA2D mra2d : mras) {
            mra2d.dwt();
            threshold.visit(mra2d);
            mra2d.idwt();
            double[][] mraReconImage = mra2d.getFilteredData();
            FileOps.data2Image(mraReconImage, FileOps.imgWriteDir + "2d_soft_thresh_" +
                    Integer.toString(fileInd));
            double varianceMRA = new Variance().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(ArrayMath.deepCopy(goldStdImage), mraReconImage)));
            System.out.format("MRA %d Mean variance: %1.3f \n", fileInd, varianceMRA);
            Assert.assertEquals(0.0, varianceMRA, EPS);
            fileInd++;
        }
    }

    @Ignore
    public void SoftThresholdingTest3D() {
        double[][][] image = FileOps.nifti2Data(FileOps.nifti3D);
        FilterBank fb = Wavelets.getFarras();
        final double SIGMA = 350;
        Threshold threshold = new Threshold(Threshold.ThreshMeth.SOFT, SIGMA);
        final double[][][] goldStdImage  = FileOps.nifti2Data(FileOps.imgSrcDir + "/brain_soft_dwt_matlab.nii");
        MRA3D mra = new MRA3D(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        MRA3D mrau = new MRA3DU(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        ArrayList<MRA3D> mras = new ArrayList<MRA3D>();
        mras.add(mra);
        mras.add(mrau);
        System.out.println("-- Test MRA 3D Soft Thresholding ");
        int fileInd = 0;
        for (MRA3D mra3d : mras) {
            mra3d.dwt();
            threshold.visit(mra3d);
            mra3d.idwt();
            double[][][] mraReconImage = mra3d.getFilteredData();
            FileOps.data2Nifti(mraReconImage, FileOps.imgWriteDir + "3d_soft_thresh_" +
                    Integer.toString(fileInd));
            double varianceMRA = new Variance().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(ArrayMath.deepCopy(goldStdImage), mraReconImage)));
            System.out.format("MRA %d Mean variance: %1.3f \n", fileInd, varianceMRA);
            Assert.assertEquals(0.0, varianceMRA, EPS);
            fileInd++;
        }
    }
    
}

        

