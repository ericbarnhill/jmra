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
import org.apache.commons.math4.stat.descriptive.moment.Mean;
import org.apache.commons.math4.stat.descriptive.rank.Median;

public class TestExactRecovery{

static final double EPS = 0.00001; // 1e-5

    @Test
    public void MRA2DTest() {
        double[][] image = FileOps.image2Array(FileOps.image2D);
        FilterBank fb = Wavelets.getFarras();
        MRA2D mra = new MRA2D(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        mra.idwt();
        double[][] reconImage = mra.getFilteredData();
        double meanDifference = new Mean().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(image, reconImage)));
        Assert.assertEquals(0.0, meanDifference, EPS);
        System.out.println("-- Test MRA 2D Exact Recovery");
        System.out.format("%s: %1.3f \n", "---- Mean difference: ", meanDifference);
    }

    @Test
    public void MRA3DTest() {
        double[][][] image = FileOps.nifti2Data(FileOps.nifti3D);
        FilterBank fb = Wavelets.getFarras();
        MRA3D mra = new MRA3D(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        mra.idwt();
        double[][][] reconImage = mra.getFilteredData();
        FileOps.data2Image(reconImage, FileOps.imgDir+"/exact_recovery_recon_mra_3d.tif");
        double meanDifference = new Mean().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(image, reconImage)));
        Assert.assertEquals(0.0, meanDifference, 0.317);
        System.out.println("-- Test MRA 3D Exact Recovery. MATLAB gold std.: 0.3162");
        System.out.format("%s: %1.3f \n", "---- Mean difference: ", meanDifference);
    }

    @Test
    public void MRA2DUTest() {
        double[][] image = FileOps.image2Array(FileOps.image2D);
        FilterBank fb = Wavelets.getFarras();
        MRA2DU mra = new MRA2DU(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        mra.idwt();
        double[][] reconImage = mra.getFilteredData();
        double meanDifference = new Mean().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(image, reconImage)));
        Assert.assertEquals(0.0, meanDifference, EPS);
        System.out.println("-- Test Undecimated MRA 2D Exact Recovery");
        System.out.format("%s: %1.3f \n", "---- Mean difference: ", meanDifference);
    }

    @Test
    public void MRA3DUTest() {
        double[][][] image = FileOps.nifti2Data(FileOps.nifti3D);
        FilterBank fb = Wavelets.getFarras();
        MRA3DU mra = new MRA3DU(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        mra.dwt();
        mra.idwt();
        double[][][] reconImage = mra.getFilteredData();
        FileOps.data2Image(reconImage, FileOps.imgDir+"/exact_recovery_recon_umra_3d.tif");
        double meanDifference = new Mean().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(image, reconImage)));
        Assert.assertEquals(0.0, meanDifference, EPS);
        System.out.println("-- Test Undecimated MRA 3D Exact Recovery");
        System.out.format("%s: %1.3f \n", "---- Mean difference: ", meanDifference);
    }
    
    @Test
    public void MRA2DDTTest() {
        double[][] image = FileOps.image2Array(FileOps.image2D);
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        DualTree2DCplx dt2d = new DualTree2DCplx(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        dt2d.dwt();
        dt2d.idwt();
        double[][] reconImage = dt2d.getFilteredData();
        double meanDifference = new Mean().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(image, reconImage)));
        FileOps.data2Image(reconImage, FileOps.imgDir+"/exact_recovery_recon_dt_2d.tif");
        Assert.assertEquals(0.0, meanDifference, EPS);
        System.out.println("-- Test Dual-Tree MRA 2D Exact Recovery");
        System.out.format("%s: %1.3f \n", "---- Mean difference: ", meanDifference);
    }

    @Ignore
    public void MRA3DDTTest() {
        double[][][] image = FileOps.nifti2Data(FileOps.nifti3D);
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        DualTree3DCplx dt3d = new DualTree3DCplx(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU);
        dt3d.dwt();
        dt3d.idwt();
        double[][][] reconImage = dt3d.getFilteredData();
        FileOps.data2Image(reconImage, FileOps.imgDir+"/exact_recovery_recon_dt_3d.tif");
        double meanDifference = new Mean().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(image, reconImage)));
        Assert.assertEquals(0.0, meanDifference, EPS);
        System.out.println("-- Test Dual-Tree MRA 3D Exact Recovery");
        System.out.format("%s: %1.3f \n", "---- Mean difference: ", meanDifference);
    }
    
    @Test
    public void MRA2DDTUTest() {
        double[][] image = FileOps.image2Array(FileOps.image2D);
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        DualTree2DCplx dt2d = new DualTree2DCplx(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU, true);
        dt2d.dwt();
        dt2d.idwt();
        double[][] reconImage = dt2d.getFilteredData();
        double meanDifference = new Mean().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(image, reconImage)));
        Assert.assertEquals(0.0, meanDifference, EPS);
        System.out.println("-- Test Dual-Tree MRA 2D Undec Exact Recovery");
        System.out.format("%s: %1.3f \n", "---- Mean difference: ", meanDifference);
    }

    @Ignore
    public void MRA3DDTUTest() {
        double[][][] image = FileOps.nifti2Data(FileOps.nifti3D);
        DTFilterBank fb = Wavelets.getFarrasKingsbury();
        DualTree3DCplx dt3d = new DualTree3DCplx(image, fb, 3, ConvolverFactory.ConvolutionType.FDCPU, true);
        dt3d.dwt();
        dt3d.idwt();
        double[][][] reconImage = dt3d.getFilteredData();
        double meanDifference = new Mean().evaluate(
                ArrayMath.vectorize(
                    ArrayMath.subtract(image, reconImage)));
        Assert.assertEquals(0.0, meanDifference, EPS);
        System.out.println("-- Test Dual-Tree MRA 3D Exact Recovery");
        System.out.format("%s: %1.3f \n", "---- Mean difference: ", meanDifference);
    }
}

        

