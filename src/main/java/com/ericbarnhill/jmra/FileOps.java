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

import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import java.security.CodeSource;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import com.ericbarnhill.niftijio.*;
import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jmra.dualTree.*;

/** File operation utilities to easily handle ImageJ and NIfTI formats. */
public class FileOps {

    public static String imgSrcDir = getImagePath() + "goldStd/";
    public static String imgWriteDir = getImagePath();
    public static String image2D = imgSrcDir + "boat_noisy.tif";
    public static String nifti3D = imgSrcDir + "brain_noisy.nii";

    public static String getJarPath() {
        String jarDir = "";
        try {
            CodeSource codeSource = FileOps.class.getProtectionDomain().getCodeSource();
            File jarFile = new File(codeSource.getLocation().toURI().getPath());
            jarDir = jarFile.getParentFile().getPath();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return jarDir;
    }

    public static String getImagePath() {
        String jarDir = getJarPath();
        String[] splitPath = jarDir.split(File.separator);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < splitPath.length - 1; i++) {
            sb.append(splitPath[i]);
            sb.append(File.separator);
        }
        sb.append(new String("src/test/img/" + File.separator));
        return sb.toString();
    }

    public static double[][] image2Array(String path) {
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

    /** 
     * Dumps 2D wavelet decomposition as a series of images.
     * These images will be automatically numbered.
     * Path provided should include directory and filename stem
     * but not file extension.
     */
    public static void decomposition2Image(MRA2D m, String path) {
        ArrayList<double[][]> decomposition = m.getDecomposition();
        for (int i = 0; i < decomposition.size(); i++) {
            String numberedPath = path + String.format("_%03d.tif", i);
            data2Image(decomposition.get(i), numberedPath);
        }
    }

    /** 
     * Dumps 2D dual-tree wavelet decomposition as a series of images.
     * These images will be automatically numbered.
     * Path provided should include directory and filename stem
     * but not file extension.
     */
    public static void dualtreeDecomposition2Image(DualTree2DCplx m, String path) {
        for (int i = 0; i < m.trees.size(); i++) {
            String treePath = path + String.format("_tree_%d", i);
            decomposition2Image((MRA2D)m.trees.get(i), treePath);
        }
    }

    /** 
     * Dumps 3D wavelet decomposition as a series of images.
     * These images will be automatically numbered.
     * Path provided should include directory and filename stem
     * but not file extension.
     */
    public static void decomposition2Nifti(MRA3D m, String path) {
        ArrayList<double[][][]> decomposition = m.getDecomposition();
        for (int i = 0; i < decomposition.size(); i++) {
            String numberedPath = path + String.format("_%03d.nii", i);
            data2Nifti(decomposition.get(i), numberedPath);
        }
    }

    /** 
     * Dumps 2D dual-tree wavelet decomposition as a series of images.
     * These images will be automatically numbered.
     * Path provided should include directory and filename stem
     * but not file extension.
     */
    public static void dualtreeDecomposition2Nifti(DualTree3DCplx m, String path) {
        for (int i = 0; i < m.trees.size(); i++) {
            String treePath = path + String.format("_tree_%d", i);
            decomposition2Nifti((MRA3D)m.trees.get(i), treePath);
        }
    }

    public static void data2Image(double[][] array, String path) {
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
        fs.saveAsTiff(path+".tif");
    }

    public static void data2Image(double[][][] data, String path) {
        int w = data.length;
        int h = data[0].length;
        int d = data[0][0].length;
        ImageStack is = new ImageStack(w,h);
        for (int k = 0; k < d; k++) {
            FloatProcessor fp = new FloatProcessor(w,h);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                            fp.putPixelValue(i,j,data[i][j][k]);
                }
            }
            is.addSlice(fp);
        }
        ImagePlus ip = new ImagePlus("", is);
        FileSaver fs = new FileSaver(ip);
        fs.saveAsTiff(path+".tif");
    }

    public static void data2Nifti(double[][][] data, String path) {
        NiftiVolume nv = new NiftiVolume(ArrayMath.convert3dto4d(data));
        try {
            nv.write(path + ".nii");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double[][][] nifti2Data(String path) {
        NiftiVolume nv = null;
        try {
            nv  = NiftiVolume.read(path);
        } catch (IOException e) {
            e.printStackTrace();
        }
        double[][][][] niftiArray = nv.data.toArray();
        double[][][] image = ArrayMath.convert4dto3d(niftiArray);
        return image;
    }

    public static void main(String [] args) {
        System.out.println("Dummy main method to enable maven shade plugin to bundle all dependencies as a fat JAR.");
    }

}
