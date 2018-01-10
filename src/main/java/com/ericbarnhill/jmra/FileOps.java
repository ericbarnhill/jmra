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

/** Utilities for test suite. */
public class FileOps {

    public static String imgDir = getImagePath();
    public static String image2D = imgDir + "lena.tif";
    public static String nifti3D = imgDir + "test.nii";

    public static String getJarPath() {
        String jarDir = "";
        try {
            CodeSource codeSource = FileOps.class.getProtectionDomain().getCodeSource();
            File jarFile = new File(codeSource.getLocation().toURI().getPath());
            jarDir = jarFile.getParentFile().getPath();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(jarDir);
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
        sb.append(new String("src/test/img" + File.separator));
        System.out.println(sb.toString());
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
        fs.saveAsTiff(path);
    }

    static public void data2Nifti(double[][][] data, String path) {
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
        fs.saveAsTiff(path);
    }

}
