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

import com.ericbarnhill.arrayMath.ArrayMath;
import com.ericbarnhill.jvcl.*;
import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.*;
import ij.io.Opener;
import ij.io.FileSaver;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

/** 
 * 3D multi-resolution analysis with serial data processing.
 * (Less intensive demands on RAM.)
 */
public class MRA3DSerial extends MRA3D {

    public ArrayList<File> waveletTempFiles;

    public MRA3DSerial() {
        super();
        initializeWaveletTempFiles();
    }
    
    public MRA3DSerial(double[][][] origData, boolean[][][] maskData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, filterBank, decompLvls, convType);
        initializeWaveletTempFiles();
    }
        
    public MRA3DSerial(double[][][] origData, FilterBank filterBank, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, filterBank, decompLvls, convType);
        initializeWaveletTempFiles();
    }
        
    public MRA3DSerial(double[][][] origData, boolean[][][] maskData, int decompLvls, ConvolverFactory.ConvolutionType convType) {
        super(origData, maskData, decompLvls, convType);
        initializeWaveletTempFiles();
    }

    private void initializeWaveletTempFiles() {
        waveletTempFiles = new ArrayList<File>();
        try {
            for (int n = 0; n < stride * decompLvls; n++) {
                File f = File.createTempFile("JMRA" + Integer.toString(n), ".tmp");
                f.deleteOnExit();
                waveletTempFiles.add(f);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public double[][][] getData(int index) {
        return Serializer.loadData(waveletTempFiles.get(index));
    }

    @Override
    public void setData(int index, double[][][] data) {
        Serializer.writeData(waveletTempFiles.get(index), data);
    }

}
