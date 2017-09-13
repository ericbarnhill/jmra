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
