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
    public void decompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl);
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl; // starting point
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            ArrayList<double[][][]> loAndHi = getDecomposition(localIndex, ind, decompLvl, dimLvl, localStride);
            double[][][] lo = loAndHi.get(0);
            Serializer.writeData(lo, waveletTempFiles.get(ind));

            double[][][] hi = loAndHi.get(1);
            Serializer.writeData(hi, waveletTempFiles.get(ind + localPair));
        }
        if (dimLvl < dimLvls - 1) {
            decompose(decompLvl, dimLvl+1);
       }
    }

    @Override
    public ArrayList<double[][][]> getDecomposition(int localIndex, int ind, int decompLvl, int dimLvl, int localStride) {
        double[][][] x = new double[0][][];
        if (dimLvl == 0) {
            if (decompLvl == 0) {
                x = ArrayMath.deepCopy(paddedData);
            } else {
                x = Serializer.loadData(waveletTempFiles.get(localIndex - stride));
            }
        } else {
            x = Serializer.loadData(waveletTempFiles.get(ind));
        }
        switch(localStride) {
            case 4: 
                x = ArrayMath.shiftDim(x, 1);
                break;
            case 2:
                x = ArrayMath.shiftDim(x, 2);
                break;
        }
        double[][][] lo = AFB(x, fb.af.lo, decompLvl);
        double[][][] hi = AFB(x, fb.af.hi, decompLvl);
        switch(localStride) {
            case 4: 
                lo = ArrayMath.shiftDim(lo, 2);
                hi = ArrayMath.shiftDim(hi, 2);
                break;
            case 2:
                lo = ArrayMath.shiftDim(lo, 1);
                hi = ArrayMath.shiftDim(hi, 1);
                break;
        }
        ArrayList<double[][][]> loAndHi = new ArrayList<double[][][]>();
        loAndHi.add(lo);
        loAndHi.add(hi);
        return loAndHi;
    }

    public void recompose(int decompLvl, int dimLvl) {
        int localStride = (int)Math.pow(2, dimLvls - dimLvl);
        int localPair = localStride / 2;
        int localIndex = stride*decompLvl;
        for (int ind = localIndex; ind < localIndex+stride; ind += localStride) { 
            double[][][] y = getRecomposition(localPair, ind, decompLvl, dimLvl, localStride);
            Serializer.writeData(y, waveletTempFiles.get(ind));
        }
        if (dimLvl > 0) {
            recompose(decompLvl, dimLvl-1);
        }
        if (decompLvl > 0) {
            Serializer.writeData(Serializer.loadData(waveletTempFiles.get(stride*(decompLvl-1))), waveletTempFiles.get(stride*decompLvl));
        }
    }

    @Override
    public double[][][] getRecomposition(int localPair, int ind, int decompLvl, int dimLvl, int localStride) {
        double[][][] lo = Serializer.loadData(waveletTempFiles.get(ind));
        double[][][] hi = Serializer.loadData(waveletTempFiles.get(ind + localPair));
        switch (localStride) {
            case 4: 
                lo = ArrayMath.shiftDim(lo, 1);
                hi = ArrayMath.shiftDim(hi, 1);
                break;
            case 2:
                lo = ArrayMath.shiftDim(lo, 2);
                hi = ArrayMath.shiftDim(hi, 2);
                break;
        }
        double[][][] y = SFB(lo, hi, fb.sf.lo, fb.sf.hi, decompLvl);
        switch (localStride) {
            case 4: 
                y = ArrayMath.shiftDim(y, 2);
                break;
            case 2:
                y = ArrayMath.shiftDim(y, 1);
                break;
        }
        return y;
    }

    @Override
    public double[][][] getFilteredData() {
        //DEBUGGING 
        System.out.println("in wrong get filtered data method");
        return ArrayMath.stripBorderPadding(Serializer.loadData(waveletTempFiles.get(0)), (wPad-w)/2, (hPad-h)/2, (dPad-d)/2);
    }

}
