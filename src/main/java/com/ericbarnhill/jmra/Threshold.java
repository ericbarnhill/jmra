package com.ericbarnhill.jmra;

import org.apache.commons.lang3.ArrayUtils;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import com.ericbarnhill.arrayMath.ArrayMath;
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Median;
import org.apache.commons.numbers.complex.Complex;
import org.apache.commons.numbers.complex.ComplexUtils;
import com.ericbarnhill.jmra.dualTree.*;

public class Threshold {

    public  enum ThreshMeth {
        HARD, SOFT, NNG
    }

    public  enum NoiseEstMeth {
        VISU_SHRINK, SURE_SHRINK, BAYES_SHRINK
    }

    ThreshMeth threshMeth;
    NoiseEstMeth noiseEstMeth;
    
    public Threshold(ThreshMeth threshMeth, NoiseEstMeth noiseEstMeth) {
        this.threshMeth = threshMeth;
        this.noiseEstMeth = noiseEstMeth;
    }

    public void visit(MRA1D mra1d) {
        // loop through each subband, pass method
        for (int i = 0; i < mra1d.waveletData.size(); i++) {
            // avoid scaling datas
            if (i % mra1d.stride != 0) {
                int level = (int)Math.floor(i / mra1d.stride);
                boolean[] maskDownsampled = ArrayMath.decimate(mra1d.paddedMask, (int)Math.pow(2,level));
                mra1d.waveletData.set(i, threshold(mra1d.waveletData.get(i), maskDownsampled));
            }
        }
    }

    public void visit(MRA1DU mra1du) {
        // loop through each subband, pass method
        for (int i = 0; i < mra1du.waveletData.size(); i++) {
            // avoid scaling datas
            if (i % mra1du.stride != 0) {
                int level = (int)Math.floor(i / mra1du.stride);
                mra1du.waveletData.set(i, threshold(mra1du.waveletData.get(i), mra1du.maskData));
            }
        }
    }

    public void visit (MRA2D mra2d) {
        // loop through each subband, pass method
        for (int i = 0; i < mra2d.waveletData.size(); i++) {
            // avoid scaling datas
            if (i % mra2d.stride != 0) {
                int level = (int)Math.floor(i / mra2d.stride);
                int decimFac = (int)Math.pow(2, level+1);
                boolean[][] maskDownsampled = ArrayMath.decimate(mra2d.paddedMask, decimFac);
                double[] waveletVec = ArrayMath.vectorize(mra2d.waveletData.get(i));
                boolean[] maskVec = ArrayMath.vectorize(maskDownsampled);
                mra2d.waveletData.set(i, 
                    ArrayMath.devectorize(
                        threshold(
                            waveletVec, maskVec)
                        ,mra2d.w/decimFac)
                    );
            }
        }
    }

    public void visit(MRA2DU mra2du) {
        // loop through each subband, pass method
        for (int i = 0; i < mra2du.waveletData.size(); i++) {
            // avoid scaling datas
            if (i % mra2du.stride != 0) {
                int level = (int)Math.floor(i / mra2du.stride);
                double[] waveletVec = ArrayMath.vectorize(mra2du.waveletData.get(i));
                boolean[] maskVec = ArrayMath.vectorize(mra2du.maskData);
                mra2du.waveletData.set(i, 
                    ArrayMath.devectorize(
                        threshold(
                            waveletVec, maskVec)
                        ,mra2du.w)
                    );
            }
        }
    }

    public void visit(MRA3D mra3d) {
        // loop through each subband, pass method
        for (int i = 0; i < mra3d.waveletData.size(); i++) {
            // avoid scaling datas
            if (i % mra3d.stride != 0) {
                int level = (int)Math.floor(i / mra3d.stride);
                int decimFac = (int)Math.pow(2, level+1);
                boolean[][][] maskDownsampled = ArrayMath.decimate(mra3d.paddedMask, decimFac);
                double[] waveletVec = ArrayMath.vectorize(mra3d.waveletData.get(i));
                boolean[] maskVec = ArrayMath.vectorize(maskDownsampled);
                mra3d.waveletData.set(i, 
                    ArrayMath.devectorize(
                        threshold(
                            waveletVec, maskVec)
                        ,mra3d.wPad/decimFac, mra3d.hPad/decimFac)
                    );
            }
        }
    }

    public void visit(MRA3DU mra3du) {
        // loop through each subband, pass method
        for (int i = 0; i < mra3du.waveletData.size(); i++) {
            // avoid scaling datas
            if (i % mra3du.stride != 0) {
                int level = (int)Math.floor(i / mra3du.stride);
                double[] waveletVec = ArrayMath.vectorize(mra3du.waveletData.get(i));
                boolean[] maskVec = ArrayMath.vectorize(mra3du.maskData);
                mra3du.waveletData.set(i, 
                    ArrayMath.devectorize(
                        threshold(
                            waveletVec, maskVec)
                        ,mra3du.w, mra3du.h)
                    );
            }
        }
    }

    /*
    public void visit(DualTree1D dt1d) {
        for (int i = 0; i < dt1d.tree1R.waveletData.size(); i++) {
            if (i % dt1d.tree1R.stride != 0) {
                Complex[] tree1 = ComplexUtils.split2Complex(dt1d.tree1R.waveletData.get(i), dt1d.tree1I.waveletData.get(i));
                Complex[] tree2 = ComplexUtils.split2Complex(dt1d.tree2R.waveletData.get(i), dt1d.tree2I.waveletData.get(i));
                boolean[] mask = dt1d.tree1R.maskData;
                tree1 = threshold(tree1, mask);
                tree2 = threshold(tree2, mask);
                dt1d.tree1R.waveletData.set(i, ComplexUtils.complex2Real(tree1));
                dt1d.tree1I.waveletData.set(i, ComplexUtils.complex2Imaginary(tree1));
                dt1d.tree2R.waveletData.set(i, ComplexUtils.complex2Real(tree2));
                dt1d.tree2I.waveletData.set(i, ComplexUtils.complex2Real(tree2));
            }
        }
    }
     */           

    public  double[] threshold(double[] data, boolean[] mask) {
        double sigma = estimateSigma(data, mask);
        double[] thresholdedPixels = applyThresh(data, sigma);
        return thresholdedPixels;
    }

    // estimate the noise using the real component of the tree
    public Complex[] threshold(Complex[] data, boolean[] mask) {
        double sigma = estimateSigma(ComplexUtils.complex2Real(data), mask);
        Complex[] thresholdedPixels = applyThresh(data, sigma);
        return thresholdedPixels;
    }

    public  double[] applyThresh(double[] data, double sigma) {
        int w = data.length;
        double[] threshedImage = new double[w];
        for (int i = 0; i < w; i++) {
            switch(threshMeth) {
                case HARD:
                    if (Math.abs(data[i]) > sigma) {
                        threshedImage[i] = data[i];
                    }
                    break;
                case SOFT:
                    if (Math.abs(data[i]) > sigma) {
                        threshedImage[i] = (Math.abs(data[i]) - sigma)*Math.signum(data[i]);
                     }
                    break;
                case NNG:
                    if (Math.abs(data[i]) > sigma) {
                        threshedImage[i] = data[i] - sigma*sigma/data[i];
                    }
                    break;
            }
        }
        return threshedImage;
    }

    public  Complex[] applyThresh(Complex[] data, double sigma) {
        int w = data.length;
        Complex[] threshedImage = new Complex[w];
        for (int i = 0; i < w; i++) {
            switch(threshMeth) {
                case HARD:
                    if ( data[i].abs() > sigma) {
                        threshedImage[i] = data[i];
                    }
                    break;
                case SOFT:
                    if (data[i].abs() > sigma) {
                        double magnitude = (data[i].abs()- sigma);
                        threshedImage[i] = data[i].multiply(magnitude / (magnitude + sigma));
                     }
                    break;
                case NNG:
                    if (data[i].abs() > sigma) {
                        double magnitude = (data[i].abs()- sigma*sigma/data[i].abs());
                        threshedImage[i] = data[i].multiply(magnitude / (magnitude + sigma));
                    }
                    break;
            }
        }
        return threshedImage;
    }



    public  double[] getMaskedPixels(double[] data, boolean[] mask) {
        List<Double> maskedPixelsList = new ArrayList<Double>();
        for (int i = 0; i < data.length; i++) {
            if (mask[i]) {
                maskedPixelsList.add(data[i]);
            }
        }
        return ArrayUtils.toPrimitive(maskedPixelsList.toArray(new Double[0]));
    }

    public  double estimateSigma(double[] data, boolean[] mask) {
        double[] maskedPixels = maskPixels(data, mask);
        double sigma = 0;
        int N = maskedPixels.length;
            double robustNoiseEst = noiseEst(maskedPixels);
                switch (noiseEstMeth) {
            case VISU_SHRINK:
                sigma = universalThreshold(maskedPixels, robustNoiseEst);
                break;
            case SURE_SHRINK:
                // hybrid scheme of Donoho et al
                // If very sparse, evaluate with universal threshold. Otherwise use SURE
                double[] sortedPixels = ArrayMath.square(
                                            ArrayMath.abs(
                                                ArrayMath.divide(maskedPixels,robustNoiseEst)));
                Arrays.sort(sortedPixels);
                // Sparsity criterion
                double nuNum = ArrayMath.sum(sortedPixels);
                double nuDenom = Math.pow(Math.log(N)/Math.log(2), 1.5);
                double nu = nuNum / nuDenom;
                if (nu < 1) {
                    sigma = universalThreshold(maskedPixels, robustNoiseEst);
                } else {
                    double minRisk = Double.MAX_VALUE;
                    int minRiskIndex = -1;
                    double cumSum = 0;
                    for (int n = 0; n < N; n++){
                        cumSum += sortedPixels[n];
                        double pixelRisk = (N - 2*n + cumSum) / N;
                        if (pixelRisk < minRisk) {
                            minRisk = pixelRisk;
                            minRiskIndex = n;
                        }
                    }
                    sigma = Math.sqrt(sortedPixels[minRiskIndex])*robustNoiseEst;
                    System.out.format("%1.4f %1.4f %1.4f \n", robustNoiseEst, sigma, cumSum);
                }
                break;
            case BAYES_SHRINK:
                Median m = new Median();
                double sigHat = m.evaluate(ArrayMath.abs(maskedPixels)) / 0.6745;
                double sig2y = ArrayMath.sum(ArrayMath.square(maskedPixels)) / N;
                double sigX = Math.sqrt(Math.max(sig2y - sigHat*sigHat, 0));
                if (sigX == 0) {
                    sigma = new Max().evaluate(maskedPixels);
                } else {
                    sigma = sigHat*sigHat / sigX;
                }
                break;
        }
        return sigma;
    }

    private  double universalThreshold(double[] maskedPixels, double noiseEst) {
        return noiseEst*Math.sqrt(2*Math.log(maskedPixels.length));
    }

    private  double noiseEst(double[] pixels) {
        return new Median().evaluate(ArrayMath.abs(pixels))/0.6745;
    }

    private  double[] maskPixels(double[] data, boolean[] mask) {
        int numTrue = 0;
        for (int i = 0; i < mask.length; i++) {
            if (mask[i]) numTrue++;
        }
        double[] maskedPixels = new double[numTrue];
        int maskedPixelsIndex = 0;
        for (int i = 0; i < mask.length; i++) {
            if (mask[i]) {
                maskedPixels[maskedPixelsIndex++] = data[i];
            }
        }
        return maskedPixels;
    }

}
