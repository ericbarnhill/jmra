package com.ericbarnhill.jmra;

import org.apache.commons.lang3.ArrayUtils;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import com.ericbarnhill.arrayMath.ArrayMath;
import org.apache.commons.math4.stat.descriptive.rank.Max;
import org.apache.commons.math4.stat.descriptive.rank.Median;

class Threshold {

    public static enum ThreshMeth {
        HARD, SOFT, NNG
    }

    public static enum NoiseEstMeth {
        VISU_SHRINK, SURE_SHRINK, BAYES_SHRINK
    }

    public static double[] threshold(double[] data, boolean[] mask, ThreshMeth threshMeth, NoiseEstMeth noiseEstMeth) {
        double sigma = estimateSigma(data, mask, noiseEstMeth);
        double[] thresholdedPixels = applyThresh(data, sigma, threshMeth);
        return thresholdedPixels;
    }

    public static double[] applyThresh(double[] data, double sigma, ThreshMeth threshMeth) {
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



    public static double[] getMaskedPixels(double[] data, boolean[] mask) {
        List<Double> maskedPixelsList = new ArrayList<Double>();
        for (int i = 0; i < data.length; i++) {
            if (mask[i]) {
                maskedPixelsList.add(data[i]);
            }
        }
        return ArrayUtils.toPrimitive(maskedPixelsList.toArray(new Double[0]));
    }

    public static double estimateSigma(double[] data, boolean[] mask,  NoiseEstMeth noiseEstMeth) {
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

    private static double universalThreshold(double[] maskedPixels, double noiseEst) {
        return noiseEst*Math.sqrt(2*Math.log(maskedPixels.length));
    }

    private static double noiseEst(double[] pixels) {
        return new Median().evaluate(ArrayMath.abs(pixels))/0.6745;
    }

    private static double[] maskPixels(double[] data, boolean[] mask) {
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
