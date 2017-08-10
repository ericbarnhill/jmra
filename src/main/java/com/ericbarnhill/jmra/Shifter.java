package com.ericbarnhill.jmra;

import com.ericbarnhill.arrayMath.ArrayMath;
import org.apache.commons.numbers.complex.Complex;

public class Shifter {
	
	public static double[] circShift(double[] vector, int shift) {
		
		int adjN;
		double[] result = new double[vector.length];
		for (int n = 0; n < vector.length; n++) {
			adjN = (n - shift) % vector.length;
			if (adjN < 0) adjN = vector.length+adjN;
			result[n] = vector[adjN];
		}
		return result;
	}
	
	public static Complex[] circShift(Complex[] vector, int shift) {
		
		int adjN;
		Complex[] result = new Complex[vector.length];
		for (int n = 0; n < vector.length; n++) {
			adjN = (n - shift) % vector.length;
			if (adjN < 0) adjN = vector.length+adjN;
			result[n] = vector[adjN];
		}
		return result;
	}
	
	public static double[][] circShift(double[][] f, int shift) {
		int fi = f.length;
		int fj = f[0].length;
		double[][] r = new double[fi][fj];
		for (int i = 0; i < fi; i++) {
			r[i] = circShift(f[i], shift);
		}
		return r;
	}
	
	public static double[][] circShift(double[][] f, int shift1, int shift2) {
		int fi = f.length;
		int fj = f[0].length;
		double[][] r = new double[fi][fj];
		for (int i = 0; i < fi; i++) {
			r[i] = circShift(f[i], shift1);
		}
        double[][] rShift = ArrayMath.shiftDim(r);
        for (int j = 0; j < fj; j++) {
            rShift[j] = circShift(rShift[j], shift2);
        }
		return ArrayMath.shiftDim(rShift);
	}
	
	public static Complex[][] circShift(Complex[][] f, int shift1, int shift2) {
		int fi = f.length;
		int fj = f[0].length;
		Complex[][] r = new Complex[fi][fj];
		for (int i = 0; i < fi; i++) {
			r[i] = circShift(f[i], shift1);
		}
        Complex[][] rShift = ArrayMath.shiftDim(r);
        for (int j = 0; j < fj; j++) {
            rShift[j] = circShift(rShift[j], shift2);
        }
		return ArrayMath.shiftDim(rShift);
	}
	

}
