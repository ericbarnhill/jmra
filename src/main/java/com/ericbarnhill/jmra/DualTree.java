package com.ericbarnhill.jmra;
import java.util.ArrayList; 
import com.ericbarnhill.jvcl.*; 

public class DualTree<N, B, V> {

    MRA<N, B, V> Tree1R;
    MRA<N, B, V> Tree1I;
    MRA<N, B, V> Tree2R;
    MRA<N, B, V> Tree2I;

    public DualTree(N originalData, B maskData, ArrayList<ArrayList<ArrayList<ArrayList<V>>>> filterBank, int decompLvls, ConvolverFactory.ConvolutionType convolutionType) {
        ArrayList<ArrayList<ArrayList<V>>> afBank = filterBank.get(0);
        ArrayList<ArrayList<ArrayList<V>>> sfBank = filterBank.get(1);
        ArrayList<ArrayList<V>> faf = afBank.get(0);
        ArrayList<ArrayList<V>> af = afBank.get(1);
        ArrayList<ArrayList<V>> fsf = sfBank.get(0);
        ArrayList<ArrayList<V>> sf = sfBank.get(1);
    }

        

}
