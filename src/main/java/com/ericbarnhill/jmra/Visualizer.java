package com.ericbarnhill.jmra;

import com.ericbarnhill.jmra.dualTree.*;

public class Visualizer {

    public static void dumpDecomposition(DualTree3DCplx dt3dc) {
        for (int i = 0; i < dt3dc.trees.size(); i++) {
            MRA<double[][][], boolean[][][], double[]> tree = dt3dc.trees.get(i);
            for (int j = 0; j < tree.waveletData.size(); j++) {
                double[][][] data = tree.waveletData.get(j);
                String filename = FilePaths.root + "_" + Integer.toString(i) + "_" + 
                    Integer.toString(j) + ".tif";
                FilePaths.data2File(data, filename);
            }
        }
    }
}


