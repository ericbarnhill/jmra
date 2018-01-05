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

import com.ericbarnhill.jmra.dualTree.*;

/** Various methods for visualizing results without leaving Java,
 * implemented using a Visitor pattern. This class is in progress.
 */
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


