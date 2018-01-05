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

package com.ericbarnhill.jmra.dualTree;

import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;

/** Creates dual tree filter banks. */
public class DTFilterBank extends FilterBank {

    final public ArrayList<FilterPair> faf;
    final public ArrayList<FilterPair> fsf;
    final public ArrayList<FilterPair> af;
    final public ArrayList<FilterPair> sf;

    public DTFilterBank(ArrayList<FilterPair> faf, ArrayList<FilterPair> fsf, ArrayList<FilterPair> af, ArrayList<FilterPair> sf) {
        this.faf = faf;
        this.fsf = fsf;
        this.af = af;
        this.sf = sf;
    }

    public DTFilterBank(FilterPair faf1, FilterPair faf2, FilterPair fsf1, FilterPair fsf2, FilterPair af1, FilterPair af2, FilterPair sf1, FilterPair sf2) {
        faf = new ArrayList<FilterPair>();
        faf.add(faf1);
        faf.add(faf2);
        fsf = new ArrayList<FilterPair>();
        fsf.add(fsf1);
        fsf.add(fsf2);
        af = new ArrayList<FilterPair>();
        af.add(af1);
        af.add(af2);
        sf = new ArrayList<FilterPair>();
        sf.add(sf1);
        sf.add(sf2);
    }

    public DTFilterBank(FilterPair faf1, FilterPair faf2, FilterPair faf3, FilterPair fsf1, FilterPair fsf2, FilterPair fsf3, FilterPair af1, FilterPair af2, FilterPair af3, FilterPair sf1, FilterPair sf2, FilterPair sf3) {
        faf = new ArrayList<FilterPair>();
        faf.add(faf1);
        faf.add(faf2);
        faf.add(faf3);
        fsf = new ArrayList<FilterPair>();
        fsf.add(fsf1);
        fsf.add(fsf2);
        fsf.add(fsf3);
        af = new ArrayList<FilterPair>();
        af.add(af1);
        af.add(af2);
        af.add(af3);
        sf = new ArrayList<FilterPair>();
        sf.add(sf1);
        sf.add(sf2);
        sf.add(sf3);
    }

}
