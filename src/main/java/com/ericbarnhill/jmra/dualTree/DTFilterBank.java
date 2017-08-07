package com.ericbarnhill.jmra.dualTree;

import com.ericbarnhill.jmra.filters.*;
import java.util.ArrayList;

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

}
