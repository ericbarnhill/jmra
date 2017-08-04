package com.ericbarnhill.jmra.filters;
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

}
