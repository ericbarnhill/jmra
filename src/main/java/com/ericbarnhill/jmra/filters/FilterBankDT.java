package com.ericbarnhill.jmra.filters;

import com.ericbarnhill.jmra.*;
import java.util.ArrayList;

public class FilterBankDT<V> {

    FilterPair faf1;
    FilterPair faf2;
    FilterPair af1;
    FilterPair af2;
    FilterPair fsf1;
    FilterPair fsf2;
    FilterPair sf1;
    FilterPair sf2;

    public FilterBankDT() {};

    public FilterBankDT(DTFilterSet dafp, DTFilterSet dsfp) {
        this.faf1 = dafp.ff1;
        this.faf2 = dafp.ff2;
        this.af1 = dafp.f1;
        this.af2 = dafp.f2;
        this.fsf1 = dsfp.ff1;
        this.fsf2 = dsfp.ff2;
        this.sf1 = dsfp.f1;
        this.sf2 = dsfp.f2;
    }



}
