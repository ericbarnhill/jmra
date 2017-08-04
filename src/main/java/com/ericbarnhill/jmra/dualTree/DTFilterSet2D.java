package com.ericbarnhill.jmra.filters;

public class DTFilterSet2D extends DTFilterBank {

    final public FilterPair faf1;
    final public FilterPair faf2;
    final public FilterPair fsf1;
    final public FilterPair fsf2;
    final public FilterPair af1;
    final public FilterPair af2;
    final public FilterPair sf1;
    final public FilterPair sf2;

    public DTFilterSet2D(FilterPair faf1, FilterPair faf2, FilterPair fsf1, FilterPair fsf2, FilterPair af1, FilterPair af2, FilterPair sf1, FilterPair sf2) {
        this.faf1 = faf1;
        this.faf2 = faf2;
        this.fsf1 = fsf1;
        this.fsf2 = fsf2;
        this.af1 = af1;
        this.af2 = af2;
        this.sf1 = sf1;
        this.sf2 = sf2;
    }

}
