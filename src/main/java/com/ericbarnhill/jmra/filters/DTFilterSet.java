package com.ericbarnhill.jmra.filters;

public class DTFilterSet extends FilterBank {
    // DTFilterSet extends FilterBank since they are both used on MRA objects
    // The full DTFilterBank is used with DualTree objects

    public FilterPair ff1;
    public FilterPair ff2;
    public FilterPair f1;
    public FilterPair f2;
    
    public DTFilterSet(FilterPair ff1, FilterPair ff2, FilterPair f1, FilterPair f2) {
        this.ff1 = ff1;
        this.ff2 = ff2;
        this.f1 = f1;
        this.f2 = f2;
    }

}
