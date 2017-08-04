package com.ericbarnhill.jmra.filters;

public class FilterBank {
    
    public FilterPair af;
    public FilterPair sf;

    public FilterBank() {} // to add fields manually
    public FilterBank(FilterPair af, FilterPair sf) {
        this.af = af;
        this.sf = sf;
    }

}
