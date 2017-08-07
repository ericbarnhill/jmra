package com.ericbarnhill.jmra.filters;

import com.ericbarnhill.jmra.*;
import com.ericbarnhill.jmra.dualTree.*;

public class Wavelets {


    public static final FilterBank getFarras() {
        FilterPair af = new FilterPair(afLoNoTree, afHiNoTree);
        FilterPair sf = new FilterPair(sfLoNoTree, sfHiNoTree);
        FilterBank farras = new FilterBank(af, sf);
        return farras;
    }

    public static final DTFilterBank getFarrasKingsbury() {
        FilterPair faf1 = new FilterPair(firstAFLoTree1, firstAFHiTree1);
        FilterPair faf2 = new FilterPair(firstAFLoTree2, firstAFHiTree2);
        FilterPair fsf1 = new FilterPair(firstSFLoTree1, firstSFHiTree1);
        FilterPair fsf2 = new FilterPair(firstSFLoTree2, firstSFHiTree2);
        FilterPair af1 = new FilterPair(afLoTree1, afHiTree1);
        FilterPair af2 = new FilterPair(afLoTree2, afHiTree2);
        FilterPair sf1 = new FilterPair(sfLoTree1, sfHiTree1);
        FilterPair sf2 = new FilterPair(sfLoTree2, sfHiTree2);
        DTFilterBank farrasKingsbury = new DTFilterBank(faf1, faf2, fsf1, fsf2, af1, af2, sf1, sf2);
        return farrasKingsbury;
    }

		
	public static final double[] afLoNoTree = new double[] {
		0, 0,
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0.01122679215254, 0.01122679215254
	};
	
	public static final double[] afHiNoTree = new double[] {
		-0.01122679215254, 0.01122679215254,
		0.08838834764832, 0.08838834764832,
		-0.69587998903400, 0.69587998903400,
		-0.08838834764832, -0.08838834764832,
		0, 0
	};
	
	public static final double[] sfLoNoTree = new double[] {
		0.01122679215254, 0.01122679215254,
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0, 0
	};

	public static final double[] sfHiNoTree = new double[] {
		0, 0, 
		-0.08838834764832,-0.08838834764832,
		0.69587998903400,-0.69587998903400,
		0.08838834764832, 0.08838834764832,
		0.01122679215254,-0.01122679215254
	};	
	
	public static final double[] firstAFLoTree1 = new double[] {
		0, 
		-0.08838834764832, 0.08838834764832, 
		0.69587998903400, 0.69587998903400, 
		0.08838834764832, -0.08838834764832,
		0.01122679215254, 0.01122679215254, 
		0 
	};
	
	public static final double[] firstAFHiTree1 = new double[] {
		0,  
		-0.01122679215254, 0.01122679215254,
		0.08838834764832, 0.08838834764832, 
		-0.69587998903400, 0.69587998903400,
		-0.08838834764832, -0.08838834764832,
		0		 
	};

	public static final double[] firstSFLoTree1 = new double[] {
		0,
		0.01122679215254, 0.01122679215254, 
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0
	};
	
	public static final double[] firstSFHiTree1 = new double[] {
		0,
		-0.08838834764832, -0.08838834764832,
		0.69587998903400, -0.69587998903400,
		0.08838834764832, 0.08838834764832, 
		0.01122679215254, -0.01122679215254,
		0		 
	};
	
	public static final double[] firstAFLoTree2 = new double[] {
		0.01122679215254, 0.01122679215254,
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0, 0
	};
	
	public static final double[] firstAFHiTree2 = new double[] {
		0, 0, 
		-0.08838834764832, -0.08838834764832, 
		0.69587998903400, -0.69587998903400,
		0.08838834764832, 0.08838834764832,
		0.01122679215254, -0.01122679215254		
	};
	
	public static final double[] firstSFLoTree2 = new double[] {
		0, 0,
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0.01122679215254, 0.01122679215254		
	};
	
	public static final double[] firstSFHiTree2 = new double[] {
		-0.01122679215254, 0.01122679215254,
		0.08838834764832, 0.08838834764832,
		-0.69587998903400, .69587998903400,
		-0.08838834764832, -0.08838834764832, 
		0, 0		
	};
	
	public static final double[] afLoTree1 = new double[] {
		0.03516384000000, 0,
		-0.08832942000000, 0.23389032000000,
		0.76027237000000, 0.58751830000000,
		0, -0.11430184000000,
		0, 0
	};
	
	public static final double[] afHiTree1 = new double[] {
		0,0,
		-0.11430184000000, 0,
		0.58751830000000, -0.7602723700000,
		0.23389032000000, 0.08832942000000,
		0, -0.03516384000000
	};
	
	public static final double[] afLoTree2 = new double[] {
		0,0,
		-0.11430184000000, 0,
		0.58751830000000, 0.7602723700000,
		0.23389032000000, -0.08832942000000,
		0, 0.03516384000000
	};
	
	public static final double[] afHiTree2 = new double[] {
		-0.03516384000000, 0,
		0.08832942000000, 0.23389032000000,
		-0.76027237000000, 0.58751830000000,
		0, -0.11430184000000,
		0, 0
	};
	
	public static final double[] sfLoTree1 = new double[] {
		0,0,
		-0.11430184000000, 0,
		0.58751830000000, 0.7602723700000,
		0.23389032000000, -0.08832942000000,
		0, 0.03516384000000
	};
	
	public static final double[] sfHiTree1 = new double[] {
		-0.03516384000000, 0,
		0.08832942000000, 0.23389032000000,
		-0.76027237000000, 0.58751830000000,
		0, -0.11430184000000,
		0, 0
	};
	
	public static final double[] sfLoTree2 = new double[] {
		0.03516384000000, 0,
		-0.08832942000000, 0.23389032000000,
		0.76027237000000, 0.58751830000000,
		0, -0.11430184000000,
		0, 0
	};
	
	public static final double[] sfHiTree2 = new double[] {
		0,0,
		-0.11430184000000, 0,
		0.58751830000000, -0.7602723700000,
		0.23389032000000, 0.08832942000000,
		0, -0.03516384000000
	};
	
	
}
