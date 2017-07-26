package com.ericbarnhill.jmra;

public class Wavelets {
		
	static final double[] afLoNoTree = new double[] {
		0, 0,
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0.01122679215254, 0.01122679215254
	};
	
	static final double[] afHiNoTree = new double[] {
		-0.01122679215254, 0.01122679215254,
		0.08838834764832, 0.08838834764832,
		-0.69587998903400, 0.69587998903400,
		-0.08838834764832, -0.08838834764832,
		0, 0
	};
	
	static final double[] sfLoNoTree = new double[] {
		0.01122679215254, 0.01122679215254,
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0, 0
	};

	static final double[] sfHiNoTree = new double[] {
		0, 0, 
		-0.08838834764832,-0.08838834764832,
		0.69587998903400,-0.69587998903400,
		0.08838834764832, 0.08838834764832,
		0.01122679215254,-0.01122679215254
	};	
	
	static final double[] firstAFLoTree1 = new double[] {
		0, 
		-0.08838834764832, 0.08838834764832, 
		0.69587998903400, 0.69587998903400, 
		0.08838834764832, -0.08838834764832,
		0.01122679215254, 0.01122679215254, 
		0 
	};
	
	static final double[] firstAFHiTree1 = new double[] {
		0,  
		-0.01122679215254, 0.01122679215254,
		0.08838834764832, 0.08838834764832, 
		-0.69587998903400, 0.69587998903400,
		-0.08838834764832, -0.08838834764832,
		0		 
	};

	static final double[] firstSFLoTree1 = new double[] {
		0,
		0.01122679215254, 0.01122679215254, 
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0
	};
	
	static final double[] firstSFHiTree1 = new double[] {
		0,
		-0.08838834764832, -0.08838834764832,
		0.69587998903400, -0.69587998903400,
		0.08838834764832, 0.08838834764832, 
		0.01122679215254, -0.01122679215254,
		0		 
	};
	
	static final double[] firstAFLoTree2 = new double[] {
		0.01122679215254, 0.01122679215254,
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0, 0
	};
	
	static final double[] firstAFHiTree2 = new double[] {
		0, 0, 
		-0.08838834764832, -0.08838834764832, 
		0.69587998903400, -0.69587998903400,
		0.08838834764832, 0.08838834764832,
		0.01122679215254, -0.01122679215254		
	};
	
	static final double[] firstSFLoTree2 = new double[] {
		0, 0,
		-0.08838834764832, 0.08838834764832,
		0.69587998903400, 0.69587998903400,
		0.08838834764832, -0.08838834764832,
		0.01122679215254, 0.01122679215254		
	};
	
	static final double[] firstSFHiTree2 = new double[] {
		-0.01122679215254, 0.01122679215254,
		0.08838834764832, 0.08838834764832,
		-0.69587998903400, .69587998903400,
		-0.08838834764832, -0.08838834764832, 
		0, 0		
	};
	
	static final double[] afLoTree1 = new double[] {
		0.03516384000000, 0,
		-0.08832942000000, 0.23389032000000,
		0.76027237000000, 0.58751830000000,
		0, -0.11430184000000,
		0, 0
	};
	
	static final double[] afHiTree1 = new double[] {
		0,0,
		-0.11430184000000, 0,
		0.58751830000000, -0.7602723700000,
		0.23389032000000, 0.08832942000000,
		0, -0.03516384000000
	};
	
	static final double[] afLoTree2 = new double[] {
		0,0,
		-0.11430184000000, 0,
		0.58751830000000, 0.7602723700000,
		0.23389032000000, -0.08832942000000,
		0, 0.03516384000000
	};
	
	static final double[] afHiTree2 = new double[] {
		-0.03516384000000, 0,
		0.08832942000000, 0.23389032000000,
		-0.76027237000000, 0.58751830000000,
		0, -0.11430184000000,
		0, 0
	};
	
	static final double[] sfLoTree1 = new double[] {
		0,0,
		-0.11430184000000, 0,
		0.58751830000000, 0.7602723700000,
		0.23389032000000, -0.08832942000000,
		0, 0.03516384000000
	};
	
	static final double[] sfHiTree1 = new double[] {
		-0.03516384000000, 0,
		0.08832942000000, 0.23389032000000,
		-0.76027237000000, 0.58751830000000,
		0, -0.11430184000000,
		0, 0
	};
	
	static final double[] sfLoTree2 = new double[] {
		0.03516384000000, 0,
		-0.08832942000000, 0.23389032000000,
		0.76027237000000, 0.58751830000000,
		0, -0.11430184000000,
		0, 0
	};
	
	static final double[] sfHiTree2 = new double[] {
		0,0,
		-0.11430184000000, 0,
		0.58751830000000, -0.7602723700000,
		0.23389032000000, 0.08832942000000,
		0, -0.03516384000000
	};
	
	
}
