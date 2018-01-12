# jmra
Java Multi-Resolution Analysis - Fast Wavelet Processing in java

(c) Eric Barnhill 2018. GNU public license.

This library aims to be the first library focused on high-performance, high-throughput wavelet analysis, particularly directed at medical imaging applications with large data sets. The library currently supports:

- 1D, 2D and 3D signals / images
- Both real- and complex- valued images
- Decimated and undecimated wavelet transforms
- Real, complex, and complex dual-tree wavelet analyses
- ImageJ TIFF and NIfTI data formats

The library also fuses object-oriented programming with the natural mathematical hierarchy of various wavelet-family transforms. Different transforms, and different-dimension data, are extensions of the fundamental MRA template. The DualTree template incorporates MRAs for trees. Thresholding is done with visitor classes.

Separate branches contain development for N-dimensional wavelet analysis and shearlet amalysis.

JMRA methods can be called from other Java code, or from MATLAB. Results can be visualizsed easily with ImageJ or MATLAB.

A JUnit testing suite and API docs are provided.
