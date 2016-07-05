EasyGC
======

This is my high throughput, command-line GC-MS analysis pipeline. It is built on top of a modified version of the PyMS python library (http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115) 

This pipeline makes it easy to analyse a large set of GC-MS runs. The input is a directory full of runs (in CDF or JDX format), and the output is a matrix of peaks aligned across all your samples, including their retention times, TIC areas, and mini mass spec per called peak. You can tweak the way peaks are called, filtered and aligned with a whole range of parameters. The key thing, though, is that this pipeline makes it VERY easy to quantitate a lot of peaks across a lot of samples with minimal fuss. It is especially useful for population-level analyses where relative peak size or presence matters most, rather than extremely accurate peak identification and quantitation.

prerequisites
-------------
- Python 2.7
- matplotlib
- netCDF
- pycdf 0.6-3  (this will only work on Linux)
- scipy.ndimage package
 
(see PyMS_UserGuide.pdf for dependencies installation instructions). 

The PyMS lib in this repository (/pyms) should be used instead of the original PyMS python library, as it has a few extra features and bug fixes.


TO DO
-----
- include dependencies in repository
- produce excel file output with error checking results shown by cell colour.
