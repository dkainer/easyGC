EasyGC
======

This is my high throughput, command-line GC-MS analysis pipeline. It is built on top of the PyMS python library (http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115)

The pipeline requires the installation of Python 2.7 and all PyMS dependencies (see PyMS_UserGuide.pdf for dependencies instructions). The PyMS lib in this repository (/lib) should be used instead of the original PyMS library, as it has a few extra features and bug fixes.

The pipeline functions are found in GCMSalign.py. Currently parameters and paths are 'hardcoded' and functions are called from the python console. This will change in the near future with a command line interface being implemented.

TO DO
-----
- Create a command-line interface to the pipeline functions 
- include dependencies in repository
- produce excel file output with error checking results shown by cell colour.
