EasyGC
======

This is a high throughput, command-line GC-MS analysis pipeline. It is built on top of a modified version of the PyMS python library (http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115) 

This pipeline makes it easy to analyse a large set of GC-MS runs. The input is a directory full of runs (in CDF or JDX format), and the output is a matrix of peaks aligned across all your samples, including their retention times, TIC areas, and mini mass spec per called peak. You can tweak the way peaks are called, filtered and aligned with a whole range of parameters. The key thing, though, is that this pipeline makes it VERY easy to quantitate a lot of peaks across a lot of samples with minimal fuss. It is especially useful for population-level analyses where relative peak size or presence matters most, rather than extremely accurate peak identification and quantitation.

**version: 0.0.1**: very early. There are currently bugs in the multithreading on linux.

prerequisites
-------------
- Python 2.7
- matplotlib
- netCDF
- pycdf 0.6-3  (this will only work on Linux)
- scipy.ndimage package
 
(see PyMS_UserGuide.pdf for dependencies installation instructions). 

The PyMS lib in this repository (/pyms) should be used instead of the original PyMS python library, as it has a few extra features and bug fixes.

manual
------
usage: python easyGC.py [command] [options]

**peakcall**

run the peak caller on a directory of GC-MS files. This will produce a <i>.expr</i> file output for each GC-MS run which is a binary format of called peaks and mass specs. The .expr files are used later by the aligner. 
```
usage: easyGC peakcall [-h] -i INDIR -f FTYPE [-TS TRIMSTART] [-TE TRIMEND]
                       [-W WINDOW] [-S SCANS] [-N MINIONS] [-R MININTENSITY]
                       [-M NOISEMULT] [-I TOPIONS]
```


```
optional arguments:
  -h, --help            show this help message and exit
  -i INDIR, --indir INDIR
                        directory containing your GC-MS files to be processed
                        (default: None)
  -f FTYPE, --ftype FTYPE
                        CDF or JDX. This is the type of input GC-MS files you
                        have. CDF is not supported on Windows (default: None)
  -TS TRIMSTART, --trimstart TRIMSTART
                        time in minutes (X.XX) in the chromatogram from where
                        the analysis should begin. Helps to cut out junk at
                        the start (default: 0.0)
  -TE TRIMEND, --trimend TRIMEND
                        time in minutes (X.XX) in the chromatogram where the
                        analysis should end. Helps to cut out junk at the end
                        (default: 20.0)
  -W WINDOW, --window WINDOW
                        peak calling: width (in scans) of window over which
                        local ion maxima are detected. Should be similar to
                        the width off your peaks. (default: 9)
  -S SCANS, --scans SCANS
                        peak calling: distance (in scans) at which locally
                        apexing ions can be combined into one peak (default:
                        3)
  -N MINIONS, --minions MINIONS
                        peak calling: min number of apexing ions with
                        intensity above a threshold required for a peak to be
                        called. Higher = less peaks called (default: 4)
  -R MININTENSITY, --minintensity MININTENSITY
                        peak calling: min intensity (percent) of an ion
                        relative to max peak intensity for that ion to be
                        included in the peak (default: 5)
  -M NOISEMULT, --noisemult NOISEMULT
                        peak calling: total peak intensity must be at least
                        this multiple of the base noise level to be called.
                        Higher multiple means fewer peaks called (default: 4)
  -I TOPIONS, --topions TOPIONS
                        from the list of most important ions in a peak, how
                        many should be outputted as a mini mass-spec?
                        (default: 10)
```

**align**

run the peak aligner on a directory of .expr files that were produced by [peakcall]() . This will produce three CSV files as ouput:

- aligned_rt.csv
- aligned_area.csv
- aligned_ions.csv

This command is especially useful if you are not happy with the aligned output from 'pipeline' and want to tweak the alignment parameters to see how they affaect your output matrix, without having to call peaks all over again. 

```
usage: easyGC align -e EXPRDIR [-D DISTANCE] [-G GAP] [-C MINCOMMON] [-T THREADS]
```
```
optional arguments:
  -h, --help            show this help message and exit
  -e EXPRDIR, --exprdir EXPRDIR
                        the path to the .expr files from a previous peak
                        calling run. These will be aligned. (default: None)
  -D DISTANCE, --distance DISTANCE
                        local alignment: distance in retention time (seconds)
                        over which the local peak aligner should search for
                        similar peaks to this one (default: 2.5)
  -G GAP, --gap GAP     local alignment: gap penalty. Lower G results in more
                        peaks in the output. Higher G result in fewer output
                        peaks but possibly some peaks contain multiple merged
                        peaks (default: 0.4)
  -C MINCOMMON, --mincommon MINCOMMON
                        local alignment: minimum number of samples that an
                        aligned peak must be called in for it to be outputted
                        (default: 1)
  -T THREADS, --threads THREADS
                        number of threads to use. Currently only multithreaded
                        on linux! (default: 1)

```

**pipeline**

run the whole shebang on a directory of GC-MS runs, inclduing peak calling through to aligned output.
```
usage: easyGC pipeline [-h] -i INDIR -f FTYPE [-TS TRIMSTART] [-TE TRIMEND]
                       [-W WINDOW] [-S SCANS] [-N MINIONS] [-R MININTENSITY]
                       [-M NOISEMULT] [-I TOPIONS] [-D DISTANCE] [-G GAP]
                       [-C MINCOMMON] [-T THREADS]
```


TO DO
-----
- include dependencies in repository
- produce excel file output with error checking results shown by cell colour.
