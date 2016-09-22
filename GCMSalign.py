# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 16:06:16 2013

@author: dkainer

@note: I Modified the PyMS Alignment class in pyms.Peak.List.DPA.class by adding
        a function write_ion_areas_csv(self, ms_file_name, minutes=True). This function
        outputs the top 20 ions and ion areas for each peak in the alignment matrix.

@note: Added a JCAMP-DX reader. This can read JDX files that are output from OpenChrom
"""

import sys, os, datetime, traceback, glob
import pyms.GCMS.IO.JCAMP.Function
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import *
from pyms.Baseline.TopHat import tophat
from pyms.Noise.Analysis import window_analyzer
from pyms.Peak.Class import Peak
from pyms.Peak.Function import peak_sum_area, peak_top_ion_areas

from pyms.Deconvolution.BillerBiemann.Function import BillerBiemann, rel_threshold, num_ions_threshold, sum_maxima

from pyms.Experiment.Class import Experiment
from pyms.Experiment.IO import store_expr, load_expr
from pyms.Peak.List.DPA.Class import PairwiseAlignment
from pyms.Peak.List.DPA.Function import align_with_tree, exprl2alignment

from pyms.Display.Function import *
from pyms.Display.Class import *

import multiprocessing as mp

#window = 9  # width of window over which local ion maxima are detected
#scans = 3  # distance at which locally apexing ions can be combined into one peak
#n = 4  # min number of ions with intensity above a threshold
#r = 5  # min percentage of mass intensity relative to max peak intensity
#top_ions = 20
#noise_mult = 4  # peak intensity must be at least this multiple of noise level to be called
#Dw = 2.5  # Within state (local) alignment rt modlation [s] (i think this is tolerance of RT shift?)
#Gw = 0.45  # Within state local alignment gap penalty. Lower G is preferable as higher G favours peak-mixing
#Db = 0.50  # Between state (global) rt modulation [s]
#Gb = 0.40  # Between state (global) gap penalty


outdir = ""
exprdir = ""

# this is the peak detection only pipeline
def detect(args):

    # define path to data files
    global outdir
    outdir = os.path.join(args.indir, 'easyGC_out')
    print outdir

    global exprdir
    exprdir = os.path.join(outdir,"expr")
    print exprdir

    if os.path.isdir(outdir) == False:
        os.mkdir(outdir)

    os.chdir(args.indir)
    runlist = []
    for file in glob.glob("*."+args.ftype):
        print file
        runlist.append(file)

    print runlist
    exprlist = detect_peaks(runlist, args)

# loop over all runs and store the peaks as 'experiments'
# within replicates alignment parameters
def detect_peaks(runs, args):
    expr_list = []
    pl_list = []
    if os.path.isdir(exprdir) == False:
        os.mkdir(exprdir)

    from sys import platform as _platform
    if _platform == "win" or _platform == "win32":
        for run in runs:
            try:
                pl, rid = detect_one_run(run, args)
                expr = store_as_expr(rid, pl, args)
                expr_list.append(expr)
            except:
                print "run failed: ", run
                traceback.print_exc()

    if _platform == "linux" or _platform == "linux2":
        pool = mp.Pool(processes=args.threads)
        results = [pool.apply_async(detect_one_run, args=(r,args,)) for r in runs]
        pl_list = [p.get() for p in results]

        for pl in pl_list:
            if len(pl[0]) > 0:
                expr = store_as_expr(pl[1], pl[0], args)
                expr_list.append(expr)

        try:
            pool.terminate()
        except Exception, e:
            print str(e)
            sys.exit()

    return expr_list


def detect_one_run(run, args):
    infile = os.path.join(args.indir, run)
    print "processing GC-MS file:", infile

   # sys.stdout("processing GCSM run:", run)
    # load the input GC-MS file
    try:
        if args.ftype == 'CDF':
            from pyms.GCMS.IO.ANDI.Function import ANDI_reader
            data = ANDI_reader(infile)
        elif args.ftype == 'JDX':
            #data = JCAMP_reader(in_file)
            data = pyms.GCMS.IO.JCAMP.Function.JCAMP_OpenChrom_reader(infile)
        else:
            raise ValueError('can only load ANDI (CDF) or JDX files!')
    except:
        print "Failure to load input file ", infile
    else:
        data.trim(args.trimstart+"m",args.trimend+"m")
        # get TIC. Would prefer to get from smoothed IM but API is faulty!
        tic = data.get_tic()
        # integer mass
        im = build_intensity_matrix_i(data)

        # would be nice to do noise_mult*noise_level using the noise level AFTER smoothing,
        # but i can't seem to get the TIC for the smoothed IM.
        peak_list = call_peaks(im, tic, True, args)
        return peak_list, run


# this is the alignment-only pipeline
def align(args):
    global exprdir
    exprdir = args.exprdir
    print exprdir

    exprlist = load_expr_list()
    multi_align_local(exprlist, args.distance, args.gap, args.mincommon, tofile=True, transposed=args.transposed)




# this is the full pipeline
def detect_and_align(args, chunked=False, numchunks=1):

    # define path to data files
    global outdir
    outdir = os.path.join(args.indir, 'easyGC_out')
    print outdir

    global exprdir
    exprdir = os.path.join(outdir,"expr")
    print exprdir

    if os.path.isdir(outdir) == False:
        os.mkdir(outdir)

    os.chdir(args.indir)
    runlist = []
    for file in glob.glob("*."+args.ftype):
        print file
        runlist.append(file)

    print runlist
    exprlist = detect_peaks(runlist, args)

    multi_align_local(exprlist, args.distance, args.gap, args.mincommon, tofile=True)



"""
Some useful Savitzky-Golay info: The width of the windows and the order of the polynomial can be changed so it
fits the noise level and complexity of the raw data, e.g., high noise levels require a wider
window (default=7), while high complexity in data requires a higher order polynomial (default=2).

The Durbin-Watson Classifier in OpenChrom can be used to determine the best width and order to use. A DW of 2.0 is optimal
"""
def call_peaks(im, tic, smooth, args):
    print "calling peaks"
    if smooth:
        print "Smoothing IM first..."
        im.crop_mass(args.lowmass, args.highmass)
        print "cropped masses..."
        # get the size of the intensity matrix
        n_scan, n_mz = im.get_size()
        print "# masses in intensity matrix: ", n_mz
        # smooth data
        for ii in range(n_mz):
            ic = im.get_ic_at_index(ii)
            #print "got ic for mass ", ii
            # ic1 = savitzky_golay(ic)
            ic_smooth = savitzky_golay(ic, window=args.window, degree=2)
            #print "savitky golay ran "
            ic_base = tophat(ic_smooth, struct="1.5m")
            #print "tophat ran "
            im.set_ic_at_index(ii, ic_base)
            #print "smoothed mass ", ii
        print "smoothed IM..."
        # noise level calc
        tic1 = savitzky_golay(tic)
        tic2 = tophat(tic1, struct="1.5m")
        noise_level = window_analyzer(tic2)
        print "Noise level in TIC: ", noise_level


    # get the list of Peak objects using BB peak detection / deconv
    pl = BillerBiemann(im, args.window, args.scans)
    print "Initial number of Peaks found:", len(pl)


    # filter down the peaks.
    #   - First: remove any masses from each peak that have intensity less than r percent of the max intensity in that peak
    #   - Second: remove any peak where there are less than n ions with intensity above the cutoff
    pl2 = rel_threshold(pl, percent=args.minintensity)
    pl3 = num_ions_threshold(pl2, n=args.minions, cutoff=noise_level * args.noisemult)
    print "Peaks remaining after filtering:", len(pl3)

    for peak in pl3:
        # peak.null_mass(73)
        peak.null_mass(207)     # column bleed
        peak.null_mass(84)      # solvent tailing

        area = peak_sum_area(im, peak)  # get the TIC area for this peak
        peak.set_area(area)
        area_dict = peak_top_ion_areas(im, peak, args.topions)  # get top n ion areas for this peak
        peak.set_ion_areas(area_dict)

    return pl3

# creates an Experiment from a GCMS run and  writes it to a .expr file
def store_as_expr(run, peak_list, args):
    # create an experiment
    print "creating expression file for run ",run
    expr = Experiment(run, peak_list)
    print "created expression file "
    # set time range for all experiments
    expr.sele_rt_range([args.trimstart+"m", args.trimend+"m"])

    filename, fileext = os.path.splitext(run)
    print "storing expression file at", os.path.join(exprdir,(filename+".expr"))
    store_expr( os.path.join(exprdir,(filename+".expr")), expr )

    return expr


def load_expr_list_from_runlist(runs):
    el = []
    for run in runs:
        file_name = os.path.join(expr_dir, run + ".expr")
        expr = load_expr(file_name)
        el.append(expr)

    return el


def load_expr_list():
    # loads expr list from a directory of exprs
    global exprdir
    el = []

    for file in glob.glob( os.path.join(exprdir, '*.expr') ):
        print file
        expr = load_expr(file)
        el.append(expr)

    return el




def load_run(infile):

    try:
        if args.ftype == 'CDF':
            from pyms.GCMS.IO.ANDI.Function import ANDI_reader
            data = ANDI_reader(infile)
        elif args.ftype == 'JDX':
            #data = JCAMP_reader(in_file)
            data = pyms.GCMS.IO.JCAMP.Function.JCAMP_OpenChrom_reader(infile)
        else:
            raise ValueError('can only load ANDI (CDF) or JDX files!')
    except:
        print "Failure to load input file ", infile
    else:
        data.trim("4.0m", "20.0m")
        # get TIC. Would prefer to get from smoothed IM but API is faulty!
        tic = data.get_tic()
        # integer mass
        return build_intensity_matrix_i(data), tic
        # return build_intensity_matrix(data), tic

# takes a list of Experiment objects and does a pairwise alignment. Optionally writes RTs, Areas and Ions to file.
def multi_align_local(exprlist, Dw, Gw, min_common=1, tofile=True, transposed=True):
    global exprdir
    global outdir
    out_prefix = (datetime.datetime.now()).strftime("%Y-%m-%d-%H%M")

    print "locally aligning peaks from ", exprdir
    F1 = exprl2alignment(exprlist)
    T1 = PairwiseAlignment(F1, Dw, Gw)
    A1 = align_with_tree(T1, min_peaks=min_common)

    if tofile == True:
        ci_list = A1.common_ion()
        if transposed == True:
            A1.write_transposed_output(os.path.join(exprdir, out_prefix + "_aligned_rt.xlsx"))
        else:
            A1.write_excel(os.path.join(exprdir, out_prefix + "_aligned_rt.xlsx"))

        A1.write_csv_dk( os.path.join(exprdir, out_prefix + '_aligned_rt.csv'), os.path.join(exprdir, out_prefix + '_aligned_area.csv') )
        #A1.write_common_ion_csv( os.path.join(exprdir, out_prefix + '_area_common_ion.csv'), ci_list)
        #A1.write_ion_areas_csv( os.path.join(exprdir, out_prefix +  '_aligned_ions.csv') )
    return A1


# takes a list of local alignments and aligns them globally
def multi_align_global(alignments_list, Db, Gb, min_common=1, tofile=True):
    print "globally aligning local alignments from list: ", alignments_list
    T1 = PairwiseAlignment(alignments_list, Db, Gb)
    A1 = align_with_tree(T1, min_peaks=min_common)

    if tofile == True:
        A1.write_csv(expr_dir + output_prefix + 'aligned_rt.csv', expr_dir + output_prefix + 'aligned_area.csv')
        # A1.write_common_ion_csv(expr_dir+output_prefix+'area_common_ion.csv', common_ion_list)
        A1.write_ion_areas_csv(expr_dir + output_prefix + 'aligned_ions.csv')
    return A1


def identify_peak(expr, peaknum):
    pl = expr.get_peak_list()
    ms = pl[peaknum].get_mass_spectrum()
    ms_search_massbank(ms)


#        for peak in pl:
#            #ions = peak.get_ion_areas()
#            #print peak.get_rt()/60, ions
#            ms = peak.get_mass_spectrum()
#            print peak.get_rt()/60, output_golm_ms(ms)

def write_peak_ions_csv(expr):
    pl = expr.get_peak_list()

    for peak in pl:
        rt = peak.get_rt()
        # ions = peak.get_ion_areas()
        ions = peak.peak_top_ion_areas()
        print rt / 60, ions


def list_peaks(expr, peak_list=None):
    if expr != None:
        peak_list = expr.get_peak_list()
    for peak in peak_list:
        print peak.get_rt() / 60, peak.get_area()


def peak_area_range(expr_list):
    peakareas = list()
    for expr in expr_list:
        for peak in expr.get_peak_list():
            peakareas.append(peak.get_area())

    peakareas.sort()
    print peakareas
    return peakareas


from itertools import chain, izip


def output_golm_ms(ms):
    a1 = ms.mass_list
    a2 = (ms.mass_spec).tolist()
    a3 = list(chain.from_iterable(izip(a1, a2)))
    return ' '.join(str(e) for e in a3)


# returns a list of filenames prefix<START> to prefix<STOP>
import random
def generate_runlist(prefix, first, last, randomize=False):
    l = []
    for i in range(first, last + 1):
        l.append(prefix + '%03i' % (i,))

    if randomize == True:
        random.shuffle(l)
    return l


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i + n]





#import osa
#import numpy
#
#
#def ms_search_massbank(ms):
#    cl = osa.Client('http://www.massbank.jp/api/services/MassBankAPI?wsdl')
#    masses = array(ms.mass_list)
#    ints = array(ms.mass_spec)
#    ints = ints / max(ints)
#    # now take the top 20
#    mzs = map(str, ms.mass_list)
#    inte = map(str, ms.mass_spec)
#    unit = ''
#    tol = ''
#    cutoff = ''
#    inst = ['all']
#    ion = 'Positive'
#    res = cl.service.searchSpectrum(mzs, inte, unit, tol, cutoff, inst, ion, 10)
#    print res


#from pyms.Gapfill.Function import *
#
#
#def fillgaps():
#    m = file2matrix(expr_dir + output_prefix + "area_common_ion.csv")
#    sample_list = mp_finder(m)
#    for s in sample_list:
#        missing_peak_finder(s, os.path.join(base_path, "CDF/" + s.get_name() + ".CDF"))
#    write_filled_csv(sample_list, expr_dir + output_prefix + "area_common_ion.csv",
#                     expr_dir + output_prefix + "area_gapfilled.csv")
