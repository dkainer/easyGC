import sys, os, argparse, glob
from GCMSalign import *

parser = argparse.ArgumentParser(prog='easyGC', description="easyGC calls peaks, deconvolutes them, quantitates them and then aligns them across all your GC-MS samples...magically!", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(title='use one of these subcommands')

peak_parser = subparsers.add_parser('peakcall', help='peakcall help', description="run the peak caller on a directory of GC-MS files", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
peak_parser.set_defaults(which='peakcall')
peak_parser.add_argument('-i', '--indir', required=True, help='[REQUIRED] absolute path to the directory containing your raw GC-MS files', type=str)
peak_parser.add_argument('-f', '--ftype', required=True, choices=['CDF','JDX'], help='[REQUIRED] This is the type of input GC-MS files you have. CDF is not supported on Windows', type=str)
peak_parser.add_argument('-TS','--trimstart', required=True, help='[REQUIRED] time in minutes (X.XX) in the chromatogram from where the analysis should begin. Helps to cut out junk at the start', type=str)
peak_parser.add_argument('-TE','--trimend', required=True, help='[REQUIRED] time in minutes (X.XX) in the chromatogram where the analysis should end. Helps to cut out junk at the end', type=str)
peak_parser.add_argument('-LM','--lowmass', required=False, default=50, help='The lowest mass found across your samples. No masses below this will be used in analysis', type=int)
peak_parser.add_argument('-HM','--highmass', required=False, default=250, help='The highest mass found across your samples. No masses above this will be used in analysis', type=int)
peak_parser.add_argument('-W', '--window', required=False, default=9, help='peak calling: width (in scans) of window over which local ion maxima are detected. Should be similar to the width off your peaks.', type=int)
peak_parser.add_argument('-S', '--scans', required=False, default=3, help='peak calling: distance (in scans) at which locally apexing ions can be combined into one peak', type=int)
peak_parser.add_argument('-N', '--minions', required=False, default=3, help='peak calling: min number of apexing ions with intensity above a threshold required for a peak to be called. Higher = less peaks called', type=int)
peak_parser.add_argument('-R', '--minintensity', required=False, default=5, help='peak calling: min intensity (percent) of an ion relative to max peak intensity for that ion to be included in the peak', type=int)
peak_parser.add_argument('-M', '--noisemult', required=False, default=3.0, help='peak calling: total peak intensity must be at least this multiple of the base noise level to be called. Higher multiple means fewer peaks called', type=float)
peak_parser.add_argument('-I', '--topions', required=False, default=5, help='from the list of most important ions in a peak, how many should be outputted as a mini mass-spec?', type=int)
peak_parser.add_argument('-T', '--threads', required=False, default=1, help='Peak calling is linearly sped up by using more threads. Currently only multithreaded on linux!', type=int)

align_parser = subparsers.add_parser('align', help='align help', description="run the peak aligner on a directory of .expr files", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
align_parser.set_defaults(which='align')
align_parser.add_argument('-e', '--exprdir', required=True, help='[REQUIRED] absolute path to the directory containing the .expr files from a previous peak calling run.', type=str)
align_parser.add_argument('-D', '--distance', required=False, default=2.5, help='local alignment: distance in retention time (seconds) over which the local peak aligner should search for similar peaks to this one', type=float)
align_parser.add_argument('-G', '--gap', required=False, default=0.40, help='local alignment: gap penalty. Lower G results in more peaks in the output. Higher G result in fewer output peaks but possibly some peaks contain multiple merged peaks', type=float)
align_parser.add_argument('-C', '--mincommon', required=False, default=1, help='local alignment: minimum number of samples that an aligned peak must be called in for it to be outputted', type=int)
align_parser.add_argument('-TR','--transposed', required=False, default=True, help='the output matrix will show compounds as columns (True) or as rows (False)', type=bool)

args = parser.parse_args()
print args

if args.which == 'peakcall':
    print "running peak detector and quantification"
    detect(args)

elif args.which == 'align':
    print "running peak aligner"
    align(args)



else:
    print "invalid command"
