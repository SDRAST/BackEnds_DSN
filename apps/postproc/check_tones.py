"""
check_tones - see whether the tones are in the right place
"""
#import IPython
#IPython.version_info = IPython.release.version.split('.')
#IPython.version_info.append('')

import dateutil
import datetime
import glob
import logging
import math
import numpy
import sys
import time
import warnings

from os import chdir, getcwd, symlink, makedirs
from os.path import basename, exists
from pylab import *

from Data_Reduction.DSN.WVSR.SpecData import read_FFT_file
from local_dirs import wvsr_fft_dir
from support.lists import unique
from support.logs import initiate_option_parser, init_logging
from support.logs import get_loglevel, set_loglevel

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
warnings.filterwarnings('error')

def make_FFTdir_name(dss, year, doy):
  """
  """
  src = "%02d" % (year - 2000)
  if dss == 14:
    src += "g"
  elif dss == 43:
    src += "c"
  elif dss == 63:
    src += "m"
  else:
    raise Exception("%s is not a valid 70-m antenna", dss)
  src += "%03d/" % doy
  return src

    
if __name__ == "__main__":
  examples = """
Examples
========
  run WVSR2SDFITS.py --console_loglevel=debug --DSS=14 \
                          --project=AUTO_EGG --date=2016/237
"""  
  p = initiate_option_parser(__doc__, examples)
  p.usage='WVSR2SDFITS.py [options]'
  p.add_argument('-a', '--activity',
               dest = 'activity',
               type = str,
               default = None,
               help = "Activity code")
  p.add_argument('--date',
               dest = 'date',
               type = str,
               default = "2016/237",
               help = 'Date of observation as YEAR/DOY string')
  p.add_argument('-D', '--dss',
               dest = 'dss',
               type = int,
               default = 14,
               help = 'DSN station number')
  p.add_argument('-s', '--scans',
               dest = 'scans',
               type = str,
               default = "",
               help = 'comma separated scan numbers')

  args = p.parse_args(sys.argv[1:])
  
  mylogger = logging.getLogger()
  mylogger.setLevel(logging.DEBUG)
  init_logging(mylogger,
                 loglevel = get_loglevel(args.file_loglevel),
                 consolevel = get_loglevel(args.console_loglevel),
                 logname = args.logpath+"WVSR2SDFITS.log")
                 
  mylogger.debug("WVSR2SDITS args: %s", args)
  yearstr, doystr = args.date.split("/")
  year = int(yearstr)
  doy = int(doystr)
  fftdir = wvsr_fft_dir+make_FFTdir_name(args.dss, year, doy)

  IFs = ["IF1","IF2"]
  
  if args.scans == "":
    allfiles = glob.glob(fftdir+yearstr[2:]+"-"+doystr+"*s*_d"+str(args.dss)+"*.fft")
    scans = []
    for fname in allfiles:
      scan = int(fname.split('-')[-1].split('_')[0][1:])
      scans.append(scan)
    scans = unique(scans)
  else:
    scans = args.scans.split(',')
    scans.sort()
  for scan in scans:
    scanstr = "%04d" % int(scan)
    scanfiles = glob.glob(fftdir+yearstr[2:]+"-"+doystr+"*s"+scanstr+"_d"+str(args.dss)+"*.fft")
    scanfiles.sort()
    fig, ax = subplots(2,2, figsize=(16,16), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    extra_peak = False
    for scanfile in scanfiles:
      if data:
        pass
      else:
        mylogger.warning("no data for scan %d", scan)
        continue
      data = read_FFT_file(scanfile)
      rowidx = scanfiles.index(scanfile)
      mylogger.debug("processing scan %s file %s for row %d",
                     scan, basename(scanfile), rowidx)
      for IF in IFs:
        IFidx = IFs.index(IF)
        indices = data[IF+'-ps'].argsort()[::-1]
        result = ax[rowidx,IFidx].hist(data['freq'][indices[:500]]/1e6 % 1, bins=100)
        spikes = where(result[0] > 30)[0]
        if len(spikes) > 1:
          mylogger.info("Scan %s has histogram peaks at %s in row %d",
                        scan, result[1][spikes], rowidx)
          extra_peak = True
          specfig = figure()
          specax = specfig.add_subplot(1,1,1)
          specax.semilogy(data['freq']/1e6, data[IF+'-ps'])
          specax.set_title("Scan "+str(scan)+" subch "+str(rowidx+1)+" "+IF)
          specax.set_xlabel("Frequency offset (MHz)")
          specfig.savefig("sc"+str(scan)+"ch"+str(rowidx+1)+IF+".png")
          specfig.show()
        ax[rowidx,IFidx].grid(True)
        ax[rowidx,IFidx].set_title(IF)
        if IFidx == 0:
          ax[rowidx,IFidx].set_ylabel("subch "+str(rowidx+1))
        if rowidx == 1:
          ax[rowidx,IFidx].set_xlabel("Relative offset (MHz)")
    fig.suptitle("Scan "+str(scan))
    if extra_peak and args.scans == "":
      # processing all scans to find extra peaks
      fig.savefig("hist-scan"+str(scan)+".png")
      fig.hold(True)
      fig.show()
    elif args.scans:
      # plotting selected scans
      fig.savefig("hist-scan"+str(scan)+".png")
      fig.show()
    else:
      # no peak found while processing all scans
      fig.hold(False)
      mylogger.debug("Plot for scan %s can be discarded", scan)


