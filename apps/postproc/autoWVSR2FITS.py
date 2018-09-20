"""
Find and translate WVSR FFT files to SDFITS
"""
import glob
import logging
import shlex
import time

from multiprocessing import cpu_count
from os.path import basename, exists
from subprocess32 import PIPE, Popen

from Automation import activity_project, AutoDirOperator
from local_dirs import fits_dir, wvsr_fft_dir
from support.logs import set_module_loggers, initiate_option_parser
from support.logs import init_logging, get_loglevel

logger = logging.getLogger(__name__)

class WVSR_to_FITS(object):
  """
  Class to convert SAO .h5 files to FITS.
  
  It can create a file (default: stdout) which can be executed as a
  shell script so user can check intended actions.
  """
  def __init__(self, projects=None):
    """
    initialize .h5 finder
    """
    mylogger = logging.getLogger(logger.name+"WVSR_to_FITS")
    if projects:
      self.projects = projects
    else:
      mylogger.error("__init__: no projects specified")
      raise RuntimeError("a list of one or more projects is required")
    self.logger = mylogger
    self.basepath = wvsr_fft_dir
  
  def prepare_jobs(self):
    """
    """
    fftdirs = glob.glob(self.basepath+"*[cgm]*")
    max_processes = cpu_count()-8
    if len(fftdirs) > max_processes:
      fftdirs = fftdirs[:max_processes]
      self.logger.warning("prepare_jobs: will only do %d directories",
                          max_processes)
    self.commands = []
    for fftdir in fftdirs:
      dir_parts = fftdir.split('/')
      session_code = dir_parts[-1]
      year = 2000+int(session_code[:2])
      cmplx_code = session_code[2]
      doy = int(session_code[3:])
      fft_files = glob.glob(fftdir+"/*.fft")
      if fft_files:
        # we need to get the station and the activity
        name_parts = basename(fft_files[0]).split('_')
        dss = int(name_parts[1][1:])
        activity = name_parts[2]
        if activity == "EGG":
          activity = "EGG0"
        project = activity_project(activity)
        if project in self.projects:
          session_postfix = "dss%2d/%4d/%03d/" % (dss, year, doy)
          self.logger.debug("prepare_jobs: doing %s for %s (%s)", session_postfix, activity, project)
        else:
          self.logger.debug("prepare_jobs: skipping %s", basename(fftdir))
          continue
      else:
        self.logger.warning("prepare_jobs: %s has no FFT files", basename(fftdir))
        continue
      # was this already done?
      fits_files = glob.glob(fits_dir+session_postfix)
      if fits_files:
        self.logger.warning("prepare_jobs: already converted: %s", fits_files)
        continue
      else:
        command = "/usr/bin/python WVSR2SDFITS.py"
        command += " --activity=%s --date=%4d/%03d --dss=%d" % \
                   (activity, year, doy, dss)
        self.commands.append(command)

  def launch_jobs(self):
    """
    """
    self.process = {}
    for index in range(len(self.commands)):
      try:
        self.process[index] = Popen(shlex.split(self.commands[index]),
                                    stdout=PIPE, stderr=PIPE)
      except Exception, details:
        self.logger.error("launch_jobs: command %d failed; %s", index, details)
      
  
  def check_outputs(self):
    self.output = {}
    self.error = {}
    while self.process.keys():
      for index in self.process.keys():
        self.logger.debug("check_outputs: for %d", index)
        if self.process[index].poll():
          try:
            self.output[index], self.error[index] = \
                                              self.process[index].communicate()
            self.logger.info("check_outputs: %d finished", index)
          except ValueError, details:
            self.logger.error("check_outputs: %s", str(details))
          self.process.pop(index)
      time.sleep(1)
    return self.output, self.error

if __name__ == "__main__":
  examples = """
  autoWVSR2FITS.py: 
    loop through all YYdDDD sub-dirs of /data/post_processing/auto
  """
  
  p = initiate_option_parser("loops through HDF5 sub-directories",examples)
  p.usage = "autoWVSR2FITS.py [options]"
  # the above defines options -h, --help, --file_loglevel, -l, --logfilepath,
  # and --module_loglevels
  p.add_argument('-d', '--date',
                 dest = 'date',
                 type = str,
                 default = None,
                 help = 'last date to process (default: yesterday)')
  p.add_argument('-n', '--no_dry_run',
                 dest = 'dry_run',
                 action = 'store_false',
                 default = True,
                 help = 'last date to process (default: yesterday)')
  p.add_argument('-p', '--projects',
                 dest = 'projects_str',
                 type = str,
                 default = 'FREGGS,ISM_RRL',
                 help = 'comma-separated list of projects')
  
  args = p.parse_args()
  logging.basicConfig(level=logging.INFO)
  mylogger = init_logging(logging.getLogger(),
                          loglevel   = get_loglevel(args.file_loglevel),
                          consolevel = get_loglevel(args.console_loglevel),
                          logname    = args.logpath+"autoWVSR2FITS.log")
  mylogger.debug(" Handlers: %s", mylogger.handlers)
  loggers = set_module_loggers(eval(args.modloglevels))
  do = WVSR_to_FITS(projects=args.projects_str.split(','))
  do.prepare_jobs()
  do.launch_jobs()
  print do.check_outputs()
  


