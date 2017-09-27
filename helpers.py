"""
"""
import glob
import logging
import re

from os.path import basename, splitext

from Automation import get_session_dirs
from Automation.sources import get_all_source_data
from Automation.WVSR import parse_scan_files
from Data_Reduction.DSN.WVSR.SpecData import get_channel_IDs
from DatesTimes import WVSR_script_time_to_timestamp
from MonitorControl.Configurations.DSN_standard import standard_equipment
from support.lists import unique


logger = logging.getLogger(__name__)

class WVSRmetadataCollector:
  """
  Class to assemble data for a session data file
  
  Public attributes::
    chan_names   - list of channel IDs
    datadir      - session directory in 'project_data' directory
    date         - datetime.datetime for observation date
    doy          - day of year of observation
    dss          - station where observations were done
    equip        - telescope equipment description
    logger       - logging.Logger object
    project      - name of AUTO project (not real name)
    real_obs_dir - project's path to session directory
    scaninfo     - metadata for every scan
    scankeys     - sorted list of scan numbers
    sourcedata   - details of sources and calibrators
    wvsrdir      - location of WVSR files
    wvsrlogs     - list of log file names
    year         - year of observation
  """
  def __init__(self, project, dss, year, doy):
    """
    initiate a WVSR configuration description collector
    
    Mainly this creates two dicts::
      wvsr_cfg - configuration for one WVSR
      fft_meta - details on the processing all the scans
    """
    self.project = project
    self.dss = dss
    self.year = year
    self.doy = doy
    self.logger = logging.getLogger(logger.name+".WVSRmetadataCollector")
    
    # get all the directories involved
    auto_obsdir, self.real_obsdir, self.project_dir, self.datadir, \
        self.wvsrdir, self.fftdir = get_session_dirs( project, dss, year, doy)
    # get the high-level metadata
    get_all_source_data(self.project_dir)
    self.get_scan_info()
    self.get_WVSR_names()
    self.get_metadata_from_WVSR_logs()
    self.equip = {}
    # specify station equipment for each IF
    for wvsr in self.wvsrnames:
      self.equip[wvsr] = {}
      self.logger.debug('__init__: %s config: %s', wvsr, self.wvsr_cfg[wvsr])
      for IF in self.wvsr_cfg[wvsr]['channels']:
        self.logger.debug('__init__: %s IF %s config: %s', wvsr, IF, self.wvsr_cfg[wvsr][IF])
        band = self.wvsr_cfg[wvsr][IF]['IF_source'].split('_')[1]
        self.equip[wvsr][IF] = standard_equipment(dss, band)
        #self.equip[wvsr][IF] = station_configuration(self.equip, project, 
        #                                             dss, year, doy,'X')
    # get the parameters for each FFT configuration
    self.get_metadata_from_FFT_logs()
  
  def get_scan_info(self):
    """
    get the metadata for channels and scans recorded
    
    Example::
      In [12]: collector.scaninfo[6]
      Out[12]: 
      {'chan_id 1': '16-237-001-s0006_d14_RCP_LCP.wvsr.fft',
       'chan_id 2': '16-237-002-s0006_d14_RCP_LCP.wvsr.fft',
       'end': datetime.datetime(2016, 8, 24, 9, 42),
       'source': 'w5w-fregg51-ref',
       'start': datetime.datetime(2016, 8, 24, 9, 40, 30)}
    """
    self.scaninfo = parse_scan_files(self.real_obsdir)
    self.scankeys = self.scaninfo.keys()
    self.scankeys.sort()
    # get IF channel IDs
    self.chan_names = get_channel_IDs(self.scaninfo)
    self.logger.debug("get_scan_info: channel keys: %s", self.chan_names)
  
  def get_WVSR_names(self):
    """
    find the recorders used
    
    creates a dict of WVSR identifiers based on simple name like 'wvsr2'::
      In [17]: collector.wvsrnames
      Out[17]: {'wvsr2': 'WVSR2_10.16-237-084500'}
    """
    self.wvsrlogs = glob.glob(self.wvsrdir+"/WVSR*")
    self.wvsrnames = {}
    for f in self.wvsrlogs:
      logname = basename(f)
      self.wvsrnames[logname[:5].lower()] = logname
    self.logger.debug("get_WVSR_names: %s", self.wvsrnames)
  
  def get_metadata_from_WVSR_logs(self):
    """
    parse the WVSR FFT logs for metadata
    
    Each WVSR typically has two digital IF channels with orthogonal pols.
    'wvsr_cfg' is a dict.  Example::
      In [3]: collector.wvsr_cfg
      Out[3]: 
      {'wvsr2': {1: {'DSP': 'dsp2',
                     'IF_source': '14_X_RCP',
                     'chan_id 1': {'bandwidth': 8000000,
                                   'bits': 8,
                                   'sfro': 209386000},
                     'chan_id 2': {'bandwidth': 8000000,
                                   'bits': 8,
                                   'sfro': 484824000},
                     'expid': '16_237_DSS-14_ARCP',
                     'pol': 'RCP',
                     'rf_to_if_lo': 8100,
                     'subchannels': ['chan_id 1', 'chan_id 2']},
                 2: {'DSP': 'dsp1',
                     'IF_source': '14_X_LCP',
                     'chan_id 1': {'bandwidth': 8000000, 'bits': 8, 'sfro': 209386000},
                     'chan_id 2': {'bandwidth': 8000000, 'bits': 8, 'sfro': 484824000},
                     'expid': '16_237_DSS-14_ALCP',
                     'pol': 'LCP',
                     'rf_to_if_lo': 8100,
                     'subchannels': ['chan_id 1', 'chan_id 2']},
                 'channels': [1, 2],
                 'from log': 'WVSR2_10.16-237-084500',
                 'user': 'cjn,'}}
    """
    self.wvsr_cfg = {}
    for logname in self.wvsrlogs:
      wvsr = basename(logname)[:5].lower()
      self.logger.debug("get_metadata_from_WVSR: WVSR ID: %s", wvsr)
      self.parse_WVSR_log(logname)
      self.logger.debug("get_metadata_from_WVSR: from log %s",
                                               self.wvsr_cfg[wvsr]['from log'])
      self.logger.debug("get_metadata_from_WVSR: user: %s",
                                                   self.wvsr_cfg[wvsr]['user'])
      self.wvsr_cfg[wvsr]["channels"] = []
      for key in self.wvsr_cfg[wvsr].keys():
        if type(key) == int:
          self.wvsr_cfg[wvsr]["channels"].append(key)
          self.wvsr_cfg[wvsr][key]['subchannels'] = []
          for subkey in self.wvsr_cfg[wvsr][key].keys():
            if subkey[:7] == 'chan_id':
              self.wvsr_cfg[wvsr][key]['subchannels'].append(subkey)
          self.wvsr_cfg[wvsr][key]['subchannels'].sort()
      self.logger.debug("get_metadata_from_WVSR: %s", self.wvsr_cfg)

  def parse_WVSR_log(self, logname):
    """
    Extracts metadata from the WVSR log
    
    Most metadata are not consider time sensitive so that the last value found
    is the value that is returned.  The exception is attenuator data which may
    change during a track.
    
    In quite a few cases we simply exec() a statement this is in the line or
    is constructed from parts of the line.  A typical line looks like this:
    16/237 08:45:15 wvsr2 ATT[2]: att = 15, des_amp = -10, cur_amp = -10.07, max_amp = 0, min_amp = -50
    """
    logfile = open(logname,'r')
    logtext = logfile.read()
    logfile.close()
    lines = logtext.split('\n')
    # WVSR and user from the first line
    line = lines[0].split()
    wvsrID, user = line[2], line[-3]
    self.logger.debug("parse_WVSR_log: for %s", wvsrID)
    EXPID = {}
    RF_TO_IF_LO = {}
    IFS = {}
    CHAN = {}
    ATT = {}
    # parse the rest of log file
    for line in lines[1:]:
      if re.search("EXPID\\[[12]\\]", line): # experiment ID
        # parse a log file line like::
        #   16/237 08:45:13 wvsr2 EXPID[1]: EVT 301 PROGRESS: \
        #              exp 16_237_DSS-14_ARCP created by Client 3625 cjn 4-2053
        # take the variable name from the 4th item leaving off the colon
        # and add item to EXPID dict
        parts = line.split()
        self.logger.debug("parse_WVSR_FFT_logs: EXPID line parts: %s", parts)
        exec(parts[3][:-1]+"= "+"'"+str(parts[8])+"'")
      if re.search('RF_TO_IF_LO\[',line): # receiver or first LO
        # parse a log file line like::
        #   16/237 08:45:15 wvsr2 RF_TO_IF_LO[2]: RF_TO_IF_LO: value = 8100
        # add to RF_TO_IF_LO for the IF channel to dict
        parts = line.split()
        exec(parts[3][:-1]+"= "+str(parts[-1]))
      if re.search("CHAN\\[.\\]: dsp", line): # digital signal processor module
        # parse a lof file line like::
        #   16/237 08:45:17 wvsr2 \
        #          CHAN[2]: dsp1:1 chan_id = 001, bandwidth = 8000000, bits = 8
        # add subchannel and its parameters to CHAN dict
        parts = line.split()
        chan_str = parts[3][:-1]
        subchan_str = "'"+parts[5]+" ' + str(int("+parts[7][:-1]+"))"
        self.logger.debug("parse_WVSR_FFT_logs: for CHAN doing %s with %s",
                     chan_str, subchan_str)
        # create subchannel item in CHAN if needed
        if CHAN == {}:
          exec(chan_str+"={}") # creates a dict for each CHAN
        try:
          eval(chan_str+".has_key("+subchan_str+")")
        except:
          exec(chan_str+"={}")
          exec(chan_str+"["+subchan_str+"]={}")
        else:
          exec(chan_str+"["+subchan_str+"]={}")
        # add assigned DSP
        exec(chan_str+"['DSP']='"+parts[4].split(":")[0]+"'")
        # add DDCLO, add bandwidth and bits/sample
        exec(chan_str+"["+subchan_str+"]['"+parts[8]+"']="
             +str(int(parts[10][:-1])))
        exec(chan_str+"["+subchan_str+"]['"+parts[11]+"']="
             +str(int(parts[13])))
      if re.search("rsp_ddclo", line): # subchannel offset LO
        # create subchannel offset
        parts = line.split()
        self.logger.debug("parse_WVSR_FFT_logs: line parts: %s", parts)
        dsp, subchan = parts[8].split(':')
        offset = int(float(parts[9])) # don't need fractional Hz
        # assign to correct channel
        for key in CHAN.keys():
          if CHAN[key]['DSP'] == dsp:
            if CHAN[key].has_key("chan_id "+subchan):
              CHAN[key]["chan_id "+subchan]["sfro"] = offset
            else:
              self.logger.info(
                  "parse_WVSR_FFT_logs: no subchannel ID %s for IF channel %s",
                     subchan,key)
      if re.search("IFS\\[.\\]:.*PROGRESS", line): # signal source for IF chanl
        # IF switch inputs
        self.logger.debug("parse_WVSR_FFT_logs: IFS from line: %s", line)
        parts = line.split()
        self.logger.debug("parse_WVSR_FFT_logs: line parts: %s", parts)
        # exclude the colon from the signal source
        exec(parts[3][:-1]+"= '"+str(parts[10][:-1])+"'")
      if re.search("ATT\\[.\\]:.*cur_amp", line): # attenuator and power
        self.logger.debug("parse_WVSR_FFT_logs: ATT line parts: %s", parts)
        parts = line.split()
        UT = WVSR_script_time_to_timestamp(*parts[:2])
        try:
          # this executes something like" "ATT[2]['att']=15"
          #self.logger.debug("parse_WVSR_FFT_logs: trying %s",
          #                  parts[3][:-1]+"['"+parts[4]+"']"+"="+parts[6].strip(','))
          #exec(parts[3][:-1]+"['"+parts[4]+"']"+"="+parts[6].strip(','))
          self.logger.debug("parse_WVSR_FFT_logs: trying %s",
                            parts[3][:-1]+"["+str(UT)+"]['"+parts[4]+"']"+"="+parts[6].strip(','))
          exec(parts[3][:-1]+"["+str(UT)+"]['"+parts[4]+"']"+"="+parts[6].strip(','))
        except KeyError,details:
          # ATT[2] is not defines as a dict
          self.logger.debug("parse_WVSR_FFT_logs: failed because %s", details)
          # this executes something like:'ATT[2]={1472028315: {}}'
          self.logger.debug("parse_WVSR_FFT_logs: trying %s",
                            parts[3][:-1]+"={"+str(UT)+": {}}")
          exec(parts[3][:-1]+"={"+str(UT)+": {}}")
          # this executes something like "ATT[2][1472028315]['att']=15"
          self.logger.debug("parse_WVSR_FFT_logs: trying %s",
                            parts[3][:-1]+"["+str(UT)+"]['"+parts[4]+"']"+"="+parts[6].strip(','))
          exec(parts[3][:-1]+"["+str(UT)+"]['"+parts[4]+"']"+"="+parts[6].strip(','))
          # now that the dict exists add the power
          exec(parts[3][:-1]+"["+str(UT)+"]['"+parts[10]+"']"+"="+parts[12].strip(',')) 
    self.logger.info("parse_WVSR_FFT_logs: EXPID = %s", EXPID)
    self.logger.info("parse_WVSR_FFT_logs: RF_TO_IF_LO = %s", RF_TO_IF_LO)
    self.logger.info("parse_WVSR_FFT_logs: IFS = %s", IFS)
    self.logger.info("parse_WVSR_FFT_logs: CHAN = %s", CHAN)
    self.logger.info("parse_WVSR_FFT_logs: ATT = %s", ATT)
    # add to WVSR configuration
    self.wvsr_cfg[wvsrID] = {'user': user,
                             'from log': basename(logname)}
    # now build the dict
    #   CHAN keys will also apply to ATT
    for key in CHAN.keys():
      self.wvsr_cfg[wvsrID][key] = {}
      self.wvsr_cfg[wvsrID][key]['pol'] = EXPID[key][-3:]
      self.wvsr_cfg[wvsrID][key]['expid'] = EXPID[key]
      self.wvsr_cfg[wvsrID][key]['rf_to_if_lo'] = RF_TO_IF_LO[key]
      self.wvsr_cfg[wvsrID][key]['IF_source'] = IFS[key]
      subchans = CHAN[key].keys()
      for subch in subchans:
        self.wvsr_cfg[wvsrID][key][subch] = {}
        if subch[:7] == 'chan_id':
          for param in CHAN[key][subch].keys():
            self.wvsr_cfg[wvsrID][key][subch][param] = CHAN[key][subch][param]
        else:
          self.wvsr_cfg[wvsrID][key][subch] = CHAN[key][subch]
      self.wvsr_cfg[wvsrID][key]['attenuation'] = {}
      for time in ATT[key].keys():
        self.wvsr_cfg[wvsrID][key]['attenuation'][time] = {}
        for param in ATT[key][time].keys():
          self.wvsr_cfg[wvsrID][key]['attenuation'][time][param] = ATT[key][time][param]
    self.logger.debug("parse_WVSR_log: %s results: %s", wvsrID,
                       self.wvsr_cfg[wvsrID])

  def get_metadata_from_FFT_logs(self):
    """"
    Gets metadata about the FFT processing of WVSR files
    
    This must be invoked after 'get_metadata_from_WVSR_logs'.
    
    This assumes that two IFs are combined into one signal with full
    polarization information.  Then the subchannels are the same for both.
    
    The problem with this right now is that there is no clean way to handle
    simultaneous WVSRs so this is going to have to assume that 'wvsrnames' has
    only one item for now.
    
    Returns a dict like::
      In [2]: collector.fft_meta[3]
      Out[2]: 
      {1: {'chan_id 1': {'datafile':  '/data2/16g237/16-237-001-s0003_d14_RCP.wvsr',
                         'first rec': '2016 237:09:34:31.000',
                         'last_rec':  '2016 237:09:35:59.000',
                         'n_recs': 90,
                         'n_secs': 90},
           'chan_id 2': {'datafile':  '/data2/16g237/16-237-002-s0003_d14_RCP.wvsr',
                         'first rec': '2016 237:09:34:31.000',
                         'last_rec':  '2016 237:09:35:59.000',
                         'n_recs': 90,
                         'n_secs': 90}},
       2: {'chan_id 1': {'datafile':  '/data2/16g237/16-237-001-s0003_d14_LCP.wvsr',
                         'first rec': '2016 237:09:34:31.000',
                         'last_rec':  '2016 237:09:35:59.000',
                         'n_recs': 90,
                         'n_secs': 90},
           'chan_id 2': {'datafile':  '/data2/16g237/16-237-002-s0003_d14_LCP.wvsr',
                         'first rec': '2016 237:09:34:31.000',
                         'last_rec':  '2016 237:09:35:59.000',
                         'n_recs': 90,
                         'n_secs': 90},
           'n_freqs':   131072,
           'n_samples': 8000000}}
    """
    fftlogs = glob.glob(self.fftdir+"*.log")
    fftlogs.sort()
    scans = []
    subchls = []
    # get the scan numbers
    for log in fftlogs:
      # parse file name
      logname = splitext(basename(log))[0]
      YR, DOY, scan, subch = self.parse_fft_log_name(logname)
      scans.append(scan)
      subchls.append(subch)
    scannums = unique(scans)
    self.logger.debug("get_metadata_from_FFT_logs: scannums: %s", scannums)
    subchannels = unique(subchls)
    self.logger.debug("get_metadata_from_FFT_logs: subchannels: %s", subchannels)
    # now parse the FFT log files
    self.fft_meta = {}
    for scan in scannums:
      self.fft_meta[scan] = {} # first level dict on scan
      # create a dict for each IF
      for channel in self.wvsr_cfg[self.wvsrnames.keys()[0]]["channels"]:
        self.fft_meta[scan][channel] = {} # second level dict on IF channel
      # each IF has the same number of subchannels
      for subch in subchannels:
        # parse the FFT log
        logname = "%2d-%03d-%03d-s%04d_d%2d_RCP_LCP.wvsr.log" % \
                  (self.year-2000, self.doy, subch, scan, self.dss)
        logfile = open(self.fftdir+logname)
        lines = logfile.readlines()
        logfile.close()
        subch_key = "chan_id "+str(subch)
        self.logger.debug("get_metadata_from_FFT_logs: for %s", subch_key)
        extracted = self.parse_fft_log(lines)
        self.logger.debug("get_metadata_from_FFT_logs: extracted %s",
                          extracted)
        # move to metadata structure
        for ch in extracted.keys(): # these keys are channels
          self.fft_meta[scan][ch][subch_key] = {} # third level dict on subch
          for key in extracted[ch].keys(): # these keys are variable names
            #self.logger.debug("get_metadata_from_FFT_logs: ch %d sub %s key '%s'",
            #                  ch, subch_key, key)
            if key == 'n_samples' or key == 'n_freqs':
              #self.fft_meta[scan][ch][key] = extracted[ch][key]
              self.fft_meta[scan][key] = extracted[ch][key]
            else:
              self.fft_meta[scan][ch][subch_key][key] = extracted[ch][key]
              #self.fft_meta[scan][subch_key][key] = extracted[ch][key]
          
  def parse_fft_log_name(self, logname):
    """
    Gets metadata from FFT log file name
    
    The filename is of the form::
      YY-DOY-sub-sSCAN_dSS_RCP_LCP.wvsr.log
    where 'sub' is the 3-digit sub-channel number, 'SCAN' is the 4-digit scan
    number, 'SS' is the 2-digit station number. 
    """
    name_parts = logname.split('-')
    YR = int(name_parts[0])
    if 2000+YR != self.year:
      self.logger.error("get_metadata_from_FFT_logs: %s is for wrong year",
                        logname)
      return False
    DOY = int(name_parts[1])
    if DOY != self.doy:
      self.logger.error("get_metadata_from_FFT_logs: %s is for wrong DOY",
                        logname)
      return False
    subch = int(name_parts[2])
    lastparts = name_parts[3].split("_")
    scan = int(lastparts[0][1:])
    dss = int(lastparts[1][1:])
    if dss != self.dss:
      self.logger.error("get_metadata_from_FFT_logs: %s is for wrong DSS",
                        logname)
      return False
    return YR, DOY, scan, subch
      
  def parse_fft_log(self, lines):
    """
    Get metadata from one scan's log file.
        
    Each log has a preamble section with details about the processing followed
    by reports for each FFT.
    
    The output is a dict like this example for scan 1::
      In [2]: collector.fft_meta[1]
      Out[2]: 
      {1: {'chan_id 2': {'datafile': '/data2/16g237/16-237-002-s0001_d14_RCP.wvsr',
                         'first rec': '2016 237:09:30:31.000',
                         'last_rec': '2016 237:09:32:00.000',
                         'n_recs': 91,
                         'n_secs': 91}},
       2: {'chan_id 2': {'datafile': '/data2/16g237/16-237-002-s0001_d14_LCP.wvsr',
                         'first rec': '2016 237:09:30:31.000',
                         'last_rec': '2016 237:09:32:00.000',
                         'n_freqs': 131072,
                         'n_recs': 91,
                         'n_samples': 8000000,
                         'n_secs': 91}}}
        where the first key is the IF channel.
    
    We don't know yet how to handle FFT logs if there are two or more WVSRs so
    'wvsr_cfg' for one and only one has a similar structure with WVSR ID in
    place of scan number.
    """
    in_preamble = True
    first_part = True
    found = {}
    for line in lines:
      parts = line.strip().split()
      if parts == []:
        continue
      if in_preamble:
        if line[:10] == "Input file": # this starts a channel section
          datafile = parts[2].lstrip('<').rstrip('>')
          VSRtype = line[-1]
        if parts[0] == "Channel" and len(parts) > 2:
          channel = int(parts[1][:-1]) # strip off colon
          found[channel] = {}
          found[channel]["datafile"] = datafile
        if re.search("First", line):
          found[channel]["first rec"] = ' '.join(parts[2:])
        if re.search("Last", line):
          found[channel]["last_rec"] = ' '.join(parts[2:])
        if line[:4] == "nrec":
          found[channel]["n_recs"] = int(line[6:10])
          found[channel]["n_secs"] = int(line[22:])
        if parts[0] == 'Unequal':
          in_preamble = False
      else:
        if re.search("ns.*npts.*power_two", line):
          found[channel]["n_samples"] = int(parts[2][:-1])
          found[channel]["n_freqs"] = int(parts[5])
          break
    return found
    
    def make_backend(self):
      """
      """
      pass

