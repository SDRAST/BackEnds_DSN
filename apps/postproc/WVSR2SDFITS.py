"""
WVSR2SDFITS - converts WVSR spectrometer FFT files to DSN SDFITS

A DSN SDFITS file has a separate extension for each backend as described in the
station configuration.

Steps::
  FITSfile.__init__()
    make_prihdu()                      # standard SDFITS HDU with required items
    add_site_data()                    # telescope location
  
  FITSfile.make_WVSR_table()
    start extension header             # puts in the required items
    add_site_data()                    # repeats what is in the primary HDU
    make_basic_columns()               # SCAN, OBJECT, DATE-OBS, etc.
    get_hardware_metadata()            # BACKEND, MAXIS1, FREQRES
    add_time_dependent_columns()       # LST, AZIMUTH, ELEVATIO
    add_IF_dependent_columns()         # TSYS
    make_data_axis()                   # DATA
    make_offset_columns()              # BEAMxOFF
    add_data()
"""
import IPython
IPython.version_info = IPython.release.version.split('.')
IPython.version_info.append('')

#from pylab import * #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import astropy.io.fits as pyfits
import astropy.units as u
import dateutil
import glob
import logging
import math
import numpy
import sys
import time
import warnings

from astropy.coordinates import FK5, SkyCoord
from astropy.time import Time
from os import chdir, getcwd, symlink, makedirs
from os.path import basename, exists
from scipy.optimize import curve_fit

from Astronomy import calendar_date
from Astronomy.redshift import V_LSR
from Automation import get_session_dirs, get_start_time
from Automation.NMClogs import NMClogManager,NMC_log_server 
from Automation.sources import get_all_source_data
from Automation.tones import tone_chnl_nums
from Automation.WVSR import get_Naudet_FFT_dir, make_datadir_name
from Data_Reduction import get_obs_dirs, get_obs_session, select_data_files
from Data_Reduction import get_num_chans, reduce_spectrum_channels
from Data_Reduction.DSN.WVSR.SpecData import read_FFT_file
from Data_Reduction.FITS.DSNFITS import FITSfile
from Data_Reduction.tipping import airmass
from DatesTimes import datetime_to_UnixTime
from local_dirs import sci_proj_path
#from Data_Reduction.DSN.SAO import parse_filename
from Math.multigauss import gaussian
from MonitorControl.BackEnds.DSN.helpers import WVSRmetadataCollector
from MonitorControl.Configurations.coordinates import DSS
from MonitorControl.Configurations.DSN_standard import standard_equipment
from MonitorControl.Configurations.GDSCC.WVSR import station_configuration
#from MonitorControl.SDFITS import FITSfile
from support import mkdir_if_needed
from support.logs import initiate_option_parser, init_logging
from support.logs import get_loglevel, set_loglevel

import Math.multigauss as multigauss

obsmode = 'LINEPSSW'
veldef = 'RADI-OBS'
num_tonespec_chls = 16

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
warnings.filterwarnings('error')

nanarray = numpy.array([numpy.nan, numpy.nan]) # used for tone fit pars
blank_array = numpy.array(256*[numpy.nan])

def biased_scaled_gaussian(x, bias, amp, offset, std):
  """
  scales an unscaled gaussian
  """
  return bias + amp*gaussian(x, offset, std)

def biased_scaled_sinc(x, bias, amp, offset, dx):
  """
  """
  return bias + \
              amp*(numpy.sin(math.pi*(x-offset)/dx)/(math.pi*(x-offset)/dx))**2


class FITSfile_from_WVSR(FITSfile):
  """
  subclass for creating FITS file from WVSR data
  
  WVSR spectra may contain tones from a phase cal generator which are typically
  spaced by 1 MHz (but could be 4 MHz). This is a feature not normally found in
  astronomical spectra but can be useful in calibration.  These tones are
  extracted from the high-resolution raw FFTs and saved in a separate table.
  """
  def __init__(self, tel):
    """
    Initialize a FITSfile_from_WVSR class
    
    This creates a table ('tbhead' and 'cols') for every connected backend and,
    if PCG tones were on, also for the PCG tones ('tonehead' and 'tonecols').
    """
    mylogger = logging.getLogger(logger.name+'.FITSfile_from_WVSR')
    FITSfile.__init__(self, tel)
    self.logger = mylogger
  
  def make_WVSR_table(self, config, collector, logserver, key, tones=True,
                      project="GBRA", observer="UNKNOWN"):
    """
    Create extension for one Backend instance from WVSR data

    The DATA column axes are::
      * CRVAL1 - data axis 1 reference value (frequency-like)
      * CRPIX1 - data axis 1 reference pixel
      * CDELT1 - data axis 1 value step
      * CRVAL2 - data axis 2 reference value (RA-like)
      * CRPIX2 - data axis 2 reference pixel
      * CDELT2 - data axis 2 value step
      * CRVAL3 - data axis 3 reference value (declination-like)
      * CRPIX3 - data axis 3 reference pixel
      * CDELT3 - data axis 3 value step
      * CRVAL4 - data axis 4 reference value (polarization code)
      * CRPIX4 - data axis 4 reference pixel
      * CDELT4 - data axis 4 value step
    The WVSR FFT post-processed data are averaged over some number of seconds
    and so there is no time axis.  Also, there is no beam axis.
    
    @param config : station configuration information ('equipment')
    @type  config : dict of dicts

    @param collector : metadata collector
    @type  collector : dict of WVSRmetadataCollector objects

    @param logserver : provides NMC data for requested time
    @type  logserver : NMC_log_server object
    
    @param key :
    @type  key : str
    
    @param tones : should tone data be extracted into a separate extension
    @type  tones : bool
    
    @param project : project ID override or default
    @type  project : str

    @param observer : team leader override or default
    @type  observer : str

    @return: pyfits.BinTableHDU instance
    """
    #  Use the first scan to get the axis length
    self.collector = collector
    self.logserver = logserver
    scans = self.collector.fft_meta.keys()
    scans.sort()
    numscans = len(scans)
    
    # the number of data cube dimensions is four following the NRAO convention.
    first_scan = scans[0]
    nchans = self.collector.fft_meta[first_scan]['n_freqs']
    # set backend attribute
    config[key]['Backend'].num_chan = nchans
    # compute number of spectrum channels for science
    subchannel_names = self.collector.wvsr_cfg[cfg_key][1]['subchannels']
    num_subchans = len(subchannel_names)
    anysubch = subchannel_names[0] # any subchannel
    
    # get the frequency and bandwidth
    obsfreq = self.collector.wvsr_cfg[cfg_key][1]['rf_to_if_lo']*1e6 \
                             + self.collector.wvsr_cfg[cfg_key][1][anysubch]['sfro']
    bandwidth = self.collector.wvsr_cfg[cfg_key][1]['chan_id 1']['bandwidth']
    
    # number of channels in the Stokes spectra
    num_Stokes_chan = get_num_chans(obsfreq, bandwidth, 0.05)
    
    # shape of data cube (spec chls, RA, dec, Stokes)
    dims = (num_Stokes_chan, 1, 1, 4)
    self.logger.debug("make_WVSR_table: science data shape is %s", dims)
    nchan = dims[0]
    nlong = dims[1]
    nlat  = dims[2]
    npols = dims[3]    
    
    # start the extension header for astronomy data
    self.exthead = self.make_basic_header()
    self.fe_key = config[key]['FrontEnd'].keys()[0] # only one FE
    self.rx_key = config[key]['Receiver'].keys()[0] # only one receiver
    self.exthead['FRONTEND'] = (config[key]['FrontEnd'][self.fe_key].name,
                          "front end ID")
    self.exthead['RECEIVER'] = (config[key]['Receiver'][self.rx_key].name,
                          "receiver ID")
    # get_hardware_metadata takes care of BACKEND
    self.exthead['projid']  = self.collector.project
    self.exthead['observer'] = self.collector.wvsr_cfg[key]['user']
    
    # add the site data (again)
    self.add_site_data(self.exthead)
    
    # make the basic columns that are always needed
    # 
    self.make_basic_columns()

    # add the backend data
    BE = config[key]['Backend'].wvsrs[key]
    self.get_hardware_metadata(BE)

    # things like LST, wind, etc.
    self.add_time_dependent_columns(1)
    
    # Add columns describing the data matrix
    #   Note that refpix defaults to 0
    axis = 1; self.make_data_axis(self.exthead, self.columns, axis,
                                  num_Stokes_chan, 'FREQ-OBS', 'D', unit='Hz',
                                comment="channel frequency in telescope frame")
    axis +=1; self.make_data_axis(self.exthead, self.columns, axis,
                                  nlong, 'RA---GLS', 'D', unit='deg',
                                  comment="RA J2000.0")
    axis +=1; self.make_data_axis(self.exthead, self.columns, axis,
                                  nlat, 'DEC--GLS','D', unit='deg',
                                  comment="decl. J2000") 
    #   Stokes axis ; get the polarizations from the spectrometer input signals
    axis+=1; self.make_data_axis(self.exthead, self.columns, axis,
                                 npols, 'STOKES',  'I',
                                 comment="polarization code: 1,2,3,4")
    # Make the data column
    fmt_multiplier = self.exthead['MAXIS1']*self.exthead['MAXIS2']* \
                     self.exthead['MAXIS3']*self.exthead['MAXIS4']
    self.logger.debug("make_WVSR_table: format multiplier = %d", fmt_multiplier)
    dimsval = "("+str(self.exthead['MAXIS1'])+"," \
                 +str(self.exthead['MAXIS2'])+"," \
                 +str(self.exthead['MAXIS3'])+"," \
                 +str(self.exthead['MAXIS4'])+")"
    self.logger.debug("make_WVSR_table: computed scan shape: %s", dimsval)
    data_format = str(fmt_multiplier)+"E"
    self.logger.debug("make_WVSR_table: data_format = %s", data_format)
    self.columns += pyfits.Column(name='DATA',
                             format=data_format, dim=dimsval)

    # add column for system temperatures
    self.columns.add_col(pyfits.Column(name='TSYS', format='2E', unit="K",
                               dim="(1,1,1,2)"))
    
    if tones:
      # start extension header for tone data
      self.make_tone_header()
      self.make_tone_columns()
      self.tonehead['backend'] = self.exthead['backend']
      self.tonehead['maxis1']  = self.exthead['maxis1']
      self.tonehead['freqres'] = self.exthead['freqres']
      # Now describe the tone data structure
      num_subch = len(subchannel_names)
      total_num_tones = 0
      for subch in subchannel_names: # count up total tones in both subchannels
        # there is one tone every MHz
        num_tones = \
            int(self.collector.wvsr_cfg[cfg_key][1][anysubch]['bandwidth']/1e6)
        total_num_tones += num_tones
      #   we take 256 channels centered on the hi-res channel nearest the tone
      toneaxis = 1; self.make_data_axis(self.tonehead, self.tonecols,
                                        toneaxis, num_tonespec_chls,
                                        'FREQ-OBS', 'D', unit='Hz',
                                comment="channel frequency in telescope frame")
      #     a Stokes spectrometer always has two input signals
      # IF1 is Stokes code -1 (RR) and IF2 is -2 (LL)
      toneaxis+=1; self.make_data_axis(self.tonehead, self.tonecols, 
                                       toneaxis, 2, 
                                       'STOKES',  'I',
                                       comment="polarization code: -2, -1")
      # Make the tone data column (MAXIS was set by make_data_axis)
      fmt_multiplier = self.tonehead['MAXIS1']*self.tonehead['MAXIS2']
      self.logger.debug("make_WVSR_table: tone format multiplier = %d",
                      fmt_multiplier)
      dimsval = "("+str(self.tonehead['MAXIS1']) + "," + \
                    str(self.tonehead['MAXIS2'])+")"
      self.logger.debug("make_WVSR_table: computed scan shape: %s", dimsval)
      data_format = str(fmt_multiplier)+"E"
      self.logger.debug("make_WVSR_table: data_format = %s", data_format)
      self.tonecols.add_col(pyfits.Column(name='TONES',
                             format=data_format, dim=dimsval))
      # columns for tone fits
      self.tonecols.add_col(pyfits.Column(name='BASELINE', format='4E',
                                          dim="(2,1,1,2)"))
      self.tonecols.add_col(pyfits.Column(name='TONEAMP',  format='4E',
                                          dim="(2,1,1,2)"))
      self.tonecols.add_col(pyfits.Column(name='TONEOFST', format='4E',
                                          dim="(2,1,1,2)"))
      self.tonecols.add_col(pyfits.Column(name='TONEWIDT', format='4E',
                                          dim="(2,1,1,2)"))
      # add IF monitor spectra to main table
      self.columns.add_col(pyfits.Column(name='IFSPECTR', format='2048E',
                                         dim="(1024,1,1,2)"))
      self.logger.debug("make_WVSR_table: columns: %s", self.columns.names)
      # create a table extension for tone data
      toneFITSrec = pyfits.FITS_rec.from_columns(self.tonecols,
                                                nrows=numscans*total_num_tones)
      tonetabhdu = pyfits.BinTableHDU(data=toneFITSrec, header=self.tonehead,
                                    name="TONES PCG")
    
    # create a table extension for astronomy data
    FITSrec = pyfits.FITS_rec.from_columns(self.columns,
                                           nrows=numscans*num_subchans)
    tabhdu =  pyfits.BinTableHDU(data=FITSrec, header=self.exthead,
                                      name="SINGLE DISH")
    if tones:
      # fill in the data and tone tables rows
      tabhdu, tonetabhdu = self.add_data(tabhdu, tonetabhdu, key)
      self.tables[key+"-pcg"] = tonetabhdu
    else:
      tabhdu, dummy = self.add_data(tabhdu, None, key)
    self.tables[key] = tabhdu
  
  def make_tone_header(self):
    """
    """
    #    I don't trust 'tbhead' because assigment doesn't necessarily make a
    # stand-alone copy.
    self.tonehead = self.make_basic_header()
    # suspect extensions are ordered alphabetically by name; TONES after SINGLE
    self.tonehead['extname'] = ("TONES_PCG", "phase cal tone extension")
    self.tonehead['FRONTEND'] = (config[key]['FrontEnd'][self.fe_key].name,
                                 "front end ID")
    self.tonehead['RECEIVER'] = (config[key]['Receiver'][self.rx_key].name,
                                 "receiver ID")
    self.tonehead['projid']  = self.collector.project
    self.tonehead['observer'] = self.collector.wvsr_cfg[key]['user']
     
  def make_tone_columns(self, numrecs=1):
    """
    Make the minimum set of columns needed by SDFITS

    This makse the REQUIRED columns for an SDFITS binary table::
    * SCAN     - scan number
    * CYCLE    - subscan number; increments by 1 for every row
    * DATE-OBS - ISO format date and time
    * EXPOSURE - integration time, sec
    * TIME     - a required FITS keyword
    * BANDWIDT - spectrometer bandwidth, Hz
    * SIDEBAND - lower or upper with respect to the LO
    * OBSFREQ  - center frequency of the receiver
    * FOFFREF1 - frequency offset for frequency switching

    To these are usually appended the various beam offsets available at DSN
    antennas.  See 'make_offset_columns'.

    @param numrecs : minimum of the number of records in each scan
    @type  numrecs : 1
    """
    # create empty column data arrays
    # create required columns.

    self.tonecols = pyfits.ColDefs([
      pyfits.Column(name='SCAN',     format='1I'             ),
      pyfits.Column(name='CYCLE',    format='1I'             ),
      pyfits.Column(name='DATE-OBS', format='16A'            ),
      pyfits.Column(name='EXPOSURE', format='1E', unit='s'   ),
      pyfits.Column(name='TIME',     format='1E', unit='s'   ),
      pyfits.Column(name='UNIXtime', format='1E', unit='s'   ),
      pyfits.Column(name='BANDWIDT', format='1E', unit='Hz'  ),
      pyfits.Column(name='SIDEBAND', format='1A'             ),
      pyfits.Column(name='OBSFREQ',  format='1D', unit='Hz'  )])
      
    # frequency switching offset
    self.tonecols += pyfits.ColDefs(
                     [pyfits.Column(name='FOFFREF1',  format='1E', unit='Hz')])

  def get_hardware_metadata(self, BE):
    """
    Initializes columns for the backends.

    Creates SDFITS columns for the backends.  Each distinct backend gets its
    own extension so BACKEND doesn't need a column.
    
    """
    self.exthead['backend'] = BE.name
    self.exthead['maxis1'] = (BE.fft_cfg['n_freqs'], "length of DATA axis 1")
    subchannel = BE.IF[1].subchannel['chan_id 1'] # any one will do
    self.exthead['freqres'] =  float(subchannel['bandwidth'])/subchannel['nchans']
  
  def add_data(self, tabhdu, tonetabhdu, cfg_key, numrecs=1):
    """
    Takes data header from WVSR post-processed FFT files and puts them into the
    SDFITS table

    Header information comes from various attributes of the 'collector'
    object::
      equip    - station configuration
      fft_meta - information about the post-processing
      fftdir   - the location of the post-processed files
      wvsr_cfg - WVSR configuration metadata
    If there are two IFs being merged to create full Stokes, then the WVSR
    configurations are the same. If there is only one IF, it will be "1".
        
    @param tabhdu : the table being filled
    @type  tabhdu : BinTableHDU object
    
    @param cfg_key : key pointing to the metadata for this back end
    @type  cfg_key : str
    """
    self.logger.debug("add_data: with %d records/scan", numrecs)
    
    scans = self.collector.fft_meta.keys()
    scans.sort()
    numscans = len(scans) # scans observed
    self.logger.debug("add_data: %d scans: %s", numscans, scans)
    
    # both IFs have the same subchannels
    subchannels = self.collector.wvsr_cfg[cfg_key][1]['subchannels']
    # create frame to return antenna RA and dec to J2000
    fk5_2000 = FK5(equinox=Time(2000, format='jyear', scale='utc'))
    
    # receiver providing the signal to the spectrometer
    rx_key = self.collector.equip[cfg_key][1]['Receiver'].keys()[0] # only one
    rx = self.collector.equip[cfg_key][1]['Receiver'][rx_key]
    self.logger.debug("add_data: receiver is %s", rx)
    
    # original number of channels
    num_chan = self.collector.equip[cfg_key][1]['Backend'].num_chan

    self.IFpower = {}
    data_row_index = 0
    if tonetabhdu:
      tone_row_index = 0
    for scan in scans:    # self.collector.scankeys:
      # dataset header is keyed on the index of the scan in the set of scans
      try:
        scan_index = scans.index(scan)
      except ValueError,details:
        self.logger.warning("add _data: scan %s not found; skipped",
                            scan)
        continue
      # use date of first record; see doc string for explanation of extra index
      year, month, day = calendar_date(self.collector.year, self.collector.doy)
      # UNIX time at midnight
      midnight = time.mktime(dateutil.parser.parse(
                          tabhdu.data[data_row_index]['DATE-OBS']).timetuple())
      subch_tone_idx = 0 # collect tone data from both subchannels
      for subch in subchannels:
        self.logger.debug(
                     "add_data: processing scan %d %s data row %d tone row %d",
                     scan, subch, data_row_index, tone_row_index)
        sub_idx = subchannels.index(subch)
        tabhdu.data[data_row_index]['SCAN'] = scan   # int
        tabhdu.data[data_row_index]['DATE-OBS'] = \
                                           "%4d/%02d/%02d" % (year, month, day)
        # each subchannel has its own cycle
        tabhdu.data[data_row_index]['CYCLE'] = sub_idx+1
        self.logger.debug("add_data: CYCLE = %d",
                          tabhdu.data[data_row_index]['CYCLE'])
        try:
          datafile = self.collector.scaninfo[scan]['subch '+str(sub_idx+1)]
        except KeyError, details:
          self.logger.warning("add_data: could not find subch %d for scan %d",
                              sub_idx+1, scan)
          continue
        # this returns a structured array with 131072 spectrum channels
        thisdata = read_FFT_file(fftdir+datafile)
        if type(thisdata) == numpy.ndarray:
          # In [26]: data.dtype.names
          # Out[26]: 
          # ('freq', 'IF1-ps', 'IF2-ps', 'IF1-phase', 'IF2-phase',
          #   'I',     'Q',      'U',      'V',           'P',
          #  'count', 'index')
          starttime = self.collector.scaninfo[scan]['start'] # datetime.datetime
          endtime = self.collector.scaninfo[scan]['end']
          self.logger.debug("add_data: for scan %d %s between %s and %s",
                            scan, subch,starttime,endtime)
          startUXtime = datetime_to_UnixTime(starttime)
          endUXtime = datetime_to_UnixTime(endtime)
          tabhdu.data[data_row_index]['UNIXtime'] = startUXtime
          tabhdu.data[data_row_index]['TIME'] = \
                               tabhdu.data[data_row_index]['UNIXtime']-midnight
          if self.collector.scaninfo[scan]['source'][-4:] == "-ref":
            tabhdu.data[data_row_index]['OBJECT'] = \
                                   self.collector.scaninfo[scan]['source'][:-4]
            tabhdu.data[data_row_index]['SIG'] = False
          else:
            tabhdu.data[data_row_index]['OBJECT'] = \
                                        self.collector.scaninfo[scan]['source']
            tabhdu.data[data_row_index]['SIG'] = True
          self.logger.debug("add_data: source is %s", 
                            tabhdu.data[data_row_index]['OBJECT'])
          response = self.logserver[cfg_key].get_azel(startUXtime, endUXtime)
          self.logger.debug("add_data: response = %s", response)
          if response:
            az,el = response
            tabhdu.data[data_row_index]['AZIMUTH'] = az
            tabhdu.data[data_row_index]['ELEVATIO'] = el
          else:
            tabhdu.data[data_row_index]['AZIMUTH'] = numpy.nan
            tabhdu.data[data_row_index]['ELEVATIO'] = numpy.nan
          tabhdu.data[data_row_index]['OBSMODE'] = obsmode
          # same exposure for all channels
          tabhdu.data[data_row_index]['EXPOSURE'] = thisdata[0]['count']
        else:
          # bad data file
          # this is probably the end of the recording
          self.logger.warning(
          "add_data: read_FFT_file return not a numpy array for scan %d %s",
            scan, subch)
          continue
        # same frequency and bandwidth for both IFs
        tabhdu.data[data_row_index]['BANDWIDT'] = \
                  self.collector.wvsr_cfg[cfg_key][1]['chan_id 1']['bandwidth']
        obsfreq = self.collector.wvsr_cfg[cfg_key][1]['rf_to_if_lo']*1e6 \
                           + self.collector.wvsr_cfg[cfg_key][1][subch]['sfro']
        tabhdu.data[data_row_index]['OBSFREQ'] = obsfreq
        tabhdu.data[data_row_index]['RESTFREQ'] = obsfreq # is this always true?
        self.logger.debug("add_data: OBJECT is '%s'",
                          tabhdu.data[data_row_index]['OBJECT'])
        sourcename = tabhdu.data[data_row_index]['OBJECT'].replace('_',' ')
        self.logger.debug("add_data: OBJECT is '%s'",
                          tabhdu.data[data_row_index]['OBJECT'])
        tabhdu.data[data_row_index]['VELOCITY'] = \
                                     self.collector.sources[sourcename]['Vlsr']
        tabhdu.data[data_row_index]['VELDEF'] = veldef
        weather = self.logserver[cfg_key].get_weather(startUXtime)
        self.logger.debug("add_data: weather at %s is %s", startUXtime, weather)
        if weather:
          tabhdu.data[data_row_index]['TAMBIENT'] = weather[0]
          tabhdu.data[data_row_index]['PRESSURE'] = weather[1]
          tabhdu.data[data_row_index]['HUMIDITY'] = weather[2]
          tabhdu.data[data_row_index]['WINDSPEE'] = weather[3]
          tabhdu.data[data_row_index]['WINDDIRE'] = weather[4]
        else:
          self.logger.debug("add_data: weather not available for %f", startUXtime)
          tabhdu.data[data_row_index]['TAMBIENT'] = numpy.nan
          tabhdu.data[data_row_index]['PRESSURE'] = numpy.nan
          tabhdu.data[data_row_index]['HUMIDITY'] = numpy.nan
          tabhdu.data[data_row_index]['WINDSPEE'] = numpy.nan
          tabhdu.data[data_row_index]['WINDDIRE'] = numpy.nan
        
        datacubeshape = tabhdu.data[data_row_index]['DATA'].shape
        # the frequency axis is first in FITS/FORTRAN order and last (of four)
        # in PYTHON/C order
        num_Stokes_chan = datacubeshape[3]
        # the sign of CDELT1 depends on the sideband of the last SSB mixer in
        # the chain
        # GBT SDFITS wants SIDEBAND
        if rx['IFmode'] == 'U':
          tabhdu.data[data_row_index]['SIDEBAND'] = +1
        elif rx['IFmode'] == 'L':
          tabhdu.data[data_row_index]['SIDEBAND'] = -1
        else:
          self.logger.error("add_data: IF mode %s is invalid; default to USB",
                            rx['IFmode'])
          tabhdu.data[data_row_index]['SIDEBAND'] = +1
             
        # second and third data axes (coordinates)
        RA, dec = \
          self.logserver[cfg_key].get_RAdec(startUXtime)
        self.logger.debug("add_data: apparent RA,dec = %f,%f", RA, dec)
        c = SkyCoord(RA, dec, unit=(u.deg, u.deg),
                     frame=FK5(equinox=Time('J'+str(year), scale='utc')))
        c2000 = c.transform_to(fk5_2000)
        self.logger.debug("add_data: RA,dec = %f,%f", c.ra.hour, c.dec.deg)
        tabhdu.data[data_row_index]['CRVAL2'] = c2000.ra.hour # hours
        tabhdu.data[data_row_index]['CRVAL3'] = c2000.dec.deg    # deg
        tabhdu.data[data_row_index]['EQUINOX'] = 2000
        # get the radial velocity of the LSR
        tabhdu.data[data_row_index]['VFRAME'] = \
                V_LSR(c2000.ra.hour, c2000.dec.deg, self.tel.number, starttime)
        tabhdu.data[data_row_index]['RVSYS'] = \
                                      tabhdu.data[data_row_index]['VELOCITY'] \
                                    + tabhdu.data[data_row_index]['VFRAME']
    
        # fourth data axis (polarization)
        tabhdu.data[data_row_index]['CRVAL4'] = -1
        tabhdu.data[data_row_index]['CDELT4'] = -1 # for I,Q,U,V (-1,-2,-3,-4)
        
        # initialize power averages
        self.IFpower[data_row_index] = {'IF1': thisdata['IF1-ps'],
                                        'IF2': thisdata['IF2-ps']}
        # fit the tones to the data using the original resolution
        bandwidth = tabhdu.data[data_row_index]['BANDWIDT']
        if tonetabhdu:
          tone_offsets, tone_chnls = tone_chnl_nums(num_chan, obsfreq,
                                                    bandwidth)
          self.logger.debug("add_data: tone offsets: %s", tone_offsets)
          self.logger.debug("add_data: tone channels: %s", tone_chnls)
          tone_indices = list(tone_chnls) # define the channels to be selected
          num_tones = len(tone_indices)
        offset, center_tone = math.modf(obsfreq/1e6) # tones every MHz
        # the objective is to fit the position of one (e.g. central) tone
        # and one std with all the other tones at a fixed distance from the
        # central tone and one std for all tones.  The individual amplitudes
        # may vary. 'multigauss.other_pars' has the fixed parameters.
        for IF in ['IF1', 'IF2']:
          self.logger.debug("add_data: processing %s", IF)
          IFidx = int(IF[-1])-1
          IFpwr = self.IFpower[data_row_index][IF]
          
          if tonetabhdu:
            for tone in tone_chnls:
              tonetabhdu.data[tone_row_index]['SCAN'] = scan   # int
              self.logger.debug("add_data: processing tone(%d) is %d",
                                subch_tone_idx, tone)
              # CYCLE increments by 1 for each row in the scan
              tonetabhdu.data[tone_row_index]['CYCLE'] = subch_tone_idx + 1
              # all the following are the same for each tone in the subchannel
              tonetabhdu.data[tone_row_index]['DATE-OBS'] = \
                                           "%4d/%02d/%02d" % (year, month, day)
              tonetabhdu.data[tone_row_index]['UNIXtime'] = startUXtime
              tonetabhdu.data[tone_row_index]['TIME'] = \
                           tonetabhdu.data[tone_row_index]['UNIXtime']-midnight
              tonetabhdu.data[tone_row_index]['EXPOSURE'] = thisdata[0]['count']
              tonetabhdu.data[tone_row_index]['BANDWIDT'] = \
                                     num_tonespec_chls*tabhdu.header['FREQRES']
              tonetabhdu.data[tone_row_index]['OBSFREQ'] = obsfreq + \
                                          (tone-65536)*tabhdu.header['FREQRES']
              tonetabhdu.data[tone_row_index]['CDELT1'] = \
                                                       tabhdu.header['FREQRES']
              self.logger.debug("add_data: freq step is %f",
                              tonetabhdu.data[tone_row_index]['CDELT1'])
              halfband = num_tonespec_chls/2 # half the number of tone band chls
              tonetabhdu.data[tone_row_index]['CRPIX1'] = halfband
              # 
              # the frequency of the channel nearest to the tone
              nearest_chnl_freq = obsfreq + \
                   (tone-65536)*tonetabhdu.data[tone_row_index]['CDELT1']
              tonefreq_kHz = (nearest_chnl_freq)/1000
              tonefreq = 1000*round(tonefreq_kHz)
              # center of the 16 channel spectrum extract
              tonetabhdu.data[tone_row_index]['CRVAL1'] = nearest_chnl_freq
              self.logger.debug("add_data: tone frequency is %f", tonefreq)
              # get the data around the tone.
              tonetabhdu.data[tone_row_index]['TONES'][IFidx,:] = \
                                thisdata[IF+'-ps'][tone-halfband:tone+halfband]
              # now fit the tone (freq in MHz)
              toneband_channels = numpy.arange(num_tonespec_chls)-halfband
              self.logger.debug("add_data: tone band channel numbers: %s",
                                                             toneband_channels)
              x = (tonetabhdu.data[tone_row_index]['CRVAL1'] + \
                                  toneband_channels * \
                                  tonetabhdu.data[tone_row_index]['CDELT1'])/1e6
              y = tonetabhdu.data[tone_row_index]['TONES'][IFidx,:]
              self.logger.debug("add_data: x = %s", x)
              self.logger.debug("add_data: y = %s", y)
              if not numpy.any(y):
                # y is all zeros.  Give up
                tonetabhdu.data[tone_row_index]['BASELINE'][IFidx,0,0,:] = \
                                                                       nanarray
                tonetabhdu.data[tone_row_index]['TONEAMP'][IFidx,0,0,:]  = \
                                                                       nanarray
                tonetabhdu.data[tone_row_index]['TONEOFST'][IFidx,0,0,:] = \
                                                                       nanarray
                tonetabhdu.data[tone_row_index]['TONEWIDT'][IFidx,0,0,:] = \
                                                                       nanarray
                continue
              est_bias = numpy.append(y[:5], y[-6:]).mean()
              initial_guess = (est_bias, y[halfband], tonefreq/1e6,
                             tabhdu.header['FREQRES']/1e6)
              self.logger.debug("add_data: initial_guess = %s", initial_guess)
              popt, pcov = curve_fit(biased_scaled_sinc, x, y,
                                   p0=(initial_guess))
              self.logger.debug("add_data: pars = %s", popt)
              self.logger.debug("add_data: covars = %s", pcov)
              bias, amp, offset, std = popt
              dbias, damp, doffset, dstd = numpy.sqrt(numpy.diag(pcov))
              tonetabhdu.data[tone_row_index]['BASELINE'][IFidx,0,0,:] = \
                                                     numpy.array([bias, dbias])
              tonetabhdu.data[tone_row_index]['TONEAMP'][IFidx,0,0,:] = \
                                                       numpy.array([amp, damp])
              tonetabhdu.data[tone_row_index]['TONEOFST'][IFidx,0,0,:] = \
                                                 numpy.array([offset, doffset])
              tonetabhdu.data[tone_row_index]['TONEWIDT'][IFidx,0,0,:] = \
                                                       numpy.array([std, dstd])
              # remove the tones from the IF power data
              #    we need the map the channels for 'x' into the 
              IFpwr[tone + toneband_channels] -= \
                                biased_scaled_sinc(x, bias, amp, offset, std)
            
              if IF == "IF2": # should be [list of IFs][-1] !!!!!!!!!!!!!!!!!!!
                # because both IF go into the same row but different positions on
                # the POL axis
                tone_row_index += 1
                subch_tone_idx += 1
              # end of tone loop
            # save smoothed IF power spectra
            # axis specs not needed; simplest set for relative frequenies is::
            newspec, newrefval, newrefpix, newdelta = \
                        reduce_spectrum_channels(IFpwr, 0, 0, 0,
                                                 num_chan=1024)
            tabhdu.data[data_row_index]['IFSPECTR'][IFidx, 0, 0,:] = newspec
            # end of tones section
          # compute the average power
          tabhdu.data[data_row_index]['TSYS'][IFidx,0,0,0] = IFpwr.mean()
          tsys_col_idx = tabhdu.data.columns.names.index('TSYS')
          tabhdu.data.columns['TSYS'].unit = "count"
          # end of IF loop
        # the data in dataset is keyed on scan number
        refval = tabhdu.data[data_row_index]['OBSFREQ']
        refpix = num_chan/2
        delta  = bandwidth/num_chan
        self.logger.debug("add_data: loading DATA")
        I, newrefval, newrefpix, newdelta = \
                 reduce_spectrum_channels(thisdata['I'], refval, refpix, delta,
                                          num_chan=num_Stokes_chan)
        Q, newrefval, newrefpix, newdelta = \
                 reduce_spectrum_channels(thisdata['Q'], refval, refpix, delta,
                                          num_chan=num_Stokes_chan)
        U, newrefval, newrefpix, newdelta = \
                 reduce_spectrum_channels(thisdata['U'], refval, refpix, delta,
                                          num_chan=num_Stokes_chan)
        V, newrefval, newrefpix, newdelta = \
                 reduce_spectrum_channels(thisdata['V'], refval, refpix, delta,
                                           num_chan=num_Stokes_chan)
        tabhdu.data[data_row_index]['DATA'][0, 0, 0,:] = I
        tabhdu.data[data_row_index]['DATA'][1, 0, 0,:] = Q
        tabhdu.data[data_row_index]['DATA'][2, 0, 0,:] = U
        tabhdu.data[data_row_index]['DATA'][3, 0, 0,:] = V
        tabhdu.data[data_row_index]['CRVAL1'] = newrefval + delta/2
        tabhdu.data[data_row_index]['CRPIX1'] = newrefpix
        tabhdu.data[data_row_index]['CDELT1'] = newdelta
        self.logger.info("add_data: finished row %d scan %d cycle %d",
                         data_row_index, tabhdu.data[data_row_index]['SCAN'],
                                         tabhdu.data[data_row_index]['CYCLE'])
        data_row_index += 1
        # end subch loop
      # end scan loop
    return tabhdu, tonetabhdu
    
if __name__ == "__main__":
  examples = """
Examples
========
  run WVSR2SDFITS.py --console_loglevel=debug --DSS=14 \
                          --project=AUTO_EGG --date=2016/237
"""  
  p = initiate_option_parser(__doc__, examples)
  p.usage='WVSR2SDFITS.py [options]'
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
  p.add_argument('-p', '--project',
               dest = 'project',
               type = str,
               default = 'AUTO_EGG',
               help = "Project code")
  args = p.parse_args(sys.argv[1:])
  
  mylogger = logging.getLogger()
  mylogger.setLevel(logging.DEBUG)
  init_logging(mylogger,
                 loglevel = get_loglevel(args.file_loglevel),
                 consolevel = get_loglevel(args.console_loglevel),
                 logname = args.logpath+"WVSR2SDFITS.log")
                 
  mylogger.critical(" WVSR2SDFITS started")
  mylogger.debug("WVSR2SDITS args: %s", args)
  if args.project[:4] != 'AUTO':
    raise RuntimeError("project name must start with AUTO")
  yearstr, doystr = args.date.split("/")
  year = int(yearstr)
  doy = int(doystr)
  
  # note that this does not handle recording sessions with multiple antennas
  obsdir, realdir, project_dir, datadir, wvsrdir, fftdir = \
                            get_session_dirs(args.project, args.dss, year, doy)

  # get data for all the sources and verifiers used by the project
  sourcedata = get_all_source_data(project_dir)
  
  # get the start and end time of the session
  timesfiles = glob.glob(obsdir+"times-*")
  mylogger.debug("WVSR2SDITS found %s", timesfiles)
  if len(timesfiles) < 1:
    raise RuntimeError("WVSR2SDITS: no times file is %s" % obsdir)
  elif len(timesfiles) > 1:
    raise RuntimeError("WVSR2SDITS: can only handle one timesfile; %s has %d" \
                       % len(timesfiles))
  starttime = get_start_time(timesfiles[0])

  # get a metadata manager for the WVSR log for this session
  collector = WVSRmetadataCollector(args.project, args.dss, year, doy, starttime)
  collector.sources = sourcedata
  mylogger.debug("__init__: equipment: %s", collector.equip)
  # add Backend to the standard equipment
  #   at this point, the collector only knows about backend equipment so
  # collector.equip.keys(), collector.wvsrnames.keys() and
  # collector.wvsr_cfg.keys() are all the same.
  config = {}
  NMClogmanager = {}
  NMClogserver  = {}
  rxband = {}
  for wvsr in collector.wvsrnames:
    # how many BINTABLEs do we need?
    if collector.wvsr_cfg[wvsr].has_key(1) and \
       collector.wvsr_cfg[wvsr].has_key(2):
      # are the two IFs the same except for polarization?
      dss1,band1,pol1 =  collector.wvsr_cfg[wvsr][1]['IF_source'].split('_')
      dss2,band2,pol2 =  collector.wvsr_cfg[wvsr][2]['IF_source'].split('_')
      if dss1 == dss2 and band1 == band2 and pol1[0] != pol2[0]:
        # yes; treat as one signal source with two pols
        dss = int(dss1)
        key = wvsr
        rxband[key] = band1
      else:
        # we can handle different bands but not antennas
        if dss1 == dss2:
          for IF in collector.wvsr_cfg[wvsr].keys():
            key = wvsr+"-"+str(IF)
            dss,band,pol = collector.wvsr_cfg[wvsr][IF]['IF_source'].split('_')
            rxband[key] = band
        else:
          raise RuntimeError("WVSR2SDFITS can only handle one antenna")
      config[key] = station_configuration(None, args.project, dss, year, doy,
                                          starttime, rxband[key])
      # get a manager for the NMC log for this session.
      NMClogmanager[key] = NMClogManager(project=args.project, station=dss,
                                         year=year, DOY=doy, starttime=starttime,
                                         use_portal=False)
      mylogger.debug("WVSR2SDITS NMC metadata available: %s",
                        NMClogmanager[wvsr].metadata.keys())
      NMClogserver[key] = NMC_log_server(args.project, dss, year, doy)
    else:
      raise RuntimeError("single IF case not yet coded")
  # Now we can start a FITSfile
  tel = config[config.keys()[0]]['Telescope'][args.dss]
  ff = FITSfile_from_WVSR(tel)
  
  # make an extension for every backend and start table headers and columns.
  # Because the extensions are for different back ends for the same observing
  # session, the number of scans in each is the same
  
  # The binary tables are constructed from datasets
  for cfg_key in config.keys():
    # each needs a separate table
    if cfg_key[-2] == '-':
      # strip off the IF
      IF = int(cfg_key[-1])
      wvsr = cfg_key[:-2]
    else:
      wvsr =  cfg_key
    ff.make_WVSR_table(config, collector, NMClogserver, cfg_key)
  
  hdulist = pyfits.HDUList([ff.prihdu]+ff.tables.values())
  parts = realdir.split('/')
  parts.remove('projects')
  parts[3] = 'RA_data'; parts[4] = 'FITS'
  fitspath = "/".join(parts)
  mkdir_if_needed(fitspath)
  fname = fitspath + "WVSR" + "_%4d-%03d-%s.fits" % (year, doy, starttime)
  mylogger.debug("WVSR2SDFITS writing %s", fname)
  hdulist.writeto(fname, clobber=True)
  mylogger.critical(" WVSR2SDFITS ended")
