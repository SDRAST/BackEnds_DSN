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
print "IPython workaround done"
import astropy.io.fits as pyfits
import astropy.units as u
import dateutil
import datetime
import glob
import logging
import math
import numpy
import sys
import time
import traceback
import warnings

from astropy.coordinates import FK5, SkyCoord
from astropy.time import Time
from os import chdir, getcwd, symlink, makedirs
from os.path import basename, exists
from scipy.optimize import curve_fit

from Astronomy import calendar_date
from Astronomy.DSN_coordinates import DSS
from Astronomy.redshift import V_LSR
from Automation import activity_project, get_session_dirs
from Automation import get_start_time, get_end_time
from Automation.NMClogs import NMClogManager,NMC_log_server 
from Automation.sources import get_all_source_data
from Automation.tones import tone_chnl_nums
from Automation.WVSR import get_Naudet_FFT_dir, make_datadir_name
from Data_Reduction import get_obs_dirs, get_obs_session, select_data_files
from Data_Reduction import get_num_chans, reduce_spectrum_channels
from Data_Reduction.DSN.WVSR.collector import WVSRmetadataCollector
from Data_Reduction.DSN.WVSR.SpecData import get_one_spectrum,read_FFT_file
from Data_Reduction.FITS.DSNFITS import FITSfile, get_row
from Data_Reduction.tipping import airmass
from DatesTimes import datetime_to_UnixTime
from local_dirs import sci_proj_path
from Math.multigauss import gaussian
#from MonitorControl.Configurations.DSN_standard import standard_equipment
from MonitorControl.Configurations.GDSCC.WVSR import station_configuration
from support import mkdir_if_needed
from support.lists import unique
from support.logs import initiate_option_parser, init_logging
from support.logs import get_loglevel, set_loglevel

import Math.multigauss as multigauss

obsmode = 'LINEPSSW'
veldef = 'RADI-OBS'

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
  
  Public attributes::
    collector - WVSRmetadataCollector object
    columns   - columns in the SINGLE DISH table
    exthead   - header for the SINGLE DISH table
    IFpower   - reduced resolution spectrum for monitoring IF
    logger    - logging.Logger object
    logserver - NMC_log_server object to fetch NMC metadata from file
    oe_end    - end of antenna dwell time for a scan from SOE
    oe_source - source for a scan from SOE
    oe_start  - start of antenna dwell time for a scan from SOE
    tables    - dict of pyfits.BinTableHDU objects
    tonehead  - header for the PCG TONES table
    
  Methods::
    add_data              -
    get_dwell_times       -
    get_hardware_metadata -
    make_tone_columns     -
    make_tone_header      -
    
  Inherited from FITSfile::
    add_site_data              -
    add_time_dependent_columns - 
    make_basic_columns         -
    make_data_axis             -
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
  
  def make_WVSR_table(self, config, collector, logserver, key,
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
    
    @param key : configuration key (index for 'config')
    @type  key : str
    
    @param project : project ID override or default
    @type  project : str

    @param observer : team leader override or default
    @type  observer : str

    @return: pyfits.BinTableHDU instance
    """
    #  Use the first scan to get the axis length
    self.collector = collector
    self.logserver = logserver
    self.scans = self.collector.fft_meta.keys()
    self.scans.sort()
    
    # the number of data cube dimensions is four following the NRAO convention.
    first_scan = self.scans[0]
    nchans = self.collector.fft_meta[first_scan]['n_freqs']
    # set backend attribute
    config[key]['Backend'].num_chan = nchans
    # compute number of spectrum channels for science
    subchannel_names = self.collector.wvsr_cfg[cfg_key][1]['subchannels']
    self.logger.debug("make_WVSR_table: subchannel_names: %s",subchannel_names)
    # initialize a SINGLE DISH record array
    FITSrec = self.init_singledish_table(config, collector)
    # provisionally, initialize a TONES PCG record array
    toneFITSrec = self.init_pcg_table()
    # now add the data using the configuration data given by 'key'
    FITSrec, toneFITSrec = self.add_data(FITSrec, toneFITSrec, key)
    # make a table from the FITS records
    self.tables[key] = self.make_table_HDU("SINGLE DISH",
                                           FITSrec, self.columns, self.exthead)
    try:
      # make a table from the tone record array
      self.tables[key+"-pcg"] = self.make_table_HDU("TONES PCG",
                                                    toneFITSrec,
                                                  self.tonecols, self.tonehead)
    except:
      self.logger.info("make_WVSR_table: no PCG tones record array")
  
  def make_table_HDU(self, extname, FITSrecords, columns, header):
    """
    Converts FITS record array into a binary table HDU
    
    Empty rows are removed
    
    @param extname : extension type, e.g. SINGLE DISK
    @type  extname : str
    
    @param FITSrecords : a FITS record array
    @type  FITSrecords : pyfits.fitsrec.FITS_rec object
    
    @param columns : column definitions
    @type  columns : pyfits.column.ColDefs objects
    
    @param header : table header
    @type  header : pyfits.Header object
    """
    # get the number of rows used
    nrows = len(FITSrecords['SCAN'].nonzero()[0])
    # create a new FITS record array with the right number of rows
    newFITSrec = pyfits.FITS_rec.from_columns(columns, nrows=nrows)
    # copy rows to the new record array
    for row in range(nrows):
      newFITSrec[row] = FITSrecords[row]
    # create the HDU
    tabhdu =  pyfits.BinTableHDU(data=newFITSrec, header=header, name=extname)
    return tabhdu
    
  def init_singledish_table(self, config, collector):
    """
    """
    subchannel_names = self.collector.wvsr_cfg[cfg_key][1]['subchannels']
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
    self.make_basic_columns()

    # add the backend data
    BE = config[key]['Backend'].wvsrs[key]
    self.get_hardware_metadata(BE)

    # things like LST, wind, etc.
    self.add_time_dependent_columns(1)
    
    # Add columns describing the data matrix
    num_subchans = len(subchannel_names)
    anysubch = subchannel_names[0] # any subchannel
    #     shape of data cube (num spec chls, num RA, num dec, num Stokes pars)
    #     num spec chls depends on the resolution in km/s
    #     to compute it use any subchannel
    obsfreq = self.collector.wvsr_cfg[cfg_key][1]['rf_to_if_lo']*1e6 \
                             + self.collector.wvsr_cfg[cfg_key][1][anysubch]['sfro']
    self.bandwidth = self.collector.wvsr_cfg[cfg_key][1]['chan_id 1']['bandwidth']
    dims = (get_num_chans(obsfreq, self.bandwidth, 0.05), 1, 1, 4)
    self.logger.debug("init_singledish_table: science data shape is %s", dims)
    #     Note that refpix defaults to 0
    axis = 1; self.make_data_axis(self.exthead, self.columns, axis,
                                  dims[0], 'FREQ-OBS', 'D', unit='Hz',
                                comment="channel frequency in telescope frame")
    axis +=1; self.make_data_axis(self.exthead, self.columns, axis,
                                  dims[1], 'RA---GLS', 'D', unit='deg',
                                  comment="RA J2000.0")
    axis +=1; self.make_data_axis(self.exthead, self.columns, axis,
                                  dims[2], 'DEC--GLS','D', unit='deg',
                                  comment="decl. J2000") 
    #   Stokes axis ; get the polarizations from the spectrometer input signals
    axis+=1; self.make_data_axis(self.exthead, self.columns, axis,
                                 dims[3], 'STOKES',  'I',
                                 comment="polarization code: 1,2,3,4")
    # Make the data column
    fmt_multiplier = self.exthead['MAXIS1']*self.exthead['MAXIS2']* \
                     self.exthead['MAXIS3']*self.exthead['MAXIS4']
    self.logger.debug("init_singledish_table: format multiplier = %d", fmt_multiplier)
    dimsval = "("+str(self.exthead['MAXIS1'])+"," \
                 +str(self.exthead['MAXIS2'])+"," \
                 +str(self.exthead['MAXIS3'])+"," \
                 +str(self.exthead['MAXIS4'])+")"
    self.logger.debug("init_singledish_table: computed scan shape: %s", dimsval)
    data_format = str(fmt_multiplier)+"E"
    self.logger.debug("init_singledish_table: data_format = %s", data_format)
    self.columns += pyfits.Column(name='DATA',
                             format=data_format, dim=dimsval)

    # add column for system temperatures
    self.columns.add_col(pyfits.Column(name='TSYS', format='2E', unit="K",
                               dim="(1,1,1,2)"))
    # add column for IF spectra
    self.columns.add_col(pyfits.Column(name='IFSPECTR', format='2048E',
                                         dim="(1024,1,1,2)"))
    self.logger.debug("init_singledish_table: columns: %s", self.columns.names)
    
    # create a structured numpy record for the data
    n_rows = len(self.scans)*num_subchans
    
    FITSrec = pyfits.FITS_rec.from_columns(self.columns, nrows=n_rows)
    return FITSrec
      
  def init_pcg_table(self, num_tonespec_chls=16):
    """
    """
    subchannel_names = self.collector.wvsr_cfg[cfg_key][1]['subchannels']
    self.num_tonespec_chls = num_tonespec_chls
    # start extension header for tone data
    self.make_tone_header()
    self.make_tone_columns()
    self.tonehead['backend'] = self.exthead['backend']
    self.tonehead['maxis1']  = self.exthead['maxis1']
    self.tonehead['freqres'] = self.exthead['freqres']
    # Now describe the tone data structure
    num_subch = len(subchannel_names)
    total_num_tones = 0
    anysubch = subchannel_names[0] # any subchannel
    for subch in subchannel_names: # count up total tones in both subchannels
      # there is one tone every MHz
      num_tones = \
            int(self.collector.wvsr_cfg[cfg_key][1][subch]['bandwidth']/1e6)
      total_num_tones += num_tones
    self.logger.debug("init_pcg_table: total number of tones is %d",
                      total_num_tones)
    #   we take 16 channels centered on the hi-res channel nearest the tone
    toneaxis = 1; self.make_data_axis(self.tonehead, self.tonecols,
                                      toneaxis, self.num_tonespec_chls,
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
    self.logger.debug("init_pcg_table: tone format multiplier = %d",
                      fmt_multiplier)
    dimsval = "("+str(self.tonehead['MAXIS1']) + "," + \
                                               str(self.tonehead['MAXIS2'])+")"
    self.logger.debug("init_pcg_table: computed scan shape: %s", dimsval)
    data_format = str(fmt_multiplier)+"E"
    self.logger.debug("init_pcg_table: data_format = %s", data_format)
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
    self.logger.debug("init_pcg_table: columns: %s", self.columns.names)
    n_rows = len(self.scans)*total_num_tones
    self.logger.debug("init_pcg_table: making %s rows", n_rows)
    toneFITSrec = pyfits.FITS_rec.from_columns(self.tonecols, nrows=n_rows)
    return toneFITSrec
 
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
  
  def add_data(self, FITSrec, toneFITSrec, cfg_key):
    """
    Takes data header from WVSR post-processed FFT files and puts them into the
    SDFITS table

    Metadata Sources
    ================
    Header information comes from various attributes of the 'collector'
    object::
      equip    - station configuration
      fft_meta - information about the post-processing
      fftdir   - the location of the post-processed files
      wvsr_cfg - WVSR configuration metadata
    If there are two IFs being merged to create full Stokes, then the WVSR
    configurations are the same. If there is only one IF, it will be "1".
    
    Psuedo-code Summary
    ===================
    data_row_index = 0 # one row for each scan, cycle
    tone_row_index = 0 # one row for each scan, cycle and tone
    for each scan:
      isubch_tone_idx = 0 # index for the tones in this scan (all subchannels)
      for each subchannel:
      add some data to columns of FITS record array (SINGLE DISH)
      read in the data from FFT file
      if the data are not valid:
        skip this file
      add more data to columns of FITS record array (SINGLE DISH)
      compute DATA cube axes parameters
      compute expected number of tones
      for each IF:
        get the spectrum
        check for tone rails in the spectrum
        if more than one rail:
          skip this IF
        add the tone data to the TONES PCG record array
          get some metadata
          for each tone:
            set the CYCLE number
            fill in other columns for this row
            compute the frequency of tone, rounded to kHz, in Hz
            get x and y data from IF spectrum
            fit the data
            store the fit results in columns of row
            remove fitted tone from IF power
            increment subchannel tone index (subch_tone_idx += 1)
            increment tone_row_index (tone_row_index += 1)
      average adjacent channels to get desired resolution
      put STOKES data into DATA cube
      increment data row index (data_row_index += 1)
        
    @param FITSrec : FITS record array for SINGLE DISH
    @type  FITSrec : pyfits.fitsrec.FITS_rec object
    
    @param toneFITSrec : FITS record array for TONES PCG
    @type  toneFITSrec : pyfits.fitsrec.FITS_rec object
    
    @param cfg_key : key pointing to the metadata for this back end
    @type  cfg_key : str
    """
    scans = self.collector.fft_meta.keys()
    scans.sort()
    numscans = len(scans) # scans observed
    self.logger.debug("add_data: %d scans: %s", numscans, scans)
    
    # both IFs have the same subchannels so use IF 1.
    subchannels = self.collector.wvsr_cfg[cfg_key][1]['subchannels']
    self.logger.debug("add_data: subchannels: %s", subchannels)
    anysubch = subchannels[0] # any subchannel
    # create frame to return antenna RA and dec to J2000
    fk5_2000 = FK5(equinox=Time(2000, format='jyear', scale='utc'))
    
    # receiver providing the signal to the spectrometer
    rx_key = self.collector.equip[cfg_key][1]['Receiver'].keys()[0] # only one
    rx = self.collector.equip[cfg_key][1]['Receiver'][rx_key]
    self.logger.debug("add_data: receiver is %s", rx)
    
    # original number of channels
    num_chan = self.collector.equip[cfg_key][1]['Backend'].num_chan

    bad_tones = []
    #data_row_index = 0 # one row for every scan and every cycle (subchannel)
    #tone_row_index = 0 #
    have_tones = False
    #IFpower = {}
    # fill up the rows scan-by-scan
    for scan in scans:
      # dataset header is keyed on the index of the scan in the set of scans
      try:
        scan_index = scans.index(scan)
      except ValueError,details:
        self.logger.warning("add _data: scan %s not found; skipped",
                            scan)
        continue
      # use date of first record; see doc string for explanation of extra index
      year, month, day = calendar_date(self.collector.year, self.collector.doy)
      date_obs = "%4d/%02d/%02d" % (year, month, day)
      fyear = self.collector.year + self.collector.doy/365.25
      self.logger.debug("add_data: for scan %d, %s is J%f", scan, date_obs, fyear)
      #subch_tone_idx = 0 # collect tone data from both subchannels
      # add a CYCLE row for every WVSR subchannel
      for subch in subchannels:
        sub_idx = subchannels.index(subch)
        try:
          datafile = self.collector.scaninfo[scan]['subch '+str(sub_idx+1)]
        except KeyError, details:
          self.logger.warning("add_data: could not find subch %d for scan %d",
                              sub_idx+1, scan)
          continue
        # this returns a structured array with 131072 spectrum channels
        thisdata = read_FFT_file(fftdir+datafile)
        # check on some reasons for discarding these data
        if type(thisdata) != numpy.ndarray:
          # bad data file
          # this is probably the end of the recording
          self.logger.warning(
          "add_data: read_FFT_file return not a numpy array for scan %d %s",
            scan, subch)
          continue
        cycle = sub_idx + 1
        data_row_index = get_row("SINGLE DISH", scans, scan=scan,
                                  num_cycles=len(subchannels), cycle=cycle)
        #self.logger.debug(
        #             "add_data: processing scan %d %s data row %d tone row %d",
        #             scan, subch, data_row_index, tone_row_index)
        self.logger.debug("add_data: processing scan %d subch %s data row %d",
                          scan, subch, data_row_index)
        FITSrec[data_row_index]['SCAN'] = scan   # int
        FITSrec[data_row_index]['DATE-OBS'] = date_obs         
        # UNIX time at midnight
        midnight = time.mktime(dateutil.parser.parse(date_obs).timetuple())
        # each subchannel has its own cycle
        FITSrec[data_row_index]['CYCLE'] = cycle
        self.logger.debug("add_data: CYCLE = %d", 
                          FITSrec[data_row_index]['CYCLE'])
        # In [26]: data.dtype.names
        # Out[26]: 
        # ('freq', 'IF1-ps', 'IF2-ps', 'IF1-phase', 'IF2-phase',
        #   'I',     'Q',      'U',      'V',           'P',
        #  'count', 'index')
        if self.collector.scaninfo[scan] == {}:
          # skip this scan
          self.logger.warning("add_data: no scan info for %d", scan)
          continue
        try:
          starttime = self.collector.scaninfo[scan]['start']
        except KeyError:
          # incomplete scan info
          continue
        # process the data
        endtime = self.collector.scaninfo[scan]['end']
        self.logger.debug("add_data: for scan %d subch %s between %s and %s",
                          scan, subch, starttime, endtime)
        startUXtime = datetime_to_UnixTime(starttime)
        endUXtime = datetime_to_UnixTime(endtime)
          # IFs used for this scan? One (single pol) or two (both pols)
        if numpy.any(thisdata['IF2-ps'] != 0):
          IFs = ['IF1', 'IF2']
        else:
          IFs = ['IF1']
        # put data in row
        FITSrec[data_row_index]['UNIXtime'] = startUXtime
        if startUXtime == 0.0:
          FITSrec[data_row_index]['CYCLE'] = 0
        FITSrec[data_row_index]['TIME'] = \
                                   FITSrec[data_row_index]['UNIXtime']-midnight
        if self.collector.scaninfo[scan]['source'][-4:] == "-ref":
          FITSrec[data_row_index]['OBJECT'] = \
                                   self.collector.scaninfo[scan]['source'][:-4]
          FITSrec[data_row_index]['SIG'] = False
        else:
          FITSrec[data_row_index]['OBJECT'] = \
                                        self.collector.scaninfo[scan]['source']
          FITSrec[data_row_index]['SIG'] = True
        self.logger.debug("add_data: source is %s", 
                          FITSrec[data_row_index]['OBJECT'])
        response = self.logserver.get_azel(startUXtime, endUXtime)
        self.logger.debug("add_data: response = %s", response)
        if response:
          az,el = response
          FITSrec[data_row_index]['AZIMUTH'] = az
          FITSrec[data_row_index]['ELEVATIO'] = el
        else:
          FITSrec[data_row_index]['AZIMUTH'] = numpy.nan
          FITSrec[data_row_index]['ELEVATIO'] = numpy.nan
        FITSrec[data_row_index]['OBSMODE'] = obsmode
        # same exposure for all channels
        FITSrec[data_row_index]['EXPOSURE'] = thisdata[0]['count']
        
        # same frequency and bandwidth for both IFs so use IF 1
        obsfreq = self.collector.wvsr_cfg[cfg_key][1]['rf_to_if_lo']*1e6 \
                             + self.collector.wvsr_cfg[cfg_key][1][subch]['sfro']

        FITSrec[data_row_index]['BANDWIDT'] = \
                  self.collector.wvsr_cfg[cfg_key][1]['chan_id 1']['bandwidth']
        # add the data to the columns of this row
        FITSrec[data_row_index]['OBSFREQ'] = obsfreq
        FITSrec[data_row_index]['RESTFREQ'] = obsfreq # is this always true?
        self.logger.debug("add_data: OBJECT is '%s'",
                          FITSrec[data_row_index]['OBJECT'])
        sourcename = FITSrec[data_row_index]['OBJECT'].replace('_',' ')
        self.logger.debug("add_data: OBJECT is '%s'",
                          FITSrec[data_row_index]['OBJECT'])
        FITSrec[data_row_index]['VELOCITY'] = \
                                     self.collector.sources[sourcename]['Vlsr']
        FITSrec[data_row_index]['VELDEF'] = veldef
        weather = self.logserver.get_weather(startUXtime)
        self.logger.debug("add_data: weather at %s is %s", startUXtime, weather)
        if weather:
          FITSrec[data_row_index]['TAMBIENT'] = weather[0]
          FITSrec[data_row_index]['PRESSURE'] = weather[1]
          FITSrec[data_row_index]['HUMIDITY'] = weather[2]
          FITSrec[data_row_index]['WINDSPEE'] = weather[3]
          FITSrec[data_row_index]['WINDDIRE'] = weather[4]
        else:
          self.logger.debug("add_data: weather not available for %f", startUXtime)
          FITSrec[data_row_index]['TAMBIENT'] = numpy.nan
          FITSrec[data_row_index]['PRESSURE'] = numpy.nan
          FITSrec[data_row_index]['HUMIDITY'] = numpy.nan
          FITSrec[data_row_index]['WINDSPEE'] = numpy.nan
          FITSrec[data_row_index]['WINDDIRE'] = numpy.nan
        # GBT SDFITS wants SIDEBAND
        if rx['IFmode'] == 'U':
          FITSrec[data_row_index]['SIDEBAND'] = +1
        elif rx['IFmode'] == 'L':
          FITSrec[data_row_index]['SIDEBAND'] = -1
        else:
          self.logger.error("add_data: IF mode %s is invalid; default to USB",
                            rx['IFmode'])
          FITSrec[data_row_index]['SIDEBAND'] = +1
        
        datacubeshape = FITSrec[data_row_index]['DATA'].shape
        # the frequency axis is first in FITS/FORTRAN order and last (of four)
        # in PYTHON/C order
        num_Stokes_chan = datacubeshape[3]
        
        # the sign of CDELT1 depends on the sideband of the last SSB mixer in
        # the chain
             
        # second and third data axes (coordinates)
        RA, dec = self.logserver.get_RAdec(startUXtime)
        self.logger.debug("add_data: apparent RA,dec = %f,%f", RA, dec)
        c = SkyCoord(RA, dec, unit=(u.deg, u.deg),
                     frame=FK5(equinox=Time('J'+str(fyear), scale='utc')))
        self.logger.debug("add_data: RA,dec = %f,%f", c.ra.hour, c.dec.deg)
        c2000 = c.transform_to(fk5_2000)
        self.logger.debug("add_data: precessed RA,dec = %f,%f",
                          c2000.ra.hour, c2000.dec.deg)
        FITSrec[data_row_index]['CRVAL2'] = c2000.ra.hour # hours
        FITSrec[data_row_index]['CRVAL3'] = c2000.dec.deg    # deg
        FITSrec[data_row_index]['EQUINOX'] = 2000
        # get the radial velocity of the LSR
        FITSrec[data_row_index]['VFRAME'] = \
                V_LSR(c2000.ra.hour, c2000.dec.deg, self.tel.number, starttime)
        FITSrec[data_row_index]['RVSYS'] = \
                                      FITSrec[data_row_index]['VELOCITY'] \
                                    + FITSrec[data_row_index]['VFRAME']
    
        # fourth data axis (polarization)
        FITSrec[data_row_index]['CRVAL4'] = -1
        FITSrec[data_row_index]['CDELT4'] = -1 # for I,Q,U,V (-1,-2,-3,-4)
        
        # initialize power averages
        #IFpower[data_row_index] = {}

        # fit the tones to the data using the original resolution
        bandwidth = FITSrec[data_row_index]['BANDWIDT']
        tone_offsets, tone_chnls = tone_chnl_nums(num_chan, obsfreq, bandwidth)
        self.logger.debug("add_data: tone offsets: %s", tone_offsets)
        self.logger.debug("add_data: tone channels: %s", tone_chnls)
        tone_indices = list(tone_chnls) # define the channels to be selected
        num_tones = len(tone_indices)
        offset, center_tone = math.modf(obsfreq/1e6) # tones every MHz
        
        # the objective is to fit the position of one (e.g. central) tone
        # and one stdev with all the other tones at a fixed distance from the
        # central tone and one stdev for all tones.  The individual amplitudes
        # may vary. 'multigauss.other_pars' has the fixed parameters.
        for IF in IFs:
          self.logger.debug("add_data: processing %s", IF)
          IFidx = IFs.index(IF)
          # this is the full spectrum IF power
          IFpwr = thisdata[IF+"-ps"]
          #IFpower[data_row_index][IF] = IFpwr
          # a rail is a set of evenly spaced tones
          rails = self.check_tones(thisdata, IF, threshold=30)
          self.logger.debug("add_data: %s has %d tone rails: %s",
                            IF, len(rails), rails)
          if len(rails) > 1:
            # bad tone rails present; skip this dataset
            FITSrec[data_row_index]['CYCLE'] = 0
            bad_tones.append(scan)
            self.logger.warning(
                           "add_data: scan %d subch %s row %d has extra tones",
                           scan, subch, data_row_index) 
            continue
          elif len(rails) == 0:
            # no tones; make compressed IF spectra
            newspec, newrefval, newrefpix, newdelta = \
                        reduce_spectrum_channels(IFpwr, 0, 0, 0, num_chan=1024)
            FITSrec[data_row_index]['IFSPECTR'][IFidx, 0, 0,:] = newspec
          else:
            have_tones = True
            # accept only one rail
            toneFITSrec, newspec = self.add_tone_data(thisdata,
                    collector, cfg_key, toneFITSrec, scan,
                    subch, IFidx, IFs, midnight, IFpwr)
            # newspec is the IF spectrum compressed to 1024 channels
            FITSrec[data_row_index]['IFSPECTR'][IFidx, 0, 0,:] = newspec
            
          # compute the average power as proxy for TSYS
          FITSrec[data_row_index]['TSYS'][IFidx,0,0,0] = IFpwr.mean()
        FITSrec.columns['TSYS'].unit = "count"

        # the data in dataset is keyed on scan number
        refval = FITSrec[data_row_index]['OBSFREQ']
        refpix = num_chan/2
        delta  = self.bandwidth/num_chan
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
        FITSrec[data_row_index]['DATA'][0, 0, 0,:] = I
        FITSrec[data_row_index]['DATA'][1, 0, 0,:] = Q
        FITSrec[data_row_index]['DATA'][2, 0, 0,:] = U
        FITSrec[data_row_index]['DATA'][3, 0, 0,:] = V
        FITSrec[data_row_index]['CRVAL1'] = newrefval + delta/2
        FITSrec[data_row_index]['CRPIX1'] = newrefpix
        FITSrec[data_row_index]['CDELT1'] = newdelta
        self.logger.info("add_data: finished row %d scan %d cycle %d",
                         data_row_index, FITSrec[data_row_index]['SCAN'],
                                         FITSrec[data_row_index]['CYCLE'])
        #data_row_index += 1
        # end subch loop
      # end scan loop
    if unique(bad_tones):
      self.exthead.add_comment("bad tones in scans %s" % str(unique(bad_tones)))
    if have_tones:
      return FITSrec, toneFITSrec
    else:
      return FITSrec, None
  
  def add_tone_data(self, thisdata, collector, cfg_key, toneFITSrec, scan,
                    subch, IFidx, IFs, midnight, IFpwr):
    """
    adds data to a row in the PCG TONES table
    
    Both IFs have the same subchannels and go into the same data cube along
    the polarization axis (4th axis).
    
    CYCLE is redefined in the TONES PCG table so each tone in each subchannel
    has a different CYCLE.  Analogous with the SINGLE DISH table, different
    CYCLE values indicate different center frequencies
    
    @param IFpwr : full raw spectrum for this scan, subchannel and IF
    @type  IFpwr : numpy.array
    """
    obsfreq = self.collector.wvsr_cfg[cfg_key][1]['rf_to_if_lo']*1e6 \
                           + self.collector.wvsr_cfg[cfg_key][1][subch]['sfro']
    bandwidth = self.collector.wvsr_cfg[cfg_key][1][subch]['bandwidth']
    self.logger.debug("add_tone_data: OBSFREQ = %f", obsfreq)
    self.logger.debug("add_tone_data: BANDWIDT = %f", bandwidth)
    num_chan = self.collector.equip[cfg_key][1]['Backend'].num_chan
    year, month, day = calendar_date(self.collector.year, self.collector.doy)
    starttime = self.collector.scaninfo[scan]['start']
    startUXtime = datetime_to_UnixTime(starttime)
    tone_offsets, tone_chnls = tone_chnl_nums(num_chan, obsfreq, bandwidth)
    self.logger.debug("add_tone_data: tone offsets: %s", tone_offsets)
    self.logger.debug("add_tone_data: tone channels: %s", tone_chnls)
    scans = self.collector.fft_meta.keys()
    subchannels = self.collector.wvsr_cfg[cfg_key][1]['subchannels']
    num_subchans = len(subchannels)
    subch_idx = subchannels.index(subch)
    subch_num = subch_idx + 1
    num_tones = len(tone_offsets)
    for tone in tone_chnls:
      tone_idx = list(tone_chnls).index(tone)
      subch_tone_idx = subch_idx*num_subchans + tone_idx
      tone_row_index = get_row("TONES PCG", scans,      scan=scan,
                               num_cycles=num_subchans, cycle=subch_num,
                               num_tones=num_tones,     tone=tone_idx)
      self.logger.debug("add_tone_data: doing row %d",  tone_row_index)
      toneFITSrec[tone_row_index]['SCAN'] = scan   # int
      self.logger.debug("add_tone_data: tone(%d) for %s is at %d",
                        subch_tone_idx, IFs[IFidx], tone)
      # CYCLE increments by 1 for each row in the scan
      toneFITSrec[tone_row_index]['CYCLE'] = subch_tone_idx + 1
       # all the following are the same for each tone in the subchannel
      toneFITSrec[tone_row_index]['DATE-OBS'] = \
                                           "%4d/%02d/%02d" % (year, month, day)
      toneFITSrec[tone_row_index]['UNIXtime'] = startUXtime
      toneFITSrec[tone_row_index]['TIME'] = \
                               toneFITSrec[tone_row_index]['UNIXtime']-midnight
      toneFITSrec[tone_row_index]['EXPOSURE'] = thisdata[0]['count']
      toneFITSrec[tone_row_index]['BANDWIDT'] = \
                                 self.num_tonespec_chls*self.exthead['FREQRES']
      toneFITSrec[tone_row_index]['OBSFREQ'] = obsfreq + \
                                           (tone-65536)*self.exthead['FREQRES']
      toneFITSrec[tone_row_index]['CDELT1'] = self.exthead['FREQRES']
      self.logger.debug("add_tone_data: freq step is %f",
                        toneFITSrec[tone_row_index]['CDELT1'])
      # half the number of tone band chls
      halfband = self.num_tonespec_chls/2
      toneFITSrec[tone_row_index]['CRPIX1'] = halfband

      # the frequency of the channel nearest to the tone
      nearest_chnl_freq = obsfreq + \
                             (tone-65536)*toneFITSrec[tone_row_index]['CDELT1']
      tonefreq_kHz = (nearest_chnl_freq)/1000
      tonefreq = 1000*round(tonefreq_kHz)
      # center of the 16 channel spectrum extract
      toneFITSrec[tone_row_index]['CRVAL1'] = nearest_chnl_freq
      self.logger.debug("add_tone_data: tone frequency is %f", tonefreq)
      # get the data around the tone.
      toneFITSrec[tone_row_index]['TONES'][IFidx,:] = \
                        thisdata[IFs[IFidx]+'-ps'][tone-halfband:tone+halfband]
      # tone spectrum channel number relative to the tone channel
      tone_spec_chls = numpy.arange(self.num_tonespec_chls)-halfband
      self.logger.debug("add_tone_data: tone band channel numbers: %s",
                        tone_spec_chls)
      # now fit the tone (freq in MHz)
      #   the objective is to fit the position of one (e.g. central) tone
      #   and one std with all the other tones at a fixed distance from the
      #   central tone and one std for all tones.  The individual amplitudes
      #   may vary. 'multigauss.other_pars' has the fixed parameters.
      x = (toneFITSrec[tone_row_index]['CRVAL1'] + \
                 tone_spec_chls * toneFITSrec[tone_row_index]['CDELT1'])/1e6
      y = toneFITSrec[tone_row_index]['TONES'][IFidx,:]
      self.logger.debug("add_tone_data: x = %s", x)
      self.logger.debug("add_tone_data: y = %s", y)
      if not numpy.any(y):
      # y is all zeros.  Give up
        toneFITSrec[tone_row_index]['BASELINE'][IFidx,0,0,:] = nanarray
        toneFITSrec[tone_row_index]['TONEAMP'][IFidx,0,0,:]  = nanarray
        toneFITSrec[tone_row_index]['TONEOFST'][IFidx,0,0,:] = nanarray
        toneFITSrec[tone_row_index]['TONEWIDT'][IFidx,0,0,:] = nanarray
        continue
      est_bias = numpy.append(y[:5], y[-6:]).mean()
      initial_guess = (est_bias, y[halfband], tonefreq/1e6,
                       self.exthead['FREQRES']/1e6)
      self.logger.debug("add_tone_data: initial_guess = %s", initial_guess)
      popt, pcov = curve_fit(biased_scaled_sinc, x, y, p0=(initial_guess))
      self.logger.debug("add_tone_data: pars = %s", popt)
      self.logger.debug("add_tone_data: covars = %s", pcov)
      bias, amp, offset, std = popt
      dbias, damp, doffset, dstd = numpy.sqrt(numpy.diag(pcov))
      toneFITSrec[tone_row_index]['BASELINE'][IFidx,0,0,:] = \
                                                     numpy.array([bias, dbias])
      toneFITSrec[tone_row_index]['TONEAMP'][IFidx,0,0,:] = \
                                                       numpy.array([amp, damp])
      toneFITSrec[tone_row_index]['TONEOFST'][IFidx,0,0,:] = \
                                                 numpy.array([offset, doffset])
      toneFITSrec[tone_row_index]['TONEWIDT'][IFidx,0,0,:] = \
                                                       numpy.array([std, dstd])
      # remove the tones from the IF power data
      #    map the frequencies in 'x' into the original spectrum channels
      #    tone is the frequency of the tone
      IFpwr[tone + tone_spec_chls] -= \
                                  biased_scaled_sinc(x, bias, amp, offset, std)
      # save smoothed IF power spectra
      newspec, newrefval, newrefpix, newdelta = \
                        reduce_spectrum_channels(IFpwr, 0, 0, 0, num_chan=1024)
    return toneFITSrec, newspec
    
  def get_dwell_times(self, activity_dir):
    """
    returns the times when the antenna is tracking a positionin the sky
    
    That is, it is not slewing.
    
    @param activity_dir : where the .oe file is located
    @type  activity_dir : str
    """
    # get the .oe file contents
    files = glob.glob(activity_dir+"*.oe")
    files.sort()
    if files:
      fname = files[-1]
      self.logger.debug("get_dwell_times: %s", basename(fname))
      f = open(fname)
      lines = f.readlines()
      f.close()
    else:
      print "No .oe file in", activity_dir
      raise RuntimeError("no .oe files")
    # parse the .oe file contents
    self.oe_source = {}
    self.oe_start = {}
    self.oe_end = {}
    for line in lines[8:]:
      parts = line.strip().split()
      scan = int(parts[0])
      mylogger.debug("get_dwell_times: scan %d line has %d parts", scan, len(parts))
      # later files have the azimuth and elevation at the end of a scan.
      if len(parts) == 8 or len(parts) == 10:
        self.oe_source[scan] = parts[1]
        self.oe_start[scan] = parts[5]
        self.oe_end[scan] = parts[6]
      elif len(parts) == 9 or len(parts) == 11:
        # source name has space
        self.oe_source[scan] = "_".join(parts[1:3])
        self.oe_start[scan] = parts[6]
        self.oe_end[scan] = parts[6]
      else:
        raise RuntimeError("OE file scan %d line has %d parts",
                           scan, len(parts))
    return True
  
  def check_recording_times(self, activity_dir):
    """
    check that recording times and antenna dwell times agree
    
    @param activity_dir : where the .oe file is located
    @type  activity_dir : str
    """
    self.get_dwell_times(activity_dir)
    scans = self.collector.scaninfo.keys()
    scans.sort()
    tol = datetime.timedelta(seconds=6)
    tablekeys = self.tables.keys()
    tablekeys.sort()
    # different tables have different number of rows but the same scans
    for key in tablekeys:
      table = self.tables[key]
      self.logger.debug("check_recording_times: table %s", key)
      goodscans = []
      badscans = []
      for scan in scans:
        rows = numpy.where(table.data['SCAN'] == scan)[0] # where returns tuple
        if len(rows):
          # there are two rows for each cycle
          obsdate = table.data['DATE-OBS'][rows[0]]
          # these are datetime objects
          on_start = datetime.datetime.strptime(obsdate+' '+self.oe_start[scan],
                                              "%Y/%m/%d %H:%M:%S")
          on_end   = datetime.datetime.strptime(obsdate+' '+self.oe_end[scan],
                                              "%Y/%m/%d %H:%M:%S")
          try:
            rec_start = self.collector.scaninfo[scan]['start']
            rec_end   = self.collector.scaninfo[scan]['end']
          except KeyError:
            continue
          if rec_start >= on_start-tol and rec_end <= on_end+tol:
            # scan accepted
            goodscans.append(scan)
            self.logger.debug(
                         "check_recording_times: scan %d times are compatible",
                         scan)
          else:
            # insufficient time overlap
            badscans.append(scan)
            self.logger.warning(
                           "check_recording_times: scan %d times incompatible",
                           scan)
            # set CYCLE to 0 in each table
            table.data['CYCLE'][rows] = 0
        else:
          # no rows with this scan number
          continue
      table.header.add_comment("wrong recording time in scans %s"
                                                       % str(unique(badscans)))
      self.logger.info("check_recording_times: good scans: %s", goodscans)
  
  def check_tones(self, data, IF, threshold=30):
    """
    find scans with bad tones in them
    
    This works with FFT files produced with the WVSR postproc Stokes program
    
    This is quite ad hoc and will be affected by things such as the number of
    raw spectrum channels.
    
    @param data : table of FFT data from a WVSR
    @type  data : numpy structured array
    """
    indices = data[IF+'-ps'].argsort()[::-1]
    result = numpy.histogram(data['freq'][indices[:500]]/1e6 % 1, bins=100)
    rails = numpy.where(result[0] > threshold)[0]
    self.logger.info("check_tones: %s has tone rail(s) at %s",
                       IF, result[1][rails])
    return rails
    
if __name__ == "__main__":
  examples = """
Examples
========
  run WVSR2SDFITS.py --console_loglevel=debug --dss=14 \
                          --activity=EGG0 --date=2016/237
"""  
  p = initiate_option_parser(__doc__, examples)
  p.usage='WVSR2SDFITS.py [options]'
  p.add_argument('-a', '--activity',
               dest = 'activity',
               type = str,
               default = None,
               help = "Project code")
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
  args = p.parse_args()
  
  mylogger = logging.getLogger()
  init_logging(mylogger,
                 loglevel = get_loglevel(args.file_loglevel),
                 consolevel = get_loglevel(args.console_loglevel),
                 logname = args.logpath+"WVSR2SDFITS.log")
                 
  mylogger.critical(" WVSR2SDFITS started")
  mylogger.debug("WVSR2SDITS args: %s", args)
  yearstr, doystr = args.date.split("/")
  year = int(yearstr)
  doy = int(doystr)
  project = activity_project(args.activity)
  # note that this does not handle recording sessions with multiple antennas
  obsdir, realdir, activity_dir, project_dir, datadir, wvsrdir, fftdir = \
                           get_session_dirs(args.activity, args.dss, year, doy)

  # do essential directories exists
  if not exists(activity_dir):
    print("There are no observing files for this session")
    raise RuntimeError("There are no observing files for this session")
  if not exists(wvsrdir):
    print("There is no WVSR data")
    raise RuntimeError("There is no WVSR data")
  if not exists(fftdir):
    print("Post-processing to FFTs was not done")
    raise RuntimeError("Post-processing to FFTs was not done")
  # get data for all the sources and verifiers used by the project
  activity_root = '/'.join(activity_dir.split('/')[:7])+'/'
  sourcedata = get_all_source_data(activity_root)
  
  # get the start and end time of the session
  timesfiles = glob.glob(activity_dir+"times-*")
  mylogger.debug("WVSR2SDITS found %s", timesfiles)
  if len(timesfiles) < 1:
    print("WVSR2SDITS: no times file in %s" % obsdir)
    raise RuntimeError("WVSR2SDITS: no times file in %s" % obsdir)
  elif len(timesfiles) > 1:
    print("WVSR2SDITS: can only handle one timesfile; %s has %d" \
                       % (activity_dir, len(timesfiles)))
    raise RuntimeError("WVSR2SDITS: can only handle one timesfile; %s has %d" \
                       % (activity_dir, len(timesfiles)))
  starttime = get_start_time(timesfiles[0])
  endtime = get_end_time(timesfiles[0])

  # get a manager for the NMC log for this session.
  NMClogmanager = NMClogManager(station=args.dss,
                                year=year, DOY=doy, time=(starttime, endtime),
                                use_portal=False)
  mylogger.debug("WVSR2SDITS NMC metadata available: %s",
                 NMClogmanager.metadata.keys())
  # the log server provides timed information from the NMC log
  NMClogserver = NMC_log_server(args.dss, year, doy)
  
  # get a metadata manager for the WVSR log for this session
  try:
    collector = WVSRmetadataCollector(args.activity, args.dss,
                                      year, doy, starttime)
  except RuntimeError as e:
    sys.exit("".join(traceback.format_exception(*sys.exc_info())))
    
  collector.sources = sourcedata
  mylogger.debug("__init__: equipment: %s", collector.equip)
  
  # add Backend to the standard equipment
  #   at this point, the collector only knows about backend equipment so
  #   collector.equip.keys(), collector.wvsrnames.keys() and
  #   collector.wvsr_cfg.keys() are all the same.
  config = {}
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
          print("WVSR2SDFITS can only handle one antenna")
          raise RuntimeError("WVSR2SDFITS can only handle one antenna")
      config[key] = station_configuration(None, project, dss, year, doy,
                                   starttime, rxband[key], collector=collector)
    else:
      print("single IF case not yet coded")
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
  ff.check_recording_times(activity_dir) 
  
  hdulist = pyfits.HDUList([ff.prihdu]+ff.tables.values())
  parts = realdir.split('/')
  parts.remove('projects')
  parts[3] = 'RA_data'; parts[4] = 'FITS'
  fitspath = "/".join(parts)
  mkdir_if_needed(fitspath)
  fname = fitspath + "WVSR" + "_%4d-%03d-%s.fits" % (year, doy, starttime)
  mylogger.info("WVSR2SDFITS writing %s", fname)
  hdulist.writeto(fname, overwrite=True)
  print("WVSR2SDFITS wrote %s" % fname)
  mylogger.critical(" WVSR2SDFITS ended")
