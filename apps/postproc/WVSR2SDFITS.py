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
    make_data_axis()                   # SPECTRUM
    make_offset_columns()              # BEAMxOFF
    add_data()
"""
import IPython
IPython.version_info = IPython.release.version.split('.')
IPython.version_info.append('')

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
#from Automation.WVSR import parse_scan_files, parse_WVSR_FFT_logs
from Data_Reduction import get_obs_dirs, get_obs_session, select_data_files
from Data_Reduction import get_num_chans, reduce_spectrum_channels
from Data_Reduction.DSN.WVSR.SpecData import read_FFT_file
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
estimate_Tsys = True
num_tonespec_chls = 16

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
warnings.filterwarnings('error')

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

class FITSfile(object):
  """
  A FITS file object having primary header and binary table extensions.
  
  The header describes the type of FITS file it is and where and when it was
  created.
  
  Each extension consists of a header followed by a column-oriented table
  (columns all have the same type of data).  The cells of the table can
  contain simple values or data arrays.  The array shape is defined in the
  header.  Columns in which all cells have the same value are considered
  'virtual' may be replaced and may be replaced by a keyword, value pair
  in the extension header.

  Public attributes::
    logger - a logging.Logger obhect
    prihdu - a pyfits.PrimaryHDU object
    tables - pyfits.BinTableHDU objects
  """
  def __init__(self, tel):
    """
    Initialize an SDFITS file with no extensions.

    Call this first to initiate an SDFITS file.  It creates the primary HDU
    suitable for a binary table extension.

    Notes
    =====

    The default length of axis 4 (STOKES) is 1 for compatibility with a
    "one spectrum per row" convention.  If this is not to be followed then
    subsequently called methods (probably Observatory.FITS_init_backends)
    must change MAXIS4.

    TDIMxx is calculated from MAXIS1 through MAXISn at the end.
    """
    self.logger = logging.getLogger(logger.name+".FITSfile")
    self.tel = tel
    self.logger.debug(" creating for %s", self.tel.name)
    self.make_prihdu()
    self.add_site_data(self.prihdu.header)
    #self.logger.debug(" Initial header: %s", self.prihdu.header)
    self.tables = {}

  def make_prihdu(self):
    """
    """
    self.prihdu = pyfits.PrimaryHDU()
    self.prihdu.header['BLOCKED'] = 'T'
    self.prihdu.header['DATE'] = time.strftime("%Y/%m/%d",time.gmtime())
    self.prihdu.header['ORIGIN'] = 'FITSfile.__init__'
    
  def add_site_data(self, hdr):
    """
    Adds telescope data to header

    This may move to the file header when the INHERIT keyword is recognized.

    @param hdr : the header to be modified
    @type  hdr : pyfits Header instance
    """
    hdr['telescop'] = self.tel.name
    hdr['sitelong'] = (self.tel['longitude'], "degrees west of Greenwich")
    hdr['sitelat']  = (self.tel['latitude'],  "degrees")
    hdr['siteelev'] = (self.tel['elevation'], "meters")
    hdr['obsgeo-x'] = (self.tel['geo-x'],     "meters")
    hdr['obsgeo-y'] = (self.tel['geo-y'],     "meters")
    hdr['obsgeo-z'] = (self.tel['geo-z'],     "meters")

  def make_basic_header(self):
    """
    Starts a header with the required values
    
    This includes values that are applicable to DSN SDFITS as well as SDFITS
    """
    head  = pyfits.Header() ## tbhdu.header
    head['extname'] = ("SINGLE DISH",               "required keyword value")
    head['nmatrix'] = (1,                           "one DATA column")
    head['veldef']  = ('FREQ-OBS',                  "raw receiver frequency")
    head['TIMESYS'] = ('UTC', "DSN standard time")
    return head
    
  def make_basic_columns(self, numrecs=1):
    """
    Make the minimum set of columns needed by SDFITS

    This make the REQUIRED columns for an SDFITS binary table::
    * SCAN     - scan number
    * CYCLE    - subscan number; increments by 1 for every row
    * DATE-OBS - ISO format date and time
    * OBJECT   - source name
    * OBSMODE  - observing mode
    * SIG      - True if on-source
    * CAL      - True if noise diode is on
    * TCAL     - effective noise diode temperature
    * EXPOSURE - integration time, sec
    * TIME     - a required FITS keyword
    * BANDWIDT - spectrometer bandwidth, Hz
    * SIDEBAND - lower or upper with respect to the LO
    * RESTFREQ - frequency of spectral line in the local frame
    * OBSFREQ  - center frequency of the receiver
    * VELDEF   - definition of the reference frame
    * RVSYS    - radial velocity of source w.r.t. telescope
    * VFRAME   - radial velocity of rest frame w.r.t. telescope
    * VELOCITY - radial velocity of source w.r.t. rest frame
    * EQUINOX  - epoch for source coordinates    
    * FOFFREF1 - frequency offset for frequency switching

    To these are usually appended the various beam offsets available at DSN
    antennas.  See 'make_offset_columns'.

    @param numrecs : minimum of the number of records in each scan
    @type  numrecs : 1
    """
    # create empty column data arrays
    # create required columns.

    cols = pyfits.ColDefs([
      pyfits.Column(name='SCAN',     format='1I'),
      pyfits.Column(name='CYCLE',    format='1I'),  # not used
      pyfits.Column(name='DATE-OBS', format='16A'),
      pyfits.Column(name='OBJECT',   format='16A'),
      pyfits.Column(name='OBSMODE',  format='8A'),
      pyfits.Column(name='SIG',      format='1L'),
      pyfits.Column(name='CAL',      format='1L'),
      pyfits.Column(name='TCAL',     format='1E'),
      pyfits.Column(name='EXPOSURE', format='1E', unit='s'),
      pyfits.Column(name='TIME',     format='1E', unit='s'),
      pyfits.Column(name='BANDWIDT', format='1E', unit='Hz'),
      pyfits.Column(name='SIDEBAND', format='1A'),
      pyfits.Column(name='RESTFREQ', format='1D', unit='Hz'),
      pyfits.Column(name='OBSFREQ',  format='1D', unit='Hz')])
    # Velocity data
    cols += pyfits.ColDefs(
                [pyfits.Column(name='VELDEF',   format='8A'),
                 pyfits.Column(name='RVSYS',    format='1E', unit='m/s'),
                 pyfits.Column(name='VFRAME',   format='1E', unit='m/s'),
                 pyfits.Column(name='VELOCITY', format='1E', unit='m/s'),
                 pyfits.Column(name='EQUINOX',  format='1E')])

    # frequency switching offset
    cols += pyfits.ColDefs(
                     [pyfits.Column(name='FOFFREF1',  format='1E', unit='Hz')])
    
    cols += self.make_offset_columns()
    return cols

  def make_offset_columns(self, numrecs=1, 
                          Aoff=False,Xoff=True,Eoff=True,
                          equatorial=False,
                          galactic=False,
                          refpos=False):
    """
    beam offsets such as those used to make maps

    Notes
    =====
    
    Pointing
    ~~~~~~~~
    These are not pointing offsets used by the antenna control computer to
    align the beam with the source.
    
    Columns
    ~~~~~~~
    Unused columns, that is, all zeros, should not be included so I think this
    method needs some arguments.  Do we need the REF offsets?
    These are the offsets used by DSN antennas::
      BEAMAOFF - azimuth offset
      BEAMXOFF - cross-elevation offset
      BEAMEOFF - elevation offset
      BEAMHOFF - hour angle offset
      BEAMCOFF - cross-declination offset
      BEAMDOFF - declination offset
    On 2016 Dec 22 these were added::
      BEAMLOFF - galactic longitude offset
      BEAMBOFF - galactic latitude offset
      BEAMGOFF - galactic cross-latitude offset\end{verbatim}
    """
    # always required
    if numrecs > 1:
      cols = pyfits.ColDefs([pyfits.Column(name='BEAMXOFF',
                                           format=str(numrecs)+'E',
                                          dim="(1,1,1,1,"+str(numrecs)+",1)"),
                             pyfits.Column(name='BEAMEOFF',
                                           format=str(numrecs)+'E',
                                          dim="(1,1,1,1,"+str(numrecs)+",1)")])
    else:
      cols = pyfits.ColDefs([pyfits.Column(name='BEAMXOFF',
                                           format='E', unit='deg'),
                             pyfits.Column(name='BEAMEOFF',
                                           format='E', unit='deg')])
                                           
    # the following are traditional columns.  Use above format for dimensioned
    # columns
    if Aoff: # if needed
      cols.add_col(pyfits.Column(name='BEAMAOFF', format='E', unit='deg'))
    if equatorial: # equatorial offsets
      cols.add_col(pyfits.Column(name='BEAMHOFF', format='E', unit='deg'))
      cols.add_col(pyfits.Column(name='BEAMCOFF', format='E', unit='deg'))
      cols.add_col(pyfits.Column(name='BEAMDOFF', format='E', unit='deg'))
    if galactic:# find this code
      pass 
    if refpos: # reference position for position switching
      cols.add_col(pyfits.Column(name='REF_HOFF', format='E', unit='deg'))
      cols.add_col(pyfits.Column(name='REF_DOFF', format='E', unit='deg'))
    return cols

  def add_time_dependent_columns(self, numrecs, cols):
    """
    create columns for the time-dependent metadata
    """
    if numrecs == 1:
      time_dim = "(1,)"
    else:
      # index order is spectrum, RA, dec, pol, time, beam
      time_dim = "(1,1,1,1,"+str(numrecs)+",1)" # FORTRAN order to C reverses
      # for a single beam the last index is dropped
    cols.add_col(pyfits.Column(name='LST',
                               format=str(numrecs)+'D',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='UNIXtime',
                               format=str(numrecs)+'D',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='AZIMUTH',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='ELEVATIO',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='TAMBIENT',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='PRESSURE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='HUMIDITY',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='WINDSPEE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    cols.add_col(pyfits.Column(name='WINDDIRE',
                               format=str(numrecs)+'E',
                               dim=time_dim))
    return cols
    
  def make_data_axis(self, header, columns,
                             daxis,
                             length,
                             dtype,
                             dformat,
                             unit=None, comment=None):
    """
    create header item and columns for a data axis

    @param tbhdu : extension (table HDU)
    @type  tbhdu : pyfits.TableHDU objetc
    
    @param daxis : axis number
    @type  daxis : int
    
    @param length : length of the axis (number of 'pixels')
    @type  length : int

    @param dtype : axis data type
    @type  dtype : str

    @param dformat : "D", "E", "I"
    @type  dformat : str
    
    @param unit : unit of measurement (defaults to SI or None)
    @type  unit : str
    """
    if comment == None:
      comment = "SPECTRUM axis "+str(daxis)
    header['ctype'+str(daxis)] = (dtype, comment)
    if unit:
      newcols = pyfits.ColDefs(
            [pyfits.Column(name='CRVAL'+str(daxis), format='1'+dformat, unit=unit),
             pyfits.Column(name='CRPIX'+str(daxis), format='1I'),
             pyfits.Column(name='CDELT'+str(daxis), format='1'+dformat, unit=unit)])
    else:
      newcols = pyfits.ColDefs(
            [pyfits.Column(name='CRVAL'+str(daxis), format='1'+dformat),
             pyfits.Column(name='CRPIX'+str(daxis), format='1I'),
             pyfits.Column(name='CDELT'+str(daxis), format='1'+dformat)])
    for col in newcols:
      columns.add_col(col)
    header['maxis'+str(daxis)] = length
    self.logger.debug("make_data_axis: MAXIS%d = %d", daxis, length)
    return header, columns

class FITSfile_from_WVSR(FITSfile):
  """
  subclass for creating FITS file from WVSR data
  
  WVSR spectra may contain tones from a phase cal generator which are typically
  spaced by 1 MHz (but could be 4 MHz). This is a feature not normally found in
  astronomical spectra but can be useful in calibration.  These tones are
  extracted from the high-resolution raw FFTs and saved in a separate table.
  """
  def __init__(self, tel, estimate_Tsys=True):
    """
    Initialize a FITSfile_from_WVSR class
    
    This creates a table ('tbhead' and 'cols') for every connected backend and,
    if PCG tones were on, also for the PCG tones ('tonehead' and 'tonecols').
    """
    mylogger = logging.getLogger(logger.name+'.FITSfile_from_WVSR')
    FITSfile.__init__(self, tel)
    self.logger = mylogger
    self.estimate_Tsys = estimate_Tsys
  
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
    #self.collector.equip[key][1]['Backend'].num_chan = nchans
    config[key]['Backend'].num_chan = nchans
    # compute number of spectrum channels for science
    subchannel_names = self.collector.wvsr_cfg[cfg_key][1]['subchannels'] # 
    num_subchans = len(subchannel_names)
    anysubch = subchannel_names[0] # any subchannel
    obsfreq = self.collector.wvsr_cfg[cfg_key][1]['rf_to_if_lo']*1e6 \
                             + self.collector.wvsr_cfg[cfg_key][1][anysubch]['sfro']
    bandwidth = self.collector.wvsr_cfg['wvsr2'][1]['chan_id 1']['bandwidth']
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
    tbhead = self.make_basic_header()
    #fe_key = self.collector.equip[key][1]['FrontEnd'].keys()[0] # only one FE
    self.fe_key = config[key]['FrontEnd'].keys()[0] # only one FE
    #rx_key = self.collector.equip[key][1]['Receiver'].keys()[0] # only one receiver
    self.rx_key = config[key]['Receiver'].keys()[0] # only one receiver
    tbhead['FRONTEND'] = (config[key]['FrontEnd'][self.fe_key].name,
                          "front end ID")
    tbhead['RECEIVER'] = (config[key]['Receiver'][self.rx_key].name,
                          "receiver ID")
    # get_hardware_metadata takes care of BACKEND
    tbhead['projid']  = self.collector.project
    tbhead['observer'] = self.collector.wvsr_cfg[key]['user']
    
    # start extension header for tone data
    #    I don't trust 'tbhead' because assigment doesn't necessarily make a
    # stand-alone copy.
    tonehead = self.make_basic_header()
    # suspect extensions are ordered alphabetically by name; TONES after SINGLE
    tonehead['extname'] = ("TONES_PCG", "phase cal tone extension") # override
    #tonehead['FRONTEND'] = (self.collector.equip[key][1]['FrontEnd'][fe_key].name, "front end ID")
    tonehead['FRONTEND'] = (config[key]['FrontEnd'][self.fe_key].name, "front end ID")
    #tonehead['RECEIVER'] = (self.collector.equip[key][1]['Receiver'][rx_key].name, "receiver ID")
    tonehead['RECEIVER'] = (config[key]['Receiver'][self.rx_key].name, "receiver ID")
    tonehead['projid']  = self.collector.project
    tonehead['observer'] = self.collector.wvsr_cfg[key]['user']
    
    # add the site data (again)
    self.add_site_data(tbhead)
    
    # make the basic columns that are always needed
    cols = self.make_basic_columns()
    tonecols = self.make_tone_columns()

    # add the backend data
    BE = config[key]['Backend'].wvsrs[key]
    tbhead, cols = self.get_hardware_metadata(BE, tbhead, cols)
    tonehead, tonecols = self.get_hardware_metadata(BE, tonehead, tonecols)

    # things like LST, wind, etc.
    cols = self.add_time_dependent_columns(1, cols)
    # we have to skip the next one if there is no system temperature information
    #if self.estimate_Tsys:
    #  cols = self.add_IF_dependent_columns(1, cols)
    
    # Add columns describing the data matrix
    #   Note that refpix defaults to 0
    axis = 1; tbhead, cols = self.make_data_axis(tbhead, cols,
                                                 axis,
                                                 num_Stokes_chan,
                                                 'FREQ-OBS', 'D',
                                                 unit='Hz',
                                comment="channel frequency in telescope frame")

    axis +=1; tbhead, cols = self.make_data_axis(tbhead, cols, axis,
                                                 nlong,
                                                 'RA---GLS', 'D',
                                                 unit='deg',
                                                 comment="RA J2000.0")
    axis +=1; tbhead, cols = self.make_data_axis(tbhead, cols, axis,
                                                 nlat,
                                                 'DEC--GLS','D',
                                                 unit='deg',
                                                 comment="decl. J2000") 
    #   Stokes axis
    #     get the polarizations from the spectrometer input signals
    axis+=1; tbhead, cols = self.make_data_axis(tbhead, cols, axis,
                                                npols, 'STOKES',  'I',
                                         comment="polarization code: 1,2,3,4")
    
    # Make the data column
    fmt_multiplier = tbhead['MAXIS1']*tbhead['MAXIS2']*tbhead['MAXIS3'] \
                    *tbhead['MAXIS4']
    self.logger.debug("make_WVSR_table: format multiplier = %d", fmt_multiplier)
    dimsval = "("+str(tbhead['MAXIS1'])+"," \
                 +str(tbhead['MAXIS2'])+"," \
                 +str(tbhead['MAXIS3'])+"," \
                 +str(tbhead['MAXIS4'])+")"
    self.logger.debug("make_WVSR_table: computed scan shape: %s", dimsval)
    data_format = str(fmt_multiplier)+"E"
    self.logger.debug("make_WVSR_table: data_format = %s", data_format)
    data_col = pyfits.Column(name='SPECTRUM',
                             format=data_format, dim=dimsval)
    cols.add_col(data_col)
    
    # Now describe the tone data structure
    num_subch = len(subchannel_names)
    total_num_tones = 0
    for subch in subchannel_names: # count up the tones
      # there is one tone every MHz
      num_tones = \
            int(self.collector.wvsr_cfg[cfg_key][1][anysubch]['bandwidth']/1e6)
      total_num_tones += num_tones
            
    #self.logger.debug("make_WVSR_table: %d tones in all subchannels", num_tones)
    #   we take 256 channels centered on the hi-res channel nearest the tone
    toneaxis = 1; tonehead, tonecols = self.make_data_axis(tonehead, tonecols,
                                                           toneaxis,
                                                           num_tonespec_chls,
                                                           'FREQ-OBS', 'D',
                                                           unit='Hz',
                                comment="channel frequency in telescope frame")
    #     a Stokes spectrometer always has two input signals
    toneaxis+=1; tonetonehead, tonecols = self.make_data_axis(tonehead, tonecols, 
                                                              toneaxis,
                                                              2, 
                                                              'STOKES',  'I',
                                           comment="polarization code: -2, -1")
    # Make the tone data column (MAXIS was set by make_data_axis)
    fmt_multiplier = tonehead['MAXIS1']*tonehead['MAXIS2']
    self.logger.debug("make_WVSR_table: tone format multiplier = %d",
                      fmt_multiplier)
    dimsval = "("+str(tonehead['MAXIS1']) + "," + str(tonehead['MAXIS2'])+")"
    self.logger.debug("make_WVSR_table: computed scan shape: %s", dimsval)
    data_format = str(fmt_multiplier)+"E"
    self.logger.debug("make_WVSR_table: data_format = %s", data_format)
    tonedata_col = pyfits.Column(name='TONES',
                             format=data_format, dim=dimsval)
    tonecols.add_col(tonedata_col)
    
    # columns for tone fits
    tonecols.add_col(pyfits.Column(name='BASELINE', format='4E', dim="(2,1,1,2)"))
    tonecols.add_col(pyfits.Column(name='TONEAMP',  format='4E', dim="(2,1,1,2)"))
    tonecols.add_col(pyfits.Column(name='TONEOFST', format='4E', dim="(2,1,1,2)"))
    tonecols.add_col(pyfits.Column(name='TONEWIDT', format='4E', dim="(2,1,1,2)"))
    #                           format=str(2*num_tones/num_subch)+'E',
    #                           dim="("+ str(num_tones/num_subch)+",1,1,2)"))
                               
    cols.add_col(pyfits.Column(name='IFSPECTR', format='2048E',
                                                          dim="(1024,1,1,2)"))
    # Tone-free averaged IF power 
    # IF dependent dim is (1, 1, 1, 2)
    # IF1 is Stokes code -1 (RR) and IF2 is -2 (LL)
    cols.add_col(pyfits.Column(name='avgpower', format='2E', dim="(1,1,1,2)"))

    if self.estimate_Tsys:
      # add column for estimated system temperatures if not measured
      cols.add_col(pyfits.Column(name='TSYS', format='2E', dim="(1,1,1,2)"))

    
    # create a table extension for astronomy data
    FITSrec = pyfits.FITS_rec.from_columns(cols, nrows=numscans*num_subchans)
    tabhdu =  pyfits.BinTableHDU(data=FITSrec, header=tbhead,
                                      name="SINGLE DISH")

    # create a table extension for tone data
    toneFITSrec = pyfits.FITS_rec.from_columns(tonecols,
                                                nrows=numscans*total_num_tones)
    tonetabhdu = pyfits.BinTableHDU(data=toneFITSrec, header=tonehead,
                                    name="TONES PCG")
    
    
    # fill in the data and tone tables rows
    tabhdu, tonetabhdu = self.add_data(tabhdu, tonetabhdu, key)
    
    self.tables[key] = tabhdu
    self.tables[key+"-pcg"] = tonetabhdu
    return tabhdu, tonetabhdu
    
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

    cols = pyfits.ColDefs([
      pyfits.Column(name='SCAN',     format='1I'             ),
      pyfits.Column(name='CYCLE',    format='1I'             ),
      pyfits.Column(name='DATE-OBS', format='16A'            ),
      pyfits.Column(name='EXPOSURE', format='1E', unit='s'   ),
      pyfits.Column(name='TIME',     format='1E', unit='s'   ),
      pyfits.Column(name='UNIXtime', format='1E', unit='s'   ),
      pyfits.Column(name='BANDWIDT', format='1E', unit='Hz'  ),
      pyfits.Column(name='SIDEBAND', format='1A'             ),
      pyfits.Column(name='OBSFREQ',  format='1D', unit='Hz'  )])
      #pyfits.Column(name='AZIMUTH',  format='1E', unit='deg' ),
      #pyfits.Column(name='ELEVATIO', format='1E', unit='deg' ),
      #pyfits.Column(name='TAMBIENT', format='1E', unit='C'   ),
      #pyfits.Column(name='PRESSURE', format='1E', unit='mB'  ),
      #pyfits.Column(name='HUMIDITY', format='1E', unit='%'   ),
      #pyfits.Column(name='WINDSPEE', format='1E', unit='km/h'),
      #pyfits.Column(name='WINDDIRE', format='1E', unit='deg' )])
      
    # frequency switching offset
    cols += pyfits.ColDefs(
                     [pyfits.Column(name='FOFFREF1',  format='1E', unit='Hz')])
    return cols

  def get_hardware_metadata(self, BE, hdr, cols):
    """
    Initializes columns for the backends.

    Creates SDFITS columns for the backends.  Each distinct backend gets its
    own extension so BACKEND doesn't need a column.
    
    """
    hdr['backend'] = BE.name
    hdr['maxis1'] =  (BE.fft_cfg['n_freqs'], "length of DATA axis 1")
    subchannel = BE.IF[1].subchannel['chan_id 1'] # any one will do
    hdr['freqres'] =  float(subchannel['bandwidth'])/subchannel['nchans']
    return hdr, cols
  
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
    numscans = len(scans) # scans observed
    # both IFs have the same subchannels
    subchannels = self.collector.wvsr_cfg[cfg_key][1]['subchannels']
    # create frame to return antenna RA and dec to J2000
    fk5_2000 = FK5(equinox=Time(2000, format='jyear', scale='utc'))
    self.IFpower = {}
    blank_array = numpy.array(256*[numpy.nan])
    
    data_row_index = 0
    tone_row_index = 0
    for scan in self.collector.scankeys:
      # dataset header is keyed on the index of the scan in the set of scans
      scan_index = scans.index(scan)
      # use date of first record; see doc string for explanation of extra index
      year, month, day = calendar_date(self.collector.year, self.collector.doy)
      # UNIX time at midnight
      midnight = time.mktime(dateutil.parser.parse(
                          tabhdu.data[data_row_index]['DATE-OBS']).timetuple())
      subch_tone_idx = 0 # collect tone data from both subchannels
      for subch in subchannels:
        self.logger.debug("add_data: processing scan %d %s data row %d",
                        scan, subch, data_row_index)
        self.logger.debug("add_data: processing scan %d %s tone row %d",
                        scan, subch, tone_row_index)
        sub_idx = subchannels.index(subch)
        tabhdu.data[data_row_index]['SCAN'] = scan   # int
        tabhdu.data[data_row_index]['DATE-OBS'] = \
                                           "%4d/%02d/%02d" % (year, month, day)
        # each subchannel has its own cycle
        tabhdu.data[data_row_index]['CYCLE'] = sub_idx+1
        self.logger.debug("add_data: CYCLE = %d",
                          tabhdu.data[data_row_index]['CYCLE'])
        datafile = self.collector.scaninfo[scan]['subch '+str(sub_idx+1)]
        # this returns a structured array with 131072 spectrum channels
        thisdata = read_FFT_file(fftdir+datafile)
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
        # same frequency and bandwidth for both IFs
        tabhdu.data[data_row_index]['BANDWIDT'] = \
                       self.collector.wvsr_cfg['wvsr2'][1]['chan_id 1']['bandwidth']
        obsfreq = self.collector.wvsr_cfg[cfg_key][1]['rf_to_if_lo']*1e6 \
                                + self.collector.wvsr_cfg[cfg_key][1][subch]['sfro']
        tabhdu.data[data_row_index]['OBSFREQ'] = obsfreq
        tabhdu.data[data_row_index]['RESTFREQ'] = obsfreq # is this always true?
        sourcename = tabhdu.data[data_row_index]['OBJECT'].replace('_',' ')
        tabhdu.data[data_row_index]['VELOCITY'] = \
                                          self.collector.sources[sourcename]['Vlsr']
        tabhdu.data[data_row_index]['VELDEF'] = veldef
        weather = self.logserver[cfg_key].get_weather(startUXtime)
        self.logger.debug("add_data: weather at %s is %s", startUXtime, weather)
        tabhdu.data[data_row_index]['TAMBIENT'] = weather[0]
        tabhdu.data[data_row_index]['PRESSURE'] = weather[1]
        tabhdu.data[data_row_index]['HUMIDITY'] = weather[2]
        tabhdu.data[data_row_index]['WINDSPEE'] = weather[3]
        tabhdu.data[data_row_index]['WINDDIRE'] = weather[4]
      
        # for first data axis (frequency)
        num_chan = self.collector.equip[cfg_key][1]['Backend'].num_chan
        tabhdu.data[data_row_index]['CRVAL1'] = \
                                         tabhdu.data[data_row_index]['OBSFREQ']
        datacubeshape = tabhdu.data[data_row_index]['SPECTRUM'].shape
        num_Stokes_chan = datacubeshape[3]
        tabhdu.data[data_row_index]['CDELT1'] = \
                        tabhdu.data[data_row_index]['BANDWIDT']/num_Stokes_chan
        tabhdu.data[data_row_index]['CRPIX1'] = num_Stokes_chan/2 # at middle
        # GBT SDFITS wants SIDEBAND
        rx_key = self.collector.equip[cfg_key][1]['Receiver'].keys()[0] # only one
        rx = self.collector.equip[cfg_key][1]['Receiver'][rx_key]
        self.logger.debug("add_data: receiver is %s", rx)
        tabhdu.data[data_row_index]['SIDEBAND'] = rx['IFmode']
        tabhdu.data[data_row_index]['CDELT1'] = \
                        tabhdu.data[data_row_index]['BANDWIDT']/num_Stokes_chan
      
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
        # fit the tones to the data.
        bandwidth = tabhdu.data[data_row_index]['BANDWIDT']
        tone_offsets, tone_chnls = tone_chnl_nums(num_chan, obsfreq, bandwidth)
        self.logger.debug("add_data: tone offsets: %s", tone_offsets)
        self.logger.debug("add_data: tone channels: %s", tone_chnls)
        offset, center_tone = math.modf(obsfreq/1e6) # tones every MHz
        # the objective is to fit the position of one (e.g. central) tone
        # and one std with all the other tones at a fixed distance from the
        # central tone and one std for all tones.  The individual amplitudes
        # may vary. 'multigauss.other_pars' has the fixed parameters.
        for IF in ['IF1', 'IF2']:
          IFidx = int(IF[-1])-1
          IFpwr = self.IFpower[data_row_index][IF]
          tone_indices = list(tone_chnls) # define the channels to be selected
          num_tones = len(tone_indices)
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
            #tonetabhdu.data[tone_row_index]['CRVAL1'] = obsfreq
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
            tonetabhdu.data[tone_row_index]['BASELINE'] = \
                                                     numpy.array([bias, dbias])
            tonetabhdu.data[tone_row_index]['TONEAMP'] = \
                                                       numpy.array([amp, damp])
            tonetabhdu.data[tone_row_index]['TONEOFST'] = \
                                                 numpy.array([offset, doffset])
            tonetabhdu.data[tone_row_index]['TONEWIDT'] = \
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
          tabhdu.data[data_row_index]['IFSPECTR'][IFidx, 0, 0,:] = \
                                 reduce_spectrum_channels(IFpwr, num_chan=1024)
          # compute the average power
          mean = IFpwr.mean()
          tabhdu.data[data_row_index]['avgpower'][IFidx,0,0,0] = mean
          # end of IF loop
        # the data in dataset in keyed on scan number
        self.logger.debug("add_data: loading SPECTRUM")
        tabhdu.data[data_row_index]['SPECTRUM'][0, 0, 0,:] = \
              reduce_spectrum_channels(thisdata['I'], num_chan=num_Stokes_chan)
        tabhdu.data[data_row_index]['SPECTRUM'][1, 0, 0,:] = \
              reduce_spectrum_channels(thisdata['Q'], num_chan=num_Stokes_chan)
        tabhdu.data[data_row_index]['SPECTRUM'][2, 0, 0,:] = \
              reduce_spectrum_channels(thisdata['U'], num_chan=num_Stokes_chan)
        tabhdu.data[data_row_index]['SPECTRUM'][3, 0, 0,:] = \
              reduce_spectrum_channels(thisdata['V'], num_chan=num_Stokes_chan)
        data_row_index += 1
        """
          # set up the fit
          y = self.IFpower[data_row_index][IF]
          initial_std = 2.0
          popt, pcov = curve_fit(MMgauss, x, y,
                                 p0=(-offset, initial_std)+tone_amps)
          position,std = popt[:2]
          amps = popt[2:]
          self.logger.debug("add_data: %s position=%f, std=%f",
                                                             IF, position, std)
          self.logger.debug("add_data: %s tone amplitudes: %s", IF, amps)
          # add the fit results to table
          tabhdu.data[data_row_index]['TONEOFST'][IFidx,0,0,0] = position
          tabhdu.data[data_row_index]['TONEWIDT'][IFidx,0,0,0] = std
          tabhdu.data[data_row_index]['TONEAMPS'][IFidx,0,0,:] = amps
          
        # extact the tone data
        for tone in tone_chnls:
          tonetabhdu.data[tone_row_index]['SCAN'] = scan   # int
          self.logger.debug("add_data: processing tone(%d) is %d",
                            tone_idx, tone)
          # CYCLE increments by 1 for each row in the scan
          tonetabhdu.data[tone_row_index]['CYCLE'] = tone_idx + 1
          tonetabhdu.data[tone_row_index]['DATE-OBS'] = \
                                           "%4d/%02d/%02d" % (year, month, day)
          tonetabhdu.data[tone_row_index]['UNIXtime'] = startUXtime
          tonetabhdu.data[tone_row_index]['TIME'] = \
                           tonetabhdu.data[tone_row_index]['UNIXtime']-midnight
          tonetabhdu.data[tone_row_index]['EXPOSURE'] = thisdata[0]['count']
          tonetabhdu.data[tone_row_index]['BANDWIDT'] = \
                                                   256*tabhdu.header['FREQRES']
          tonetabhdu.data[tone_row_index]['TAMBIENT'] = weather[0]
          tonetabhdu.data[tone_row_index]['PRESSURE'] = weather[1]
          tonetabhdu.data[tone_row_index]['HUMIDITY'] = weather[2]
          tonetabhdu.data[tone_row_index]['WINDSPEE'] = weather[3]
          tonetabhdu.data[tone_row_index]['WINDDIRE'] = weather[4]
          tonetabhdu.data[tone_row_index]['OBSFREQ'] = obsfreq + \
                                          (tone-65536)*tabhdu.header['FREQRES']
          tonetabhdu.data[tone_row_index]['CDELT1'] = tabhdu.header['FREQRES']
          tonetabhdu.data[tone_row_index]['TONES'][0,:] = \
                                          thisdata['IF1-ps'][tone-128:tone+128]
          tonetabhdu.data[tone_row_index]['TONES'][1,:] = \
                                          thisdata['IF2-ps'][tone-128:tone+128]
          tone_idx += 1
          tone_row_index += 1
        # after all the tones have been removed, compute the average
        for IF in ['IF1', 'IF2']:
          #   first take out the tones
          IFidx = int(IF[-1])-1
          offset = tabhdu.data[data_row_index]['TONEOFST'][IFidx,0,0,0]
          std = tabhdu.data[data_row_index]['TONEWIDT'][IFidx,0,0,0]
          amps = list(tabhdu.data[data_row_index]['TONEAMPS'][IFidx,0,0,:])
          self.logger.debug("add_data: std, offset: %f, %f", std, offset)
          self.logger.debug("add_data: TONEAMPS: %s type %s", amps, type(amps))
          IFpwr = self.IFpower[data_row_index][IF] - \
                                                 MMgauss(x, offset, std, *amps)
          # save smoothed IF power spectra
          tabhdu.data[data_row_index]['IFSPECTR'][IFidx, 0, 0,:] = \
                                 reduce_spectrum_channels(IFpwr, num_chan=1024)
          # compute the average power
          mean = IFpwr.mean()
          tabhdu.data[data_row_index]['avgpower'][IFidx,0,0,0] = mean
        """
    # when all the scans are done, optionally estimate system temperatures
    if estimate_Tsys:
      tabhdu = self.fit_mean_power_to_airmass(config, cfg_key, tabhdu)

    return tabhdu, tonetabhdu

  def fit_mean_power_to_airmass(self, config, cfg_key, tabhdu):
    """
    fits the mean power data vs airmass to the radiative transfer equation
    """
    def opacity_fitting(x, a, tau):
      x_rad = numpy.deg2rad(x)
      x_sec = 1/numpy.sin(x_rad)
      return a + tau*x_sec
      
    # get the system temperature
    fe = config[cfg_key]['FrontEnd'][self.fe_key]
    
    # subch 1 is in the even numbered rows:  
    # subch 2 is in the odd numbered rows    
    # the first two rows in a set of four are on source 
    # the second two rows in a set of four are off source
    # subch1 on is  [0::4] (subch 0, pos 0)
    # subch2 on is  [1::4] (subch 1, pos 0)
    # subch1 off is [2::4] (subch 0, pos 1)
    # subch2 off is [3::4] (subch 1, pos 1)
    for IFidx in [0,1]: # IF1 and IF2
      pol = ["R","L"][IFidx]
      Tvac = fe.Tsys_vacuum(pol=pol)
      msg = "estimated zenith IF%d Tsys in vacuum= %6.2f" % (IFidx,Tvac)
      tabhdu.header.add_history(msg)
      for subch in [0,1]: # 1, 2
        mean_power = tabhdu.data['avgpower'][:,IFidx,0,0,0][subch::2]
        elevation  = tabhdu.data['ELEVATIO'][subch::2]
        # get an array of the 'nan' values in the elevations
        good_elev = ~numpy.isnan(elevation)
        # extract the good values
        elv = elevation[good_elev]
        pwr = mean_power[good_elev]
        # fit the data
        popt, pcov = curve_fit(opacity_fitting, elv, pwr, p0=[0, 0])
        intercept, slope = popt[0], popt[1]
        self.logger.debug("fit_mean_power_to_airmass: intercept, slope: %f, %f",
                        intercept, slope)
        msg="IF%d, subch%d gain=%9.3e counts, gain_slope=%9.3e counts/airmass"\
              % (IFidx+1, subch+1, intercept, slope)
        tabhdu.header.add_history(msg)
        gain = Tvac/intercept
        K_per_am = gain*slope
        tabhdu.data['TSYS'][:,IFidx,0,0,0][subch::2] = gain * mean_power
    return tabhdu
    
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
  p.add_argument('-D', '--DSS',
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
                 logname = args.logpath+"SAO2SDFITS.log")
  mylogger.debug("WVSR2SDITS args: %s", args)
  
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
  mylogger.debug(" found %s", timesfiles)
  if len(timesfiles) < 1:
    raise RuntimeError("WVSR2SDITS: no times file is %s" % obsdir)
  elif len(timesfiles) > 1:
    raise RuntimeError("WVSR2SDITS: can only handle one timesfile; %s has %d" \
                       % len(timesfiles))
  starttime = get_start_time(timesfiles[0])

  # get a metadata manager for the WVSR log for this session
  collector = WVSRmetadataCollector(args.project, args.dss, year, doy)
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
                                          rxband[key])
      # get a manager for the NMC log for this session.
      NMClogmanager[key] = NMClogManager(station=dss, year=year, DOY=doy,
                                starttime=starttime, use_portal=False)
      mylogger.debug("WVSR2SDITS NMC metadata available: %s",
                        NMClogmanager[wvsr].metadata.keys())
      NMClogserver[key] = NMC_log_server(args.project, dss, year, doy)
    else:
      raise RuntimeError("single IF case not yet coded")
  # Now we can start a FITSfile
  tel = config[config.keys()[0]]['Telescope'][args.dss]
  ff = FITSfile_from_WVSR(tel, estimate_Tsys)
  
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
  parts = datadir.split('/')
  parts[3] = 'RA_data'; parts[4] = 'FITS'
  fitspath = "/".join(parts)
  mkdir_if_needed(fitspath)
  fname = fitspath + "WVSR" + "_%4d-%03d-%s.fits" % (year, doy, starttime)
  hdulist.writeto(fname, clobber=True)

