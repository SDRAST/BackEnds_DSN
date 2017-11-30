"""
gets data from FFT files and metadata from scripts

Significant variables::
  Directories::
    projdir - path to the project link in the DSAO/Science directory
    obsdir  - path to an observing session
    wvsrdir - path to the WVSR data and log files
"""
import IPython
IPython.version_info = IPython.release.version.split('.')
IPython.version_info.append('')

import glob
import logging
import numpy
import cPickle as pickle
import sys
import time

from math import ceil, log
from os import chdir, getcwd, symlink, makedirs
from os.path import basename, exists
from sys import stdout

from Automation import get_project_code, get_real_project
from Automation.sources import get_all_source_data
from Automation.WVSR import get_Naudet_FFT_dir, make_datadir_name
from Automation.WVSR import parse_scan_files, parse_WVSR_FFT_logs
from Data_Reduction import reduce_spectrum_channels
from Data_Reduction.DSN.WVSR.SpecData import get_channel_IDs, get_FFT_data
from Data_Reduction.DSN.WVSR.SpecData import read_FFT_file
from Math.multigauss import MMgauss
from support.logs import init_logging, get_loglevel

from Automation.tones import fit_tones, tone_chnl_nums, get_tone_data_chans

logger = logging.getLogger(__name__)

class Something(object):
  def __init__self):
    pass

  def plot_Tsys(self):
    """
    """
    # Output the run of system temperatures as a diagnostic:
    if self.dataset == {}:
      self.get_datasets()
    ds0 = self.dataset[0]
    fig = figure()
    for dskey in self.dataset.keys():
      ds = self.dataset[dskey]
      self.logger.debug("plot_Tsys: dataset %d scan keys: %s", dskey, ds.scan_keys)
      for key in ds.scan_keys:
        index = ds.scan_keys.index(key)
        plottimes = UnixTime_to_MPL(ds.header['time'][index])
        color_idx = 0
        for pol in [0,1]:
          for beam in [0,1]:
            tsys = ds.header['TSYS'][index][pol,:len(plottimes),beam]
            plot_date(plottimes[:len(tsys)], tsys, linestyle='-', 
                      color=plotcolors[color_idx], marker=plotsymbols[dskey])
            if key == ds.scan_keys[0]:
              plot_date(plottimes[0], tsys[0],
                        color=plotcolors[color_idx],
                        marker=plotsymbols[dskey],
                        label=clean_TeX(basename(ds.file))+", Pol "+str(pol) \
                        +" Beam "+str(beam))
            else:
              plot_date(plottimes[0], tsys[0], color=plotcolors[color_idx],
                        marker=plotsymbols[dskey])
            color_idx += 1
    ylabel(r"T$_{sys}$ (K)")
    legend(loc='best', fontsize='xx-small', numpoints=1)
    grid()
    titlestr = clean_TeX(str(ds0.year)+"/"+str(ds0.doy))
    title(titlestr)
    fig.autofmt_xdate()

if __name__ == "__main__":
  
  import argparse
  parser = argparse.ArgumentParser(
                      description=
                         'Convert WVSR FFT files and metadata to Python dict.')
  parser.add_argument('-d', '--doy',
                      dest='doy',
                      type=int,
                      default=None,
                      help='day of year')
  parser.add_argument('-D', '--dss',
                      dest='dss',
                      type=int,
                      default=14,
                      help='DSN station')
  parser.add_argument('--debug',
                      dest = 'loglevel',
                      type = str,
                      default = 'warning',
                      help = "Set debugging level. Default: warning")
  parser.add_argument('-p', '--project',
                      dest='project',
                      type=str,
                      default='AUTO_EGG',
                      help='observing project ID')
  parser.add_argument('--pickle',
                      dest='pkl',
                      action='store_true',
                      default=True,
                      help='save data and metadata')
  parser.add_argument('-y', '--year',
                      dest='year',
                      type=int,
                      default=2016,
                      help='year of observation')
  opts  = parser.parse_args()

  logging.basicConfig(level=opts.loglevel)
  mylogger = logging.getLogger()
  mylogger.setLevel(get_loglevel(opts.loglevel))
  logger = mylogger
  
  if opts.doy:
    pass
  else:
    logger.error(" doy is a required argument")
    sys.exit()
    
  dss = "dss%02d" % opts.dss
  dss_name = "DSS-%02d" % opts.dss
  
  # get data for all the sources and verifiers used by the project
  projdir = "/usr/local/projects/DSAO/Science/"+opts.project+"/"
  sourcedata, dt = get_all_source_data(projdir, opts.year, opts.doy)
  #                -------------------
  logger.debug(" date/time of Doppler calculation is %s", dt.ctime())
  
  obs_suffix = "%s/%4d/%03d/" % (dss, opts.year, opts.doy)
  obsdir = projdir+obs_suffix
  mylogger.debug(" changing to %s", obsdir)
  fftdir = get_Naudet_FFT_dir(obsdir)
  
  # return here when done
  progdir = getcwd() 
  # go to the observation directory
  chdir(obsdir)
  
  # get the scans recorded
  scaninfo = parse_scan_files(obsdir)
  #          ----------------
  scankeys = scaninfo.keys()
  scankeys.sort()
  
  # define data directory for data structure files (load or dump)
  project = get_real_project(opts.project)
  datadir = "/usr/local/project_data/"+project+"/"+obs_suffix
  if exists(datadir) == False:
    makedirs(datadir, mode=0775)
    mylogger.info(" created %s", datadir)
  
  # now we have enough information to start populating a session
  # data structure
  data = {"ORIGIN": "DSAO",
          "DATE": dt.strftime("%Y/%m/%d"),
          "PROJID": project,
          "TELESCOP": dss_name}

  # get IF channel IDs
  chan_keys = get_channel_IDs(scaninfo)
  channames = []
  for chan in chan_keys:
    channames.append(chan)

  # find the recorders used
  wvsrdir = "/data2/" + make_datadir_name(obsdir)
  wvsrlogs = glob.glob(wvsrdir+"/WVSR*")
  extensions = []
  wvsrnames = {}
  for f in wvsrlogs:
    logname = basename(f)
    wvsrnames[logname[:5].lower()] = logname
  extnum = 0
  for wvsr in wvsrnames:
    data["extension "+str(extnum)] = {'EXTNAME': wvsr}
    extnum += 1
  
  # parse the WVSR FFT logs for metadata
  self.wvsr_cfg = parse_WVSR_FFT_logs(wvsrlogs)
  #                    -------------------
  
  # add the metadata to the data structure
  data['OBSERVER'] = user
  data['BACKEND'] = wvsr_id
  
  # make extension-like dicts for each digital IF channel
  wvsrkeys = wvsrnames.keys()
  extrow = {} # each extension has a row counter
  metadata = {}
  for wvsr in wvsrkeys:
    # if there are multiple WVSRs; usually only one
    wvsrindex = wvsrkeys.index(wvsr)
    extension = "extension "+str(wvsrindex) # the name of the extension
    metadata[extension] = {}
    # order of spectra
    metadata[extension]["labels"] = ["I",  "Q",  "U",       "V",
                                     "IF1","IF2","IF1-phase","IF2-phase"]
    # each WVSR typically has two digital IF channels with orthogonal pols.
    polkeys = pol.keys()
    polkeys.sort()
    for key in polkeys:
      # there is metadata for each polarization into the WVSR
      metakeys = pol[key].keys()
      metakeys.sort()
      for mkey in metakeys:
        # extract the data for each digital IF channel for the extension head
        if mkey[:7] == 'chan_id':
          mindex = channames.index(mkey)
          # these data are for a specified digital IF channel
          metadata[extension][mindex] = {}
          metadata[extension][mindex]['expid'] = pol[key]['expid']
          #data[extension][row]['pol'] = pol[key]['pol']
          metadata[extension][mindex]['OBSFREQ'] = pol[key]['rf_to_if_lo']*1e6 \
                                                 + pol[key][mkey]['sfro']
          metadata[extension][mindex]['BANDWIDT'] = pol[key][mkey]['bandwidth']
  
  # now we can build up the scans from scaninfo.  Note that what are called
  # scans here are really CYCLES.  We could merge the extensions into one and
  # have separate rows for each scan with two CYCLES each.
  for wvsr in wvsrkeys:
    wvsrindex = wvsrkeys.index(wvsr)
    extension = "extension "+str(wvsrindex)
    row = 0
    for scan in scankeys:
      scanmeta = scaninfo[scan]
      for chan in chan_keys:
        mindex = channames.index(chan)
        data[extension][row] = {'CYCLE': mindex}
        data[extension][row]['SCAN'] = scan
        data[extension][row]['OBSFREQ'] = metadata[extension][mindex]['OBSFREQ']
        data[extension][row]['BANDWIDT'] = metadata[extension][mindex]['BANDWIDT']
        data[extension][row]['datafile'] = scaninfo[scan][chan]
        source = scaninfo[scan]['source'].replace('_',' ')
        data[extension][row]['OBJECT'] = source
        if source[-4:] == "-ref":
          data[extension][row]['RA'] = sourcedata[source[:-4]]['RA']
          data[extension][row]['dec'] = sourcedata[source[:-4]]['dec']
          data[extension][row]['VELOCITY'] = sourcedata[source[:-4]]['Vlsr']
        else:
          data[extension][row]['RA'] = sourcedata[source]['RA']
          data[extension][row]['dec'] = sourcedata[source]['dec']
          data[extension][row]['VELOCITY'] = sourcedata[source]['Vlsr']
        data[extension][row]['start'] = scaninfo[scan]['start']
        data[extension][row]['end'] = scaninfo[scan]['end']
        row += 1    
  
    # get the data; remember that each subchannel corresponds to a CYCLE
    num_rows = len(data[extension]) - 1 # because of the extension name
    for row in range(num_rows):
      logger.debug(" getting data for row %s", row)
      # for each row
      datafile = data[extension][row]['datafile']
      subchan = data[extension][row]['CYCLE']
      scan = data[extension][row]['SCAN']
      thisdata = read_FFT_file(fftdir+datafile)
      #          -------------
      if type(thisdata) != numpy.ndarray:
        # no data returned
        logger.warning(" no data for extension %s scan %d", subchan, scan) 
        continue
      if data[extension].has_key('freq'):
        # we already have freq data from previous scan
        pass
      else:
        # first scan with freq data
        data[extension]['freq'] = thisdata['freq']
      
      # more metadata for this row
      data[extension][row]['EXPOSURE'] = thisdata['count']
      
      # see if there are tones.
      num_chans = len(thisdata['IF1-ps'])
      tone_channels = tone_chnl_nums(num_chans,
                                     data[extension][row]['OBSFREQ']/1e6,
                                     data[extension][row]['BANDWIDT']/1e6)
      logger.info(" Tone channels: %s", tone_channels)
      tone_data_chans = get_tone_data_chans(tone_channels, 16)
      data[extension][row]['tonedata'] = numpy.empty((len(tone_data_chans),9)) 
      data[extension][row]['tonedata'][:,0] = thisdata['freq'][tone_data_chans]
      data[extension][row]['tonedata'][:,1] = thisdata['I'][tone_data_chans]
      data[extension][row]['tonedata'][:,2] = thisdata['Q'][tone_data_chans]
      data[extension][row]['tonedata'][:,3] = thisdata['U'][tone_data_chans]
      data[extension][row]['tonedata'][:,4] = thisdata['V'][tone_data_chans]
      data[extension][row]['tonedata'][:,5] = thisdata['IF1-ps'][tone_data_chans]
      data[extension][row]['tonedata'][:,6] = thisdata['IF2-ps'][tone_data_chans]
      data[extension][row]['tonedata'][:,7] = thisdata['IF1-phase'][tone_data_chans]
      data[extension][row]['tonedata'][:,8] = thisdata['IF2-phase'][tone_data_chans]
      
      # fit the tones
      freqs = data[extension]['freq']/1e6
      data[extension][row]['tonemeta'] = numpy.empty((tone_channels.shape[0]+2,6)) 
      x, y, position, std, amps, center_tone, offsets = \
                            fit_tones(data, extension, row, which='I')
      data[extension][row]['tonemeta'][:,0] = \
                               numpy.append(numpy.array([position]+[std]),amps)
      tonefree = thisdata['I'] - MMgauss(freqs, position, std, *amps)
      I = reduce_spectrum_channels(tonefree, 
                                   linefreq=data[extension][row]['OBSFREQ'],
                                   bandwidth=data[extension][row]['BANDWIDT'],
                                   max_vel_width=0.05)
 
      x, y, position, std, amps, center_tone, offsets = \
                            fit_tones(data, extension, row, which='Q')
      data[extension][row]['tonemeta'][:,1] = \
                               numpy.append(numpy.array([position]+[std]),amps)
      tonefree = thisdata['Q'] - MMgauss(freqs, position, std, *amps)
      Q = reduce_spectrum_channels(tonefree, 
                                   linefreq=data[extension][row]['OBSFREQ'],
                                   bandwidth=data[extension][row]['BANDWIDT'],
                                   max_vel_width=0.05)
                                   
      x, y, position, std, amps, center_tone, offsets = \
                            fit_tones(data, extension, row, which='U')
      data[extension][row]['tonemeta'][:,2] = \
                               numpy.append(numpy.array([position]+[std]),amps)
      tonefree = thisdata['U'] - MMgauss(freqs, position, std, *amps)
      U = reduce_spectrum_channels(tonefree, 
                                   linefreq=data[extension][row]['OBSFREQ'],
                                   bandwidth=data[extension][row]['BANDWIDT'],
                                   max_vel_width=0.05)
                                   
      x, y, position, std, amps, center_tone, offsets = \
                            fit_tones(data, extension, row, which='V')
      tonefree = thisdata['V'] - MMgauss(freqs, position, std, *amps)
      data[extension][row]['tonemeta'][:,3] = \
                               numpy.append(numpy.array([position]+[std]),amps)
      V = reduce_spectrum_channels(tonefree, 
                                   linefreq=data[extension][row]['OBSFREQ'],
                                   bandwidth=data[extension][row]['BANDWIDT'],
                                   max_vel_width=0.05)
                                   
      x, y, position, std, amps, center_tone, offsets = \
                            fit_tones(data, extension, row, which='IF1')
      tonefree = thisdata['IF1-ps'] - MMgauss(freqs, position, std, *amps)
      data[extension][row]['tonemeta'][:,4] = \
                               numpy.append(numpy.array([position]+[std]),amps)
      IF1 = reduce_spectrum_channels(tonefree, 
                                   linefreq=data[extension][row]['OBSFREQ'],
                                   bandwidth=data[extension][row]['BANDWIDT'],
                                   max_vel_width=0.05)
                                   
      x, y, position, std, amps, center_tone, offsets = \
                            fit_tones(data, extension, row, which='IF2')
      tonefree = thisdata['IF2-ps'] - MMgauss(freqs, position, std, *amps)
      data[extension][row]['tonemeta'][:,5] = \
                               numpy.append(numpy.array([position]+[std]),amps)
      IF2 = reduce_spectrum_channels(tonefree, 
                                   linefreq=data[extension][row]['OBSFREQ'],
                                   bandwidth=data[extension][row]['BANDWIDT'],
                                   max_vel_width=0.05)
      
      data[extension][row]['data'] = numpy.empty((len(I),1,1,4))
      data[extension][row]['data'][:,0,0,0] = I
      data[extension][row]['data'][:,0,0,1] = Q
      data[extension][row]['data'][:,0,0,2] = U
      data[extension][row]['data'][:,0,0,3] = V
      
      # this adds the IF monitor spectra
      data[extension][row]['IF-mon']   = numpy.empty((512,2))
      IF1 = reduce_spectrum_channels(thisdata['IF1-ps'], num_chan=512)
      IF2 = reduce_spectrum_channels(thisdata['IF2-ps'], num_chan=512)
      data[extension][row]['IF-mon'][:,0] = IF1
      data[extension][row]['IF-mon'][:,1] = IF2
      logger.debug(" got data for %s row %d", extension, row)

  
  # save the data
  outfile = open(datadir+"IQUV.pkl","w")
  pickle.dump(data, outfile)
  outfile.close()
  
  # return to original directory
  chdir(progdir)
