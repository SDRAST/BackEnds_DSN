"""
"""
import logging

from MonitorControl import IF, Spectrum
from MonitorControl.BackEnds import Backend, get_freq_array
#from MonitorControl.BackEnds.DSN.helpers import WVSRmetadataCollector
from Data_Reduction.DSN.WVSR.collector import WVSRmetadataCollector

logger = logging.getLogger(__name__)

class WVSRbackend(Backend):
  """
  """
  def __init__(self, name, collector, inputs=None, output_names=None):
    """
    initiate a backend consisting of WVSRs
    
    @param collector : gets metadata from WVSR log files
    @type  collector : WVSRmetadataCollector object
    
    @param inputs : signal sources
    @type  inputs : Port objects
    
    @param output_names : names of output ports (not required)
    """
    mylogger = logging.getLogger(logger.name+".WVSRbackend")
    Backend.__init__(self, name, inputs=inputs, output_names=output_names)
    self.collector = collector
    self.logger = mylogger
    self.logger.info("__init__: %s input channels: %s", self, self.inputs)
    self.logger.info("__init__: %s output channels: %s", self, self.outputs)
    self.wvsrs = {}
    for name in collector.wvsrnames:
      self.wvsrs[name] = WVSRbackend.WVSR(self, name, inputs=inputs,
                                          output_names = output_names)
        
  class WVSR(Backend.DSProc):
    """
    """
    def __init__(self, parent, name, inputs=None, output_names=None):
      """
      initialize a WVSR
      
      """
      self.name = name
      self.parent = parent
      mylogger = logging.getLogger(self.parent.name + ".WVSR")
      mylogger.debug("__init__: parent is %s", self.parent)
      Backend.DSProc.__init__(self, self.parent, self.name, inputs=inputs,
                              output_names=output_names)
      self.logger = mylogger
      self.logger.info("__init__: %s input channels: %s", self, self.inputs)
      self.logger.info("__init__: %s output channels: %s", self, self.outputs)
      self.cfg = self.parent.collector.wvsr_cfg[self.name]
      scankeys = self.parent.collector.fft_meta.keys()
      # this assumes that all FFT configurations are the same for all scans
      self.fft_cfg = self.parent.collector.fft_meta[scankeys[0]]
      self.data['log'] = self.cfg['from log']
      self.data['user'] = self.cfg['user']
      self.data['bandwidth'] = 640
      self.data['frequency'] = 320
      # define the IF channels
      self.IF = {}
      for inname in self.inputs.keys():
        self.logger.debug("__init__: processing input %s", inname)
        wvsr, IF = inname.split('.')
        IFnum = int(IF[-1])
        if self.IF.has_key(IFnum):
          self.logger.debug("__init__: %s already exists", self.IF[IFnum])
        else:
          # there is one input from each IF
          self.IF[IFnum] = WVSRbackend.WVSR.IFchannel(self, inname,
                                               inputs={inname: inputs[inname]},
                                               output_names=output_names)
    
    class IFchannel(Backend.DSProc.Channel):
      """
      Typically there are two, one for each polarization
      
      Input names are of the form 'wvsr2.IF1'.
      Output names are of the form 'X14.chan_id 1.I'.
      """
      def __init__(self, parent, name, inputs=None, output_names=None):
        """
        In an IF channel there is a one-to-one mapping of inputs to outputs
        """
        self.name = name
        self.parent = parent
        mylogger = logging.getLogger(self.parent.name + ".IFchannel")
        mylogger.debug("__init__: parent is %s", self.parent)
        Backend.DSProc.Channel.__init__(self, self.parent, self.name,
                                        inputs=inputs,
                                        output_names=output_names)
        self.logger = mylogger
        self.logger.info("__init__: %s input channels: %s", self, self.inputs)
        self.logger.info("__init__: %s output channels: %s",self, self.outputs)
        IFnum = int(self.name[-1])
        self.data['DSP'] = self.parent.cfg[IFnum]['DSP']
        self.data['expid'] = self.parent.cfg[IFnum]['expid']
        self.data['pol'] = self.parent.cfg[IFnum]['pol']
        self.data['LO'] = self.parent.cfg[IFnum]['rf_to_if_lo']
        # back propagate the LO and IF freq to the down-converter
        for key in self.inputs.keys(): # a formality; there is only one input
          self.inputs[key].source.parent['LO'] = self.data['LO']
        subchannels = self.parent.cfg[IFnum]['subchannels']
        self.subchannel = {}
        # FIX FIX FIX vvvvvvv
        #self.logger.debug("__init__: FFT cfg keys: %s",
        #                                     self.parent.fft_cfg[IFnum].keys())
        #nchans = self.parent.fft_cfg[2]['n_freqs'] # UGLY HACK FOR ERROR in 
        #                                                            helpers.py
        for inname in self.inputs.keys():
          # define the data channels; each IF shares the data channels
          for outname in self.outputs.keys():
            rx, subch, stokes = outname.split('.')
            if self.subchannel.has_key(subch):
              self.logger.debug("__init__: %s already exists", 
                                                        self.subchannel[subch])
            else:
              # create a subchannel for this IF
              self.subchannel[subch] =  WVSRbackend.WVSR.Subchannel(
                             self, subch, inputs={inname: self.inputs[inname]},
                             output_names=output_names)
            self.outputs[outname].source = \
                                        self.subchannel[subch].outputs[outname]
            #self.outputs[outname].source = self.inputs[inname]
            #self.outputs[outname].signal = Spectrum(self.inputs[inname].signal,
            #                                        num_chans=nchans)
          self.inputs[inname].destinations = self.outputs.values()
          # redefine the parent output in terms of the channel output
          # this attaches properties to the parents output.
            #self.parent.outputs[outname] = self.outputs[outname]
            #self.outputs[outname].signal['IF bandwidth'] = \
            #                                      self.parent.data['bandwidth']
            #self.parent.outputs[outname] = self.outputs[outname]
          

      #for subch in subchannels:
      #  
      #    
      #  else:
      #    
      #    # select the outputs for this subchannel
      #    for outname in output_names:
      #      rx, sub, stokes = outname.split('.')
      #      if sub == subch:
      #        outnames.append(outname)
      #    self.logger.debug("__init__: subchannel % outputs: %s", subch, outnames)
      #    
             
    class Subchannel(Backend.DSProc.Channel):
        """
        
        """
        def __init__(self, parent, name, inputs=None, output_names=None):
          """
          """
          self.name = name
          self.parent = parent
          mylogger = logging.getLogger(self.parent.name + ".Subchannel")
          Backend.DSProc.Channel.__init__(self, self.parent, self.name,
                                      inputs=inputs, output_names=output_names)
          self.logger = mylogger
          self.logger.info("__init__: %s input channels: %s", 
                           self, self.inputs)
          self.logger.info("__init__: %s output channels: %s",
                           self, self.outputs)
          cfg = self.parent.parent.cfg
          fft_cfg = self.parent.parent.fft_cfg
          IFnum = int(self.parent.name[-1])
          self.data['LO'] = cfg[IFnum][self.name]['sfro']
          self.data['bandwidth'] = cfg[IFnum][self.name]['bandwidth']
          self.data['nchans'] = fft_cfg['n_freqs'] 
          for outname in self.outputs.keys():
            stokes = outname[-1]
            for inname in self.inputs.keys():
              self.outputs[outname].source = self.inputs[inname]
              self.outputs[outname].signal = Spectrum(
                                           self.outputs[outname].source.signal,
                                           name=stokes,
                                           num_chans=self.data['nchans'])
              #self.outputs[outname].signal['IF frequency'] = \
              #                                    self.parent.data['frequency']
              #self.outputs[outname].signal['IF bandwidth'] = \
              #                                    self.parent.data['bandwidth']
              #self.parent.outputs[outname] = self.outputs[outname]
              pass

