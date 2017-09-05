"""
"""
import logging

from MonitorControl.BackEnds import Backend, get_freq_array
from MonitorControl.BackEnds.DSN.helpers import WVSRmetadataCollector

logger = logging.getLogger(__name__)

class WVSRbackend(Backend):
  """
  """
  def __init__(self, name, collector, inputs=None, output_names=None):
    """
    """
    mylogger = logging.getLogger(logger.name+".WVSRbackend")
    Backend.__init__(self, name, inputs=inputs, output_names=output_names)
    self.collector = collector
    self.logger = mylogger
    self.logger.info("__init__: %s input channels: %s", self, self.inputs)
    self.logger.info("__init__: %s output channels: %s", self, self.outputs)
    self.wvsrs = {}
    for name in collector.wvsrnames: # there is only one WVSR
      self.wvsrs[name] = WVSRbackend.WVSR(self, name, inputs=inputs,
                                          output_names = output_names)
        
  class WVSR(Backend.DSProc):
    """
    """
    def __init__(self, parent, name, inputs=None, output_names=None):
      """
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
      collector = self.parent.collector
      self.IF = {}
      for IF in collector.wvsr_cfg[name]['channels']:
        self.IF[IF] = WVSRbackend.WVSR.IFchannel(self, self.name,
                                                 inputs=inputs,
                                                 output_names=output_names)
    
    class IFchannel(Backend.DSProc.Channel):
      """
      """
      def __init__(self, parent, name, inputs=None, output_names=None):
        """
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
        self.logger.info("__init__: %s output channels: %s",
                         self, self.outputs)
        collector = self.parent.parent.collector
        IFnumber = int(self.name[-1])
        self.subchannel = {}
        self.logger.debug("__init__: %s has keys %s",
                          collector.wvsr_cfg[self.parent.name][IFnumber],
                         collector.wvsr_cfg[self.parent.name][IFnumber].keys())
        for subch in \
                 collector.wvsr_cfg[self.parent.name][IFnumber]['subchannels']:
          self.subchannel[subch] =  WVSRbackend.WVSR.IFchannel.Subchannel(
                         self, subch, inputs=inputs, output_names=output_names)
             
      class Subchannel(Backend.DSProc.Channel.DataChl):
        """
        """
        def __init__(self, parent, name, inputs=None, output_names=None):
          """
          """
          self.name = name
          self.parent = parent
          mylogger = logging.getLogger(self.parent.name + ".Subchannel")
          Backend.DSProc.Channel.DataChl.__init__(self, self.parent, self.name,
                                      inputs=inputs, output_names=output_names)
          self.logger = mylogger
          self.logger.info("__init__: %s input channels: %s", 
                           self, self.inputs)
          self.logger.info("__init__: %s output channels: %s",
                           self, self.outputs)
          collector = self.parent.parent.parent.collector
          
      self.DC[inname].outputs[outname].source = self.DC[inname].inputs[inname]
      self.DC[inname].outputs[outname].source.destinations.append(
                                             self.DC[inname].outputs[outname])
      self.DC[inname].outputs[outname].signal = IF(
                           self.DC[inname].outputs[outname].source.signal,'U')
      self.DC[inname].outputs[outname].signal['IF frequency'] = self.data['frequency']
      self.DC[inname].outputs[outname].signal['IF bandwidth'] = self.data['bandwidth']
      self.outputs[outname] = self.DC[inname].outputs[outname]

