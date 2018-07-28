"""
what happens here


metadata gathered

indices initialized
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
          increment tone_row_index
    average adjacent channels to get desired resolution
    put STOKES data into DATA cube
    increment data row index (data_row_index += 1)
    
Since each subchannel has its own tones, each scan will num_tones x num_subch rows
    
    
"""
