#!/usr/bin/env python

# Wrapper script for using elPrep, elPrep sfm, and elPrep gnupar   
# Can be used as an entrypoint for Docker

import sys
import os
import subprocess
import elprep_im
import elprep_sfm
import elprep_sfm_gnupar

# actual script

def elprep_entrypoint():
  if len(sys.argv) >= 2:
    script_selector = sys.argv[1]
    if (script_selector == "sfm"):
      elprep_sfm.elprep_sfm(sys.argv[1:])
    elif (script_selector == "sfm-gnupar"):
      elprep_sfm_gnupar.elprep_sfm_gnupar(sys.argv[1:])
    else:
      elprep_im.elprep_im(sys.argv)
  else:
    subprocess.call(["elprep"])

elprep_entrypoint()    
