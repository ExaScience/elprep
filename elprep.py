#!/usr/bin/env python

import sys
import elprep_io_wrapper

def elprep():
  elprep_io_wrapper.cmd_wrap_io(["elprep"], sys.argv[1], sys.argv[2], sys.argv[3:])

elprep()
