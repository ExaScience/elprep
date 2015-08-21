#!/usr/bin/env python

import elprep_io_wrapper

def elprep_im(argv):
    elprep_io_wrapper.cmd_wrap_io(["elprep"], argv[1], argv[2], argv[3:])
