#!/usr/bin/env python

# Example script for using elPrep split/merge tools.  

import sys
import subprocess
import os
import time
import elprep_io_wrapper

# actual script

def elprep_sfm ():
  # set up directories for intermediate results
  file_in = sys.argv[1]
  file_out = sys.argv[2]
  input = os.path.basename(file_in)
  output_prefix, output_extension = os.path.splitext(input)
  stamp = str(time.time())
  split_dir = os.path.join(os.getcwd(), "temp-" + stamp + os.sep)
  result_dir = os.path.join(os.getcwd(), "temp-processed-" + stamp + os.sep)
  os.mkdir(split_dir)
  os.mkdir(result_dir)
  # split command
  nr_of_threads_opt = elprep_io_wrapper.cmd_option("--nr-of-threads", sys.argv)
  cmd_opts = sys.argv[3:]
  elprep_io_wrapper.cmd_wrap_input(["elprep", "split"], file_in, split_dir, ["--output-prefix", output_prefix, "--output-type", "sam"] + nr_of_threads_opt)
  # run filter command for split files
  for root, dirs, files in os.walk(split_dir):
    for file in files:
      ext = os.path.splitext(file)[1]
      if (ext == ".sam" or ext == ".bam" or ext == ".cram"):
        ffile = os.path.join(root, file)
        processed_file = os.path.join(result_dir, os.path.basename(file))
        elprep_io_wrapper.cmd_wrap_io(["elprep"], ffile, processed_file, cmd_opts + ["--split-file"])
        os.remove(ffile)
  os.rmdir(split_dir)
  # merge command
  elprep_io_wrapper.cmd_wrap_output(["elprep", "merge"], result_dir, file_out, nr_of_threads_opt)
  # remove directories for intermediate results
  for root, dirs, files in os.walk(result_dir):
    for file in files:
      ffile = os.path.join(root, file)
      os.remove(ffile)
  os.rmdir(result_dir)

elprep_sfm()
