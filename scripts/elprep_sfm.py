#!/usr/bin/env python

# Example script for using elPrep split/merge tools.  

import subprocess
import os
import time
import elprep_io_wrapper
import multiprocessing

# actual script

def elprep_sfm (argv):
  # set up directories for intermediate results
  file_in = argv[1]
  file_out = argv[2]
  stamp = str(time.time())
  split_dir = os.path.join(os.getcwd(), "temp-" + stamp + os.sep)
  result_dir = os.path.join(os.getcwd(), "temp-processed-" + stamp + os.sep)
  os.mkdir(split_dir)
  os.mkdir(result_dir)
  # split command
  nr_of_threads_opt = elprep_io_wrapper.cmd_option("--nr-of-threads", argv) or ["--nr-of-threads", str(multiprocessing.cpu_count())]

  if os.path.isdir(file_in):
    output_prefix, output_extension = os.path.splitext(os.path.basename(os.listdir(file_in)[0]))
  else: # file
    output_prefix, output_extension = os.path.splitext(os.path.basename(file_in))
  if output_extension == '': # e.g. /dev/stdin
    output_extension = '.sam'

  intermediate_files_opt = elprep_io_wrapper.cmd_option("--intermediate-files-output-type", argv) 
  if intermediate_files_opt:
    intermediate_files_output_type = intermediate_files_opt[1]
  else:
    intermediate_files_output_type = output_extension[1:]

  intermediate_files_op_opt = elprep_io_wrapper.cmd_option("--intermediate-files-output-prefix", argv) 
  if intermediate_files_op_opt:
    output_prefix = intermediate_files_op_opt[1]

  given_cmd_opts = elprep_io_wrapper.remove_cmd_option(argv[3:], "--intermediate-files-output-type")

  given_cmd_opts = elprep_io_wrapper.remove_cmd_option(given_cmd_opts, "--intermediate-files-output-prefix")

  cmd_opts = given_cmd_opts
  nr_of_threads_opt_given = elprep_io_wrapper.cmd_option("--nr-of-threads", cmd_opts)
  if not nr_of_threads_opt_given:
    cmd_opts = cmd_opts + nr_of_threads_opt # so we pass --nr-of-threads to the elprep filter command explicitly

  single_end_opt = elprep_io_wrapper.flg_option("--single-end", argv)

  elprep_io_wrapper.cmd_wrap_input(["elprep", "split"], file_in, split_dir, ["--output-prefix", output_prefix, "--output-type", intermediate_files_output_type] + nr_of_threads_opt + single_end_opt)

  if single_end_opt:
    splits_path = split_dir
    cmd_opts = elprep_io_wrapper.remove_flg_option(cmd_opts, "--single-end")
  else:
    splits_path = os.path.join(split_dir, "splits" + os.sep)

  # run filter command for split files
  for root, dirs, files in os.walk(splits_path):
    for file in files:
      ext = os.path.splitext(file)[1]
      if (ext == ".sam" or ext == ".bam" or ext == ".cram"):
        ffile = os.path.join(root, file)
        processed_file = os.path.join(result_dir, os.path.basename(file))
        elprep_io_wrapper.cmd_wrap_io(["elprep", "filter"], ffile, processed_file, cmd_opts)
        os.remove(ffile)

  # command for spread file
  # commands for split files and spread file are the same, but the files are stored in different folders 
  # we keep them in seperate folders for backwards compatibility with elPrep v2.61
  # in a future release this will be simplified, so all split files will be in the same folder
  if not single_end_opt:
    os.rmdir(splits_path)
    spread_file = os.path.join(split_dir, output_prefix + "-spread." + intermediate_files_output_type)
    spread_out_file = os.path.join(result_dir, output_prefix + "-spread." + intermediate_files_output_type)
    elprep_io_wrapper.cmd_wrap_io(["elprep", "filter"], spread_file, spread_out_file , cmd_opts)
    os.remove(spread_file)
  os.rmdir(split_dir)

  # merge command
  elprep_io_wrapper.cmd_wrap_output(["elprep", "merge"], result_dir, file_out, nr_of_threads_opt + single_end_opt)
  # remove directories for intermediate results
  for root, dirs, files in os.walk(result_dir):
    for file in files:
      ffile = os.path.join(root, file)
      os.remove(ffile)
  os.rmdir(result_dir)
