#!/usr/bin/env python

# Example script for using elPrep split/merge tools.  

import subprocess
import os
import time
import elprep_io_wrapper
import operator
import multiprocessing

# actual script

def append_cmd(cmd1, cmd2):
  return cmd1 + " " + cmd2

def elprep_sfm_gnupar (argv):
  # set up directories for intermediate results
  file_in = argv[1]
  file_out = argv[2]
  stamp = str(time.time())
  split_dir = os.path.join(os.getcwd(), "temp-" + stamp + os.sep)
  result_dir = os.path.join(os.getcwd(), "temp-processed-" + stamp + os.sep)
  os.mkdir(split_dir)
  os.mkdir(result_dir)
  nr_of_jobs_opt = elprep_io_wrapper.cmd_option('--nr-of-jobs', argv) or ['--nr-of-jobs', 1]
  nr_of_threads_opt = elprep_io_wrapper.cmd_option('--nr-of-threads', argv) or ['--nr-of-threads', str((multiprocessing.cpu_count() / int(nr_of_jobs_opt[1])) or 1)]
  # split command
  nr_of_split_merge_threads_opt = ['--nr-of-threads', str(int(nr_of_jobs_opt[1]) * int(nr_of_threads_opt[1]))]
 
  if os.path.isdir(file_in):
    output_prefix, output_extension = os.path.splitext(os.path.basename(os.listdir(file_in)[0]))
  else: # file
    output_prefix, output_extension = os.path.splitext(os.path.basename(file_in))
  if output_extension == '': # e.g. /dev/stdin
    output_extension = '.sam'

  intermediate_files_opt = elprep_io_wrapper.cmd_option('--intermediate-files-output-type', argv) 
  if intermediate_files_opt:
    intermediate_files_output_type = intermediate_files_opt[1]
  else:
    intermediate_files_output_type = output_extension[1:] 

  intermediate_files_op_opt = elprep_io_wrapper.cmd_option('--intermediate-files-output-prefix', argv) 
  if intermediate_files_op_opt:
    output_prefix = intermediate_files_op_opt[1]

  single_end_opt = elprep_io_wrapper.flg_option("--single-end", argv)

  elprep_io_wrapper.cmd_wrap_input(["elprep", "split"], file_in, split_dir, ["--output-prefix", output_prefix, "--output-type", intermediate_files_output_type] + nr_of_split_merge_threads_opt + single_end_opt)

  # gnu parallel command
  read_group_string = elprep_io_wrapper.cmd_option('--replace-read-group', argv)
  given_cmd_opts = elprep_io_wrapper.remove_cmd_option(argv[3:], '--nr-of-jobs')
  given_cmd_opts = elprep_io_wrapper.remove_cmd_option(given_cmd_opts, '--intermediate-files-output-type')
  given_cmd_opts = elprep_io_wrapper.remove_cmd_option(given_cmd_opts, '--intermediate-files-output-prefix')
  if single_end_opt:
    splits_path = split_dir
    given_cmd_opts = elprep_io_wrapper.remove_flg_option(given_cmd_opts, "--single-end")
  else:
    splits_path = os.path.join(split_dir, "splits" + os.sep)

  given_cmd_opts = elprep_io_wrapper.remove_cmd_option(given_cmd_opts, '--nr-of-threads')
  cmd_opts = given_cmd_opts + nr_of_threads_opt
  if read_group_string:
    cmd_opts = elprep_io_wrapper.remove_cmd_option(cmd_opts, '--replace-read-group') 
    cmd_opts = cmd_opts + ['--replace-read-group', '\"' + read_group_string[1] + '\"']
  cmd_list = ["elprep filter"]
  elprep_cmd = '\'' + reduce(append_cmd, cmd_list + ['{}', result_dir + '{/.}.' + intermediate_files_output_type ] + cmd_opts) + '\''
  gnu_cmd = 'parallel --gnu -j ' + str(nr_of_jobs_opt[1]) + ' ' + elprep_cmd + ' ::: ' + splits_path + '*.' + intermediate_files_output_type
  ret = subprocess.check_call(gnu_cmd, shell=True)
  if ret != 0: raise SystemExit, ret

  # command for spread file

  if not single_end_opt:
    spread_file = os.path.join(split_dir, output_prefix + "-spread." + intermediate_files_output_type)
    spread_out_file = os.path.join(result_dir, output_prefix + "-spread." + intermediate_files_output_type)
    elprep_io_wrapper.cmd_wrap_io(["elprep", "filter"], spread_file, spread_out_file , given_cmd_opts + nr_of_split_merge_threads_opt)
    os.remove(spread_file)

  # merge command
  elprep_io_wrapper.cmd_wrap_output(["elprep", "merge"], result_dir, file_out, nr_of_split_merge_threads_opt + single_end_opt)
  # remove directories for intermediate results
  for root, dirs, files in os.walk(splits_path):
    for file in files:
      ffile = os.path.join(root, file)
      os.remove(ffile)
  if not single_end_opt:
    os.rmdir(splits_path)
  os.rmdir(split_dir)
  for root, dirs, files in os.walk(result_dir):
    for file in files:
      ffile = os.path.join(root, file)
      os.remove(ffile)
  os.rmdir(result_dir)
