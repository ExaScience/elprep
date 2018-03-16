#!/usr/bin/env python

import sys
import subprocess
import os
import multiprocessing

# utilities

def member(el,l):
  i = iter(l)
  for val in i:
    if val == el:
      rest = list(i)
      return [el] + rest

def flg_option(opt, cmdl):
  val = member(opt, cmdl)
  return [opt] if val else []

def cmd_option(opt, cmdl):
  val = member(opt, cmdl)
  return [opt,val[1]] if val else []

def remove_flg_option(cmd_l, name):
  nl = []
  i = iter(cmd_l)
  val = next(i, None)
  while(val):
    if not(val == name):
       nl = nl + [val]
    val = next(i, None)
  return nl

def remove_cmd_option(cmd_l, name):
  "removes a command with its option value from a list of commands"
  nl = []
  i = iter(cmd_l)
  val = next(i, None)
  while(val):
    if not(val == name):
       nl = nl + [val]
       val = next(i, None)
    else:
       next(i, None)
       val = next(i, None)   
  return nl

def cmd_wrap_input (cmd_list, file_in, file_out, cmd_opts):
  input = os.path.basename(file_in)
  input_prefix, input_extension = os.path.splitext(input)
  nr_of_threads_opt = cmd_option("--nr-of-threads", cmd_opts)
  nr_of_threads = str(nr_of_threads_opt[1]) if nr_of_threads_opt else str(multiprocessing.cpu_count()) 
  if (input_extension == ".bam" or input_extension == ".cram"):
    p1 = subprocess.Popen(["samtools", "view", "-h", "-@", nr_of_threads, file_in], bufsize=-1, stdout=subprocess.PIPE)          
    p2 = subprocess.Popen(cmd_list + ["/dev/stdin", file_out] + cmd_opts, bufsize=-1, stdin=p1.stdout)
    p2.communicate()
    if p2.returncode != 0: raise SystemExit, p2.returncode
  else:
    ret = subprocess.call(cmd_list + [file_in, file_out] + cmd_opts)
    if ret != 0: raise SystemExit, ret

def cmd_wrap_output (cmd_list, file_in, file_out, cmd_opts):
  output = os.path.basename(file_out)
  output_prefix, output_extension = os.path.splitext(output)
  nr_of_threads_opt = cmd_option("--nr-of-threads", cmd_opts)
  nr_of_threads = str(nr_of_threads_opt[1]) if nr_of_threads_opt else str(multiprocessing.cpu_count()) 
  if output_extension == ".bam":
    p1 = subprocess.Popen(cmd_list + [file_in, "/dev/stdout"] + cmd_opts, bufsize=-1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["samtools", "view", "-bS", "-@", nr_of_threads, "-o", file_out, "-"], bufsize=-1, stdin=p1.stdout) 
    p2.communicate()
    if p2.returncode != 0: raise SystemExit, p2.returncode
  if output_extension == ".cram":
    reference_t_opt = cmd_option("--reference-t", cmd_opts)
    reference_bigT_opt = cmd_option("--reference-T", cmd_opts)
    if not(reference_t_opt) and not(reference_bigT_opt): 
        print("Converting to .cram. Need to pass reference-t or reference-T")
        return
    opt_to_delete = "--reference-t" if reference_t_opt else "--reference-T"
    p1 = subprocess.Popen(cmd_list + [file_in, "/dev/stdout"] + remove_cmd_option(cmd_opts, opt_to_delete), bufsize=-1, stdout=subprocess.PIPE)
    t_opt = ["-t", reference_t_opt[1]] if reference_t_opt else ["-T", reference_bigT_opt[1]]
    p2 = subprocess.Popen(["samtools", "view", "-C", "-@", nr_of_threads] + t_opt + ["-o", file_out, "-"], bufsize=-1, stdin=p1.stdout)
    p2.communicate()
    if p2.returncode != 0: raise SystemExit, p2.returncode
  else:
    ret = subprocess.call(cmd_list + [file_in, file_out] + cmd_opts)
    if ret != 0: raise SystemExit, ret

def cmd_wrap_io(cmd_list, file_in, file_out, cmd_opts):
  input = os.path.basename(file_in)
  output = os.path.basename(file_out)
  output_prefix, output_extension = os.path.splitext(output)
  input_prefix, input_extension = os.path.splitext(input)
  nr_of_threads_opt = cmd_option("--nr-of-threads", cmd_opts)
  nr_of_threads = str(nr_of_threads_opt[1]) if nr_of_threads_opt else str(multiprocessing.cpu_count())  
  if (input_extension == ".bam" or input_extension == ".cram"):
    p1 = subprocess.Popen(["samtools", "view", "-h", "-@", nr_of_threads, file_in], bufsize=-1, stdout=subprocess.PIPE)          
    if (output_extension == ".bam"):  
      p2 = subprocess.Popen(cmd_list + ["/dev/stdin", "/dev/stdout"] + cmd_opts, bufsize=-1, stdin=p1.stdout, stdout=subprocess.PIPE)
      p3 = subprocess.Popen(["samtools", "view", "-bS", "-@", nr_of_threads, "-o", file_out, "-"], bufsize=-1, stdin=p2.stdout) 
      p3.communicate()
      if p3.returncode != 0: raise SystemExit, p3.returncode
    elif (output_extension == ".cram"):
      reference_t_opt = cmd_option("--reference-t", cmd_opts)
      reference_bigT_opt = cmd_option("--reference-T", cmd_opts)
      if not(reference_t_opt) and not(reference_bigT_opt): 
          print("Converting to .cram. Need to pass reference-t or reference-T")
          return
      opt_to_delete = "--reference-t" if reference_t_opt else "--reference-T"
      p2 = subprocess.Popen(cmd_list + ["/dev/stdin", "/dev/stdout"] + remove_cmd_option(cmd_opts, opt_to_delete), bufsize=-1, stdin=p1.stdout, stdout=subprocess.PIPE)
      t_opt = ["-t", reference_t_opt[1]] if reference_t_opt else ["-T", reference_bigT_opt[1]]
      p3 = subprocess.Popen(["samtools", "view", "-C", "-@", nr_of_threads] + t_opt + ["-o", file_out, "-"], bufsize=-1, stdin=p2.stdout)
      p3.communicate()
      if p3.returncode != 0: raise SystemExit, p3.returncode
    else: # output_extension == ".sam"
      p2 = subprocess.Popen(cmd_list + ["/dev/stdin", file_out] + cmd_opts, bufsize=-1, stdin=p1.stdout)
      p2.communicate()
      if p2.returncode != 0: raise SystemExit, p2.returncode
  else: # input_extension == ".sam"
    if (output_extension == ".bam"):
      p1 = subprocess.Popen(cmd_list + [file_in, "/dev/stdout"] + cmd_opts, bufsize=-1, stdout=subprocess.PIPE)
      p2 = subprocess.Popen(["samtools", "view", "-bS", "-@", nr_of_threads, "-o", file_out, "-"], bufsize=-1, stdin=p1.stdout)
      p2.communicate()
      if p2.returncode != 0: raise SystemExit, p2.returncode
    elif (output_extension == ".cram"):
      reference_t_opt = cmd_option("--reference-t", cmd_opts)
      reference_bigT_opt = cmd_option("--reference-T", cmd_opts)
      if not(reference_t_opt) and not(reference_bigT_opt): 
          print("Converting to .cram. Need to pass reference-t or reference-T")
          return
      opt_to_delete = "--reference-t" if reference_t_opt else "--reference-T"
      p1 = subprocess.Popen(cmd_list + [file_in, "/dev/stdout"] + remove_cmd_option(cmd_opts, opt_to_delete), bufsize=-1, stdout=subprocess.PIPE)
      t_opt = ["-t", reference_t_opt[1]] if reference_t_opt else ["-T", reference_bigT_opt[1]]
      p2 = subprocess.Popen(["samtools", "view", "-C", "-@", nr_of_threads] + t_opt + ["-o", file_out, "-"], bufsize=-1, stdin=p1.stdout)
      p2.communicate()
      if p2.returncode != 0: raise SystemExit, p2.returncode
    else: # output_extension == ".sam"    
      ret = subprocess.call(cmd_list + [file_in, file_out] + cmd_opts)
      if ret != 0: raise SystemExit, ret
