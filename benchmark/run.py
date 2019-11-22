#! /usr/bin/env python

import functools
from glob import glob
import math
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pickle
import subprocess
import vtk

# Settings
numruns = 20
inpath = "./input/surf_1"
stickingcoefficient = 0.2
numofrays = 8e7
numofthreads = 1
cmd = "../build/rti_m --infile {} --outfile {} --sticking-coefficient {} --triangles --number-of-rays {} --max-threads {}" # format string

# Prepare
infiles = glob(inpath + '**/*', recursive=True)
infiles = filter(lambda ss: os.path.isfile(ss), infiles)
infiles = map(lambda ss: [ss, os.path.dirname(ss), os.path.basename(ss)], infiles)

# Define functions
def run_processes():
  # Compile
  subprocess.run(["(cd ../build; cmake --build . --target rti)"], shell=True)
  #
  def handle(triple):
    global cmd
    inname, dirname, basename = triple
    #print(f'{inname} {dirname} {basename}\n')
    # reduced path
    intermediatePath = os.path.relpath(dirname, inpath)
    extbasename = "result--" + basename
    #print(f'intermediatePath == {intermediatePath}')
    listoutpaths = []
    for runnumber in range(numruns):
      outpath = os.path.join("./output", intermediatePath, basename + '.d', str(runnumber))
      os.makedirs(outpath, exist_ok=True)
      outname = os.path.join(outpath, extbasename);
      concretecmd = cmd.format(inname, outname, stickingcoefficient, numofrays, numofthreads)
      #print(f'cmd == {cmd}')
      # run command. We do not save standard output.
      subprocess.run(concretecmd, shell=True)
      # Collect outpaths
      listoutpaths.append(outname)
    return listoutpaths
  result = []
  for triple in infiles:
    tt = handle(triple)
    result.extend(tt)
  ## Note: the results is a list of out-file-paths without the '.vtp' suffix
  return result

def extract_data(pp):
  #print(pp)

  #def handle(outfile):
  #  vtpfilepath = outfile + '.vtp'
  #  vtkreader = vtk.vtkXMLPolyDataReader()
  #  #print(vtpfilepath)
  #  vtkreader.SetFileName(vtpfilepath)
  #  vtkreader.Update()
  #  vtkrunningtime = vtkreader.GetOutput().GetFieldData().GetAbstractArray("running-time[ns]")
  #  assert vtkrunningtime.GetNumberOfTuples() == 1, "Error: expected one value for the running time"
  #  runningtime = vtkrunningtime.GetValue(0) # get first value in vtk array
  #  #print(f'running time: {runningtime}')

  def handle(triple, printtoterminal=False):
    dirpath, maxidxstr, filename = triple
    maxidx = int(maxidxstr)
    numelems = maxidx + 1
    #print(f'{dirpath}, {maxidx}, {filename}')
    sum = 0
    sumsqrs = 0
    for idx in range(numelems):
      vtpfilepath = os.path.join(dirpath, str(idx), filename) + '.vtp'
      vtkreader = vtk.vtkXMLPolyDataReader()
      #print(vtpfilepath)
      vtkreader.SetFileName(vtpfilepath)
      vtkreader.Update()
      vtkrunningtime = vtkreader.GetOutput().GetFieldData().GetAbstractArray("running-time[ns]")
      assert vtkrunningtime.GetNumberOfTuples() == 1, "Error: expected one value for the running time"
      runningtime = int(vtkrunningtime.GetValue(0)) # get first value in vtk array
      sum += runningtime
      sumsqrs += runningtime ** 2
      #print(f'running time: {runningtime}')
    avg = sum / numelems
    variance = (sumsqrs / numelems - avg ** 2) * numelems / (numelems - 1)
    ci95delta = 1.96 * math.sqrt(variance / numelems)
    ci95 = [avg + ci95delta, avg - ci95delta]
    #print(f'dirpath: {dirpath} avg: {avg}, 0.95 confidence interval: {ci95}')
    if printtoterminal: print(f'{dirpath}, {avg}, {ci95[0]}, {ci95[1]}')
    return [dirpath, avg, ci95[0], ci95[1]]

  def split(aa):
    head1, tail1 = os.path.split(aa)
    assert head1 != aa
    head2, tail2 = os.path.split(head1)
    #print(f'head2: {head2}, tail2: {tail2}, tail1: {tail1}')
    return [head2, tail2, tail1]

  bb = list(map(split, pp))
  #print(f'bb: {bb}')
  # filter the list bb such that only the highest values in the second column remain
  col = [ee for ee in bb if int(ee[1]) == max([int(ii[1]) for ii in bb])]
  #print(f'cc: {cc}')
  printtoterminal = True
  if printtoterminal: print("# filename, average running time [ns], time confidence interval low [ns], time confidence interval high [ns]")
  acc = []
  for triple in col:
    acc.append( handle(triple, printtoterminal) )
  return acc

def plot(list_lt, list_m, list_p, outfile, xtiks):
  print('in plot()')

  plt.close('all')
  plt.figure(0)
  plt.figure(figsize=(9,5))
  plt.rcParams.update({'font.size': 22})

  #print(f'list_lt\n{list_lt}')
  # Note: in numpy an array holds data of a particular type. That is, the data
  # is converted to string if there are any strings in the input.
  list_lt = np.array(list_lt)
  #print(f'list_lt\n{list_lt}')
  list_lt = list_lt[list_lt[:,0].argsort()] # sort by first column
  #print(f'list_lt\n{list_lt}')
  list_m = np.array(list_m)
  list_m = list_m[list_m[:,0].argsort()] # sort by first column
  assert len(list_p) == 1
  list_p = list_p[0]
  print(f'list_lt\n{list_lt}')
  print(f'list_m\n{list_m}')
  print(f'list_p\n{list_p}')
  #errorbars_lt = [[e[1]-e[2], e[3]-e[1]] for e in list_lt]
  #errorbars_lt = list(zip([e[0] for e in errorbars_lt], [e[1] for e in errorbars_lt]))
  errorbars_lt, errorbars_m = map(lambda ll: [(e[1]-e[2]) for e in ll], [list_lt, list_m])
  #print(f'errorbars_lt\n{errorbars_lt}')
  plt.errorbar(range(len(list_lt)), list_lt[:,1], yerr=errorbars_lt, fmt='-r', clip_on = False, label="LT")
  plt.errorbar(range(len(list_m)), list_m[:,1], yerr=errorbars_m, fmt='-b', clip_on = False, label="M")
  plt.hlines(list_p[1], xmin=0, xmax=6, colors='g', label="Orgnl")
  #plt.plot(list_m[:,1], '-b^', clip_on = False, label="LT")

  plt.legend(loc=0, prop={'size': 20})
  plt.xlabel('percentage of vertices of the original mesh remaining')
  plt.ylabel('running time of ray tracing [sec]')
  plt.xticks(np.arange(7), xtiks)
  plt.savefig(outfile, bbox_inches='tight', dpi=1200, format='eps')


############
### Main ###
############

# pickle files for saving data
picklefile1='picklefile1'
picklefile2='picklefile2'

####
#tt = run_processes()
### write list tt to file
#with open(picklefile1, 'wb') as fp:
#  pickle.dump(tt, fp)
#
### read list tt from file
#with open(picklefile1, 'rb') as fp:
#  tt = pickle.load(fp)
#  # process
#  data = extract_data(tt)
####
#with open(picklefile2, 'wb') as fp:
#  pickle.dump(data, fp)

with open(picklefile2, 'rb') as fp:
  data = pickle.load(fp)

  toSec = 10**-9
  # convert nano-seconds to seconds
  data = [[e[0], toSec*e[1], toSec*e[2], toSec*e[3]] for e in data]

  #print(f'data\n{data}')
  surf_1_data = [ee for ee in data if "Surface_1" in ee[0]]
  surf_2_data = [ee for ee in data if "Surface_2" in ee[0]]
  #print(f'data1:\n{surf_1_data}')
  #print(f'data2:\n{surf_2_data}')

  surf_2_lt_data = [ ee for ee in surf_2_data if "_lt_" in ee[0]]
  surf_2_m_data = [ ee for ee in surf_2_data if "_m_" in ee[0]]
  surf_2_p_data = [ ee for ee in surf_2_data if "_p_" in ee[0]]

  surf_1_lt_data = [ee for ee in surf_1_data if "_lt_" in ee[0]]
  surf_1_m_data = [ee for ee in surf_1_data if "_m_" in ee[0]]
  surf_1_p_data = [ee for ee in surf_1_data if "_p_" in ee[0]]

  #print(f'surf_1_lt_data\n{surf_1_lt_data}')
  #print(f'surf_1_m_data\n{surf_1_m_data}')
  #print(f'surf_1_p_data\n{surf_1_p_data}')

  # A function that takes a list of lists and filters the first column for digits
  func1 = lambda lol: [ [int(''.join(list(filter(str.isdigit, ee[0]))))] + ee[1:] for ee in lol ]
  # Apply the function to each data array. Effectively removing the strings
  surf_2_lt_data, surf_2_m_data, surf_2_p_data = map(lambda e1: func1(e1), [surf_2_lt_data, surf_2_m_data, surf_2_p_data])
  surf_1_lt_data, surf_1_m_data, surf_1_p_data = map(lambda e1: func1(e1), [surf_1_lt_data, surf_1_m_data, surf_1_p_data])

  #print('after filtering')
  #print(f'surf 2 lt data:\n{surf_2_lt_data}')
  #print(f'surf 2 m data:\n{surf_2_m_data}')
  #print(f'surf 2 p data:\n{surf_2_p_data}')
  #print(f'surf_1_lt_data\n{surf_1_lt_data}')
  #print(f'surf_1_m_data\n{surf_1_m_data}')
  #print(f'surf_1_p_data\n{surf_1_p_data}')

  if surf_2_lt_data and surf_2_m_data and surf_2_p_data:
    plot(surf_2_lt_data, surf_2_m_data, surf_2_p_data, "figure_surf_2.eps", ('72%','47%','33%', '21%','12%','5%','3%'))
  if surf_1_lt_data and surf_1_m_data and surf_1_p_data:
    plot(surf_1_lt_data, surf_1_m_data, surf_1_p_data, "figure_surf_1.eps", ('60%','42%','36%', '21%','11%','4%','2%'))
