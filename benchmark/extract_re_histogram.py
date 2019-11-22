#! /usr/bin/env python

import matplotlib.pyplot as pyplot
import numpy
import sys
import vtk

if len(sys.argv) <= 1:
  print('No input file given.')
  sys.exit()

infiles = sys.argv[1:]
#print(f'infiles\n{infiles}')
colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black']
coloridx = 0
for infile in infiles:
  vtkreader = vtk.vtkXMLPolyDataReader()
  vtkreader.SetFileName(infile)
  vtkreader.Update()
  vtknumcells = vtkreader.GetNumberOfCells()
  vtkrerr = vtkreader.GetOutput().GetCellData().GetArray("relative-error")
  vtkhitcnt = vtkreader.GetOutput().GetCellData().GetArray("hitcnt")
  rerr = []
  hitcnt = []
  for idx in range(vtknumcells):
    vale = vtkrerr.GetValue(idx)
    valh = vtkhitcnt.GetValue(idx)
    rerr.append(vale)
    hitcnt.append(valh)
  # filter elements which have no hits at all
  rerr = [e1 for e1, e2 in filter(lambda e: e[1]!=0, zip(rerr,hitcnt))]
  numbins = 64
  n, bins, patches = pyplot.hist(rerr, numbins, facecolor=colors[coloridx % len(colors)], alpha=0.5)
  coloridx += 1

pyplot.show()
