#!/usr/bin/env python

##########################################################################################
##
## CNV_CopywriterGetModeCorrectionFactor.py
##
## Calculates the location parameter using Mode, then infer correction factor.
##
##########################################################################################

import numpy as np
import os
import sys
from scipy import stats

def myMode(a):
        sarr = np.sort(a)
        kde = stats.gaussian_kde(sarr)
        y = kde.evaluate(sarr)
        peak = sarr[y.argmax()]
        return peak

def CalculateMode(name):
    CopyWriterPath = os.getcwd()+"/"+name+"/results/Copywriter/"
    sub ="SegmentsChromosome"
    files = os.listdir(CopyWriterPath)
    fileArray = [s for s in files if sub in s]
    ModeLocations = []
    print("Collect Segment data")
    for file in fileArray:
          logRR = np.loadtxt(CopyWriterPath+file)
          if logRR.size>1:
              myPeak = myMode(logRR)
              ModeLocations.append(myPeak)
          if logRR.size==1:
              logRR = logRR.tolist()
              ModeLocations.append(logRR)
    ModeLocations = np.array(ModeLocations)
    print("Calculate Mode correction for Sample "+name)
    myPeak = myMode(ModeLocations)
    print("Writing mode correction estimate to "+CopyWriterPath+"/"+name+".KDEestimate.txt")
    OutFile = CopyWriterPath+"/"+name+".KDEestimate.txt"
    f = open(OutFile, 'w')
    f.write("name\tShift\n")
    f.write(name+"\t"+str(myPeak)+"\n")
    f.close()

name = sys.argv[1]
CalculateMode(name)