import sys, getopt
import os
import subprocess
from subprocess import Popen, PIPE
import ROOT
from ROOT import *
import re
import cmath

mypdg = 0

if(len(sys.argv)) < 2:
  print "Must have at least one argument (how many files to loop over)"
  sys.exit()
if(len(sys.argv)) == 3:
  mypdg = int(sys.argv[2])
  if mypdg == 13:
    print "Sweet! Lookin' for muons"
  elif mypdg == 111:
    print "Okay, pions are more your thing"
  elif mypdg == 2212:
    print "Protons! Those are cool, i guess"
  else:
    print mypdg
    print "We don't have any of that."
    sys.exit()

outf = open("pathlist.txt",'w')

filelim = int(sys.argv[1])
filenames = ''
cmd2 = 'samweb list-definition-files prodgenie_bnb_nu_uboone_mcc7_ana'
filenames = os.popen(cmd2).read()
files = filenames.splitlines()

firstfile = list(files)[0]
cmd3 = 'samweb locate-file %s' %(firstfile)
firstfile_loc = os.popen(cmd3).read()
# Now, truncate the dumb samweb format extra stuff off of the path
filepath = re.split(':|\(',firstfile_loc)[1]+'/'

ct = 0
for File in files:
  ct += 1
  outf.write(str(filepath+File+"\n"))
  if ct == filelim:
    break
outf.close()

n_evt = 0
n_pass = 0

# Load up our plotter
gROOT.ProcessLine('.L anatree_looper.C+')

gROOT.ProcessLine('loop(%i)'%(mypdg))

gROOT.ProcessLine('draw()')
