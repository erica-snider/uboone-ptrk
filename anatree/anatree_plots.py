import sys, getopt
import os
import subprocess
from subprocess import Popen, PIPE
import ROOT
from ROOT import *
import re
import cmath

hdisttovert_kalman = ROOT.TH1D("Distance to Vertex","Distance to Vertex (Kalman); cm, presumably; entries",200,0,100)
hdisttovert_pandora = ROOT.TH1D("Distance to Vertex","Distance to Vertex (Pandora); cm, presumably; entries",200,0,100)
hclosestapproach_kalman = ROOT.TH1D("Closest approach","Closest Approach (Kalman); cm",200,0,100)
hclosestapproach_pandora = ROOT.TH1D("Closest approach","Closest Approach (Pandora); cm",200,0,100)
htrackangle_pandora = ROOT.TH1D("Track Angle","Track Angle (Pandora); radians", 200, 0, 6.5)
htrackangle_kalman = ROOT.TH1D("Track Angle", "Track Angle (Kalman); radians", 200, 0, 6.5)

suffix = ""

if(len(sys.argv)) < 2:
  print "Must have at least one argument (how many files to loop over)"
  sys.exit()
if(len(sys.argv)) == 3:
  mypdg = int(sys.argv[2])
  if mypdg == 13:
    print "Sweet! Lookin' for muons"
    suffix = "_mu"
  elif mypdg == 111:
    print "Okay, pions are more your thing"
    suffix = "_pi0"
  elif mypdg == 2212:
    print "Protons! Those are cool, i guess"
    suffix = "_p"
  else:
    print mypdg
    print "We don't have any of that."
    sys.exit()

filelim = int(sys.argv[1])
filepaths= []
filenames = ''
cmd2 = 'samweb list-definition-files prodgenie_bnb_nu_uboone_mcc7_ana'
filenames = os.popen(cmd2).read()
files = filenames.splitlines()

if len(files) > 0:
  firstfile = list(files)[0]
  cmd3 = 'samweb locate-file %s' %(firstfile)
  firstfile_loc = os.popen(cmd3).read()
  # Now, truncate the dumb samweb format extra stuff off of the path
  filepath = re.split(':|\(',firstfile_loc)[1]+'/'

ct = 0
for File in files:
  ct += 1
  filepaths.append(str(filepath+File))
  if ct == filelim:
    break

n_evt = 0
n_pass = 0

for myfile in filepaths:
  infile = ROOT.TFile(myfile)
  tree = infile.Get("analysistree/anatree")

  for event in tree:
    n_evt += 1
    # we want to check if an event has at least one proton, muon and pi0
    n_p = 0
    n_mu = 0
    n_pi = 0
    for part in range(0,event.geant_list_size):
      if event.pdg[part] == 111:
        n_pi += 1
      if event.pdg[part] == 13:
        n_mu += 1
      if event.pdg[part] == 2212:
        n_p += 1
   
    if n_p and n_mu and n_pi:
      n_pass += 1
    else:
      continue


    for i in range(0,event.mcevts_truth):
      for j in range(0,event.ntracks_trackkalmanhit):
        if len(sys.argv) == 3:
	  if event.trkpdgtruth_trackkalmanhit[j] != mypdg:
            continue
        dx = event.trkstartx_trackkalmanhit[j] - event.nuvtxx_truth[i]
        dy = event.trkstarty_trackkalmanhit[j] - event.nuvtxy_truth[i]
        dz = event.trkstartz_trackkalmanhit[j] - event.nuvtxz_truth[i]
        hdisttovert_kalman.Fill(sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2)))

        vec_start = ROOT.TVector3(dx, dy, dz)
        vec_end = ROOT.TVector3(event.trkendx_trackkalmanhit[j] - event.nuvtxx_truth[i], event.trkendy_trackkalmanhit[j] - event.nuvtxy_truth[i], event.trkendz_trackkalmanhit[j] - event.nuvtxz_truth[i])
        start_end = ROOT.TVector3(event.trkstartx_trackkalmanhit[j] - event.trkendx_trackkalmanhit[j], event.trkstarty_trackkalmanhit[j] - event.trkendy_trackkalmanhit[j], event.trkstartz_trackkalmanhit[j] - event.trkendz_trackkalmanhit[j])
        proj = ( vec_start.Dot(start_end) / (start_end.Dot(start_end)) ) * start_end
        perp = proj - vec_start

        theta = asin( (vec_start.Cross(vec_end)).Mag() / (vec_start.Mag() * vec_end.Mag()) )
        alpha = asin( (vec_start.Cross(perp)).Mag() / (vec_start.Mag() * perp.Mag()) )
        beta = asin( (vec_end.Cross(perp)).Mag() / (vec_end.Mag() * perp.Mag()) )

        if (theta >= alpha) and (theta >= beta):
          clap = perp.Mag()
        else:
          clap = min(vec_start.Mag(), vec_end.Mag())

        hclosestapproach_kalman.Fill(clap)
        htrackangle_kalman.Fill(theta)

      for j in range(0,event.ntracks_pandoraNuKHit):
        if len(sys.argv) == 3:
          if event.trkpdgtruth_pandoraNuKHit[j] != mypdg:
            continue
        dx = event.trkstartx_pandoraNuKHit[j] - event.nuvtxx_truth[i]
        dy = event.trkstarty_pandoraNuKHit[j] - event.nuvtxy_truth[i]
        dz = event.trkstartz_pandoraNuKHit[j] - event.nuvtxz_truth[i]
        hdisttovert_pandora.Fill(sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2)))
        
        vec_start = ROOT.TVector3(dx, dy, dz)
        vec_end = ROOT.TVector3(event.trkendx_pandoraNuKHit[j] - event.nuvtxx_truth[i], event.trkendy_pandoraNuKHit[j] - event.nuvtxy_truth[i], event.trkendz_pandoraNuKHit[j] - event.nuvtxz_truth[i])
        start_end = ROOT.TVector3(event.trkstartx_pandoraNuKHit[j] - event.trkendx_pandoraNuKHit[j], event.trkstarty_pandoraNuKHit[j] - event.trkendy_pandoraNuKHit[j], event.trkstartz_pandoraNuKHit[j] - event.trkendz_pandoraNuKHit[j])
        proj = ( vec_start.Dot(start_end) / (start_end.Dot(start_end)) ) * start_end
        perp = proj - vec_start

        theta = asin( (vec_start.Cross(vec_end)).Mag() / (vec_start.Mag() * vec_end.Mag()) )
        alpha = asin( (vec_start.Cross(perp)).Mag() / (vec_start.Mag() * perp.Mag()) )
        beta = asin( (vec_end.Cross(perp)).Mag() / (vec_end.Mag() * perp.Mag()) )

        if (theta >= alpha) and (theta >= beta):
          clap = perp.Mag()
        else:
          clap = min(vec_start.Mag(), vec_end.Mag())

        hclosestapproach_pandora.Fill(clap)
        htrackangle_pandora.Fill(theta)

print n_pass
print n_evt


canv = ROOT.TCanvas()
hdisttovert_pandora.Draw()
canv.SaveAs("myplot_pandora"+suffix+".eps")
hdisttovert_kalman.Draw()
canv.SaveAs("myplot_kalman"+suffix+".eps")
hclosestapproach_pandora.Draw()
canv.SaveAs("clap_pandora"+suffix+".eps")
hclosestapproach_kalman.Draw()
canv.SaveAs("clap_kalman"+suffix+".eps")
htrackangle_pandora.Draw()
canv.SaveAs("theta_pandora"+suffix+".eps")
htrackangle_kalman.Draw()
canv.SaveAs("theta_kalman"+suffix+".eps")


