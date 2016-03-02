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


filelim = 200
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
    


    for i in range(0,event.mcevts_truth):
      for j in range(0,event.ntracks_trackkalmanhit):
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

      for j in range(0,event.ntracks_pandoraNuKHit):
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



canv = ROOT.TCanvas()
hdisttovert_pandora.Draw()
canv.SaveAs("myplot_pandora.eps")
hdisttovert_kalman.Draw()
canv.SaveAs("myplot_kalman.eps")
hclosestapproach_pandora.Draw()
canv.SaveAs("clap_pandora.eps")
hclosestapproach_kalman.Draw()
canv.SaveAs("clap_kalman.eps")


