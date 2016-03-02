#define short_track_study_cxx
#include "short_track_study.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TProfile.h>
#include "TF1.h"
#include "TF2.h"
#include <stdio.h>
#include <stdlib.h>
#include "algorithm"
#include <iostream>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include "TLegend.h"
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TGraphErrors.h>

void short_track_study::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout<<nentries<<"\n";
	Int_t j=0;
	Long64_t fltrdEvents[100000]={-999};
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      Int_t mu=0,prot=0,pi0=0;
      for(Int_t i=0;i<geant_list_size;i++) if(pdg[i]==13) mu++;
      for(Int_t i=0;i<geant_list_size;i++) if(pdg[i]==2212) prot++;
      for(Int_t i=0;i<geant_list_size;i++) if(pdg[i]==111) pi0++;
      if (mu>=1 && prot>=1 && pi0>=1) {
      fltrdEvents[j]=jentry; 
      std::cout<<"filtered event num is: "<<fltrdEvents[j]<<"\n";
      j++;}
   }
   for(Int_t k=0;k<j;k++){
   Long64_t kentry = LoadTree(fltrdEvents[k]);
   if (kentry < 0) break;
   std::cout<<"opening filtered event "<<fltrdEvents[k]<<"  "<<kentry<<"\n";
   
    for(Int_t n = 0; n < mcevts_truth && n < 10; n++) {
    std::cout<<"nuvtxx_truth "<<nuvtxx_truth<<"nuvtxy_truth "<<nuvtxy_truth<<"nuvtxz_truth "<<nuvtxz_truth<<"\n";
							}
}
}
