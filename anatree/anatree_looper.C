/*
 * AnaTree Looper
 * author: Davio Cianci with Claire Savard
 *
 * This module is called by anatree_plotter.py because python is really slow and bad at looping through large numbers
 * of things.
 * loop() is called once, and it reads a .txt file containing a list of paths to anatree hist files. We then loop through
 * all events in each file and fill several histograms.
 * draw() is called once, and it prints out a brief summary and saves the plots to a folder called "plots"
 *
 * */

#include "anatree_looper.h"

TH1D* hdisttovert_kalman = new TH1D("Distance to Vertex","Distance to Vertex (Kalman); cm, presumably; entries",200,0,100);
TH1D* hdisttovert_pandora = new TH1D("Distance to Vertex","Distance to Vertex (Pandora); cm, presumably; entries",200,0,100);
TH1D* hclosestapproach_kalman = new TH1D("Closest approach","Closest Approach (Kalman); cm",200,0,100);
TH1D* hclosestapproach_pandora = new TH1D("Closest approach","Closest Approach (Pandora); cm",200,0,100);
TH1D* htrackangle_pandora = new TH1D("Track Angle","Track Angle (Pandora); radians", 200, 0, 6.5);
TH1D* htrackangle_kalman = new TH1D("Track Angle", "Track Angle (Kalman); radians", 200, 0, 6.5);
TH1F* htracklen_pandora = new TH1F("Track Length", "Track Length (Pandora); cm", 1000, 0, 500);
TH1F* htracklen_kalman = new TH1F("Track Length", "Track Length (Kalman); cm", 1000, 0, 500);

std::vector <std::string> paths;

#define kmax 25000
int pdg[kmax], status[kmax];
short ntracks_trackkalmanhit, ntracks_pandoraNuKHit;
int trkpdgtruth_trackkalmanhit[kmax], trkpdgtruth_pandoraNuKHit[kmax], mcevts_truth, geant_list_size;
int trkpidpdg_trackkalmanhit[kmax], trkpidpdg_pandoraNuKHit[kmax];
float trkstartx_trackkalmanhit[kmax], trkstarty_trackkalmanhit[kmax], trkstartz_trackkalmanhit[kmax], trkstartx_pandoraNuKHit[kmax], trkstarty_pandoraNuKHit[kmax], trkstartz_pandoraNuKHit[kmax];
float trkendx_trackkalmanhit[kmax], trkendy_trackkalmanhit[kmax], trkendz_trackkalmanhit[kmax], trkendx_pandoraNuKHit[kmax], trkendy_pandoraNuKHit[kmax], trkendz_pandoraNuKHit[kmax];
float nuvtxx_truth[10], nuvtxy_truth[10], nuvtxz_truth[10];
float trklen_trackkalmanhit[kmax], trklen_pandoraNuKHit[kmax];
bool vtx_truth = true;
bool yz_only = true;

TVector3 vec_start, vec_end, start_end, proj, perp;

std::string s_suffix = "";
int n_evt = 0;
int n_pass = 0;
int n_protons = 0; int n_protons_kalman = 0; int n_protons_pandora = 0;

/*if vtx_truth == true {
	nuvtxx = nuvtxx_truth;
	nuvtxy = nuvtxy_truth;
	nuvtxx = nuvtxz_truth;
} else {
	nuvtxx = 
	nuvtxy = 
	nuvtxz = 
}*/ 
	

void loop(int mypdg){
	if(mypdg == 13) s_suffix = "_mu";
	if(mypdg == 111) s_suffix = "_pi0";
	if(mypdg == 2212) s_suffix = "_p";

	// Read in the pathlists and fill vector
	paths.resize(0);
	std::string thisline;
	ifstream infile;
	infile.open("pathlist.txt");
	while(std::getline(infile,thisline)){
		paths.push_back(thisline);
	}
	infile.close();


	for(unsigned int fs = 0; fs < paths.size(); fs++){
		if(fs % 25 == 0) std::cout << fs << "/" << paths.size() << " files scanned" << std::endl;

		TFile* infile = new TFile(paths[fs].c_str());
		TTree* tree = (TTree*)(infile->Get("analysistree/anatree"));

		// set branches...
		tree->SetBranchAddress("geant_list_size", &geant_list_size);
		tree->SetBranchAddress("pdg", pdg);
		tree->SetBranchAddress("status",status);
		tree->SetBranchAddress("mcevts_truth", &mcevts_truth);
		tree->SetBranchAddress("ntracks_trackkalmanhit", &ntracks_trackkalmanhit);
		tree->SetBranchAddress("ntracks_pandoraNuKHit", &ntracks_pandoraNuKHit);
		tree->SetBranchAddress("trkpdgtruth_trackkalmanhit", trkpdgtruth_trackkalmanhit);
		tree->SetBranchAddress("trkpdgtruth_pandoraNuKHit", trkpdgtruth_pandoraNuKHit);
		tree->SetBranchAddress("trkstartx_trackkalmanhit", trkstartx_trackkalmanhit);
		tree->SetBranchAddress("trkstarty_trackkalmanhit", trkstarty_trackkalmanhit);
		tree->SetBranchAddress("trkstartz_trackkalmanhit", trkstartz_trackkalmanhit);
		tree->SetBranchAddress("trkstartx_pandoraNuKHit", trkstartx_pandoraNuKHit);
		tree->SetBranchAddress("trkstarty_pandoraNuKHit", trkstarty_pandoraNuKHit);
		tree->SetBranchAddress("trkstartz_pandoraNuKHit", trkstartz_pandoraNuKHit);
		tree->SetBranchAddress("trkendx_trackkalmanhit", trkendx_trackkalmanhit);
		tree->SetBranchAddress("trkendy_trackkalmanhit", trkendy_trackkalmanhit);
		tree->SetBranchAddress("trkendz_trackkalmanhit", trkendz_trackkalmanhit);
		tree->SetBranchAddress("trkendx_pandoraNuKHit", trkendx_pandoraNuKHit);
		tree->SetBranchAddress("trkendy_pandoraNuKHit", trkendy_pandoraNuKHit);
		tree->SetBranchAddress("trkendz_pandoraNuKHit", trkendz_pandoraNuKHit);
		tree->SetBranchAddress("nuvtxx_truth", nuvtxx_truth);
		tree->SetBranchAddress("nuvtxy_truth", nuvtxy_truth);
		tree->SetBranchAddress("nuvtxz_truth", nuvtxz_truth);
		tree->SetBranchAddress("trklen_trackkalmanhit", trklen_trackkalmanhit);
		tree->SetBranchAddress("trklen_pandoraNuKHit", trklen_pandoraNuKHit);
		tree->SetBranchAddress("trkpidpdg_pandoraNuKHit", trkpidpdg_pandoraNuKHit);
		tree->SetBranchAddress("trkpidpdg_trackkalmanhit", trkpidpdg_trackkalmanhit);

		for(int g = 0; g < tree->GetEntries(); g++){
			tree->GetEntry(g);
			n_evt ++;

			// Our filtah
			int n_p = 0; int n_mu = 0; int n_pi = 0;
			for(int part = 0; part < geant_list_size; part++){
				if(status[part] != 1) continue;
				if(pdg[part] == 111) 	n_pi++;
				if(pdg[part] == 13)	n_mu++;
				if(pdg[part] == 2212) 	n_p++;
			}
			if(n_p && n_pi && n_mu)
				n_pass ++;
				else
				continue;
			n_protons += n_p;

	    	for(int i = 0; i < mcevts_truth; i++){
	     		for(int j = 0; j < ntracks_trackkalmanhit; j++){
					if(trkpidpdg_trackkalmanhit[j] == 2212)  n_protons_kalman++;
					if(mypdg != 0)
					  if( trkpidpdg_trackkalmanhit[j] != mypdg)
					    continue;

				double clap;
				double dx_start = 0;
				double dx_end = 0;

				if (yz_only == false) {
					dx_end = trkendx_trackkalmanhit[j] - nuvtxx_truth[i];
					dx_start = trkstartx_trackkalmanhit[j] - nuvtxx_truth[i];
				}
					
				double dy_start = trkstarty_trackkalmanhit[j] - nuvtxy_truth[i];
	        		double dz_start = trkstartz_trackkalmanhit[j] - nuvtxz_truth[i];
				double dy_end = trkendy_trackkalmanhit[j] - nuvtxy_truth[i];
				double dz_end = trkendz_trackkalmanhit[j] - nuvtxz_truth[i];
	        		hdisttovert_kalman->Fill(sqrt(pow(dx_start,2)+pow(dy_start,2)+pow(dz_start,2)));

	        		vec_start =  TVector3(dx_start, dy_start, dz_start);
	        		vec_end =  TVector3(dx_end, dy_end, dz_end);
	        		start_end =  TVector3(trkstartx_trackkalmanhit[j] - trkendx_trackkalmanhit[j], trkstarty_trackkalmanhit[j] - trkendy_trackkalmanhit[j], trkstartz_trackkalmanhit[j] - trkendz_trackkalmanhit[j]);
	        		proj = ( vec_start.Dot(start_end) / (start_end.Dot(start_end)) ) * start_end;
	        		perp = proj - vec_start;

	        		double theta = asin( (vec_start.Cross(vec_end)).Mag() / (vec_start.Mag() * vec_end.Mag()) );
	        		double alpha = asin( (vec_start.Cross(perp)).Mag() / (vec_start.Mag() * perp.Mag()) );
	        		double beta = asin( (vec_end.Cross(perp)).Mag() / (vec_end.Mag() * perp.Mag()) );

	        		if (theta >= alpha && theta >= beta)
	        			clap = perp.Mag();
	        		else
	          			clap = min(vec_start.Mag(), vec_end.Mag());

	        		hclosestapproach_kalman->Fill(clap);
	        		htrackangle_kalman->Fill(theta);
				htracklen_kalman->Fill(trklen_trackkalmanhit[j]);
				}
	      		for(int j = 0; j < ntracks_pandoraNuKHit; j++){
					if(trkpidpdg_pandoraNuKHit[j] == 2212) n_protons_pandora++;
					if(mypdg != 0)
					  if( trkpidpdg_pandoraNuKHit[j] != mypdg)
					    continue;

				double clap;
				double dx_start = 0;
				double dx_end = 0;
				if (yz_only == false) {
					dx_end = trkendx_pandoraNuKHit[j] - nuvtxx_truth[i];
					dx_start = trkstartx_pandoraNuKHit[j] - nuvtxx_truth[i];
				}
					
				double dy_start = trkstarty_pandoraNuKHit[j] - nuvtxy_truth[i];
	        		double dz_start = trkstartz_pandoraNuKHit[j] - nuvtxz_truth[i];
				double dy_end = trkendy_pandoraNuKHit[j] - nuvtxy_truth[i];
				double dz_end = trkendz_pandoraNuKHit[j] - nuvtxz_truth[i];
	        		double dy = trkstarty_pandoraNuKHit[j] - nuvtxy_truth[i];
	        		double dz = trkstartz_pandoraNuKHit[j] - nuvtxz_truth[i];
	        		hdisttovert_pandora->Fill(sqrt(pow(dx_start,2)+pow(dy_start,2)+pow(dz_start,2)));

	        		vec_start = TVector3(dx_start, dy_start, dz_start);
	        		vec_end = TVector3(dx_end, dy_end, dz_end);
	        		start_end = TVector3(trkstartx_pandoraNuKHit[j] - trkendx_pandoraNuKHit[j], trkstarty_pandoraNuKHit[j] - trkendy_pandoraNuKHit[j], trkstartz_pandoraNuKHit[j] - trkendz_pandoraNuKHit[j]);
	        		proj = ( vec_start.Dot(start_end) / (start_end.Dot(start_end)) ) * start_end;
	        		perp = proj - vec_start;

	        		double theta = asin( (vec_start.Cross(vec_end)).Mag() / (vec_start.Mag() * vec_end.Mag()) );
	        		double alpha = asin( (vec_start.Cross(perp)).Mag() / (vec_start.Mag() * perp.Mag()) );
	        		double beta = asin( (vec_end.Cross(perp)).Mag() / (vec_end.Mag() * perp.Mag()) );

	        		if (theta >= alpha && theta >= beta)
	          			clap = perp.Mag();
	        		else
	          			clap = min(vec_start.Mag(), vec_end.Mag());

	        		hclosestapproach_pandora->Fill(clap);
	        		htrackangle_pandora->Fill(theta);
				htracklen_pandora->Fill(trklen_pandoraNuKHit[j]);
				}
			}
		}
		infile->Close();
	}
}

bool draw(){

	std::cout << n_pass << "/" << n_evt << " events passed our lovely little filter." << std::endl;
	std::cout << "We had " << n_protons << " protons, of which Kalman reconstructed " << n_protons_kalman << " tracks and pandora did " << n_protons_pandora << " tracks." << std::endl;
	std::cout << "Pandora reco efficiency: " << float(n_protons_pandora)/float(n_protons) << " Kalman reco efficiency: " << float(n_protons_kalman)/float(n_protons) << "\n\n\n\n\n\n\n" << std::endl;

	TCanvas* canv = new TCanvas();
	hdisttovert_pandora->Draw();
	canv->SaveAs(("plots/myplot_pandora"+s_suffix+".eps").c_str());
	hdisttovert_kalman->Draw();
	canv->SaveAs(("plots/myplot_kalman"+s_suffix+".eps").c_str());
	hclosestapproach_pandora->Draw();
	canv->SaveAs(("plots/clap_pandora"+s_suffix+".eps").c_str());
	hclosestapproach_kalman->Draw();
	canv->SaveAs(("plots/clap_kalman"+s_suffix+".eps").c_str());
	htrackangle_pandora->Draw();
	canv->SaveAs(("plots/theta_pandora"+s_suffix+".eps").c_str());
	htrackangle_kalman->Draw();
	canv->SaveAs(("plots/theta_kalman"+s_suffix+".eps").c_str());
	htracklen_kalman->Draw();
	canv->SaveAs(("plots/trklen_kalman"+s_suffix+".eps").c_str());
	htracklen_pandora->Draw();
	canv->SaveAs(("plots/trklen_pandora"+s_suffix+".eps").c_str());

	return false;
}
