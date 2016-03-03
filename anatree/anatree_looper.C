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

TH1D* hdisttovert_kalman = new TH1D("Distance to Vertex K","Distance to Vertex (Kalman); cm, presumably; entries",200,0,100);
TH1D* hdisttovert_pandora = new TH1D("Distance to Vertex P","Distance to Vertex (Pandora); cm, presumably; entries",200,0,100);
TH1D* hclosestapproach_kalman = new TH1D("Closest approach K","Closest Approach (Kalman); cm",200,0,100);
TH1D* hclosestapproach_pandora = new TH1D("Closest approach P","Closest Approach (Pandora); cm",200,0,100);
//TH1D* hzoomclap_kalman = new TH1D("
TH1D* htrackangle_pandora = new TH1D("Track Angle P","Track Angle (Pandora); radians", 200, 0, 6.5);
TH1D* htrackangle_kalman = new TH1D("Track Angle K", "Track Angle (Kalman); radians", 200, 0, 6.5);
TH1F* htracklen_pandora = new TH1F("Track Length P", "Track Length (Pandora); cm", 1000, 0, 500);
TH1F* htracklen_kalman = new TH1F("Track Length K", "Track Length (Kalman); cm", 1000, 0, 500);
TH1F* htracklenshort_pandora = new TH1F("Short Track Length P", "Short Track Length (Pandora); cm", 50, 0, 50);
TH1F* htracklenshort_kalman = new TH1F("Short Track Length K", "Short Track Length (Kalman); cm", 50, 0, 50);
TH1D* hprotondistance = new TH1D("Proton Distance","Proton Distance;cm;entries",200,0,10);

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
float StartPointx[kmax], StartPointy[kmax], StartPointz[kmax], EndPointx[kmax], EndPointy[kmax], EndPointz[kmax];
bool vtx_truth = true;
bool yz_only = true;

TVector3 vec_start, vec_end, start_end, proj, perp;

std::string s_suffix = "";
int n_evt = 0;
int n_pass = 0;
int n_protons = 0; int n_protons_kalman = 0; int n_protons_pandora = 0;
long cnttruepdg_kalman, cnttruepdg_pandora;
int cnttrashpdg_kalman, cnttrashpdg_pandora;
int cntbothval_kalman, cntbothval_pandora;
int cntdiffer_kalman, cntdiffer_pandora;
float pdist = 0.;

/*if vtx_truth == true {
	nuvtxx = nuvtxx_truth;
	nuvtxy = nuvtxy_truth;
	nuvtxx = nuvtxz_truth;
} else {
	nuvtxx =
	nuvtxy =
	nuvtxz =
}*/

bool closeEnough(std::vector<TVector3> pis, std::vector<TVector3> mus, std::vector<TVector3> ps){
	float buffer = .5;

	for(unsigned int i = 0; i < pis.size(); i++){
		for(unsigned int j = 0; j < pis.size(); j++){
			if((pis[i]- mus[j]).Mag() < buffer){
				for(unsigned int k = 0; k < ps.size(); k++){
					if((pis[i] - ps[k]).Mag() < buffer && (mus[j] - ps[k]).Mag() < buffer)
						return true;
				}
			}
		}
	}
	return false;
}

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
		tree->SetBranchAddress("StartPointx", StartPointx);
		tree->SetBranchAddress("StartPointy", StartPointy);
		tree->SetBranchAddress("StartPointz", StartPointz);
		tree->SetBranchAddress("EndPointx", EndPointx);
		tree->SetBranchAddress("EndPointy", EndPointy);
		tree->SetBranchAddress("EndPointz", EndPointz);

		for(int g = 0; g < tree->GetEntries(); g++){
			tree->GetEntry(g);
			n_evt ++;

			// Our filtah
			std::vector<TVector3> pis; std::vector<TVector3> mus; std::vector<TVector3> ps;
			int n_p = 0; int n_pi = 0; int n_mu = 0;
			for(int part = 0; part < geant_list_size; part++){

				if(status[part] != 1)
					continue;
				if(pdg[part] == 2212){
					pdist = sqrt(pow(StartPointx[part] - EndPointx[part],2) + pow(StartPointy[part] - EndPointy[part],2) + pow(StartPointz[part] - EndPointz[part],2));
					hprotondistance->Fill(pdist);
					if(pdist < 1.) continue;
				}
				if(pdg[part] == 111){
					pis.push_back(TVector3(StartPointx[part],StartPointy[part],StartPointz[part]));
					n_pi ++;
				}
				if(pdg[part] == 13){
					mus.push_back(TVector3(StartPointx[part],StartPointy[part],StartPointz[part]));
					n_mu ++;
				}
				if(pdg[part] == 2212){
					ps.push_back(TVector3(StartPointx[part],StartPointy[part],StartPointz[part]));
					n_p ++;
				}
			}
			if(n_p && n_mu && n_pi){
				// Now, let's loop through
				if(closeEnough(pis,mus,ps))
					n_pass++;
				else
					continue;
			}
			else continue;

			n_protons += n_p;
			//std::cout << "           n_protons: " << n_protons << std::endl;

	    	for(int i = 0; i < mcevts_truth; i++){
	     		for(int j = 0; j < ntracks_trackkalmanhit; j++){
					if(trkpidpdg_trackkalmanhit[j] == 2212)  n_protons_kalman++;
					if(mypdg != 0)
					  if( trkpidpdg_trackkalmanhit[j] != mypdg)
					    continue;

				//std::cout << "---------------Kalman---------------------" << std::endl;
				//std::cout << "true pdg: " << trkpdgtruth_trackkalmanhit[j] << "       pid pdg: " << trkpidpdg_trackkalmanhit[j] << std::endl;
				cnttruepdg_kalman++;
				if (trkpdgtruth_trackkalmanhit[j] == -99999) {
				      cnttrashpdg_kalman++;
				}
				else if (trkpdgtruth_trackkalmanhit[j] != -99999) {
				      cntbothval_kalman++;
				      if (abs(trkpdgtruth_trackkalmanhit[j]) != trkpidpdg_trackkalmanhit[j]) {
					    cntdiffer_kalman++;
				      }
				}

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
				double angfromvtx = asin( (vec_start.Cross(-1*start_end)).Mag() / (vec_start.Mag() * start_end.Mag()) );

	        		if (theta >= alpha && theta >= beta)
	        			clap = perp.Mag();
	        		else
	          			clap = min(vec_start.Mag(), vec_end.Mag());

	        		hclosestapproach_kalman->Fill(clap);
	        		htrackangle_kalman->Fill(angfromvtx);
				htracklen_kalman->Fill(trklen_trackkalmanhit[j]);
				htracklenshort_kalman->Fill(trklen_trackkalmanhit[j]);
				}
	      		for(int j = 0; j < ntracks_pandoraNuKHit; j++){
					if(trkpidpdg_pandoraNuKHit[j] == 2212) n_protons_pandora++;
					if(mypdg != 0)
					  if( trkpidpdg_pandoraNuKHit[j] != mypdg)
					    continue;
				//std::cout << "---------------Pandora---------------------" << std::endl;
				//std::cout << "true pdg: " << trkpdgtruth_pandoraNuKHit[j] << "       pid pdg: " << trkpidpdg_pandoraNuKHit[j] << std::endl;
				cnttruepdg_pandora++;
				if (trkpdgtruth_pandoraNuKHit[j] == -99999) {
				      cnttrashpdg_pandora++;
				}
				else if (trkpdgtruth_pandoraNuKHit[j] != -99999) {
				      cntbothval_pandora++;
				      if (abs(trkpdgtruth_pandoraNuKHit[j]) != trkpidpdg_pandoraNuKHit[j]) {
					    cntdiffer_pandora++;
				      }
				}

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
				double angfromvtx = asin( (vec_start.Cross(-1*start_end)).Mag() / (vec_start.Mag() * start_end.Mag()) );

	        		if (theta >= alpha && theta >= beta)
	          			clap = perp.Mag();
	        		else
	          			clap = min(vec_start.Mag(), vec_end.Mag());

	        		hclosestapproach_pandora->Fill(clap);
	        		htrackangle_pandora->Fill(angfromvtx);
				htracklen_pandora->Fill(trklen_pandoraNuKHit[j]);
				htracklenshort_pandora->Fill(trklen_pandoraNuKHit[j]);
				}
			}
			/*std::cout << "---Kalman---" << std::endl;
			std::cout << "Total amount of particles: " << cnttruepdg_kalman << std::endl;
			std::cout << "Number of trash true pdg's: " << cnttrashpdg_kalman << "        Number of not trash pdg's: " << cntbothval_kalman << std::endl;
			std::cout << "Number of times true and pid pdg differed when true pdg is not trash: " << cntdiffer_kalman << std::endl << std::endl;
			std::cout << "---Pandora---" << std::endl;
			std::cout << "Total amount of particles: " << cnttruepdg_pandora << std::endl;
			std::cout << "Number of trash true pdg's: " << cnttrashpdg_pandora << "        Number of not trash pdg's: " << cntbothval_pandora << std::endl;
			std::cout << "Number of times true and pid pdg differed when true pdg is not trash: " << cntdiffer_pandora << std::endl;
			std::cout << "\n\n\n\n";*/
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
	canv->SaveAs(("plots/myplot_pandora"+s_suffix+".pdf").c_str());
	hdisttovert_kalman->Draw();
	canv->SaveAs(("plots/myplot_kalman"+s_suffix+".pdf").c_str());
	hclosestapproach_pandora->Draw();
	canv->SaveAs(("plots/clap_pandora"+s_suffix+".pdf").c_str());
	hclosestapproach_kalman->Draw();
	canv->SaveAs(("plots/clap_kalman"+s_suffix+".pdf").c_str());
	canv->SetLogy(true);
	htrackangle_pandora->Draw();
	canv->SaveAs(("plots/theta_pandora"+s_suffix+".pdf").c_str());
	htrackangle_kalman->Draw();
	canv->SetLogy(false);
	canv->SaveAs(("plots/theta_kalman"+s_suffix+".pdf").c_str());
	htracklen_kalman->Draw();
	canv->SaveAs(("plots/trklen_kalman"+s_suffix+".pdf").c_str());
	htracklen_pandora->Draw();
	canv->SaveAs(("plots/trklen_pandora"+s_suffix+".eps").c_str());
	htracklenshort_kalman->Draw();
	canv->SaveAs(("plots/trklenshort_kalman"+s_suffix+".pdf").c_str());
	htracklenshort_pandora->Draw();
	canv->SaveAs(("plots/trklenshort_pandora"+s_suffix+".pdf").c_str());
	hprotondistance->Draw();
	canv->SaveAs(("plots/protondist"+s_suffix+".eps").c_str());

	return false;
}
