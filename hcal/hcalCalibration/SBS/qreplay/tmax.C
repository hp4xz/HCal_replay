//SSeeds 9.1.23 Applies ADCt corrections and builds agg histo for sbs.hcal.adc.tmax calibration

#include <vector>
#include <iostream>
#include <algorithm>
#include <vector>
#include <filesystem>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

const int total_bins = 900;
const int lower_lim = -150;
const int upper_lim = 150;
const int total_bins_short = 90;
const int lower_lim_short = -15;
const int upper_lim_short = 15;

const double xmin = -10.;
const double xmax = 10.;

const double aFWHM = 5; //full width half max approximate adct cluster1 - cluster2 dist

double skewedGaussPlusConstant(double *x, double *p) {
    double gaussPart = p[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - p[1])/p[2], 2));
    double skewPart = 1 + TMath::Erf(p[3] * (x[0] - p[1]) / p[2]);
    return gaussPart * skewPart + p[4];
}

double computeHWHMLeft(TF1 *fitFunc, double x_max, double max) {
    double halfMax = max / 2.0;
    double a = x_max;
    double b = x_max - aFWHM;

    while (fabs(fitFunc->Eval((a + b) / 2) - halfMax) > 1e-6) {  // 1e-6 is precision, can be adjusted
        double c = (a + b) / 2;
        if (fitFunc->Eval(c) > halfMax) {
            a = c;
        } else {
            b = c;
        }
    }
    return x_max - (a + b) / 2;  // Difference in x
}

double computeHWHMLeftZero(TF1 *fitFunc, double x_max, double max) {
    double halfMax = max / 2.0;
    double a = x_max;
    double b = x_max - aFWHM;

    while (fabs(fitFunc->Eval((a + b) / 2) - halfMax) > 1e-6) {  // 1e-6 is precision, can be adjusted
        double c = (a + b) / 2;
        if (fitFunc->Eval(c) > halfMax) {
            a = c;
        } else {
            b = c;
        }
    }
   return 0.0 - (a + b) / 2;  // Difference from x=0
}

double computeHWHMRightZero(TF1 *fitFunc, double x_max, double max) {
    double halfMax = max / 2.0;
    double a = x_max;
    double b = x_max + aFWHM;  // Arbitrary; adjust based on your data

    while (fabs(fitFunc->Eval((a + b) / 2) - halfMax) > 1e-6) {
        double c = (a + b) / 2;
        if (fitFunc->Eval(c) > halfMax) {  // This is the corrected condition
            a = c;
        } else {
            b = c;
        }
    }
    return (a + b) / 2 - 0.0;  // Difference from x=0
}

std::pair<double, double> computeMaxValue(TF1 *fitFunc, double mean, double range) {
    double step = 0.1;  // Adjust this as per the histogram resolution
    double x_start = mean - range;
    double x_end = mean + range;

    double maxVal = fitFunc->Eval(mean);
    double x_max = mean;
    for (double x = x_start; x <= x_end; x += step) {
        double val = fitFunc->Eval(x);
        if (val > maxVal) {
            maxVal = val;
            x_max = x;
        }
    }
    return std::make_pair(maxVal, x_max);
}

void overlayWithGaussianFits(TFile *file1, TH2D* hist, TCanvas* canvas, TH1D* hist2) {
  if (!hist || !canvas || !hist2) {
    std::cerr << "Null histogram or canvas pointer!" << std::endl;
    return;
  }

  // Use the provided canvas
  canvas->cd();

  // Create TGraphErrors for standard and red points
  TGraphErrors* graph = new TGraphErrors();

  for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
    // Project the 2D histogram onto a 1D histogram along the y-axis for the current x-bin
    TH1D* projY = hist->ProjectionY("_py", binX, binX);

    //Skip the fit if less than 100 entries exist in the bin
    if( projY->GetEntries()<100 )
      continue;

    //hist2 = projY;
	
    hist2 = (TH1D*)projY->Clone(Form("hist2_%d",binX));
    file1->cd();
    hist2->Write();

    //declare some dynamic fit variables
    Double_t FWHM = aFWHM;
    Int_t binMax = projY->GetMaximumBin();
    Double_t binCenter = projY->GetBinCenter( binMax );
    Double_t fitLowerLim = binCenter - FWHM;
    Double_t fitUpperLim = binCenter + FWHM;
	
    // Fit a Gaussian to this projection
    TF1 *gaussFit = new TF1("gausFit", "gaus");
    projY->Fit(gaussFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode
    
    if (gaussFit) { // Ensure the fit was successful

      // Determine the x-value as the center of the current bin
      double xCenter = hist->GetXaxis()->GetBinCenter(binX);

      double mean = gaussFit->GetParameter(1);
      double sigma = gaussFit->GetParameter(2);

      // Get the errors on the mean and sigma
      double meanError = gaussFit->GetParError(1);
      double sigmaError = gaussFit->GetParError(2);

      // Fill the graph
      Int_t binXval = hist->GetXaxis()->GetBinCenter(binX);
      int pointIdx = graph->GetN();
      graph->SetPoint(pointIdx, xCenter, mean);
      graph->SetPointError(pointIdx, 0, sigma);
    }
    delete projY; // Clean up
  }

  // Draw the original histogram and standard graph
  hist->Draw("COLZ");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("P SAME");

  // Update the canvas
  canvas->Update();

  //Add a legend
  auto legend = new TLegend(0.63,0.7,0.89,0.89);
  legend->SetTextSize(0.03);
  legend->AddEntry(graph,"Gaus fit, mean and sigma","l");
  legend->Draw();

  // Update the canvas
  canvas->Update();
}

//for monitoring - sends current coefficient set to console
void sendOffsetsToConsole(const double arr[hcal::maxHCalChan], std::string type_var) {
  std::cout << std::endl << Form( "%s =",type_var.c_str() ) << std::endl;
  for( int r=0; r<hcal::maxHCalRows; r++ ){
    for ( int c=0; c<hcal::maxHCalCols; c++){
      int i = r*hcal::maxHCalCols+c;
      std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

//Main (select indicates to the script which cluster block comparison is the highest with sufficient stats)
void tmax( std::string experiment = "gmn", 
	   int config = 4, 
	   int pass = 0, 
	   int set = 1,
	   int max_runs = 10,
	   int run_begin = 0,
	   int run_end = 0,
	   int select = 6,
	   int nsig = 4 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  // Get the path to the new offsets
  std::string patch_path;
  if( config == 4 && set == 1 )
    patch_path = "../timing/parameters/adctsetoffsets_gmn_conf4_pass0_0_to_0_exclude_0_to_0.txt";
  else if( config == 7 && set == 1 )
    patch_path = "../timing/parameters/adctsetoffsets_gmn_conf7_pass0_0_to_0_exclude_0_to_0.txt";
  else if( config == 14 && set == 1 )
    patch_path = "../timing/parameters/adctsetoffsets_gmn_conf14_pass1_13239_to_13260_exclude_0_to_0.txt";
  else if( config == 14 && set == 2 )
    patch_path = "../timing/parameters/adctsetoffsets_gmn_conf14_pass1_13261_to_13407_exclude_0_to_0.txt";
  else if( config == 11 && set == 1 )
    patch_path = "../timing/parameters/adctsetoffsets_gmn_conf11_pass1_0_to_0_exclude_12450_to_12860.txt";
  else if( config == 11 && set == 2 )
    patch_path = "../timing/parameters/adctsetoffsets_gmn_conf11_pass1_12450_to_12860_exclude_0_to_0.txt";
  else if( config == 9 && set == 1 )
    patch_path = "../timing/parameters/adctsetoffsets_gmn_conf9_pass1_0_to_0_exclude_0_to_0.txt";
  else if( config == 8 && set == 1 )
    patch_path = "../timing/parameters/adctsetoffsets_gmn_conf8_pass1_0_to_0_exclude_0_to_0.txt";
  else{
    std::cerr << "ERROR: config/set combination does not exist. Reconfigure and rerun." << std::endl;
    return;
  }

  // Get the date
  string date = util::getDate();

  //Set up path variables and output files
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + Form("/hcal_calibrations/qreplay/hcalsolo_adct_%s_sbs%d_set%d.root",experiment.c_str(),config,set);
  //std::string fout_path = "fit_test.root";

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );
  
  // Gets database directory from system, must reconfigure if not set
  // std::string db_path = gSystem->Getenv("DB_DIR");
  // std::string hcal_db_path = db_path + "/db_sbs.hcal.dat";
  std::string hcal_db_path = "/w/halla-scshelf2102/sbs/seeds/alt_sbsreplay/SBS-replay/DB/db_sbs.hcal.dat";

  //get database variable name
  std::string db_variable = "sbs.hcal.adc.timeoffset";
  
  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment.c_str()); //unique to my environment for now
  int nruns = -1; //Analyze all available runs for this configuration
  int verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);   
  std::cout << "Configuration parameters for data to qreplay:" << std::endl;
  std::cout << config_parameters;
  double ebeam = config_parameters.GetEbeam();
  double hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  double hcaldist = config_parameters.GetHCALdist();
  double sbsdist = config_parameters.GetSBSdist();
  double bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  //Histograms
  TH2D *ht_run = new TH2D("ht_run",
			  "adct vs run",
			  (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
			  runs[0].runnum-1,
			  runs[nruns-1].runnum+1,
			  total_bins,
			  lower_lim,
			  upper_lim);

  TH2D *ht_id = new TH2D("ht_id",
			 "adct vs id",
			 hcal::maxHCalChan,
			 0,
			 hcal::maxHCalChan,
			 total_bins,
			 lower_lim,
			 upper_lim);

  TH2D *ht_cidx = new TH2D("ht_cidx",
			   "adct diff vs cluster block index; index; ns",
			   20,
			   0,
			   20,
			   80,
			   -10,
			   10);

  TH1D *ht_allchan_allruns = new TH1D("ht_allchan_allruns",
				      "adct all channels/runs",
				      total_bins,
				      lower_lim,
				      upper_lim);

  TH1D *ht_all_diff = new TH1D("ht_all_diff",
			       "adct all channels/runs highE cluster block difference",
			       total_bins_short,
			       lower_lim_short,
			       upper_lim_short);


  TH1D *ht_all_best_diff = new TH1D("ht_all_best_diff",
				    "adct all channels/runs highE cluster best block difference",
				    total_bins_short,
				    lower_lim_short,
				    upper_lim_short);


  TH1D *ht_all_best_diff_idx = new TH1D("ht_all_best_diff_idx",
					"adct all channels/runs highE cluster best block difference index",
					total_bins_short,
					lower_lim_short,
					upper_lim_short);


  TH1D *ht_all_select_diff = new TH1D("ht_all_select_diff",
  				    "adct all channels/runs highE cluster select block difference",
  				    total_bins_short,
  				    lower_lim_short,
  				    upper_lim_short);

  //reporting/run indices
  int target_change_index;
  double config_sampling_fraction;
  double config_e_sigma_ratio;
  std::string ts_compare = "";
  bool cut_first = true;
  vector<int> used_runs;
  bool first = false;
  int good_run = 0;

  //loop over runs and exclude. should update util::ReadRunList to improve efficiency in the future
  for (int r=0; r<nruns; r++) {

    if( good_run > max_runs ) continue;

    int current_runnumber = runs[r].runnum;
    
    //deal with user-specified range only
    if( run_begin == 0 && run_end ==0 ){
      cout << "All runs selected. Proceeding with run " << current_runnumber << endl;
    }else if( current_runnumber > run_end ||
	      current_runnumber < run_begin )
      continue;

    cout << "run " << r << "/" << nruns << ": " << current_runnumber << endl;      
    
    bool ts_different = false;
    
    //Get run experimental parameters
    std::string current_timestamp = runs[r].adcg_ts; //For now, only select on energy calibration sets
    std::string current_target = runs[r].target;
    
    std::cout << "Current run number: " << current_runnumber << ", current timestamp: " << current_timestamp << ", current target: " << current_target << std::endl;
    
    if( !(current_timestamp.compare(ts_compare)==0) ){
      ts_compare = current_timestamp;
      ts_different = true;
    }
        
    
    std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    int mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A
    
    //Get run paths
    std::string rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-%s/pass%d/SBS%d/%s/rootfiles/",experiment.c_str(),pass,config,targ_uppercase.c_str());
    std::string rootfile_path = rootfile_dir + Form("*%d*",current_runnumber);

    //Get target configuration
    SBStarget target_parameters(current_target);
    int target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
    double target_length = target_parameters.GetTargLength();
    double target_rho = target_parameters.GetTargRho();
    double cell_rho = target_parameters.GetCellRho();
    double cell_diam = target_parameters.GetCellDiam();
    double cell_dEdx = target_parameters.GetCelldEdx();
    double upstream_wthick = target_parameters.GetUpstreamWallThick();
    double target_dEdx = target_parameters.GetTargdEdx();
    double M_avg = target_parameters.GetAvgMass();

    //Record old gain and offset parameters by tstamp from database (assumes one file for all timestamps)
    double old_time_offsets[hcal::maxHCalChan] = {0.};

    //switch to correct timestamp where necessary
    std::string time_timestamp = runs[r].adct_ts;

    std::cout << "Reading parameters from file " << hcal_db_path << ":" << std::endl;
    util::readDB( hcal_db_path, time_timestamp, db_variable, old_time_offsets );

    //get new adct offsets from calibrated parameters
    double new_time_offsets[hcal::maxHCalChan] = {0.};

    //if no selected timestamp is passed, we use the first set of new coefficients that appears
    string quality_timestamp = "";
    if( quality_timestamp.compare("")==0 ) quality_timestamp = "none";

    //read in the new offsets. assumes first set in file at active_patch_path is the one desired
    std::cout << "Reading parameters from file: " << patch_path << std::endl;
    util::readDB( patch_path, "none", db_variable, new_time_offsets );

    //Report offset parameters on first loop
    if( first ){
      std::cout << "Loaded set of OLD adct offset coefficients:" << std::endl;
      sendOffsetsToConsole( old_time_offsets, db_variable );
      
      std::cout << "Loaded set of NEW adct offset coefficients:" << std::endl;
      sendOffsetsToConsole( new_time_offsets, db_variable );
      
      first = false;
    }

    //Get available cuts for current config/target/field combination. Set index = 0 as no calibration is being done here and physics cuts do not vary by calibration set
    vector<calcut> cut;
    util::ReadCutList(struct_dir,experiment,config,0,pass,current_target,mag,verb,cut);
    if( cut_first ){
      std::cout << cut[0];
      cut_first = false;
    }

    // Setting up chain and branch addresses
    C = new TChain("T");
    C->Add(rootfile_path.c_str());

    C->SetBranchStatus("*",0);    
    double BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
    double BBtr_vz[hcal::maxTracks];
    double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;	
    double HCALx, HCALy, HCALe;
    double pblkrow, pblkcol, nblk, nclus;
    int Ncid, Ncblkid;
    double cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan], cblkagain[hcal::maxHCalChan];
    double cid[hcal::maxHCalClus], crow[hcal::maxHCalClus], ccol[hcal::maxHCalClus], ce[hcal::maxHCalClus], cx[hcal::maxHCalClus], cy[hcal::maxHCalClus], catime[hcal::maxHCalClus];
    double HODOtmean;

    // Speed up processing by switching on only useful branches
    C->SetBranchStatus( "*", 0 );
    C->SetBranchStatus( "sbs.hcal.x", 1 );
    C->SetBranchStatus( "sbs.hcal.y", 1 );
    C->SetBranchStatus( "sbs.hcal.e", 1 );
    C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
    C->SetBranchStatus( "sbs.hcal.colblk", 1 );
    C->SetBranchStatus( "sbs.hcal.nblk", 1 );
    C->SetBranchStatus( "sbs.hcal.nclus", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    //C->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.row", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.col", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.x", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.y", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
    C->SetBranchStatus( "bb.tr.n", 1 );
    C->SetBranchStatus( "bb.tr.px", 1 );
    C->SetBranchStatus( "bb.tr.py", 1 );
    C->SetBranchStatus( "bb.tr.pz", 1 );    
    C->SetBranchStatus( "bb.tr.vz", 1 );
    C->SetBranchStatus( "bb.tr.p", 1 );
    C->SetBranchStatus( "bb.ps.e", 1 );
    C->SetBranchStatus( "bb.ps.x", 1 );
    C->SetBranchStatus( "bb.ps.y", 1 );
    C->SetBranchStatus( "bb.sh.e", 1 );
    C->SetBranchStatus( "bb.sh.x", 1 );
    C->SetBranchStatus( "bb.sh.y", 1 );
    C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    C->SetBranchStatus( "bb.gem.track.nhits", 1 );
    C->SetBranchStatus( "bb.etot_over_p", 1 );
    C->SetBranchStatus( "Ndata.sbs.hcal.clus.id", 1 ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus
    C->SetBranchStatus( "Ndata.sbs.hcal.clus_blk.id", 1 );

    // Linking memory
    C->SetBranchAddress( "sbs.hcal.x", &HCALx );
    C->SetBranchAddress( "sbs.hcal.y", &HCALy );
    C->SetBranchAddress( "sbs.hcal.e", &HCALe );
    C->SetBranchAddress( "sbs.hcal.rowblk", &pblkrow );
    C->SetBranchAddress( "sbs.hcal.colblk", &pblkcol );
    C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
    C->SetBranchAddress( "sbs.hcal.nclus", &nclus ); // Total number of clusters
    C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
    C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
    C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime ); // Array of block TDC times
    C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times
    //C->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain ); // Array of block ADC gain coeff
    C->SetBranchAddress( "sbs.hcal.clus.id", cid );
    C->SetBranchAddress( "sbs.hcal.clus.row", crow );
    C->SetBranchAddress( "sbs.hcal.clus.col", ccol );
    C->SetBranchAddress( "sbs.hcal.clus.e", ce );
    C->SetBranchAddress( "sbs.hcal.clus.x", cx );
    C->SetBranchAddress( "sbs.hcal.clus.y", cy );
    C->SetBranchAddress( "sbs.hcal.clus.atime", catime );
    C->SetBranchAddress( "bb.tr.n", &BBtr_n );
    C->SetBranchAddress( "bb.tr.px", BBtr_px );
    C->SetBranchAddress( "bb.tr.py", BBtr_py );
    C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
    C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
    C->SetBranchAddress( "bb.tr.p", BBtr_p );
    C->SetBranchAddress( "bb.ps.e", &BBps_e );
    C->SetBranchAddress( "bb.ps.x", &BBps_x );
    C->SetBranchAddress( "bb.ps.y", &BBps_y );
    C->SetBranchAddress( "bb.sh.e", &BBsh_e );
    C->SetBranchAddress( "bb.sh.x", &BBsh_x );
    C->SetBranchAddress( "bb.sh.y", &BBsh_y ); 
    C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Ncid ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus
    C->SetBranchAddress( "Ndata.sbs.hcal.clus_blk.id", &Ncblkid );
    
    //globalcut enables
    C->SetBranchStatus( "bb.tr.tg_th", 1 );
    C->SetBranchStatus( "bb.tr.tg_ph", 1 );
    
    //Use TTreeFormula to avoid looping over data an additional time
    TCut GCut = cut[0].gcut.c_str();
    
    //Add globalcut and elastic cuts for reporting
    if( target_change_index != target_index ){
      std::cout << target_parameters;
      target_change_index = target_index;
    }
    
    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );
    
    // Set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta_rad, hcalaxes );
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal::HCalvoff*hcalaxes[0];
    double BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    double Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;
    
    long nevent = 0, nevents = C->GetEntries(); 
    int treenum = 0, currenttreenum = 0;
    
    cout << "Analyzing run " << current_runnumber << " with " << nevents << ". This is good run number " << good_run+1 << endl;
    good_run++;

    //Main loop over events in run
    while (C->GetEntry(nevent++)) {
      
      cout << "Analyzing " << current_target << " run " << current_runnumber << ": " <<  nevent << "/" << nevents << "  \r";
      cout.flush();
      
      ///////
      //Single-loop elastic globalcut method. Save pass/fail for output tree.
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;
      
      if( failedglobal ) 
	continue;
      
      ///////
      //Physics calculations
      //correct beam energy with vertex information, primary track
      double ebeam_c = ebeam - ( (BBtr_vz[0]+target_length/2.0) * target_rho * target_dEdx + upstream_wthick * cell_rho * cell_dEdx );
      
      TVector3 vertex( 0., 0., BBtr_vz[0] );
      
      //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      double precon = BBtr_p[0] + Eloss_outgoing;
      
      //set up four-momenta with some empty for various calculation methods
      TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
      TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
      TLorentzVector ptarg; //target momentum
      ptarg.SetPxPyPzE( 0., 0., 0., M_avg );
      TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
      TVector3 qv = q.Vect();
      TLorentzVector pN; //N' momentum
      
      //simple calculations for e' and N'
      double etheta = acos( pe.Pz() / pe.E() );
      double ephi = atan2( pe.Py(), pe.Px() );
      double pcent = ebeam_c/( 1. + ( ebeam_c/M_avg )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
      double phNexp = ephi + hcal::PI;
      double Q2, W2;
      
      //e' p reconstruction with track angles (not momentum)
      double nu = pbeam.E() - pcent;
      double pNexp = sqrt( pow(nu, 2.) + 2. * M_avg * nu );
      double thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
      TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
      pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
      Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
      W2 = pow( M_avg, 2.0 ) + 2.0*M_avg * (ebeam_c-pe.E()) - Q2;
      
      //Calculate h-arm quantities
      vector<double> xyhcalexp; util::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
      double dx = HCALx - xyhcalexp[0];
      double dy = HCALy - xyhcalexp[1];
      
      ///////
      //BigBite/SBS Acceptance matching.
      bool failedaccmatch = 
	xyhcalexp[1] > hcal::posHCalYf ||
	xyhcalexp[1] < hcal::posHCalYi ||
	xyhcalexp[0] > hcal::posHCalXf ||
	xyhcalexp[0] < hcal::posHCalXi;
      
      if( failedaccmatch ) continue;
      
      double pblkid = (double)cblkid[0]-1;
      
      double hcalatime = cblkatime[0];
      //double hcalatime_bc = hcalatime-BBsh_atime; //enable if bbcal better resolution than hodo
      double hcalatime_bc = hcalatime-HODOtmean;
      double hcalatime_c = hcalatime + old_time_offsets[(int)pblkid] - new_time_offsets[(int)pblkid];
      //double hcalatime_fc = hcalatime-BBsh_atime; //enable if bbcal better resolution than hodo
      double hcalatime_fc = hcalatime_c-HODOtmean;
      
      //Assuming the tdc calib hasn't changed, which is true on SBS14 and SBS11 from pass0/1 to pass2
      double hcaltime = cblktime[0];
      double hcaltime_bc = hcaltime-HODOtmean;
      
      double time_before;
      double time_after;
      
      time_before = hcalatime;
      time_after = hcalatime_c;
      
      ht_run->Fill(current_runnumber,time_after);
      ht_id->Fill(pblkid,time_after);
      ht_allchan_allruns->Fill(time_after);
      
      //loop over hcal primary cluster blocks to return cluster time differences
      double diff_mark = 1000.;
      double diff_mark_idx = -1;
      for( int c=1; c<Ncblkid; c++ ){
	double id = cblkid[c]-1;
	double atime = cblkatime[c];
	double atime_c = atime + old_time_offsets[(int)id] - new_time_offsets[(int)id];
	
	//time difference between primary cluster and each additional cluster
	double diff = hcalatime - atime;
	double diff_c = hcalatime_c - atime_c;
	ht_all_diff->Fill(diff_c);

	ht_cidx->Fill(c,diff_c);
	if( c==select )
	  ht_all_select_diff->Fill(diff_c);

	//check to see if the magnitude of the time diff is closer to zero than the last block
	if( abs(diff_c) < diff_mark && c!=0 ){
	  diff_mark_idx = c;
	  diff_mark = diff_c;
	}
      }
      
      if( diff_mark_idx!=-1) 
	ht_all_best_diff->Fill(diff_mark);
      ht_all_best_diff_idx->Fill(diff_mark_idx);

    }//end loop over event
    
    // getting ready for the next run
    C->Reset();
        
  }//endloop over runs
  
  gStyle->SetOptStat(0);

  //set up blank TH1D for widest tdiff histogram on canvas 3
  TH1D *hslice;

  TCanvas* c2 = new TCanvas("c2", "Fit overlay", 1600, 800);

  overlayWithGaussianFits(fout,ht_cidx,c2,hslice);

  string plotdir2 = Form("../quality_plots/%s/conf%d/tmax_bdiff_sbs%d_set%d_timevblk.png",experiment.c_str(),config,config,set);
  c2->Write();
  c2->SaveAs(plotdir2.c_str());

  TCanvas* c1 = new TCanvas("c1", "Gaussian Fit", 1600, 800);
  c1->SetGridx(1);

  // Define the Gaussian fit function and constrain its range
  TF1 *fitFunction = new TF1("fitFunction", skewedGaussPlusConstant, xmin, xmax, 5);

  // Set initial parameter guesses
  fitFunction->SetParameter(0, 50000);     // Initial amplitude guess
  fitFunction->SetParameter(1, 0);    // Initial mean guess
  fitFunction->SetParameter(2, 3);    // Initial sigma guess
  fitFunction->SetParameter(3, 1);     // Initial skewness guess
  fitFunction->SetParameter(4, 500);     // Initial y-intercept guess for the linear

  ht_all_diff->Fit(fitFunction, "Q", "", xmin, xmax);  // "Q" option for quiet mode

  auto [maxVal, x_max] = computeMaxValue(fitFunction, fitFunction->GetParameter(1), aFWHM);

  double hwhmLeftDifference = computeHWHMLeft(fitFunction, x_max, maxVal);

  double diffFromZeroLeft = computeHWHMLeftZero(fitFunction, x_max, maxVal);
  double diffFromZeroRight = computeHWHMRightZero(fitFunction, x_max, maxVal);

  // Access the parameters of the fit
  double amplitude = fitFunction->GetParameter(0);
  double mean = fitFunction->GetParameter(1);
  double sigma = fitFunction->GetParameter(2);

  // Draw the histogram (with the fit overlaid)
  ht_all_diff->SetTitle(Form("ADCt PClus PBlk-Blk (all diff) SBS-%d, Set %d",config,set));
  ht_all_diff->GetXaxis()->SetTitle("ns");
  ht_all_diff->Draw();

  // Draw the line representing the HWHM
  double y1 = fitFunction->Eval(x_max - hwhmLeftDifference);  // y-location of the HWHM
  double y2 = fitFunction->Eval(x_max);                        // y-location of the maximum
  
  TLine *line1 = new TLine(x_max, y1, x_max, y2);
  line1->SetLineColor(kBlue);
  line1->Draw();
  
  TLine *line2 = new TLine(x_max - hwhmLeftDifference, y1, x_max, y1);
  line2->SetLineColor(kBlue);
  line2->SetLineStyle(2);  // Dashed line
  line2->Draw();

  double yOffset = maxVal * 0.0125;  // Adjust this value to control the vertical offset

  double xcut = nsig*diffFromZeroLeft;

  // Draw lines which represent difference between zero and HWHM right/left
  TLine *lineLeft = new TLine(-diffFromZeroLeft, y1-yOffset, 0., y1-yOffset);
  lineLeft->SetLineColor(kGreen);
  lineLeft->SetLineStyle(2);  // Dashed line
  lineLeft->SetLineWidth(2);
  lineLeft->Draw();

  TLine *lineRight = new TLine(0., y1-yOffset, diffFromZeroRight, y1-yOffset);
  lineRight->SetLineColor(kMagenta);
  lineRight->SetLineStyle(2);  // Dashed line
  lineRight->SetLineWidth(2);
  lineRight->Draw();

  //draw line to show cut
  TLine *lineCutLeft = new TLine(-xcut, yOffset-2*yOffset, -xcut, yOffset+2*yOffset);
  lineCutLeft->SetLineColor(kGreen);
  //lineCutLeft->SetLineStyle(2);  // Dashed line
  lineCutLeft->SetLineWidth(2);
  lineCutLeft->Draw();

  //draw line to show cut
  TLine *lineCutRight = new TLine(xcut, yOffset-2*yOffset, xcut, yOffset+2*yOffset);
  lineCutRight->SetLineColor(kGreen);
  //lineCutRight->SetLineStyle(2);  // Dashed line
  lineCutRight->SetLineWidth(2);
  lineCutRight->Draw();

  // Create a legend to report the fit parameters
  TLegend* leg = new TLegend(0.55, 0.6, 0.88, 0.88);
  //leg->AddEntry(ht_allchan_allruns, Form("ADCt SBS-%d, Set %d",config,set), "");
  leg->AddEntry(fitFunction, "Skewed Gaussian Fit", "l");
  leg->AddEntry(line1, "Half Max", "l");
  leg->AddEntry(line2, "HWHM Left", "l");
  leg->AddEntry(lineLeft, "HWHM Difference from x=0 (Left)", "l");
  leg->AddEntry(lineRight, "HWHM Difference from x=0 (Right)", "l");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry((TObject*)0, Form("HWHM Left Magnitude (ns): %0.2f", hwhmLeftDifference ), "");
  leg->AddEntry((TObject*)0, Form("HWHM Left from x=0 Magnitude (ns): %0.2f", diffFromZeroLeft ), "");
  leg->AddEntry((TObject*)0, Form("HWHM Right from x=0 Magnitude (ns): %0.2f", diffFromZeroRight ), "");
  leg->AddEntry((TObject*)0, Form("tmax (%d HWHM Left): %0.2f",nsig,nsig*hwhmLeftDifference ), "");
  leg->AddEntry((TObject*)0, Form("tmax (%d HWHM from x=0 Left): %0.2f",nsig,nsig*diffFromZeroLeft ), "");
  leg->AddEntry((TObject*)0, Form("tmax (%d HWHM from x=0 Right): %0.2f",nsig,nsig*diffFromZeroRight ), "");
  leg->Draw();

  // Save the canvas to a file (if needed)
  string plotdir = Form("../quality_plots/%s/conf%d/tmax_bdiff_sbs%d_set%d.png",experiment.c_str(),config,config,set);
  c1->Write();
  c1->SaveAs(plotdir.c_str());

  //To tune pars
  cout << "Skewed gaussian fit parameters:" << endl;
  for( int i=0; i<5; ++i )
    cout << "par" << i << ": " << fitFunction->GetParameter(i) << endl;


  fout->Write();

  cout << endl << "tmax.C complete. Output file written to " << fout_path << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    

}
