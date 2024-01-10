//SSeeds 8.2.23 Script to use gain coefficients from a passed calibration set timestamp with data and database parameters from a given experiment/configuration to investigate effects of these coefficents/offsets on a broad set of data. Will use to evaluate the need for separate calibrations per kinematic

#include <vector>
#include <iostream>
#include <algorithm>
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

const Double_t minEventPerCell = 100; //minimum elastic events detected in cell to be considered for calibration
const Double_t highDelta = 0.01; //max allowed discrepancy factor between energy deposited in block and expected
const Int_t hcal_first_channel = 0; //first channel in hcal index
const Int_t bins_magnet_percent = 21; //SBS field at 2100A, so divide by this factor to get percent
const Double_t start_magnet_percent = 0.;
const Double_t end_magnet_percent = 105.;
const Int_t bins_dxdy = 250;
const Double_t start_dx = -4.;
const Double_t end_dx = 3.;
const Double_t start_dy = -1.25;
const Double_t end_dy = 1.25;
const Double_t bins_E = 400;
const Double_t start_E = 0.;
//const Double_t end_E = 1.0;
//const Double_t end_E = 1.8; //for now comment in for SBS11
const Double_t end_E = 1.4; //for now comment in for SBS7
const Double_t bins_SF = 400;
const Double_t start_SF = 0.;
const Double_t end_SF = 0.4;
const Double_t bins_pC = 400;
const Double_t start_pC = 0.;
const Double_t end_pC = 1.0/0.0017; //estimated range from cosmic gain coefficient
const Int_t linecount = 25;
const Int_t atimeNsig = 6;
const Int_t total_bins = 320;
const Int_t tdc_lower_lim = -140;
const Int_t tdc_upper_lim = 20;
const Int_t lower_lim = -60;
const Int_t upper_lim = 100;
const Double_t E_approx_FWHM = 0.3;
const Double_t pC_approx_FWHM = 0.3/0.0017;
const Int_t samples[2] = {161,80}; //Sample channels to check timing alignment over runs within calibrations set
//const Double_t test_coeff = 0.0017;


void overlayWithGaussianFits(TH2D* hist, TCanvas* canvas, bool pCopt, const std::vector<int>& redBins) {
  if (!hist || !canvas) {
    std::cerr << "Null histogram or canvas pointer!" << std::endl;
    return;
  }

  // Use the provided canvas
  canvas->cd();

  // Create TGraphErrors for standard and red points
  TGraphErrors* graph = new TGraphErrors();
  TGraphErrors* redGraph = new TGraphErrors();

  // //get x value of first bin
  // Int_t xstart = hist->GetBinCenter(1);
  // //

  for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
    // Project the 2D histogram onto a 1D histogram along the y-axis for the current x-bin
    TH1D* projY = hist->ProjectionY("_py", binX, binX);

    //Skip the fit if less than 100 entries exist in the bin
    if( projY->GetEntries()<100 )
      continue;
	
    //declare some dynamic fit variables
    Double_t FWHM = E_approx_FWHM;
    if(pCopt) FWHM = pC_approx_FWHM;
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

      // Determine which graph to fill based on whether the bin is in redBins
      

      TGraphErrors* currentGraph = graph;
      Int_t binXval = hist->GetXaxis()->GetBinCenter(binX);

      //cout << "binXval: " << binXval << endl;

      if (std::find(redBins.begin(), redBins.end(), binXval) != redBins.end()) {
	currentGraph = redGraph;
      }

      int pointIdx = currentGraph->GetN();
      currentGraph->SetPoint(pointIdx, xCenter, mean);
      currentGraph->SetPointError(pointIdx, 0, sigma);
    }
    delete projY; // Clean up
  }

  // Draw the original histogram, standard graph, and red graph
  hist->Draw("COLZ");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("P SAME");

  redGraph->SetMarkerStyle(20);
  redGraph->SetMarkerColor(kRed);
  redGraph->SetLineColor(kRed);
  redGraph->Draw("P SAME");

  // Update the canvas
  canvas->Update();
}

//Main <experiment> <configuration> <replay-pass> <target-option> <selected-energy-timestamp> <replacement-experiment> <replacement-configuration> <replacement-pass> <replacement-database-timestamp> <verbose-option>;
//Timestamp GMn SBS4 set 2: -------[ 2021-10-24 04:30:00 ]
//For the first set from a given kinematic, pass "none" to avoid any runs from additional sets

void qreplay_standalone( const char *experiment = "gmn", 
			 Int_t config = 8, 
			 Int_t pass = 1,  
			 const char *sts = "", //if entire set should be analyzed with new gain parameters, pass the argument ""
			 const char *rexperiment = "gmn", 
			 Int_t rconfig = 8, 
			 Int_t rpass = 1, 
			 const char *qts = "",
			 bool h2only = false,
			 bool verbose = true ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  // Get the date
  string date = util::getDate();

  // outfile path
  std::string h2opt = "";
  if( h2only )
    h2opt = "_lh2only";
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string qr_path = outdir_path + Form("/hcal_calibrations/qreplay/qreplay%s_from_%s_%d_%d_to_%s_%d_%d.root",h2opt.c_str(),rexperiment,rconfig,rpass,experiment,config,pass);

  //Set up path variables and output files
  TFile *fout = new TFile( qr_path.c_str(), "RECREATE" );

  std::string qr_report_path = Form("report_path/qreplay%s_from_%s_%d_%d_to_%s_%d_%d.txt",h2opt.c_str(),experiment,config,pass,rexperiment,rconfig,rpass);
  std::string db_path = gSystem->Getenv("DB_DIR");
  std::string hcal_db_path = db_path + "/db_sbs.hcal.dat";
  //Using trial calibrated set
  //std::string r_adcgain_path = Form("../energy/parameters/adcgaincoeff_%s%s_conf%d_pass%d.txt",rexperiment,h2opt.c_str(),rconfig,rpass);
  std::string r_adcgain_path = Form("../energy/parameters/adcgaincoeff_%s_lh2only_conf%d_pass%d.txt",rexperiment,rconfig,rpass); //check only calibration constants from lh2 for comparisons
  //Use base calibrated set, not trial replacement calibrated set
  std::string r_adct_path = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string r_tdcoffset_path = Form("../timing/parameters/tdcoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string db_gain_variable = "sbs.hcal.adc.gain";
  std::string db_adcg_variable = "sbs.hcal.adc.gain";
  std::string db_adct_variable = "sbs.hcal.adc.timeoffset";
  std::string db_tdcoffset_variable = "sbs.hcal.tdc.offset";
  std::string db_tdccalib_variable = "sbs.hcal.tdc.calib";
  
  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available

  //define structures for holding cut report data and indices	    
  reportset report_set[hcal::gNmag];
  Int_t report_set_size = 0;
  Int_t report_set_idx = 0;
  
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);   
  std::cout << "Configuration parameters for data to qreplay:" << std::endl;
  std::cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  //Add quality plots   
  TH1D * hE_old_nocut = new TH1D("hE_old_nocut",
				 "HCal Primary Cluster E, Old Coeff, No Elastic Cuts; GeV",
				 bins_E,
				 start_E,
				 end_E);

  TH1D * hE_old = new TH1D("hE_old",
			   "HCal Primary Cluster E, Old Coeff, Elastic Cuts; GeV",
			   bins_E,
			   start_E,
			   end_E);

  TH1D * hpC_nocut = new TH1D("hpC_old_nocut",
				 "HCal Primary Cluster pC, No Coeff, No Elastic Cuts; GeV",
				 bins_pC,
				 start_pC,
				 end_pC);

  TH1D * hpC = new TH1D("hpC_old",
			   "HCal Primary Cluster pC, No Coeff, Elastic Cuts; GeV",
			   bins_pC,
			   start_pC,
			   end_pC);

  TH1D * hE_new_nocut = new TH1D("hE_new_nocut",
				 "HCal Primary Cluster E, New Coeff, No Elastic Cuts; GeV",
				 bins_E,
				 start_E,
				 end_E);

  TH1D * hE_new = new TH1D("hE_new",
			   "HCal Primary Cluster E, New Coeff, Elastic Cuts; GeV",
			   bins_E,
			   start_E,
			   end_E);
    
  TH1D * hSF_old_nocut = new TH1D("hSF_old_nocut",
				 "HCal Primary Cluster SF, Old Coeff, No Elastic Cuts; GeV",
				 bins_SF,
				 start_SF,
				 end_SF);

  TH1D * hSF_old = new TH1D("hSF_old",
			   "HCal Primary Cluster SF, Old Coeff, Elastic Cuts; GeV",
			   bins_SF,
			   start_SF,
			   end_SF);

  TH1D * hSF_new_nocut = new TH1D("hSF_new_nocut",
				 "HCal Primary Cluster SF, New Coeff, No Elastic Cuts; GeV",
				 bins_SF,
				 start_SF,
				 end_SF);

  TH1D * hSF_new = new TH1D("hSF_new",
			   "HCal Primary Cluster SF, New Coeff, Elastic Cuts; GeV",
			   bins_SF,
			   start_SF,
			   end_SF);

  TH2D *hdxvmag_h = new TH2D("hdxvmag_h",
			     "dx vs magnetic field lh2",
			     bins_magnet_percent,
			     start_magnet_percent,
			     end_magnet_percent,
			     bins_dxdy,
			     start_dx,
			     end_dx);
    
  TH2D *hdyvmag_h = new TH2D("hdyvmag_h",
			     "dy vs magnetic field lh2",
			     bins_magnet_percent,
			     start_magnet_percent,
			     end_magnet_percent,
			     bins_dxdy,
			     start_dy,
			     end_dy);
    
  TH2D *hdxvmag_d = new TH2D("hdxvmag_d",
			     "dx vs magnetic field LD2",
			     bins_magnet_percent,
			     start_magnet_percent,
			     end_magnet_percent,
			     bins_dxdy,
			     start_dx,
			     end_dx);
    
  TH2D *hdyvmag_d = new TH2D("hdyvmag_d",
			     "dy vs magnetic field LD2",
			     bins_magnet_percent,
			     start_magnet_percent,
			     end_magnet_percent,
			     bins_dxdy,
			     start_dy,
			     end_dy);
    
  TH2D *hdxvmag_h3 = new TH2D("hdxvmag_h3",
			      "dx vs magnetic field He3",
			      bins_magnet_percent,
			      start_magnet_percent,
			      end_magnet_percent,
			      bins_dxdy,
			      start_dx,
			      end_dx);
    
  TH2D *hdyvmag_h3 = new TH2D("hdyvmag_h3",
			      "dy vs magnetic field He3",
			      bins_magnet_percent,
			      start_magnet_percent,
			      end_magnet_percent,
			      bins_dxdy,
			      start_dy,
			      end_dy);
   
  TH2D *hEvX = new TH2D("hEvX",
			"HCal energy vs HCal X",
			hcal::maxHCalRows,
			hcal::posHCalXi,
			hcal::posHCalXf,
			bins_E,
			start_E,
			end_E);
   
  TH2D *hEvY = new TH2D("hEvY",
			"HCal energy vs HCal Y",
			hcal::maxHCalRows,
			hcal::posHCalYi,
			hcal::posHCalYf,
			bins_E,
			start_E,
			end_E);

  TH2D *hSFvID = new TH2D("hSFvID",
			  "sampling fraction vs ID",
			  hcal::maxHCalChan,
			  hcal_first_channel,
			  hcal::maxHCalChan,
			  bins_SF,
			  start_SF,
			  end_SF);

  TH2D *hSFvrow = new TH2D("hSFvrow",
			   "sampling fraction vs row",
			   hcal::maxHCalRows,
			   hcal_first_channel,
			   hcal::maxHCalRows,
			   bins_SF,
			   start_SF,
			   end_SF);

  TH2D *hSFvcol = new TH2D("hSFvcol",
			   "sampling fraction vs col",
			   hcal::maxHCalCols,
			   hcal_first_channel,
			   hcal::maxHCalCols,
			   bins_SF,
			   start_SF,
			   end_SF);

  TH2D *hap_hodocorr_ID = new TH2D("hap_hodocorr_ID",
				   "adct hodoscope corrected vs ID",
				   hcal::maxHCalChan,
				   hcal_first_channel,
				   hcal::maxHCalChan,
				   total_bins,
				   lower_lim,
				   upper_lim);

  //Assumes that the runlist is ordered in ascending run numbers
  TH2D *hadct_samp1_run = new TH2D(Form("hadct_s%d_run",samples[0]),
				   Form("adct vs run channel %d",samples[0]),
				   (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				   runs[0].runnum-1,
				   runs[nruns-1].runnum+1,
				   total_bins,
				   lower_lim,
				   upper_lim);

  TH2D *hadct_samp2_run = new TH2D(Form("hadct_s%d_run",samples[1]),
				   Form("adct vs run channel %d",samples[1]),
				   (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				   runs[0].runnum-1,
				   runs[nruns-1].runnum+1,
				   total_bins,
				   lower_lim,
				   upper_lim);


  TH2D *hE_run_nocut = new TH2D("hE_run_nocut",
				"Corrected E vs run all channels no elastic cut",
				(runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				runs[0].runnum-1,
				runs[nruns-1].runnum+1,
				bins_E,
				start_E,
				end_E);


  TH2D *hE_old_run = new TH2D("hE_old_run",
			  "E (old coeff) vs run all channels",
			  (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
			  runs[0].runnum-1,
			  runs[nruns-1].runnum+1,
			  bins_E,
			  start_E,
			  end_E);



  TH2D *hE_run = new TH2D("hE_run",
			  "Corrected E vs run all channels",
			  (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
			  runs[0].runnum-1,
			  runs[nruns-1].runnum+1,
			  bins_E,
			  start_E,
			  end_E);


  TH2D *hpC_run = new TH2D("hpC_run",
			  "pC vs run all channels",
			  (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
			  runs[0].runnum-1,
			  runs[nruns-1].runnum+1,
			  bins_pC,
			  start_pC,
			  end_pC);

    
  TH1D *hadctblk = new TH1D("hadctblk",
			    "adct all channels all blocks primary cluster",
			    total_bins,
			    lower_lim,
			    upper_lim);

    
  TH1D *hadct_nocut = new TH1D("hadct_nocut",
			 "adct all channels primary block primary cluster no elastic cut",
			 total_bins,
			 lower_lim,
			 upper_lim);

  TH1D *hadct = new TH1D("hadct",
			 "adct all channels primary block primary cluster",
			 total_bins,
			 lower_lim,
			 upper_lim);

  TH2D *htp_hodocorr_ID = new TH2D("htp_hodocorr_ID",
				   "tdc hodoscope corrected vs ID",
				   hcal::maxHCalChan,
				   hcal_first_channel,
				   hcal::maxHCalChan,
				   total_bins,
				   tdc_lower_lim,
				   tdc_upper_lim);

  
  //cut quality plots
  TH1D * hfglobal = new TH1D("hfglobal",
			     "failed global cut, 1=fail",
			     3,
			     -1,
			     2);

  TH1D * hfactivearea = new TH1D("hfactivearea",
				 "failed active area cut, 1=fail",
				 3,
				 -1,
				 2);
  TH1D * hfcoin = new TH1D("hfcoin",
			   "failed adct coin cut, 1=fail",
			   3,
			   -1,
			   2);
  TH1D * hfdy = new TH1D("hfdy",
			 "failed dy cut, 1=fail",
			 3,
			 -1,
			 2);
  TH1D * hpid = new TH1D("hpid",
			 "Particle",
			 6,
			 -2,
			 4);
  TH1D * hfW2 = new TH1D("hfW2",
			 "failed W2 cut, 1=fail",
			 3,
			 -1,
			 2);


  //reporting indices
  Int_t target_change_index;
  Double_t config_sampling_fraction;
  Double_t config_e_sigma_ratio;
  std::string ts_compare = "";
  bool first = true;
  bool cut_first = true;
  std::vector<Int_t> lh2runs;
  
  //TEST
  vector<int> elastics_per_run;
  vector<int> elastic_runs;
  int total_elastics_allruns = 0;
  
  //Main loop over runs
  for (Int_t r=0; r<nruns; r++) {
    bool ts_different = false;

    //Get run experimental parameters
    std::string current_timestamp = runs[r].adcg_ts; //For now, only select on energy calibration sets
    Int_t current_runnumber = runs[r].runnum;
    std::string current_target = runs[r].target;

    //Will use these to differentiate the tgraph later on
    if( current_target.compare("lh2")==0 )
      lh2runs.push_back(current_runnumber);

    if(verbose) 
      std::cout << "Current run number: " << current_runnumber << ", current timestamp: " << current_timestamp << ", current target: " << current_target << std::endl;

    if( !(current_timestamp.compare(ts_compare)==0) ){
      ts_compare = current_timestamp;
      ts_different = true;
    }

    //Skip any runs outside of select timestamp calibration set, where it is passed
    std::string select_timestamp = sts;
    std::string quality_timestamp = qts;
    if( !(current_timestamp.compare(select_timestamp)==0) && !(select_timestamp.compare("")==0) )
      continue;

    std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    Int_t mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A

    //Get run paths
    std::string rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-%s/pass%d/SBS%d/%s/rootfiles/",experiment,pass,config,targ_uppercase.c_str());
    std::string rootfile_path = rootfile_dir + Form("*%d*",current_runnumber);

    //Get target configuration
    SBStarget target_parameters(current_target);
    Int_t target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
    Double_t target_length = target_parameters.GetTargLength();
    Double_t target_rho = target_parameters.GetTargRho();
    Double_t cell_rho = target_parameters.GetCellRho();
    Double_t cell_diam = target_parameters.GetCellDiam();
    Double_t cell_dEdx = target_parameters.GetCelldEdx();
    Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
    Double_t target_dEdx = target_parameters.GetTargdEdx();
    Double_t M_avg = target_parameters.GetAvgMass();
        
    //check target and continue on deuterium of only calibrating with hydrogen
    if( h2only && target_index!=1 )
      continue;

    //Record old gain and offset parameters by tstamp from database (assumes one file for all timestamps)
    Double_t old_adcg_coeff[hcal::maxHCalChan] = {0.};
    Double_t old_adct_offsets[hcal::maxHCalChan] = {0.};
    Double_t old_tdc_offsets[hcal::maxHCalChan] = {0.};
    Double_t old_tdc_calib = 1;

    if(verbose) std::cout << "Reading parameters from file: " << hcal_db_path << std::endl;
    util::readDB( hcal_db_path, runs[r].adcg_ts, db_adcg_variable, old_adcg_coeff );
    util::readDB( hcal_db_path, runs[r].adct_ts, db_adct_variable, old_adct_offsets );
    util::readDB( hcal_db_path, runs[r].tdc_ts, db_tdcoffset_variable, old_tdc_offsets );
    util::readDB( hcal_db_path, runs[r].tdcc_ts, db_tdccalib_variable, old_tdc_calib );

    //get new adct offsets from calibrated parameters
    Double_t r_adcg_coeff[hcal::maxHCalChan] = {0.};
    Double_t r_adct_offsets[hcal::maxHCalChan] = {0.};
    Double_t r_tdc_offsets[hcal::maxHCalChan] = {0.};

    //if no selected timestamp is passed, we use the first set of new coefficients that appears
    if( quality_timestamp.compare("")==0 ) quality_timestamp = "none";
    if(verbose) std::cout << "Reading parameters from file: " << r_adcgain_path << std::endl;
    util::readDB( r_adcgain_path, quality_timestamp, db_adcg_variable, r_adcg_coeff );
    if(verbose) std::cout << "Reading parameters from file: " << r_adct_path << std::endl;
    util::readDB( r_adct_path, "none", db_adct_variable, r_adct_offsets );
    if(verbose) std::cout << "Reading parameters from file: " << r_tdcoffset_path << std::endl;
    util::readDB( r_tdcoffset_path, "none", db_tdcoffset_variable, r_tdc_offsets );
    
    //Report offset parameters on first loop
    if( verbose && first ){
      std::cout << "Current set of old ADCt offset coefficients:" << std::endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
	for ( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  std::cout << old_adct_offsets[i] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << "Current set of replacement ADCt offset coefficients:" << std::endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
	for ( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  std::cout << r_adct_offsets[i] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << "Current set of old TDC offset coefficients:" << std::endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
	for ( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  std::cout << old_tdc_offsets[i] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << "Current set of replacement TDC offset coefficients:" << std::endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
	for ( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  std::cout << r_tdc_offsets[i] << " ";
	}
	std::cout << std::endl;
      }
      first = false;
    }

    //Report gain parameters on each database change as necessary
    if( verbose && ts_different ){
      std::cout << std::endl << std::endl << "Current set of old ADC gain coefficients:" << std::endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
	for ( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  std::cout << old_adcg_coeff[i] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << "Current set of replacement ADC gain coefficients:" << std::endl;
      for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
	for ( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  std::cout << r_adcg_coeff[i] << " ";
	}
	std::cout << std::endl;
      }
    }

    //Get available cuts for current config/target/field combination. Set index = 0 as no calibration is being done here and physics cuts do not vary by calibration set
    vector<calcut> cut;
    util::ReadCutList(struct_dir,experiment,config,0,pass,current_target,mag,verb,cut);
    if( cut_first ){
      std::cout << cut[0];
      cut_first = false;
    }

    //TEST
    elastic_runs.push_back(current_runnumber);

    // Setting up chain and branch addresses
    C = new TChain("T");
    C->Add(rootfile_path.c_str());

    C->SetBranchStatus("*",0);    
    Double_t BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
    Double_t BBtr_vz[hcal::maxTracks];
    Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;	
    Double_t HCALx, HCALy, HCALe;
    Double_t pblkrow, pblkcol, nblk, nclus;
    Int_t Ncid;
    Double_t cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan], cblkagain[hcal::maxHCalChan];
    Double_t cid[hcal::maxHCalClus], crow[hcal::maxHCalClus], ccol[hcal::maxHCalClus], ce[hcal::maxHCalClus], cx[hcal::maxHCalClus], cy[hcal::maxHCalClus], catime[hcal::maxHCalClus];
    Double_t HODOtmean;

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
    Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;

    long nevent = 0, nevents = C->GetEntries(); 
    Int_t treenum = 0, currenttreenum = 0;
    
    //TEST
    int total_elastics = 0;

    //Main loop over events in run
    while (C->GetEntry(nevent++)) {

      cout << "Analyzing " << current_target << " run " << current_runnumber << ": " <<  nevent << "/" << nevents << " \r";
      cout.flush();


      ///////
      //Single-loop elastic globalcut method. Save pass/fail for output tree.
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;
	  
      // if( failedglobal ) 
      // 	continue;

      ///////
      //HCal Active Area Cut
      bool failedactivearea = 
	pblkrow==0 || 
	pblkrow==23 || 
	pblkcol==0 || 
	pblkcol==11;

      // if( failedactivearea ) 
      // 	continue; //All events with primary cluster element on edge blocks cut

      ///////
      //HCal primary cluster coincidence time cut (using adctime while hcal tdc suspect, new offsets)
      Int_t pblkid = cblkid[0]-1; //define primary block, primary cluster ID

      Double_t natime = cblkatime[0]+old_adct_offsets[pblkid]-r_adct_offsets[pblkid]; //new atime
      Double_t atime0 = cut[0].atime0; //observed elastic peak in adc time
      Double_t atimesig = cut[0].atime_sig; //observed width of elastic peak in adc time
      Double_t natime_hodo = natime - HODOtmean;

      Double_t ntime = cblktime[0]+old_tdc_offsets[pblkid]-r_tdc_offsets[pblkid];
      Double_t ntime_hodo = ntime - HODOtmean;

      bool failedcoin = abs(natime-atime0)>atimeNsig*atimesig;

      // if( failedcoin ) 
      // 	continue; //All events where adctime outside of reasonable window cut

      ///////
      //Physics calculations
      //correct beam energy with vertex information, primary track
      Double_t ebeam_c = ebeam - ( (BBtr_vz[0]+target_length/2.0) * target_rho * target_dEdx + upstream_wthick * cell_rho * cell_dEdx );

      TVector3 vertex( 0., 0., BBtr_vz[0] );

      //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      Double_t precon = BBtr_p[0] + Eloss_outgoing;

      //set up four-momenta with some empty for various calculation methods
      TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
      TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
      TLorentzVector ptarg; //target momentum
      ptarg.SetPxPyPzE( 0., 0., 0., M_avg );
      TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
      TVector3 qv = q.Vect();
      TLorentzVector pN; //N' momentum
      
      //simple calculations for e' and N'
      Double_t etheta = acos( pe.Pz() / pe.E() );
      Double_t ephi = atan2( pe.Py(), pe.Px() );
      Double_t pcent = ebeam_c/( 1. + ( ebeam_c/M_avg )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
      Double_t phNexp = ephi + hcal::PI;
      Double_t Q2, W2;

      //e' p reconstruction with track angles (not momentum)
      Double_t nu = pbeam.E() - pcent;
      Double_t pNexp = sqrt( pow(nu, 2.) + 2. * M_avg * nu );
      Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
      TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
      pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
      Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
      W2 = pow( M_avg, 2.0 ) + 2.0*M_avg * (ebeam_c-pe.E()) - Q2;

      //Calculate h-arm quantities
      vector<Double_t> xyhcalexp; util::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
      Double_t KE_p = nu; //For elastics total predicted by earm
      Double_t SFrac = HCALe/KE_p; //Measured
      Double_t dx = HCALx - xyhcalexp[0];
      Double_t dy = HCALy - xyhcalexp[1];
      TVector3 neutdir = (hcalpos - vertex).Unit();
      Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + hcal::sbsdipolegap/2.0) );
      TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
      Double_t thetapq_p = acos( protdir.Dot( qv.Unit() ) );
      Double_t thetapq_n = acos( neutdir.Dot( qv.Unit() ) );
	  
      if( current_target.compare("lh2")==0 ){
	hdxvmag_h->Fill( mag, dx );
	hdyvmag_h->Fill( mag, dy );
      }
      if( current_target.compare("ld2")==0 ){
	hdxvmag_d->Fill( mag, dx );
	hdyvmag_d->Fill( mag, dy );
      }
      if( current_target.compare("he3")==0 ){
	hdxvmag_h3->Fill( mag, dx );
	hdyvmag_h3->Fill( mag, dy );
      }

      ///////
      //dy elastic cut
      Double_t dy0 = cut[0].dy0;
      Double_t dysig = cut[0].dy_sig;
      bool faileddy = abs(dy-dy0)>3*dysig;

      // if( faileddy ) 
      // 	continue;

      ///////
      //PID
      Int_t pid = -1;
      util::checkPID(current_target,cut[0].dx0_p,cut[0].dx0_n,cut[0].dx_sig_p,cut[0].dx_sig_n,dx,pid);

      ///////
      //dx elastic cut
      // if( pid==-1 ) 
      // 	continue;

      ////////////////////////////
      //Primary W2 cut on elastics
      Double_t W2_mean = cut[0].W2_mean;
      Double_t W2_sig = cut[0].W2_sig;
      bool failedW2 = fabs(W2-W2_mean)>W2_sig;

      // if( failedW2 ) 
      // 	continue; //Observed mean W2 cut on elastic peak
      
      ///////
      //Get corrected primary cluster energy
      Double_t clusE = 0.0;
      Double_t cluspC = 0.0;
      Double_t clusE_corr = 0.0;
      bool badclus = false;

      //Loop over blocks in primary cluster
      for( Int_t blk = 0; blk<nblk; blk++ ){
	Int_t blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	Double_t blke = cblke[blk];
	Double_t blke_corr = cblke[blk]/old_adcg_coeff[blkid]*r_adcg_coeff[blkid];
	//Double_t blke_corr = cblke[blk]/old_adcg_coeff[blkid]*test_coeff;
	Double_t blkpC = cblke[blk]/old_adcg_coeff[blkid];

	Double_t blknatime = cblkatime[blk]+old_adct_offsets[blkid]-r_adct_offsets[blkid]; //new atime

	////////
	//Cluster block coincidence adc time check, 5sig estimate
	if( abs(blknatime-atime0)>atimeNsig*atimesig ){
	  badclus = 1;
	  //continue;
	}
	
	clusE += blke;
	cluspC += blkpC;
	clusE_corr += blke_corr;
   
      }

      //Assume that no cut placed on adc time in the cluster loop
      if( clusE != HCALe )
	cout << "ERROR: sbs.hcal.e ne primary cluster sum E" << endl;

      Double_t sampling_fraction = clusE/KE_p;
      Double_t sampling_fraction_corr = clusE_corr/KE_p;

      hE_old_nocut->Fill( clusE );
      hE_new_nocut->Fill( clusE_corr );
      hpC_nocut->Fill( cluspC );
      hSF_old_nocut->Fill( sampling_fraction );
      hSF_new_nocut->Fill( sampling_fraction_corr );
      hfglobal->Fill( failedglobal );
      hfactivearea->Fill( failedactivearea );
      hfcoin->Fill( failedcoin );
      hfdy->Fill( faileddy );
      hpid->Fill( pid );
      hfW2->Fill( failedW2 );
      hE_run_nocut->Fill( current_runnumber, clusE_corr);
      hadct_nocut->Fill( natime );

      if( failedglobal || failedactivearea || failedcoin || faileddy || pid==-1 || failedW2 )
	continue;

      //TEST
      total_elastics++;
      total_elastics_allruns++;

      if(pblkid==samples[0])
	hadct_samp1_run->Fill(current_runnumber,natime);
      if(pblkid==samples[1])
	hadct_samp2_run->Fill(current_runnumber,natime);

      hE_old->Fill( clusE );
      hE_new->Fill( clusE_corr );
      hpC->Fill( cluspC );
      hSF_old->Fill( sampling_fraction );
      hSF_new->Fill( sampling_fraction_corr );
      hEvX->Fill( HCALx, clusE );
      hEvY->Fill( HCALy, clusE );
      hSFvID->Fill( cblkid[0], sampling_fraction );
      hSFvrow->Fill( crow[0], sampling_fraction );
      hSFvcol->Fill( ccol[0], sampling_fraction );
      hE_run->Fill( current_runnumber, clusE_corr);
      hE_old_run->Fill( current_runnumber, clusE);
      hpC_run->Fill( current_runnumber, cluspC);
      hap_hodocorr_ID->Fill( pblkid, natime_hodo );
      htp_hodocorr_ID->Fill( pblkid, ntime_hodo );
      hadct->Fill( natime );
     

    }//endloop over events

    //TEST
    elastics_per_run.push_back(total_elastics);

    // getting ready for the next run
    C->Reset();

  }//endloop over runs

  cout << "lh2 runs are as follows: ";
  for( int i=0; i<lh2runs.size(); i++ )
    cout << lh2runs[i] << " ";
  cout << endl;

  //overlay energy vs run histo
  TCanvas *c1 = new TCanvas("c1","E vs Run",1200,1200);
  //gStyle->SetOptStat(0);
  c1->cd();
  overlayWithGaussianFits(hE_run,c1,false,lh2runs);
  c1->Write();

  //overlay pC vs run histo
  TCanvas *c2 = new TCanvas("c2","pC vs Run",1200,1200);
  //gStyle->SetOptStat(0);
  c2->cd();
  overlayWithGaussianFits(hpC_run,c2,true,lh2runs);
  c2->Write();

  //overlay old coeff energy vs run histo
  TCanvas *c3 = new TCanvas("c3","E (old coeff) vs Run",1200,1200);
  //gStyle->SetOptStat(0);
  c3->cd();
  overlayWithGaussianFits(hE_old_run,c3,false,lh2runs);
  c3->Write();

  cout << "Ended loop over runs. Report printed to working directory." << endl;

  if( verbose ){
    std::cout << "LH2 only option:                     " << h2only << std::endl;
    std::cout << "Test experiment parameter:           " << experiment << std::endl;
    std::cout << "Test configuration parameter:        " << config << std::endl;
    std::cout << "Test pass parameter:                 " << pass << std::endl;
    std::cout << "Selected test data timestamp:        " << sts << std::endl;
    std::cout << "Replacement experiment parameter:    " << rexperiment << std::endl;
    std::cout << "Replacement configuration parameter: " << rconfig << std::endl;
    std::cout << "Replacement pass parameter:          " << rpass << std::endl;
    std::cout << "Selected replacement data timestamp: " << qts << std::endl;
  }

  //write report file
  ofstream qreplay_report;
  qreplay_report.open( qr_report_path );
  qreplay_report << "Report file for qreplay path/file " << qr_path << std::endl;
  qreplay_report << "Quality replay from qreplay_standalone.C performed " << date.c_str() << std::endl << std::endl;
  qreplay_report << "LH2 only option:                     " << h2only << std::endl;
  qreplay_report << "Test experiment parameter:           " << experiment << std::endl;
  qreplay_report << "Test configuration parameter:        " << config << std::endl;
  qreplay_report << "Test pass parameter:                 " << pass << std::endl;
  qreplay_report << "Selected test data timestamp:        " << sts << std::endl;
  qreplay_report << "Replacement experiment parameter:    " << rexperiment << std::endl;
  qreplay_report << "Replacement configuration parameter: " << rconfig << std::endl;
  qreplay_report << "Replacement pass parameter:          " << rpass << std::endl;
  qreplay_report << "Selected replacement data timestamp: " << qts << std::endl << std::endl;
  qreplay_report << std::endl << std::endl;
  qreplay_report.close();

  fout->Write();

  st->Stop();

  //TEST
  if( elastics_per_run.size() != elastic_runs.size() ){
    cout << "TEST ERROR: vector size mismatch" << endl;
    return;
  }

  //TEST
  for( int i=0; i<elastic_runs.size(); i++ ){
    cout << "run:" << elastic_runs[i] << "  events: " << elastics_per_run[i] << endl;

  }
  cout << "Total elastics: " << total_elastics_allruns << endl;

  cout << endl << "Output file written to " << qr_path << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    

}
