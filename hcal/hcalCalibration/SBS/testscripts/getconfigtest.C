#include "../include/sbs.h"
#include <vector>
#include <iostream>
#include <algorithm>

const Int_t hcal_first_chan = 0; //first channel in hcal index
const Int_t total_tdc_bins = 320;
const Int_t lower_tdc_lim = -120;
const Int_t upper_tdc_lim = 40;
const Int_t total_adct_bins = 320;
const Int_t lower_adct_lim = -60;
const Int_t upper_adct_lim = 100;
const Double_t bins_SFE = 400;
const Double_t start_SFE = 0.;
const Double_t end_SFE = 1.;

void getconfigtest(const char *experiment = "gmn",Int_t runN=13453,Int_t config = 8,Int_t pass=1){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  // Get the date
  string date = util::getDate();

  // // Get outdir path and declare outfile
  // std::string outdir_path = gSystem->Getenv("OUT_DIR");
  // std::string sbsofflinecheck_path = outdir_path + Form("/hcal_calibrations/qreplay/sbsoffline_check_run%d_pass%d.root",runN,pass);

  // //Set up path variables and output files
  // TFile *fout = new TFile( sbsofflinecheck_path.c_str(), "RECREATE" );

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs;
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
  
  //define structures for holding calibration data and indices
  calset gain_coeff[hcal::gNstamp];
  Int_t Ncal_set_size = 0;
  Int_t Ncal_set_idx = 0;

  //define structures for holding cut report data and indices	    
  reportset report_set[hcal::gNmag];
  Int_t report_set_size = 0;
  Int_t report_set_idx = 0;

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  std::cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  std::cout << "Setting up histograms.." << std::endl;

  TH1D *hE = new TH1D("hE_set%d",
		       "HCal Primary Cluster Energy",
		       bins_SFE,
		       start_SFE,
		       end_SFE);

  TH1D *hSF = new TH1D("hSF_set%d",
		       "HCal Primary Cluster Sampling Fraction",
		       bins_SFE,
		       start_SFE,
		       end_SFE);

  TH2D *htp_hodocorr_ID = new TH2D("htp_hodocorr_ID_set%d",
				   "TDC (hodo corrected) vs ID",
				   hcal::maxHCalChan,
				   hcal_first_chan,
				   hcal::maxHCalChan,
				   total_tdc_bins,
				   lower_tdc_lim,
				   upper_tdc_lim);

  TH2D *hap_hodocorr_ID = new TH2D("hap_hodocorr_ID_set%d",
				   "ADCt (hodo corrected) vs ID",
				   hcal::maxHCalChan,
				   hcal_first_chan,
				   hcal::maxHCalChan,
				   total_adct_bins,
				   lower_adct_lim,
				   upper_adct_lim);
  

  //Find some relevant parameters from the runlist
  Int_t mag;
  std::string current_target;
  std::string current_timestamp;
  
  for (Int_t r=0; r<nruns; r++) {
    
    Int_t current_runnumber = runs[r].runnum;

    if( current_runnumber == runN ){
      mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A
      current_timestamp = runs[r].adcg_ts;
      current_target = runs[r].target;
      break;
    }

    if(r==nruns-1){
      cout << "Did not find run number in run list, aborting.." << endl;
      return;
    }

  }

  //if necessary, convert target string to uppercase
  std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );

  //Get run path - hardcoded to my environment, should contain a test replay of a single run, min 100k events
  std::string rootfile_dir = "/volatile/halla/sbs/seeds/";
  std::string rootfile_path = rootfile_dir + Form("*%d*",runN);

  //Get target configuration
  SBStarget target_parameters(current_target);
  std::cout << target_parameters;
  Int_t target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
  Double_t target_length = target_parameters.GetTargLength();
  Double_t target_rho = target_parameters.GetTargRho();
  Double_t cell_rho = target_parameters.GetCellRho();
  Double_t cell_diam = target_parameters.GetCellDiam();
  Double_t cell_dEdx = target_parameters.GetCelldEdx();
  Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
  Double_t target_dEdx = target_parameters.GetTargdEdx();
  Double_t M_avg = target_parameters.GetAvgMass();

  //Get the relevant index. Note: will need to be updated if more than two calibration sets exist on a single configuration
  Ncal_set_idx = 0;
  std::string active_timestamp;
  util::tsCompare( current_timestamp, config_parameters.GetSBSTimestamp(), active_timestamp );
  if( active_timestamp.compare( current_timestamp )!=0 )
    Ncal_set_idx = 1;

  cout << current_timestamp << " " << config_parameters.GetSBSTimestamp() << endl;


  cout << Ncal_set_idx << endl;

  // //Get available cuts for current config/target/field combination. Use first element (0) of cut
  // vector<calcut> cut;
  // util::ReadCutList(struct_dir,experiment,config,Ncal_set_idx,pass,current_target,mag,verb,cut);
  // std::cout << cut[0];


}
