//SSeeds 8.16.23 Careful reconstruction of tdc_tw.C to exclude by-block parameter extraction then add adct linear timewalk calibrations, tdc linear timewalk calibrations, and traditional tdc timewalk fit calibrations
// Note: last two arguments should not be configured by user on this version of the script unless changes to code to enable their use is performed

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

// Time vs E extraction constraints
const Int_t first_hcal_chan = 0;
const Int_t tdc_bins = 200;
const Int_t adct_bins = 200;
const Double_t tdc_lower_lim = -50.;
const Double_t tdc_upper_lim = 50.;
const Double_t adct_lower_lim = -50.;
const Double_t adct_upper_lim = 50.;
const Int_t E_bins = 1000;
const Double_t E_lower_lim = 0.;
const Double_t E_upper_lim = 0.5;

// expo timewalk fit parameters - may need to be tuned at different hadron momenta
const Double_t tw_fit_llim = 0.02;
const Double_t tw_fit_ulim = 0.50;
const Double_t tw_p0_llim = 4.5;
const Double_t tw_p0_ulim = 9.;
const Double_t tw_p1_llim = 8.;
const Double_t tw_p1_ulim = 12.5;
const Double_t tw_p2_llim = -1.0;
const Double_t tw_p2_ulim = 0.0;

// traditional fit tw fit
const Double_t trad_p0_llim = -2.2;
const Double_t trad_p0_ulim = -1.8;
const Double_t trad_p1_llim = 0.75;
const Double_t trad_p1_ulim = 1.25;
const Double_t trad_p2_llim = 0.49;
const Double_t trad_p2_ulim = 0.51;

//pol1 tdc tw fit parameters
const Double_t p0_set = 7.;
const Double_t p1_set = -18.;
const Double_t p0_ll = 5.;
const Double_t p0_ul = 9.;
const Double_t p1_ll = -20.;
const Double_t p1_ul = -16.;
const Double_t fit_ul = 0.50;
const Double_t fit_ll = 0.01;

//pol1 adct tw fit parameters
const Double_t ap0_set = 0.6;
const Double_t ap1_set = -5.2;
const Double_t ap0_ll = 1.;
const Double_t ap0_ul = 3.;
const Double_t ap1_ll = -7.;
const Double_t ap1_ul = -4.;
const Double_t afit_ul = 0.50;
const Double_t afit_ll = 0.01;

//pol1 adct tw fit parameters

//timewalk tdc set
const Int_t twset = 0;

//trad fit parameters
const Double_t minE = 0.0001; //GeV - added to prevent division by zero on missing hcal clusters

//additional cut constraints
const Int_t fit_event_min = 60;
const Int_t adct_Nsig = 4;

//for testing
const bool verbose = false; //set to true if full list of loaded db params desired sent to console
void coutParams( Double_t params[hcal::maxHCalChan] ){
  Int_t col = 0;
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){
    std::cout << params[i] << " ";
    if(col==11){
      col=0;
      std::cout << std::endl;
      continue;
    }
    col++;
  }
  std::cout << std::endl << std::endl;
}

//Main <experiment> <configuration> <quasi-replay> <replay-pass> <number-of-calibration-sets> <parameter-override-option> <override-timestamp>; qreplay should only be performed after new offsets obtained
void timewalk( const char *experiment = "gmn", Int_t config = 8, bool qreplay = true, Int_t pass = 0 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  std::string date = util::getDate();
  
  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string tw_rootpath = outdir_path + Form("/hcal_calibrations/pass%d/timing/timewalk_%s_conf%d_qr%d_pass%d.root",pass,experiment,config,(Int_t)qreplay,pass);

  // offset paths and variables
  std::string db_path = gSystem->Getenv("DB_DIR");
  std::string new_tdcoffset_path = Form("../timing/parameters/tdcoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string new_adcgain_path = Form("../energy/parameters/adcgaincoeff_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string new_adctoffset_path = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string new_tw_path = Form("../timing/parameters/timewalk_%s_conf%d_pass%d.txt",experiment,config,pass);
  
  //std::string new_tw_path = "test.txt";

  std::string old_db_path = db_path + "/db_sbs.hcal.dat";
  std::string db_tdcoffset_variable = "sbs.hcal.tdc.offset";
  std::string db_tdccalib_variable = "sbs.hcal.tdc.calib";

  // trial timewalk calibration variables
  std::string db_tdctw_trad_variable = "sbs.hcal.tdc.tradtw"; //a+b/E^n fit tdc v E, p0=b, p1=n, a asymptotes to signal
  std::string db_tdctw_expo_variable = "sbs.hcal.tdc.expotw"; //expo fit tdc v E, p0, p1
  std::string db_tdctw_pol_variable = "sbs.hcal.tdc.poltw"; //pol1 fit slope tdc v E
  std::string db_adcttw_pol_variable = "sbs.hcal.adct.poltw"; //pol1 fit slope adct v E

  std::string db_gain_variable = "sbs.hcal.adc.gain";
  std::string db_adctoffset_variable = "sbs.hcal.adc.timeoffset";

  // Declare output analysis file
  TFile *fout = new TFile( tw_rootpath.c_str(), "RECREATE" );

  //TFile *fout = new TFile( "test.root", "RECREATE" );

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
  
  //Set up structure to hold calibration parameters and initialize size/index
  calset tw_cal[hcal::gNstamp];
  Int_t Ncal_set_size = 0;
  Int_t Ncal_set_idx = 0;
  
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  
  //create output tree
  TTree *P = new TTree("P","Analysis Tree");

  // Calibration set
  Int_t set_out;

  // Timing
  Double_t pblkid_out; //hcal cluster, primary block id
  Double_t tdc_out; //raw hcal cluster tdc time, tree
  Double_t tdc_corr_out; //hcal cluster tdc time, aligned, hodo corrected
  Double_t tdc_tradtw_out; //corr hcal cluster tdc time, qreplay trad fit tw corrections primary cluster primary block
  Double_t tdc_expotw_out; //corr hcal cluster tdc time, qreplay expo fit tw corrections primary cluster primary block
  Double_t tdc_poltw_out; //corr hcal cluster tdc time, qreplay linear fit tw corrections primary cluster primary block
  Double_t atime_out; //hcal cluster adc time, tree
  Double_t atime_corr_out; //hcal cluster adc time, aligned, hodo corrected
  Double_t atime_poltw_out; //corr hcal cluster adc time, qreplay linear fit tw corrections primary cluster primary block
  Double_t hodotmean_out; //hodoscope cluster mean tdc time

  // Physics
  Double_t dx_out; //hcal cluster dx
  Double_t dy_out; //hcal cluster dy
  Double_t W2_out; //Invariant mass squared W2
  Double_t Q2_out; //Inverse momentum squared Q2
  Double_t hcale_out; //hcal cluster energy
  Double_t ep_out; //track reconstructed e' momentum
  Int_t nblk_out; //total blocks in primary cluster
  Int_t mag_out; //sbs magnetic field strength (percent)
  Int_t run_out; //run number
  Int_t tar_out; //target, LH2 or LD2

  // Set output tree branches
  P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
  P->Branch( "tdc", &tdc_out, "tdc/D" );
  P->Branch( "tdc_corr", &tdc_corr_out, "tdc_corr/D" );
  P->Branch( "tdc_tradtw", &tdc_tradtw_out, "tdc_tradtw/D" );
  P->Branch( "tdc_expotw", &tdc_expotw_out, "tdc_expotw/D" );
  P->Branch( "tdc_poltw", &tdc_poltw_out, "tdc_poltw/D" );
  P->Branch( "atime", &atime_out, "atime/D" );
  P->Branch( "atime_corr", &atime_corr_out, "atime_corr/D" );
  P->Branch( "atime_poltw", &atime_poltw_out, "atime_poltw/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "nblk", &nblk_out, "nblk_out/I" );
  P->Branch( "mag", &mag_out, "mag_out/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "tar", &tar_out, "tar_out/I" );
  P->Branch( "set", &set_out, "set/I" );

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  //Declare histograms for modification later where additional calibration sets are necessary
  TH2D *htdcvE[hcal::gNstamp];
  TH2D *hadctvE[hcal::gNstamp];
  TH1D *htesttdc[hcal::gNstamp];
  TH1D *htdc[hcal::gNstamp];
  TH1D *hadct[hcal::gNstamp];
  TH1D *htdc_twtrad[hcal::gNstamp];
  TH1D *htdc_twexpo[hcal::gNstamp];
  TH1D *htdc_twpol[hcal::gNstamp];
  TH1D *hadct_twpol[hcal::gNstamp];

  for( Int_t set=0; set<hcal::gNstamp; set++ ){

    htdcvE[set] = new TH2D(Form("htdcvE_set%d",set),
				    "null",
				    E_bins,
				    E_lower_lim,
				    E_upper_lim,
				    tdc_bins,
				    tdc_lower_lim,
				    tdc_upper_lim);

    hadctvE[set] = new TH2D(Form("hadctvE_set%d",set),
				     "null",
				     E_bins,
				     E_lower_lim,
				     E_upper_lim,
				     adct_bins,
				     adct_lower_lim,
				     adct_upper_lim);

    htesttdc[set] = new TH1D(Form("htesttdc_set%d",set),
			 "null",
			 tdc_bins,
			 tdc_lower_lim,
			 tdc_upper_lim);

    htdc[set] = new TH1D(Form("htdc_set%d",set),
			 "null",
			 tdc_bins,
			 tdc_lower_lim,
			 tdc_upper_lim);

    hadct[set] = new TH1D(Form("hadct_set%d",set),
			  "null",
			  adct_bins,
			  adct_lower_lim,
			  adct_upper_lim);

    htdc_twtrad[set] = new TH1D(Form("htdc_twtrad_set%d",set),
				"null",
				tdc_bins,
				tdc_lower_lim,
				tdc_upper_lim);

    htdc_twexpo[set] = new TH1D(Form("htdc_twexpo_set%d",set),
				"null",
				tdc_bins,
				tdc_lower_lim,
				tdc_upper_lim);

    htdc_twpol[set] = new TH1D(Form("htdc_twpol_set%d",set),
			       "null",
			       tdc_bins,
			       tdc_lower_lim,
			       tdc_upper_lim);

    hadct_twpol[set] = new TH1D(Form("hadct_twpol_set%d",set),
				"null",
				adct_bins,
				adct_lower_lim,
				adct_upper_lim);

  }

  //minimize reporting
  Int_t target_change_index;
  Int_t total_analyzed_events;
  Int_t missing_tdc_events;
  std::string target_marker = "";
  Int_t mag_marker = -1;

  cout << "Setup complete. Beginning loop over runs.." << endl;

  //Main loop over runs
  for (Int_t r=0; r<nruns; r++) {

    // if(r>2)
    //   break;

    //Get run experimental parameters
    std::string current_gain_timestamp = runs[r].adcg_ts;
    std::string current_adctoffset_timestamp = runs[r].adct_ts;
    std::string current_tdcoffset_timestamp = runs[r].tdc_ts;
    std::string current_tdccalib_timestamp = runs[r].tdcc_ts;
    Int_t current_runnumber = runs[r].runnum;

    std::string current_target = runs[r].target;
    std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    Int_t mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A

    //Get run paths
    std::string rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-%s/pass%d/SBS%d/%s/rootfiles/",experiment,pass,config,targ_uppercase.c_str());
    std::string rootfile_path = rootfile_dir + Form("*%d*",current_runnumber);

    //Get target configuration
    SBStarget target_parameters(current_target);
    Int_t target_index = target_parameters.GetTargIndex();
    Double_t target_length = target_parameters.GetTargLength();
    Double_t target_rho = target_parameters.GetTargRho();
    Double_t cell_rho = target_parameters.GetCellRho();
    Double_t cell_diam = target_parameters.GetCellDiam();
    Double_t cell_dEdx = target_parameters.GetCelldEdx();
    Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
    Double_t target_dEdx = target_parameters.GetTargdEdx();
    Double_t M_avg = target_parameters.GetAvgMass();
    
    //send to console if the target changes
    if( target_change_index != target_index ){
      cout << endl << target_parameters;
      target_change_index = target_index;
    }

    ////Timestamp juggling. There is probably a more efficient and clear way to do this.

    //Record new offset/coefficient parameters. Repopulate array on each run for changing timestamps   
    Double_t new_tdc_offsets[hcal::maxHCalChan] = {0.};
    Double_t new_adct_offsets[hcal::maxHCalChan] = {0.};
    Double_t new_adcg_coeff[hcal::maxHCalChan] = {1.};
    vector<Double_t> new_tdctw_tradpar;
    vector<Double_t> new_tdctw_expopar;
    vector<Double_t> new_tdctw_polpar;
    vector<Double_t> new_adcttw_polpar;

    //Get timestamp for new tdc offsets for this run
    std::string tdc_active_timestamp = "-------[ 0000-00-00 00:00:00 ]";
    util::tsCompare( current_tdcoffset_timestamp, config_parameters.GetSBSTimestamp(), tdc_active_timestamp );
    util::tsCompare( tdc_active_timestamp, current_tdccalib_timestamp, tdc_active_timestamp );
    util::readDB( new_tdcoffset_path, tdc_active_timestamp, db_tdcoffset_variable, new_tdc_offsets );

    //get timestamp for new adct offsets for this run
    std::string adct_active_timestamp = "-------[ 0000-00-00 00:00:00 ]";
    util::tsCompare( current_adctoffset_timestamp, config_parameters.GetSBSTimestamp(), adct_active_timestamp );
    util::readDB( new_adctoffset_path, adct_active_timestamp, db_adctoffset_variable, new_adct_offsets );

    //get timestamp for new adc gain coefficients for this run
    std::string active_timestamp = "-------[ 0000-00-00 00:00:00 ]";
    util::tsCompare( current_gain_timestamp, config_parameters.GetSBSTimestamp(), active_timestamp );
    util::readDB( new_adcgain_path, active_timestamp, db_gain_variable, new_adcg_coeff );

    if( qreplay ){
      std::string single_cal_timestamp = config_parameters.GetSBSTimestamp();

      // util::readDB( new_tw_path, active_timestamp, db_tdctw_trad_variable, new_tdctw_tradpar );
      // util::readDB( new_tw_path, active_timestamp, db_tdctw_expo_variable, new_tdctw_expopar );
      // util::readDB( new_tw_path, active_timestamp, db_tdctw_pol_variable, new_tdctw_polpar );
      // util::readDB( new_tw_path, active_timestamp, db_adcttw_pol_variable, new_adcttw_polpar );

      util::readDB( new_tw_path, single_cal_timestamp, db_tdctw_trad_variable, new_tdctw_tradpar );
      util::readDB( new_tw_path, single_cal_timestamp, db_tdctw_expo_variable, new_tdctw_expopar );
      util::readDB( new_tw_path, single_cal_timestamp, db_tdctw_pol_variable, new_tdctw_polpar );
      util::readDB( new_tw_path, single_cal_timestamp, db_adcttw_pol_variable, new_adcttw_polpar );
    }

    if(verbose){
      cout << "TDC offsets from path: " << new_tdcoffset_path << endl;
      coutParams(new_tdc_offsets);
      cout << "ADCt offsets from path: " << new_adctoffset_path << endl;
      coutParams(new_adct_offsets);
      cout << "ADC gain pars from path: " << new_adcgain_path << endl;
      coutParams(new_adcg_coeff);
    }

    //quasi-replay, get new timewalk parameters. Where param_override is not false, use passed params
    // if( qreplay && param_override.compare("false")==0 ){
    //   util::readDB( new_tw_path, active_timestamp, db_tdctw_expo_variable, new_tdctw_expopar );      
    // }else if( qreplay ){
    //   util::readDB( param_override, param_ts, db_tdctw_expo_variable, new_tdctw_expopar );  
    // }

    //////////
    //Select correct calibration set. Timewalk corrections will go with energy, active_timestamp left unaltered.
    bool new_ts = true;
    bool new_offset_ts = true;
    //Look through available sets for duplicate timestamps
    for( Int_t cs=0; cs<Ncal_set_size; cs++ ){ 
      if( tw_cal[cs].timestamp == active_timestamp ){
	new_ts = false;
	Ncal_set_idx = cs;
      }
    }
    
    //if a new timestamp exists, expand the set and update index/old-offsets/old-calib
    if( new_ts || Ncal_set_size==0){

      cout << endl << "Adding new calibration set. Active (energy calibration) timestamp:" << active_timestamp << endl << endl;

      Ncal_set_idx = Ncal_set_size;

      if( Ncal_set_size == hcal::gNstamp-1 ){
	cout << "Error: Allowable calibration sets exceed maximum of " << hcal::gNstamp << ". Check timestamps and adjust limits as necessary." << endl;
	return;
      }

      //update timestamp and old coefficients for corrections. Don't use active timestamp, use db timestamps for old
      tw_cal[Ncal_set_size].timestamp = active_timestamp.c_str();
      util::readDB( old_db_path, current_gain_timestamp, db_gain_variable, tw_cal[Ncal_set_size].old_param ); //gain
      util::readDB( old_db_path, current_tdcoffset_timestamp, db_tdcoffset_variable, tw_cal[Ncal_set_size].old_paramB );  //tdcoffset
      util::readDB( old_db_path, current_adctoffset_timestamp, db_adctoffset_variable, tw_cal[Ncal_set_size].old_paramC );  //adctoffset  
      util::readDB( old_db_path, current_tdccalib_timestamp, db_tdccalib_variable, tw_cal[Ncal_set_size].tdc_calib );  //tdc calib

      if(verbose){
	cout << "OLD ADC gain pars from path: " << old_db_path << endl;
	coutParams(tw_cal[Ncal_set_size].old_param);
	cout << "OLD TDC offsets pars from path: " << old_db_path << endl;
	coutParams(tw_cal[Ncal_set_size].old_paramB);
	cout << "OLD ADCt offsets from path: " << old_db_path << endl;
	coutParams(tw_cal[Ncal_set_size].old_paramC);
	cout << "OLD TDC calib from path: " << old_db_path << endl;
	cout << tw_cal[Ncal_set_size].tdc_calib << endl;
      }


      //build analysis histograms
      htdcvE[Ncal_set_size]->SetTitle(Form("TDC vs E primary block only, set: %d;MeV;ns",Ncal_set_size));
      hadctvE[Ncal_set_size]->SetTitle(Form("ADCt vs E primary block only, set: %d;MeV;ns",Ncal_set_size));
      htdc[Ncal_set_size]->SetTitle(Form("TDC hodo/aligned (no tw), set: %d; ns",Ncal_set_size));
      htesttdc[Ncal_set_size]->SetTitle(Form("TEST TDC hodo/aligned (no tw), set: %d; ns",Ncal_set_size));
      hadct[Ncal_set_size]->SetTitle(Form("ADCt hodo/aligned (no tw), set: %d; ns",Ncal_set_size));
      htdc_twexpo[Ncal_set_size]->SetTitle(Form("TDC timewalk expo corrected, set: %d; ns",Ncal_set_size));
      htdc_twpol[Ncal_set_size]->SetTitle(Form("TDC timewalk pol1 corrected, set: %d; ns",Ncal_set_size));
      hadct_twpol[Ncal_set_size]->SetTitle(Form("ADCt timewalk pol1 corrected, set: %d; ns",Ncal_set_size));

      Ncal_set_size++;
    }     

    //Get available cuts for current config/target/field combination. Use first element (0) of cut
    vector<calcut> cut;
    util::ReadCutList(struct_dir,experiment,config,Ncal_set_idx,pass,current_target,mag,verb,cut);
    if( mag_marker!=mag || target_marker.compare(current_target)!=0 ){
      cout << endl << cut[0];
      mag_marker = mag;
      target_marker = current_target;
    }
    Double_t atime0 = cut[0].atime0; //observed elastic peak in adc time
    Double_t atimesig = cut[0].atime_sig; //observed width of elastic peak in adc time


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
    C->SetBranchStatus( "bb.tr.tg_th", 1 );
    C->SetBranchStatus( "bb.tr.tg_ph", 1 );
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

    //Use TTreeFormula to avoid looping over data an additional time
    TCut GCut = cut[0].gcut.c_str();
    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

    // Set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta_rad, hcalaxes );
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal::HCalvoff*hcalaxes[0];
    Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;

    long nevent = 0, nevents = C->GetEntries(); 
    Int_t treenum = 0, currenttreenum = 0;

    cout << endl << endl << "Looping over run " << r << "/" << nruns << endl;

    //Main loop over events in run
    while(C->GetEntry(nevent++)) {

      cout << "Analyzing run " << current_runnumber << ": " <<  nevent << "/" << nevents << ". Elastics so far: " << total_analyzed_events << " \r";
      cout.flush();

      ///////
      //Single-loop elastic globalcut method. Save pass/fail for output tree.
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
    	treenum = currenttreenum; 
    	GlobalCut->UpdateFormulaLeaves();
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;
	  
      if( failedglobal ) continue;

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
      TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary cluster
      Double_t dx = HCALx - xyhcalexp[0];
      Double_t dy = HCALy - xyhcalexp[1];

      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      Double_t W2mean = cut[0].W2_mean;
      Double_t W2sig = cut[0].W2_sig;
      bool failedw2 = fabs(W2-W2mean)>W2sig; //Observed mean W2 cut on elastic peak

      if( failedw2 ) continue;

      ///////
      //BigBite/SBS Acceptance matching.
      bool failedaccmatch = 
    	xyhcalexp[1] > hcal::posHCalYf ||
    	xyhcalexp[1] < hcal::posHCalYi ||
    	xyhcalexp[0] > hcal::posHCalXf ||
    	xyhcalexp[0] < hcal::posHCalXi;
      
      if( failedaccmatch ) continue;

      //////
      //Cut on bad tdc and keep track of missing primary cluster events
      total_analyzed_events++;
      bool missedtdc = cblktime[0]==-1000 ||
      		       cblktime[0]>10000;
      if( missedtdc )
      	missing_tdc_events++;
	
      ///////
      //Address timewalk with primary cluster block tdc time vs energy
      //primary block id and primary block energy used for both tdc and adct
      Double_t pblkid = cblkid[0]-1; //define primary block, primary cluster ID

      Double_t pblke = cblke[0]/tw_cal[Ncal_set_idx].old_param[(int)pblkid]*new_adcg_coeff[(int)pblkid]; //corrected with new gain parameters

      //ADC time 
      Double_t pblkatime = cblkatime[0];
      //align atime with new params, paramC is adct offset
      Double_t hcaladct = pblkatime+tw_cal[Ncal_set_idx].old_paramC[(int)pblkid]-new_adct_offsets[(int)pblkid]; 
      //correct the block adc time with hodoscope to mitigate trigger jitter
      Double_t hcaladct_hodo = hcaladct - HODOtmean;
      //Double_t atime0 = cut[0].atime0; //observed elastic peak in adc time
      //Double_t atimesig = cut[0].atime_sig; //observed width of elastic peak in adc time
      bool failedcoin = abs(hcaladct-atime0)>adct_Nsig*atimesig; //may use this later

      //TDC time
      Double_t pblktime = cblktime[0];
      //align block tdc time with new offset parameters, paramB is tdc offset
      Double_t pblktime_aligned = pblktime + 
    	tw_cal[Ncal_set_idx].old_paramB[(int)pblkid] * tw_cal[Ncal_set_idx].tdc_calib - 
    	new_tdc_offsets[(int)pblkid] * tw_cal[Ncal_set_idx].tdc_calib;
      //correct the block tdc time with hodoscope to mitigate trigger jitter
      Double_t pblktime_hodo = pblktime_aligned - HODOtmean;

      //correct the primary block times with loaded timewalk corrections on qreplay
      Double_t pblktime_tw_trad = 0.;
      Double_t pblktime_tw_expo = 0.;
      Double_t pblktime_tw_pol = 0.;
      Double_t hcaladct_tw_pol = 0.;
      if( qreplay ){
    	Double_t trad_corr = new_tdctw_tradpar[0]/pow( (pblke+minE), new_tdctw_tradpar[1] );
    	pblktime_tw_trad = pblktime_hodo - trad_corr;
    	htdc_twtrad[Ncal_set_idx]->Fill( pblktime_tw_trad );

    	Double_t expo_corr = new_tdctw_expopar[0]*exp( -new_tdctw_expopar[1]*pblke );
    	pblktime_tw_expo = pblktime_hodo - expo_corr;
    	htdc_twexpo[Ncal_set_idx]->Fill( pblktime_tw_expo );

    	Double_t tdc_pol_corr = new_tdctw_polpar[0]*pblke;
    	pblktime_tw_pol = pblktime_hodo - tdc_pol_corr;
    	htdc_twpol[Ncal_set_idx]->Fill( pblktime_tw_pol );

    	Double_t adct_pol_corr = new_adcttw_polpar[0]*pblke;
    	hcaladct_tw_pol = hcaladct_hodo - adct_pol_corr;
    	hadct_twpol[Ncal_set_idx]->Fill( hcaladct_tw_pol );
      }

      //For clarity, reset all values to -1000 if the tdc was missed on this event
      if( missedtdc ){
	pblktime = -1000;
	pblktime_aligned = -1000;
	pblktime_hodo = -1000;
	pblktime_tw_trad = -1000;
	pblktime_tw_expo = -1000;
	pblktime_tw_pol = -1000;
      }

      //build tdc vs E from primary block to extract timewalk parameters
      htdcvE[Ncal_set_idx]->Fill( pblke, pblktime_hodo );

      //build tdc vs E from primary block to extract timewalk parameters
      hadctvE[Ncal_set_idx]->Fill( pblke, hcaladct_hodo );

      //build tdc general histos for quality check
      htdc[Ncal_set_idx]->Fill( pblktime_hodo );
      htesttdc[Ncal_set_idx]->Fill( cblkatime[0] );
      hadct[Ncal_set_idx]->Fill( hcaladct_hodo );

      //fill analysis tree
      pblkid_out = pblkid;
      tdc_out = pblktime;
      tdc_corr_out = pblktime_hodo;
      tdc_tradtw_out = pblktime_tw_trad;
      tdc_expotw_out = pblktime_tw_expo;
      tdc_poltw_out = pblktime_tw_pol;
      atime_out = pblkatime;
      atime_corr_out = hcaladct_hodo;
      atime_poltw_out = hcaladct_tw_pol;
      hcale_out = HCALe;
      dx_out = dx;
      dy_out = dy;
      W2_out = W2;
      Q2_out = Q2;
      nblk_out = nblk;
      mag_out = mag;
      run_out = current_runnumber;
      tar_out = target_index;
      hodotmean_out = HODOtmean;
      
      P->Fill();
      
    }//end loop over event
    
    // getting ready for the next run
    C->Reset();
    
  }//endloop over runs
  
  // set up arrays for timewalk fits
  Double_t tdcvE_tradP1 = 0.;
  Double_t tdcvE_tradP2 = 0.;
  Double_t tdcvE_expoP0 = 0.;
  Double_t tdcvE_expoP1 = 0.;
  Double_t tdcvE_pol = 0.;
  Double_t adctvE_pol = 0.;

  // Declare canvases
  TCanvas *c1 = new TCanvas("c1","TDC_tradfit",1600,1200);
  TCanvas *c2 = new TCanvas("c2","TDC_expofit",1600,1200);
  TCanvas *c3 = new TCanvas("c3","TDC_linearfit",1600,1200);
  TCanvas *c4 = new TCanvas("c4","ADCt_linearfit",1600,1200);

  // set up clones of time v E for various fits
  TH2D *htdc_tradclone = (TH2D*)htdcvE[twset]->Clone("htdc_tradclone");
  TH2D *htdc_expoclone = (TH2D*)htdcvE[twset]->Clone("htdc_expoclone");
  TH2D *htdc_polclone = (TH2D*)htdcvE[twset]->Clone("htdc_polclone");
  TH2D *hadct_polclone = (TH2D*)hadctvE[twset]->Clone("hadct_polclone");

  ////
  //traditional tdc fit

  std::cout << "working on traditional tdc fit.." << std::endl;
  
  c1->cd();

  //Fit the TDC vs E plot with traditional timewalk fit
  TF1 *fittdc_trad = new TF1( "fittdc_trad", util::g_tradtwfit, tw_fit_llim, tw_fit_ulim, 3 );

  fittdc_trad->SetParameters(( trad_p0_ulim + trad_p0_llim )/2,
  			     ( trad_p1_ulim + trad_p1_llim )/2,
  			     ( trad_p2_ulim + trad_p2_llim )/2);

  fittdc_trad->SetParLimits(0,trad_p0_llim,trad_p0_ulim);
  fittdc_trad->SetParLimits(1,trad_p1_llim,trad_p1_ulim);
  fittdc_trad->SetParLimits(2,trad_p2_llim,trad_p2_ulim);

  //fittdc_trad->SetParLimits(2,0.3,1.8);

  htdc_tradclone->Fit("fittdc_trad","WQ","",tw_fit_llim,tw_fit_ulim);

  tdcvE_tradP1 = fittdc_trad->GetParameter(1);
  tdcvE_tradP2 = fittdc_trad->GetParameter(2);

  htdc_tradclone->SetTitle(Form("tdc tradfit P0:%0.2f P1:%0.2f P2:%0.2f",fittdc_trad->GetParameter(0),tdcvE_tradP1,tdcvE_tradP2));  
  htdc_tradclone->Draw("colz");

  c1->Write();

  ////
  //exponential tdc fit

  c2->cd();
  
  cout << "working on exponential tdc fit.." << endl;
  
  //Fit the TDC vs E plot with exponential fit
  TF1 *fittdc_expo = new TF1( "fittdc_expo", util::g_twfit, tw_fit_llim, tw_fit_ulim, 3 );

  fittdc_expo->SetParameters((tw_p0_ulim+tw_p0_llim)/2,
			     (tw_p1_ulim+tw_p1_llim)/2,
			     (tw_p2_ulim+tw_p2_llim)/2);
  
  fittdc_expo->SetParLimits(0,tw_p0_llim,tw_p0_ulim);
  fittdc_expo->SetParLimits(1,tw_p1_llim,tw_p1_ulim);
  fittdc_expo->SetParLimits(2,tw_p2_llim,tw_p2_ulim);
  
  htdc_expoclone->Fit("fittdc_expo","WQ","",tw_fit_llim,tw_fit_ulim);
  
  tdcvE_expoP0 = fittdc_expo->GetParameter(0);
  tdcvE_expoP1 = fittdc_expo->GetParameter(1);
  
  htdc_expoclone->SetTitle(Form("tdc expofit P0:%0.2f P1:%0.2f P2:%0.2f",tdcvE_expoP0,tdcvE_expoP1,fittdc_expo->GetParameter(2)));
  htdc_expoclone->Draw("colz");
  
  c2->Write();
    
  ////
  //linear tdc fit

  c3->cd();
  
  cout << "working on linear tdc fit.." << endl;
  
  //Fit the TDC vs E plots with pol1
  TF1 *fittdc_pol = new TF1( "fittdc_pol", util::g_lfit, fit_ll, fit_ul, 2 );
  fittdc_pol->SetParameters(p0_set,p1_set);
  
  fittdc_pol->SetParLimits(0,p0_ll,p0_ul);
  fittdc_pol->SetParLimits(1,p1_ll,p1_ul);

  htdc_polclone->Fit("fittdc_pol","WQ","",fit_ll,fit_ul);
  tdcvE_pol = fittdc_pol->GetParameter(1);

  htdc_polclone->SetTitle(Form("tdc linearfit P0:%0.2f P1:%0.2f",fittdc_pol->GetParameter(0),tdcvE_pol));

  htdc_polclone->Draw("colz");

  c3->Write();

  ////
  //linear adct fit

  c4->cd();

  cout << "working on linear adct fit.." << endl;

  //Fit the ADCt vs E plots with pol1
  TF1 *fitadct_pol = new TF1( "fitadct_pol", util::g_lfit, afit_ll, afit_ul, 2 );
  fitadct_pol->SetParameters(ap0_set,ap1_set);

  fitadct_pol->SetParLimits(0,ap0_ll,ap0_ul);
  fitadct_pol->SetParLimits(1,ap1_ll,ap1_ul);

  hadct_polclone->Fit("fitadct_pol","WQ","",afit_ll,afit_ul);
  adctvE_pol = fitadct_pol->GetParameter(1);

  hadct_polclone->SetTitle(Form("adct linearfit P0:%0.2f P1:%0.2f",fitadct_pol->GetParameter(0),adctvE_pol));
  hadct_polclone->Draw("colz");

  c4->Write();

  cout << "moving to write parameters to file.." << endl;

  // Create timewalk fit parameter files
  ofstream tw;
  if( qreplay )
    new_tw_path = "/dev/null"; //Safety to prevent overwriting constants on quasi-replay

  if( !qreplay ){
    
    //check if timestamp for calibration set is newer than config timestamp. Use the newer of the two.
    std::string active_timestamp;
    util::tsCompare(sbs_timestamp,tw_cal[twset].timestamp,active_timestamp);
        
    tw.open( new_tw_path );

    tw << endl << endl << "#HCal tdc vs E fit parameters obtained " << date.c_str() << endl;

    tw << endl << active_timestamp << endl << endl;

    tw << "#Traditional tdc vs E fit for all PMTs -> y = P0 + P1/( x^P2 ). P0 normal to signal." << endl;
	      
    tw << db_tdctw_trad_variable << " = ";
    tw << tdcvE_tradP1 << " ";
    tw << tdcvE_tradP2 << endl << endl;

    tw << "#Exponential tdc vs E fit for all PMTs -> y = P0*exp(-P1*x) + P2. P2 normal to signal." << endl;
	      
    tw << db_tdctw_expo_variable << " = ";
    tw << tdcvE_expoP0 << " ";
    tw << tdcvE_expoP1 << endl << endl;

    tw << endl << "#First order polynomial tdc vs E fit for all PMTs -> P0 + P1*x. P0 normal to signal." << endl;
      
    tw << db_tdctw_pol_variable << " = ";
    tw << tdcvE_pol << endl << endl;

    tw << endl << "#First order polynomial adc time vs E fit for all PMTs -> P0 + P1*x. P0 normal to signal." << endl;
      
    tw << db_adcttw_pol_variable << " = ";
    tw << adctvE_pol << endl << endl;
    
    tw << endl << endl; //needed for now to prevent bug in util::readDB
    
  }

  tw.close();

  //clean up
  for( Int_t unused_set_index=Ncal_set_size; unused_set_index<hcal::gNstamp; unused_set_index++ ){
    delete htdcvE[unused_set_index];
    delete hadctvE[unused_set_index];
    delete htdc[unused_set_index];
    delete hadct[unused_set_index];
    delete htdc_twtrad[unused_set_index];
    delete htdc_twexpo[unused_set_index];
    delete htdc_twpol[unused_set_index];
    delete hadct_twpol[unused_set_index];
  }

  cout << endl << "Writing parameters to " << new_tw_path << endl << endl;
  cout << "Writing output root file to " << tw_rootpath << endl << endl;

  fout->Write();

  st->Stop();

  // Send missing primary cluster tdc events to console for reference
  std::cout << "Missing TDC signals on all elastics analyzed: " << missing_tdc_events << "/" << total_analyzed_events << std::endl;

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;    

}

