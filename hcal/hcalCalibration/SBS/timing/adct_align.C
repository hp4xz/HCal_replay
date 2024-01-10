//SSeeds 5.8.23 - Script adapted from tcal.C at https://github.com/sebastianseeds/HCal_replay to align adct elastic peaks across all HCal block ids
// 5.23.23 update - Added adct all-channels distribution and fit for adct sigma extraction and calibration of sbs.hcal.tmax
// 7.31.23 update - Added tracking of timing peaks for two sample channels over runs to monitor timing shifts evident in GMn. Shifted target ADCt value to zero. Staged shift to timing difference between HCal - BBCal atime to improve timing resolution for tighter elastic selection and added branch to output analysis tree.
//NOTE: requires $DB_DIR and $OUT_DIR paths set correctly in current environment. Script assumes hardware differences on adctoffset timestamps and outputs alignment constants for all timestamps within the configuration provided.

#include <vector>
#include <iostream>
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

// ADCt offset extraction constraints
const Int_t first_hcal_chan = 0;
const Int_t total_bins = 320;
const Int_t lower_lim = -60;
const Int_t upper_lim = 100;
// const Int_t total_bins = 600;
// const Int_t lower_lim = -150;
// const Int_t upper_lim = 150;
const Int_t fit_event_min = 50; //Minimum events in histogram to fit
const Double_t observed_adct_sigma = 4.0; //rough estimate
const Double_t Nsig = 2;
const Double_t ADCt_target = 0.; //set the adc time target (ns)
const Double_t tmax_factor = 6.; //number of sigma to consider for tmax cut on cluster elements
const Int_t linecount = 12;
const Int_t samples[2] = {161,80}; //Sample channels to check timing alignment over runs within calibrations set

bool wonkset = false; //setup for run range exceptions

bool tempfix = false; //temporary adjustment of binmaxX for large variation on small data sets
bool channelCorrect = true; //temporary adjustment of single adct values on channel to fit in gen window
Int_t corrChan = 264; //Channel to correct
const Double_t tempChanOffset = 330.; //temporary offset for corrChan to return to general window

void overlayWithGaussianFits(TH2D *h2d, TCanvas *canvas) {
  // switch to canvas location
  canvas->cd();

  // Create a TGraphErrors to hold the means and standard deviations
  TGraphErrors *graphErrors = new TGraphErrors(h2d->GetNbinsX());

  // Loop over the bins in the TH2D
  for (int i = 1; i <= h2d->GetNbinsX(); ++i) {
    // Get the projection of this bin along the y-axis
    TH1D *projY = h2d->ProjectionY("_py", i, i);

    //Skip the fit if less than 100 entries exist in the bin
    if( projY->GetEntries()<100 )
      continue;

    //declare some dynamic fit variables
    Int_t binMax = projY->GetMaximumBin();
    Double_t binCenter = projY->GetBinCenter( binMax );
    Double_t fitLowerLim = binCenter - 2*observed_adct_sigma;
    Double_t fitUpperLim = binCenter + 2*observed_adct_sigma;

    // Fit a Gaussian to this projection
    TF1 *gausFit = new TF1("gausFit", "gaus");
    projY->Fit(gausFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode

    // Get the mean and standard deviation from the fit
    double mean = gausFit->GetParameter(1);
    double stdDev = gausFit->GetParameter(2);

    // Get the error on the mean and standard deviation
    double meanError = gausFit->GetParError(1);
    double stdDevError = gausFit->GetParError(2);

    // Set this point in the TGraphErrors
    graphErrors->SetPoint(i - 1, h2d->GetXaxis()->GetBinCenter(i), mean);
    graphErrors->SetPointError(i - 1, 0, stdDev); // Assuming no error on x
  }

  // Draw the TH2D
  h2d->Draw("COLZ");

  // Customize the appearance of the TGraphErrors
  graphErrors->SetMarkerStyle(20);
  graphErrors->SetMarkerColor(kBlack);
  graphErrors->SetLineColor(kBlack);

  // Draw the TGraphErrors on the same canvas
  graphErrors->Draw("P SAME");

  // Update the canvas
  gPad->Update();
}

//Main <experiment> <configuration> <quasi-replay> <replay-pass> <target_option>; qreplay should only be performed after new offsets obtained
void adct_align( const char *experiment = "gmn", Int_t config=4, bool qreplay = true, Int_t pass = 0, bool h2only = false ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  std::string date = util::getDate();

  // outfile path
  std::string h2opt = "";
  if( h2only )
    h2opt = "_lh2only";

  std::string wonkopt = "";
  if( config==11 || config==14 ){
    if( wonkset )
      wonkopt = "_wonkset";
    else
      wonkopt = "_normalset";
  }

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string adctalign_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/adctalign%s_%s_conf%d_qr%d_pass%d%s.root",pass,h2opt.c_str(),experiment,config,(Int_t)qreplay,pass,wonkopt.c_str());
  std::string plotdir = Form("../quality_plots/%s/conf%d/",experiment,config);

  // offset paths and variables
  //std::string db_path = gSystem->Getenv("DB_DIR");
  std::string db_path = "/w/halla-scshelf2102/sbs/seeds/alt_sbsreplay/SBS-replay/DB";

  std::string new_adctoffset_path = Form("parameters/adctoffsets_%s_conf%d_pass%d%s.txt",experiment,config,pass,wonkopt.c_str());
  std::string new_tmax_path = Form("parameters/adcttmax_%s_conf%d_pass%d%s.txt",experiment,config,pass,wonkopt.c_str());
  std::string old_db_path = db_path + "/db_sbs.hcal.dat";
  std::string db_adctoffset_variable = "sbs.hcal.adc.timeoffset";
  std::string db_tmax_variable = "sbs.hcal.tmax";

  //std::string testpath = "adcttest.root";

  // Set up output analysis file
  TFile *fout = new TFile( adctalign_path.c_str(), "RECREATE" );

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  Int_t verb = 0; //Don't print diagnostic info by default
  
  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
 
  //Set up structure to hold calibration parameters and initialize size/index
  calset adct_cal[hcal::gNstamp];
  Int_t Ncal_set_size = 0;
  Int_t Ncal_set_idx = 0;
    
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  
  // create output tree
  TTree *P = new TTree("P","Analysis Tree");

  // Timing
  Double_t pblkid_out; //hcal cluster, primary block id
  Double_t tdc_out; //hcal cluster tdc time, tree
  Double_t atime_out; //hcal cluster adc time, tree
  Double_t hodotmean_out; //hodoscope cluster mean tdc time
  Double_t atime_hc_out; //hodoscope corrected hcal adct
  Double_t atime_bc_out; //bbcal corrected hcal adct

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

  //cluster block variables
  Double_t cblkid_out[hcal::maxHCalBlk];
  Double_t cblkatime_out[hcal::maxHCalBlk];
  Double_t cblktime_out[hcal::maxHCalBlk];
  Double_t cblke_out[hcal::maxHCalBlk];

  // Set output tree branches
  P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
  P->Branch( "tdc", &tdc_out, "tdc/D" );
  P->Branch( "atime", &atime_out, "atime/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "atime_hc", &atime_hc_out, "atime_hc/D" );
  P->Branch( "atime_bc", &atime_bc_out, "atime_bc/D" );
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "nblk", &nblk_out, "nblk_out/I" );
  P->Branch( "mag", &mag_out, "mag_out/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "tar", &tar_out, "tar_out/I" );

  P->Branch( "cblkid", &cblkid_out, Form("cblkid[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblkatime", &cblkatime_out, Form("cblkatime[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblktime", &cblktime_out, Form("cblktime[%d]/D",hcal::maxHCalBlk) );
  P->Branch( "cblke", &cblke_out, Form("cblke[%d]/D",hcal::maxHCalBlk) );

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  // Declare histograms for modification later where additional calibration sets are necessary
  TH2D *hap_hodocorr_ID[hcal::gNstamp];
  TH2D *hadct_samp1_run[hcal::gNstamp];
  TH2D *hadct_samp2_run[hcal::gNstamp];
  TH1D *hadctblk[hcal::gNstamp];
  TH1D *hadct[hcal::gNstamp];

  for( Int_t set=0; set<hcal::gNstamp; set++ ){
    hap_hodocorr_ID[set] = new TH2D(Form("hap_hodocorr_ID_set%d",set),
				    "null",
				    hcal::maxHCalChan,
				    first_hcal_chan,
				    hcal::maxHCalChan,
				    total_bins,
				    lower_lim,
				    upper_lim);

    //Assumes that the runlist is ordered in ascending run numbers
    hadct_samp1_run[set] = new TH2D(Form("hadct_s%d_run_set%d",samples[0],set),
				   "null",
				   (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				   runs[0].runnum-1,
				   runs[nruns-1].runnum+1,
				   total_bins,
				   lower_lim,
				   upper_lim);

    hadct_samp2_run[set] = new TH2D(Form("hadct_s%d_run_set%d",samples[1],set),
				   "null",
				   (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				   runs[0].runnum-1,
				   runs[nruns-1].runnum+1,
				   total_bins,
				   lower_lim,
				   upper_lim);
    
    hadctblk[set] = new TH1D(Form("hadctblk_set%d",set),
			  "null",
			  total_bins,
			  lower_lim,
			  upper_lim);

    
    hadct[set] = new TH1D(Form("hadct_set%d",set),
			  "null",
			  total_bins,
			  lower_lim,
			  upper_lim);


  }

  //minimize reporting
  Int_t target_change_index;

  //Main loop over runs
  for (Int_t r=0; r<nruns; r++) {

    //Get run experimental parameters
    std::string current_offset_timestamp = runs[r].adct_ts;
    Int_t current_runnumber = runs[r].runnum;

    //cout << wonkset << " " << current_runnumber << endl;

    //TEMPORARY FIX for pass2 to account for timing shift
    if( config==11 ){
      if( wonkset ){
    	if( current_runnumber<12450 || current_runnumber>12860 )
    	  continue;
      }else{
    	if( current_runnumber>12450 && current_runnumber<12860 )
	  //cout << current_runnumber << endl;
    	  continue;
      }
    }else if( config==14 ){
      if( wonkset ){
    	if( current_runnumber>13260 )
    	  continue;
      }else{
    	if( current_runnumber<13260 )
    	  continue;
      }
    }

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
      cout << target_parameters;
      target_change_index = target_index;
    }

    //Record new offset parameters from qreplay==false. Repopulate array on each run for changing timestamps
    Double_t new_adct_offsets[hcal::maxHCalChan] = {0.};
    if( qreplay ){ //Check which timestamp is newest (DB:adctoffset, DB:tdccalib, configstart) and get new offsets
      std::string active_timestamp;
      util::tsCompare( current_offset_timestamp, config_parameters.GetSBSTimestamp(), active_timestamp );
      util::readDB( new_adctoffset_path, active_timestamp, db_adctoffset_variable, new_adct_offsets );
    }

    //////////
    //Select correct calibration set.
    bool new_offset_ts = true;
    //Look through available sets for duplicate timestamps
    for( Int_t cs=0; cs<Ncal_set_size; cs++ ){ 

      //if timestamp found in previous or current set, update the calibration index
      if( adct_cal[cs].timestamp == current_offset_timestamp ){
	new_offset_ts = false;
	Ncal_set_idx = cs;
      }

    }

    //if a new timestamp exists, expand the set and update index/old-offsets
    if( new_offset_ts || Ncal_set_size==0){

      cout << "Adding new calibration set. Offset ts:" << current_offset_timestamp << endl << endl;

      Ncal_set_idx = Ncal_set_size;

      if( Ncal_set_size == hcal::gNstamp-1 ){
	cout << "Error: Allowable calibration sets exceed maximum of " << hcal::gNstamp << ". Check timestamps and adjust limits as necessary." << endl;
	return;
      }

      adct_cal[Ncal_set_size].timestamp = current_offset_timestamp.c_str();
      util::readDB( old_db_path, current_offset_timestamp, db_adctoffset_variable, adct_cal[Ncal_set_size].old_param ); 
      hap_hodocorr_ID[Ncal_set_size]->SetTitle(Form("ADCt Primary Block - TDC hodo, offset ts: %s;Channel;ADCt_{HCAL}-TDC_{HODO} (ns)",adct_cal[Ncal_set_size].timestamp.c_str()));
      hadct_samp1_run[Ncal_set_size]->SetTitle(Form("ADCt vs Run, offset ts: %s;Run;ADCt_{HCAL} (ns)",adct_cal[Ncal_set_size].timestamp.c_str()));
      hadct_samp2_run[Ncal_set_size]->SetTitle(Form("ADCt vs Run, offset ts: %s;Run;ADCt_{HCAL} (ns)",adct_cal[Ncal_set_size].timestamp.c_str()));
      hadct[Ncal_set_size]->SetTitle(Form("ADCt All Blocks Primary Cluster, offset ts: %s;Raw Cluster Block ADCt_{HCAL} (ns)",adct_cal[Ncal_set_size].timestamp.c_str()));
      hadctblk[Ncal_set_size]->SetTitle(Form("ADCt Primary Block Primary Cluster, offset ts: %s;Raw Cluster Block ADCt_{HCAL} (ns)",adct_cal[Ncal_set_size].timestamp.c_str()));

      Ncal_set_size++;
    }     

    //Get available cuts for current config/target/field combination. Use first element (0) of cut
    vector<calcut> cut;
    util::ReadCutList(struct_dir,experiment,config,Ncal_set_idx,pass,current_target,mag,verb,cut);

    // Setting up chain and branch addresses
    C = new TChain("T");
    C->Add(rootfile_path.c_str());

    C->SetBranchStatus("*",0);    
    Double_t BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
    Double_t BBtr_vz[hcal::maxTracks];
    Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBps_atime, BBsh_x, BBsh_y, BBsh_e;	
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
    C->SetBranchStatus( "bb.ps.atimeblk", 1 );
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
    C->SetBranchAddress( "bb.ps.atimeblk", &BBps_atime );
    C->SetBranchAddress( "bb.sh.e", &BBsh_e );
    C->SetBranchAddress( "bb.sh.x", &BBsh_x );
    C->SetBranchAddress( "bb.sh.y", &BBsh_y ); 
    C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Ncid ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

    //Globalcut enables
    C->SetBranchStatus( "bb.tr.tg_th", 1 );
    C->SetBranchStatus( "bb.tr.tg_ph", 1 );

    //Use TTreeFormula to avoid looping over data an additional time
    TCut GCut = cut[0].gcut.c_str();

    //Add globalcut and elastic cuts for reporting
    adct_cal[Ncal_set_idx].gcut = cut[0].gcut;
    adct_cal[Ncal_set_idx].minEv = fit_event_min;

    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

    // Set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta_rad, hcalaxes );
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal::HCalvoff*hcalaxes[0];
    Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;

    long nevent = 0, nevents = C->GetEntries(); 
    Int_t treenum = 0, currenttreenum = 0;

    //Main loop over events in run
    while (C->GetEntry(nevent++)) {

      cout << "Analyzing run (" << r << "/" << nruns << ") " << current_runnumber << ": " <<  nevent << "/" << nevents << " \r";
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
      TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
      Double_t dx = HCALx - xyhcalexp[0];
      Double_t dy = HCALy - xyhcalexp[1];

      ///////
      //BigBite/SBS Acceptance matching.
      bool failedaccmatch = 
	xyhcalexp[1] > hcal::posHCalYf ||
	xyhcalexp[1] < hcal::posHCalYi ||
	xyhcalexp[0] > hcal::posHCalXf ||
	xyhcalexp[0] < hcal::posHCalXi;
	  
      if( failedaccmatch ) continue;

      Double_t pblkid = (Double_t)cblkid[0]-1;

      Double_t hcalatime = cblkatime[0];
      if( qreplay ) 
	hcalatime = hcalatime + adct_cal[Ncal_set_idx].old_param[(Int_t)pblkid] - new_adct_offsets[(Int_t)pblkid];
      hadct[Ncal_set_idx]->Fill(hcalatime);

      //Double_t adct_tc = hcalatime-BBps_atime; //Primary cluster, primary block tdc with trigger correction (ns)
      Double_t adct_tc = hcalatime-HODOtmean; //Primary cluster, primary block tdc with trigger correction (ns)

      //Fill primary timing alignment histogram. Will use hodoscope corrected time for cancellation with BBCal atime later
      //hap_hodocorr_ID[Ncal_set_idx]->Fill( pblkid, adct_tc );
      
      //Correct channel way out of window 
      if(pblkid==corrChan&&channelCorrect){
	hap_hodocorr_ID[Ncal_set_idx]->Fill( pblkid, adct_tc+tempChanOffset );
      }else{
	hap_hodocorr_ID[Ncal_set_idx]->Fill( pblkid, adct_tc );
      }

      //Fill single channel histograms for comparison over many runs
      if(pblkid==samples[0])
	hadct_samp1_run[Ncal_set_idx]->Fill(current_runnumber,hcalatime);
      if(pblkid==samples[1])
	hadct_samp2_run[Ncal_set_idx]->Fill(current_runnumber,hcalatime);

      //Fill the adc time for all blocks in the primary cluster for tmax extraction. Note that no hodoscope correction is applied as it is not included in SBS-offline and on qreplay new offsets are used
      for( Int_t blk = 0; blk<nblk; blk++ ){
	Double_t blkid = (Double_t)cblkid[blk]-1;
	Double_t blkatime = cblkatime[blk];
	if( qreplay ) 
	  blkatime = blkatime + adct_cal[Ncal_set_idx].old_param[(Int_t)blkid] - new_adct_offsets[(Int_t)blkid];
	hadctblk[Ncal_set_idx]->Fill(blkatime);
      }

      // pblkid_out = pblkid;
      // tdc_out = cblktime[0];
      // atime_out = hcalatime;
      // atime_hc_out = hcalatime-HODOtmean;
      // atime_bc_out = hcalatime-BBps_atime;
      // hcale_out = HCALe;
      // dx_out = dx;
      // dy_out = dy;
      // W2_out = W2;
      // Q2_out = Q2;
      // nblk_out = nblk;   
      // mag_out = mag;
      // run_out = current_runnumber;
      // tar_out = target_index;
      // hodotmean_out = HODOtmean;

      P->Fill();
	  
    }//end loop over event

    // getting ready for the next run
    C->Reset();

  }//endloop over runs

  cout << "Ended loop over events. Proceeding to fits.." << endl;

  // Declare canvas and graph for plots
  TCanvas *ADCt_top[Ncal_set_size];
  TCanvas *ADCt_bot[Ncal_set_size];
  TCanvas *ADCtblk_all[Ncal_set_size];
  TCanvas *ADCt_all[Ncal_set_size];
  TCanvas *samp1_run[Ncal_set_size];
  TCanvas *samp2_run[Ncal_set_size];

  //No error on cell location
  Double_t cellerr[hcal::maxHCalChan] = {0.};

  //Make arrays for adct tgraphs
  Double_t tcell[Ncal_set_size][hcal::maxHCalChan];
  Double_t tcval[Ncal_set_size][hcal::maxHCalChan];
  Double_t tcvalw[Ncal_set_size][hcal::maxHCalChan];
  Double_t tcerr[Ncal_set_size][hcal::maxHCalChan];
  TH1D *tcellslice[Ncal_set_size][hcal::maxHCalChan];

  //Get averages for empty channels
  Double_t tcval_avg[Ncal_set_size];
  Int_t tcval_Ng[Ncal_set_size];

  //Get sigma for tmax/tcut calibrations
  Double_t allblksig[Ncal_set_size];
  Double_t allsig[Ncal_set_size];

  //Define offset calibration set
  std::string offset_timestamp;

  for( Int_t set=0; set<Ncal_set_size; set++ ){
        
    //Define calibration set by newest timestamp (config/tdcoffset/tdccalib) and record constant
    offset_timestamp = adct_cal[set].timestamp;
    std::string config_timestamp = config_parameters.GetSBSTimestamp();
    std::string active_timestamp;
    util::tsCompare( offset_timestamp, config_timestamp, active_timestamp );
    
    tcval_avg[set] = 0.;
    tcval_Ng[set] = 0;

    //Set up canvases for quality checks
    ADCtblk_all[set] = new TCanvas(Form("adct_tmax_set%d",set),Form("ADCtblk_all, timestamp: %s",active_timestamp.c_str()),1600,1200); 
    ADCt_all[set] = new TCanvas(Form("adct_tcut_set%d",set),Form("ADCt_all, timestamp: %s",active_timestamp.c_str()),1600,1200); 
    ADCt_top[set] = new TCanvas(Form("adct_fit_top_set%d",set),Form("ADCt_top, timestamp: %s",active_timestamp.c_str()),1600,1200);
    ADCt_bot[set] = new TCanvas(Form("adct_fit_bot_set%d",set),Form("ADCt_bot, timestamp: %s",active_timestamp.c_str()),1600,1200);
    samp1_run[set] = new TCanvas(Form("samp1_run_set%d",set),Form("ADCt vs Run, Channel: %d, ts:%s",samples[0],active_timestamp.c_str()),1600,1200);
    samp2_run[set] = new TCanvas(Form("samp2_run_set%d",set),Form("ADCt vs Run, Channel: %d, ts:%s",samples[1],active_timestamp.c_str()),1600,1200); 
    ADCt_top[set]->Divide(12,12);
    ADCt_bot[set]->Divide(12,12);
    gStyle->SetOptStat(0);
    
    //Fit channels and overlay for run-to-run comparison
    overlayWithGaussianFits(hadct_samp1_run[set],samp1_run[set]);
    samp1_run[set]->Write();
    overlayWithGaussianFits(hadct_samp2_run[set],samp2_run[set]);
    samp2_run[set]->Write();

    //Fit all channel all block adct and extract sigma for tmax calibration
    ADCtblk_all[set]->cd();
    Double_t allblkfitl = 0.;
    Double_t allblkfith = 0.;

    Int_t adctblkN = hadctblk[set]->GetEntries();
    Int_t allblkbinmax = hadctblk[set]->GetMaximumBin();
    Double_t allblkbinmaxX = lower_lim+allblkbinmax*(upper_lim-lower_lim)/total_bins;
    Double_t allblkbinmaxY = hadctblk[set]->GetBinContent(allblkbinmax);
    allblkfitl = allblkbinmaxX - 2*observed_adct_sigma;
    allblkfith = allblkbinmaxX + 2*observed_adct_sigma;

    TF1 *gausfitallblk = new TF1("gausfitallblk",util::g_gfit_bg,allblkfitl,allblkfith,5);
    gausfitallblk->SetLineWidth(4);
    gausfitallblk->SetParameter(0,allblkbinmaxY);
    gausfitallblk->SetParameter(1,allblkbinmaxX);
    gausfitallblk->SetParLimits(1,allblkfitl,allblkfith);
    gausfitallblk->SetParameter(2,observed_adct_sigma);
    gausfitallblk->SetParLimits(2,1.,3*observed_adct_sigma);
    gausfitallblk->SetParameter(3,25);
    gausfitallblk->SetParameter(4,20);
    gausfitallblk->SetParLimits(3,0,200);
    
    hadctblk[set]->Fit("gausfitallblk","RBMQ");
    hadctblk[set]->Draw();

    Double_t allblkmean = gausfitallblk->GetParameter(1);
    allblksig[set] = gausfitallblk->GetParameter(2);
    hadctblk[set]->SetTitle(Form("All Blk Set:%d QR:%d N:%d Mean:%f Sigma:%f",set,(Int_t)qreplay,adctblkN,allblkmean,allblksig[set]));    
    //hadctblk[set]->Write();
    ADCtblk_all[set]->Write();

    //Fit all channel adct and extract sigma for tcut calibration
    ADCt_all[set]->cd();
    Double_t allfitl = 0.;
    Double_t allfith = 0.;

    Int_t adctN = hadct[set]->GetEntries();
    Int_t allbinmax = hadct[set]->GetMaximumBin();
    Double_t allbinmaxX = lower_lim+allbinmax*(upper_lim-lower_lim)/total_bins;
    Double_t allbinmaxY = hadct[set]->GetBinContent(allbinmax);
    allfitl = allbinmaxX - 2*observed_adct_sigma;
    allfith = allbinmaxX + 2*observed_adct_sigma;

    TF1 *gausfitall = new TF1("gausfitall",util::g_gfit_bg,allfitl,allfith,5);
    gausfitall->SetLineWidth(4);
    gausfitall->SetParameter(0,allbinmaxY);
    gausfitall->SetParameter(1,allbinmaxX);
    gausfitall->SetParLimits(1,allfitl,allfith);
    gausfitall->SetParameter(2,observed_adct_sigma);
    gausfitall->SetParLimits(2,1.,3*observed_adct_sigma);
    gausfitall->SetParameter(3,25);
    gausfitall->SetParameter(4,20);
    gausfitall->SetParLimits(3,0,200);
    
    hadct[set]->Fit("gausfitall","RBMQ");
    hadct[set]->Draw();

    Double_t allmean = gausfitall->GetParameter(1);
    allsig[set] = gausfitall->GetParameter(2);
    hadct[set]->SetTitle(Form("Set:%d QR:%d N:%d Mean:%f Sigma:%f",set,(Int_t)qreplay,adctN,allmean,allsig[set]));    
    //hadct[set]->Write();
    ADCt_all[set]->Write();

    for(Int_t c=0; c<hcal::maxHCalChan; c++){
      tcvalw[set][c] = 0.;
      tcerr[set][c] = 0.;

      Int_t half_chan = hcal::maxHCalChan/2;

      //Index through the canvas
      ADCt_top[set]->cd(c+1);
      if( c>=half_chan ){
	ADCt_bot[set]->cd(c+1-half_chan);
      }

      //Get slices from htDiff_ID and fit for mean vals
      Double_t tfitl = 0.;
      Double_t tfith = 0.;
      tcell[set][c] = c;
      tcellslice[set][c] = hap_hodocorr_ID[set]->ProjectionY(Form("tcellslice_%d_%d",set,c+1),c+1,c+1); 
      tcval[set][c] = adct_cal[set].old_param[c];

      Int_t sliceN = tcellslice[set][c]->GetEntries();
      Int_t binmax = tcellslice[set][c]->GetMaximumBin();
      Double_t binmaxX = lower_lim+binmax*(upper_lim-lower_lim)/total_bins;

      //TEMPORARY FIX
      if(binmaxX>1&&binmaxX<20&&tempfix)
	binmaxX=48;

      Double_t binmaxY = tcellslice[set][c]->GetBinContent(binmax);
      tfitl = binmaxX - Nsig*observed_adct_sigma;
      tfith = binmaxX + Nsig*observed_adct_sigma;

      if( sliceN<fit_event_min || binmaxX<=lower_lim || binmaxX>=upper_lim ){
	tcellslice[set][c]->Draw();
	tcellslice[set][c]->SetLineColor(kRed);
	continue;
      }

      //TF1 *gausfit = new TF1("gausfit",util::g_gfit_bg,tfitl,tfith,5);
      TF1 *gausfit = new TF1("gausfit","gaus",tfitl,tfith);
      gausfit->SetLineWidth(4);
      gausfit->SetLineColor(kGreen);
      gausfit->SetParameter(0,binmaxY);
      gausfit->SetParameter(1,binmaxX);
      gausfit->SetParLimits(1,tfitl,tfith);
      gausfit->SetParameter(2,observed_adct_sigma);
      gausfit->SetParLimits(2,1.,3*observed_adct_sigma);
      //gausfit->SetParameter(3,25);
      //gausfit->SetParLimits(3,0,200);

      tcellslice[set][c]->Fit("gausfit","RBMQ");
      tcellslice[set][c]->Draw();
      
      tcvalw[set][c] = gausfit->GetParameter(1);
      tcerr[set][c] = gausfit->GetParameter(2);

      //Catch bad fits on very small data sets      
      if( gausfit->GetParameter(1)>tfith || gausfit->GetParameter(1)<tfitl ){
	tcval[set][c] = binmaxX;
	gausfit->SetLineColor(kRed);
      }else{
	tcval[set][c] = gausfit->GetParameter(1);
      }
	
      //finish correction of channel well out of window
      if(c==corrChan&&channelCorrect)
	tcval[set][c] = gausfit->GetParameter(1) - tempChanOffset;

      //Write histogram
      tcellslice[set][c]->SetTitle(Form("Set:%d N:%d MaxX:%f Mean:%f Sigma:%f",set,sliceN,binmaxX,tcval[set][c],tcerr[set][c]));    
      tcellslice[set][c]->Write();

      tcval_avg[set] += tcval[set][c];
      tcval_Ng[set]++;

    }    
    ADCt_top[set]->Write();
    ADCt_bot[set]->Write();

    if( !qreplay ){
      std::string qplotpath_blkall = plotdir + Form("adct_tmax_set%d%s.png",set,wonkopt.c_str());
      std::string qplotpath_all = plotdir + Form("adct_tcut_set%d%s.png",set,wonkopt.c_str());
      std::string qplotpath_top = plotdir + Form("adct_fit_top_set%d%s.png",set,wonkopt.c_str());
      std::string qplotpath_bot = plotdir + Form("adct_fit_bot_set%d%s.png",set,wonkopt.c_str());
      ADCtblk_all[set]->SaveAs(qplotpath_blkall.c_str());
      //ADCt_all[set]->SaveAs(qplotpath_all.c_str());
      ADCt_top[set]->SaveAs(qplotpath_top.c_str());
      ADCt_bot[set]->SaveAs(qplotpath_bot.c_str());
    }else{
      std::string samp1plotpath = plotdir + Form("samp1_run_set%d%s.png",set,wonkopt.c_str());
      std::string samp2plotpath = plotdir + Form("samp2_run_set%d%s.png",set,wonkopt.c_str());
      samp1_run[set]->SaveAs(samp1plotpath.c_str());
      samp2_run[set]->SaveAs(samp2plotpath.c_str());
    } 	

    tcval_avg[set] /= tcval_Ng[set];
	
  }//endloop over nsets

  TCanvas *ADCt_graph[Ncal_set_size];
  TGraphErrors *gadct_c[Ncal_set_size];

  // Output text file for new ADCT offsets
  ofstream writeParFile;
  if( qreplay )
    new_adctoffset_path = "/dev/null"; //Safety to prevent overwriting constants on quasi-replay
  writeParFile.open( new_adctoffset_path );

  //Build tgraph and record adct offsets
  for( Int_t set=0; set<Ncal_set_size; set++ ){

    //Define calibration set by newest timestamp (config/tdcoffset/tdccalib) and record constant
    std::string offset_timestamp = adct_cal[set].timestamp;
    std::string config_timestamp = config_parameters.GetSBSTimestamp();
    std::string active_timestamp;
    util::tsCompare( offset_timestamp, config_timestamp, active_timestamp );
	  
    ADCt_graph[set] = new TCanvas(Form("ADCt_graph, timestamp: %s",active_timestamp.c_str()),"ADCt vs ID, means",1600,1200);
    ADCt_graph[set]->cd();

    //Make graphs with errors for reporting. All failed fits are zero here
    gadct_c[set] = new TGraphErrors( hcal::maxHCalChan, tcell[set], tcvalw[set], cellerr, tcerr[set] );
    gadct_c[set]->GetXaxis()->SetLimits(-10,290);  
    gadct_c[set]->GetYaxis()->SetLimits(lower_lim,upper_lim);
    gadct_c[set]->SetTitle(Form("ADCt_{hcal}-TDCmean_{hodo} vs Cell, timestamp: %s, qreplay: %d",active_timestamp.c_str(),(Int_t)qreplay ) );
    gadct_c[set]->GetXaxis()->SetTitle("Cell");
    gadct_c[set]->GetYaxis()->SetTitle("ADCt_{HCAL}-TDC_{MEAN,HODO}");
    gadct_c[set]->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
    gadct_c[set]->Draw();
    gadct_c[set]->Write(Form("gadct_c_s%d",set));

    // Output text file for new ADCT offsets
    if( qreplay ) continue; //Don't write if quasi-replay
    writeParFile << "#HCal ADCt offsets obtained " << date.c_str() << " for db adct offset timestamp " << offset_timestamp << endl;
    writeParFile << "#Offsets obtained from fits over ADCt distributions." << endl;
    writeParFile << active_timestamp << endl;
    writeParFile << db_adctoffset_variable << " =" << endl;

    Int_t cell = 0;

    for( Int_t r = 0; r<hcal::maxHCalRows; r++){
      for( Int_t c = 0; c<hcal::maxHCalCols; c++){
	if( tcval[set][cell]==0 || abs(tcerr[set][cell])<0.05 ){ //This will only work if most channels are closely aligned already and oldtdcoffsets aren't believable
	  writeParFile << tcval_avg[set] + adct_cal[set].old_param[cell] - ADCt_target << " ";
	}else{
	  writeParFile << tcval[set][cell] + adct_cal[set].old_param[cell] - ADCt_target << " ";
	}
	cell++;
      } // endloop over columns 
      writeParFile << endl;
    }  // endloop over rows 
    writeParFile << endl << endl;

    cell = 0;

    writeParFile << "(inverted) " << db_adctoffset_variable << " =" << endl;
    for( Int_t r = 0; r<hcal::maxHCalRows; r++){
      for( Int_t c = 0; c<hcal::maxHCalCols; c++){
	if( tcval[set][cell]==0 || abs(tcerr[set][cell])<0.05 ){ //This will only work if most channels are closely aligned already and oldtdcoffsets aren't believable
	  writeParFile << tcval_avg[set] - adct_cal[set].old_param[cell] - ADCt_target << " ";
	}else{
	  writeParFile << tcval[set][cell] - adct_cal[set].old_param[cell] - ADCt_target << " ";
	}
	cell++;
      } // endloop over columns 
      writeParFile << endl;
    }  // endloop over rows 
    writeParFile << endl << endl;
      
  } //endloop over calibration sets

  writeParFile.close();

  //Write tmax calibrations per set on qreplay
  ofstream write_tmax_file;
  if( !qreplay )
    new_tmax_path = "/dev/null";
  write_tmax_file.open( new_tmax_path );
  for( Int_t set=0; set<Ncal_set_size; set++ ){
    write_tmax_file << "#tmax parameters obtained " << date.c_str() << " for db adct offset timestamp " << offset_timestamp << endl;
    write_tmax_file << db_tmax_variable << " = " << tmax_factor*allsig[set] << endl;
  }
  write_tmax_file.close();

    //Add output report canvas
  TCanvas *c1[Ncal_set_size];
  Int_t report_height = 1800;

  for( Int_t s=0; s<Ncal_set_size; s++ ){

    c1[s] = new TCanvas(Form("adctalignreport_%d",s), Form("Configuration/Cut Information, Calibration Set %d",s), report_height, 900);
  
    // Set margin.
    c1[s]->SetLeftMargin(0.01);

    // Create a TText object.
    TText *t = new TText();

    // Set text align to the left (horizontal alignment = 1).
    t->SetTextAlign(11);
    t->SetTextSize(0.02);

    //make an array of strings
    std::string target_option = "All Available";
    if( h2only )
      target_option = "LH2";
    std::string report[linecount] = {
      "General ADCt Alignment Info",
      Form("Experiment: %s, Configuration: %d %s, Pass: %d", experiment, config, wonkopt.c_str(), pass),
      Form("Creation Date: %s", date.c_str() ),
      Form("Target(s) Used: %s", target_option.c_str() ),
      Form("Calibration Set: %s", adct_cal[s].timestamp.c_str() ),
      "",
      "Elastic Cuts",
      Form("Global Elastic Cuts: %s", adct_cal[s].gcut.c_str() ),
      "",
      "Other Cuts",
      Form("Minimum Ev per Cell : %d", adct_cal[s].minEv),
      "HCal Active Area (Projected Nucleon 1 row/col Within HCal Acceptance)"
    };
    // Loop to write the lines to the canvas.
    for( Int_t i = 0; i<linecount; i++ ) {
      // Vertical position adjusted according to line number.
      Double_t verticalPosition = 0.9 - i * 0.04;
      t->DrawTextNDC(0.1, verticalPosition, report[i].c_str());
    }

    c1[s]->Write();    
  }

  //clean up
  for( Int_t unused_set_index=Ncal_set_size; unused_set_index<hcal::gNstamp; unused_set_index++ ){
    delete hap_hodocorr_ID[unused_set_index];
  }

  fout->Write();
  st->Stop();

  cout << "Analysis and fits written to " << adctalign_path << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
