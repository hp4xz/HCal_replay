//SSeeds 9.1.23 Shortened qreplay designed to get adct and tdc vs run before new offsets and after offsets over two sets

#include <vector>
#include <iostream>
#include <algorithm>
#include <vector>
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

const Double_t minEventPerCell = 150; //minimum elastic events detected in cell to be considered for calibration
const Int_t total_bins = 500;
const Int_t lower_lim = -150;
const Int_t upper_lim = 150;
bool verbose = true;
const Double_t tFWHM = 4;

void overlayWithMeans(TH2D* hist, TCanvas* canvas, std::vector<Double_t> &means, std::vector<Double_t> &xcent) {
  if (!hist || !canvas) {
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

    //declare some dynamic fit variables
    Double_t FWHM = tFWHM;
    Int_t binMax = projY->GetMaximumBin();
    Double_t binCenter = projY->GetBinCenter( binMax );
    Double_t fitLowerLim = binCenter - FWHM;
    Double_t fitUpperLim = binCenter + FWHM;

    Double_t mean = -1000.;

    // Fit a Gaussian to this projection if sufficient stats, if no entries skip, else get arithmetic mean
    TF1 *gaussFit = new TF1("gausFit", "gaus");
    if( projY->GetEntries()==0 )
      continue;
    else if( projY->GetEntries() < minEventPerCell )
      mean = projY->GetMean(); //arithmetic mean on low statistics
    else{
      projY->Fit(gaussFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode
      mean = gaussFit->GetParameter(1);
    }

    //Double_t mean = projY->GetMean(); //Arithmetic mean here due to low statistics, visual aid only
    means.push_back(mean);

    // Determine the x-value as the center of the current bin
    Double_t xCenter = hist->GetXaxis()->GetBinCenter(binX);
    xcent.push_back(xCenter);

    int pointIdx = graph->GetN();
    graph->SetPoint(pointIdx, xCenter, mean);
    graph->SetPointError(pointIdx, 0, 0);
 
    delete projY; // Clean up
  }

  // Draw the original histogram, standard graph, and red graph
  hist->Draw("COLZ");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("P SAME");

  // Update the canvas
  canvas->Update();
}

//Main, configured for GMn pass0/1. tdcexclusive should be enabled to account for tdc difference from adct cal sets
void qreplay_timepatch( Int_t config = 11, Int_t pass = 1, bool nowonk = false, bool tdcexclusive = false ){
  
  const char *experiment = "gmn";

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  // Get the date
  string date = util::getDate();

  std::string tdce_word = "";
  if( tdcexclusive )
    tdce_word = "_tdcexclusive";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + Form("/hcal_calibrations/qreplay/timepatch_%s_sbs%d%s.root",experiment,config,tdce_word.c_str());

  //Set up path variables and output files
  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  // std::string db_path = gSystem->Getenv("DB_DIR");
  // std::string hcal_db_path = db_path + "/db_sbs.hcal.dat";

  std::string hcal_db_path = "/w/halla-scshelf2102/sbs/seeds/alt_sbsreplay/SBS-replay/DB/db_sbs.hcal.dat";

  std::string new_adct_wonkpath = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d_wonkset.txt",experiment,config,pass);
  std::string new_tdcoffset_wonkpath = Form("../timing/parameters/tdcoffsets_class_%s_conf%d_pass%d_wonkset.txt",experiment,config,pass);

  std::string new_adct_normalpath = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d_normalset.txt",experiment,config,pass);
  std::string new_tdcoffset_normalpath = Form("../timing/parameters/tdcoffsets_class_%s_conf%d_pass%d_normalset.txt",experiment,config,pass);

  std::string new_adct_nowonkpath = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string new_tdcoffset_nowonkpath = Form("../timing/parameters/tdcoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);

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


  //Histograms
  TH2D *hadct_run_before = new TH2D("hadct_run_before",
				   "ADCt-HODOtmean vs Run, Before Alignment",
				   (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				   runs[0].runnum-1,
				   runs[nruns-1].runnum+1,
				   total_bins,
				   lower_lim,
				   upper_lim);
  
  TH2D *hadct_run_after = new TH2D("hadct_run_after",
				   "ADCt-HODOtmean vs Run, After Alignment",
				   (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				   runs[0].runnum-1,
				   runs[nruns-1].runnum+1,
				   total_bins,
				   lower_lim,
				   upper_lim);

  TH2D *htdc_run_before = new TH2D("htdc_run_before",
				   "TDC-HODOtmean vs Run, Before Alignment",
				   (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				   runs[0].runnum-1,
				   runs[nruns-1].runnum+1,
				   total_bins,
				   lower_lim,
				   upper_lim);
  
  TH2D *htdc_run_after = new TH2D("htdc_run_after",
				   "TDC-HODOtmean vs Run, After Alignment",
				   (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				   runs[0].runnum-1,
				   runs[nruns-1].runnum+1,
				   total_bins,
				   lower_lim,
				   upper_lim);

  TH1D *hadct_allchan_allruns_before = new TH1D("hadct_allchan_allruns_before",
					"ADCt-HODOtmean All Runs All Channels, Before Alignment",
					total_bins,
					lower_lim,
					upper_lim);

  TH1D *hadct_allchan_allruns_after = new TH1D("hadct_allchan_allruns_after",
					"ADCt-HODOtmean All Runs All Channels, After Alignment",
					total_bins,
					lower_lim,
					upper_lim);

  TH1D *htdc_allchan_allruns_before = new TH1D("htdc_allchan_allruns_before",
					"TDC-HODOtmean All Runs All Channels, Before Alignment",
					total_bins,
					lower_lim,
					upper_lim);

  TH1D *htdc_allchan_allruns_after = new TH1D("htdc_allchan_allruns_after",
					"TDC-HODOtmean All Runs All Channels, After Alignment",
					total_bins,
					lower_lim,
					upper_lim);


  //reporting indices
  Int_t target_change_index;
  Double_t config_sampling_fraction;
  Double_t config_e_sigma_ratio;
  std::string ts_compare = "";
  bool cut_first = true;
  
  //loop over sets
  for (Int_t s=0; s<2; s++){

    //if configured to consider the configuration as a single set, skip loop over added sets
    if(nowonk && s==1)
      continue;

    bool wonkset = true;
    if(s==1)
      wonkset = false;

    bool first = true;

    //loop over runs
    for (Int_t r=0; r<nruns; r++) {

      bool ts_different = false;
      
      //Get run experimental parameters
      std::string current_timestamp = runs[r].adcg_ts; //For now, only select on energy calibration sets
      Int_t current_runnumber = runs[r].runnum;
      std::string current_target = runs[r].target;
      
      //Hell of nested ifs incoming
      //adct opt
      if( !nowonk ){
	if( !tdcexclusive ){
	  //for pass2 to account for timing shift
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
	}else{
	  //tdc opt
	  if( config==11 ){
	    if( wonkset ){
	      if( current_runnumber<12995 )
		continue;
	    }else{
	      if( current_runnumber>12995 )
		//cout << current_runnumber << endl;
		continue;
	    }
	  }else if( config==9 ){
	    if( wonkset ){
	      if( current_runnumber!=13682 )
		continue;
	    }else{
	      if( current_runnumber==13682 )
		continue;
	    }
	  }
	}
      }	

      if(verbose) 
	std::cout << "Current run number: " << current_runnumber << ", current timestamp: " << current_timestamp << ", current target: " << current_target << std::endl;

      if( !(current_timestamp.compare(ts_compare)==0) ){
	ts_compare = current_timestamp;
	ts_different = true;
      }

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

      //Record old gain and offset parameters by tstamp from database (assumes one file for all timestamps)
      Double_t old_adct_offsets[hcal::maxHCalChan] = {0.};
      Double_t old_tdc_offsets[hcal::maxHCalChan] = {0.};
      Double_t old_tdc_calib = 1;

      if(verbose) std::cout << "Reading parameters from file: " << hcal_db_path << std::endl;
      util::readDB( hcal_db_path, runs[r].adct_ts, db_adct_variable, old_adct_offsets );
      util::readDB( hcal_db_path, runs[r].tdc_ts, db_tdcoffset_variable, old_tdc_offsets );
      util::readDB( hcal_db_path, runs[r].tdcc_ts, db_tdccalib_variable, old_tdc_calib );

      //get new adct offsets from calibrated parameters
      Double_t new_adct_offsets[hcal::maxHCalChan] = {0.};
      Double_t new_tdc_offsets[hcal::maxHCalChan] = {0.};

      //if no selected timestamp is passed, we use the first set of new coefficients that appears
      string quality_timestamp = "";
      if( quality_timestamp.compare("")==0 ) quality_timestamp = "none";

      //select between the two offset sets on both sbs11 and 14 for adct quality plots and sbs11/9 for tdc
      if( !tdcexclusive ){
	if(wonkset&&!nowonk){
	  if(verbose) std::cout << "Reading parameters from file: " << new_adct_wonkpath << std::endl;
	  util::readDB( new_adct_wonkpath, "none", db_adct_variable, new_adct_offsets );
	  if(verbose) std::cout << "Reading parameters from file: " << new_tdcoffset_wonkpath << std::endl;
	  util::readDB( new_tdcoffset_wonkpath, "none", db_tdcoffset_variable, new_tdc_offsets );
	}else if(!wonkset&&!nowonk){
	  if(verbose) std::cout << "Reading parameters from file: " << new_adct_normalpath << std::endl;
	  util::readDB( new_adct_normalpath, "none", db_adct_variable, new_adct_offsets );
	  if(verbose) std::cout << "Reading parameters from file: " << new_tdcoffset_normalpath << std::endl;
	  util::readDB( new_tdcoffset_normalpath, "none", db_tdcoffset_variable, new_tdc_offsets );
	}else if(nowonk){
	  if(verbose) std::cout << "Reading parameters from file: " << new_adct_nowonkpath << std::endl;
	  util::readDB( new_adct_nowonkpath, "none", db_adct_variable, new_adct_offsets );
	  if(verbose) std::cout << "Reading parameters from file: " << new_tdcoffset_nowonkpath << std::endl;
	  util::readDB( new_tdcoffset_nowonkpath, "none", db_tdcoffset_variable, new_tdc_offsets );

	}else{
	  cout << "ERROR: configured to read no database parameters. Option should be logically excluded. Check script." << endl;
	  return;
	}
      }else{
	if(wonkset&&!nowonk){
	  if(verbose) std::cout << "Reading parameters from file: " << new_adct_nowonkpath << std::endl;
	  util::readDB( new_adct_nowonkpath, "none", db_adct_variable, new_adct_offsets );
	  if(verbose) std::cout << "Reading parameters from file: " << new_tdcoffset_wonkpath << std::endl;
	  util::readDB( new_tdcoffset_wonkpath, "none", db_tdcoffset_variable, new_tdc_offsets );
	}else if(!wonkset&&!nowonk){
	  if(verbose) std::cout << "Reading parameters from file: " << new_adct_nowonkpath << std::endl;
	  util::readDB( new_adct_nowonkpath, "none", db_adct_variable, new_adct_offsets );
	  if(verbose) std::cout << "Reading parameters from file: " << new_tdcoffset_normalpath << std::endl;
	  util::readDB( new_tdcoffset_normalpath, "none", db_tdcoffset_variable, new_tdc_offsets );
	}else if(nowonk){
	  if(verbose) std::cout << "Reading parameters from file: " << new_adct_nowonkpath << std::endl;
	  util::readDB( new_adct_nowonkpath, "none", db_adct_variable, new_adct_offsets );
	  if(verbose) std::cout << "Reading parameters from file: " << new_tdcoffset_nowonkpath << std::endl;
	  util::readDB( new_tdcoffset_nowonkpath, "none", db_tdcoffset_variable, new_tdc_offsets );

	}else{
	  cout << "ERROR: configured to read no database parameters. Option should be logically excluded. Check script." << endl;
	  return;
	}
      }


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
	    std::cout << new_adct_offsets[i] << " ";
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
	    std::cout << new_tdc_offsets[i] << " ";
	  }
	  std::cout << std::endl;
	}
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

      std::string wonkyness = std::to_string(wonkset);
      if(nowonk)
	wonkyness = "nowonk";

      //Main loop over events in run
      while (C->GetEntry(nevent++)) {

	cout << "Wonkset = " << wonkyness << ". Run index " << r << "/" << nruns << ". Analyzing " << current_target << " run " << current_runnumber << ": " <<  nevent << "/" << nevents << " \r";
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
	Double_t hcalatime_bc = hcalatime-HODOtmean;
	Double_t hcalatime_c = hcalatime + old_adct_offsets[(Int_t)pblkid] - new_adct_offsets[(Int_t)pblkid];
	Double_t hcalatime_fc = hcalatime_c-HODOtmean;

	//Assuming the tdc calib hasn't changed, which is true on SBS14 and SBS11 from pass0/1 to pass2
	Double_t hcaltime = cblktime[0];
	Double_t hcaltime_bc = hcaltime-HODOtmean;
	Double_t hcaltime_c = hcaltime +
		       old_tdc_offsets[(Int_t)pblkid] * old_tdc_calib - 
		       new_tdc_offsets[(Int_t)pblkid] * old_tdc_calib;
	Double_t hcaltime_fc = hcaltime_c-HODOtmean;

	//Double_t adct_tc = hcalatime-BBps_atime; //Primary cluster, primary block tdc with trigger correction (ns)
	Double_t adct_tc = hcalatime-HODOtmean; //Primary cluster, primary block tdc with trigger correction (ns)

	hadct_run_before->Fill(current_runnumber,hcalatime_bc);
	hadct_run_after->Fill(current_runnumber,hcalatime_fc);
	hadct_allchan_allruns_before->Fill(hcalatime_bc);
	hadct_allchan_allruns_after->Fill(hcalatime_fc);

	htdc_run_before->Fill(current_runnumber,hcaltime_bc);
	htdc_run_after->Fill(current_runnumber,hcaltime_fc);
	htdc_allchan_allruns_before->Fill(hcaltime_bc);
	htdc_allchan_allruns_after->Fill(hcaltime_fc);

      }//end loop over event

      // getting ready for the next run
      C->Reset();

    }//endloop over runs

  }//endloop over sets

  //overlay adct vs run histo before new coeff
  TCanvas *c1 = new TCanvas("c1","ADCt vs Run, Before",1600,800);
  std::vector<Double_t> c1y;
  std::vector<Double_t> c1x;
  //gStyle->SetOptStat(0);
  //c1->cd();
  overlayWithMeans(hadct_run_before,c1,c1y,c1x);
 
  c1->SetGridy();
  c1->Write();

  //overlay adct vs run histo after new coeff
  TCanvas *c2 = new TCanvas("c2","ADCt vs Run, After",1600,800);
  vector<Double_t> c2y;
  vector<Double_t> c2x;
  //gStyle->SetOptStat(0);
  //c2->cd();
  overlayWithMeans(hadct_run_after,c2,c2y,c2x);
 
  c2->SetGridy();
  c2->Write();

  //overlay adct vs run histo before new coeff
  TCanvas *c3 = new TCanvas("c3","TDC vs Run, Before",1600,800);
  vector<Double_t> c3y;
  vector<Double_t> c3x;
  //gStyle->SetOptStat(0);
  //c3->cd();
  overlayWithMeans(htdc_run_before,c3,c3y,c3x);


  c3->SetGridy();
  c3->Write();

  //overlay adct vs run histo before new coeff
  TCanvas *c4 = new TCanvas("c4","TDC vs Run, After",1600,800);
  vector<Double_t> c4y;
  vector<Double_t> c4x;
  //gStyle->SetOptStat(0);
  //c4->cd();
  overlayWithMeans(htdc_run_after,c4,c4y,c4x);

  c4->SetGridy();
  c4->Write();

  TCanvas *c5 = new TCanvas("c5","ADCt amean before/after",1600,800);
  //c5->cd();
  TGraph *g5bef = new TGraph();
  for( int i=0; i<c1y.size(); i++ ){
    int pIdx = g5bef->GetN();
    g5bef->SetPoint(pIdx, c1x[i], c1y[i]);
  }
  TGraph *g5af = new TGraph();
  for( int i=0; i<c2y.size(); i++ ){
    c2x[i]+=0.1;
    int pIdx = g5af->GetN();
    g5af->SetPoint(pIdx, c2x[i], c2y[i]);
  }
  g5bef->SetMarkerStyle(20);
  g5bef->SetMarkerColor(kBlack);
  g5bef->SetLineColor(kRed);
  g5bef->GetYaxis()->SetRangeUser(-80,80);
  g5bef->Draw("P");

  g5af->SetMarkerStyle(20);
  g5af->SetMarkerColor(kBlack);
  g5af->SetLineColor(kBlack);
  g5af->Draw("P SAME");

  c5->SetGridy();
  c5->Write();


  TCanvas *c6 = new TCanvas("c6","TDC amean before/after",1600,800);
  c6->cd();
  TGraph *g6bef = new TGraph();
  for( int i=0; i<c3y.size(); i++ ){
    int pIdx = g6bef->GetN();
    g6bef->SetPoint(pIdx, c3x[i], c3y[i]);
  }
  TGraph *g6af = new TGraph();
  for( int i=0; i<c4y.size(); i++ ){
    c4x[i]+=0.1;
    int pIdx = g6af->GetN();
    g6af->SetPoint(pIdx, c4x[i], c4y[i]);
  }
  g6bef->SetMarkerStyle(20);
  g6bef->SetMarkerColor(kBlack);
  g6bef->SetMarkerColor(kRed);
  g6bef->GetYaxis()->SetRangeUser(-80,80);
  g6bef->Draw("P");

  g6af->SetMarkerStyle(20);
  g6af->SetMarkerColor(kBlack);
  g6af->SetMarkerColor(kBlack);
  g6af->Draw("P SAME");

  c6->SetGridy();
  c6->Write();

  fout->Write();

  st->Stop();

  cout << endl << "Timepatch complete. Output file written to " << fout_path << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    

}
