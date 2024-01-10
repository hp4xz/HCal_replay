//SSeeds 9.1.23 Shortened qreplay for timing only where set ranges and plots are configured by user. Intended for use with output offset parameters via adct_setalign.C and tdc_setalign.C

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

const double minEventPerCell = 300; //minimum elastic events detected in cell to be considered for calibration
const int total_bins = 500;
const int lower_lim = -150;
const int upper_lim = 150;
bool verbose = true;
const double tFWHM = 4;

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

//extracts the run ranges and exceptions from the name of the patch files
std::vector<int> extractLastFourIntegers(const std::string& s) {
  std::vector<int> allIntegers, lastFourIntegers;
  int num = 0;
  bool foundDigit = false;

  for (char ch : s) {
    if (isdigit(ch)) {
      num = num * 10 + (ch - '0');
      foundDigit = true;
    } else if (foundDigit) {
      allIntegers.push_back(num);
      num = 0;
      foundDigit = false;
    }
  }

  // Handle the case where the string ends with a number
  if (foundDigit) {
    allIntegers.push_back(num);
  }

  // Copy the last four numbers if there are at least four numbers
  if (allIntegers.size() >= 4) {
    lastFourIntegers.insert(lastFourIntegers.end(), allIntegers.end() - 4, allIntegers.end());
  } else {
    lastFourIntegers = allIntegers;
  }

  return lastFourIntegers;
}

void overlayWithMeans(TH2D* hist, TCanvas* canvas, std::vector<double> &means, std::vector<double> &xcent) {
  if (!hist || !canvas) {
    std::cerr << "Null histogram or canvas pointer!" << std::endl;
    return;
  }

  //Make histogram to store all tdc slices
  TH1D *hist_all = hist->ProjectionY("hist_all", 1, hist->GetNbinsX(), "e");
  //Get fit parameters
  vector<Double_t> fit_vec = util::fitGaussianAndGetFineParams(hist_all, tFWHM);
  double win_low = fit_vec[1] - 5*fit_vec[2];
  double win_high = fit_vec[1] + 5*fit_vec[2];

  // Use the provided canvas
  canvas->cd();

  // Create TGraphErrors for standard and red points
  TGraphErrors* graph = new TGraphErrors();

  for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
    // Project the 2D histogram onto a 1D histogram along the y-axis for the current x-bin
    TH1D* projY = hist->ProjectionY("_py", binX, binX);

    //declare some dynamic fit variables
    double FWHM = tFWHM;
    int binMax = projY->GetMaximumBin();
    double binCenter = projY->GetBinCenter( binMax );
    double fitLowerLim = binCenter - FWHM;
    double fitUpperLim = binCenter + FWHM;

    double mean = -1000.;

    // Fit a Gaussian to this projection if sufficient stats, if no entries skip, else get arithmetic mean
    TF1 *gaussFit = new TF1("gausFit", "gaus");

    // Calculate the integral (total number of entries) in the range
    double sliceN = projY->Integral(1, projY->GetNbinsX());
    if( sliceN<minEventPerCell ) //continue if too sparse to get distribution
      continue;
    else{
      projY->Fit(gaussFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode
      mean = gaussFit->GetParameter(1);
    }

    //double mean = projY->GetMean(); //Arithmetic mean here due to low statistics, visual aid only
    means.push_back(mean);

    // Determine the x-value as the center of the current bin
    double xCenter = hist->GetXaxis()->GetBinCenter(binX);
    xcent.push_back(xCenter);

    int pointIdx = graph->GetN();
    graph->SetPoint(pointIdx, xCenter, mean);
    graph->SetPointError(pointIdx, 0, 0);
 
    delete projY; // Clean up
  }

  // Draw the original histogram, standard graph, and red graph
  hist->GetYaxis()->SetRangeUser(win_low,win_high);
  hist->Draw("COLZ");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("P SAME");

  // Update the canvas
  canvas->Update();
}

//Main. Currently configured to take a maximum of four timing offset sets. leave as none if not needed, build patch list from first to fourth in run order. pass absolute paths
void qreplay_settimepatch( std::string experiment = "gen", 
			   int config = 2, 
			   int pass = 1, 
			   bool best_clus = false,
			   std::string type = "adct", 
			   std::string first_patch_path = "none",
			   std::string second_patch_path = "none",
			   std::string third_patch_path = "none",
			   std::string fourth_patch_path = "none" ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  // Get the date
  string date = util::getDate();
  
  // consolidate pass0/1
  if(pass==0)
    pass=1;

  //added for ease of use
  if( second_patch_path.compare("0")==0 )
    second_patch_path = "none";
  if( third_patch_path.compare("0")==0 )
    third_patch_path = "none";
  if( fourth_patch_path.compare("0")==0 )
    fourth_patch_path = "none";

  //try to catch user error
  if( first_patch_path.compare("none")==0 ||
      ( second_patch_path.compare("none")==0 && 
	( third_patch_path.compare("none")!=0 || 
	  fourth_patch_path.compare("none")!=0 ) ) ||
      ( third_patch_path.compare("none")==0 && 
	fourth_patch_path.compare("none")!=0 ) ){
    std::cout << "ERROR: Reorder patch paths filling from first to fourth and leaving unused as none, then retry." << std::endl;
    return;
  }

  //Set up path variables and output files
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + Form("/hcal_calibrations/qreplay/%s_settimepatch_fix_%s_sbs%d.root",type.c_str(),experiment.c_str(),config);
  //std::string fout_path = "test.root";
  std::string plotdir = Form("../quality_plots/%s/conf%d_settime/%s/",experiment.c_str(),config,type.c_str());

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );
  
  //Paths to database. These should point to database files which produced the data analyzed.
  std::string db_path;
  std::string exper = experiment;
  if( exper.compare("gmn")==0 ){
    if( pass==0 || pass==1 )
      db_path = hcal::gmn_p1_db_path;
    else if( pass==2 )
      db_path = hcal::gmn_p2_db_path;
    else{
      cout << "ERROR: Enter a valid pass for gmn." << endl;
      return;
    }
  }else if( exper.compare("gen")==0 ){
    db_path = hcal::gmn_p2_db_path;
  }else{
    cout << "ERROR: Enter a valid experiment (gmn or gen)" << endl;
    return;
  }

  std::string hcal_db_path = db_path + "/db_sbs.hcal.dat";

  //Get run ranges from paths
  std::vector<std::vector<int>> patch_ranges;
  std::vector<std::string> patches;

  patches.push_back( first_patch_path );
  patch_ranges.push_back( extractLastFourIntegers(first_patch_path) );

  //Get additional run paths and ranges as necessary
  if( second_patch_path.compare("none") != 0 ){
    patches.push_back( second_patch_path );
    patch_ranges.push_back( extractLastFourIntegers(patches[1]) );
  }
  if( third_patch_path.compare("none") != 0 ){
    patches.push_back( third_patch_path );
    patch_ranges.push_back( extractLastFourIntegers(patches[2]) );
  }
  if( fourth_patch_path.compare("none") != 0 ){
    patches.push_back( fourth_patch_path );
    patch_ranges.push_back( extractLastFourIntegers(patches[3]) );
  }

  //set the type for this run
  bool adct = type.compare("adct")==0 || type.compare("ADCt")==0 || type.compare("ADCt")==0;
  bool tdc = type.compare("tdc")==0 || type.compare("TDC")==0;
  if( adct == tdc ){
    std::cout << "ERROR: provide either adct or tdc for type argument." << std::endl;
    return;
  }

  //get database variable name
  std::string db_variable;
  if(adct)
    db_variable = "sbs.hcal.adc.timeoffset";
  else
    db_variable = "sbs.hcal.tdc.offset";
  std::string db_tdccalib_variable = "sbs.hcal.tdc.calib";
  
  std::string db_inv_variable;
  if(adct)
    db_inv_variable = "(inverted) sbs.hcal.adc.timeoffset";
  else
    db_inv_variable = "sbs.hcal.tdc.offset";

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s_p%d/",experiment.c_str(),pass); //unique to my environment for now
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
  TH2D *ht_run_before = new TH2D("ht_run_before",
				 Form( "%s-HODOtmean vs Run, Before Alignment",type.c_str() ),
				 (runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				 runs[0].runnum-1,
				 runs[nruns-1].runnum+1,
				 total_bins,
				 lower_lim,
				 upper_lim);
  
  TH2D *ht_run_after = new TH2D("ht_run_after",
				Form( "%s-HODOtmean vs Run, After Alignment",type.c_str() ),
				(runs[nruns-1].runnum+1) - (runs[0].runnum-1),
				runs[0].runnum-1,
				runs[nruns-1].runnum+1,
				total_bins,
				lower_lim,
				upper_lim);

  TH2D *ht_id_before = new TH2D("ht_id_before",
				Form( "%s-HODOtmean vs ID, Before Alignment",type.c_str() ),
				hcal::maxHCalChan,
				0,
				hcal::maxHCalChan,
				total_bins,
				lower_lim,
				upper_lim);
  
  TH2D *ht_id_after = new TH2D("ht_id_after",
			       Form( "%s-HODOtmean vs ID, After Alignment",type.c_str() ),
			       hcal::maxHCalChan,
			       0,
			       hcal::maxHCalChan,
			       total_bins,
			       lower_lim,
			       upper_lim);

  TH2D *ht_row_before = new TH2D("ht_row_before",
				Form( "%s-HODOtmean vs row, Before Alignment",type.c_str() ),
				hcal::maxHCalRows,
				0,
				hcal::maxHCalRows,
				total_bins,
				lower_lim,
				upper_lim);
  
  TH2D *ht_row_after = new TH2D("ht_row_after",
			       Form( "%s-HODOtmean vs row, After Alignment",type.c_str() ),
			       hcal::maxHCalRows,
			       0,
			       hcal::maxHCalRows,
			       total_bins,
			       lower_lim,
			       upper_lim);

  TH1D *ht_allchan_allruns_before = new TH1D("ht_allchan_allruns_before",
					     Form( "%s-HODOtmean All Runs All Channels, Before Alignment", type.c_str() ),
					     total_bins,
					     lower_lim,
					     upper_lim);

  TH1D *ht_allchan_allruns_after = new TH1D("hadct_allchan_allruns_after",
					    Form( "%s-HODOtmean All Runs All Channels, After Alignment", type.c_str() ),
					    total_bins,
					    lower_lim,
					    upper_lim);


  //reporting/run indices
  int target_change_index;
  double config_sampling_fraction;
  double config_e_sigma_ratio;
  std::string ts_compare = "";
  bool cut_first = true;
  vector<int> used_runs;

  //loop over patches
  for (size_t s=0; s<patch_ranges.size(); s++){

    std::string active_patch_path = patches[s];

    std::filesystem::path patch_file(active_patch_path);
    
    std::cout << std::endl << "Switching to patch " << patch_file.filename() << std::endl;

    bool first = true;

    //loop over runs and exclude. should update util::ReadRunList to improve efficiency in the future
    for (int r=0; r<nruns; r++) {
      int current_runnumber = runs[r].runnum;

      //first deal with exclusions
      int patch_starting_run = patch_ranges[s][0];
      int patch_ending_run = patch_ranges[s][1];
      int patch_starting_exclude_run = patch_ranges[s][2];
      int patch_ending_exclude_run = patch_ranges[s][3];

      //zero in these spots indicate no ranges/exceptions are to be made
      bool runall = patch_starting_run == 0;
      bool noexceptions = patch_starting_exclude_run == 0;

      if( !runall &&
	  (current_runnumber < patch_starting_run ||
	   current_runnumber > patch_ending_run) )
	continue;

      if( !noexceptions && 
	  (current_runnumber >= patch_starting_exclude_run &&
	   current_runnumber <= patch_ending_exclude_run) )
	continue;

      //error out where the patches overlap
      for( size_t r=0; r<used_runs.size(); ++r ){
	if( current_runnumber == used_runs[r] ){
	  std::cout << "ERROR: patch ranges overlap. Check input file paths and logic." << std::endl;
	  return 0;
	}
	used_runs.push_back(current_runnumber);
      }

      bool ts_different = false;

      //Get run experimental parameters
      std::string current_timestamp = runs[r].adcg_ts; //For now, only select on energy calibration sets
      std::string current_target = runs[r].target;
    
      if(verbose) 
	std::cout << "Current run number: " << current_runnumber << ", current timestamp: " << current_timestamp << ", current target: " << current_target << std::endl;

      if( !(current_timestamp.compare(ts_compare)==0) ){
	ts_compare = current_timestamp;
	ts_different = true;
      }

      std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
      int mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A

      //Get run paths
      // std::string rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-%s/pass%d/SBS%d/%s/rootfiles/",experiment.c_str(),pass,config,targ_uppercase.c_str());
      // if( experiment.compare("gmn")==0 && pass==2 )
      // 	rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2/QA_FINAL/SBS%d/%s/rootfiles/",config,targ_uppercase.c_str());

      std::string rootfile_dir = runs[r].rootfile_dir;
      std::string rootfile_path = rootfile_dir + Form("/*%d*",current_runnumber);

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
      double old_tdc_calib = 1;

      //switch to correct timestamp where necessary
      std::string time_timestamp;
      if(adct)
	time_timestamp = runs[r].adct_ts;
      else
	time_timestamp = runs[r].tdc_ts;

      if(verbose) std::cout << "Reading parameters from file " << hcal_db_path << ":" << std::endl;
      util::readDB( hcal_db_path, time_timestamp, db_variable, old_time_offsets );
      //util::readDB( hcal_db_path, runs[r].tdcc_ts, db_tdccalib_variable, old_tdc_calib );
      util::readDB( hcal_db_path, time_timestamp, db_tdccalib_variable, old_tdc_calib );

      //get new adct offsets from calibrated parameters
      double new_time_offsets[hcal::maxHCalChan] = {0.};

      //if no selected timestamp is passed, we use the first set of new coefficients that appears
      string quality_timestamp = "";
      if( quality_timestamp.compare("")==0 ) quality_timestamp = "none";

      //read in the new offsets. assumes first set in file at active_patch_path is the one desired
      if(verbose) std::cout << "Reading parameters from file: " << active_patch_path << std::endl;
      util::readDB( active_patch_path, "none", db_inv_variable, new_time_offsets );

      //Report offset parameters on first loop
      if( verbose && first ){
	std::cout << Form( "Loaded set of old %s offset coefficients:", type.c_str() ) << std::endl;
	sendOffsetsToConsole( old_time_offsets, db_variable );

	std::cout << Form( "Loaded set of new %s offset coefficients:", type.c_str() ) << std::endl;
	sendOffsetsToConsole( new_time_offsets, db_inv_variable );

	first = false;
      }

      //Get available cuts for current config/target/field combination. Set index = 0 as no calibration is being done here and physics cuts do not vary by calibration set
      //vector<calcut> cut;
      vector<caldiag> cut;
      //util::ReadCutList(struct_dir,experiment,config,0,pass,current_target,mag,verb,cut);
      util::ReadDiagnosticCutList(struct_dir,experiment,config,current_target,mag,verb,cut);
      if( cut_first ){
	std::cout << cut[0];
	cut_first = false;
      }

      std::string gcut = cut[0].gcut;
      Double_t W2min = cut[0].W2_min;
      Double_t W2max = cut[0].W2_max;
      Double_t dx0_p = cut[0].dx0_p;
      Double_t dx0_n = cut[0].dx0_n;
      Double_t dy0 = cut[0].dy0;
      Double_t dxsig_p = cut[0].dx_sig_p;
      Double_t dxsig_n = cut[0].dx_sig_n;
      Double_t dysig = cut[0].dy_sig;
      Double_t atime0 = cut[0].atime0;
      Double_t atimesig = cut[0].atime_sig;
      Int_t min_ev = cut[0].min_ev;

      // Setting up chain and branch addresses
      C = new TChain("T");
      
      if(verbose)
	cout << "Loading rootfile from path: " << rootfile_path << endl;
      
      C->Add(rootfile_path.c_str());

      C->SetBranchStatus("*",0);    
      double BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
      double BBtr_vz[hcal::maxTracks];
      double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, BBsh_atime, BBsh_nclus;	
      double HCALx, HCALy, HCALe;
      double pblkrow, pblkcol, nblk, nclus;
      int Ncid;
      double cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan], cblkagain[hcal::maxHCalChan];
      double cid[hcal::maxHCalClus], crow[hcal::maxHCalClus], ccol[hcal::maxHCalClus], ce[hcal::maxHCalClus], cx[hcal::maxHCalClus], cy[hcal::maxHCalClus], catime[hcal::maxHCalClus], ctdctime[hcal::maxHCalClus];
      double hodotmean[hcal::maxHCalClus];

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
      C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
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
      C->SetBranchStatus( "bb.sh.atimeblk", 1 );
      C->SetBranchStatus( "bb.sh.nclus", 1 );
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
      C->SetBranchAddress( "sbs.hcal.clus.tdctime", ctdctime );
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
      C->SetBranchAddress( "bb.sh.atimeblk", &BBsh_atime ); 
      C->SetBranchAddress( "bb.sh.nclus", &BBsh_nclus ); 
      C->SetBranchAddress( "bb.hodotdc.clus.tmean", hodotmean );
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
      double BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
      double Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;

      long nevent = 0, nevents = C->GetEntries(); 
      int treenum = 0, currenttreenum = 0;

      //Main loop over events in run
      while (C->GetEntry(nevent++)) {

	cout << "Patch set index " << s << "/" << patch_ranges.size() << " with potentially many runs. Analyzing " << current_target << " run " << current_runnumber << ": " <<  nevent << "/" << nevents << " \r";
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

	//e' p reconstruction with track angles (not momentum)
	double nu = pbeam.E() - pcent;
	double pNexp = sqrt( pow(nu, 2.) + 2. * M_avg * nu );
	double thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	double Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
	double W2 = pow( M_avg, 2.0 ) + 2.0*M_avg * (ebeam_c-pe.E()) - Q2;

	/////////////////////
	//W2 elastic cut
	bool failedW2 = W2>W2max || W2<W2min;
	if(failedW2)
	  continue;

	//Calculate h-arm quantities
	vector<double> xyhcalexp; util::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
	double dx = HCALx - xyhcalexp[0];
	double dy = HCALy - xyhcalexp[1];

	
	//Best Cluster analysis using intime algorithm for simplicity on diagnosis

	//Set up best cluster index and energy for sorting
	Double_t best_cluster[2] = {-1,0.};

	//loop through all clusters and select without HCal position information
	for( int c=0; c<Ncid; c++ ){
	
	  //calculate h-arm physics quantities per cluster
	  double atime = catime[c];
	  double atime_hodo = atime - hodotmean[0];
	  double cluse = ce[c];
	
	  //wide cut around HCal/Hodo time
	  bool passedTime = abs(atime_hodo-atime0)<3*atimesig;

	  if( passedTime && cluse>best_cluster[1] ){
	    best_cluster[0] = c;
	    best_cluster[1] = cluse;
	  }
	}//endloop over cluster elements

	//If no cluster passes atime cut, use default
	if(best_cluster[0]==-1)
	  best_cluster[0]=0;

	//Switch between best clusters for systematic analysis
	Int_t cidx_best = best_cluster[0];
      
	//Calculations from the best cluster
	Double_t dx_bestcluster = cx[cidx_best] - xyhcalexp[0];
	Double_t dy_bestcluster = cy[cidx_best] - xyhcalexp[1];
	Double_t hatime_bestcluster = catime[cidx_best];
	Double_t htime_bestcluster = ctdctime[cidx_best];
	Double_t hcoin_bestcluster = catime[cidx_best] - BBsh_atime;
	Double_t ce_bestcluster = ce[cidx_best];
	Double_t x_bestcluster = cx[cidx_best];
	Double_t y_bestcluster = cy[cidx_best];
	Double_t row_bestcluster = crow[cidx_best];
	Double_t col_bestcluster = ccol[cidx_best];
	Double_t id_bestcluster = cid[cidx_best];

	//Cut on dxdy spots bools
	bool is_proton = util::Nspotcheck(dy, dx, dy0, dx0_p, 3*dysig, 3*dxsig_p);
	bool is_neutron = util::Nspotcheck(dy, dx, dy0, dx0_n, 3*dysig, 3*dxsig_n);
	bool is_proton_bc = util::Nspotcheck(dy_bestcluster, dx_bestcluster, dy0, dx0_p, 3*dysig, 3*dxsig_p);
	bool is_neutron_bc = util::Nspotcheck(dy_bestcluster, dx_bestcluster, dy0, dx0_n, 3*dysig, 3*dxsig_n);

	double pblkid;
	double pblkrow;
	double hcalatime;
	double hcaltime;
	bool pass_spot = false;
	if(best_clus){
	  if(is_proton_bc || is_neutron_bc)
	    pass_spot = true;
	  pblkid = (double)cid[cidx_best];
	  pblkrow = (double)crow[cidx_best];
	  hcalatime = catime[cidx_best];
	  hcaltime = ctdctime[cidx_best];
	}else{
	  if(is_proton || is_neutron)
	    pass_spot = true;
	  pblkid = (double)cid[0];
	  pblkrow = (double)crow[0];
	  hcalatime = catime[0];
	  hcaltime = ctdctime[0];
	}

	//double hcalatime = cblkatime[0];
	double hcalatime_bc = hcalatime-hodotmean[0];
	double hcalatime_c = hcalatime - 
	  old_time_offsets[(int)pblkid] + 
	  new_time_offsets[(int)pblkid];
	double hcalatime_fc = hcalatime_c-hodotmean[0];
	//Assuming the tdc calib hasn't changed, which is true on SBS14 and SBS11 from pass0/1 to pass2
	//double hcaltime = cblktime[0];
	double hcaltime_bc = hcaltime-hodotmean[0];
	double hcaltime_c = hcaltime +
	  old_time_offsets[(int)pblkid] * old_tdc_calib - 
	  new_time_offsets[(int)pblkid] * old_tdc_calib;
	double hcaltime_fc = hcaltime_c-hodotmean[0];
	
	double time_before;
	double time_after;
	if(adct){
	  time_before = hcalatime_bc;
	  time_after = hcalatime_fc;
	}else{
	  time_before = hcaltime_bc;
	  time_after = hcaltime_fc;
	} 

	//Cut on nucleon spot
	if( !pass_spot )
	  continue;

	ht_id_before->Fill(pblkid,time_before);
	ht_id_after->Fill(pblkid,time_after);
	ht_row_before->Fill(pblkrow,time_before);
	ht_row_after->Fill(pblkrow,time_after);
	ht_run_before->Fill(current_runnumber,time_before);
	ht_run_after->Fill(current_runnumber,time_after);
	ht_allchan_allruns_before->Fill(time_before);
	ht_allchan_allruns_after->Fill(time_after);

      }//end loop over event

      // getting ready for the next run
      C->Reset();

    }//endloop over runs

  }//endloop over sets
  
  gStyle->SetPalette(55);

  //overlay adct vs run histo before new coeff
  TCanvas *c1 = new TCanvas("c1",Form("%s vs Run, Before",type.c_str()),1600,800);
  std::vector<double> c1y;
  std::vector<double> c1x;
  //gStyle->SetOptStat(0);
  //c1->cd();
  overlayWithMeans(ht_run_before,c1,c1y,c1x);

  c1->SetGridy();
  c1->Write();
  std::string c1_path = plotdir + Form("%svRun_before.png",type.c_str());
  c1->SaveAs(c1_path.c_str());

  //overlay adct vs run histo after new coeff
  TCanvas *c2 = new TCanvas("c2",Form("%s vs Run, After",type.c_str()),1600,800);
  vector<double> c2y;
  vector<double> c2x;
  //gStyle->SetOptStat(0);
  //c2->cd();
  overlayWithMeans(ht_run_after,c2,c2y,c2x);

  c2->SetGridy();
  c2->Write();
  std::string c2_path = plotdir + Form("%svRun_after.png",type.c_str());
  c2->SaveAs(c2_path.c_str());

  //overlay adct vs run histo before new coeff
  TCanvas *c3 = new TCanvas("c3",Form("%s vs ID, Before",type.c_str()),1600,800);
  std::vector<double> c3y;
  std::vector<double> c3x;
  //gStyle->SetOptStat(0);
  //c3->cd();
  overlayWithMeans(ht_id_before,c3,c3y,c3x);

  c3->SetGridy();
  c3->Write();
  std::string c3_path = plotdir + Form("%svID_before.png",type.c_str());
  c3->SaveAs(c3_path.c_str());

  //overlay adct vs run histo after new coeff
  TCanvas *c4 = new TCanvas("c4",Form("%s vs ID, After",type.c_str()),1600,800);
  vector<double> c4y;
  vector<double> c4x;
  //gStyle->SetOptStat(0);
  //c4->cd();
  overlayWithMeans(ht_id_after,c4,c4y,c4x);

  c4->SetGridy();
  c4->Write();
  std::string c4_path = plotdir + Form("%svID_after.png",type.c_str());
  c4->SaveAs(c4_path.c_str());

  //overlay adct vs row histo before new coeff
  TCanvas *c3r = new TCanvas("c3r",Form("%s vs Row, Before",type.c_str()),1600,800);
  std::vector<double> c3ry;
  std::vector<double> c3rx;
  //gStyle->SetOptStat(0);
  //c3r->cd();
  overlayWithMeans(ht_row_before,c3r,c3ry,c3rx);

  c3r->SetGridy();
  c3r->Write();
  std::string c3r_path = plotdir + Form("%svRow_before.png",type.c_str());
  c3r->SaveAs(c3r_path.c_str());

  //overlay adct vs row histo after new coeff
  TCanvas *c4r = new TCanvas("c4r",Form("%s vs Row, After",type.c_str()),1600,800);
  vector<double> c4ry;
  vector<double> c4rx;
  //gStyle->SetOptStat(0);
  //c4r->cd();
  overlayWithMeans(ht_row_after,c4r,c4ry,c4rx);

  c4r->SetGridy();
  c4r->Write();
  std::string c4r_path = plotdir + Form("%svRow_after.png",type.c_str());
  c4r->SaveAs(c4r_path.c_str());

  //set up before tgraph
  int Nbef = c1x.size();
  double *befx = &c1x[0];
  double *befy = &c1y[0];
  TGraph *gbef = new TGraph(Nbef,befx,befy);

  //set up after tgraph
  int Naf = c2x.size();
  double *afx = &c2x[0];
  double *afy = &c2y[0];
  TGraph *gaf = new TGraph(Naf,afx,afy);

  gbef->SetMarkerStyle(20);
  gbef->SetMarkerColor(kRed);
  //gbef->SetLineColor(kRed);
  gbef->GetYaxis()->SetRangeUser(-80,80);
  //gbef->Draw("P");

  gaf->SetMarkerStyle(20);
  gaf->SetMarkerColor(kBlack);
  //gaf->SetLineColor(kBlack);
  //gaf->Draw("P SAME");

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gbef);
  mg->Add(gaf);

  TCanvas *c5 = new TCanvas("c5",Form("%s mpv vs run before/after",type.c_str()),1600,800);
  c5->cd();

  c5->SetGridy();

  mg->Draw("AP");
  mg->SetTitle(Form("%s elastic signal peak vs run number; run number; ns",type.c_str()));
  
  TLegend *leg = new TLegend(0.75,0.75,0.89,0.89);
  leg->AddEntry(gbef, "Before Alignment", "p");
  leg->AddEntry(gaf, "After Alignment", "p");
  leg->Draw();

  c5->Update();
  c5->Write();
  std::string c5_path = plotdir + Form("%sMPVvRun_comp.png",type.c_str());
  c5->SaveAs(c5_path.c_str());

  //set up before id tgraph
  int Nidbef = c3x.size();
  double *idbefx = &c3x[0];
  double *idbefy = &c3y[0];
  TGraph *gidbef = new TGraph(Nidbef,idbefx,idbefy);

  //set up after id tgraph
  int Nidaf = c4x.size();
  double *idafx = &c4x[0];
  double *idafy = &c4y[0];
  TGraph *gidaf = new TGraph(Nidaf,idafx,idafy);

  gidbef->SetMarkerStyle(20);
  gidbef->SetMarkerColor(kRed);
  //gbef->SetLineColor(kRed);
  gidbef->GetYaxis()->SetRangeUser(-80,80);
  //gbef->Draw("P");

  gidaf->SetMarkerStyle(20);
  gidaf->SetMarkerColor(kBlack);
  //gaf->SetLineColor(kBlack);
  //gaf->Draw("P SAME");

  TMultiGraph *idmg = new TMultiGraph();
  idmg->Add(gidbef);
  idmg->Add(gidaf);

  TCanvas *c6 = new TCanvas("c6",Form("%s mpv vs ID before/after",type.c_str()),1600,800);
  c6->cd();

  c6->SetGridy();

  idmg->Draw("AP");
  idmg->SetTitle(Form("%s elastic signal peak vs ID; channel; ns",type.c_str()));
  
  TLegend *idleg = new TLegend(0.75,0.75,0.89,0.89);
  idleg->AddEntry(gidbef, "Before Alignment", "1p");
  idleg->AddEntry(gidaf, "After Alignment", "1p");
  idleg->Draw();

  c6->Update();
  c6->Write();
  std::string c6_path = plotdir + Form("%sMPVvID_comp.png",type.c_str());
  c6->SaveAs(c6_path.c_str());

  TCanvas *c7 = new TCanvas("c7",Form("%s All Channels All Runs before/after",type.c_str()),1600,800);
  c7->cd();

  c7->SetGridx();

  ht_allchan_allruns_before->SetLineColor(kRed);
  ht_allchan_allruns_before->SetLineWidth(2);
  ht_allchan_allruns_before->SetLineStyle(2);
  ht_allchan_allruns_before->Draw();
  ht_allchan_allruns_before->SetTitle(Form("%s-HODOtmean All Runs All Channels; ns",type.c_str()));
  
  ht_allchan_allruns_after->SetLineColor(kBlack);
  ht_allchan_allruns_after->SetLineWidth(2);
  ht_allchan_allruns_after->Draw("same");  

  TLegend *allleg = new TLegend(0.75,0.75,0.89,0.89);
  allleg->AddEntry(ht_allchan_allruns_before, "Before Alignment", "l");
  allleg->AddEntry(ht_allchan_allruns_after, "After Alignment", "l");
  allleg->Draw();

  gStyle->SetOptStat(0);

  c7->Update();
  c7->Write();
  std::string c7_path = plotdir + Form("%sAll_comp.png",type.c_str());
  c7->SaveAs(c7_path.c_str());
  
  fout->Write();

  st->Stop();

  cout << endl << "Timepatch complete on " << type << ". Output file written to " << fout_path << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    

}
