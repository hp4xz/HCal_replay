//SSeeds 9.1.23 - script adapted from tdc_align.C to allow for run ranges passed as arguments. Also reduces overhead where possible. Quality replay is expected to be done with qreplay_timepatch.C.

#include <vector>
#include <iostream>
#include <algorithm>
#include <vector>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

// TDC offset extraction constraints
const int first_hcal_chan = 0;
const int total_bins = 400;
const int lower_lim = -160;
const int upper_lim = 40;
const int fit_event_min = 50;
const double observed_tdc_sigma = 2.5; //rough estimate
const double TDC_target = 0.; //Target value for peak tdc.
const int linecount = 14;

//for small data sets, want to set bad fit signal mean to maximum value. true if there is a unique maximum, else false.
bool hasUniqueMaximum(TH1D* hist) {
    // Check if histogram is valid
    if (!hist) {
        std::cerr << "Invalid histogram pointer!" << std::endl;
        return false;
    }

    // Get the number of bins in the histogram
    int n_bins = hist->GetNbinsX();

    // Find the maximum value in the histogram
    double max_val = hist->GetMaximum();

    // Count how many times the maximum value appears
    int max_count = 0;
    for (int i = 1; i <= n_bins; i++) {
        if (hist->GetBinContent(i) == max_val) {
            max_count++;
        }
    }

    // If the maximum value appears more than once, it's not unique
    return max_count == 1;
}

void sliceFits(TH1D* hist, 
	       TCanvas* canvas, 
	       int canvasIndex, 
	       double &mpv, 
	       double &sigma, 
	       TFile* outputFile) {
  
  // Check if histogram is valid
  if (!hist) {
    std::cerr << "Invalid histogram pointer!" << std::endl;
    return;
  }
  if (!outputFile) {
    std::cerr << "Invalid TFile pointer!" << std::endl;
    return;
  }

  canvas->cd(canvasIndex);
  //Get slices from htDiff_ID and fit for mpv vals
  double fitl = 0.;
  double fith = 0.; 

  bool unique_maximum = hasUniqueMaximum(hist);
  int N = hist->GetEntries();
  int binmax = hist->GetMaximumBin();
  double binmaxX = hist->GetBinCenter(binmax);
  double binmaxY = hist->GetBinContent(binmax);
  fitl = binmaxX - 2*observed_tdc_sigma;
  fith = binmaxX + 3*observed_tdc_sigma;
  
  if( N<fit_event_min || binmaxX<=lower_lim || binmaxX>=upper_lim ){
    hist->SetLineColor(kRed);
    hist->SetTitle(Form( "Insufficient Statistics, max X unique = %d, max X %0.2f", unique_maximum, binmaxX ) );
    if(unique_maximum)
      mpv = binmaxX;
    hist->Draw();
    hist->Write();
    return;
  }

  // Define the skewed gaussian function for fitting
  TF1 *sgfit = new TF1("sgfit",util::g_sgfit_bg,fitl,fith,5);

  //Set and constrain the fit parameters
  sgfit->SetLineWidth(4);
  sgfit->SetLineColor(kGreen);
  sgfit->SetParameter(0,binmaxY);
  sgfit->SetParLimits(0,0,1.5*binmaxY);
  sgfit->SetParameter(1,binmaxX);
  sgfit->SetParLimits(1,fitl,fith);
  sgfit->SetParameter(2,observed_tdc_sigma);
  sgfit->SetParLimits(2,1.,3*observed_tdc_sigma);
  sgfit->SetParameter(3,1.1);
  sgfit->SetParLimits(3,0.3,8);
  sgfit->SetParameter(4,25);
  sgfit->SetParLimits(4,5,30);
  
  // Fit the histogram
  hist->Fit("sgfit","RBMQ");

  //Catch bad fits on very small data sets, assuming 0.01 is very small wrt sigma of fit    
  if( sgfit->GetMaximumX()>fith-0.01 || sgfit->GetMaximumX()<fitl+0.01 ){
    mpv = binmaxX;
    hist->SetLineColor(kOrange);
  }else{
    mpv = sgfit->GetMaximumX(); //for a skewed gaussian, parameter 1 isn't the MPV
  }

  sigma = sgfit->GetParameter(2);

  hist->SetTitle(Form("N:%d MaxX:%f MPV:%f Sigma:%f",N,binmaxX,mpv,sigma));

  // Check if the output file is valid and write the fitted histogram to it
  hist->Write();

  canvas->Update();

  delete sgfit;
}

//Main. If no arguments are passed to run_b and run_e, alignment over all runs will proceed.
void tdc_setalign( const char *experiment = "gmn", int config = 7, int pass = 2, int run_b = 0, int run_e = 0, int run_exclude_b = 0, int run_exclude_e = 0, bool lowstatopt = false, bool best_clus = false ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  std::string date = util::getDate();

  //switch if all runs or no exclusions should be considered
  bool runall = false;
  if( run_b==0 && run_e==0 )
    runall = true;
  bool noexclude = false;
  if( run_exclude_b==0 && run_exclude_e==0 )
    noexclude = true;

  //Gets database directory from system vars - USER change if necessary
  //std::string db_path = gSystem->Getenv("DB_DIR");
  std::string db_path = "/work/halla/sbs/seeds/alt_sbsreplay/SBS-replay/DB";

  std::string run_exclude_word = "";
  if( (run_exclude_b != run_exclude_e && run_exclude_b == 0) ||
      (run_exclude_b != run_exclude_e && run_exclude_e == 0) ||
      run_exclude_e < run_exclude_b ||
      (run_b != run_e && run_b == 0) ||
      (run_b != run_e && run_e == 0) ||
      run_e < run_b ){
    std::cout << "ERROR: Run or exclude range configured improperly. Reconfigure and retry." << std::endl;
    return;
  }

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string tdc_align_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/tdcsetalign_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.root",pass,experiment,config,pass,run_b,run_e,run_exclude_b,run_exclude_e);

  std::string new_tdcoffset_path = Form("parameters/tdcsetoffsets_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.txt",experiment,config,pass,run_b,run_e,run_exclude_b,run_exclude_e);

  std::string old_db_path = db_path + "/db_sbs.hcal.dat";
  std::string db_tdcoffset_variable = "sbs.hcal.tdc.offset";
  std::string db_tdccalib_variable = "sbs.hcal.tdc.calib";

  // Declare output analysis file
  TFile *fout = new TFile( tdc_align_path.c_str(), "RECREATE" );

  std::string plotdir = Form("../quality_plots/%s/conf%d/",experiment,config);

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  int nruns = -1; //Analyze all available runs for this configuration
  int verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
  
  //Set up structure to hold calibration parameters and initialize size/index
  calset tdc_cal;
  int Ncal_set_size = 0;
  int Ncal_set_idx = 0;
  
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  cout << config_parameters;
  double ebeam = config_parameters.GetEbeam();
  double hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  double hcaldist = config_parameters.GetHCALdist();
  double sbsdist = config_parameters.GetSBSdist();
  double bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  //Declare histograms for modification later where additional calibration sets are necessary
  TH2D *  htp_hodocorr_ID = new TH2D("htp_hodocorr_ID",
				     "tdc-hodotmean vs hcal id; chan; ns",
				     hcal::maxHCalChan,
				     first_hcal_chan,
				     hcal::maxHCalChan,
				     total_bins,
				     lower_lim,
				     upper_lim);
  
  TH1D *htp_allhodocorr = new TH1D("htp_allhodocorr",
				   "tdc-hodotmean all channels; ns",
				   total_bins,
				   lower_lim,
				   upper_lim);
    

  //minimize reporting
  int target_change_index;
  std::string set_timestamp = "";

  //Main loop over runs
  for (int r=0; r<nruns; r++) {

    //Get run experimental parameters
    std::string current_offset_timestamp = runs[r].tdc_ts;
    std::string current_calib_timestamp = runs[r].tdcc_ts;
    if( r==0 )
      set_timestamp = current_offset_timestamp;
    int current_runnumber = runs[r].runnum;

    //Address run range
    if( !runall && (current_runnumber < run_b || current_runnumber > run_e) )
      continue;

    //Address exclusion range
    if( !noexclude && (current_runnumber >= run_exclude_b && current_runnumber <= run_exclude_e) )
      continue;

    std::string current_target = runs[r].target;
    std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    int mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A

    //Get run paths
    std::string rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-%s/pass%d/SBS%d/%s/rootfiles/",experiment,pass,config,targ_uppercase.c_str());
    std::string rootfile_path = rootfile_dir + Form("*%d*",current_runnumber);

    //Get target configuration
    SBStarget target_parameters(current_target);
    int target_index = target_parameters.GetTargIndex();
    double target_length = target_parameters.GetTargLength();
    double target_rho = target_parameters.GetTargRho();
    double cell_rho = target_parameters.GetCellRho();
    double cell_diam = target_parameters.GetCellDiam();
    double cell_dEdx = target_parameters.GetCelldEdx();
    double upstream_wthick = target_parameters.GetUpstreamWallThick();
    double target_dEdx = target_parameters.GetTargdEdx();
    double M_avg = target_parameters.GetAvgMass();

    //send to console if the target changes
    if( target_change_index != target_index ){
      cout << target_parameters;
      target_change_index = target_index;
    }

    tdc_cal.timestamp = current_offset_timestamp.c_str();
    tdc_cal.calib_ts = current_calib_timestamp.c_str();
    util::readDB( old_db_path, current_offset_timestamp, db_tdcoffset_variable, tdc_cal.old_param ); 
    util::readDB( old_db_path, current_calib_timestamp, db_tdccalib_variable, tdc_cal.tdc_calib );
  
    //Get available cuts for current config/target/field combination. Use first element (0) of cut
    //vector<calcut> cut;
    vector<caldiag> cut;
    //util::ReadCutList(struct_dir,experiment,config,0,pass,current_target,mag,verb,cut);
    util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ,mag,verb,cut);  

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

    // Setting up chain and branch addresses
    C = new TChain("T");
    C->Add(rootfile_path.c_str());

    C->SetBranchStatus("*",0);    
    double BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
    double BBtr_vz[hcal::maxTracks];
    double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, BBsh_atime;	
    double HCALx, HCALy, HCALe;
    double pblkrow, pblkcol, nblk, nclus;
    int Ncid;
    double cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan], cblkagain[hcal::maxHCalChan];
    double cid[hcal::maxHCalClus], crow[hcal::maxHCalClus], ccol[hcal::maxHCalClus], ce[hcal::maxHCalClus], cx[hcal::maxHCalClus], cy[hcal::maxHCalClus], catime[hcal::maxHCalClus];
    double HODOtmean[hcal::maxHCalClus];

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
    C->SetBranchStatus( "bb.sh.atimeblk", 1 );
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
    C->SetBranchAddress( "bb.sh.atimeblk", &BBsh_atime ); 
    C->SetBranchAddress( "bb.hodotdc.clus.tmean", HODOtmean );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Ncid ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

    //Globalcut enables
    C->SetBranchStatus( "bb.tr.tg_th", 1 );
    C->SetBranchStatus( "bb.tr.tg_ph", 1 );

    //Use TTreeFormula to avoid looping over data an additional time
    TCut GCut = cut[0].gcut.c_str();

    //Add globalcut and elastic cuts for reporting
    tdc_cal.gcut = cut[0].gcut;
    tdc_cal.minEv = fit_event_min;

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

      double pblkid = (double)cblkid[0]-1;
      double hcaltdc = cblktime[0];

      double tdc_tc = hcaltdc-HODOtmean; //Primary cluster, primary block tdc with trigger correction (ns)
      
      //Best Cluster analysis using intime algorithm for simplicity on diagnosis

      //Set up best cluster index and energy for sorting
      Double_t best_cluster[2] = {-1,0.};

      //loop through all clusters and select without HCal position information
      for( int c=0; c<Nhcalcid; c++ ){
	
	//calculate h-arm physics quantities per cluster
	double atime = catime[c];
	double atime_hodo = atime - HODOtmean[0];
	double ce = ce[c];
	
	//wide cut around HCal/Hodo time
	bool passedTime = abs(atime_hodo-atime0)<3*atimesig;

	if( passedTime && ce>best_cluster[1] ){
	  best_cluster[0] = c;
	  best_cluster[1] = ce;
	}
      }//endloop over cluster elements
      
      //If no cluster passes atime cut, use default
      if(best_cluster[0]==-1)
	best_cluster[0]=0;

      //Switch between best clusters for systematic analysis
      Int_t cidx_best = best_cluster[0];
      
      //Calculations from the best cluster
      Double_t dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
      Double_t dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
      Double_t hatime_bestcluster = hcalcatime[cidx_best];
      Double_t htime_bestcluster = hcalctdctime[cidx_best];
      Double_t hcoin_bestcluster = hcalcatime[cidx_best] - atimeSH;
      Double_t ce_bestcluster = hcalce[cidx_best];
      Double_t x_bestcluster = hcalcx[cidx_best];
      Double_t y_bestcluster = hcalcy[cidx_best];
      Double_t row_bestcluster = hcalcrow[cidx_best];
      Double_t col_bestcluster = hcalccol[cidx_best];
      Double_t id_bestcluster = hcalcid[cidx_best];

      //prevent corruption of later logic which selects histograms to fit based on total number of events
      if(tdc_tc<-800)
	continue;

      htp_hodocorr_ID->Fill( pblkid, tdc_tc );
      htp_allhodocorr->Fill( tdc_tc );
      
    }//end loop over event
    
    // getting ready for the next run
    C->Reset();
    
  }//endloop over runs
  
  cout << "Ended loop over events. Proceeding to fits.." << endl;

  // Declare canvas and graph for plots
  TCanvas *TDC_top = new TCanvas("tdc_fit_top",Form("TDC_top, set TS: %s",set_timestamp.c_str()),1600,1200);
  TCanvas *TDC_bot = new TCanvas("tdc_fit_bot",Form("TDC_bot, set TS: %s",set_timestamp.c_str()),1600,1200);
  TCanvas *TDC_allchan = new TCanvas("tdc_fit_allchan",Form("TDC_fit_allchan, set TS: %s",set_timestamp.c_str()),1600,1200);;
  TDC_top->Divide(12,12);
  TDC_bot->Divide(12,12);
  gStyle->SetOptStat(0);

  // Set up arrays and vars for tgraph
  Double_t calib_const = tdc_cal.tdc_calib;

  //No error on cell location
  double cellerr[hcal::maxHCalChan] = {0.};

  //Make arrays for tdc tgraphs
  double tcell[hcal::maxHCalChan];
  double tcval[hcal::maxHCalChan];
  double tcerr[hcal::maxHCalChan];
  TH1D *tcellslice[hcal::maxHCalChan];
  double alltcval = 0.;
  double alltcerr = 0.;
  
  sliceFits(htp_allhodocorr,TDC_allchan,1,alltcval,alltcerr,fout);

  //loop over top half of hcal for tdc fits by channel
  for(int c=0; c<hcal::maxHCalChan; c++){
    //initialize the graph arrays
    tcerr[c] = 0.;
    tcval[c] = 0.; 
    tcell[c] = c;
    
    double mpv = -1.;
    double sigma = 0.;

    int half_chan = hcal::maxHCalChan/2;
      
    //Get slices from htp_hodocorr_ID and fit for most probable vals
    tcellslice[c] = htp_hodocorr_ID->ProjectionY(Form("tcellslice_%d",c+1),c+1,c+1);
      
    if( c < half_chan )
      sliceFits(tcellslice[c],TDC_top,c+1,mpv,sigma,fout);
    else
      sliceFits(tcellslice[c],TDC_bot,c+1-half_chan,mpv,sigma,fout);

    tcval[c] = mpv;
    tcerr[c] = sigma;

  }    
  TDC_top->Write();
  TDC_bot->Write();
  TDC_allchan->Write();

  // Output text file for new TDC offsets
  ofstream writeParFile;
  writeParFile.open( new_tdcoffset_path );

  TCanvas *TDC_graph = new TCanvas(Form("TDC_graph, offset TS: %s",set_timestamp.c_str()),"TDC vs ID, means",1600,1200);
  TGraphErrors *gtdc_c = new TGraphErrors( hcal::maxHCalChan, tcell, tcval, cellerr, tcerr );;

  TDC_graph->cd();

  //Make graphs with errors for reporting. All failed fits are zero here
  gtdc_c->GetXaxis()->SetLimits(-10,290);  
  gtdc_c->GetYaxis()->SetLimits(lower_lim,upper_lim);
  gtdc_c->SetTitle(Form("TDC_{hcal}-TDCmean_{hodo} vs Cell, set TS: %s",set_timestamp.c_str()));
  gtdc_c->GetXaxis()->SetTitle("Cell");
  gtdc_c->GetYaxis()->SetTitle("TDC_{HCAL}-TDC_{MEAN,HODO}");
  gtdc_c->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  gtdc_c->Draw();
  gtdc_c->Write("gtdc_c");  

  //Write to output text file
  writeParFile << "#HCal TDC offsets obtained " << date.c_str() << " for db tdc offset timestamp " << set_timestamp <<  endl;
  writeParFile << "#Offsets obtained from fits over TDC distributions." << endl;
  writeParFile << set_timestamp << endl;
  writeParFile << db_tdcoffset_variable << " =" << endl;
  
  int cell = 0;
  
  for( int r = 0; r<hcal::maxHCalRows; r++){
    for( int c = 0; c<hcal::maxHCalCols; c++){
      
      if( lowstatopt )
	writeParFile << alltcval/calib_const + tdc_cal.old_param[cell] - TDC_target/calib_const << " "; //Applies offset implied by distribution over all channels to each individually. Will not produce reasonable results if the offsets are not already reasonably aligned from a broader data set.
      else if(  tcval[cell]==-1000 || abs(tcerr[cell])<0.05  ){ //This will only work if most channels are closely aligned already and oldtdcoffsets aren't believable
	writeParFile << alltcval/calib_const + tdc_cal.old_param[cell] - TDC_target/calib_const << " ";
      }else{
	writeParFile << tcval[cell]/calib_const + tdc_cal.old_param[cell] - TDC_target/calib_const << " ";
      }
      cell++;
    } // endloop over columns 
    writeParFile << endl;
  } // endloop over rows
  writeParFile << endl << endl;
  
  //Add output report canvas
  TCanvas *reportcanvas = new TCanvas("tdcalignreport", "Configuration/Cut Information, Calibration", 1800, 900);
  
  // Set margin.
  reportcanvas->SetLeftMargin(0.01);
  
  // Create a TText object.
  TText *t = new TText();
  
  // Set text align to the left (horizontal alignment = 1).
  t->SetTextAlign(11);
  t->SetTextSize(0.02);
  
  //make an array of strings
  std::string target_option = "All Available";
  std::string report[linecount] = {
    "General Set TDC Alignment Info",
    Form( "Experiment: %s, Configuration: %d, Pass: %d", experiment, config, pass ),
    Form( "Creation Date: %s", date.c_str() ),
    Form( "Run range %d - %d", run_b, run_e ),
    Form( "Exclusion range %d - %d", run_exclude_b, run_exclude_e ),

    Form( "Target(s) Used: %s", target_option.c_str() ),
    Form( "Calibration Set: %s", tdc_cal.timestamp.c_str() ),
    "",
    "Elastic Cuts",
    Form( "Global Elastic Cuts: %s", tdc_cal.gcut.c_str() ),
    "",
    "Other Cuts",
    Form( "Minimum Ev per Cell : %d", tdc_cal.minEv ),
    "HCal Acceptance Match (Projected Nucleon Within HCal Acceptance)"
  };
  // Loop to write the lines to the canvas.
  for( int i = 0; i<linecount; i++ ) {
    // Vertical position adjusted according to line number.
    double verticalPosition = 0.9 - i * 0.04;
    t->DrawTextNDC(0.1, verticalPosition, report[i].c_str());
  }

  reportcanvas->Write();    

  writeParFile.close();
  fout->Write();
  st->Stop();

  cout << "Analysis and fits written to " << tdc_align_path << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
