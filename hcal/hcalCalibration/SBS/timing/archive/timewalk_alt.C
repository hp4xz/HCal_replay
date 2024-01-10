//sseeds 8.16.23 Alternative version of timewalk.C with a more fundamental scripting structure (via tdc_align.C) with aim to fix irregular seg faults.

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

// TDC offset extraction constraints
const Int_t first_hcal_chan = 0;
const Int_t total_bins = 320;
const Int_t lower_lim = -120;
const Int_t upper_lim = 40;
const Int_t fit_event_min = 50;
const Double_t observed_tdc_sigma = 2.5; //rough estimate
const Double_t TDC_target = 0.; //Target value for peak tdc.
const Int_t linecount = 12;

//Main <experiment> <configuration> <quasi-replay> <replay-pass> <target_option>; qreplay should only be performed after new offsets obtained
void timewalk_alt( const char *experiment = "gmn", Int_t config = 7, bool qreplay = false, Int_t pass = 0, bool h2only = true ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  std::string date = util::getDate();

  // Set up universal paths and variables
  std::string h2opt = "";
  if( h2only )
    h2opt = "_lh2only";
  std::string db_path = gSystem->Getenv("DB_DIR");
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string tdc_align_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/timewalkalt%s_%s_conf%d_qr%d_pass%d.root",pass,h2opt.c_str(),experiment,config,(Int_t)qreplay,pass);

  std::string new_tdcoffset_path = Form("parameters/timewalkalt_%s_conf%d_pass%d.txt",experiment,config,pass);
  //std::string new_tdcoffset_path = "test.txt";
  std::string old_db_path = db_path + "/db_sbs.hcal.dat";
  std::string db_tdcoffset_variable = "sbs.hcal.tdc.offset";
  std::string db_tdccalib_variable = "sbs.hcal.tdc.calib";

  // Declare output analysis file
  //TFile *fout = new TFile( tdc_align_path.c_str(), "RECREATE" );

  TFile *fout = new TFile("test.root","RECREATE");

  std::string plotdir = Form("../quality_plots/%s/conf%d/",experiment,config);

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
  
  //Set up structure to hold calibration parameters and initialize size/index
  calset tdc_cal[hcal::gNstamp];
  Int_t Ncal_set_size = 0;
  Int_t Ncal_set_idx = 0;
  
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  
  // create output tree
  TTree *P = new TTree("P","Analysis Tree");

  // Timing
  Double_t pblkid_out; //hcal cluster, primary block id
  Double_t tdc_out; //hcal cluster tdc time, tree
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

  //cluster block variables
  Double_t cblkid_out[hcal::maxHCalBlk];
  Double_t cblkatime_out[hcal::maxHCalBlk];
  Double_t cblktime_out[hcal::maxHCalBlk];
  Double_t cblke_out[hcal::maxHCalBlk];

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

  //Declare histograms for modification later where additional calibration sets are necessary
  TH2D *htp_hodocorr_ID[hcal::gNstamp];

  for( Int_t set=0; set<hcal::gNstamp; set++ ){
    htp_hodocorr_ID[set] = new TH2D(Form("htp_hodocorr_ID_set%d",set),
				    "null",
				    hcal::maxHCalChan,
				    first_hcal_chan,
				    hcal::maxHCalChan,
				    total_bins,
				    lower_lim,
				    upper_lim);
    
  }

  //minimize reporting
  Int_t target_change_index;

  //Main loop over runs
  for (Int_t r=0; r<nruns; r++) {

    //Get run experimental parameters
    std::string current_offset_timestamp = runs[r].tdc_ts;
    std::string current_calib_timestamp = runs[r].tdcc_ts;
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
    
    //check target and continue on deuterium of only calibrating with hydrogen
    if( h2only && target_index!=1 )
      continue;

    //send to console if the target changes
    if( target_change_index != target_index ){
      cout << target_parameters;
      target_change_index = target_index;
    }

    //Record new offset parameters from qreplay==false. Repopulate array on each run for changing timestamps
    Double_t new_tdc_offsets[hcal::maxHCalChan] = {0.};
    if( qreplay ){ //Check which timestamp is newest (DB:tdcoffset, DB:tdccalib, new)
      std::string active_timestamp;
      util::tsCompare( current_offset_timestamp, config_parameters.GetSBSTimestamp(), active_timestamp );
      util::tsCompare( active_timestamp, current_calib_timestamp, active_timestamp );
      util::readDB( new_tdcoffset_path, active_timestamp, db_tdcoffset_variable, new_tdc_offsets );
    }

    //////////
    //Select correct calibration set.
    bool new_calib_ts = true;
    bool new_offset_ts = true;
    //Look through available sets for duplicate timestamps
    for( Int_t cs=0; cs<Ncal_set_size; cs++ ){ 

      if( tdc_cal[cs].timestamp == current_offset_timestamp )
	new_offset_ts = false;
      if( tdc_cal[cs].calib_ts == current_calib_timestamp )
	new_calib_ts = false;
      //set the index where both exist
      if( tdc_cal[cs].calib_ts == current_calib_timestamp &&
	  tdc_cal[cs].timestamp == current_offset_timestamp )
	Ncal_set_idx = cs;
    }

    //if a new timestamp exists, expand the set and update index/old-offsets/old-calib
    if( new_calib_ts || new_offset_ts || Ncal_set_size==0){

      cout << "Adding new calibration set. Offset ts:" << current_offset_timestamp << ", calib ts:" << current_calib_timestamp << endl << endl;

      Ncal_set_idx = Ncal_set_size;

      if( Ncal_set_size == hcal::gNstamp-1 ){
	cout << "Error: Allowable calibration sets exceed maximum of " << hcal::gNstamp << ". Check timestamps and adjust limits as necessary." << endl;
	return;
      }

      tdc_cal[Ncal_set_size].timestamp = current_offset_timestamp.c_str();
      tdc_cal[Ncal_set_size].calib_ts = current_calib_timestamp.c_str();
      util::readDB( old_db_path, current_offset_timestamp, db_tdcoffset_variable, tdc_cal[Ncal_set_size].old_param ); 
      util::readDB( old_db_path, current_calib_timestamp, db_tdccalib_variable, tdc_cal[Ncal_set_size].tdc_calib );
      htp_hodocorr_ID[Ncal_set_size]->SetTitle(Form("TDC Primary Block - TDC hodo, offset ts: %s, calib ts: %s;Channel;TDC_{HCAL}-TDC_{HODO} (ns)",tdc_cal[Ncal_set_size].timestamp.c_str(),tdc_cal[Ncal_set_size].calib_ts.c_str()));

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

    //Globalcut enables
    C->SetBranchStatus( "bb.tr.tg_th", 1 );
    C->SetBranchStatus( "bb.tr.tg_ph", 1 );

    //Use TTreeFormula to avoid looping over data an additional time
    TCut GCut = cut[0].gcut.c_str();

    //Add globalcut and elastic cuts for reporting
    tdc_cal[Ncal_set_idx].gcut = cut[0].gcut;
    tdc_cal[Ncal_set_idx].minEv = fit_event_min;

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

      cout << "Analyzing run " << current_runnumber << ": " <<  nevent << "/" << nevents << " \r";
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

      Double_t pblkid = (double)cblkid[0]-1;
      Double_t hcaltdc = cblktime[0];
      //correct the tree tdc time with new offset
      if( qreplay ) 
	hcaltdc = hcaltdc + 
		       tdc_cal[Ncal_set_idx].old_param[(Int_t)pblkid] * tdc_cal[Ncal_set_idx].tdc_calib - 
		       new_tdc_offsets[(Int_t)pblkid] * tdc_cal[Ncal_set_idx].tdc_calib;

      Double_t tdc_tc = hcaltdc-HODOtmean; //Primary cluster, primary block tdc with trigger correction (ns)
      
      htp_hodocorr_ID[Ncal_set_idx]->Fill( pblkid, tdc_tc );
      
      pblkid_out = pblkid;
      tdc_out = cblktime[0];
      atime_out = cblkatime[0];
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
  
  // Declare canvas and graph for plots
  TCanvas *TDC_top[Ncal_set_size];
  TCanvas *TDC_bot[Ncal_set_size];

  // Set up arrays for tgraph
  //No error on cell location
  Double_t cellerr[hcal::maxHCalChan] = {0.};

  //Make arrays for tdc tgraphs
  Double_t tcell[Ncal_set_size][hcal::maxHCalChan];
  Double_t tcval[Ncal_set_size][hcal::maxHCalChan];
  Double_t tcvalw[Ncal_set_size][hcal::maxHCalChan];
  Double_t tcerr[Ncal_set_size][hcal::maxHCalChan];
  TH1D *tcellslice[Ncal_set_size][hcal::maxHCalChan];

  //Get averages for empty channels
  Double_t tcval_avg[Ncal_set_size];
  Int_t tcval_Ng[Ncal_set_size];

  for( Int_t set=0; set<Ncal_set_size; set++ ){
    
    //Define calibration set by newest timestamp (config/tdcoffset/tdccalib) and record constant
    Double_t calib_const = tdc_cal[set].tdc_calib;
    std::string offset_timestamp = tdc_cal[set].timestamp;
    std::string calib_timestamp = tdc_cal[set].calib_ts;
    std::string config_timestamp = config_parameters.GetSBSTimestamp();
    std::string active_timestamp;
    util::tsCompare( offset_timestamp, calib_timestamp, active_timestamp );
    util::tsCompare( active_timestamp, config_timestamp, active_timestamp );

    tcval_avg[set] = 0.;
    tcval_Ng[set] = 0;

    //Set up canvases for quality check
    TDC_top[set] = new TCanvas(Form("tdc_fit_top_set%d",set),Form("TDC_top, set TS: %s, offset TS: %s, calib TS: %s",active_timestamp.c_str(),offset_timestamp.c_str(),calib_timestamp.c_str()),1600,1200);
    TDC_bot[set] = new TCanvas(Form("tdc_fit_bot_set%d",set),Form("TDC_bot, set TS: %s, offset TS: %s, calib TS: %s",active_timestamp.c_str(),offset_timestamp.c_str(),calib_timestamp.c_str()),1600,1200);
    TDC_top[set]->Divide(12,12);
    TDC_bot[set]->Divide(12,12);
    gStyle->SetOptStat(0);
    
    for(Int_t c=0; c<hcal::maxHCalChan; c++){
      //initialize the graph arrays
      tcvalw[set][c] = 0.;
      tcerr[set][c] = 0.;

      Int_t half_chan = hcal::maxHCalChan/2;

      //Index through the canvas
      TDC_top[set]->cd(c+1);
      if( c>=half_chan )
	TDC_bot[set]->cd(c+1-half_chan);
      
      //Get slices from htDiff_ID and fit for mean vals
      Double_t tfitl = 0.;
      Double_t tfith = 0.;
      tcell[set][c] = c;
      tcellslice[set][c] = htp_hodocorr_ID[set]->ProjectionY(Form("tcellslice_%d_%d",set,c+1),c+1,c+1);
      tcval[set][c] = tdc_cal[set].old_param[c]; 
      
      tcvalw[set][c] = 0; 

      Int_t sliceN = tcellslice[set][c]->GetEntries();
      Int_t binmax = tcellslice[set][c]->GetMaximumBin();
      Double_t binmaxX = lower_lim+binmax*(upper_lim-lower_lim)/total_bins;
      Double_t binmaxY = tcellslice[set][c]->GetBinContent(binmax);
      tfitl = binmaxX - 4*observed_tdc_sigma;
      tfith = binmaxX + 4*observed_tdc_sigma;

      if( sliceN<fit_event_min || binmaxX<=lower_lim || binmaxX>=upper_lim ){
	tcellslice[set][c]->Draw();
	continue;
      }

      TF1 *sgausfit = new TF1("sgausfit",util::g_sgfit_bg,tfitl,tfith,5);
      sgausfit->SetLineWidth(4);
      sgausfit->SetParameter(0,binmaxY);
      sgausfit->SetParameter(1,binmaxX);
      sgausfit->SetParLimits(1,tfitl,tfith);
      sgausfit->SetParameter(2,observed_tdc_sigma);
      sgausfit->SetParLimits(2,1.,3*observed_tdc_sigma);
      sgausfit->SetParameter(3,1.1);
      sgausfit->SetParLimits(3,0.3,8);
      sgausfit->SetParameter(4,25);
      sgausfit->SetParLimits(4,5,30);
      tcellslice[set][c]->Fit("sgausfit","RBMQ");

      tcellslice[set][c]->Draw();
      
      tcval[set][c] = sgausfit->GetParameter(1);
      tcvalw[set][c] = sgausfit->GetParameter(1);
      tcerr[set][c] = sgausfit->GetParameter(2);
      tcellslice[set][c]->SetTitle(Form("Set:%d N:%d MaxX:%f Mean:%f Sigma:%f Alpha:%f",set,sliceN,binmaxX,tcval[set][c],tcerr[set][c],sgausfit->GetParameter(3)));    
      tcellslice[set][c]->Write();

      tcval_avg[set] += tcval[set][c];
      tcval_Ng[set]++;

    }    
    TDC_top[set]->Write();
    TDC_bot[set]->Write();

    // if( !qreplay ){
    //   std::string qplotpath_top = plotdir + Form("tdc_fit_top_set%d.png",set);
    //   std::string qplotpath_bot = plotdir + Form("tdc_fit_bot_set%d.png",set);
    //   TDC_top[set]->SaveAs(qplotpath_top.c_str());
    //   TDC_top[set]->SaveAs(qplotpath_bot.c_str());
    // } 

    tcval_avg[set] /= tcval_Ng[set];

  }//endloop over slice fits

  TCanvas *TDC_graph[Ncal_set_size];
  TGraphErrors *gtdc_c[Ncal_set_size];

  // Output text file for new TDC offsets
  ofstream writeParFile;
  if( qreplay )
    new_tdcoffset_path = "/dev/null"; //Safety to prevent overwriting constants on quasi-replay
  writeParFile.open( new_tdcoffset_path );

  // Get tgraph and record tdc offsets
  for( Int_t set=0; set<Ncal_set_size; set++ ){

    //Define calibration set by newest timestamp (config/tdcoffset/tdccalib) and record constant
    Double_t calib_const = tdc_cal[set].tdc_calib;
    std::string offset_timestamp = tdc_cal[set].timestamp;
    std::string calib_timestamp = tdc_cal[set].calib_ts;
    std::string config_timestamp = config_parameters.GetSBSTimestamp();
    std::string active_timestamp;
    util::tsCompare( offset_timestamp, calib_timestamp, active_timestamp );
    util::tsCompare( active_timestamp, config_timestamp, active_timestamp );

    TDC_graph[set] = new TCanvas(Form("TDC_graph, set TS: %s, offset TS: %s, calib TS: %s",active_timestamp.c_str(),offset_timestamp.c_str(),calib_timestamp.c_str()),"TDC vs ID, means",1600,1200);
    TDC_graph[set]->cd();

    //Make graphs with errors for reporting. All failed fits are zero here
    gtdc_c[set] = new TGraphErrors( hcal::maxHCalChan, tcell[set], tcvalw[set], cellerr, tcerr[set] );
    gtdc_c[set]->GetXaxis()->SetLimits(-10,290);  
    gtdc_c[set]->GetYaxis()->SetLimits(lower_lim,upper_lim);
    gtdc_c[set]->SetTitle(Form("TDC_{hcal}-TDCmean_{hodo} vs Cell, set TS: %s, offset TS: %s, calib TS: %s, qreplay: %d",active_timestamp.c_str(),offset_timestamp.c_str(),calib_timestamp.c_str(),(Int_t)qreplay ) );
    gtdc_c[set]->GetXaxis()->SetTitle("Cell");
    gtdc_c[set]->GetYaxis()->SetTitle("TDC_{HCAL}-TDC_{MEAN,HODO}");
    gtdc_c[set]->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
    gtdc_c[set]->Draw();
    gtdc_c[set]->Write(Form("gtdc_c_s%d",set));  

    //Write to output text file
    if( qreplay ) continue; //Don't write if quasi-replay
    writeParFile << "#HCal TDC offsets obtained " << date.c_str() << " for db tdc offset timestamp " << offset_timestamp << " and for tdc calib timestamp " << calib_timestamp << endl;
    writeParFile << "#Offsets obtained from fits over TDC distributions." << endl;
    writeParFile << active_timestamp << endl;
    writeParFile << db_tdcoffset_variable << " =" << endl;
  
    Int_t cell = 0;

    for( Int_t r = 0; r<hcal::maxHCalRows; r++){
      for( Int_t c = 0; c<hcal::maxHCalCols; c++){
	
	if( tdc_cal[set].old_param[cell]>0 && ( tcval[set][cell]==0 || abs(tcerr[set][cell])<0.05 ) ){ //This will only work if most channels are closely aligned already and oldtdcoffsets aren't believable
	  writeParFile << tcval_avg[set]/calib_const + tdc_cal[set].old_param[cell] - TDC_target/calib_const << " ";
	}else{
	  writeParFile << tcval[set][cell]/calib_const + tdc_cal[set].old_param[cell] - TDC_target/calib_const << " ";
	}
	cell++;
      } // endloop over columns 
      writeParFile << endl;
    } // endloop over rows
    writeParFile << endl << endl;
    
  } // endloop over calibration sets

  writeParFile << endl << endl;

    //Add output report canvas
  TCanvas *c1[Ncal_set_size];
  Int_t report_height = 1800;

  for( Int_t s=0; s<Ncal_set_size; s++ ){

    c1[s] = new TCanvas(Form("tdcalignreport_set%d",s), Form("Configuration/Cut Information, Calibration Set %d",s), report_height, 900);
  
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
      "General TDC Alignment Info",
      Form("Experiment: %s, Configuration: %d, Pass: %d", experiment, config, pass),
      Form("Creation Date: %s", date.c_str() ),
      Form("Target(s) Used: %s", target_option.c_str() ),
      Form("Calibration Set: %s", tdc_cal[s].timestamp.c_str() ),
      "",
      "Elastic Cuts",
      Form("Global Elastic Cuts: %s", tdc_cal[s].gcut.c_str() ),
      "",
      "Other Cuts",
      Form("Minimum Ev per Cell : %d", tdc_cal[s].minEv),
      "HCal Acceptance Match (Projected Nucleon Within HCal Acceptance)"
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
    delete htp_hodocorr_ID[unused_set_index];
  }

  writeParFile.close();
  fout->Write();
  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
