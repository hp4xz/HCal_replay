//sseeds 11.23.23 - Diagnostic script to produce histograms for later alignments. Intended to be run several times to determine need for constrained run ranges. First run with all run_* options = 0.

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

// time histogram limits
const Double_t min_t = -100.;
const Double_t max_t = 100.;
const Double_t min_tdct = -100.;
const Double_t max_tdct = 100.;
const Int_t bins_t = 400;
const Int_t intime_sig = 5;
const Int_t spot_sig = 3;

bool verb = false;

//MAIN
void timing_diagnostic( const char *experiment = "gen", Int_t kine=2, Int_t pass=1, int run_b = 0, int run_e = 0, int run_exclude_b = 0, int run_exclude_e = 0, bool all_targets=false )
{   

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  std::string date = util::getDate();

  //One set of data files in json config shared between pass0/1 per kinematic
  if( pass==0 )
    pass=1;

  //switch if all runs or no exclusions should be considered
  bool runall = false;
  if( run_b==0 && run_e==0 )
    runall = true;
  bool noexclude = false;
  if( run_exclude_b==0 && run_exclude_e==0 )
    noexclude = true;

  //run inclusion/exclusion logical integrity check
  if( (run_exclude_b != run_exclude_e && run_exclude_b == 0) ||
      (run_exclude_b != run_exclude_e && run_exclude_e == 0) ||
      run_exclude_e < run_exclude_b ||
      (run_b != run_e && run_b == 0) ||
      (run_b != run_e && run_e == 0) ||
      run_e < run_b ){
    std::cerr << "ERROR: Run or exclude range configured improperly. Reconfigure arguments and retry." << std::endl;
    return;
  }

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s_p%d/",experiment,pass);

  Int_t nruns = -1; //Analyze all available runs for this configuration

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,kine,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available

  // outfile path
  std::string out_dir = gSystem->Getenv("OUT_DIR");
  std::string out_path = out_dir + Form("/hcal_calibrations/pass%d/diagnostics/timingdiag_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.root",pass,experiment,kine,pass,run_b,run_e,run_exclude_b,run_exclude_e);

  //set up output files
  TFile *fout = new TFile( out_path.c_str(), "RECREATE" );

  //set up diagnostic histograms
  TH2D *hdxdy = new TH2D( "hdxdy", "HCal dx vs dy; dy; dx", 400, -2, 2, 500, -3, 2 );
  TH2D *hdxdy_cut = new TH2D( "hdxdy_cut", "HCal dx vs dy (spot cut); dy; dx", 400, -2, 2, 500, -3, 2 );
  TH2D *hdxdy_bc = new TH2D( "hdxdy_bc", "HCal dx vs dy (best cluster); dy; dx", 400, -2, 2, 500, -3, 2 );
  TH2D *hdxdy_cut_bc = new TH2D( "hdxdy_cut_bc", "HCal dx vs dy (best cluster, spot cut); dy; dx", 400, -2, 2, 500, -3, 2 );

  TH2D *hatimerun_ecut = new TH2D( "hatimerun_ecut", "HCal-Hodo Atime vs Run Num (Electron Arm Cuts Only); Run; ns", nruns, 0, nruns, bins_t, min_t, max_t );
  TH2D *hatimerow_ecut = new TH2D( "hatimerow_ecut", "HCal-Hodo Atime vs HCal Row (Electron Arm Cuts Only); Row; ns", hcal::maxHCalRows, 0, hcal::maxHCalRows, bins_t, min_t, max_t );
  TH2D *hatimecol_ecut = new TH2D( "hatimecol_ecut", "HCal-Hodo Atime vs HCal Col (Electron Arm Cuts Only); Col; ns", hcal::maxHCalCols, 0, hcal::maxHCalCols, bins_t, min_t, max_t );
  TH2D *hatimeid_ecut = new TH2D( "hatimeid_ecut", "HCal-Hodo Atime vs HCal ID (Electron Arm Cuts Only); ID; ns", hcal::maxHCalChan, 0, hcal::maxHCalChan, bins_t, min_t, max_t );
  TH2D *hatimex_ecut = new TH2D( "hatimex_ecut", "HCal-Hodo Atime vs HCal X (Electron Arm Cuts Only); Row; ns", 500, -3, 2, bins_t, min_t, max_t );
  TH2D *hatimey_ecut = new TH2D( "hatimey_ecut", "HCal-Hodo Atime vs HCal Y (Electron Arm Cuts Only); Row; ns", 400, -2, 2, bins_t, min_t, max_t );

  TH2D *hatimerun = new TH2D( "hatimerun", "HCal-Hodo Atime vs Run Num; Run; ns", nruns, 0, nruns, bins_t, min_t, max_t );
  TH2D *hatimerow = new TH2D( "hatimerow", "HCal-Hodo Atime vs HCal Row; Row; ns", hcal::maxHCalRows, 0, hcal::maxHCalRows, bins_t, min_t, max_t );
  TH2D *hatimecol = new TH2D( "hatimecol", "HCal-Hodo Atime vs HCal Col; Col; ns", hcal::maxHCalCols, 0, hcal::maxHCalCols, bins_t, min_t, max_t );
  TH2D *hatimeid = new TH2D( "hatimeid", "HCal-Hodo Atime vs HCal ID; ID; ns", hcal::maxHCalChan, 0, hcal::maxHCalChan, bins_t, min_t, max_t );
  TH2D *hatimex = new TH2D( "hatimex", "HCal-Hodo Atime vs HCal X; Row; ns", 500, -3, 2, bins_t, min_t, max_t );
  TH2D *hatimey = new TH2D( "hatimey", "HCal-Hodo Atime vs HCal Y; Row; ns", 400, -2, 2, bins_t, min_t, max_t );

  TH2D *hatimerun_bc = new TH2D( "hatimerun_bc", "HCal-Hodo Atime vs Run Num (intime cluster); Run; ns", nruns, 0, nruns, bins_t, min_t, max_t );
  TH2D *hatimerow_bc = new TH2D( "hatimerow_bc", "HCal-Hodo Atime vs HCal Row (intime cluster); Row; ns", hcal::maxHCalRows, 0, hcal::maxHCalRows, bins_t, min_t, max_t );
  TH2D *hatimecol_bc = new TH2D( "hatimecol_bc", "HCal-Hodo Atime vs HCal Col (intime cluster); Col; ns", hcal::maxHCalCols, 0, hcal::maxHCalCols, bins_t, min_t, max_t );
  TH2D *hatimeid_bc = new TH2D( "hatimeid_bc", "HCal-Hodo Atime vs HCal ID (intime cluster); ID; ns", hcal::maxHCalChan, 0, hcal::maxHCalChan, bins_t, min_t, max_t );
  TH2D *hatimex_bc = new TH2D( "hatimex_bc", "HCal-Hodo Atime vs HCal X (intime cluster); Row; ns", 500, -3, 2, bins_t, min_t, max_t );
  TH2D *hatimey_bc = new TH2D( "hatimey_bc", "HCal-Hodo Atime vs HCal Y (intime cluster); Row; ns", 400, -2, 2, bins_t, min_t, max_t );

  TH2D *htimerun_ecut = new TH2D( "htimerun_ecut", "HCal-Hodo Time vs Run Num (Electron Arm Cuts Only); Run; ns", nruns, 0, nruns, bins_t, min_tdct, max_tdct );
  TH2D *htimerow_ecut = new TH2D( "htimerow_ecut", "HCal-Hodo Time vs HCal Row (Electron Arm Cuts Only); Row; ns", hcal::maxHCalRows, 0, hcal::maxHCalRows, bins_t, min_tdct, max_tdct );
  TH2D *htimecol_ecut = new TH2D( "htimecol_ecut", "HCal-Hodo Time vs HCal Col (Electron Arm Cuts Only); Col; ns", hcal::maxHCalCols, 0, hcal::maxHCalCols, bins_t, min_tdct, max_tdct );
  TH2D *htimeid_ecut = new TH2D( "htimeid_ecut", "HCal-Hodo Time vs HCal ID (Electron Arm Cuts Only); ID; ns", hcal::maxHCalChan, 0, hcal::maxHCalChan, bins_t, min_tdct, max_tdct );
  TH2D *htimex_ecut = new TH2D( "htimex_ecut", "HCal-Hodo Atime vs HCal X (Electron Arm Cuts Only); Row; ns", 500, -3, 2, bins_t, min_tdct, max_tdct );
  TH2D *htimey_ecut = new TH2D( "htimey_ecut", "HCal-Hodo Atime vs HCal Y (Electron Arm Cuts Only); Row; ns", 400, -2, 2, bins_t, min_tdct, max_tdct );

  TH2D *htimerun = new TH2D( "htimerun", "HCal-Hodo Time vs Run Num; Run; ns", nruns, 0, nruns, bins_t, min_tdct, max_tdct );
  TH2D *htimerow = new TH2D( "htimerow", "HCal-Hodo Time vs HCal Row; Row; ns", hcal::maxHCalRows, 0, hcal::maxHCalRows, bins_t, min_tdct, max_tdct );
  TH2D *htimecol = new TH2D( "htimecol", "HCal-Hodo Time vs HCal Col; Col; ns", hcal::maxHCalCols, 0, hcal::maxHCalCols, bins_t, min_tdct, max_tdct );
  TH2D *htimeid = new TH2D( "htimeid", "HCal-Hodo Time vs HCal ID; ID; ns", hcal::maxHCalChan, 0, hcal::maxHCalChan, bins_t, min_tdct, max_tdct );
  TH2D *htimex = new TH2D( "htimex", "HCal-Hodo Atime vs HCal X; Row; ns", 500, -3, 2, bins_t, min_tdct, max_tdct );
  TH2D *htimey = new TH2D( "htimey", "HCal-Hodo Atime vs HCal Y; Row; ns", 400, -2, 2, bins_t, min_tdct, max_tdct );

  TH2D *htimerun_bc = new TH2D( "htimerun_bc", "HCal-Hodo Time vs Run Num (intime cluster); Run; ns", nruns, 0, nruns, bins_t, min_tdct, max_tdct );
  TH2D *htimerow_bc = new TH2D( "htimerow_bc", "HCal-Hodo Time vs HCal Row (intime cluster); Row; ns", hcal::maxHCalRows, 0, hcal::maxHCalRows, bins_t, min_tdct, max_tdct );
  TH2D *htimecol_bc = new TH2D( "htimecol_bc", "HCal-Hodo Time vs HCal Col (intime cluster); Col; ns", hcal::maxHCalCols, 0, hcal::maxHCalCols, bins_t, min_tdct, max_tdct );
  TH2D *htimeid_bc = new TH2D( "htimeid_bc", "HCal-Hodo Time vs HCal ID (intime cluster); ID; ns", hcal::maxHCalChan, 0, hcal::maxHCalChan, bins_t, min_tdct, max_tdct );
  TH2D *htimex_bc = new TH2D( "htimex_bc", "HCal-Hodo Time vs HCal X (intime cluster); Row; ns", 500, -3, 2, bins_t, min_tdct, max_tdct );
  TH2D *htimey_bc = new TH2D( "htimey_bc", "HCal-Hodo Time vs HCal Y (intime cluster); Row; ns", 400, -2, 2, bins_t, min_tdct, max_tdct );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  TChain *TS = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // new output tree vars
  Double_t dx_out;
  Double_t dy_out;
  Double_t W2_out;
  Double_t Q2_out;
  Double_t atimediff_bc_out;
  Double_t timediff_bc_out;
  Int_t run_out;
  Int_t tar_out; //0:LH2, 1:LD2
  Int_t mag_out;
  Int_t isP_out;
  Int_t isN_out;
  Int_t isPbc_out;
  Int_t isNbc_out;
  
  // relevant old output tree vars
  Double_t bb_hodotdc_clus_tmean_out;
  Double_t sbs_hcal_clus_id_out[hcal::maxClus];
  Double_t sbs_hcal_clus_e_out[hcal::maxClus];
  Double_t sbs_hcal_clus_x_out[hcal::maxClus];
  Double_t sbs_hcal_clus_y_out[hcal::maxClus];
  Double_t sbs_hcal_clus_tdctime_out[hcal::maxClus];
  Double_t sbs_hcal_clus_atime_out[hcal::maxClus];
  Double_t sbs_hcal_nclus_out;
  Double_t sbs_hcal_nblk_out;
  Double_t sbs_hcal_clus_blk_id_out[hcal::maxBlk];
  Double_t sbs_hcal_clus_blk_e_out[hcal::maxBlk];
  Double_t sbs_hcal_clus_blk_x_out[hcal::maxBlk];
  Double_t sbs_hcal_clus_blk_y_out[hcal::maxBlk];
  Double_t sbs_hcal_clus_blk_atime_out[hcal::maxBlk];
  Double_t sbs_hcal_clus_blk_tdctime_out[hcal::maxBlk];
  Int_t Ndata_sbs_hcal_clus_blk_id_out;
  Int_t Ndata_sbs_hcal_clus_id_out;

  // set new output tree branches
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "atimediff_bc", &atimediff_bc_out, "atimediff_bc/D" );
  P->Branch( "timediff_bc", &timediff_bc_out, "timediff_bc/D" );
  P->Branch( "mag", &mag_out, "mag/I" );
  P->Branch( "run", &run_out, "run/I" );
  P->Branch( "tar", &tar_out, "tar/I" );
  P->Branch( "isP", &isP_out, "isP/I" );
  P->Branch( "isN", &isN_out, "isN/I" );
  P->Branch( "isPbc", &isPbc_out, "isPbc/I" );
  P->Branch( "isNbc", &isNbc_out, "isNbc/I" );

  // set relevant old output tree branches
  P->Branch( "bb.hodotdc.clus.tmean", &bb_hodotdc_clus_tmean_out, "bb.hodotdc.clus.tmean/D" );
  P->Branch( "sbs.hcal.clus.id", &sbs_hcal_clus_id_out, "sbs.hcal.clus.id/D" );
  P->Branch( "sbs.hcal.clus.e", &sbs_hcal_clus_e_out, "sbs.hcal.clus.e/D" );
  P->Branch( "sbs.hcal.clus.x", &sbs_hcal_clus_x_out, "sbs.hcal.clus.x/D" );
  P->Branch( "sbs.hcal.clus.y", &sbs_hcal_clus_y_out, "sbs.hcal.clus.y/D" );
  P->Branch( "sbs.hcal.clus.tdctime", &sbs_hcal_clus_tdctime_out, "sbs.hcal.clus.tdctime/D" );
  P->Branch( "sbs.hcal.clus.atime", &sbs_hcal_clus_atime_out, "sbs.hcal.clus.atime/D" );
  P->Branch( "sbs.hcal.clus_blk.id", &sbs_hcal_clus_blk_id_out, "sbs.hcal.clus_blk.id/D" );
  P->Branch( "sbs.hcal.clus_blk.e", &sbs_hcal_clus_blk_e_out, "sbs.hcal.clus_blk.e/D" );
  P->Branch( "sbs.hcal.clus_blk.x", &sbs_hcal_clus_blk_x_out, "sbs.hcal.clus_blk.x/D" );
  P->Branch( "sbs.hcal.clus_blk.y", &sbs_hcal_clus_blk_y_out, "sbs.hcal.clus_blk.y/D" );
  P->Branch( "sbs.hcal.clus_blk.tdctime", &sbs_hcal_clus_blk_tdctime_out, "sbs.hcal.clus_blk.tdctime/D" );
  P->Branch( "sbs.hcal.clus_blk.atime", &sbs_hcal_clus_blk_atime_out, "sbs.hcal.clus_blk.atime/D" );
  P->Branch( "sbs.hcal.nclus", &sbs_hcal_nclus_out, "sbs.hcal.nclus/D" );
  P->Branch( "sbs.hcal.nblk", &sbs_hcal_nblk_out, "sbs.hcal.nblk/D" );
  P->Branch( "Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id_out, "Ndata.sbs.hcal.clus.id/I" );
  P->Branch( "Ndata.sbs.hcal.clus_blk.id", &Ndata_sbs_hcal_clus_blk_id_out, "Ndata.sbs.hcal.clus_blk.id/I" );

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,kine);    
  cout << config_parameters << endl;
  Double_t hcaltheta = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbthr = config_parameters.GetBBtheta_rad(); //in radians

  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";
  vector<caldiag> allcut;

  vector<string> used_targets;
  vector<int> used_fields;

  for (Int_t irun=0; irun<nruns; irun++) {

    // only analyze hydrogen data to minimize better isolate elastic selection cuts
    std::string targ = runs[irun].target;
    Int_t t=-1;
    if( targ.compare("lh2")==0 )
      t=0;
    else if( targ.compare("ld2")==0 )
      t=1;
    else if( targ.compare("he3")==0 )
      t=2;
    else{
      cout << "ERROR: Bad target input. Check and rerun." << endl;
      return;
    }
    
    if( !all_targets && t>0 )
      continue;

    //Address run range
    Int_t runnum = runs[irun].runnum;
    if( !runall && (runnum < run_b || runnum > run_e) )
      continue;

    //Address exclusion range
    if( !noexclude && (runnum >= run_exclude_b && runnum <= run_exclude_e) )
      continue;

    std::string targ_uppercase = targ; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    
    cout << "Processing run " << runnum << ".." << endl;

    // get additional information from run object
    Int_t mag = runs[irun].sbsmag / 21; //convert to percent where max field is at 2100A
    Double_t ebeam = runs[irun].ebeam; //get beam energy per run
    Double_t charge = runs[irun].charge; //get total charge per run
    
    //Read in rootfile directory and get path
    std::string rootfile_dir = runs[irun].rootfile_dir;
    std::string rootfile_path = rootfile_dir + Form("/*%d*",runnum);

    //Get target parameters
    SBStarget target_parameters(targ);
    Int_t target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
    Double_t target_length = target_parameters.GetTargLength();
    Double_t target_rho = target_parameters.GetTargRho();
    Double_t cell_rho = target_parameters.GetCellRho();
    Double_t cell_diam = target_parameters.GetCellDiam();
    Double_t cell_dEdx = target_parameters.GetCelldEdx();
    Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
    Double_t target_dEdx = target_parameters.GetTargdEdx();
    Double_t M_avg = target_parameters.GetAvgMass();
    
    //Read in hcal vertical offset per run
    Double_t hcal_v_offset = runs[irun].hcal_v_offset;
    
    //Get available cuts for current config/target/field combination. Use first element (0) of cut
    vector<caldiag> cut;

    util::ReadDiagnosticCutList(struct_dir,experiment,kine,targ,mag,verb,cut);  
    
    cout << "Cut list loaded.." << endl;

    //cout << cut[0].gcut << endl;

    cout << cut[0] << endl;

    //search the used fields and targets, where current field/targ exists, don't add to allcut
    bool found_targ = false;
    if (std::find(used_targets.begin(), used_targets.end(), targ) != used_targets.end()) {
      found_targ = true; // Found the target
    }else{
      used_targets.push_back(targ); // Didn't find the target. Add to vector.
      cout << "Found new target, adding to caldiag allcut vector.." << endl;
    }

    bool found_field = false;
    if (std::find(used_fields.begin(), used_fields.end(), mag) != used_fields.end()) {
      found_field = true; // Found the field
    }else{
      used_fields.push_back(mag); // Didn't find the field. Add to vector.
      cout << "Found new SBS field setting, adding to caldiag allcut vector.." << endl;
    }

    //If either weren't found, new cut set exists. Add to allcut vector
    if (!found_field || !found_targ)
      allcut.push_back(cut[0]);
    
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
    
    C = new TChain("T");
    C->Add(rootfile_path.c_str());
    
    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);    
    C->SetMakeClass(1); //For event header branches

    // HCal general
    Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime, hcalidx, nclus, nblk;

    C->SetBranchStatus( "sbs.hcal.idblk", 1 );
    C->SetBranchStatus( "sbs.hcal.e", 1 );
    C->SetBranchStatus( "sbs.hcal.x", 1 );
    C->SetBranchStatus( "sbs.hcal.y", 1 );
    C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
    C->SetBranchStatus( "sbs.hcal.colblk", 1 );
    C->SetBranchStatus( "sbs.hcal.tdctimeblk", 1 );
    C->SetBranchStatus( "sbs.hcal.atimeblk", 1 );
    C->SetBranchStatus( "sbs.hcal.index", 1 );
    C->SetBranchStatus( "sbs.hcal.nclus", 1 );
    C->SetBranchStatus( "sbs.hcal.nblk", 1 );

    C->SetBranchAddress( "sbs.hcal.idblk", &hcalid );
    C->SetBranchAddress( "sbs.hcal.e", &hcale );
    C->SetBranchAddress( "sbs.hcal.x", &hcalx );
    C->SetBranchAddress( "sbs.hcal.y", &hcaly );
    C->SetBranchAddress( "sbs.hcal.rowblk", &hcalr );
    C->SetBranchAddress( "sbs.hcal.colblk", &hcalc );
    C->SetBranchAddress( "sbs.hcal.tdctimeblk", &hcaltdc );
    C->SetBranchAddress( "sbs.hcal.atimeblk", &hcalatime );
    C->SetBranchAddress( "sbs.hcal.index", &hcalidx );
    C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
    C->SetBranchAddress( "sbs.hcal.nblk", &nblk );

    // HCal cluster branches
    Double_t hcalcid[hcal::maxClus], hcalce[hcal::maxClus], hcalcx[hcal::maxClus], hcalcy[hcal::maxClus], hcalctdctime[hcal::maxClus], hcalcatime[hcal::maxClus], hcalcrow[hcal::maxClus], hcalccol[hcal::maxClus];
    Int_t Nhcalcid;

    C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.x", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.y", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.row", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.col", 1 );
    C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
    C->SetBranchStatus( "Ndata.sbs.hcal.clus.id", 1 );

    C->SetBranchAddress( "sbs.hcal.clus.id", hcalcid );
    C->SetBranchAddress( "sbs.hcal.clus.e", hcalce );
    C->SetBranchAddress( "sbs.hcal.clus.x", hcalcx );
    C->SetBranchAddress( "sbs.hcal.clus.y", hcalcy );
    C->SetBranchAddress( "sbs.hcal.clus.tdctime", hcalctdctime );
    C->SetBranchAddress( "sbs.hcal.clus.atime", hcalcatime );
    C->SetBranchAddress( "sbs.hcal.clus.row", hcalcrow );
    C->SetBranchAddress( "sbs.hcal.clus.col", hcalccol );
    C->SetBranchAddress( "sbs.hcal.clus.id", hcalcid );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Nhcalcid );

    // HCal cluster blk branches
    Double_t hcalcbid[hcal::maxClus], hcalcbe[hcal::maxClus], hcalcbx[hcal::maxClus], hcalcby[hcal::maxClus], hcalcbtdctime[hcal::maxClus], hcalcbatime[hcal::maxClus];
    Int_t Nhcalcbid;
    
    C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.x", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.y", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
    C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    C->SetBranchStatus( "Ndata.sbs.hcal.clus_blk.id", 1 );
   
    C->SetBranchAddress( "sbs.hcal.clus_blk.id", hcalcbid );
    C->SetBranchAddress( "sbs.hcal.clus_blk.e", hcalcbe );
    C->SetBranchAddress( "sbs.hcal.clus_blk.x", hcalcbx );
    C->SetBranchAddress( "sbs.hcal.clus_blk.y", hcalcby );
    C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", hcalcbtdctime );
    C->SetBranchAddress( "sbs.hcal.clus_blk.atime", hcalcbatime );
    C->SetBranchAddress( "Ndata.sbs.hcal.clus_blk.id", &Nhcalcbid );

    // bbcal clus var
    Double_t eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, nclusSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;

    C->SetBranchStatus( "bb.sh.e", 1 );
    C->SetBranchStatus( "bb.sh.x", 1 );
    C->SetBranchStatus( "bb.sh.y", 1 );
    C->SetBranchStatus( "bb.sh.rowblk", 1 );
    C->SetBranchStatus( "bb.sh.colblk", 1 );
    C->SetBranchStatus( "bb.sh.idblk", 1 );
    C->SetBranchStatus( "bb.sh.atimeblk", 1 );
    C->SetBranchStatus( "bb.sh.nclus", 1 );
    C->SetBranchStatus( "bb.ps.e", 1 );
    C->SetBranchStatus( "bb.ps.rowblk", 1 );
    C->SetBranchStatus( "bb.ps.colblk", 1 );
    C->SetBranchStatus( "bb.ps.idblk", 1 );
    C->SetBranchStatus( "bb.ps.atimeblk", 1 );

    C->SetBranchAddress( "bb.sh.e", &eSH );
    C->SetBranchAddress( "bb.sh.x", &xSH );
    C->SetBranchAddress( "bb.sh.y", &ySH );
    C->SetBranchAddress( "bb.sh.rowblk", &rblkSH );
    C->SetBranchAddress( "bb.sh.colblk", &cblkSH );
    C->SetBranchAddress( "bb.sh.idblk", &idblkSH );
    C->SetBranchAddress( "bb.sh.atimeblk", &atimeSH );
    C->SetBranchAddress( "bb.sh.nclus", &nclusSH );
    C->SetBranchAddress( "bb.ps.e", &ePS );
    C->SetBranchAddress( "bb.ps.rowblk", &rblkPS );
    C->SetBranchAddress( "bb.ps.colblk", &cblkPS );
    C->SetBranchAddress( "bb.ps.idblk", &idblkPS );
    C->SetBranchAddress( "bb.ps.atimeblk", &atimePS );

    // hodoscope cluster mean time
    Int_t Nhodotmean; 
    Double_t hodotmean[hcal::maxClus];
      
    C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
    C->SetBranchStatus( "Ndata.bb.hodotdc.clus.tmean", 1 );
      
    C->SetBranchAddress( "bb.hodotdc.clus.tmean", hodotmean );
    C->SetBranchAddress( "Ndata.bb.hodotdc.clus.tmean", &Nhodotmean );

    // track branches
    Double_t ntrack, p[hcal::maxTracks],px[hcal::maxTracks],py[hcal::maxTracks],pz[hcal::maxTracks],xtr[hcal::maxTracks],ytr[hcal::maxTracks],thtr[hcal::maxTracks],phtr[hcal::maxTracks];
    Double_t vx[hcal::maxTracks],vy[hcal::maxTracks],vz[hcal::maxTracks];
    Double_t xtgt[hcal::maxTracks],ytgt[hcal::maxTracks],thtgt[hcal::maxTracks],phtgt[hcal::maxTracks];

    C->SetBranchStatus( "bb.tr.n", 1 );
    C->SetBranchStatus( "bb.tr.p", 1 );
    C->SetBranchStatus( "bb.tr.px", 1 );
    C->SetBranchStatus( "bb.tr.py", 1 );
    C->SetBranchStatus( "bb.tr.pz", 1 );    
    C->SetBranchStatus( "bb.tr.x", 1 );    
    C->SetBranchStatus( "bb.tr.y", 1 );    
    C->SetBranchStatus( "bb.tr.th", 1 );    
    C->SetBranchStatus( "bb.tr.ph", 1 );    
    C->SetBranchStatus( "bb.tr.vx", 1 );
    C->SetBranchStatus( "bb.tr.vy", 1 );
    C->SetBranchStatus( "bb.tr.vz", 1 );
    C->SetBranchStatus( "bb.tr.tg_x", 1 );
    C->SetBranchStatus( "bb.tr.tg_y", 1 );
    C->SetBranchStatus( "bb.tr.tg_th", 1 );
    C->SetBranchStatus( "bb.tr.tg_ph", 1 );
    
    C->SetBranchAddress( "bb.tr.n", &ntrack );
    C->SetBranchAddress( "bb.tr.p", p );
    C->SetBranchAddress( "bb.tr.px", px );
    C->SetBranchAddress( "bb.tr.py", py );
    C->SetBranchAddress( "bb.tr.pz", pz );    
    C->SetBranchAddress( "bb.tr.x", xtr );    
    C->SetBranchAddress( "bb.tr.y", ytr );    
    C->SetBranchAddress( "bb.tr.th", thtr );    
    C->SetBranchAddress( "bb.tr.ph", phtr );    
    C->SetBranchAddress( "bb.tr.vx", vx );
    C->SetBranchAddress( "bb.tr.vy", vy );
    C->SetBranchAddress( "bb.tr.vz", vz );
    C->SetBranchAddress( "bb.tr.tg_x", xtgt );
    C->SetBranchAddress( "bb.tr.tg_y", ytgt );
    C->SetBranchAddress( "bb.tr.tg_th", thtgt );
    C->SetBranchAddress( "bb.tr.tg_ph", phtgt );

    // ekine branches
    Double_t ekineQ2, ekineW2, ekineeps, ekinenu, ekineqx, ekineqy, ekineqz;

    C->SetBranchStatus( "e.kine.Q2", 1 );
    C->SetBranchStatus( "e.kine.W2", 1 );
    C->SetBranchStatus( "e.kine.epsilon", 1 );
    C->SetBranchStatus( "e.kine.nu", 1 );
    C->SetBranchStatus( "e.kine.q_x", 1 );
    C->SetBranchStatus( "e.kine.q_y", 1 );
    C->SetBranchStatus( "e.kine.q_z", 1 );

    C->SetBranchAddress( "e.kine.Q2", &ekineQ2 );
    C->SetBranchAddress( "e.kine.W2", &ekineW2 );
    C->SetBranchAddress( "e.kine.epsilon", &ekineeps );
    C->SetBranchAddress( "e.kine.nu", &ekinenu );
    C->SetBranchAddress( "e.kine.q_x", &ekineqx );
    C->SetBranchAddress( "e.kine.q_y", &ekineqy );
    C->SetBranchAddress( "e.kine.q_z", &ekineqz );

    // fEvtHdr branches
    UInt_t rnum, gevnum, trigbits;

    C->SetBranchStatus( "fEvtHdr.fRun", 1 );
    C->SetBranchStatus( "fEvtHdr.fEvtNum", 1 );
    C->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );

    C->SetBranchAddress( "fEvtHdr.fRun", &rnum );
    C->SetBranchAddress( "fEvtHdr.fEvtNum", &gevnum );
    C->SetBranchAddress( "fEvtHdr.fTrigBits", &trigbits );

    // other bb branches
    Double_t gemNhits, eop;

    C->SetBranchStatus( "bb.gem.track.nhits", 1 );
    C->SetBranchStatus( "bb.etot_over_p", 1 );

    C->SetBranchAddress( "bb.gem.track.nhits", &gemNhits );
    C->SetBranchAddress( "bb.etot_over_p", &eop );

    // define very wide globalcut
    TCut GCut = gcut.c_str();

    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

    // get experimental quantities by run
    // set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta, hcalaxes );
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];
    Double_t BdL = hcal::sbsmaxfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = cell_diam/2.0/sin(bbthr) * target_rho * target_dEdx;

    // set nucleon default for LH2
    std::string nucleon = "p"; 

    // event indices
    long nevent = 0, nevents = C->GetEntries(); 
    Int_t treenum = 0, currenttreenum = 0;

    while (C->GetEntry(nevent++)) {
	
      std::cout << "Processing run " << runnum << "(" << irun << " / " << nruns << "), event " << nevent << " / " << nevents << "\r";
      std::cout.flush();
	
      ///////
      //Single-loop globalcut method. Save pass/fail for output tree.
      bool failedglobal = false;

      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
      }
      failedglobal = GlobalCut->EvalInstance(0) == 0;
      if( failedglobal )
	continue;

      ///////
      //Physics calculations
      //correct beam energy with vertex information
      Double_t ebeam_c = ebeam - ( (vz[0]+target_length/2.0) * target_rho * target_dEdx + upstream_wthick * cell_rho * cell_dEdx );

      TVector3 vertex( 0., 0., vz[0] );

      //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      Double_t precon = p[0] + Eloss_outgoing;

      //set up four-momenta with some empty for various calculation methods
      TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
      TLorentzVector pe( precon*px[0]/p[0], precon*py[0]/p[0], precon*pz[0]/p[0], precon ); //e' recon plvect
      TLorentzVector ptarg; ptarg.SetPxPyPzE( 0., 0., 0., M_avg ); //target momentum
      TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
      TVector3 qv = q.Vect();
      TLorentzVector pN; //N' momentum
      
      //simple calculations for e' and N'
      Double_t etheta = acos( pe.Pz() / pe.E() );
      Double_t ephi = atan2( pe.Py(), pe.Px() );
      Double_t pcent = ebeam_c/( 1. + ( ebeam_c/M_avg )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
      Double_t phNexp = ephi + hcal::PI;

      //e' p reconstruction with track angles (not momentum)
      Double_t nu = pbeam.E() - pcent;
      Double_t pNexp = sqrt( pow(nu, 2.) + 2. * M_avg * nu );
      Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
      TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
      pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
      Double_t Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
      Double_t W2 = pow( M_avg, 2.0 ) + 2.0*M_avg * (ebeam_c-pe.E()) - Q2;

      /////////////////////
      //W2 elastic cut
      bool failedW2 = W2>W2max || W2<W2min;
      if(failedW2)
	continue;

      //Calculate h-arm quantities
      vector<Double_t> xyhcalexp; util::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      TVector3 hcalpos = hcalorigin + hcalx*hcalaxes[0] + hcaly*hcalaxes[1]; //from primary blk
      Double_t KE_p = nu; //For elastics total predicted by earm
      Double_t SFrac = hcale/KE_p; //Measured
      Double_t dx = hcalx - xyhcalexp[0];
      Double_t dy = hcaly - xyhcalexp[1];
      TVector3 neutdir = (hcalpos - vertex).Unit();
      Double_t thetapq = acos( neutdir.Dot( qv.Unit() ) );
      Double_t eoverp = (ePS + eSH) / p[0];

      //Best Cluster analysis using intime algorithm for simplicity on diagnosis

      //Set up best cluster index and energy for sorting
      Double_t best_cluster[2] = {-1,0.};

      //loop through all clusters and select without HCal position information
      for( int c=0; c<Nhcalcid; c++ ){
	
	//calculate h-arm physics quantities per cluster
	double atime = hcalcatime[c];
	double atime_hodo = atime - hodotmean[0];
	double ce = hcalce[c];
	
	//wide cut around HCal/Hodo time
	bool passedTime = abs(atime_hodo-atime0)<intime_sig*atimesig;

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
      
      hdxdy->Fill(dy,dx);
      hdxdy_bc->Fill(dy_bestcluster,dx_bestcluster);

      //Cut on dxdy spots bools
      bool is_proton = util::Nspotcheck(dy, dx, dy0, dx0_p, spot_sig*dysig, spot_sig*dxsig_p);
      bool is_neutron = util::Nspotcheck(dy, dx, dy0, dx0_n, spot_sig*dysig, spot_sig*dxsig_n);
      bool is_proton_bc = util::Nspotcheck(dy_bestcluster, dx_bestcluster, dy0, dx0_p, spot_sig*dysig, spot_sig*dxsig_p);
      bool is_neutron_bc = util::Nspotcheck(dy_bestcluster, dx_bestcluster, dy0, dx0_n, spot_sig*dysig, spot_sig*dxsig_n);

      //Time quantities
      Double_t hadiff_bestcluster = hcalcatime[cidx_best] - hodotmean[0];
      Double_t hadiff = hcalcatime[0] - hodotmean[0];
      Double_t hdiff_bestcluster = hcalctdctime[cidx_best] - hodotmean[0];
      Double_t hdiff = hcalctdctime[0] - hodotmean[0];

      hatimerun_ecut->Fill( irun, hadiff);
      hatimerow_ecut->Fill( hcalcrow[0], hadiff);
      hatimecol_ecut->Fill( hcalccol[0], hadiff);
      hatimeid_ecut->Fill( hcalcid[0], hadiff);
      hatimex_ecut->Fill( hcalcx[0], hadiff);
      hatimey_ecut->Fill( hcalcy[0], hadiff);

      htimerun_ecut->Fill( irun, hdiff);
      htimerow_ecut->Fill( hcalcrow[0], hdiff);
      htimecol_ecut->Fill( hcalccol[0], hdiff);
      htimeid_ecut->Fill( hcalcid[0], hdiff);
      htimex_ecut->Fill( hcalcx[0], hdiff);
      htimey_ecut->Fill( hcalcy[0], hdiff);

      isP_out=0;
      if( is_proton )
	isP_out=1;

      isN_out=0;
      if( is_neutron )
	isN_out=1;

      if( is_proton || is_neutron ){
	hdxdy_cut->Fill(dy,dx);
	hatimerun->Fill( irun, hadiff);
	hatimerow->Fill( hcalcrow[0], hadiff);
	hatimecol->Fill( hcalccol[0], hadiff);
	hatimeid->Fill( hcalcid[0], hadiff);
	hatimex->Fill( hcalcx[0], hadiff);
	hatimey->Fill( hcalcy[0], hadiff);

	htimerun->Fill( irun, hdiff);
	htimerow->Fill( hcalcrow[0], hdiff);
	htimecol->Fill( hcalccol[0], hdiff);
	htimeid->Fill( hcalcid[0], hdiff);
	htimex->Fill( hcalcx[0], hdiff);
	htimey->Fill( hcalcy[0], hdiff);
      }      

      isPbc_out=0;
      if( is_proton_bc )
	isPbc_out=1;

      isNbc_out=0;
      if( is_neutron_bc )
	isNbc_out=1;

      if( is_proton_bc || is_neutron_bc ){
	hdxdy_cut_bc->Fill(dy_bestcluster,dx_bestcluster);
	hatimerun_bc->Fill( irun, hadiff_bestcluster);
	hatimerow_bc->Fill( row_bestcluster, hadiff_bestcluster);
	hatimecol_bc->Fill( col_bestcluster, hadiff_bestcluster);
	hatimeid_bc->Fill( id_bestcluster, hadiff_bestcluster);
	hatimex_bc->Fill( x_bestcluster, hadiff_bestcluster);
	hatimey_bc->Fill( y_bestcluster, hadiff_bestcluster);

	htimerun_bc->Fill( irun, hdiff_bestcluster);
	htimerow_bc->Fill( row_bestcluster, hdiff_bestcluster);
	htimecol_bc->Fill( col_bestcluster, hdiff_bestcluster);
	htimeid_bc->Fill( id_bestcluster, hdiff_bestcluster);
	htimex_bc->Fill( x_bestcluster, hdiff_bestcluster);
	htimey_bc->Fill( y_bestcluster, hdiff_bestcluster);
      }

      //Fill new output tree     
      dx_out = dx;
      dy_out = dy;
      W2_out = W2;
      Q2_out = Q2;
      atimediff_bc_out = hadiff_bestcluster;
      timediff_bc_out = hdiff_bestcluster;
      mag_out = mag;
      run_out = runnum;
      tar_out = t;

      //Fill old output tree
      bb_hodotdc_clus_tmean_out = hodotmean[0];
      sbs_hcal_nclus_out = nclus;
      Ndata_sbs_hcal_clus_id_out = Nhcalcid;
      for( Int_t c=0; c<Nhcalcid; ++c ){
	Double_t cid = hcalcid[c];
	Double_t ce = hcalce[c];
	Double_t cx = hcalcx[c];
	Double_t cy = hcalcy[c];
	Double_t ctdc = hcalctdctime[c];
	Double_t catime = hcalcatime[c];

	sbs_hcal_clus_id_out[c] = cid;
	sbs_hcal_clus_e_out[c] = ce;
	sbs_hcal_clus_x_out[c] = cx;	  
	sbs_hcal_clus_y_out[c] = cy;
	sbs_hcal_clus_tdctime_out[c] = ctdc;
	sbs_hcal_clus_atime_out[c] = catime;
      }
      sbs_hcal_nblk_out = nblk;
      Ndata_sbs_hcal_clus_blk_id_out = Nhcalcbid;
      for( Int_t b=0; b<Nhcalcbid; ++b ){
	Double_t cbid = hcalcbid[b];
	Double_t cbe = hcalcbe[b];
	Double_t cbx = hcalcbx[b];
	Double_t cby = hcalcby[b];
	Double_t cbtdc = hcalcbtdctime[b];
	Double_t cbatime = hcalcbatime[b];

	sbs_hcal_clus_blk_id_out[b] = cbid;
	sbs_hcal_clus_blk_e_out[b] = cbe;
	sbs_hcal_clus_blk_x_out[b] = cbx;	  
	sbs_hcal_clus_blk_y_out[b] = cby;
	sbs_hcal_clus_blk_tdctime_out[b] = cbtdc;
	sbs_hcal_clus_blk_atime_out[b] = cbatime;

      }

      P->Fill();

    }//end event loop

    // getting ready for the next run
    C->Reset();

  }//end run loop

  //Send data to report canvas file
  std::string report_path = out_dir + Form("/hcal_calibrations/pass%d/diagnostics/timingdiag_report_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.root",pass,experiment,kine,pass,run_b,run_e,run_exclude_b,run_exclude_e);

  //Create instance of struct to hold experimental parameters
  suppset supplement;
  std::string exper = experiment;
  std::string target;
  if(all_targets)
    target = "All Available";
  else
    target = "LH2 Only";

  supplement.exper = exper;
  supplement.kine = kine;
  supplement.pass = pass;
  supplement.date = date;
  supplement.runb = run_b;
  supplement.rune = run_e;
  supplement.runexb = run_exclude_b;
  supplement.runexe = run_exclude_e;
  supplement.targ = target;
  supplement.spotsig = spot_sig;

  fout->Write();

  //Use utility function to write all reports
  util::diagnosticReport(allcut,supplement,report_path);

  std::cout << std::endl << "Diagnostic complete. Output written to " << out_path << std::endl;
  std::cout << "Report file written to " << report_path << std::endl << std::endl;

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}
