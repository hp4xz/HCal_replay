//SSeeds 9.1.23 Script to obtain adc gain coefficients (pC->GeV) from small data set for use during startup with new HV settings.

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"

const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxBlk = 20; //More than the limit on possible blocks in a cluster in hcal
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const double atimeNsig = 3; //number of sigma to cut around hcal-bbcal coin peak

double hcalheight = -0.2897; //m The height of the center of HCAL above beam. Database should have this taken care of
//double hcalheight = 0.;
double sampFrac = 0.10; //HCal sampling frac (rough, mc)

// LH2
const double lh2tarrho = 0.0723;     //g/cc, target density
const double lh2cthick = 0.02;       //cm, target cell thickness
const double lh2uwallthick = 0.0145; //cm, upstream wall thickness
const double lh2dwallthick = 0.015;  //cm, downstream wall thickness
const double lh2dEdx = 0.00574; 

const double PI = TMath::Pi();
const double M_e = 0.00051;
const double M_p = 0.938272;

// Aluminum
const double rho_al = 2.7; //g/cc
const double aldEdx = 0.0021;

//all cryotarg
const double l_tgt = 0.15;

//double oldGain[kNcell] = {0.0};
//double oldRatio[kNcell] = {0.0};

double coeff[kNcell] = {0.0};
double gOldConst[kNcell] = {0.0};

// Get today's date
string getDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);
  
  string yyyy = to_string(1900 + ltm.tm_year);
  string mm = to_string(1 + ltm.tm_mon);
  string dd = to_string(ltm.tm_mday);
  string date = mm + '_' + dd + '_' + yyyy;
  
  return date;
}

void ecal_startup( const char *configfilename = "setup_ecal_startup.cfg", int run = -1 ){
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();
  
  // Declare Chain for many root files
  TChain *C = new TChain("T");

  string constPath = "coeff.txt";
  //string ratioPath = "ratio.txt";
  //string oneBlockPath = "oneBlock.txt";

  // Declare general physics parameters to be modified by input config file
  int uni_N = 400; // Total number of bins used to measure detection uniformity (hSampFrac histos)
  const int xN = 48; //2*kNrows, total number of dispersive bins detection uni
  const int yN = 24; //2*kNcols, total number of transverse bins detection uni
  double E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  int minEventPerCell = 50; // Minimum number of scattered p in cell required to calibrate, try for 200
  int maxEventPerCell = 4000; // Maximum number of scattered p events to contribute
  double highDelta = 0.1; // Minimum M(i,j)/b(i) factor allowed 
  double HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  double HCal_th = 35.0; // Angle that the center of HCal is at  
  double opticsCorr = 1.0; //added in case optics need correction for calibration
  double HCal_div = 0.15; // Transverse width in x and y per cell
  double HCal_Xi = -2.655; // Distance from beam center to top of HCal in analysis coordinates (m)
  double HCal_Xf = 1.155; // Distance from beam center to bottom of HCal in analysis coordinates (m)
  double HCal_Yi = -0.92964; // Distance from beam center to opposite-beam side of HCal in analysis coordinates (m)
  double HCal_Yf = 0.92964; // Distance from beam center to beam side of HCal in analysis coordinates (m)
  double W2_mean = 0.93; // Mean of W2 at current kinematic
  double W2_sig = 0.10; // Width of W2 at current kinematic
  double dx0 = 0.; // Position of proton spot, x-x_expected
  double dy0 = 0.; // Position of proton spot, y-y_expected
  double dx_sig = 0.10; // Max spread of proton spot, x-x_expected
  double dy_sig = 0.10; // Max spread of proton spot, y-y_expected
  double atime0 = 50.; //peak location of hcal-bbcal atime
  double atime_sig = 5.; //sigma of hcal-bbcal atime
  double ScaleFac = 1.0; //not used

  int elasYield = 0; // Keep track of total elastics analyzed

  //For position reconstruction
  // double HCal_Xmin = HCal_Xi-HCal_div/2;
  // double HCal_Xmax = HCal_Xf+HCal_div/2;
  // double HCal_Ymin = HCal_Yi-HCal_div/2;
  // double HCal_Ymax = HCal_Yf-HCal_div/2;

  // Reading config file
  cout << "Reading config file at " << configfilename << endl;

  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
      cout << "Loading file: " << currentline << ".." << endl;
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    int ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Loading beam energy: " << E_e << endl;
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "Loading HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof() * TMath::DegToRad();	
	cout << "Loading HCal angle: " << HCal_th << endl;
      }
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_mean = sval.Atof();
	cout << "Loading W2 mean cut: " << W2_mean << endl;
      }
      if( skey == "W2_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_sig = sval.Atof();
	cout << "Loading W2 sigma cut: " << W2_sig << endl;
      }
      if( skey == "dx0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0 = sval.Atof();
	cout << "Loading x position of spot: " << dx0 << endl;
      }
      if( skey == "dy0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0 = sval.Atof();
	cout << "Loading y position of spot: " << dy0 << endl;
      }
      if( skey == "dx_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig = sval.Atof();
	cout << "Loading x sigma of spot: " << dx_sig << endl;
      }
      if( skey == "dy_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig = sval.Atof();
	cout << "Loading y sigma of spot: " << dy_sig << endl;
      }
      if( skey == "atime0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	atime0 = sval.Atof();
	cout << "Loading hcal-bbcal adct peak location: " << atime0 << endl;
      }
      if( skey == "dy_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	atime_sig = sval.Atof();
	cout << "Loading hcal-bbcal adct sigma: " << atime_sig << endl;
      }
    }
    delete tokens;
  }
  
  cout << endl << endl << "Setup parameters loaded." << endl;

  // Declare general detector parameters
  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vz[maxTracks], BBtr_chi2[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, BBsh_atime, BBsh_nclus;

  double TDCT_id[maxTdcChan], TDCT_tdc[maxTdcChan]; 
  int TDCTndata;

  double HCALx, HCALy, HCALe;
  double cblkagain[maxBlk], cblkatime[maxBlk];
  double crow, ccol, nblk, nclus;
  double cblkid[kNcell], cblke[kNcell], cblkrow[kNcell], cblkcol[kNcell];


  // Declare root tree variables and set values to memory locations in root file
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
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );

  C->SetBranchStatus( "bb.tr.chi2", 1 );
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
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.sh.atimeblk", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
  C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus ); // Total number of blocks in highest E clus
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // cluster block adc time
  C->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain ); // cluster block adc gain params
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", cblkrow ); // cluster block row
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", cblkcol ); // cluster block column
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.chi2", BBtr_chi2 );
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
  C->SetBranchAddress( "bb.sh.nclus", &BBsh_nclus );
  C->SetBranchAddress( "bb.sh.atimeblk", &BBsh_atime );
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );  


  TCut GCut = globalcut;

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );


  // TEventList *elist = new TEventList("elist","Elastic Event List");
  // C->Draw(">>elist",globalcut);
  
  // cout << "Event list populated with cut placed on elastics." << endl;
  
  // Declare outfile
  TFile *fout = new TFile( "ecal_startup_out.root", "RECREATE" );
  
  // Initialize vectors and arrays
  double gRatio[kNcell] = {0.0}; //Old ratios from first iteration
  double oneBlock[kNcell] = {0.0};
  double GCoeff[kNcell] = {0.0};
  double GCoeff_divide[kNcell] = {0.0};
  //double A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Memory is reset on each event.

  // Initialize histograms	
  TH1D *hHCALe = new TH1D( "hHCALe","HCal Cluster E", 400, 0., 1. );
  TH1D *hSampFrac = new TH1D( "hSampFrac","HCal Cluster E / Expected KE", 400, 0., 1. );
  TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *hKE_p = new TH1D( "hkE_p", "proton kinetic energy; GeV", 500, 0.0, E_e*1.5 );
  TH1D *hKE_exp = new TH1D( "hkE_exp", "proton kinetic energy * sampling fraction; GeV", 500, 0.0, E_e*1.5 );



  TH1D *hDeltaE = new TH1D( "hDeltaE","1.0-Eclus/p_rec", 100, -1.5, 1.5 );
  TH1D *hEDiffChan = new TH1D( "hEDiffChan","EDiff Events Per Channel", 288, 0, 288 );
  TH1D *hClusE = new TH1D( "hClusE","Best Cluster Energy", 100, 0.0, 2.0);
  TH2D *hClusBlk_ClusE = new TH2D( "hClusBlk_ClusE","Cluster Size vs Recon Cluster E", 8, 0, 8, 100, 2.0, 3.0 );
  TH2D *hPAngleCorr = new TH2D( "hPAngCorr","Track p vs Track ang", 100, 30, 60, 100, 0.4, 1.2 );
  TH2D *hPAngleCorr_2 = new TH2D( "hPAngCorr_2","Track p vs Track ang v2", 100, 30, 60, 100, 0.4, 1.2 );
  TH2D *hClusE_vs_X = new TH2D("hClusE_vs_X",";X-pos (m);E_dep (GeV)",500,-3,2,100,0.0,1.0);  
  TH2D *hClusE_vs_Y = new TH2D("hClusE_vs_Y",";Y-pos (m);E_dep (GeV)",200,-1,1,100,0.0,1.0);  
  TH1D *hept = new TH1D( "hept", "Targ P", 600, 0, 6 );
  TH1D *hpp = new TH1D( "hpp", "Elastic Proton Momentum", 600, 0, 6 );
  TH1D *hdx = new TH1D( "hdx", "HCal X - X Expected; m", 400, -2, 2 );
  TH1D *hdy = new TH1D( "hdy", "HCal Y - Y Expected; m", 400, -2, 2 );
  TH1D *hW2 = new TH1D( "hW2", "W2 wide cut; GeV^2", 250, 0.3, 1.5 );
  TH1D *hNBlk = new TH1D( "hNBlk", "Number of Blocks in Primary Cluster", 25, 0, 25 );
  TH1D *hW2_cuts = new TH1D( "hW2_cuts", "W2 tight cuts; GeV^2", 250, 0.3, 1.5 );
  TH1D *hQ2 = new TH1D( "hQ2", "Q2; GeV^2", 250, 0.5, 3.0 );
  TH1D *hE_ep = new TH1D( "hE_ep","e' energy;GeV", 500, 0.0, E_e*1.5 ); 
  TH1D *hE_pp = new TH1D( "hE_pp", "proton energy; GeV", 500, 0.0, E_e*1.5 );
  TH2D *hADC = new TH2D( "hADC", "HCal Int_ADC Spectra: W Cut", 288, 0, 288, 100., 0., 1. );
  TH2D *hADC_amp = new TH2D( "hADC_amp", "HCal ADC_amp Spectra: W Cut", 288, 0, 288, 100., 0., 10. );
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hXY = new TH2D("hXY",";y_{expect} (m); x_{expect} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH2D *hXY_HCAL = new TH2D("hXY_HCAL",";y_{HCAL} (m); x_{HCAL} (m)", 250, -5.0, 5.0, 250, -10, 10 );
  TH1D *hvz = new TH1D("hvz",";vertex z (m);", 250,-0.125,0.125);
  TH1D *hvz_cut = new TH1D("hvz_cut",";vertex z (m);", 250,-0.125,0.125);
  TH2D *coefficients = new TH2D("coefficients",";Channel ;GeV/mV", 288,0,288,250,0,0.025);
  TH2D *coefficients_1b = new TH2D("coefficients_1b",";Channel ;GeV/mV", 288,0,288,250,0,0.025);
  TH2D *hEDiff_vs_X = new TH2D("hEDiff_vs_X",";X-pos (m);(E_exp-E_dep)/E_exp (GeV)",xN,HCal_Xi,HCal_Xf,100,-3,3);  
  TH2D *hEDiff_vs_Y = new TH2D("hEDiff_vs_Y",";Y-pos (m);(E_exp-E_dep)/E_exp (GeV)",yN,HCal_Yi,HCal_Yf,100,-3,3);
  TH2D *hEDiff_vs_block = new TH2D("hEDiff_vs_block",";Block (cell);(E_exp-E_dep)/E_exp (GeV)",288,0,288,100,-3,3);
  TH2D *htDiff_vs_ADCint[kNcell+1];
  TH2D *hSampFrac_vs_X = new TH2D( "hSampFrac_vs_X","Sampling Fraction ;X (m) ;HCal_E / Exp_KE", xN, HCal_Xi, HCal_Xf, uni_N, 0., 1. );
  TH2D *hSampFrac_vs_Y = new TH2D( "hSampFrac_vs_Y","Sampling Fraction ;Y (m) ;HCal_E / Exp_KE", yN, HCal_Yi, HCal_Yf, uni_N, 0., 1. );
  TH1D *hE_pp_corr = new TH1D( "Deposited Elastic Proton Energy", "E_pp_corr", 500, 0.0, E_e*1.5 );

  // Set long int to keep track of total entries
  Long64_t nevent = 0, Nevents = C->GetEntries();
  int treenum = 0, currenttreenum = 0;

  cout << endl << "All parameters loaded." << endl << endl;
  cout << "Opened up tree with nentries: " << Nevents << ".." << endl << endl;

  //Set up HCal coordinate system
  vector<TVector3> hcalaxes;
  TVector3 hcal_zaxis( sin(-HCal_th), 0, cos(-HCal_th) ); // Clock-wise rotation about Y axis
  TVector3 hcal_xaxis( 0, -1, 0 ); // -Y axis of Hall coordinate system = X axis of hcal coordinate system
  TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
  hcalaxes.push_back(hcal_xaxis);
  hcalaxes.push_back(hcal_yaxis);
  hcalaxes.push_back(hcal_zaxis);

  TVector3 hcalorigin = HCal_d*hcalaxes[2] + hcalheight*hcalaxes[0];

  //Declare matrices for chi-square min calibration scheme and keep track of calibrated events
  TMatrixD Ma(kNcell,kNcell);
  TVectorD ba(kNcell);
  TVectorD bb(kNcell);
  int NEV[kNcell] = {0};
  
  //Loop over events
  cout << "Main loop over all data commencing.." << endl;
  while (C->GetEntry(nevent++)) {
    cout << "Analyzing event " <<  nevent << "/" << Nevents << " \r";
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
    //HCal Active Area Cut
    int pblkrow = (int)cblkrow[0];
    int pblkcol = (int)cblkcol[0];
    int pblkid = (int)cblkid[0]-1;
    gOldConst[pblkid] = cblkagain[0];
    bool failedactivearea = 
      pblkrow==0 || 
      pblkrow==23 || 
      pblkcol==0 || 
      pblkcol==11;

    if( failedactivearea ) 
      continue; //All events with primary cluster element on edge blocks cut

    ///////
    //HCal primary cluster coincidence time cut (using adctime while hcal tdc suspect, new offsets)
    double atimediff = cblkatime[0]-BBsh_atime;
    bool failedcoin = abs(atimediff-atime0)>atimeNsig*atime_sig;
    hDiff->Fill(atimediff);

    if( failedcoin ) 
      continue; //All events where adctime outside of reasonable window cut

    double A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev

    double ebeam_c = E_e - ( (BBtr_vz[0]+l_tgt/2.0) * lh2tarrho * lh2dEdx + lh2uwallthick * rho_al * aldEdx );


    //vertex position
    TVector3 vertex( 0., 0., BBtr_vz[0] );

    //track momentum
    double trackp = BBtr_p[0];

    //set up four-momenta and physics variables per run
    TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
    //TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' plvect
    TLorentzVector pe( trackp*BBtr_px[0]/BBtr_p[0], 
		       trackp*BBtr_py[0]/BBtr_p[0], 
		       trackp*BBtr_pz[0]/BBtr_p[0], 
		       trackp ); //e' recon plvect
    TLorentzVector ptarg;
    ptarg.SetPxPyPzE( 0., 0., 0., M_p );

    TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
    TVector3 qv = q.Vect();
    double etheta = acos( pe.Pz() / pe.E() ); 
    double ephi = atan2( pe.Py(), pe.Px() );
    double pcent = E_e/( 1. + ( E_e/M_p )*( 1.0 - cos(etheta) ) );
    double phNexp = ephi + PI;
    TLorentzVector pN; //N' momentum
    double Q2, W2, nu, thNexp, pNexp;

    //Calculate using reconstructed track angles.
    nu = pbeam.E() - pcent;
    pNexp = sqrt( pow(nu, 2.) + 2. * M_p * nu );
    thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
    TVector3 pNhat( sin(thNexp) * cos(phNexp), 
		    sin(thNexp) * sin(phNexp), 
		    cos(thNexp) );
    pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
    Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
    W2 = pow( M_p, 2.0 ) + 2.0*M_p * (pbeam.E()-pe.E()) - Q2;

    hW2->Fill(W2);
    hQ2->Fill(Q2);

    /////////////////////
    //W2 elastic cut bool
    bool failedW2 = abs(W2-W2_mean)>W2_sig;

    if(failedW2)
      continue;

    //HCal active area cut (acceptance matching). Save pass/fail for output tree.
    vector<double> xyhcalexp; 
    double sintersect = ( hcalorigin - vertex).Dot(hcalaxes[2] ) / ( pNhat.Dot(hcalaxes[2]) );
    TVector3 hcal_intersect = vertex + sintersect*pNhat; 
    double xexpect_hcal = ( hcal_intersect - hcalorigin ).Dot( hcalaxes[0] );
    double yexpect_hcal = ( hcal_intersect - hcalorigin ).Dot( hcalaxes[1] );
    xyhcalexp.push_back( xexpect_hcal );
    xyhcalexp.push_back( yexpect_hcal );

    double dx = HCALx - xyhcalexp[0];
    double dy = HCALy - xyhcalexp[1];


    //Plot 2D histo of position
    hXY_HCAL->Fill( HCALy, HCALx );

    hXY->Fill( xyhcalexp[1], xyhcalexp[0] );

    double E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy
    hE_ep->Fill( E_ep ); // Fill histogram
	
	
    //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
    double E_pp = nu+M_p; // Get energy of the proton
    double Enucleon = sqrt(pow(pNexp,2)+pow(M_p,2)); // Check on E_pp, same
    hE_pp->Fill( E_pp ); // Fill histogram

    double KE_p = nu; //For elastics
    double KE_exp = KE_p*sampFrac; //Expected E in HCal, equivalent to KE*modified_samp_frac

    //Fill some histograms
    hKE_p->Fill( KE_exp );   
    hpp->Fill( pNexp );
    hdx->Fill( dx );
    hdy->Fill( dy );
    hdxdy_HCAL->Fill( dy, dx );
      
    ///////
    //dy elastic cut
    bool faileddy = abs(dy-dy0)>3*dy_sig;

    if( faileddy ) 
      continue;

    //Check energy deposited and uniformity in HCal after elastics cuts
    double SFrac = HCALe/KE_p;
    hHCALe->Fill( HCALe );
    hSampFrac->Fill( SFrac );
    hSampFrac_vs_X->Fill( HCALx, SFrac );
    hSampFrac_vs_Y->Fill( HCALy, SFrac );
    
    //Fill vertex position histogram for cut on tracks
    hvz_cut->Fill( BBtr_vz[0] );

    //Events that pass the above cuts constitute elastics
    elasYield++;
	
    double clusE = 0.0;
    double cluspC = 0.0;
	
    hNBlk->Fill( nblk );
	
    // Get energies with simplest scheme from clusters only
    for( int blk = 0; blk<(int)nblk; blk++ ){
      int blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
      double blke = cblke[blk];
      Double_t blkpC = cblke[blk]/cblkagain[blk]; //divide off the old coefficient
      
      clusE += blke;
      cluspC += blkpC;
      A[blkid] += blkpC;
	
      NEV[blkid]++;
	    
    }

    //Build the 288x288 matrix as clearly as possible
    for(Int_t icol = 0; icol<kNcell; icol++){ //matrix column here
      ba(icol) += A[icol];

      for(Int_t irow = 0; irow<kNcell; irow++){ //matrix row here
	Ma(icol,irow) += A[icol]*A[irow]/KE_exp;
      } //endloop over matrix element rows
    } //endloop over matrix element cols

  } //endloop over events


  cout << endl << "Checking data, inverting matrix, and solving for coefficients.." << endl << endl;

  //Reject the bad cells and normalize the oneblock check
  int badcell[kNcell];
  int cellBad = 0;
  double y[kNcell] = {0.0}; // For easy TGraphErrors build
  
  for(int i=0; i<kNcell; i++){
    badcell[i] = 0;
    y[i] = i;
    
    //Do not change ADC gain coeff if insufficient events or energy dep in cell
    if( NEV[i] < minEventPerCell || Ma(i,i) < 0.1*ba(i) ){ 

      cellBad = 1;

      double elemRatio = Ma(i,i)/ba(i);

      ba(i) = 1.0;  // Set RHS vector for cell i to 1.0 
      Ma(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
      for(int j=0; j<kNcell; j++){
	if( j != i ){
	  
	  Ma(i,j) = 0.0; 
	  Ma(j,i) = 0.0;
	}
      }
      badcell[i] = 1;
      cout << "cell" << i << " bad" << endl;
      cout << "Number of events in bad cell =" << NEV[i] << endl;
      cout << "Matrix element/vector element ratio =" << elemRatio << endl;
    }  
  }

  if( cellBad==0 ) cout << "No bad cells detected.." << endl << endl;
  
  //Invert the matrix, solve for ratios
  TMatrixD M_inv = Ma.Invert();
  TVectorD Coeff = M_inv*ba; // Stays unmodified for reference

  for(int i=0; i<kNcell; i++){
    if(badcell[i]==0){
      GCoeff[i]=Coeff[i]; // The new gain coefficient is the old coefficient multiplied by the gain factor that we just solved for. We will take these coefficents as the input for another set of data to estimate the error per channel.
    }else{
      GCoeff[i]=gOldConst[i]; // If the cell is bad, use the old coefficient
      GCoeff_divide[i]=-1.0;
    }
  }

  cout << "Inversion complete. Building histograms and writing out coefficients.." << endl << endl;

  for(int i=0; i<kNcell; i++){
    coefficients->Fill(i,GCoeff[i]);
  }

  double yErr[kNcell] = {0.0};
  double constErr[kNcell] = {0.0}; //Will need to improve error here

  TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, GCoeff, yErr, constErr ); 
  ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph->GetYaxis()->SetLimits(0.0,0.25);
  ccgraph->SetTitle("Calibration Coefficients");
  ccgraph->GetXaxis()->SetTitle("Channel");
  ccgraph->GetXaxis()->SetTitle("Unitless");
  ccgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  ccgraph->Write("constants");

  TGraphErrors *ccgraph_divide = new TGraphErrors( kNcell, y, GCoeff_divide, yErr, yErr ); 
  ccgraph_divide->GetXaxis()->SetLimits(0.0,288.0);  
  ccgraph_divide->SetTitle("Calibration Coefficients / OneBlock Coeff");
  ccgraph_divide->GetXaxis()->SetTitle("Channel");
  ccgraph_divide->GetXaxis()->SetTitle("Unitless");
  ccgraph_divide->SetMarkerStyle(21);
  ccgraph_divide->Write("constants_divide");

  //Write out diagnostic histos and print to console
  //fout->Write();

  cout << "Gain Coefficients: " << endl << endl;

  int cell = 0;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      //if( GCoeff[cell]==0 ) GCoeff[cell] = gOldConst[cell];
      cout << GCoeff[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  fout->Write();

  cout << endl << "Number of events available for calibration: " << endl << endl;
  
  cell = 0;

  for( int r = 0; r<kNrows; r++){
    for( int c = 0; c<kNcols; c++){
      cout << NEV[cell] << "  ";
      cell++;
    }
    cout << endl;
  }

  cout << endl << endl;

  //Declare outfile
  ofstream cfile;

  //Write to outfiles
  cfile.open( constPath );
  cfile << "#HCal gain coeffiencts obtained " << date.c_str() << endl;
  cfile << "sbs.hcal.adc.gain = " << endl;

  cell = 0;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){
      cfile << GCoeff[cell] << "  ";
      cell++;
    }
    cfile << endl;
  }

  cfile.close();

  cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total events analyzed: " << Nevents << "." << endl << endl;
  

  cout << "Calibration complete and constants written to file at " << constPath << "." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}

//Sample setup file, comment with #
////////////////////////////////////////////////////////
//setup_HCal_Calibration.txt
////////////////////////////////////////////////////////
// #LH2
// #/adaqfs/home/a-onl/sbs/Rootfiles/gmn_replayed_11263*
// endlist
// bb.tr.n==1&&bb.ps.e>0.2&&sbs.hcal.nclus>0&&bb.sh.nclus>0
// endcut
// E_e 1.0
// HCal_d 1.0
// HCal_th 1.0
// opticsCorr 1.0
// W_mean 1.0
// W_sig 0.1
// dx0 0.
// dy0 0.
// dx_sig 0.1
// dy_sig 0.1
// atime0 50.
// atime_sig 5.
// ScaleFac 1.0
/////////////////////////////////////////////////////////
