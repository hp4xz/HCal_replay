//SSeeds - 9.1.23 Startup energy calibration script configured for ease of use.

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <unistd.h>
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
#include "TString.h"
#include "TF1.h"

//Detector constants
const Int_t ifac = 3; // Inclusion factor, number of sigma to keep around cut peaks
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const Int_t uni_N = 400; // Total number of bins used to measure detection uniformity (hSampFrac histos)
const Int_t xN = 48; //2*kNrows, total number of dispersive bins detection uni
const Int_t yN = 24; //2*kNcols, total number of transverse bins detection uni
const Double_t hcalheight = 0.;
const Double_t sampFrac = 0.10;
const Double_t dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const Double_t dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const Double_t uwallthick_LH2 = 0.0145; //cm
const Double_t dwallthick_LH2 = 0.015; //cm
const Double_t cellthick_LH2 = 0.02; //cm, this is a guess;
const Double_t Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch  
const Double_t dxlim_l = 4.0;
const Double_t dxlim_u = 3.0;
const Double_t l_tgt = 0.15; // Length of the target (m)
const Double_t rho_tgt = 0.0723; // Density of target (g/cc)
const Double_t rho_Al = 2.7; // Density of aluminum windows (g/cc)
const Double_t celldiameter = 1.6*2.54; //cm
const Double_t Ztgt = 1.0;
const Double_t Atgt = 1.0;
const Double_t Mmol_tgt = 1.008; //g/mol
const UInt_t second = 1000000; 

//Math
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;

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

//Main.
void ecal( const char *config_path, Int_t iter=0 ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = getDate();

  //Set up bool to track quasi-replay
  bool qreplay = iter==1;

  // Declare file name paths
  string constPath;
  if( qreplay ){
    constPath = "coeff_qreplay.txt";
  }else{
    constPath = "coeff.txt";
  }

  //Get gain parameters from database for pass0 where sbs.hcal.again not on tree
  Double_t gConst_iter1[kNcell] = {0.}; 
  //Set up with all coefficients 1 have no effect on iter == 0
  for( Int_t c=0; c<kNcell; c++){
    gConst_iter1[c] = 1.;
  }

  if( qreplay ){
    string inConstPath = "coeff.txt";

    cout << "Loading new gain coefficents for quasi-replay from: " << inConstPath << ".." << endl;
    ifstream inConstFile( inConstPath );
    if( !inConstFile ){
      cerr << endl << "ERROR: No input constant file present -> path to coeff_sbs<kine>.txt expected." << endl;
      return 0;
    }

    Int_t n1=0;
    Double_t d1=0;
    string line;
    bool skip_line = true;

    while( getline( inConstFile, line ) ){
      
      if( n1==( kNcell ) ) break;
      
      TString Tline = (TString)line;
      
      if( Tline.BeginsWith("sbs.hcal.adc.gain") && skip_line==true ){	
	skip_line = false;
	continue;
      }

      if( skip_line==false ){
	istringstream iss( line );
	while( iss >> d1 ){
	  gConst_iter1[n1] = d1;
	  n1++;
	}
      }
    }
    
    cout << endl << endl << "ADC gain params from iteration 0 calibration: " << endl;
    
    for( Int_t r=0; r<kNrows; r++){
      for( Int_t c=0; c<kNcols; c++){
	Int_t i = r*kNcols+c;
	cout << gConst_iter1[i] << " ";
      }
      cout << endl;
    }

    usleep( 3*second ); //Give some time for review of step

  } 

  cout << endl;

  TChain *C = new TChain("T");


  //Set up multiple chains and elists for hydrogen and deuterium
  TChain *Ch[nfset_lh2[kIdx]];
  TChain *Cd[nfset_ld2[kIdx]];

  TEventList *elist_h[nfset_lh2[kIdx]];
  Long64_t Nevents_h[nfset_lh2[kIdx]];
  Long64_t NTevents_h[nfset_lh2[kIdx]];

  TEventList *elist_d[nfset_ld2[kIdx]];
  Long64_t Nevents_d[nfset_ld2[kIdx]];
  Long64_t NTevents_d[nfset_ld2[kIdx]];

  TCut globalcut_h[nfset_lh2[kIdx]];
  TCut globalcut_d[nfset_ld2[kIdx]];

  // Declare general detector parameters - lh2
  Double_t BBtr_p_h[nfset_lh2[kIdx]][maxTracks], BBtr_px_h[nfset_lh2[kIdx]][maxTracks], BBtr_py_h[nfset_lh2[kIdx]][maxTracks], BBtr_pz_h[nfset_lh2[kIdx]][maxTracks];
  Double_t BBtr_vz_h[nfset_lh2[kIdx]][maxTracks], BBtr_chi2_h[nfset_lh2[kIdx]][maxTracks];
  Double_t BBtr_n_h[nfset_lh2[kIdx]], BBps_x_h[nfset_lh2[kIdx]], BBps_y_h[nfset_lh2[kIdx]], BBps_e_h[nfset_lh2[kIdx]], BBsh_x_h[nfset_lh2[kIdx]], BBsh_y_h[nfset_lh2[kIdx]], BBsh_e_h[nfset_lh2[kIdx]];

  Double_t HCALx_h[nfset_lh2[kIdx]], HCALy_h[nfset_lh2[kIdx]], HCALe_h[nfset_lh2[kIdx]];
  Double_t crow_h[nfset_lh2[kIdx]], ccol_h[nfset_lh2[kIdx]], nblk_h[nfset_lh2[kIdx]];
  Double_t cblkid_h[nfset_lh2[kIdx]][kNcell], cblke_h[nfset_lh2[kIdx]][kNcell], cblkatime_h[nfset_lh2[kIdx]][kNcell];
  Double_t cblkagain_h[nfset_lh2[kIdx]][kNcell];
  Double_t ekineW2_h[nfset_lh2[kIdx]];

  // Declare general detector parameters - ld2
  Double_t BBtr_p_d[nfset_ld2[kIdx]][maxTracks], BBtr_px_d[nfset_ld2[kIdx]][maxTracks], BBtr_py_d[nfset_ld2[kIdx]][maxTracks], BBtr_pz_d[nfset_ld2[kIdx]][maxTracks];
  Double_t BBtr_vz_d[nfset_ld2[kIdx]][maxTracks], BBtr_chi2_d[nfset_ld2[kIdx]][maxTracks];
  Double_t BBtr_n_d[nfset_ld2[kIdx]], BBps_x_d[nfset_ld2[kIdx]], BBps_y_d[nfset_ld2[kIdx]], BBps_e_d[nfset_ld2[kIdx]], BBsh_x_d[nfset_ld2[kIdx]], BBsh_y_d[nfset_ld2[kIdx]], BBsh_e_d[nfset_ld2[kIdx]];

  Double_t HCALx_d[nfset_ld2[kIdx]], HCALy_d[nfset_ld2[kIdx]], HCALe_d[nfset_ld2[kIdx]];
  Double_t crow_d[nfset_ld2[kIdx]], ccol_d[nfset_ld2[kIdx]], nblk_d[nfset_ld2[kIdx]];
  Double_t cblkid_d[nfset_ld2[kIdx]][kNcell], cblke_d[nfset_ld2[kIdx]][kNcell], cblkatime_d[nfset_ld2[kIdx]][kNcell];
  Double_t cblkagain_d[nfset_ld2[kIdx]][kNcell];
  Double_t ekineW2_d[nfset_ld2[kIdx]];

  string configfilename_h[nfset_lh2[kIdx]];
  for( Int_t h=0; h<nfset_lh2[kIdx]; h++){
    Int_t field = fset_lh2[kIdx][h];
    if(field==-1) cout << "Error. configfilename array out of bounds." << endl;
    configfilename_h[h] = Form("../config/GMn/SBS%d/secal_lh2_sbs%d_f%d.cfg",kine,kine,fset_lh2[kIdx][h]);
  }
  string configfilename_d[nfset_ld2[kIdx]];
  for( Int_t d=0; d<nfset_ld2[kIdx]; d++){
    Int_t field = fset_ld2[kIdx][d];
    if(field==-1) cout << "Error. configfilename array out of bounds." << endl;
    configfilename_d[d] = Form("../config/GMn/SBS%d/secal_ld2_sbs%d_f%d.cfg",kine,kine,fset_ld2[kIdx][d]);
  }

  // Declare general physics/fit parameters
  Int_t minEventPerCell = 30; // Minimum number of scattered p in cell required to calibrate
  Int_t maxEventPerCell = 4000; // Maximum number of scattered p events to contribute
  Double_t HCal_divx = 0.15875; // Transverse width in x and y per cell
  Double_t HCal_divy = 0.15494;
  // Double_t HCal_Xi = -2.355005; // Distance from beam center to top of HCal in m from MC database
  // Double_t HCal_Xf = 1.454995; // Distance from beam center to bottom of HCal in m from MC database
  // Double_t HCal_Yi = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
  // Double_t HCal_Yf = 0.92964;
  Double_t HCal_Xi = -2.3531; // Distance from beam center to top of HCal in m from MC database
  Double_t HCal_Xf = 1.45309; // Distance from beam center to bottom of HCal in m from MC database
  Double_t HCal_Yi = -0.93155; // Distance from beam center to opposite-beam side of HCal in m from MC database
  Double_t HCal_Yf = 0.93155;
  Double_t highDelta = 0.05; // Minimum M(i,j)/b(i) factor allowed 

  // Declare general physics parameters to be modified by input config file - LH2
  //Double_t E_e_h[nfset_lh2[kIdx]] = {0.}; // Energy of beam (incoming electrons from accelerator)  
  Double_t E_e_h[maxfset] = {-1.}; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d_h[maxfset] = {-1.}; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th_h[maxfset] = {-1.}; // Angle that the center of HCal is at
  Double_t BB_th_h[maxfset] = {-1.}; // Angle that the center of BBCal is at
  Double_t W2_mean_h[maxfset] = {-1.}; // Mean of W at current kinematic
  Double_t W2_sig_h[maxfset] = {-1.}; // Width of W at current kinematic
  Double_t dx0_n_h[maxfset] = {-1.}; // Position of neutron spot, x-x_expected
  Double_t dx0_p_h[maxfset] = {-1.}; // Position of proton spot, x-x_expected
  Double_t dy0_h[maxfset] = {-1.}; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n_h[maxfset] = {-1.}; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p_h[maxfset] = {-1.}; // Max spread of proton spot, x-x_expected
  Double_t dy_sig_h[maxfset] = {-1.}; // Max spread of hadron spot, y-y_expected
  Double_t atime0_h[maxfset] = {-1.}; // Expected location in ADC time of signal
  Double_t atime_sig_h[maxfset] = {-1.}; // 1 sig of atime distribution
  Int_t useAlshield_h[maxfset] = {0};

  // Declare general physics parameters to be modified by input config file - LD2
  //Double_t E_e_d[nfset_ld2[kIdx]] = {0.}; // Energy of beam (incoming electrons from accelerator)
  Double_t E_e_d[maxfset] = {-1.}; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d_d[maxfset] = {-1.}; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th_d[maxfset] = {-1.}; // Angle that the center of HCal is at
  Double_t BB_th_d[maxfset] = {-1.}; // Angle that the center of BBCal is at
  Double_t W2_mean_d[maxfset] = {-1.}; // Mean of W at current kinematic
  Double_t W2_sig_d[maxfset] = {-1.}; // Width of W at current kinematic
  Double_t dx0_n_d[maxfset] = {-1.}; // Position of neutron spot, x-x_expected
  Double_t dx0_p_d[maxfset] = {-1.}; // Position of proton spot, x-x_expected
  Double_t dy0_d[maxfset] = {-1.}; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n_d[maxfset] = {-1.}; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p_d[maxfset] = {-1.}; // Max spread of proton spot, x-x_expected
  Double_t dy_sig_d[maxfset] = {-1.}; // Max spread of hadron spot, y-y_expected
  Double_t atime0_d[maxfset] = {-1.}; // Expected location in ADC time of signal
  Double_t atime_sig_d[maxfset] = {-1.}; // 1 sig of atime distribution
  Int_t useAlshield_d[maxfset] = {0};

  Int_t elasYield = 0; // Keep track of total elastics analyzed
  Int_t badtimeblkclus = 0; // Keep track of multi-blk clusters with out of time blks
  Int_t multblkclus = 0; // Keep track of total multi-block clusters (nblk>1)

  //For position reconstruction
  Double_t HCal_Xmin = HCal_Xi-HCal_divx/2;
  Double_t HCal_Xmax = HCal_Xf+HCal_divx/2;
  Double_t HCal_Ymin = HCal_Yi-HCal_divy/2;
  Double_t HCal_Ymax = HCal_Yf+HCal_divy/2;

  // Reading all relevant config files for lh2 and setting up chain
  for( Int_t f=0; f<nfset_lh2[kIdx]; f++ ){
    Ch[f] = new TChain("T");

    ifstream configfile(configfilename_h[f]);
    TString currentline;
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){

      if( !currentline.BeginsWith("#") ){
	Ch[f]->Add(currentline);
	report << "Loading file: " << currentline << ".." << endl;
      }    
    }
    //TCut globalcut = "";
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut_h[f] += currentline;
      }    
    }
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      Int_t ntokens = tokens->GetEntries();
      if( ntokens>1 ){
	TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
	if( skey == "E_e" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  E_e_h[f] = sval.Atof();
	  report << "Loading beam energy: " << E_e_h[f] << endl;
	}
	if( skey == "HCal_d" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_d_h[f] = sval.Atof();
	  report << "Loading HCal distance: " << HCal_d_h[f] << endl;
	}
	if( skey == "HCal_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_th_h[f] = sval.Atof() * TMath::DegToRad();	
	  report << "Loading HCal angle: " << HCal_th_h[f] << endl;
	}
	if( skey == "BB_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  BB_th_h[f] = sval.Atof() * TMath::DegToRad();	
	  report << "Loading BBCal angle: " << BB_th_h[f] << endl;
	}
	if( skey == "W2_mean" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_mean_h[f] = sval.Atof();
	  report << "Loading W2 mean cut: " << W2_mean_h[f] << endl;
	}
	if( skey == "W2_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_sig_h[f] = sval.Atof();
	  report << "Loading W2 sigma cut: " << W2_sig_h[f] << endl;
	}
	if( skey == "dx0_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_n_h[f] = sval.Atof();
	  report << "Loading x position of neutron spot: " << dx0_n_h[f] << endl;
	}
	if( skey == "dx0_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_p_h[f] = sval.Atof();
	  report << "Loading y position of proton spot: " << dx0_p_h[f] << endl;
	}
	if( skey == "dy0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy0_h[f] = sval.Atof();
	  report << "Loading y position of both hadron spots: " << dy0_h[f] << endl;
	}
	if( skey == "dx_sig_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_n_h[f] = sval.Atof();
	  report << "Loading x sigma of neutron spot: " << dx_sig_n_h[f] << endl;
	}
	if( skey == "dx_sig_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_p_h[f] = sval.Atof();
	  report << "Loading x sigma of proton spot: " << dx_sig_p_h[f] << endl;
	}
	if( skey == "dy_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy_sig_h[f] = sval.Atof();
	  report << "Loading y sigma of both hadron spots: " << dy_sig_h[f] << endl;
	}
	if( skey == "atime0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime0_h[f] = sval.Atof();
	  report << "Loading ADC time mean: " << atime0_h[f] << endl;
	}
	if( skey == "atime_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime_sig_h[f] = sval.Atof();
	  report << "Loading ADC time sigma: " << atime_sig_h[f] << endl;
	}
	if( skey == "useAlshield" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  useAlshield_h[f] = sval.Atoi();
	  report << "Loading Aluminum absorber option: " << useAlshield_h[f] << endl;
	}
      }
      delete tokens;
    }

    NTevents_h[f] = Ch[f]->GetEntries();
    elist_h[f] = new TEventList(Form("elist_lh2_sbs%d_f%d",kine,fset_lh2[kIdx][f]),Form("Elastic Event List, LH2, SBS%d, Mag%d",kine,fset_lh2[kIdx][f]));
    Ch[f]->Draw(Form(">>elist_lh2_sbs%d_f%d",kine,fset_lh2[kIdx][f]),globalcut_h[f]);
    Nevents_h[f] = elist_h[f]->GetN();

    Ch[f]->SetBranchStatus( "*", 0 );
    Ch[f]->SetBranchStatus( "sbs.hcal.x", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.y", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.e", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.rowblk", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.colblk", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.nblk", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    Ch[f]->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.chi2", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.n", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.px", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.py", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.pz", 1 );    
    Ch[f]->SetBranchStatus( "bb.tr.vz", 1 );
    Ch[f]->SetBranchStatus( "bb.tr.p", 1 );
    Ch[f]->SetBranchStatus( "bb.ps.e", 1 );
    Ch[f]->SetBranchStatus( "bb.ps.x", 1 );
    Ch[f]->SetBranchStatus( "bb.ps.y", 1 );
    Ch[f]->SetBranchStatus( "bb.sh.e", 1 );
    Ch[f]->SetBranchStatus( "bb.sh.x", 1 );
    Ch[f]->SetBranchStatus( "bb.sh.y", 1 );
    Ch[f]->SetBranchStatus( "e.kine.W2", 1 );

    Ch[f]->SetBranchAddress( "sbs.hcal.x", &HCALx_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.y", &HCALy_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.e", &HCALe_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.rowblk", &crow_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.colblk", &ccol_h[f] );
    Ch[f]->SetBranchAddress( "sbs.hcal.nblk", &nblk_h[f] ); // Total number of blocks in highest E clus
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid_h[f] ); // kNcell-1 index for each block
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke_h[f] ); // Array of block energies
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime_h[f] ); // Array of block ADC times
    Ch[f]->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain_h[f] ); // Array of block ADC gain coeff
    Ch[f]->SetBranchAddress( "bb.tr.n", &BBtr_n_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.chi2", BBtr_chi2_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.px", BBtr_px_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.py", BBtr_py_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.pz", BBtr_pz_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.vz", BBtr_vz_h[f] );
    Ch[f]->SetBranchAddress( "bb.tr.p", BBtr_p_h[f] );
    Ch[f]->SetBranchAddress( "bb.ps.e", &BBps_e_h[f] );
    Ch[f]->SetBranchAddress( "bb.ps.x", &BBps_x_h[f] );
    Ch[f]->SetBranchAddress( "bb.ps.y", &BBps_y_h[f] );
    Ch[f]->SetBranchAddress( "bb.sh.e", &BBsh_e_h[f] );
    Ch[f]->SetBranchAddress( "bb.sh.x", &BBsh_x_h[f] );
    Ch[f]->SetBranchAddress( "bb.sh.y", &BBsh_y_h[f] ); 
    Ch[f]->SetBranchAddress( "e.kine.W2", &ekineW2_h[f] );


  }//end config file for over lh2 field settings

 // Reading all relevant config files for ld2
  for( Int_t f=0; f<nfset_ld2[kIdx]; f++ ){
    Cd[f] = new TChain("T");

    ifstream configfile(configfilename_d[f]);
    TString currentline;
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	Cd[f]->Add(currentline);
	report << "Loading file: " << currentline << ".." << endl;
      }    
    }
    //TCut globalcut = "";
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut_d[f] += currentline;
      }    
    }
    while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      Int_t ntokens = tokens->GetEntries();
      if( ntokens>1 ){
	TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
	if( skey == "E_e" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  E_e_d[f] = sval.Atof();
	  report << "Loading beam energy: " << E_e_d[f] << endl;
	}
	if( skey == "HCal_d" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_d_d[f] = sval.Atof();
	  report << "Loading HCal distance: " << HCal_d_d[f] << endl;
	}
	if( skey == "HCal_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  HCal_th_d[f] = sval.Atof() * TMath::DegToRad();	
	  report << "Loading HCal angle: " << HCal_th_d[f] << endl;
	}
	if( skey == "BB_th" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  BB_th_d[f] = sval.Atof() * TMath::DegToRad();	
	  report << "Loading BBCal angle: " << BB_th_d[f] << endl;
	}
	if( skey == "W2_mean" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_mean_d[f] = sval.Atof();
	  report << "Loading W2 mean cut: " << W2_mean_d[f] << endl;
	}
	if( skey == "W2_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  W2_sig_d[f] = sval.Atof();
	  report << "Loading W2 sigma cut: " << W2_sig_d[f] << endl;
	}
	if( skey == "dx0_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_n_d[f] = sval.Atof();
	  report << "Loading x position of neutron spot: " << dx0_n_d[f] << endl;
	}
	if( skey == "dx0_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx0_p_d[f] = sval.Atof();
	  report << "Loading y position of proton spot: " << dx0_p_d[f] << endl;
	}
	if( skey == "dy0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy0_d[f] = sval.Atof();
	  report << "Loading y position of both hadron spots: " << dy0_d[f] << endl;
	}
	if( skey == "dx_sig_n" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_n_d[f] = sval.Atof();
	  report << "Loading x sigma of neutron spot: " << dx_sig_n_d[f] << endl;
	}
	if( skey == "dx_sig_p" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dx_sig_p_d[f] = sval.Atof();
	  report << "Loading x sigma of proton spot: " << dx_sig_p_d[f] << endl;
	}
	if( skey == "dy_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  dy_sig_d[f] = sval.Atof();
	  report << "Loading y sigma of both hadron spots: " << dy_sig_d[f] << endl;
	}
	if( skey == "atime0" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime0_d[f] = sval.Atof();
	  report << "Loading ADC time mean: " << atime0_d[f] << endl;
	}
	if( skey == "atime_sig" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  atime_sig_d[f] = sval.Atof();
	  report << "Loading ADC time sigma: " << atime_sig_d[f] << endl;
	}
	if( skey == "useAlshield" ){
	  TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	  useAlshield_d[f] = sval.Atoi();
	  report << "Loading Aluminum absorber option: " << useAlshield_d[f] << endl;
	}
      }
      delete tokens;
    }

    NTevents_d[f] = Cd[f]->GetEntries();
    elist_d[f] = new TEventList(Form("elist_ld2_sbs%d_f%d",kine,fset_ld2[kIdx][f]),Form("Elastic Event List, LD2, SBS%d, Mag%d",kine,fset_ld2[kIdx][f]));
    Cd[f]->Draw(Form(">>elist_ld2_sbs%d_f%d",kine,fset_ld2[kIdx][f]),globalcut_d[f]);
    Nevents_d[f] = elist_d[f]->GetN();

    Cd[f]->SetBranchStatus( "*", 0 );
    Cd[f]->SetBranchStatus( "sbs.hcal.x", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.y", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.e", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.rowblk", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.colblk", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.nblk", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
    Cd[f]->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.chi2", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.n", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.px", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.py", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.pz", 1 );    
    Cd[f]->SetBranchStatus( "bb.tr.vz", 1 );
    Cd[f]->SetBranchStatus( "bb.tr.p", 1 );
    Cd[f]->SetBranchStatus( "bb.ps.e", 1 );
    Cd[f]->SetBranchStatus( "bb.ps.x", 1 );
    Cd[f]->SetBranchStatus( "bb.ps.y", 1 );
    Cd[f]->SetBranchStatus( "bb.sh.e", 1 );
    Cd[f]->SetBranchStatus( "bb.sh.x", 1 );
    Cd[f]->SetBranchStatus( "bb.sh.y", 1 );
    Cd[f]->SetBranchStatus( "e.kine.W2", 1 );

    Cd[f]->SetBranchAddress( "sbs.hcal.x", &HCALx_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.y", &HCALy_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.e", &HCALe_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.rowblk", &crow_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.colblk", &ccol_d[f] );
    Cd[f]->SetBranchAddress( "sbs.hcal.nblk", &nblk_d[f] ); // Total number of blocks in highest E clus
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid_d[f] ); // kNcell-1 index for each block
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke_d[f] ); // Array of block energies
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime_d[f] ); // Array of block ADC times
    Cd[f]->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain_d[f] ); // Array of block ADC gain coeff
    Cd[f]->SetBranchAddress( "bb.tr.n", &BBtr_n_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.chi2", BBtr_chi2_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.px", BBtr_px_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.py", BBtr_py_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.pz", BBtr_pz_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.vz", BBtr_vz_d[f] );
    Cd[f]->SetBranchAddress( "bb.tr.p", BBtr_p_d[f] );
    Cd[f]->SetBranchAddress( "bb.ps.e", &BBps_e_d[f] );
    Cd[f]->SetBranchAddress( "bb.ps.x", &BBps_x_d[f] );
    Cd[f]->SetBranchAddress( "bb.ps.y", &BBps_y_d[f] );
    Cd[f]->SetBranchAddress( "bb.sh.e", &BBsh_e_d[f] );
    Cd[f]->SetBranchAddress( "bb.sh.x", &BBsh_x_d[f] );
    Cd[f]->SetBranchAddress( "bb.sh.y", &BBsh_y_d[f] ); 
    Cd[f]->SetBranchAddress( "e.kine.W2", &ekineW2_d[f] );

  }//end config file for over ld2 field settings

  cout << endl << endl << "Setup parameters loaded and chains linked." << endl;
  report << endl << endl << "Setup parameters loaded and chains linked." << endl;

  // Create stopwatch to track processing time
  TStopwatch *sw = new TStopwatch();
  sw->Start();
  
  // Declare outfile
  TFile *fout;
  if( qreplay ){
    fout = new TFile( Form("reporting/qreplay_sbs%d.root", kine ), "RECREATE" );
  }else{
    fout = new TFile( Form("reporting/eCalEOut_sbs%d.root", kine ), "RECREATE" );
  }

  report << "Output file located here: " << Form("reporting/qreplay_sbs%d.root", kine ) << endl;

  // Initialize vectors and arrays
  Double_t GCoeff[kNcell] = {0.0};
  Double_t GCoeff_oneblock[kNcell] = {0.0};
  Double_t GCoeff_divide[kNcell] = {0.0};
  Double_t gOldConst[kNcell] = {0.0};

  //No histograms needed. Will have to replay data with new database params, then use diagnostic script to assess
  Int_t lh2Events = 0;
  for( Int_t f=0; f<nfset_lh2[kIdx]; f++ ){
    lh2Events+=Nevents_h[f];
  }
  Int_t ld2Events = 0;
  for( Int_t f=0; f<nfset_ld2[kIdx]; f++ ){
    ld2Events+=Nevents_d[f];
  }
  
  cout << "Opened many trees with " << lh2Events << " LH2 events passing global cuts." << endl;
  cout << "Opened many trees with " << ld2Events << " LD2 events passing global cuts." << endl << endl;
  cout << "Total available events to calibrate = " << lh2Events + ld2Events << "." << endl << endl;

  report << "Opened many trees with " << lh2Events << " LH2 events passing global cuts." << endl;
  report << "Opened many trees with " << ld2Events << " LD2 events passing global cuts." << endl << endl;
  report << "Total available events to calibrate = " << lh2Events + ld2Events << "." << endl << endl;

  //Declare matrices for chi-square min calibration scheme and keep track of calibrated events
  TMatrixD Ma(kNcell,kNcell);
  TMatrixD Ma_oneblock(kNcell,kNcell);
  TMatrixD Ma_err(kNcell,kNcell);
  TVectorD ba(kNcell);
  TVectorD ba_oneblock(kNcell);
  TVectorD ba_err(kNcell);
  Int_t TNEV_h = 0;
  Int_t TNEV_d = 0;
  Int_t NEV[kNcell] = {0};
  Int_t NEV_oneblock[kNcell] = {0};
  Double_t err[kNcell] = {0.};
  Double_t err_ev[kNcell] = {0.};
  Double_t err_oneblock[kNcell] = {0.};
  Double_t err_ev_oneblock[kNcell] = {0.};
  Int_t deadclus = 0;
  
  //Declare diagnostic histograms (as sparse as possible)
  TH2D *hdx_mag_h = new TH2D( "hdx_mag_h", "Delta X vs Field Setting (LH2); field (percent); x_{HCAL} - x_{exp} (m)", 21, 0, 105, 250, -4, 3 );
  TH2D *hdy_mag_h = new TH2D( "hdy_mag_h", "Delta Y vs Field Setting (LH2); field (percent); y_{HCAL} - y_{exp} (m)", 21, 0, 105, 250, -1.25, 1.25 );
  TH2D *hdx_mag_d = new TH2D( "hdx_mag_d", "Delta X vs Field Setting (LD2); field (percent); x_{HCAL} - x_{exp} (m)", 21, 0, 105, 250, -4, 3 );
  TH2D *hdy_mag_d = new TH2D( "hdy_mag_d", "Delta Y vs Field Setting (LD2); field (percent); y_{HCAL} - y_{exp} (m)", 21, 0, 105, 250, -1.25, 1.25 );
  TH1D *hHCALe = new TH1D( "hHCALe","HCal Cluster E", 400, 0., 1. );
  TH2D *hHCALeX = new TH2D( "hHCALeX","HCal Cluster E vs X", kNrows, HCal_Xi, HCal_Xf, 400, 0., 1. );
  TH2D *hHCALeY = new TH2D( "hHCALeY","HCal Cluster E vs Y", kNcols, HCal_Yi, HCal_Yf, 400, 0., 1. );
  TH1D *hHCALblke = new TH1D( "hHCALblke","HCal Cluster Seed E", 400, 0., 1. );
  TH1D *hSampFrac = new TH1D( "hSampFrac","W2 Cut HCal Cluster E / Expected KE", 400, 0., 1. );
  TH1D *hblkid = new TH1D( "hblkid","HCal Block ID", 400, -50, 350 );
  TH2D *hSampFracvID = new TH2D( "hSampFracvID","HCal Cluster E / Expected KE vs ID", kNcell, 0, kNcell, 400, 0., 1. );
  TH2D *hSampFracvRow = new TH2D( "hSampFracvRow","HCal Cluster E / Expected KE vs Row", kNrows, 0, kNrows, 400, 0., 1. );
  TH2D *hSampFracvCol = new TH2D( "hSampFracvCol","HCal Cluster E / Expected KE vs Col", kNcols, 0, kNcols, 400, 0., 1. );

  //Loop over events
  cout << "Looping over hydrogen data.." << endl;
  report << "Looping over hydrogen data.." << endl;
  if(pass0) gROOT->ProcessLine( "gErrorIgnoreLevel = 6001" ); //Suppress error output to avoid undetectable sbs.hcal.clus_blk.again for pass0

  //Loop over all hydrogen data
  for( Int_t f=0; f<nfset_lh2[kIdx]; f++ ){
    
    Int_t hfieldset = fset_lh2[kIdx][f];

    //Declare energy loss parameters for beam going through the target
    Double_t pBeam = E_e_h[f]/(1.+E_e_h[f]/M_p*(1.-cos(BB_th_h[f])));
    Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th_h[f]) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)
    
    for( Long64_t nevent = 1; nevent <Nevents_h[f]; nevent++){

      if ( nevent%10000==0 && !qreplay ){
	cout << "LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << ". Total events gathered for calibration: " << TNEV_h << " \r";
      }else if( nevent%10000==0 ){
	cout << "qreplay at LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << " \r";
      }
      cout.flush();
      
      if ( nevent%100000==0 && !qreplay){
	report << "LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << ". Block atime cut dead clusters: " << deadclus << ". Total events gathered for calibration: " << TNEV_h << endl;
      }else if( nevent%100000==0 ){
	report << "qreplay at LH2 kinematic " << kine << " at field " << hfieldset << "% , entry: " << nevent << "/" << Nevents_h[f] << " \r";
      }

      Ch[f]->GetEntry( elist_h[f]->GetEntry( nevent ) ); 


      //////////////////////////////////////////////////////////////////////
      //Make cuts that don't rely on projections first to improve efficiency
      if( crow_h[f]==0 || 
	  crow_h[f]==23 || 
	  ccol_h[f]==0 || 
	  ccol_h[f]==11 ) continue; //All events with primary cluster element on edge blocks cut
      if( abs(cblkatime_h[f][0]-atime0_h[f])>ifac*atime_sig_h[f] ) continue; //All events where adctime outside of reasonable window cut
      //////////////////////////////////////////////////////////////////////


      Double_t A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
      Double_t A_oneblock[kNcell] = {0.0}; // Array to keep track of ADC values per cell for one block clusters only. Outscope on each ev
      
      
      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz_h[f][0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e_h[f] - Eloss;
      Double_t p_corr = BBtr_p_h[f][0] - Eloss_outgoing; //Neglecting the mass of e'

      Double_t etheta = acos( BBtr_pz_h[f][0]/BBtr_p_h[f][0] );
      Double_t ephi = atan2( BBtr_py_h[f][0], BBtr_px_h[f][0] );

      TVector3 vertex(0,0,BBtr_vz_h[f][0]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_corr,E_corr); //Mass of e negligable
      TLorentzVector kprime(BBtr_px_h[f][0],BBtr_py_h[f][0],BBtr_pz_h[f][0],BBtr_p_h[f][0]);
      TLorentzVector Ptarg(0,0,0,M_p); // assume proton for both LH2 and LD2 - can refine where useful with exclusive LH2 data set. Likely better to refine first with dxdy spot cuts on both protons and neutrons.

      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

      Double_t pel = E_corr/(1.+E_corr/M_p*(1.-cos(etheta)));
      Double_t nu = E_corr - BBtr_p_h[f][0];
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu);
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz_h[f][0])/pp ); //use elastic constraint on nucleon kinematics
	
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th_h[f]),0,cos(-HCal_th_h[f]));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
      TVector3 HCAL_origin = HCal_d_h[f] * HCAL_zaxis + hcalheight * HCAL_xaxis;

      //Define intersection points for hadron vector
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );

      TVector3 HCAL_intersect = vertex + sintersect * pNhat;

      Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p_h[f][0],2) ); // Obtain the scattered electron energy
	
      Double_t p_ep = BBtr_p_h[f][0];
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz_h[f][0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
	
      Double_t W = PgammaN.M();
      Double_t W2 = ekineW2_h[f];

      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      Double_t E_pp = nu+M_p; // Get energy of the proton
      Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      Double_t KE_p = nu; //For elastics
      //Double_t KE_p = Enucleon; //For sanity check on sampling fraction (2.28.23)
      Double_t dx = HCALx_h[f] - xexpect_HCAL;
      Double_t dy = HCALy_h[f] - yexpect_HCAL;
      
      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      if( fabs(W2-W2_mean_h[f])>W2_sig_h[f] ) continue; //Observed mean W2 cut on elastic peak
      elasYield++; //Events that pass the above cuts constitute elastics
      /////////////////////////////////////////////////

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe_h[f]/KE_p;
      Double_t E_exp = KE_p*sampFrac[kIdx]/pCsigfactor[kIdx]; //Currently testing pCsigfactor
      Int_t rec_row = ( HCALx_h[f] - HCal_Xmin )/HCal_divx;
      Int_t rec_col = ( HCALy_h[f] - HCal_Ymin )/HCal_divy;
      Int_t rec_cell = rec_row*kNcols + rec_col;

      hdx_mag_h->Fill( hfieldset, dx );
      hdy_mag_h->Fill( hfieldset, dy );
      
      //////////////////////////////////////////////////////////////////
      //Cut on dx and dy.
      bool pass_y = abs(dy-dy0_h[f])<ifac*dy_sig_h[f];
      if( !pass_y ) continue;
      bool pass_p = abs(dx-dx0_p_h[f])<ifac*dx_sig_p_h[f];
      bool pass_n = abs(dx-dx0_n_h[f])<ifac*dx_sig_n_h[f]; //Redundant for LH2
      if( !pass_p && !pass_n ) continue; //Cut on both n and p spots for each event, cannot know which apriori
      //////////////////////////////////////////////////////////////////

      //Write out diagnostic histograms
      Double_t clusE = 0.0;
      Double_t clusblkE = 0.0;
      for( Int_t blk = 0; blk<(int)nblk_h[f]; blk++ ){
      	Int_t blkid = int(cblkid_h[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)

	if( qreplay ){
	  if( pass0 ){
	    clusE += cblke_h[f][blk]/gOldConst_pass0[blkid]*gConst_iter1[blkid];
	    if( blk==0 ) clusblkE += cblke_h[f][blk]/gOldConst_pass0[blkid]*gConst_iter1[blkid];
	  }
	  if( !pass0 ){
	    clusE += cblke_h[f][blk]/cblkagain_h[f][blk]*gConst_iter1[blkid];
	    if( blk==0 ) clusblkE += cblke_h[f][blk]/cblkagain_h[f][blk]*gConst_iter1[blkid];
	  }
	}else{
	  clusE += cblke_h[f][blk];
	  if( blk==0 ) clusblkE += cblke_h[f][blk];
	}
      }

      hHCALe->Fill( clusE );
      hHCALeX->Fill( HCALx_h[f], clusE );
      hHCALeY->Fill( HCALy_h[f], clusE );
      hHCALblke->Fill( clusblkE );
      hSampFrac->Fill( clusE/KE_p );
      hSampFracvID->Fill( cblkid_h[f][0], clusE/KE_p );
      hSampFracvRow->Fill( crow_h[f], clusE/KE_p );
      hSampFracvCol->Fill( ccol_h[f], clusE/KE_p );
      

      //////////////////////////////////////////////////////////////////
      //Continue if iter == 1 (second quasi-replay pass for diagnostics)
      if( qreplay ) continue;
      //////////////////////////////////////////////////////////////////

      // Get pC from primary clusters only
      bool badclus = false;
      clusE = 0.0;
      Double_t cluspC = 0.0;
      for( Int_t blk = 0; blk<(int)nblk_h[f]; blk++ ){
	Int_t blkid = int(cblkid_h[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	hblkid->Fill( blkid );

	//Populate list of old gain coefficients in case dof is bad
	if( gOldConst[blkid]==0.0 ) gOldConst[blkid]=cblkagain_h[f][blk];
	
	
	//Diagnostic cut for now
	if( abs(cblkatime_h[f][blk]-cblkatime_h[f][0])>ifac*atime_sig_h[f] ) {
	  badclus = true;
	  continue; //Cut blocks that are out of time
	}
	
	clusE += cblke_h[f][blk];

	// Get pC with old gain coefficients where pass0 extracted from db_sbs.hcal.dat
	if( pass0 ){
	  cluspC += cblke_h[f][blk]/gOldConst_pass0[blkid]; 
	  A[blkid] += cblke_h[f][blk]/gOldConst_pass0[blkid];	
	  if(nblk_h[f]==1) {
	    A_oneblock[blkid] += cblke_h[f][blk]/gOldConst_pass0[blkid];
	  }
	}else{
	  cluspC += cblke_h[f][blk]/cblkagain_h[f][blk]; 
	  A[blkid] += cblke_h[f][blk]/cblkagain_h[f][blk];	
	  if(nblk_h[f]==1) {
	    A_oneblock[blkid] += cblke_h[f][blk]/cblkagain_h[f][blk];
	  }
	}

	err_ev[blkid] += pow( ((1.+cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Getting ready for std. dev of mean

	// Simple estimation of the coefficients assuming 100% energy deposition in one block. Will check against the more sophisticated version with chi^2 reduction.
	if(nblk_h[f]==1) {
	 
	  err_ev_oneblock[blkid] += pow( ((1.+0.5*cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Probably better approximation of the error since all the KE_p should be deposited in this block.
	  NEV_oneblock[blkid]++;
	}
	NEV[blkid]++;
	TNEV_h++;
      }

      if(clusE<0.025){
	deadclus++;
	report << "Warning: Cluster block cut on atime killed cluster." << endl;
	continue;
      }

      if( badclus ) badtimeblkclus++;
      if( nblk_h[f]>1 ) multblkclus++;

      //Calculate clusE variables
      Double_t E_dev = (E_exp-clusE)/E_exp; //Energy deviation as percent of expected
      
      //Build the matrix as simply as possible
      for(Int_t icol = 0; icol<kNcell; icol++){
	ba(icol)+= A[icol];

	if(nblk_h[f]==1) ba_oneblock(icol)+=A_oneblock[icol];

	for(Int_t irow = 0; irow<kNcell; irow++){
	  Ma(icol,irow) += A[icol]*A[irow]/E_exp;
	  if(nblk_h[f]==1){	    
	    // cout << "do something to correct oneblock" << endl;
	    Ma_oneblock(icol,irow) += A_oneblock[icol]*A_oneblock[irow]/E_exp;
	  } //oneblock
	} //loop over matrix element rows
      } //loop over matrix element cols
    } //loop over events
  } //loop for lh2

  cout << "Hydrogen loop complete." << endl;
  report << "Hydrogen loop complete." << endl;

  if( !qreplay ){
    Int_t cell = 0;

    cout << endl << "Number of events available for calibration from hydrogen alone: " << endl << endl;
    report << endl << "Number of events available for calibration from hydrogen alone: " << endl << endl;

    for( Int_t r = 0; r<kNrows; r++){
      for( Int_t c = 0; c<kNcols; c++){
	cout << NEV[cell] << "  ";
	report << NEV[cell] << "  ";
	cell++;
      }
      cout << endl;
      report << endl;
    }
  
    cout << endl;
    report << endl;

    usleep( 3*second ); //Give some time for review of step
  }

  cout << "Looping over deuterium data.." << endl << endl;
  report << "Looping over deuterium data.." << endl << endl;

  //Loop over all deuterium data
  for( Int_t f=0; f<nfset_ld2[kIdx]; f++ ){

    Int_t dfieldset = fset_ld2[kIdx][f];

    //Declare energy loss parameters for beam going through the target
    Double_t pBeam = E_e_d[f]/(1.+E_e_d[f]/M_p*(1.-cos(BB_th_d[f])));
    Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th_d[f]) * rho_tgt * dEdx_tgt; //Mean energy loss of the beam prior to the scattering, approximately 1 MeV, could correct further with raster position (likely negligable)

    for( Long64_t nevent = 1; nevent <Nevents_d[f]; nevent++){

      if ( nevent%10000==0 && !qreplay ){
	cout << "LD2 kinematic " << kine << " at field " << dfieldset << "% , entry: " << nevent << "/" << Nevents_d[f] << ". Total events gathered for calibration: " << TNEV_d << " \r";
      }else if( nevent%10000==0 ){
	cout << "qreplay at LD2 kinematic " << kine << " at field " << dfieldset << "% , entry: " << nevent << "/" << Nevents_d[f] << " \r";
      }
      cout.flush();
      
      if ( nevent%100000==0 && !qreplay){
	report << "LD2 kinematic " << kine << " at field " << dfieldset << "% , entry: " << nevent << "/" << Nevents_d[f] << ". Block atime cut dead clusters: " << deadclus << ". Total events gathered for calibration: " << TNEV_d << endl;
      }else if( nevent%100000==0 ){
	report << "qreplay at LD2 kinematic " << kine << " at field " << dfieldset << "% , entry: " << nevent << "/" << Nevents_d[f] << " \r";
      }
      
      Cd[f]->GetEntry( elist_d[f]->GetEntry( nevent ) ); 

      //////////////////////////////////////////////////////////////////////
      //Make cuts that don't rely on projections first to improve efficiency
      if( crow_d[f]==0 || 
      	  crow_d[f]==23 || 
      	  ccol_d[f]==0 || 
      	  ccol_d[f]==11 ) continue; //All events with primary cluster element on edge blocks cut
      if( abs(cblkatime_d[f][0]-atime0_d[f])>3*atime_sig_d[f] ) continue; //All events where adctime outside of reasonable window cut
      //////////////////////////////////////////////////////////////////////

      Double_t A[kNcell] = {0.0}; // Array to keep track of ADC values per cell. Outscope on each ev
      Double_t A_oneblock[kNcell] = {0.0}; // Array to keep track of ADC values per cell for one block clusters only. Outscope on each ev
      
      
      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz_d[f][0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e_d[f] - Eloss;
      Double_t p_corr = BBtr_p_d[f][0] - Eloss_outgoing; //Neglecting the mass of e'

      Double_t etheta = acos( BBtr_pz_d[f][0]/BBtr_p_d[f][0] );
      Double_t ephi = atan2( BBtr_py_d[f][0], BBtr_px_d[f][0] );

      TVector3 vertex(0,0,BBtr_vz_d[f][0]); // z location of vertex in hall coordinates
      TLorentzVector Pbeam(0,0,E_corr,E_corr); //Mass of e negligable
      TLorentzVector kprime(BBtr_px_d[f][0],BBtr_py_d[f][0],BBtr_pz_d[f][0],BBtr_p_d[f][0]);
      TLorentzVector Ptarg(0,0,0,M_p); // assume proton for both LH2 and LD2 - can refine where useful with exclusive LH2 data set. Likely better to refine first with dxdy spot cuts on both protons and neutrons.

      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

      Double_t pel = E_corr/(1.+E_corr/M_p*(1.-cos(etheta)));
      Double_t nu = E_corr - BBtr_p_d[f][0];
      Double_t pp = sqrt(pow(nu,2)+2.*M_p*nu);
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz_d[f][0])/pp ); //use elastic constraint on nucleon kinematics
	
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

      //Define HCal coordinate system
      TVector3 HCAL_zaxis(sin(-HCal_th_d[f]),0,cos(-HCal_th_d[f]));
      TVector3 HCAL_xaxis(0,-1,0);
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
	
      TVector3 HCAL_origin = HCal_d_d[f] * HCAL_zaxis + hcalheight * HCAL_xaxis;

      //Define intersection points for hadron vector
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / ( pNhat.Dot( HCAL_zaxis ) );

      TVector3 HCAL_intersect = vertex + sintersect * pNhat;

      Double_t yexpect_dCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_dCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );

      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p_d[f][0],2) ); // Obtain the scattered electron energy
	
      Double_t p_ep = BBtr_p_d[f][0];
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz_d[f][0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
	
      Double_t W = PgammaN.M();
      Double_t W2 = ekineW2_d[f];

      //Use the electron kinematics to predict the proton momentum assuming elastic scattering on free proton at rest:
      Double_t E_pp = nu+M_p; // Get energy of the proton
      Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same
      Double_t KE_p = nu; //For elastics
      //Double_t KE_p = Enucleon; //For sanity check on sampling fraction (2.28.23)
      Double_t dx = HCALx_d[f] - xexpect_dCAL;
      Double_t dy = HCALy_d[f] - yexpect_dCAL;
      
      /////////////////////////////////////////////////
      //Primary W2 cut on elastics and HCal active area
      if( fabs(W2-W2_mean_d[f])>W2_sig_d[f] ) continue; //Observed mean W2 cut on elastic peak
      elasYield++; //Events that pass the above cuts constitute elastics
      /////////////////////////////////////////////////

      //Calculate/declare new variables for analysis
      Double_t SFrac = HCALe_d[f]/KE_p;
      Double_t E_exp = KE_p*sampFrac[kIdx]/pCsigfactor[kIdx]; //Currently testing pCsigfactor
      Int_t rec_row = ( HCALx_d[f] - HCal_Xmin )/HCal_divx;
      Int_t rec_col = ( HCALy_d[f] - HCal_Ymin )/HCal_divy;
      Int_t rec_cell = rec_row*kNcols + rec_col;
      
      hdx_mag_d->Fill( dfieldset, dx );
      hdy_mag_d->Fill( dfieldset, dy );

      //////////////////////////////////////////////////////////////////
      //Cut on dx and dy.
      bool pass_y = abs(dy-dy0_d[f])<ifac*dy_sig_d[f];
      if( !pass_y ) continue;
      bool pass_p = abs(dx-dx0_p_d[f])<ifac*dx_sig_p_d[f];
      bool pass_n = abs(dx-dx0_n_d[f])<ifac*dx_sig_n_d[f];
      if( !pass_p && !pass_n ) continue; //Cut on both n and p spots for each event, cannot know which apriori      
      //////////////////////////////////////////////////////////////////

      //Write out diagnostic histograms
      Double_t clusE = 0.0;
      Double_t clusblkE = 0.0;
      for( Int_t blk = 0; blk<(int)nblk_d[f]; blk++ ){
      	Int_t blkid = int(cblkid_d[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)

	if( qreplay ){
	  if( pass0 ){
	    clusE += cblke_d[f][blk]/gOldConst_pass0[blkid]*gConst_iter1[blkid];
	    if( blk==0 ) clusblkE += cblke_d[f][blk]/gOldConst_pass0[blkid]*gConst_iter1[blkid];
	  }
	  if( !pass0 ){
	    clusE += cblke_d[f][blk]/cblkagain_d[f][blk]*gConst_iter1[blkid];
	    if( blk==0 ) clusblkE += cblke_d[f][blk]/cblkagain_d[f][blk]*gConst_iter1[blkid];
	  }
	}else{
	  clusE += cblke_d[f][blk];
	  if( blk==0 ) clusblkE += cblke_d[f][blk];
	}
      }

      hHCALe->Fill( clusE );
      hHCALeX->Fill( HCALx_d[f], clusE );
      hHCALeY->Fill( HCALy_d[f], clusE );
      hHCALblke->Fill( clusblkE );
      hSampFrac->Fill( clusE/KE_p );
      hSampFracvID->Fill( cblkid_d[f][0], clusE/KE_p );
      hSampFracvRow->Fill( crow_d[f], clusE/KE_p );
      hSampFracvCol->Fill( ccol_d[f], clusE/KE_p );

      //////////////////////////////////////////////////////////////////
      //Continue if iter == 1 (second quasi-replay pass for diagnostics)
      if( qreplay ) continue;
      //////////////////////////////////////////////////////////////////

      // Get pC from primary clusters only
      bool badclus = false;
      clusE = 0.0;
      Double_t cluspC = 0.0;
      for( Int_t blk = 0; blk<(int)nblk_d[f]; blk++ ){
	Int_t blkid = int(cblkid_d[f][blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	hblkid->Fill( blkid );

	//Populate list of old gain coefficients in case dof is bad
	if( gOldConst[blkid]==0.0 ) gOldConst[blkid]=cblkagain_d[f][blk];

	//Diagnostic cut for now
	if( abs(cblkatime_d[f][blk]-cblkatime_d[f][0])>ifac*atime_sig_d[f] ) {
	  badclus = true;
	  continue; //Cut blocks that are out of time
	}
	
	clusE += cblke_d[f][blk];

	// Get pC with old gain coefficients where pass0 extracted from db_sbs.hcal.dat
	if( pass0 ){
	  cluspC += cblke_d[f][blk]/gOldConst_pass0[blkid]; 
	  A[blkid] += cblke_d[f][blk]/gOldConst_pass0[blkid];	
	  if(nblk_d[f]==1) {
	    A_oneblock[blkid] += cblke_d[f][blk]/gOldConst_pass0[blkid];
	  }
	}else{
	  cluspC += cblke_d[f][blk]/cblkagain_d[f][blk]; 
	  A[blkid] += cblke_d[f][blk]/cblkagain_d[f][blk];	
	  if(nblk_d[f]==1) {
	    A_oneblock[blkid] += cblke_d[f][blk]/cblkagain_d[f][blk];
	  }
	}

	err_ev[blkid] += pow( ((1.+cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Getting ready for std. dev of mean

	// Simple estimation of the coefficients assuming 100% energy deposition in one block. Will check against the more sophisticated version with chi^2 reduction.
	if(nblk_d[f]==1) {
	 
	  err_ev_oneblock[blkid] += pow( ((1.+0.5*cluspC)/cluspC+(0.05*KE_p)/KE_p) , 2 ); //Probably better approximation of the error since all the KE_p should be deposited in this block.
	  NEV_oneblock[blkid]++;
	}
	NEV[blkid]++;
	TNEV_d++;
      }

      if(clusE<0.025){
	deadclus++;
	report << "Warning: Cluster block cut on atime killed cluster." << endl;
	continue;
      }

      if( badclus ) badtimeblkclus++;
      if( nblk_d[f]>1 ) multblkclus++;

      //Calculate clusE variables
      Double_t E_dev = (E_exp-clusE)/E_exp; //Energy deviation as percent of expected
      
      //Build the matrix as simply as possible
      for(Int_t icol = 0; icol<kNcell; icol++){
	ba(icol)+= A[icol];

	if(nblk_d[f]==1) ba_oneblock(icol)+=A_oneblock[icol];

	for(Int_t irow = 0; irow<kNcell; irow++){
	  Ma(icol,irow) += A[icol]*A[irow]/E_exp;
	  if(nblk_d[f]==1){	    
	    // cout << "do something to correct oneblock" << endl;
	    Ma_oneblock(icol,irow) += A_oneblock[icol]*A_oneblock[irow]/E_exp;
	  } //oneblock
	} //loop over matrix element rows
      } //loop over matrix element cols
    } //loop over events
  } //loop for ld2

  cout << "Deuterium loop complete." << endl << endl;
  report << "Deuterium loop complete." << endl << endl;

  if( !qreplay ){
    cout << endl << "Checking data, inverting matrix, and solving for coefficients.." << endl << endl;
    report << endl << "Checking data, inverting matrix, and solving for coefficients.." << endl << endl;
  }

  //Reject the bad cells and normalize the oneblock check
  Int_t badcell[kNcell];
  Int_t badcell_oneblock[kNcell];
  Int_t cellBad = 0;
  Double_t y[kNcell] = {0.0}; // For easy TGraphErrors build
  
  for(Int_t i=0; i<kNcell; i++){
    if( qreplay ){
      report << "Skipping matrix element checks for iteration 1 (quasi-replay).." << endl << endl;
      break;
    }else{
      report << "Proceeding to energy matrix quality check analysis.." << endl << endl;
    }

    badcell[i] = 0;
    y[i] = i;
    
    //Do not change ADC gain coeff if insufficient events or energy dep in cell
    if( NEV[i] < minEventPerCell || Ma(i,i) < highDelta*ba(i) ){ 

      cellBad = 1;

      Double_t elemRatio = Ma(i,i)/ba(i);

      ba(i) = 1.0;  // Set RHS vector for cell i to 1.0 
      Ma(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
      for(Int_t j=0; j<kNcell; j++){
	if( j != i ){
	  
	  Ma(i,j) = 0.0;
	  Ma(j,i) = 0.0;
	}
      }
      badcell[i] = 1;
      report << "Cell " << i << " bad" << endl;
      report << "Number of events in bad cell =" << NEV[i] << endl;
      report << "Matrix element/vector element ratio =" << elemRatio << endl << endl;
    } 

    badcell_oneblock[i] = 0;

    //Do not change ADC gain coeff if insufficient events or energy dep in cell, oneblock
    if( NEV_oneblock[i] < minEventPerCell || Ma_oneblock(i,i) < highDelta*ba_oneblock(i) ){ 

      cellBad = 1;

      Double_t elemRatio = Ma_oneblock(i,i)/ba_oneblock(i);

      ba_oneblock(i) = 1.0;  // Set RHS vector for cell i to 1.0 
      Ma_oneblock(i,i) = 1.0; // Set diagonal element of Matrix M for cell i to 1.0 
      
      for(Int_t j=0; j<kNcell; j++){
	if( j != i ){
	  
	  Ma_oneblock(i,j) = 0.0; 
	  Ma_oneblock(j,i) = 0.0;
	}
      }
      badcell_oneblock[i] = 1;
      report << "Oneblock: cell" << i << " bad" << endl;
      report << "Oneblock: number of events in bad cell =" << NEV_oneblock[i] << endl;
      report << "Oneblock: matrix element/vector element ratio =" << elemRatio << endl << endl;
    }   

    //Calculate error per block (element of A)
    err[i] = sqrt(err_ev[i]/NEV[i]);
    err_oneblock[i] = sqrt(err_ev_oneblock[i]/NEV_oneblock[i]);
  }

  if( cellBad!=0 && !qreplay ) cout << "Bad cells detected. See report for details." << endl << endl;
  if( cellBad==0 && !qreplay ) report << "No bad cells detected!" << endl << endl;

  //If iteration == 0, invert the matrix, solve for ratios
  if( !qreplay ){
    TMatrixD M_inv = Ma.Invert();
    TMatrixD M_inv_oneblock = Ma_oneblock.Invert();
    TVectorD Coeff = M_inv*ba; // Stays unmodified for reference
    TVectorD Coeff_oneblock = M_inv_oneblock*ba_oneblock; // Stays unmodified for reference

    for(Int_t i=0; i<kNcell; i++){
      if( pass0 ){
	if( gOldConst[i]==0.0 ) gOldConst[i]=gOldConst_pass0[i];
      }

      if(badcell_oneblock[i]==0){
	GCoeff_oneblock[i]=Coeff_oneblock[i];
      }else{
	GCoeff_oneblock[i]=gOldConst[i];
      }

      if(badcell[i]==0&&Coeff[i]>0){ //Only update the coefficent if coefficient passes quality checks
	GCoeff[i]=Coeff[i];
	GCoeff_divide[i]=GCoeff[i]/GCoeff_oneblock[i];
      }else{
	GCoeff[i]=gOldConst[i]; // If the cell calibration is bad, use the old coefficient
	GCoeff_divide[i]=-1.0;
      }
    }

    Double_t yErr[kNcell] = {0.0};

    TGraphErrors *ccgraph = new TGraphErrors( kNcell, y, GCoeff, yErr, err ); 
    ccgraph->GetXaxis()->SetLimits(0.0,288.0);  
    ccgraph->GetYaxis()->SetLimits(0.0,0.25);
    ccgraph->SetTitle("Calibration Coefficients");
    ccgraph->GetXaxis()->SetTitle("Channel");
    ccgraph->GetYaxis()->SetTitle("Unitless");
    ccgraph->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
    ccgraph->Write("constants");

    TGraphErrors *ccgraph_oneBlock = new TGraphErrors( kNcell, y, GCoeff_oneblock, yErr, err_oneblock ); 
    ccgraph_oneBlock->GetXaxis()->SetLimits(0.0,288.0);  
    ccgraph_oneBlock->GetYaxis()->SetLimits(0.0,0.25);
    ccgraph_oneBlock->SetTitle("Calibration Coefficients One Block");
    ccgraph_oneBlock->GetXaxis()->SetTitle("Channel");
    ccgraph_oneBlock->GetYaxis()->SetTitle("Unitless");
    ccgraph_oneBlock->SetMarkerStyle(21);
    ccgraph_oneBlock->Write("constants_oneblock");

    TGraphErrors *ccgraph_divide = new TGraphErrors( kNcell, y, GCoeff_divide, yErr, yErr ); 
    ccgraph_divide->GetXaxis()->SetLimits(0.0,288.0);  
    ccgraph_divide->SetTitle("Calibration Coefficients / OneBlock Coeff");
    ccgraph_divide->GetXaxis()->SetTitle("Channel");
    ccgraph_divide->GetYaxis()->SetTitle("Unitless");
    ccgraph_divide->SetMarkerStyle(21);
    ccgraph_divide->Write("constants_divide");

    //Declare coeff outfile
    ofstream GainCoeff;

    GainCoeff.open( constPath );

    //Console/txt outs
    Int_t cell = 0;

    GainCoeff << "#HCal gain coefficients from SBS-" << kine << " obtained " << date.c_str() << endl;
    GainCoeff << "sbs.hcal.adc.gain =" << endl;

    cout << "Gain Coefficients" << endl;
    report << "Gain Coefficients" << endl;
    for( Int_t r=0; r<kNrows; r++ ){
      for( Int_t c=0; c<kNcols; c++ ){
	GainCoeff << GCoeff[cell] << "  ";
	cout << GCoeff[cell] << "  ";
	report << GCoeff[cell] << "  ";
	cell++;
      }
      GainCoeff << endl;
      cout << endl;
      report << endl;
    }

    cell = 0;

    cout << endl << "One Block Std:" << endl;
    report << endl << "One Block Std:" << endl;
    GainCoeff << endl << "#One Block Std = " << endl;

    for( Int_t r = 0; r<kNrows; r++){
      for( Int_t c = 0; c<kNcols; c++){
	GainCoeff << GCoeff_oneblock[cell] << "  ";
	cout << GCoeff_oneblock[cell] << "  ";
	report << GCoeff_oneblock[cell] << "  ";
	cell++;
      }
      GainCoeff << endl;
      cout << endl;
    }

    cell = 0;

    cout << endl << "Total Number of events available for calibration: " << endl << endl;
    report << endl << "Total Number of events available for calibration: " << endl << endl;
    GainCoeff << endl << "#Number of events available for calibration = " << endl;

    for( Int_t r = 0; r<kNrows; r++){
      for( Int_t c = 0; c<kNcols; c++){
	GainCoeff << NEV[cell] << "  ";
	cout << NEV[cell] << "  ";
	report << NEV[cell] << "  ";
	cell++;
      }
      GainCoeff << endl;
      cout << endl;
      report << endl;
    }

    GainCoeff.close();
  } //end if iter==0

  fout->Write();
  
  if( !qreplay ){
    report << endl << endl << "Total blocks out of time with primary block / Multi-block clusters: " << badtimeblkclus << "/" << multblkclus << endl;
    
    cout << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total event bins available for calibration = LH2 + LD2: " << TNEV_h << " + " << TNEV_d << " = " << TNEV_h + TNEV_d << ". Total number of events analyzed: " << lh2Events + ld2Events << "." << endl;
    
    report << endl << endl << "Elastic yield for analyzed runs: " << elasYield << ". Total event bins available for calibration = LH2 + LD2: " << TNEV_h << " + " << TNEV_d << " = " << TNEV_h + TNEV_d << ". Total number of events analyzed: " << lh2Events + ld2Events << "." << endl;
  }    

  if( qreplay ){
    cout << endl << "Diagnostic iteration complete. Histograms written to file." << endl << endl;
    report << endl <<  "Diagnostic iteration complete. Histograms written to file." << endl << endl;
  }else{
    cout << endl <<  "Calibration complete. Constants and histograms written to file." << endl << endl;
    report << endl <<  "Calibration complete. Constants and histograms written to file." << endl << endl;
  }

  st->Stop();

  // Send time efficiency report to console/report
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  report << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  report.close();
}
