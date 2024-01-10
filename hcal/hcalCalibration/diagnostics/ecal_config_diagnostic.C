//sseeds 2.22.23
//Short script to:
// evaluate dx, dy, and W2 limits for hcal energy calibration by kinematic
// produce sampling fraction and cluster energy plots to quickly evaluate calibration performance
#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include "TLorentzVector.h"
#include "TCut.h"
#include "TLatex.h"
#include "TFile.h"

using namespace std;

//global variables
const Int_t kNcell = 288;
const Int_t maxTracks = 16; 
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.0051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;
const Double_t sampFrac = 0.0795;

const Double_t hcalheight = -0.2897;
const Double_t HCal_divy = 0.15494; // Horizontal length of all HCAL blocks in m from MC database
const Double_t HCal_divx = 0.15875; // Vertical length of all HCAL blocks in m from MC database
const Double_t HCal_Xi = -2.355005; // Distance from beam center to top of HCal in m from MC database
const Double_t HCal_Xf = 1.454995; // Distance from beam center to bottom of HCal in m from MC database
const Double_t HCal_Yi = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
const Double_t HCal_Yf = 0.92964; // Distance from beam center to beam side of HCal in m from MC database
const Int_t xN = 48; //2*kNrows, total number of dispersive bins detection uni
const Int_t yN = 24; //2*kNcols, total number of transverse bins detection uni
const Int_t uni_N = 400; // Total number of bins used to measure detection uniformity (hSampFrac histos)

void ecal_config_diagnostic( const char *experiment = "GEn", Int_t kine=3, const char *tar = "lh2", Int_t field=100, Int_t iter=0 ){//main

  std::string TAR = tar;

  if( kine==-1 || TAR.compare("")==0 || field==-1 ){
    cout << "Error: Input parameters out of bounds. Please execute with the format:" << endl;
    cout << "  root -l \"ecal_diagnostic.C( <kine>, <tar>, <mag> )\" " << endl;
    cout << "  ..where kine = { 4, 7, 8, 9, 11, 14 }" << endl;
    cout << "  ..where tar = { \"lh2\", \"ld2\" }" << endl;
    cout << "  ..where field = { 0, 30, 50, 70, 85, 100 }" << endl;
    cout << "  ..where iter = { 0, 1} (0 for full diagnostic and 1 for completed calibration review plots)" << endl;
    cout << "  ..and note that some of these combinations will not have associated data. Beware!" << endl;
    return;
  }

  TChain *C = new TChain("T");

  cout << "Selected Kinematic " << kine << "." << endl;

  // Declare file name paths
  TString configfilename;
  TString outputfilename;
  if( iter==0 ){
    configfilename = Form( "../config/%s/SBS%d/secal_%s_sbs%d_f%d.cfg", experiment, kine, tar, kine, field );
    outputfilename = Form( "output/%s/SBS%d/ecal_diag_%s_sbs%d_f%d.root", experiment, kine, tar, kine, field );
  }else{
    configfilename = Form( "../config/%s/SBS%d/calReview_%s_sbs%d_f%d.cfg", experiment, kine, tar, kine, field );
    outputfilename = Form( "output/%s/SBS%d/ecal_calReview_sbs%d.root", experiment, kine, kine );
  }

  Double_t E_e = -1000.; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d = -1000.; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = -1000.; // Angle that the center of HCal is at
  Double_t BB_th = -1000.; // Angle that the center of BBCal is at
  Double_t W2_mean = -1000.; // Mean of W at current kinematic
  Double_t W2_sig = -1000.; // Width of W at current kinematic
  Double_t dx0_n = -1000.; // Position of neutron spot, x-x_expected
  Double_t dx0_p = -1000.; // Position of proton spot, x-x_expected
  Double_t dy0 = -1000.; // Position of hadron spot, y-y_expected
  Double_t dx_sig_n = -1000.; // Max spread of neutron spot, x-x_expected
  Double_t dx_sig_p = -1000.; // Max spread of proton spot, x-x_expected
  Double_t dy_sig = -1000.; // Max spread of hadron spot, y-y_expected
  Double_t atime0 = -1000.; // Expected location in ADC time of signal
  Double_t atime_sig = -1000.; // 1 sig of atime distribution
  Int_t useAlshield = -1000;

 // Reading config file
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
    Int_t ntokens = tokens->GetEntries();
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
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_th = sval.Atof() * TMath::DegToRad();	
	cout << "Loading BBCal angle: " << BB_th << endl;
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
      if( skey == "dx0_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_n = sval.Atof();
	cout << "Loading x position of neutron spot: " << dx0_n << endl;
      }
      if( skey == "dx0_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_p = sval.Atof();
	cout << "Loading y position of proton spot: " << dx0_p << endl;
      }
      if( skey == "dy0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0 = sval.Atof();
	cout << "Loading y position of both hadron spots: " << dy0 << endl;
      }
      if( skey == "dx_sig_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_n = sval.Atof();
	cout << "Loading x sigma of neutron spot: " << dx_sig_n << endl;
      }
      if( skey == "dx_sig_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_p = sval.Atof();
	cout << "Loading x sigma of proton spot: " << dx_sig_p << endl;
      }
      if( skey == "dy_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig = sval.Atof();
	cout << "Loading y sigma of both hadron spots: " << dy_sig << endl;
      }
      if( skey == "atime0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	atime0 = sval.Atof();
	cout << "Loading ADC time mean: " << atime0 << endl;
      }
      if( skey == "atime_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	atime_sig = sval.Atof();
	cout << "Loading ADC time sigma: " << atime_sig << endl;
      }
      if( skey == "useAlshield" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	useAlshield = sval.Atoi();
	cout << "Loading Aluminum absorber option: " << useAlshield << endl;
      }
    }
    delete tokens;
  }

  if( kine == -1000 ||
      E_e == -1000 || 
      HCal_d == -1000 ||
      HCal_th == -1000 ||
      BB_th == -1000 ||
      W2_mean == -1000 ||
      W2_sig == -1000 ||
      dx0_n == -1000 ||
      dx0_p == -1000 ||
      dy0 == -1000 ||
      dx_sig_n == -1000 ||
      dx_sig_p == -1000 ||
      dy_sig == -1000 ||
      atime0 == -1000 ||
      atime_sig == -1000 ||
      useAlshield == -1000 ){
    cout << "Error: Setup parameters not fully loaded. Check config file." << endl;
    return;
  }

  cout << "Total entries from files: " << C->GetEntries() << endl;
  
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vz[maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;
  Double_t HCALx, HCALy, HCALe, ekineW2;
  Double_t cblkatime[kNcell];

  C ->SetBranchStatus("*",0);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );


  C->SetBranchAddress("sbs.hcal.x", &HCALx);
  C->SetBranchAddress("sbs.hcal.y", &HCALy);
  C->SetBranchAddress("sbs.hcal.e", &HCALe);
  C->SetBranchAddress("bb.tr.n", &BBtr_n);
  C->SetBranchAddress("bb.tr.px", BBtr_px);
  C->SetBranchAddress("bb.tr.py", BBtr_py);
  C->SetBranchAddress("bb.tr.pz", BBtr_pz);
  C->SetBranchAddress("bb.tr.p", BBtr_p);
  C->SetBranchAddress("bb.tr.vz", BBtr_vz);
  C->SetBranchAddress("bb.ps.e", &BBps_e);
  C->SetBranchAddress("bb.ps.x", &BBps_x);
  C->SetBranchAddress("bb.ps.y", &BBps_y);
  C->SetBranchAddress("bb.sh.e", &BBsh_e);
  C->SetBranchAddress("bb.sh.x", &BBsh_x);
  C->SetBranchAddress("bb.sh.y", &BBsh_y);
  C->SetBranchAddress("e.kine.W2", &ekineW2);
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times

  TFile *fout = new TFile( outputfilename, "RECREATE" );

  TH1D *hdx_HCAL = new TH1D( "hdx_HCAL ", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3 );
  TH1D *hdx_HCAL_final = new TH1D( "hdx_HCAL_final ", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3 );
  TH1D *hdy_HCAL = new TH1D( "hdy_HCAL","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25 );
  TH1D *hW2 = new TH1D( "hW2", " ;GeV2  ", 200,0,5 );
  TH1D *hQ2 = new TH1D( "hQ2", " ;GeV2  ", 200,0,15 );
  TH1D *hHCALe = new TH1D( "hHCALe","HCal Cluster E", 400, 0., 1. );
  TH1D *hhadke = new TH1D( "hhadke","Hadron KE", 400, 0., 10. );
  TH1D *hSampFrac = new TH1D( "hSampFrac","W2 Cut HCal Cluster E / Expected KE", 400, 0., 1. );
  TH2D *hSampFrac_vs_X = new TH2D( "hSampFrac_vs_X","Sampling Fraction ;X (m) ;HCal_E / Exp_KE", xN, HCal_Xi, HCal_Xf, uni_N, 0., 1. );
  TH2D *hSampFrac_vs_Y = new TH2D( "hSampFrac_vs_Y","Sampling Fraction ;Y (m) ;HCal_E / Exp_KE", yN, HCal_Yi, HCal_Yf, uni_N, 0., 1. );
  TH1D *hatime = new TH1D( "hatime", "Cluster adc time, primary block; ns", 220, -50, 170);

  //Long64_t mark = 0;

  Long64_t Nevents = elist ->GetN();
  
  cout << Nevents << " passed global cut." << endl;

  for( Long64_t nevent = 1; nevent < Nevents; nevent++){

    if (nevent%10000==0) cout << "Analyzing entry = " << nevent << " / " << Nevents << " \r";
    cout.flush();
    
    C->GetEntry(elist ->GetEntry(nevent)); 
    
    Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0]);
    Double_t ephi = atan2(BBtr_py[0],BBtr_px[0]);    
    TVector3 vertex(0,0,BBtr_vz[0]);
    TLorentzVector Pbeam(0,0,E_e,E_e);
    TLorentzVector kprime(BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0]);
    TLorentzVector Ptarg(0, 0, 0, M_p);

    TLorentzVector q = Pbeam - kprime;
    TLorentzVector PgammaN = Ptarg +q; //should go through and write this out. Momentum of virtual photon

    Double_t pel = E_e/ (1. +E_e/M_p*(1.-cos(etheta)));//momentum of elastically scattered electron 
    Double_t nu = E_e -BBtr_p[0]; //kinetic energy of the elasticlly scattered electron 
    Double_t pp = sqrt(pow(nu,2)+2 *M_p*nu); 
    Double_t phinucleon = ephi + PI; //coplanar 
    Double_t thetanucleon = acos((E_e - BBtr_pz[0])/pp);
    TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );    
    TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
    TVector3 HCAL_xaxis(0,-1,0);
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
    TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;
    Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis) / (pNhat.Dot( HCAL_zaxis ) );
    TVector3 HCAL_intersect = vertex + sintersect * pNhat; 

    Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis);
    Double_t xexpect_HCAL =  (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis);

    Double_t Q2 = 2*E_e *BBtr_p[0]*( 1-BBtr_pz[0]/BBtr_p[0]);
    Double_t KE_p = nu; //For elastics
    Double_t SFrac = HCALe/KE_p;

    //Double_t W = PgammaN.M();
    //Double_t W2 = pow(W,2);

    Double_t W2 = ekineW2;

    hW2->Fill(W2);
    hQ2->Fill(Q2);

    Double_t dx = HCALx - xexpect_HCAL;
    Double_t dy = HCALy - yexpect_HCAL;

    if (abs(W2 - W2_mean)>W2_sig ) continue; 

    hatime->Fill( cblkatime[0] );
    hdy_HCAL->Fill(dy);
    hdx_HCAL->Fill(dx);

    if(abs(dy-dy0)>dy_sig) continue;

    hdx_HCAL_final->Fill(dx);
    hSampFrac->Fill( SFrac );
    hSampFrac_vs_X->Fill( HCALx, SFrac );
    hSampFrac_vs_Y->Fill( HCALy, SFrac );
    hHCALe->Fill( HCALe );
    hhadke->Fill( HCALe/sampFrac );

  }

  //cout << "mark = " << mark;

  cout << endl;

  fout->Write();
  fout->Close();

}// end main







