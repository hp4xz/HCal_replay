//sseeds 11.11.23 - script to plot E and SF comparison results from short calibration run

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

const Double_t SF_hard_llim = 0.02; //Hard lower limit on hcal SF for fitting
const Double_t SF_hard_ulim = 0.30; //Hard lower limit on hcal SF for fitting
const Double_t E_hard_llim = 0.02; //Hard lower limit on hcal E for fitting
const Double_t E_hard_ulim = 1.6; //Hard lower limit on hcal E for fitting (0.8 for all except sbs11: 1.6)

const int calset = 1; //change for kinematics with multiple HV settings

void plot_mc_data_ESF( int kine = 4, int pass = 2, int qr = 0 ){
  
  std::string fout_path = Form("ESF_datamc_k%d_p%d_qr%d_out.root",kine,pass,qr);

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);
  TCanvas *c1 = new TCanvas("c1","HCal Energy MC Comparison",2100,800);
  c1->Divide(2,1);
  c1->cd(1);
  c1->SetGrid();

  //MC volatile location
  string MC_dir = "/lustre19/expphy/volatile/halla/sbs/seeds/hcal_calibrations/MC";
  string MC_ext = MC_dir + Form("/hcalE_idx_mc_dig_gmn_conf%d.root",kine);

  //Open MC file, get energy histogram
  cout << "Opening MC E file here: " << MC_ext << endl;
  TFile *fEmc = TFile::Open(MC_ext.c_str());
  TH1D *hEmc = (TH1D*)fEmc->Get("hE_mc");;

  //get max values of histogram
  Double_t hEmc_max = hEmc->GetMaximum();
  Int_t hEmc_binMax = hEmc->GetMaximumBin();
  Double_t hEmc_binCenter = hEmc->GetBinCenter( hEmc_binMax );
  double fitEmcs = 0.5*hEmc_binCenter; //sig/E ~ 0.5, 2 sigma spread gets fit region
  double fitEmcl = hEmc_binCenter - fitEmcs;
  double fitEmch = hEmc_binCenter + fitEmcs;    

  // Fit a Gaussian to this hist
  TF1 *hEmc_gfit = new TF1("hEmc_gfit", "gaus");
  hEmc->Fit( "hEmc_gfit", "QN", "", fitEmcl, fitEmch ); // "Q" for quiet mode, "N" for no draw
  Double_t hEmc_fit_max = hEmc_gfit->GetMaximum(fitEmcl,fitEmch);
  double Emc_mean = hEmc_gfit->GetParameter(1);
  double Emc_sig = hEmc_gfit->GetParameter(2);

  //data volatile location
  string data_dir = "/lustre19/expphy/volatile/halla/sbs/seeds/hcal_calibrations/pass2/energy";
  string data_ext = data_dir + Form("/ecal_lh2only_gmn_conf%d_qr%d_pass%d.root",kine,qr,pass);

  //open data file, get energy histogram
  cout << "Opening E file here: " << data_ext << endl;
  TFile *fE = TFile::Open(data_ext.c_str());
  TH1D *hE = (TH1D*)fE->Get(Form("hE_set%d",calset));

  //get max value of fit
  Int_t hE_binMax = hE->GetMaximumBin();
  Double_t hE_binCenter = hE->GetBinCenter( hE_binMax );
  double fitEs = 0.5*hE_binCenter; //sig/E ~ 0.5, 2 sigma spread gets fit region
  Double_t fitEl = hE_binCenter - fitEs;
  Double_t fitEh = hE_binCenter + fitEs;

  // Fit a Gaussian to this hist
  TF1 *hE_gfit = new TF1( "hE_gfit", "gaus" );
  hE->Fit( "hE_gfit", "QN", "", fitEl, fitEh ); // "Q" for quiet mode, "N" for no draw
  Double_t hE_fit_max = hE_gfit->GetMaximum(fitEl,fitEh);
  double E_mean = hE_gfit->GetParameter(1);
  double E_sig = hE_gfit->GetParameter(2);

  //get data to MC ratios
  double data_mc_mean_ratio = E_mean/Emc_mean;
  double data_mc_sig_ratio = E_sig/Emc_sig;
  
  //get scale factor, scale, and get max value
  Double_t scale_factor = hEmc_fit_max / hE_fit_max;
  hE->Scale(scale_factor);
  Double_t hE_max = hE->GetMaximum();

  Double_t yrangescale = hEmc_max*1.1;
  if( hEmc_max < hE_max )
    yrangescale = hE_max*1.1;

  hEmc->SetLineWidth(1);
  hEmc->SetLineColor(kGreen);
  hEmc->SetFillStyle(3004);
  hEmc->SetFillColor(kGreen);
  hEmc->SetTitle(Form("HCal Energy (pass%d, qr%d)",pass,qr));
  hEmc->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hEmc->GetYaxis()->SetRangeUser(0,yrangescale);
  hEmc->GetXaxis()->SetTitle("GeV");
  hEmc->Draw("E");

  hE->SetLineWidth(1);
  hE->SetLineColor(kBlack);
  hE->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hE->GetXaxis()->SetTitle("GeV");
  hE->Draw("hist same");

  //Add a legend
  auto legend = new TLegend(0.63,0.7,0.89,0.89);
  legend->SetTextSize(0.03);
  legend->SetHeader(Form("SBS%d, Cal Set%d",kine,calset));
  legend->AddEntry(hEmc,"MC","l");
  legend->AddEntry(hE,"Data (Scaled)","l");
  legend->AddEntry( (TObject*)0, "", "");
  legend->AddEntry( (TObject*)0, Form("Data/MC mean ratio: %0.2f",data_mc_mean_ratio), "");
  //legend->AddEntry( (TObject*)0, Form("Data/MC sig ratio: %0.2f",data_mc_sig_ratio), "");
  legend->Draw();

  c1->Modified();

  //Now do sampling fraction
  c1->cd(2);

  //get SF mc plot
  TH1D *hSFmc = (TH1D*)fEmc->Get("hSF_mc");;

  //get max values of histogram
  Double_t hSFmc_max = hSFmc->GetMaximum();
  Int_t hSFmc_binMax = hSFmc->GetMaximumBin();
  Double_t hSFmc_binCenter = hSFmc->GetBinCenter( hSFmc_binMax );
  double fitSFmcs = 0.5*hSFmc_binCenter; //sig/E ~ 0.5, 2 sigma spread gets fit region
  double fitSFmcl = hSFmc_binCenter - fitSFmcs;
  double fitSFmch = hSFmc_binCenter + fitSFmcs;    

  // Fit a Gaussian to this hist
  TF1 *hSFmc_gfit = new TF1("hSFmc_gfit", "gaus");
  hSFmc->Fit( "hSFmc_gfit", "QN", "", fitSFmcl, fitSFmch ); // "Q" for quiet mode, "N" for no draw
  Double_t hSFmc_fit_max = hSFmc_gfit->GetMaximum(fitSFmcl,fitSFmch);
  double SFmc_mean = hSFmc_gfit->GetParameter(1);
  double SFmc_sig = hSFmc_gfit->GetParameter(2);
  
  //get SF data plot
  TH1D *hSF = (TH1D*)fE->Get(Form("hSF_set%d",calset));

  //get max value of fit
  Int_t hSF_binMax = hSF->GetMaximumBin();
  Double_t hSF_binCenter = hSF->GetBinCenter( hSF_binMax );
  double fitSFs = 0.5*hSF_binCenter; //sig/E ~ 0.5, 2 sigma spread gets fit region
  Double_t fitSFl = hSF_binCenter - fitSFs;
  Double_t fitSFh = hSF_binCenter + fitSFs;
  
  // Fit a Gaussian to this hist
  TF1 *hSF_gfit = new TF1( "hSF_gfit", "gaus" );
  hSF->Fit( "hSF_gfit", "QN", "", fitSFl, fitSFh ); // "Q" for quiet mode, "N" for no draw
  Double_t hSF_fit_max = hSF_gfit->GetMaximum(fitSFl,fitSFh);
  double SF_mean = hSF_gfit->GetParameter(1);
  double SF_sig = hSF_gfit->GetParameter(2);
  
  //get data to MC ratios
  double data_mc_mean_ratio_SF = SF_mean/SFmc_mean;
  double data_mc_sig_ratio_SF = SF_sig/SFmc_sig;

  //get scale factor, scale, and get max value
  Double_t scale_factor_SF = hSFmc_fit_max / hSF_fit_max;
  hSF->Scale(scale_factor_SF);
  Double_t hSF_max = hSF->GetMaximum();

  Double_t yrangescale_SF = hSFmc_max*1.1;
  if( hSFmc_max < hSF_max )
    yrangescale_SF = hSF_max*1.1;

  hSFmc->SetLineWidth(1);
  hSFmc->SetLineColor(kGreen);
  hSFmc->SetFillStyle(3004);
  hSFmc->SetFillColor(kGreen);
  hSFmc->SetTitle(Form("HCal Sampling Fraction (pass%d, qr%d)",pass,qr));
  hSFmc->GetXaxis()->SetRangeUser(0,SF_hard_ulim);
  hSFmc->GetYaxis()->SetRangeUser(0,yrangescale_SF);
  hSFmc->GetXaxis()->SetTitle("GeV");
  hSFmc->Draw("E");

  hSF->SetLineWidth(2);
  hSF->SetLineColor(kBlack);
  hSF->GetXaxis()->SetRangeUser(0,SF_hard_ulim);
  hSF->GetXaxis()->SetTitle("GeV");
  hSF->Draw("hist same");

  //Add a legend
  auto legend_SF = new TLegend(0.63,0.7,0.89,0.89);
  legend_SF->SetTextSize(0.03);
  legend_SF->SetHeader(Form("SBS%d, Cal Set%d",kine,calset));
  legend_SF->AddEntry(hSFmc,"MC","l");
  legend_SF->AddEntry(hSF,"Data (Scaled)","l");
  legend_SF->AddEntry( (TObject*)0, Form("Data/MC mean ratio: %0.2f",data_mc_mean_ratio_SF), "");
  //legend_SF->AddEntry( (TObject*)0, Form("Data/MC sig ratio: %0.2f",data_mc_sig_ratio_SF), "");
  legend_SF->Draw();

  // c2->Modified();
  // c2->Write();

  c1->Update();
  c1->Write();

  fout->Write();

  cout << "Analysis written to " << fout_path << endl;


}
