//sseeds 11.12.23 - Script to generate quality plots for timing diagnostics

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

//const Double_t atime_fwhm = 5.;
const Double_t time_fwhm = 4.;

bool verb = false;

//Main.
void timing_diagnostic_plots( const char *experiment = "gmn", Int_t kine = 8, Int_t pass = 2 ){
  
  // set style
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);

  // Get information from .csv files for run numbers
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,kine,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available

  vector<caldiag> cut;
  util::ReadDiagnosticCutList(struct_dir,experiment,kine,runs[0].targ,runs[0].mag,verb,cut); //needed only for min ev
  Int_t min_epc = cut[0].min_ev;
  Double_t atime_fwhm = cut[0].atime_sig;

  //get vector of runnumbers
  vector<Double_t> runnumbers;
  for (Int_t irun=0; irun<nruns; ++irun){
    runnumbers.push_back(runs[irun].runnum);
  }

  // outfile path
  std::string out_dir = gSystem->Getenv("OUT_DIR");
  std::string in_path = out_dir + Form("/hcal_calibrations/pass%d/diagnostics/timingdiag_%s_sbs%d_pass%d.root",pass,experiment,kine,pass);

  // open file, get histograms
  TFile *f1 = TFile::Open(in_path.c_str());

  TH2D *hatimerun = (TH2D*)f1->Get("hatimerun");
  TH2D *hatimerow = (TH2D*)f1->Get("hatimerow");
  TH2D *hatimecol = (TH2D*)f1->Get("hatimecol");
  TH2D *hatimeid = (TH2D*)f1->Get("hatimeid");
  TH2D *hatimex = (TH2D*)f1->Get("hatimex");
  TH2D *hatimey = (TH2D*)f1->Get("hatimey");

  TH2D *htimerun = (TH2D*)f1->Get("htimerun");
  TH2D *htimerow = (TH2D*)f1->Get("htimerow");
  TH2D *htimecol = (TH2D*)f1->Get("htimecol");
  TH2D *htimeid = (TH2D*)f1->Get("htimeid");
  TH2D *htimex = (TH2D*)f1->Get("htimex");
  TH2D *htimey = (TH2D*)f1->Get("htimey");

  TH2D *hatimerun_bc = (TH2D*)f1->Get("hatimerun_bc");
  TH2D *hatimerow_bc = (TH2D*)f1->Get("hatimerow_bc");
  TH2D *hatimecol_bc = (TH2D*)f1->Get("hatimecol_bc");
  TH2D *hatimeid_bc = (TH2D*)f1->Get("hatimeid_bc");
  TH2D *hatimex_bc = (TH2D*)f1->Get("hatimex_bc");
  TH2D *hatimey_bc = (TH2D*)f1->Get("hatimey_bc");

  TH2D *htimerun_bc = (TH2D*)f1->Get("htimerun_bc");
  TH2D *htimerow_bc = (TH2D*)f1->Get("htimerow_bc");
  TH2D *htimecol_bc = (TH2D*)f1->Get("htimecol_bc");
  TH2D *htimeid_bc = (TH2D*)f1->Get("htimeid_bc");
  TH2D *htimex_bc = (TH2D*)f1->Get("htimex_bc");
  TH2D *htimey_bc = (TH2D*)f1->Get("htimey_bc");

  // Default plot output directory
  std::string plotdir = Form("../quality_plots/diagnostics/%s/",experiment);

  ////////////////////
  // ADCtime vs run
  TCanvas *c1 = new TCanvas("c1","Time Diff vs Run Number",1600,1200);
  c1->Divide(1,2);
  c1->SetGrid();
    
  c1->cd(1);

  vector<Double_t> c1cell;
  vector<Double_t> c1mean;
  vector<Double_t> c1err;
    
  Int_t c1binsX = hatimerun->GetNbinsX();

  util::sliceHisto(hatimerun, c1binsX, atime_fwhm, min_epc, c1cell, c1mean, c1err );
    
  // Convert vectors to arrays
  Double_t* c1x = &c1cell[0];  
  Double_t* c1y = &c1mean[0];
  Double_t* c1ey = &c1err[0];

  TAxis *xaxis = hatimerun->GetXaxis();

  // Set the bin labels
  for (int i = 1; i <= xaxis->GetNbins(); i++)
    if (i-1 < runnumbers.size())
      xaxis->SetBinLabel(i, std::to_string((int)runnumbers[i-1]).c_str());
  
  hatimerun->GetYaxis()->SetRangeUser(-10.,10.);
  hatimerun->GetYaxis()->SetTitle("ns");
  hatimerun->GetXaxis()->SetTitle("");
  hatimerun->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c1g = new TGraphErrors( c1cell.size(), c1x, c1y, 0, c1ey );
    
  c1g->SetTitle(Form("HCal adctime vs runnumber, SBS%d",kine));
  c1g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c1g->SetMarkerColor(kBlack);
  c1g->SetMarkerSize(1);
  c1g->SetLineColor(kBlack);
  c1g->SetLineWidth(3);
  c1g->Draw("same P");

  c1->Update();
  
  c1->cd(2);

  vector<Double_t> c1cell_bc;
  vector<Double_t> c1mean_bc;
  vector<Double_t> c1err_bc;
    
  Int_t c1binsX_bc = hatimerun_bc->GetNbinsX();

  util::sliceHisto(hatimerun_bc, c1binsX_bc, atime_fwhm, min_epc, c1cell_bc, c1mean_bc, c1err_bc );
    
  //Convert vectors to arrays
  Double_t* c1x_bc = &c1cell_bc[0];  
  Double_t* c1y_bc = &c1mean_bc[0];
  Double_t* c1ey_bc = &c1err_bc[0];

  TAxis *xaxis_bc = hatimerun_bc->GetXaxis();

  // Set the bin labels
  for (int i = 1; i <= xaxis_bc->GetNbins(); i++)
    if (i-1 < runnumbers.size())
      xaxis_bc->SetBinLabel(i, std::to_string((int)runnumbers[i-1]).c_str());

  hatimerun_bc->GetYaxis()->SetRangeUser(-10.,10.);
  hatimerun_bc->GetYaxis()->SetTitle("ns");
  hatimerun_bc->GetXaxis()->SetTitle("");
  hatimerun_bc->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c1g_bc = new TGraphErrors( c1cell_bc.size(), c1x_bc, c1y_bc, 0, c1ey_bc );
    
  c1g_bc->SetTitle(Form("HCal adctime vs runnumber best cluster, SBS%d",kine));
  c1g_bc->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c1g_bc->SetMarkerColor(kBlack);
  c1g_bc->SetMarkerSize(1);
  c1g_bc->SetLineColor(kBlack);
  c1g_bc->SetLineWidth(3);
  c1g_bc->Draw("same P");

  c1->Update();

  std::string hatimerun_path = plotdir + "atimerun.png";
  c1->SaveAs(hatimerun_path.c_str());
 
  ////////////////////
  // ADCtime vs Row
  TCanvas *c2 = new TCanvas("c2","Time Diff vs Row",1600,1200);
  c2->Divide(1,2);
  c2->SetGrid();
    
  c2->cd(1);

  vector<Double_t> c2cell;
  vector<Double_t> c2mean;
  vector<Double_t> c2err;
    
  Int_t c2binsX = hatimerow->GetNbinsX();

  util::sliceHisto(hatimerow, c2binsX, atime_fwhm, min_epc, c2cell, c2mean, c2err );
    
  // Convert vectors to arrays
  Double_t* c2x = &c2cell[0];  
  Double_t* c2y = &c2mean[0];
  Double_t* c2ey = &c2err[0];
  
  hatimerow->GetYaxis()->SetRangeUser(-10.,10.);
  hatimerow->GetYaxis()->SetTitle("ns");
  hatimerow->GetXaxis()->SetTitle("");
  hatimerow->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c2g = new TGraphErrors( c2cell.size(), c2x, c2y, 0, c2ey );
    
  c2g->SetTitle(Form("HCal adctime vs row, SBS%d",kine));
  c2g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c2g->SetMarkerColor(kBlack);
  c2g->SetMarkerSize(1);
  c2g->SetLineColor(kBlack);
  c2g->SetLineWidth(3);
  c2g->Draw("same P");

  c2->Update();
  
  c2->cd(2);

  vector<Double_t> c2cell_bc;
  vector<Double_t> c2mean_bc;
  vector<Double_t> c2err_bc;
    
  Int_t c2binsX_bc = hatimerow_bc->GetNbinsX();

  util::sliceHisto(hatimerow_bc, c2binsX_bc, atime_fwhm, min_epc, c2cell_bc, c2mean_bc, c2err_bc );
    
  //Convert vectors to arrays
  Double_t* c2x_bc = &c2cell_bc[0];  
  Double_t* c2y_bc = &c2mean_bc[0];
  Double_t* c2ey_bc = &c2err_bc[0];

  hatimerow_bc->GetYaxis()->SetRangeUser(-10.,10.);
  hatimerow_bc->GetYaxis()->SetTitle("ns");
  hatimerow_bc->GetXaxis()->SetTitle("");
  hatimerow_bc->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c2g_bc = new TGraphErrors( c2cell_bc.size(), c2x_bc, c2y_bc, 0, c2ey_bc );
    
  c2g_bc->SetTitle(Form("HCal adctime vs row best cluster, SBS%d",kine));
  c2g_bc->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c2g_bc->SetMarkerColor(kBlack);
  c2g_bc->SetMarkerSize(1);
  c2g_bc->SetLineColor(kBlack);
  c2g_bc->SetLineWidth(3);
  c2g_bc->Draw("same P");

  c2->Update();

  std::string hatimerow_path = plotdir + "atimerow.png";
  c2->SaveAs(hatimerow_path.c_str());

 
  ////////////////////
  // ADCtime vs Col
  TCanvas *c3 = new TCanvas("c3","Time Diff vs Col",1600,1200);
  c3->Divide(1,2);
  c3->SetGrid();
    
  c3->cd(1);

  vector<Double_t> c3cell;
  vector<Double_t> c3mean;
  vector<Double_t> c3err;
    
  Int_t c3binsX = hatimecol->GetNbinsX();

  util::sliceHisto(hatimecol, c3binsX, atime_fwhm, min_epc, c3cell, c3mean, c3err );
    
  // Convert vectors to arrays
  Double_t* c3x = &c3cell[0];  
  Double_t* c3y = &c3mean[0];
  Double_t* c3ey = &c3err[0];
  
  hatimecol->GetYaxis()->SetRangeUser(-10.,10.);
  hatimecol->GetYaxis()->SetTitle("ns");
  hatimecol->GetXaxis()->SetTitle("");
  hatimecol->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c3g = new TGraphErrors( c3cell.size(), c3x, c3y, 0, c3ey );
    
  c3g->SetTitle(Form("HCal adctime vs col, SBS%d",kine));
  c3g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c3g->SetMarkerColor(kBlack);
  c3g->SetMarkerSize(1);
  c3g->SetLineColor(kBlack);
  c3g->SetLineWidth(3);
  c3g->Draw("same P");

  c3->Update();
  
  c3->cd(2);

  vector<Double_t> c3cell_bc;
  vector<Double_t> c3mean_bc;
  vector<Double_t> c3err_bc;
    
  Int_t c3binsX_bc = hatimecol_bc->GetNbinsX();

  util::sliceHisto(hatimecol_bc, c3binsX_bc, atime_fwhm, min_epc, c3cell_bc, c3mean_bc, c3err_bc );
    
  //Convert vectors to arrays
  Double_t* c3x_bc = &c3cell_bc[0];  
  Double_t* c3y_bc = &c3mean_bc[0];
  Double_t* c3ey_bc = &c3err_bc[0];

  hatimecol_bc->GetYaxis()->SetRangeUser(-10.,10.);
  hatimecol_bc->GetYaxis()->SetTitle("ns");
  hatimecol_bc->GetXaxis()->SetTitle("");
  hatimecol_bc->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c3g_bc = new TGraphErrors( c3cell_bc.size(), c3x_bc, c3y_bc, 0, c3ey_bc );
    
  c3g_bc->SetTitle(Form("HCal adctime vs col best cluster, SBS%d",kine));
  c3g_bc->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c3g_bc->SetMarkerColor(kBlack);
  c3g_bc->SetMarkerSize(1);
  c3g_bc->SetLineColor(kBlack);
  c3g_bc->SetLineWidth(3);
  c3g_bc->Draw("same P");

  c3->Update();

  std::string hatimecol_path = plotdir + "atimecol.png";
  c3->SaveAs(hatimecol_path.c_str());

  ////////////////////
  // ADCtime vs ID
  TCanvas *c4 = new TCanvas("c4","Time Diff vs Id",1600,1200);
  c4->Divide(1,2);
  c4->SetGrid();
    
  c4->cd(1);

  vector<Double_t> c4cell;
  vector<Double_t> c4mean;
  vector<Double_t> c4err;
    
  Int_t c4binsX = hatimeid->GetNbinsX();

  util::sliceHisto(hatimeid, c4binsX, atime_fwhm, min_epc, c4cell, c4mean, c4err );
    
  // Convert vectors to arrays
  Double_t* c4x = &c4cell[0];  
  Double_t* c4y = &c4mean[0];
  Double_t* c4ey = &c4err[0];
  
  hatimeid->GetYaxis()->SetRangeUser(-10.,10.);
  hatimeid->GetYaxis()->SetTitle("ns");
  hatimeid->GetXaxis()->SetTitle("");
  hatimeid->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c4g = new TGraphErrors( c4cell.size(), c4x, c4y, 0, c4ey );
    
  c4g->SetTitle(Form("HCal adctime vs id, SBS%d",kine));
  c4g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c4g->SetMarkerColor(kBlack);
  c4g->SetMarkerSize(1);
  c4g->SetLineColor(kBlack);
  c4g->SetLineWidth(3);
  c4g->Draw("same P");

  c4->Update();
  
  c4->cd(2);

  vector<Double_t> c4cell_bc;
  vector<Double_t> c4mean_bc;
  vector<Double_t> c4err_bc;
    
  Int_t c4binsX_bc = hatimeid_bc->GetNbinsX();

  util::sliceHisto(hatimeid_bc, c4binsX_bc, atime_fwhm, min_epc, c4cell_bc, c4mean_bc, c4err_bc );
    
  //Convert vectors to arrays
  Double_t* c4x_bc = &c4cell_bc[0];  
  Double_t* c4y_bc = &c4mean_bc[0];
  Double_t* c4ey_bc = &c4err_bc[0];

  hatimeid_bc->GetYaxis()->SetRangeUser(-10.,10.);
  hatimeid_bc->GetYaxis()->SetTitle("ns");
  hatimeid_bc->GetXaxis()->SetTitle("");
  hatimeid_bc->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c4g_bc = new TGraphErrors( c4cell_bc.size(), c4x_bc, c4y_bc, 0, c4ey_bc );
    
  c4g_bc->SetTitle(Form("HCal adctime vs id best cluster, SBS%d",kine));
  c4g_bc->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c4g_bc->SetMarkerColor(kBlack);
  c4g_bc->SetMarkerSize(1);
  c4g_bc->SetLineColor(kBlack);
  c4g_bc->SetLineWidth(3);
  c4g_bc->Draw("same P");

  c4->Update();

  std::string hatimeid_path = plotdir + "atimeid.png";
  c4->SaveAs(hatimeid_path.c_str());

  ////////////////////
  // TDC time vs ID
  TCanvas *c5 = new TCanvas("c5","Time Diff vs Id",1600,1200);
  c5->Divide(1,2);
  c5->SetGrid();
    
  c5->cd(1);

  vector<Double_t> c5cell;
  vector<Double_t> c5mean;
  vector<Double_t> c5err;
    
  Int_t c5binsX = htimeid->GetNbinsX();

  util::sliceHisto(htimeid, c5binsX, time_fwhm, min_epc, c5cell, c5mean, c5err );
    
  // Convert vectors to arrays
  Double_t* c5x = &c5cell[0];  
  Double_t* c5y = &c5mean[0];
  Double_t* c5ey = &c5err[0];
  
  htimeid->GetYaxis()->SetRangeUser(-10.,10.);
  htimeid->GetYaxis()->SetTitle("ns");
  htimeid->GetXaxis()->SetTitle("");
  htimeid->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c5g = new TGraphErrors( c5cell.size(), c5x, c5y, 0, c5ey );
    
  c5g->SetTitle(Form("HCal tdc time vs id, SBS%d",kine));
  c5g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c5g->SetMarkerColor(kBlack);
  c5g->SetMarkerSize(1);
  c5g->SetLineColor(kBlack);
  c5g->SetLineWidth(3);
  c5g->Draw("same P");

  c5->Update();
  
  c5->cd(2);

  vector<Double_t> c5cell_bc;
  vector<Double_t> c5mean_bc;
  vector<Double_t> c5err_bc;
    
  Int_t c5binsX_bc = htimeid_bc->GetNbinsX();

  util::sliceHisto(htimeid_bc, c5binsX_bc, time_fwhm, min_epc, c5cell_bc, c5mean_bc, c5err_bc );
    
  //Convert vectors to arrays
  Double_t* c5x_bc = &c5cell_bc[0];  
  Double_t* c5y_bc = &c5mean_bc[0];
  Double_t* c5ey_bc = &c5err_bc[0];

  htimeid_bc->GetYaxis()->SetRangeUser(-10.,10.);
  htimeid_bc->GetYaxis()->SetTitle("ns");
  htimeid_bc->GetXaxis()->SetTitle("");
  htimeid_bc->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c5g_bc = new TGraphErrors( c5cell_bc.size(), c5x_bc, c5y_bc, 0, c5ey_bc );
    
  c5g_bc->SetTitle(Form("HCal tdc time vs id best cluster, SBS%d",kine));
  c5g_bc->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c5g_bc->SetMarkerColor(kBlack);
  c5g_bc->SetMarkerSize(1);
  c5g_bc->SetLineColor(kBlack);
  c5g_bc->SetLineWidth(3);
  c5g_bc->Draw("same P");

  c5->Update();

  std::string htimeid_path = plotdir + "timeid.png";
  c5->SaveAs(htimeid_path.c_str());

  ////////////////////
  // ADCtime vs X
  TCanvas *c6 = new TCanvas("c6","Time Diff vs X",1600,1200);
  c6->Divide(1,2);
  c6->SetGrid();
    
  c6->cd(1);

  vector<Double_t> c6cell;
  vector<Double_t> c6mean;
  vector<Double_t> c6err;
    
  Int_t c6binsX = hatimex->GetNbinsX();

  util::sliceHisto(hatimex, c6binsX, atime_fwhm, min_epc, c6cell, c6mean, c6err );
    
  // Convert vectors to arrays
  Double_t* c6x = &c6cell[0];  
  Double_t* c6y = &c6mean[0];
  Double_t* c6ey = &c6err[0];
  
  hatimex->GetYaxis()->SetRangeUser(-10.,10.);
  hatimex->GetYaxis()->SetTitle("ns");
  hatimex->GetXaxis()->SetTitle("");
  hatimex->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c6g = new TGraphErrors( c6cell.size(), c6x, c6y, 0, c6ey );
    
  c6g->SetTitle(Form("HCal adctime vs X, SBS%d",kine));
  c6g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c6g->SetMarkerColor(kBlack);
  c6g->SetMarkerSize(1);
  c6g->SetLineColor(kBlack);
  c6g->SetLineWidth(3);
  c6g->Draw("same P");

  c6->Update();
  
  c6->cd(2);

  vector<Double_t> c6cell_bc;
  vector<Double_t> c6mean_bc;
  vector<Double_t> c6err_bc;
    
  Int_t c6binsX_bc = hatimex_bc->GetNbinsX();

  util::sliceHisto(hatimex_bc, c6binsX_bc, atime_fwhm, min_epc, c6cell_bc, c6mean_bc, c6err_bc );
    
  //Convert vectors to arrays
  Double_t* c6x_bc = &c6cell_bc[0];  
  Double_t* c6y_bc = &c6mean_bc[0];
  Double_t* c6ey_bc = &c6err_bc[0];

  hatimex_bc->GetYaxis()->SetRangeUser(-10.,10.);
  hatimex_bc->GetYaxis()->SetTitle("ns");
  hatimex_bc->GetXaxis()->SetTitle("");
  hatimex_bc->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c6g_bc = new TGraphErrors( c6cell_bc.size(), c6x_bc, c6y_bc, 0, c6ey_bc );
    
  c6g_bc->SetTitle(Form("HCal adctime vs X best cluster, SBS%d",kine));
  c6g_bc->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c6g_bc->SetMarkerColor(kBlack);
  c6g_bc->SetMarkerSize(1);
  c6g_bc->SetLineColor(kBlack);
  c6g_bc->SetLineWidth(3);
  c6g_bc->Draw("same P");

  c6->Update();

  std::string hatimex_path = plotdir + "atimex.png";
  c6->SaveAs(hatimex_path.c_str());

  ////////////////////
  // ADCtime vs Y
  TCanvas *c7 = new TCanvas("c7","Time Diff vs Y",1600,1200);
  c7->Divide(1,2);
  c7->SetGrid();
    
  c7->cd(1);

  vector<Double_t> c7cell;
  vector<Double_t> c7mean;
  vector<Double_t> c7err;
    
  Int_t c7binsX = hatimey->GetNbinsX();

  util::sliceHisto(hatimey, c7binsX, atime_fwhm, min_epc, c7cell, c7mean, c7err );
    
  // Convert vectors to arrays
  Double_t* c7x = &c7cell[0];  
  Double_t* c7y = &c7mean[0];
  Double_t* c7ey = &c7err[0];
  
  hatimey->GetYaxis()->SetRangeUser(-10.,10.);
  hatimey->GetYaxis()->SetTitle("ns");
  hatimey->GetXaxis()->SetTitle("");
  hatimey->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c7g = new TGraphErrors( c7cell.size(), c7x, c7y, 0, c7ey );
    
  c7g->SetTitle(Form("HCal adctime vs Y, SBS%d",kine));
  c7g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c7g->SetMarkerColor(kBlack);
  c7g->SetMarkerSize(1);
  c7g->SetLineColor(kBlack);
  c7g->SetLineWidth(3);
  c7g->Draw("same P");

  c7->Update();
  
  c7->cd(2);

  vector<Double_t> c7cell_bc;
  vector<Double_t> c7mean_bc;
  vector<Double_t> c7err_bc;
    
  Int_t c7binsX_bc = hatimey_bc->GetNbinsX();

  util::sliceHisto(hatimey_bc, c7binsX_bc, atime_fwhm, min_epc, c7cell_bc, c7mean_bc, c7err_bc );
    
  //Convert vectors to arrays
  Double_t* c7x_bc = &c7cell_bc[0];  
  Double_t* c7y_bc = &c7mean_bc[0];
  Double_t* c7ey_bc = &c7err_bc[0];

  hatimey_bc->GetYaxis()->SetRangeUser(-10.,10.);
  hatimey_bc->GetYaxis()->SetTitle("ns");
  hatimey_bc->GetXaxis()->SetTitle("");
  hatimey_bc->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c7g_bc = new TGraphErrors( c7cell_bc.size(), c7x_bc, c7y_bc, 0, c7ey_bc );
    
  c7g_bc->SetTitle(Form("HCal adctime vs Y best cluster, SBS%d",kine));
  c7g_bc->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c7g_bc->SetMarkerColor(kBlack);
  c7g_bc->SetMarkerSize(1);
  c7g_bc->SetLineColor(kBlack);
  c7g_bc->SetLineWidth(3);
  c7g_bc->Draw("same P");

  c7->Update();

  std::string hatimey_path = plotdir + "atimey.png";
  c7->SaveAs(hatimey_path.c_str());

}
  
