//sseeds 8.2.23 script to compare output of qreplay hcal energy distributions against MC

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

//const Double_t E_hard_ulim = 0.9;

//const Double_t E_hard_ulim = 1.95; //GMn kine 11
//const Double_t E_hard_ulim = 1.35; //GMn kine 7
const Double_t E_hard_ulim = 0.95; //all GMn kine other than 11/7
const Double_t SF_hard_ulim = 0.39;
const Double_t E_approx_FWHM = 0.3;
const Double_t SFvrowcol_yrange = 0.12;

void resizeAndSaveCanvas(const char* filename, const char* canvasName, int config) {
    // Open the TFile
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Get the canvas
    TCanvas* canvas = (TCanvas*)file->Get(canvasName);
    if (!canvas) {
        std::cerr << "Error: Canvas " << canvasName << " not found in " << filename << std::endl;
        file->Close();
        return;
    }

    // Resize the canvas
    canvas->SetCanvasSize(1600, 500);  // For example, set size to 800x600

    // Save the canvas as a .png
    canvas->SaveAs(Form("EvRun_8to%d.png",config));

    // Clean up
    file->Close();
}

void overlayWithGaussianFits(TH2D* hist, TCanvas* canvas, const std::vector<int>& redBins) {
  if (!hist || !canvas) {
    std::cerr << "Null histogram or canvas pointer!" << std::endl;
    return;
  }

  // Use the provided canvas
  canvas->cd();

  // Create TGraphErrors for standard and red points
  TGraphErrors* graph = new TGraphErrors();
  TGraphErrors* redGraph = new TGraphErrors();

  // //get x value of first bin
  // Int_t xstart = hist->GetBinCenter(1);
  // //

  for (int binX = 1; binX <= hist->GetNbinsX(); ++binX) {
    // Project the 2D histogram onto a 1D histogram along the y-axis for the current x-bin
    TH1D* projY = hist->ProjectionY("_py", binX, binX);

    //Skip the fit if less than 100 entries exist in the bin
    if( projY->GetEntries()<100 )
      continue;
	
    //declare some dynamic fit variables
    Double_t FWHM = E_approx_FWHM;
    Int_t binMax = projY->GetMaximumBin();
    Double_t binCenter = projY->GetBinCenter( binMax );
    Double_t fitLowerLim = binCenter - FWHM;
    Double_t fitUpperLim = binCenter + FWHM;
	
    // Fit a Gaussian to this projection
    TF1 *gaussFit = new TF1("gausFit", "gaus");
    projY->Fit(gaussFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode
    
    if (gaussFit) { // Ensure the fit was successful

      // Determine the x-value as the center of the current bin
      double xCenter = hist->GetXaxis()->GetBinCenter(binX);

      double mean = gaussFit->GetParameter(1);
      double sigma = gaussFit->GetParameter(2);

      // Get the errors on the mean and sigma
      double meanError = gaussFit->GetParError(1);
      double sigmaError = gaussFit->GetParError(2);

      // Determine which graph to fill based on whether the bin is in redBins
      

      TGraphErrors* currentGraph = graph;
      Int_t binXval = hist->GetXaxis()->GetBinCenter(binX);

      //cout << "binXval: " << binXval << endl;

      if (std::find(redBins.begin(), redBins.end(), binXval) != redBins.end()) {
	currentGraph = redGraph;
      }

      int pointIdx = currentGraph->GetN();
      currentGraph->SetPoint(pointIdx, xCenter, mean);
      currentGraph->SetPointError(pointIdx, 0, sigma);
    }
    delete projY; // Clean up
  }

  // Draw the original histogram, standard graph, and red graph
  hist->Draw("COLZ");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("P SAME");

  redGraph->SetMarkerStyle(20);
  redGraph->SetMarkerColor(kRed);
  redGraph->SetLineColor(kRed);
  redGraph->Draw("P SAME");

  // Update the canvas
  canvas->Update();

  //Add a legend
  auto legend = new TLegend(0.63,0.7,0.89,0.89);
  legend->SetTextSize(0.03);
  legend->AddEntry(graph,"LD2 Runs","l");
  legend->AddEntry(redGraph,"LH2 Runs","l");
  legend->Draw();

  // Update the canvas
  canvas->Update();
}

void SFoverlayWithGaussianFits(TH2D* hist, TCanvas* canvas, Double_t ymax, string axis) {
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

    //Skip the fit if less than 100 entries exist in the bin
    if( projY->GetEntries()<100 )
      continue;
	
    //declare some dynamic fit variables
    Double_t FWHM = E_approx_FWHM;
    Int_t binMax = projY->GetMaximumBin();
    Double_t binCenter = projY->GetBinCenter( binMax );
    Double_t fitLowerLim = binCenter - FWHM;
    Double_t fitUpperLim = binCenter + FWHM;
	
    // Fit a Gaussian to this projection
    TF1 *gaussFit = new TF1("gausFit", "gaus");
    projY->Fit(gaussFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode
    
    if (gaussFit) { // Ensure the fit was successful

      // Determine the x-value as the center of the current bin
      double xCenter = hist->GetXaxis()->GetBinCenter(binX);

      double mean = gaussFit->GetParameter(1);
      double sigma = gaussFit->GetParameter(2);

      // Get the errors on the mean and sigma
      double meanError = gaussFit->GetParError(1);
      double sigmaError = gaussFit->GetParError(2);

      TGraphErrors* currentGraph = graph;

      int pointIdx = currentGraph->GetN();
      currentGraph->SetPoint(pointIdx, xCenter, mean);
      currentGraph->SetPointError(pointIdx, 0, sigma);
    }
    delete projY; // Clean up
  }

  // Draw the original histogram, standard graph, and red graph
  hist->GetYaxis()->SetRangeUser(0,ymax);
  hist->GetYaxis()->SetTitle("E_{clus}/T_{N}");
  hist->GetXaxis()->SetTitle(axis.c_str());
  hist->Draw("COLZ");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("P SAME");

  // Update the canvas
  canvas->Update();

  // //Add a legend
  // auto legend = new TLegend(0.63,0.7,0.89,0.89);
  // legend->SetTextSize(0.03);
  // legend->AddEntry(graph,"LD2 Runs","l");
  // legend->AddEntry(redGraph,"LH2 Runs","l");
  // legend->Draw();

  // Update the canvas
  //canvas->Update();
}


//Compares the energy spectra from qreplay_standalone.C root file and mc data
void hcale_sbs8_compare( const char *experiment = "gmn", Int_t config = 4, Int_t pass = 0, Double_t FWHM = 0.03, bool lh2only = true ){

  //hardcoded for my environment for now
  std::string outdir_path = "/lustre19/expphy/volatile/halla/sbs/seeds";

  std::string lh2_marker = "";
  if( lh2only )
    lh2_marker = "_lh2only";

  //set up relative paths
  std::string plotdir = Form("../quality_plots/%s/conf%d%s/",experiment,config,"_s4comp");
  std::string qr_path = outdir_path + Form("/hcal_calibrations/qreplay/qreplay%s_from_gmn_4_0_to_%s_%d_%d.root",lh2_marker.c_str(),experiment,config,pass); //only looking at sbs8 compared with each other kinematic's full data set (LH2+LD2)
  std::string digmc_path = outdir_path + Form("/hcal_calibrations/MC/hcalE_idx_mc_dig_%s_conf%d.root",experiment,config);
  std::string energy_dmc_path = plotdir + "qr_energy_dmc.png";
  std::string sf_dmc_path = plotdir + "qr_sf_dmc.png";
  std::string erun_path = plotdir + "qr_e_run.png";
  std::string sf_x_path = plotdir + "qr_sf_x.png";
  std::string sf_y_path = plotdir + "qr_sf_y.png";

  std::string fout_path = outdir_path + Form("/hcal_calibrations/qreplay/ecomp_from_gmn_4_0_to_%s_%d_%d.root",experiment,config,pass);

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  TCanvas *c1 = new TCanvas("c1","mc/data energy",1600,1200);
  c1->cd();
  gStyle->SetOptStat(0);

  TFile *qr_file = new TFile(qr_path.c_str(), "READ"); // Open the ROOT file in read mode


  ////
  //Get E and SF comparison plots

  TH1D *hE = (TH1D*)qr_file->Get("hE_new");
  Int_t hE_entries = hE->GetEntries();

  //get max value of fit to peak for hE
  Int_t hE_binMax = hE->GetMaximumBin();
  Double_t hE_binCenter = hE->GetBinCenter( hE_binMax );
  Double_t hE_fitLL = hE_binCenter - FWHM;
  Double_t hE_fitUL = hE_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hE_gfit = new TF1("hE_gfit", "gaus");
  hE->Fit("hE_gfit", "QN", "", hE_fitLL, hE_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hE_fit_max = hE_gfit->GetMaximum(hE_fitLL,hE_fitUL);

  //cout << "hE fit max " << hE_fit_max << endl;

  TFile *digmc_file = new TFile(digmc_path.c_str(), "READ"); // Open the ROOT file in read mode
  TH1D *hEmc = (TH1D*)digmc_file->Get("hE_mc");
  Int_t hEmc_entries = hEmc->GetEntries();
  TH1D *hEmc_nocut =  (TH1D*)digmc_file->Get("hE_mc_nocut");
  Int_t hEmc_nocut_entries = hEmc_nocut->GetEntries();

  //get max value of fit to peak for hEmc
  Double_t hEmc_max = hEmc->GetMaximum();

  Int_t hEmc_binMax = hEmc->GetMaximumBin();
  Double_t hEmc_binCenter = hEmc->GetBinCenter( hEmc_binMax );
  Double_t hEmc_fitLL = hEmc_binCenter - FWHM;
  Double_t hEmc_fitUL = hEmc_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hEmc_gfit = new TF1("hEmc_gfit", "gaus");
  hEmc->Fit("hEmc_gfit", "QN", "", hEmc_fitLL, hEmc_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hEmc_fit_max = hEmc_gfit->GetMaximum(hEmc_fitLL,hEmc_fitUL);

  //get max value of fit to peak for hEmc_nocut
  Double_t hEmc_nocut_max = hEmc_nocut->GetMaximum();

  Int_t hEmc_nocut_binMax = hEmc_nocut->GetMaximumBin();
  Double_t hEmc_nocut_binCenter = hEmc_nocut->GetBinCenter( hEmc_nocut_binMax );
  Double_t hEmc_nocut_fitLL = hEmc_nocut_binCenter - FWHM;
  Double_t hEmc_nocut_fitUL = hEmc_nocut_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hEmc_nocut_gfit = new TF1("hEmc_nocut_gfit", "gaus");
  hEmc_nocut->Fit("hEmc_nocut_gfit", "QN", "", hEmc_nocut_fitLL, hEmc_nocut_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hEmc_nocut_fit_max = hEmc_nocut_gfit->GetMaximum(hEmc_nocut_fitLL,hEmc_nocut_fitUL);

  // Obtain scale factor and scale other histograms
  Double_t escale_factor = hEmc_fit_max / hE_fit_max;
  hE->Scale(escale_factor);
  Double_t escale_nocut_factor = hEmc_fit_max / hEmc_nocut_fit_max;
  hEmc_nocut->Scale(escale_nocut_factor);
  Double_t hE_max = hE->GetMaximum();

  //cout << hEmc_max << ":" << hE_max << endl;

  Double_t eyrangescale = hEmc_max*1.1;
  if( hEmc_max < hE_max )
    eyrangescale = hE_max*1.1;

  hEmc->SetLineWidth(1);
  hEmc->SetLineColor(kGreen);
  hEmc->SetFillStyle(3004);
  hEmc->SetFillColor(kGreen);
  hEmc->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hEmc->GetYaxis()->SetRangeUser(0,eyrangescale);
  hEmc->GetXaxis()->SetTitle("GeV");
  hEmc->Draw("hist");

  hEmc_nocut->SetLineWidth(1);
  hEmc_nocut->SetLineColor(kCyan);
  hEmc_nocut->SetFillStyle(3005);
  hEmc_nocut->SetFillColor(kCyan);
  hEmc_nocut->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hEmc_nocut->GetYaxis()->SetRangeUser(0,eyrangescale);
  hEmc_nocut->GetXaxis()->SetTitle("GeV");
  hEmc_nocut->Draw("hist same");
    
  hE->SetLineWidth(2);
  hE->SetLineColor(kBlack);
  hE->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hE->GetXaxis()->SetTitle("GeV");
  hE->Draw("hist same");
  
  //Add a legend
  auto elegend = new TLegend(0.33,0.7,0.89,0.89);
  elegend->SetTextSize(0.03);
  elegend->SetHeader(Form("SBS%d, E",config));
  elegend->AddEntry(hEmc,Form("MC, %d entries",hEmc_entries),"l");
  elegend->AddEntry(hEmc_nocut,Form("MC no cuts (scale: %0.2f), %d entries",escale_nocut_factor,hEmc_nocut_entries),"l");
  elegend->AddEntry(hE,Form("Data (scale: %0.2f), %d entries",escale_factor,hE_entries),"l");
  elegend->Draw();

  c1->SaveAs(energy_dmc_path.c_str());
  //c1->Write();

  //Sampling fraction
  TCanvas *c2 = new TCanvas("c2","mc/data sampling fraction",1600,1200);
  c2->cd();
  gStyle->SetOptStat(0);

  //get mc sf
  //TFile *fSFmc = TFile::Open(digmc_path.c_str());
  TH1D *hSFmc = (TH1D*)digmc_file->Get("hSF_mc");
  Int_t hSFmc_entries = hSFmc->GetEntries();
  TH1D *hSFmc_nocut = (TH1D*)digmc_file->Get("hSF_mc_nocut");
  Int_t hSFmc_nocut_entries = hSFmc_nocut->GetEntries();

  //get max value of hSFmc histogram
  Double_t hSFmc_max = hSFmc->GetMaximum();

  Int_t hSFmc_binMax = hSFmc->GetMaximumBin();
  Double_t hSFmc_binCenter = hSFmc->GetBinCenter( hSFmc_binMax );
  Double_t hSFmc_fitLL = hSFmc_binCenter - FWHM;
  Double_t hSFmc_fitUL = hSFmc_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hSFmc_gfit = new TF1("hSFmc_gfit", "gaus");
  hSFmc->Fit("hSFmc_gfit", "QN", "", hSFmc_fitLL, hSFmc_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hSFmc_fit_max = hSFmc_gfit->GetMaximum(hSFmc_fitLL,hSFmc_fitUL);

  //get max value of hSFmc_nocut histogram
  Double_t hSFmc_nocut_max = hSFmc_nocut->GetMaximum();

  Int_t hSFmc_nocut_binMax = hSFmc_nocut->GetMaximumBin();
  Double_t hSFmc_nocut_binCenter = hSFmc_nocut->GetBinCenter( hSFmc_nocut_binMax );
  Double_t hSFmc_nocut_fitLL = hSFmc_nocut_binCenter - FWHM;
  Double_t hSFmc_nocut_fitUL = hSFmc_nocut_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hSFmc_nocut_gfit = new TF1("hSFmc_nocut_gfit", "gaus");
  hSFmc_nocut->Fit("hSFmc_nocut_gfit", "QN", "", hSFmc_nocut_fitLL, hSFmc_nocut_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hSFmc_nocut_fit_max = hSFmc_nocut_gfit->GetMaximum(hSFmc_nocut_fitLL,hSFmc_nocut_fitUL);

  //get data sf
  //TFile *qr_file = new TFile(qr_path.c_str(), "READ"); // Open the ROOT file in read mode
  TH1D *hSF = (TH1D*)qr_file->Get("hSF_new");
  Int_t hSF_entries = hSF->GetEntries();

  //get max value of fit
  Int_t hSF_binMax = hSF->GetMaximumBin();
  Double_t hSF_binCenter = hSF->GetBinCenter( hSF_binMax );
  Double_t hSF_fitLL = hSF_binCenter - FWHM;
  Double_t hSF_fitUL = hSF_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hSF_gfit = new TF1("hSF_gfit", "gaus");
  hSF->Fit("hSF_gfit", "QN", "", hSF_fitLL, hSF_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hSF_fit_max = hSF_gfit->GetMaximum(hSF_fitLL,hSF_fitUL);

  //get scale factor
  Double_t sfscale_factor = hSFmc_fit_max / hSF_fit_max;
  hSF->Scale(sfscale_factor);
  Double_t sfscale_nocut_factor = hSFmc_fit_max / hSFmc_nocut_fit_max;
  hSFmc_nocut->Scale(sfscale_nocut_factor);
  Double_t hSF_max = hSF->GetMaximum();

  Double_t sfyrangescale = hSFmc_max*1.1;
  if( hSFmc_max < hSF_max )
    sfyrangescale = hSF_max*1.1;

  hSFmc->SetLineWidth(1);
  hSFmc->SetLineColor(kOrange);
  hSFmc->SetFillStyle(3004);
  hSFmc->SetFillColor(kOrange);
  hSFmc->GetXaxis()->SetRangeUser(0,SF_hard_ulim);
  hSFmc->GetYaxis()->SetRangeUser(0,sfyrangescale);
  hSFmc->GetXaxis()->SetTitle("E_{clus}/KE");
  hSFmc->Draw("hist");

  hSFmc_nocut->SetLineWidth(1);
  hSFmc_nocut->SetLineColor(kRed);
  hSFmc_nocut->SetFillStyle(3005);
  hSFmc_nocut->SetFillColor(kRed);
  hSFmc_nocut->GetXaxis()->SetRangeUser(0,SF_hard_ulim);
  hSFmc_nocut->GetYaxis()->SetRangeUser(0,sfyrangescale);
  hSFmc_nocut->GetXaxis()->SetTitle("E_{clus}/KE");
  hSFmc_nocut->Draw("hist same");

  hSF->SetLineWidth(2);
  hSF->SetLineColor(kBlack);
  hSF->GetXaxis()->SetRangeUser(0,SF_hard_ulim);
  hSF->GetXaxis()->SetTitle("E_{clus}/KE");
  hSF->Draw("hist same");

  //Add a legend
  auto sflegend = new TLegend(0.33,0.7,0.89,0.89);
  sflegend->SetTextSize(0.03);
  sflegend->SetHeader(Form("SBS%d, SF",config));
  sflegend->AddEntry(hSFmc,Form("MC, %d entries",hSFmc_entries),"l");
  sflegend->AddEntry(hSFmc_nocut,Form("MC no cuts (scale: %0.2f), %d entries",sfscale_nocut_factor,hSFmc_nocut_entries),"l");
  sflegend->AddEntry(hSF,Form("Data (scale: %0.2f), %d entries",sfscale_factor,hSF_entries),"l");
  //sflegend->AddEntry(hSFmc,"MC","l");
  //sflegend->AddEntry(hSF,"Data (Scaled)","l");
  sflegend->Draw();

  c2->SaveAs(sf_dmc_path.c_str());
  //c2->Write();

  ////
  //Get canvases with E vs run, resize, and print
  //get list of lh2runs in vector
  std::string struct_dir = Form("../config/%s/","gmn");
  Int_t nruns = -1;
  Int_t verb = 0;
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,"gmn",nruns,config,pass,verb,runs);
  std::vector<Int_t> lh2runs;

  for (Int_t r=0; r<nruns; r++) {
    Int_t current_runnumber = runs[r].runnum;
    std::string current_target = runs[r].target;

    if( current_target.compare("lh2")==0 )
      lh2runs.push_back(current_runnumber);
  }  

  std::string qrall_path = outdir_path + Form("/hcal_calibrations/qreplay/qreplay_from_gmn_4_0_to_%s_%d_%d.root",experiment,config,pass); //only looking at sbs8 compared with each other kinematic's full data set (LH2+LD2)
  TFile *qrall_file = new TFile(qrall_path.c_str(), "READ"); // Open the ROOT file in read mode

  //resizeAndSaveCanvas( qrall_path.c_str(), "c1", config );
  TH2D *hErun = (TH2D*)qrall_file->Get("hE_run");
  
  TCanvas *c3 = new TCanvas("c3","E vs Run, LH2 and LD2",1600,600);
  c3->cd();

  overlayWithGaussianFits(hErun,c3,lh2runs);

  c3->SaveAs(erun_path.c_str());

  ////
  //Get SF vs row and SF vs col and overlay

  TH2D *hSFvrow = (TH2D*)qr_file->Get("hSFvrow");
  //hSFvrow->GetYaxis()->SetRangeUser(0,0.12);

  TCanvas *c4 = new TCanvas("c3","E vs X",1600,600);
  c4->cd();

  SFoverlayWithGaussianFits(hSFvrow,c4,SFvrowcol_yrange,"row");

  c4->SaveAs(sf_x_path.c_str());

  TH2D *hSFvcol = (TH2D*)qr_file->Get("hSFvcol");
  //hSFvcol->GetYaxis()->SetRangeUser(0,0.12);

  TCanvas *c5 = new TCanvas("c3","E vs Y",1600,600);
  c5->cd();

  SFoverlayWithGaussianFits(hSFvcol,c5,SFvrowcol_yrange,"col");

  c5->SaveAs(sf_y_path.c_str());


  fout->Write();

}
