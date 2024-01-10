//sseeds 9.1.23 script to read in two adct or tdc vs run histos and plot together with tgraph overlay

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

const Int_t tfitmin = 50;
const Double_t fitl = 0.;
const Double_t fitu = 0.2;

//main
void timevrun_patchwork( const char *type = "adct", int kine = 14 ) {
  
  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);

  // Open the ROOT files
  string file_dir = "/lustre19/expphy/volatile/halla/sbs/seeds/hcal_calibrations/pass1/timing";
  string f1_path = file_dir + Form("/%salign_gmn_conf%d_qr1_pass1_normalset.root";
  string f2_path = file_dir + Form("/%salign_gmn_conf%d_qr1_pass1_wonkset.root";

  TFile *f1 = TFile::Open(f1_path.c_str());
  TFile *f2 = TFile::Open(f2_path.c_str());

  // Get the TH2D histograms
  TH2D *h1 = (TH2D*)f1->Get("name_of_hist1");  // Replace with the correct name of your histogram
  TH2D *h2 = (TH2D*)f2->Get("name_of_hist2");  // Replace with the correct name of your histogram

  // Create a new histogram to store the sum
  TH2D *hsum = (TH2D*)h1->Clone("hsum");
  hsum->Add(h2);

  // Create a canvas for drawing
  TCanvas *c1 = new TCanvas("c1", Form("%s vs run, sbs%d",type,kine), 1600, 600);
  c1->cd();

  // Plot the combined TH2D
  hsum->Draw("COLZ"); // Use COLZ option to draw a color plot

  // Create vectors to store the means and standard deviations
  vector<double> bin_centers, means, std_devs, std_devs_errors;

  // Fit each bin in x with a Gaussian and extract mean and std deviation
  for (int i = 1; i <= hsum->GetNbinsX(); i++) {
    TH1D *proj = hsum->ProjectionY("_py", i, i); // Project to 1D histogram along y for the bin
    TF1 *gauss = new TF1("gauss", "gaus", proj->GetMean() - 2*proj->GetStdDev(), proj->GetMean() + 2*proj->GetStdDev());
    proj->Fit(gauss, "Q"); // Quiet mode
        
    bin_centers.push_back(hsum->GetXaxis()->GetBinCenter(i));
    means.push_back(gauss->GetParameter(1)); // Mean (parameter 1 of the Gaussian)
    std_devs.push_back(gauss->GetParameter(2)); // Std deviation (parameter 2 of the Gaussian)
    std_devs_errors.push_back(gauss->GetParError(2)); // Error on std deviation
  }

  // Convert vectors to arrays (needed for TGraphErrors)
  double *x = &bin_centers[0];
  double *y = &means[0];
  double *yerrors = &std_devs[0];
  double *xerrors = new double[bin_centers.size()](); // Assuming no error on x-axis

  // Create the TGraphErrors
  TGraphErrors *graph = new TGraphErrors(bin_centers.size(), x, y, xerrors, yerrors);

  // Plot the TGraphErrors on the same canvas
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kRed);
  graph->Draw("P SAME");  // Draws on the same canvas

  c1->SaveAs("result.png");  // Save the canvas as a png image

  // Cleanup
  delete c1;
  f1->Close();
  f2->Close();
}

