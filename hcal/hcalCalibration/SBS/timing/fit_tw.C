//sseeds 9.23.23 macro to fit various timewalk plots from output of timewalk.C

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
#include <fstream>
#include <string>

//y and x axis fit limits
const Double_t adct_ymin = -8.;
const Double_t adct_ymax = 9.;
const Double_t tdc_ymin = -3.;
const Double_t tdc_ymax = 7.;
const Double_t tw_fit_llim = 0.01;
const Double_t tw_fit_ulim = 0.50;

// expo timewalk fit parameters
const Double_t tw_p0_llim = 6.0;
const Double_t tw_p0_ulim = 7.5;
const Double_t tw_p1_llim = 9.;
const Double_t tw_p1_ulim = 11.;
const Double_t tw_p2_llim = -4;
const Double_t tw_p2_ulim = 4;

// traditional fit tw fit
// const Double_t trad_p0_llim = -2.1;
// const Double_t trad_p0_ulim = -1.9;
// const Double_t trad_p1_llim = 1.0;
// const Double_t trad_p1_ulim = 1.4;
// const Double_t trad_p2_llim = 0.5;
// const Double_t trad_p2_ulim = 0.5;

//pol1 tdc tw fit parameters
// const Double_t p0_set = 4.;
// const Double_t p1_set = -18.;
// const Double_t p0_ll = 3.9;
// const Double_t p0_ul = 4.1;
// const Double_t p1_ll = -19.;
// const Double_t p1_ul = -16.;
const Double_t fit_ul = 0.50;
const Double_t fit_ll = 0.01;

//pol1 adct tw fit parameters
// const Double_t ap0_set = 2.0;
// const Double_t ap1_set = -5.2;
// const Double_t ap0_ll = 1.;
// const Double_t ap0_ul = 3.;
// const Double_t ap1_ll = -7.;
// const Double_t ap1_ul = -4.;
const Double_t afit_ul = 0.50;
const Double_t afit_ll = 0.01;

// Usage:
// TH2D* myClonedHist = CloneTH2DWithinYRange(originalHist, yMinValue, yMaxValue);

void TradFitTH2DWithinYRange(TH2D* origHist, double yMin, double yMax, vector<double> &fitParams) {
    // Create a new TH2D with the same binning in X, but limited Y range
    TH2D* clonedHist = new TH2D("clonedHist", origHist->GetTitle(),
                                origHist->GetXaxis()->GetNbins(),
                                origHist->GetXaxis()->GetXmin(),
                                origHist->GetXaxis()->GetXmax(),
                                origHist->GetYaxis()->FindBin(yMax) - origHist->GetYaxis()->FindBin(yMin) + 1,
                                yMin, yMax);

    // Loop over the X and Y bins of the original TH2D
    for (int i = 1; i <= origHist->GetXaxis()->GetNbins(); i++) {
        for (int j = origHist->GetYaxis()->FindBin(yMin); j <= origHist->GetYaxis()->FindBin(yMax); j++) {
            double binContent = origHist->GetBinContent(i, j);
            clonedHist->SetBinContent(i, j - origHist->GetYaxis()->FindBin(yMin) + 1, binContent);
        }
    }

    //Fit the TDC vs E plot with traditional timewalk fit
    TF1 *fittdc_trad = new TF1( "fittdc_trad", util::g_tradtwfit, tw_fit_llim, tw_fit_ulim, 3 );
    
    // fittdc_trad->SetParameters(( trad_p0_ulim + trad_p0_llim )/2,
    // 			       ( trad_p1_ulim + trad_p1_llim )/2,
    // 			       ( trad_p2_ulim + trad_p2_llim )/2);

    fittdc_trad->FixParameter(2,0.5); //obtain better results for majority of data with this par fixed

    clonedHist->Fit("fittdc_trad","MQRN0","",tw_fit_llim,tw_fit_ulim);

    fitParams.push_back(fittdc_trad->GetParameter(0));
    fitParams.push_back(fittdc_trad->GetParameter(1));
    fitParams.push_back(fittdc_trad->GetParameter(2));

}

void ExpoFitTH2DWithinYRange(TH2D* origHist, double yMin, double yMax, vector<double> &fitParams) {
    // Create a new TH2D with the same binning in X, but limited Y range
    TH2D* clonedHist = new TH2D("clonedHist", origHist->GetTitle(),
                                origHist->GetXaxis()->GetNbins(),
                                origHist->GetXaxis()->GetXmin(),
                                origHist->GetXaxis()->GetXmax(),
                                origHist->GetYaxis()->FindBin(yMax) - origHist->GetYaxis()->FindBin(yMin) + 1,
                                yMin, yMax);

    // Loop over the X and Y bins of the original TH2D
    for (int i = 1; i <= origHist->GetXaxis()->GetNbins(); i++) {
        for (int j = origHist->GetYaxis()->FindBin(yMin); j <= origHist->GetYaxis()->FindBin(yMax); j++) {
            double binContent = origHist->GetBinContent(i, j);
            clonedHist->SetBinContent(i, j - origHist->GetYaxis()->FindBin(yMin) + 1, binContent);
        }
    }

    //Fit the TDC vs E plot with exponential fit
    TF1 *fittdc_expo = new TF1( "fittdc_expo", util::g_twfit, tw_fit_llim, tw_fit_ulim, 3 );

    fittdc_expo->SetParameters((tw_p0_ulim+tw_p0_llim)/2,
    			       (tw_p1_ulim+tw_p1_llim)/2,
    			       (tw_p2_ulim+tw_p2_llim)/2);
  
    // fittdc_expo->SetParLimits(0,tw_p0_llim,tw_p0_ulim);
    // fittdc_expo->SetParLimits(1,tw_p1_llim,tw_p1_ulim);
    fittdc_expo->SetParLimits(2,tw_p2_llim,tw_p2_ulim);
  
    clonedHist->Fit("fittdc_expo","MQRN0","",tw_fit_llim,tw_fit_ulim);

    fitParams.push_back(fittdc_expo->GetParameter(0));
    fitParams.push_back(fittdc_expo->GetParameter(1));
    fitParams.push_back(fittdc_expo->GetParameter(2));

}

void tdcPolFitTH2DWithinYRange(TH2D* origHist, double yMin, double yMax, vector<double> &fitParams) {
    // Create a new TH2D with the same binning in X, but limited Y range
    TH2D* clonedHist = new TH2D("clonedHist", origHist->GetTitle(),
                                origHist->GetXaxis()->GetNbins(),
                                origHist->GetXaxis()->GetXmin(),
                                origHist->GetXaxis()->GetXmax(),
                                origHist->GetYaxis()->FindBin(yMax) - origHist->GetYaxis()->FindBin(yMin) + 1,
                                yMin, yMax);

    // Loop over the X and Y bins of the original TH2D
    for (int i = 1; i <= origHist->GetXaxis()->GetNbins(); i++) {
        for (int j = origHist->GetYaxis()->FindBin(yMin); j <= origHist->GetYaxis()->FindBin(yMax); j++) {
            double binContent = origHist->GetBinContent(i, j);
            clonedHist->SetBinContent(i, j - origHist->GetYaxis()->FindBin(yMin) + 1, binContent);
        }
    }

    //Fit the TDC vs E plots with pol1
    TF1 *fittdc_pol = new TF1( "fittdc_pol", util::g_lfit, fit_ll, fit_ul, 2 );
    //fittdc_pol->SetParameters(p0_set,p1_set);
  
    // fittdc_pol->SetParLimits(0,p0_ll,p0_ul);
    // fittdc_pol->SetParLimits(1,p1_ll,p1_ul);

    clonedHist->Fit("fittdc_pol","MQRN0","",fit_ll,fit_ul);

    fitParams.push_back(fittdc_pol->GetParameter(0));
    fitParams.push_back(fittdc_pol->GetParameter(1));

}

void adctPolFitTH2DWithinYRange(TH2D* origHist, double yMin, double yMax, vector<double> &fitParams) {
    // Create a new TH2D with the same binning in X, but limited Y range
    TH2D* clonedHist = new TH2D("clonedHist", origHist->GetTitle(),
                                origHist->GetXaxis()->GetNbins(),
                                origHist->GetXaxis()->GetXmin(),
                                origHist->GetXaxis()->GetXmax(),
                                origHist->GetYaxis()->FindBin(yMax) - origHist->GetYaxis()->FindBin(yMin) + 1,
                                yMin, yMax);

    // Loop over the X and Y bins of the original TH2D
    for (int i = 1; i <= origHist->GetXaxis()->GetNbins(); i++) {
        for (int j = origHist->GetYaxis()->FindBin(yMin); j <= origHist->GetYaxis()->FindBin(yMax); j++) {
            double binContent = origHist->GetBinContent(i, j);
            clonedHist->SetBinContent(i, j - origHist->GetYaxis()->FindBin(yMin) + 1, binContent);
        }
    }

    //Fit the ADCt vs E plots with pol1
    TF1 *fitadct_pol = new TF1( "fitadct_pol", util::g_lfit, afit_ll, afit_ul, 2 );
    //fitadct_pol->SetParameters(ap0_set,ap1_set);

    // fitadct_pol->SetParLimits(0,ap0_ll,ap0_ul);
    // fitadct_pol->SetParLimits(1,ap1_ll,ap1_ul);

    clonedHist->Fit("fitadct_pol","MQRN0","",afit_ll,afit_ul);

    fitParams.push_back(fitadct_pol->GetParameter(0));
    fitParams.push_back(fitadct_pol->GetParameter(1));

    fitadct_pol->Delete();
    clonedHist->Delete();

}


// Usage:
// replaceLineWithPrefix("path/to/file.txt", "prefixToSearch", "replacementLine");
void replaceLineWithPrefix(const std::string& filePath, const std::string& searchString, const std::string& replacementString) {
    std::ifstream inFile(filePath);  // Open file for reading
    if (!inFile.is_open()) {
        std::cerr << "Unable to open file for reading: " << filePath << std::endl;
        return;
    }

    std::string tempFilePath = filePath + ".tmp";
    std::ofstream outFile(tempFilePath);  // Create a temporary file for writing
    if (!outFile.is_open()) {
        std::cerr << "Unable to open temp file for writing: " << tempFilePath << std::endl;
        return;
    }

    std::string line;
    while (getline(inFile, line)) {
        if (line.find(searchString) == 0) {  // If the line starts with the searchString
            outFile << replacementString << std::endl;  // Write the replacement string to the temp file
        } else {
            outFile << line << std::endl;  // Write the original line to the temp file
        }
    }

    inFile.close();
    outFile.close();

    // Replace the original file with the temporary file
    if (std::rename(tempFilePath.c_str(), filePath.c_str()) != 0) {
        std::cerr << "Error renaming temp file to original." << std::endl;
        std::remove(tempFilePath.c_str());  // Cleanup: delete the temp file
    }
}


//main
void fit_tw( string experiment = "gmn", Int_t config = 9, bool overwrite_parameters = true, Int_t epass = -1 ){

  int pass = 0;
  if( experiment.compare("gmn")==0 && ( config==11 || config==14 || config==8 || config==9 ) )
    pass = 1;
  if( epass!=-1 )
    pass = epass;

  //if overwrite parameters is chosen, add safety check
  std::string response = "n";

  if( overwrite_parameters ){
    // Prompt the user for input
    std::cout << "User has chosen to overwrite staged timewalk calibrations. Do you wish to proceed? (yes/no): ";
    std::cin >> response;

    // Convert the response to lowercase
    std::transform(response.begin(), response.end(), response.begin(), ::tolower);

    // Check the user's response
    if (response != "yes" && response != "y") {
      std::cout << "Exiting the program." << std::endl;
      return 1;
    }
  }

  //hardcoded for my environment for now
  std::string outdir_path = "/lustre19/expphy/volatile/halla/sbs/seeds";

  //paths to output plots
  std::string plotdir = Form("../quality_plots/%s/conf%d%s/",experiment.c_str(),config,"_comp");
  std::string twfit_tdctrad_path = plotdir + "twfit_tdctrad.png";
  std::string twfit_tdcexpo_path = plotdir + "twfit_tdcexpo.png";
  std::string twfit_tdcpol_path = plotdir + "twfit_tdcpol.png";
  std::string twfit_adctpol_path = plotdir + "twfit_adctpol.png";

  //path to timewalk.C output file
  std::string tw_plot_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/timewalk_%s_conf%d_qr0_pass%d.root",pass,experiment.c_str(),config,pass); //only looking at sbs8 compared with each other kinematic's full data set (LH2+LD2)

  //path to output analysis file
  std::string fout_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/timewalk_fits_%s_%d_%d.root",pass,experiment.c_str(),config,pass);

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  cout << "output file path: " << fout_path << endl;

  //path to timewalk parameters
  std::string new_tw_path = Form("../timing/parameters/timewalk_%s_conf%d_pass%d.txt",experiment.c_str(),config,pass);

  //Begin analysis by clearing the fit overlay from the TH2Ds from the original analysis file
  TFile *tw_plot_file = new TFile(tw_plot_path.c_str(), "READ"); // Open the ROOT file in read mode
  TH2D *h1 = (TH2D*)tw_plot_file->Get("htdc_tradclone");
  TH2D *htdc_trad = (TH2D*)h1->Clone("htdc_trad");
  htdc_trad->GetListOfFunctions()->Clear();

  TH2D *h2 = (TH2D*)tw_plot_file->Get("htdc_expoclone");
  TH2D *htdc_expo = (TH2D*)h2->Clone("htdc_expo");
  htdc_expo->GetListOfFunctions()->Clear();

  TH2D *h3 = (TH2D*)tw_plot_file->Get("htdc_polclone");
  TH2D *htdc_pol = (TH2D*)h3->Clone("htdc_pol");
  htdc_pol->GetListOfFunctions()->Clear();

  TH2D *h4 = (TH2D*)tw_plot_file->Get("hadct_polclone");
  TH2D *hadct_pol = (TH2D*)h4->Clone("hadct_pol");
  hadct_pol->GetListOfFunctions()->Clear();

  vector<double> trad_tdcfit_params;
  vector<double> expo_tdcfit_params;
  vector<double> linear_tdcfit_params;
  vector<double> linear_adctfit_params;

  // Declare canvases
  TCanvas *c1 = new TCanvas("c1","TDC_tradfit",1600,1200);
  TCanvas *c2 = new TCanvas("c2","TDC_expofit",1600,1200);
  TCanvas *c3 = new TCanvas("c3","TDC_linearfit",1600,1200);

  gStyle->SetPalette(53);

  ////
  //traditional tdc fit

  std::cout << "working on traditional tdc fit.." << std::endl;
  
  c1->cd();

  //Fit the TDC vs E plot with traditional timewalk fit
  TradFitTH2DWithinYRange(htdc_trad,tdc_ymin,tdc_ymax,trad_tdcfit_params);

  TF1 *f1 = new TF1("f1", util::g_tradtwfit, tw_fit_llim, tw_fit_ulim, 3);

  for( int i=0; i<trad_tdcfit_params.size(); i++)
    f1->FixParameter(i,trad_tdcfit_params[i]);

  htdc_trad->SetTitle(Form("tdc tradfit P0:%0.2f P1:%0.2f P2:%0.2f",trad_tdcfit_params[0],trad_tdcfit_params[1],trad_tdcfit_params[2]));
  htdc_trad->Draw("colz");
  f1->SetLineColor(kGreen);
  f1->SetLineWidth(3);
  f1->Draw("same");

  gPad->Update();

  c1->Write();
  c1->SaveAs(twfit_tdctrad_path.c_str());

  ////
  //exponential tdc fit

  c2->cd();

  cout << "working on exponential tdc fit.." << endl;

  //Fit the TDC vs E plot with traditional timewalk fit
  ExpoFitTH2DWithinYRange(htdc_expo,tdc_ymin,tdc_ymax,expo_tdcfit_params);

  TF1 *f2 = new TF1("f2", util::g_twfit, tw_fit_llim, tw_fit_ulim, 3 );

  for( int i=0; i<expo_tdcfit_params.size(); i++)
    f2->FixParameter(i,expo_tdcfit_params[i]);

  htdc_expo->SetTitle(Form("tdc expofit P0:%0.2f P1:%0.2f P2:%0.2f",expo_tdcfit_params[0],expo_tdcfit_params[1],expo_tdcfit_params[2]));
  htdc_expo->Draw("colz");
  f2->SetLineColor(kGreen);
  f2->SetLineWidth(3);
  f2->Draw("same");

  gPad->Update();
  
  c2->Write();
  c2->SaveAs(twfit_tdcexpo_path.c_str());

  ////
  //linear tdc fit

  c3->cd();
  
  cout << "working on linear tdc fit.." << endl;
  
  //Fit the TDC vs E plot with traditional timewalk fit
  tdcPolFitTH2DWithinYRange(htdc_pol,tdc_ymin,tdc_ymax,linear_tdcfit_params);

  TF1 *f3 = new TF1("f3", util::g_lfit, fit_ll, fit_ul, 2 );

  for( int i=0; i<linear_tdcfit_params.size(); i++)
    f3->FixParameter(i,linear_tdcfit_params[i]);

  htdc_pol->SetTitle(Form("tdc polfit P0:%0.2f P1:%0.2f",linear_tdcfit_params[0],linear_tdcfit_params[1]));
  htdc_pol->Draw("colz");
  f3->SetLineColor(kGreen);
  f3->SetLineWidth(3);
  f3->Draw("same");

  gPad->Update();

  c3->Write();
  c3->SaveAs(twfit_tdcpol_path.c_str());

  ////
  //linear adct fit

  cout << "working on linear adct fit.." << endl;

  //Fit the TDC vs E plot with traditional timewalk fit
  adctPolFitTH2DWithinYRange(hadct_pol,adct_ymin,adct_ymax,linear_adctfit_params);

  TCanvas *c4 = new TCanvas("c4","ADCt_linearfit",1600,1200);
  c4->cd();

  TF1 *f4 = new TF1("f4", util::g_lfit, afit_ll, afit_ul, 2 );

  for( int i=0; i<linear_adctfit_params.size(); i++)
    f4->FixParameter(i,linear_adctfit_params[i]);

  hadct_pol->SetTitle(Form("adct polfit P0:%0.2f P1:%0.2f",linear_adctfit_params[0],linear_adctfit_params[1]));
  hadct_pol->Draw("colz");
  f4->SetLineColor(kGreen);
  f4->SetLineWidth(3);
  f4->Draw("same");

  gPad->Update();

  c4->Write();
  c4->SaveAs(twfit_adctpol_path.c_str());

  // Create timewalk fit parameter files
  std::string db_tdctw_trad_variable = "sbs.hcal.tdc.tradtw"; //a+b/E^n fit tdc v E, p0=b, p1=n, a asymptotes to signal
  std::string db_tdctw_expo_variable = "sbs.hcal.tdc.expotw"; //expo fit tdc v E, p0, p1
  std::string db_tdctw_pol_variable = "sbs.hcal.tdc.poltw"; //pol1 fit slope tdc v E
  std::string db_adcttw_pol_variable = "sbs.hcal.adct.poltw"; //pol1 fit slope adct v E

  ofstream tw;
  if( !overwrite_parameters )
    new_tw_path = "/dev/null"; //Safety to prevent overwriting constants on quasi-replay

  if( overwrite_parameters ){

    //overwrite tdc traditional fit line
    std::string new_trad_param_line = db_tdctw_trad_variable + Form(" = %0.2f %0.2f",trad_tdcfit_params[1],trad_tdcfit_params[2]);
    replaceLineWithPrefix( new_tw_path.c_str(), db_tdctw_trad_variable.c_str(), new_trad_param_line.c_str());
    
    //overwrite tdc exponential fit line
    std::string new_expo_param_line = db_tdctw_expo_variable + Form(" = %0.2f %0.2f",expo_tdcfit_params[0],expo_tdcfit_params[1]);
    replaceLineWithPrefix( new_tw_path.c_str(), db_tdctw_expo_variable.c_str(), new_expo_param_line.c_str());
    
    //overwrite tdc linear fit line
    std::string new_tdcpol_param_line = db_tdctw_pol_variable + Form(" = %0.2f",linear_tdcfit_params[1]);
    replaceLineWithPrefix( new_tw_path.c_str(), db_tdctw_pol_variable.c_str(), new_tdcpol_param_line.c_str());
    
    //overwrite adt linear fit line
    std::string new_adctpol_param_line = db_adcttw_pol_variable + Form(" = %0.2f",linear_adctfit_params[1]);
    replaceLineWithPrefix( new_tw_path.c_str(), db_adcttw_pol_variable.c_str(), new_adctpol_param_line.c_str());
    
    tw.close();

  }

  fout->Write();

}
