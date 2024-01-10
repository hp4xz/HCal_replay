//sseeds 9.23.23 macro to fit various timewalk plots

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

const Double_t E_hard_ulim = 0.9;

//SBS4

// // expo timewalk fit parameters - may need to be tuned at different hadron momenta
// const Double_t tw_fit_llim = 0.01;
// const Double_t tw_fit_ulim = 0.50;
// const Double_t tw_p0_llim = 5.5;
// const Double_t tw_p0_ulim = 6.5;
// const Double_t tw_p1_llim = 9.;
// const Double_t tw_p1_ulim = 10.;
// const Double_t tw_p2_llim = -1.3;
// const Double_t tw_p2_ulim = -1.1;

// // traditional fit tw fit
// const Double_t trad_p0_llim = -2.1;
// const Double_t trad_p0_ulim = -1.9;
// const Double_t trad_p1_llim = 1.0;
// const Double_t trad_p1_ulim = 1.4;
// const Double_t trad_p2_llim = 0.5;
// const Double_t trad_p2_ulim = 0.5;

// //pol1 tdc tw fit parameters
// const Double_t p0_set = 4.;
// const Double_t p1_set = -18.;
// const Double_t p0_ll = 3.9;
// const Double_t p0_ul = 4.1;
// const Double_t p1_ll = -19.;
// const Double_t p1_ul = -16.;
// const Double_t fit_ul = 0.50;
// const Double_t fit_ll = 0.01;

// //pol1 adct tw fit parameters
// const Double_t ap0_set = 2.0;
// const Double_t ap1_set = -5.2;
// const Double_t ap0_ll = 1.;
// const Double_t ap0_ul = 3.;
// const Double_t ap1_ll = -7.;
// const Double_t ap1_ul = -4.;
// const Double_t afit_ul = 0.50;
// const Double_t afit_ll = 0.01;

//SBS7
const Double_t adct_ymin = -10.;
const Double_t adct_ymax = 10.;
// const Double_t tdc_ymin = -10.;
// const Double_t tdc_ymax = 10.;
const Double_t tdc_ymin = -3.;
const Double_t tdc_ymax = 7.;
const Double_t tw_fit_llim = 0.01;
const Double_t tw_fit_ulim = 0.50;

// expo timewalk fit parameters - may need to be tuned at different hadron momenta
const Double_t tw_p0_llim = 6.0;
const Double_t tw_p0_ulim = 7.5;
const Double_t tw_p1_llim = 9.;
const Double_t tw_p1_ulim = 11.;
const Double_t tw_p2_llim = -4;
const Double_t tw_p2_ulim = 4;

// traditional fit tw fit
const Double_t trad_p0_llim = -2.1;
const Double_t trad_p0_ulim = -1.9;
const Double_t trad_p1_llim = 1.0;
const Double_t trad_p1_ulim = 1.4;
const Double_t trad_p2_llim = 0.5;
const Double_t trad_p2_ulim = 0.5;

//pol1 tdc tw fit parameters
const Double_t p0_set = 4.;
const Double_t p1_set = -18.;
const Double_t p0_ll = 3.9;
const Double_t p0_ul = 4.1;
const Double_t p1_ll = -19.;
const Double_t p1_ul = -16.;
const Double_t fit_ul = 0.50;
const Double_t fit_ll = 0.01;

//pol1 adct tw fit parameters
const Double_t ap0_set = 2.0;
const Double_t ap1_set = -5.2;
const Double_t ap0_ll = 1.;
const Double_t ap0_ul = 3.;
const Double_t ap1_ll = -7.;
const Double_t ap1_ul = -4.;
const Double_t afit_ul = 0.50;
const Double_t afit_ll = 0.01;


#include <iostream>
#include <fstream>
#include <string>


// Usage:
// TH2D* myClonedHist = CloneTH2DWithinYRange(originalHist, yMinValue, yMaxValue);

TH2D* CloneTH2DWithinYRange(TH2D* origHist, double yMin, double yMax) {
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

    return clonedHist;
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
void fit_tw( string experiment = "gmn", Int_t config = 11, bool overwrite_parameters = false, Int_t epass = -1 ){

  int pass = 0;
  if( experiment.compare("gmn")==0 && ( config==11 || config==14 || config==8 || config==9 ) )
    pass = 1;
  if( epass!=-1 )
    pass = epass;

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

  //Begin analysis

  TFile *tw_plot_file = new TFile(tw_plot_path.c_str(), "READ"); // Open the ROOT file in read mode
  TH2D *htdc_trad = (TH2D*)tw_plot_file->Get("htdc_tradclone");
  TH2D *htdc_expo = (TH2D*)tw_plot_file->Get("htdc_expoclone");
  TH2D *htdc_pol = (TH2D*)tw_plot_file->Get("htdc_polclone");
  TH2D *hadct_pol = (TH2D*)tw_plot_file->Get("hadct_polclone");

  // clone into new histograms with better fitting range in y
  TH2D *htdc_tradclone = CloneTH2DWithinYRange( htdc_trad, tdc_ymin, tdc_ymax );
  TH2D *htdc_expoclone = CloneTH2DWithinYRange( htdc_expo, tdc_ymin, tdc_ymax );
  TH2D *htdc_polclone = CloneTH2DWithinYRange( htdc_pol, tdc_ymin, tdc_ymax );
  TH2D *hadct_polclone = CloneTH2DWithinYRange( hadct_pol, adct_ymin, adct_ymax );

  //TH2D *htdc_tradplotclone = (TH2D*) htdc_trad->Clone("htdc_tradplotclone");


  // TFile *tw_plot_file = new TFile(tw_plot_path.c_str(), "READ"); // Open the ROOT file in read mode
  // TH2D *htdc_tradclone = (TH2D*)tw_plot_file->Get("htdc_tradclone");
  // TH2D *htdc_expoclone = (TH2D*)tw_plot_file->Get("htdc_expoclone");
  // TH2D *htdc_polclone = (TH2D*)tw_plot_file->Get("htdc_polclone");
  // TH2D *hadct_polclone = (TH2D*)tw_plot_file->Get("hadct_polclone");


  // set up arrays for timewalk fits
  Double_t tdcvE_tradP0 = 0.;
  Double_t tdcvE_tradP1 = 0.;
  Double_t tdcvE_tradP2 = 0.;
  Double_t tdcvE_expoP0 = 0.;
  Double_t tdcvE_expoP1 = 0.;
  Double_t tdcvE_expoP2 = 0.;
  Double_t tdcvE_pol = 0.;
  Double_t tdcvE_pol_offset = 0.;
  Double_t adctvE_pol = 0.;
  Double_t adctvE_pol_offset = 0.;

  // Declare canvases
  TCanvas *c1 = new TCanvas("c1","TDC_tradfit",1600,1200);
  TCanvas *c2 = new TCanvas("c2","TDC_expofit",1600,1200);
  TCanvas *c3 = new TCanvas("c3","TDC_linearfit",1600,1200);
  TCanvas *c4 = new TCanvas("c4","ADCt_linearfit",1600,1200);

  gStyle->SetPalette(53);

  ////
  //traditional tdc fit

  std::cout << "working on traditional tdc fit.." << std::endl;
  
  c1->cd();

  //Fit the TDC vs E plot with traditional timewalk fit
  TF1 *fittdc_trad = new TF1( "fittdc_trad", util::g_tradtwfit, tw_fit_llim, tw_fit_ulim, 3 );

  fittdc_trad->SetParameters(( trad_p0_ulim + trad_p0_llim )/2,
  			     ( trad_p1_ulim + trad_p1_llim )/2,
  			     ( trad_p2_ulim + trad_p2_llim )/2);

  // fittdc_trad->SetParLimits(0,trad_p0_llim,trad_p0_ulim);
  // fittdc_trad->SetParLimits(1,trad_p1_llim,trad_p1_ulim);
  // fittdc_trad->SetParLimits(2,trad_p2_llim,trad_p2_ulim);

  //fittdc_trad->SetParLimits(2,0.3,1.8);

  fittdc_trad->FixParameter(2,0.5); //obtain better results for majority of data with this par fixed

  htdc_tradclone->Fit("fittdc_trad","MQRN0","",tw_fit_llim,tw_fit_ulim);

  tdcvE_tradP0 = fittdc_trad->GetParameter(0);
  tdcvE_tradP1 = fittdc_trad->GetParameter(1);
  tdcvE_tradP2 = fittdc_trad->GetParameter(2);

  htdc_tradclone->SetTitle(Form("tdc tradfit P0:%0.2f P1:%0.2f P2:%0.2f",tdcvE_tradP0,tdcvE_tradP1,tdcvE_tradP2));
  htdc_tradclone->Draw("colz");

  // TF1 *fittdc_tradclone = new TF1( "fittdc_trad", util::g_tradtwfit, tw_fit_llim, tw_fit_ulim, 3 );

  // fittdc_trad->FixParameter(0,tdcvE_tradP0);
  // fittdc_trad->FixParameter(1,tdcvE_tradP1);
  // fittdc_trad->FixParameter(2,tdcvE_tradP2);

  // htdc_tradplotclone->Draw("colz");
  //fittdc_tradclone->Draw("same");

  c1->Write();
  c1->SaveAs(twfit_tdctrad_path.c_str());

  ////
  //exponential tdc fit

  c2->cd();

  cout << "working on exponential tdc fit.." << endl;
  
  //Fit the TDC vs E plot with exponential fit
  TF1 *fittdc_expo = new TF1( "fittdc_expo", util::g_twfit, tw_fit_llim, tw_fit_ulim, 3 );

  fittdc_expo->SetParameters((tw_p0_ulim+tw_p0_llim)/2,
			     (tw_p1_ulim+tw_p1_llim)/2,
			     (tw_p2_ulim+tw_p2_llim)/2);
  
  // fittdc_expo->SetParLimits(0,tw_p0_llim,tw_p0_ulim);
  // fittdc_expo->SetParLimits(1,tw_p1_llim,tw_p1_ulim);
  fittdc_expo->SetParLimits(2,tw_p2_llim,tw_p2_ulim);
  
  htdc_expoclone->Fit("fittdc_expo","MQR","",tw_fit_llim,tw_fit_ulim);
  
  tdcvE_expoP0 = fittdc_expo->GetParameter(0);
  tdcvE_expoP1 = fittdc_expo->GetParameter(1);
  
  htdc_expoclone->SetTitle(Form("tdc expofit P0:%0.2f P1:%0.2f P2:%0.2f",tdcvE_expoP0,tdcvE_expoP1,fittdc_expo->GetParameter(2)));
  htdc_expoclone->Draw("colz");
  
  c2->Write();
  c2->SaveAs(twfit_tdcexpo_path.c_str());

  ////
  //linear tdc fit

  c3->cd();
  
  cout << "working on linear tdc fit.." << endl;
  
  //Fit the TDC vs E plots with pol1
  TF1 *fittdc_pol = new TF1( "fittdc_pol", util::g_lfit, fit_ll, fit_ul, 2 );
  fittdc_pol->SetParameters(p0_set,p1_set);
  
  // fittdc_pol->SetParLimits(0,p0_ll,p0_ul);
  // fittdc_pol->SetParLimits(1,p1_ll,p1_ul);

  htdc_polclone->Fit("fittdc_pol","MQR","",fit_ll,fit_ul);
  tdcvE_pol = fittdc_pol->GetParameter(1);

  htdc_polclone->SetTitle(Form("tdc linearfit P0:%0.2f P1:%0.2f",fittdc_pol->GetParameter(0),tdcvE_pol));

  htdc_polclone->Draw("colz");

  c3->Write();
  c3->SaveAs(twfit_tdcpol_path.c_str());

  ////
  //linear adct fit

  c4->cd();

  cout << "working on linear adct fit.." << endl;

  //Fit the ADCt vs E plots with pol1
  TF1 *fitadct_pol = new TF1( "fitadct_pol", util::g_lfit, afit_ll, afit_ul, 2 );
  fitadct_pol->SetParameters(ap0_set,ap1_set);

  // fitadct_pol->SetParLimits(0,ap0_ll,ap0_ul);
  // fitadct_pol->SetParLimits(1,ap1_ll,ap1_ul);

  hadct_polclone->Fit("fitadct_pol","MQR","",afit_ll,afit_ul);
  adctvE_pol = fitadct_pol->GetParameter(1);

  hadct_polclone->SetTitle(Form("adct linearfit P0:%0.2f P1:%0.2f",fitadct_pol->GetParameter(0),adctvE_pol));
  hadct_polclone->Draw("colz");

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

    // tw.open( new_tw_path );

    //overwrite tdc traditional fit line
    std::string new_trad_param_line = db_tdctw_trad_variable + Form(" = %0.1f %0.1f",tdcvE_tradP1,tdcvE_tradP2);
    replaceLineWithPrefix( new_tw_path.c_str(), db_tdctw_trad_variable.c_str(), new_trad_param_line.c_str());
    
    //overwrite tdc exponential fit line
    std::string new_expo_param_line = db_tdctw_expo_variable + Form(" = %0.1f %0.1f",tdcvE_expoP0,tdcvE_expoP1);
    replaceLineWithPrefix( new_tw_path.c_str(), db_tdctw_expo_variable.c_str(), new_expo_param_line.c_str());
    
    //overwrite tdc linear fit line
    std::string new_tdcpol_param_line = db_tdctw_pol_variable + Form(" = %0.1f",tdcvE_pol);
    replaceLineWithPrefix( new_tw_path.c_str(), db_tdctw_pol_variable.c_str(), new_tdcpol_param_line.c_str());
    
    //overwrite traditional fit line
    std::string new_adctpol_param_line = db_adcttw_pol_variable + Form(" = %0.1f",adctvE_pol);
    replaceLineWithPrefix( new_tw_path.c_str(), db_adcttw_pol_variable.c_str(), new_adctpol_param_line.c_str());



    // //check if timestamp for calibration set is newer than config timestamp. Use the newer of the two.
    // std::string active_timestamp;
    // util::tsCompare(sbs_timestamp,tw_cal[twset].timestamp,active_timestamp);
        
    // tw.open( new_tw_path );

    // tw << endl << endl << "#HCal tdc vs E fit parameters obtained " << date.c_str() << endl;

    // tw << endl << active_timestamp << endl << endl;

    // tw << "#Traditional tdc vs E fit for all PMTs -> y = P0 + P1/( x^P2 ). P0 normal to signal." << endl;
	      
    // tw << db_tdctw_trad_variable << " = ";
    // tw << tdcvE_tradP1 << " ";
    // tw << tdcvE_tradP2 << endl << endl;

    // tw << "#Exponential tdc vs E fit for all PMTs -> y = P0*exp(-P1*x) + P2. P2 normal to signal." << endl;
	      
    // tw << db_tdctw_expo_variable << " = ";
    // tw << tdcvE_expoP0 << " ";
    // tw << tdcvE_expoP1 << endl << endl;

    // tw << endl << "#First order polynomial tdc vs E fit for all PMTs -> P0 + P1*x. P0 normal to signal." << endl;
      
    // tw << db_tdctw_pol_variable << " = ";
    // tw << tdcvE_pol << endl << endl;

    // tw << endl << "#First order polynomial adc time vs E fit for all PMTs -> P0 + P1*x. P0 normal to signal." << endl;
      
    // tw << db_adcttw_pol_variable << " = ";
    // tw << adctvE_pol << endl << endl;
    
    // tw << endl << endl; //needed for now to prevent bug in util::readDB
    
    tw.close();

  }

  fout->Write();

}
