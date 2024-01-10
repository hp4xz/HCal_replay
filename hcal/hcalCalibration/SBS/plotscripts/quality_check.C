//sseeds 5.23.23 - Script to generate supplemental quality plots for timing alignments and energy calibrations

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

const Double_t hcalE_sigma_factor = 0.05; //Rough correspondance between energy spectra width and beam energy (beamE * hcalE_sigma_factor = hcalE distribution sigma. Used for fit limits only.
const Double_t Y_range_factor = 20; //plot +/- about average value for comparison. Adjust depending on variance across HCal ID
const Double_t SF_hard_llim = 0.02; //Hard lower limit on hcal SF for fitting
const Double_t SF_hard_ulim = 0.30; //Hard lower limit on hcal SF for fitting
const Double_t E_hard_llim = 0.02; //Hard lower limit on hcal E for fitting
const Double_t E_hard_ulim = 1.6; //Hard lower limit on hcal E for fitting (0.8 for all except sbs11: 1.6)
const Double_t tdc_fit_llim = -120.;
const Double_t tdc_fit_ulim = 50.;
const Int_t fit_min = 50;
const Double_t FWHM = 0.05; //approx fullwidthhalfmax hcal E dist
const Double_t adct_FWHM = 5; //approx fullwidthhalfmax hcal adc time

void fitAggHistograms(const char* file1, const char* file2, const char* adct_fitcomp_path, int set ) {
    // Open the ROOT files and get the histograms
    TFile f1(file1, "READ");
    TFile f2(file2, "READ");
    
    TH1D* hist1 = (TH1D*) f1.Get(Form("hadct_set%d",set)); // replace name_of_hist1 with the actual name
    TH1D* hist2 = (TH1D*) f2.Get(Form("hadct_set%d",set)); // replace name_of_hist2 with the actual name
    
    // Remove existing fits
    hist1->GetListOfFunctions()->Clear();
    hist2->GetListOfFunctions()->Clear();
    
    // Name the resulting histogram comparison
    hist1->SetTitle("ADC time alignment comparison, all channels");

    // Fit the histograms dynamically
    double max1 = hist1->GetBinCenter(hist1->GetMaximumBin());
    double max2 = hist2->GetBinCenter(hist2->GetMaximumBin());
    
    hist1->Fit("gaus", "Q", "", max1 - adct_FWHM, max1 + adct_FWHM);
    hist2->Fit("gaus", "Q", "", max2 - adct_FWHM, max2 + adct_FWHM);
    
    // Extract fit results
    TF1* fit1 = hist1->GetFunction("gaus");
    TF1* fit2 = hist2->GetFunction("gaus");

    double mean1 = fit1->GetParameter(1);
    double sigma1 = fit1->GetParameter(2);

    double mean2 = fit2->GetParameter(1);
    double sigma2 = fit2->GetParameter(2);

    // Create a canvas and plot
    TCanvas* canvas = new TCanvas("canvas", "Comparison", 1600, 800);
    
    hist1->SetLineColor(kBlack);
    hist2->SetLineColor(kGreen);
    
    hist1->Draw();
    hist2->Draw("SAME");
    
    // Create the legend
    TLegend* leg = new TLegend(0.7, 0.8, 0.89, 0.89);
    leg->AddEntry(hist1, Form("Before Alignment: Mean = %.2f, #sigma = %.2f", mean1, sigma1), "l");
    leg->AddEntry(hist2, Form("After Alignment: Mean = %.2f, #sigma = %.2f", mean2, sigma2), "l");
    leg->Draw();

    canvas->SaveAs(adct_fitcomp_path);
    
    f1.Close();
    f2.Close();
}

//Gets the maximum bin on a user passed range
Int_t get_max_bin_on_range( TH1D *h, Double_t xlow, Double_t xhigh ){

  Int_t bin1 = h->FindBin(xlow);  // convert xlow to bin number
  Int_t bin2 = h->FindBin(xhigh); // convert xhigh to bin number

  // Initial maximum and its position
  Double_t maxVal = h->GetBinContent(bin1);
  Int_t maxBin = bin1;

  // Loop over bins in the range to find the maximum
  for (Int_t i = bin1; i <= bin2; i++) {
    Double_t binContent = h->GetBinContent(i);
    if (binContent > maxVal) {
      maxVal = binContent;
      maxBin = i;
    }
  }

  return maxBin;
}

void slice_histo( TH2D *h2, Int_t Nslices, Double_t fwhm, vector<Double_t> &cell, vector<Double_t> &mean, vector<Double_t> &err ){
  
  TH1D *cellslice[Nslices];

  for( Int_t i=0; i<Nslices; i++ ){
    
    Double_t cellval = (Double_t)i+0.5;
    Double_t meanval = 0.;
    Double_t errval = 0.;

    cellslice[i] = h2->ProjectionY(Form("cellslice_%d",i+1),i+1,i+1);
    Int_t sliceN = cellslice[i]->GetEntries();
    if( sliceN<fit_min )
      continue;

    //get reasonable limits on fit
    Double_t min = cellslice[i]->GetXaxis()->GetXmin();
    Double_t max = cellslice[i]->GetXaxis()->GetXmax();
    Int_t totalbins = cellslice[i]->GetNbinsX();
    Int_t binmax = cellslice[i]->GetMaximumBin();
    Double_t binmaxX = min+binmax*(max-min)/totalbins;
    Double_t llim = binmaxX-fwhm/2;
    Double_t ulim = binmaxX+fwhm/2;
    if( llim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
      binmax = get_max_bin_on_range( cellslice[i], E_hard_llim, max );
      binmaxX = min+binmax*(max-min)/totalbins;
      llim = E_hard_llim;
      ulim = binmaxX + ( binmaxX - llim );
    }
    
    TF1 *fit1;
    cellslice[i]->Fit("gaus","Q","",llim,ulim);
    fit1 = cellslice[i]->GetFunction("gaus");
    meanval = fit1->GetParameter(1);
    errval = fit1->GetParameter(2);

    cell.push_back(cellval);
    mean.push_back(meanval);
    err.push_back(errval);
  }

}

//Generates multigraph with timing alignment comparison
void generate_mg( std::string qr0path, std::string qr1path, std::string type, TMultiGraph *mg, Int_t set_N, Double_t avg ){

  //declare necessary variables
  Double_t x,y;
  avg = 0.;
  
  //open files and get tgraphs
  TFile *f = TFile::Open(qr0path.c_str());
  TFile *fqr = TFile::Open(qr1path.c_str()); 
  
  TGraphErrors *g = (TGraphErrors*)f->Get(Form("g%s_c_s%d",type.c_str(),set_N));
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){ //remove points where not enough data to fit (y=0) and get Y average
    x=-1;
    y=-1;
    g->GetPoint(i,x,y);
    if( y==0 ) g->RemovePoint( i );
    avg += y;

    if( x==hcal::maxHCalChan-1 ) break; //avoid adding points to mg
  }
  avg /= hcal::maxHCalChan;
  g->SetTitle("Unaligned");
  g->SetMarkerColor(kGray);
  g->SetMarkerStyle(33);
  g->SetMarkerSize(1);
  g->SetLineColor(kGray);
  g->SetLineWidth(1);
  mg->Add(g);

  TGraphErrors *gqr = (TGraphErrors*)fqr->Get(Form("g%s_c_s%d",type.c_str(),set_N));
  for( Int_t i=0; i<hcal::maxHCalChan; i++ ){ //remove points and offset values by 0.5 in X for clarity
    x=-1;
    y=-1;
    gqr->GetPoint(i,x,y);
    if( y!=0 ){
      gqr->SetPoint( i,x+0.5,y );
    }else{
      gqr->RemovePoint( i );
    }

    if( x==hcal::maxHCalChan-1 ) break; //avoid adding points to mg
  }
  //options for second tgrapherrors
  gqr->SetTitle("Aligned");
  gqr->SetMarkerColor(kBlack);
  gqr->SetMarkerStyle(34);
  gqr->SetMarkerSize(1);
  gqr->SetLineColor(kBlack);
  gqr->SetLineWidth(1);
  mg->Add(gqr);

}

//Main. Script to generate all necessary plots for calibration quality checks. See overleaf for quality check template
void quality_check( const char *experiment = "gmn", Int_t config = 9, Int_t pass = 1, bool h2only = true ){
  
  // set style
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);

  std::string h2opt = "";
  if( h2only )
    h2opt = "_lh2only";

  // Default plot output directory
  std::string plotdir = Form("../quality_plots/%s/conf%d%s/",experiment,config,h2opt.c_str());

  // Calibration analysis files output directory
  std::string outdir_path = gSystem->Getenv("OUT_DIR");

  // Get beam energy by configuration
  SBSconfig config_parameters(experiment,config);    
  Double_t beamE = config_parameters.GetEbeam();

  ////////////////////////
  //ADCt alignment plots//
  ////////////////////////

  //get paths to calibration output files and type
  std::string adct_type = "adct";
  std::string adct_param = "sbs.hcal.adc.timeoffset";
  std::string adct_align_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_%s_conf%d_qr0_pass%d.root",pass,adct_type.c_str(),experiment,config,pass);
  std::string adct_align_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_%s_conf%d_qr1_pass%d.root",pass,adct_type.c_str(),experiment,config,pass);
  std::string adct_param_path = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
  std::string adc_gain_path = Form("../energy/parameters/adcgaincoeff_%s%s_conf%d_pass%d.txt",experiment,h2opt.c_str(),config,pass);
  

  //determine how many calibration sets exist
  Int_t adct_align_Nsets = util::countSets( adct_param_path, adct_param );

  ///////////////////////////////////////
  // Get ADCt alignment comparison canvas
  TCanvas *c1 = new TCanvas("c1","adct comparison",1600/adct_align_Nsets,1200);
  
  // divide the canvas
  c1->Divide(1,adct_align_Nsets);

  //Generate adct multigraphs by calibration set
  for( Int_t set=0; set<adct_align_Nsets; set++ ){
    
    c1->cd(set+1);

    TMultiGraph *adctmg = new TMultiGraph();
    Double_t adct_avg;

    generate_mg( adct_align_qr0_path, adct_align_qr1_path, adct_type, adctmg, set, adct_avg );

    //options for mg
    adctmg->SetTitle(Form("HCal %s Timing Alignment Comparison, %s Config %d, Calibration Set %d",adct_type.c_str(),experiment,config,set));
    adctmg->GetXaxis()->SetTitle("HCal Channel");
    adctmg->GetYaxis()->SetTitle(Form("HCal %s (ns)",adct_type.c_str()));
    adctmg->Draw("AP");
    
    adctmg->GetYaxis()->SetRangeUser(adct_avg-Y_range_factor,adct_avg+Y_range_factor);
    
    c1->Modified();
    
    //c1->BuildLegend(); //must happen after mg is drawn

  }

  //ADC time aggregate plot comparisons

  std::string adct_fitcomp_path = plotdir + "adct_all_fitcomp.png";
  c1->SaveAs(adct_fitcomp_path.c_str());

  //Generate adct aggregate plots by calibration set
  for( Int_t set=0; set<adct_align_Nsets; set++ )
    fitAggHistograms(adct_align_qr0_path.c_str(),adct_align_qr1_path.c_str(),adct_fitcomp_path.c_str(),set);

  /////////////
  //TDC plots//
  /////////////

  //get paths to calibration output files and type
  std::string tdc_type = "tdc";
  std::string tdc_param = "sbs.hcal.tdc.offset";
  std::string tdc_align_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_%s_conf%d_qr0_pass%d.root",pass,tdc_type.c_str(),experiment,config,pass);
  std::string tdc_align_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/%salign_%s_conf%d_qr1_pass%d.root",pass,tdc_type.c_str(),experiment,config,pass);
  std::string tdc_param_path = Form("../timing/parameters/tdcoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);

  //determine how many calibration sets exist
  Int_t tdc_align_Nsets = util::countSets( tdc_param_path, tdc_param );

  ///////////////////////////
  // Get TDC alignment canvas
  TCanvas *c2 = new TCanvas("c2","tdc alignment",1600/tdc_align_Nsets,1200);
  
  // divide the canvas
  c2->Divide(1,tdc_align_Nsets); 

  //Generate tdc multigraphs by calibration set
  for( Int_t set=0; set<tdc_align_Nsets; set++ ){
    
    c2->cd(set+1);

    TMultiGraph *tdcmg = new TMultiGraph();
    Double_t tdc_avg;

    generate_mg( tdc_align_qr0_path, tdc_align_qr1_path, tdc_type, tdcmg, set, tdc_avg );

    //options for mg
    tdcmg->SetTitle(Form("HCal %s Timing Alignment Comparison, %s Config %d, Calibration Set %d",tdc_type.c_str(),experiment,config,set));
    tdcmg->GetXaxis()->SetTitle("HCal Channel");
    tdcmg->GetYaxis()->SetTitle(Form("HCal %s (ns)",tdc_type.c_str()));
    tdcmg->Draw("AP");
    
    tdcmg->GetYaxis()->SetRangeUser(tdc_avg-Y_range_factor,tdc_avg+Y_range_factor);
    
    c2->Modified();
    
    //c2->BuildLegend(); //must happen after mg is drawn

  }

  std::string tdc_fitcomp_path = plotdir + "tdc_fitcomp.png";
  c2->SaveAs(tdc_fitcomp_path.c_str());

  ////////////////
  //Energy Plots//
  ////////////////
  
  //get paths to calibration output files and type
  std::string adc_gain_param = "sbs.hcal.adc.gain";
  // std::string hcal_energy_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal_class_test_%s_conf%d_qr0_pass%d.root",pass,experiment,config,pass);
  // std::string hcal_energy_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal_class_test_%s_conf%d_qr1_pass%d.root",pass,experiment,config,pass);
  // std::string hcal_energy_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal%s_%s_conf%d_qr0_pass%d.root",pass,h2opt.c_str(),experiment,config,pass);
  // std::string hcal_energy_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal%s_%s_conf%d_qr1_pass%d.root",pass,h2opt.c_str(),experiment,config,pass);
  std::string hcal_energy_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal_lh2only_%s_conf%d_qr0_pass%d.root",pass,experiment,config,pass);
  std::string hcal_energy_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/energy/ecal_lh2only_%s_conf%d_qr1_pass%d.root",pass,experiment,config,pass);
  //std::string hcal_energy_qr1_path = outdir_path + Form("/hcal_calibrations/qreplay/qreplay%s_from_gmn_8_1_to_%s_%d_%d.root",h2opt.c_str(),experiment,config,pass);

  //determine how many calibration sets exist
  Int_t hcal_energy_Nsets = util::countSets( adc_gain_path, adc_gain_param );

  ///////////////////////////////
  // Get sampling fraction canvas
  TCanvas *c3 = new TCanvas("c3","sampling fraction",1600/hcal_energy_Nsets,1200);

  // divide the canvas
  c3->Divide(1,hcal_energy_Nsets);

  //Generate adct multigraphs by calibration set
  Double_t obs_fwhm = beamE*hcalE_sigma_factor;
  for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
    c3->cd(set+1);

    //Get pre-ecal samp frac
    TFile *f1 = TFile::Open(hcal_energy_qr0_path.c_str());
    TH1D *h1= (TH1D*)f1->Get(Form("hSF_set%d",set));
    Double_t x1min = h1->GetXaxis()->GetXmin();
    Double_t x1max = h1->GetXaxis()->GetXmax();
    Int_t x1totalbins = h1->GetNbinsX();
    Int_t x1maxX = (SF_hard_ulim-x1min)*x1totalbins/(x1max-x1min);
    h1->GetXaxis()->SetRange(0,x1maxX);
    h1->SetLineWidth(3);
    h1->SetLineColor(kRed);
    h1->Draw("hist");

    //Get reasonable limits on fit
    Int_t x1binmax = h1->GetMaximumBin();
    Double_t x1binmaxX = x1min+x1binmax*(x1max-x1min)/x1totalbins;
    Double_t x1llim = x1binmaxX-obs_fwhm/2;
    Double_t x1ulim = x1binmaxX+obs_fwhm/2;
    if( x1llim < SF_hard_llim ){ //if the lower limit is below minimum, recalculate
      x1binmax = get_max_bin_on_range( h1, SF_hard_llim, x1max );
      x1binmaxX = x1min+x1binmax*(x1max-x1min)/x1totalbins;
      x1llim = SF_hard_llim;
      x1ulim = x1binmaxX + ( x1binmaxX - x1llim );
    }

    //Get post-ecal samp frac
    TFile *f2 = TFile::Open(hcal_energy_qr1_path.c_str());
    TH1D *h2 = (TH1D*)f2->Get(Form("hSF_set%d",set));
    Double_t x2min = h2->GetXaxis()->GetXmin();
    Double_t x2max = h2->GetXaxis()->GetXmax();
    Int_t x2totalbins = h2->GetNbinsX();
    h2->SetLineWidth(2);
    h2->SetLineColor(kBlack);
    h2->Draw("hist same");

    //Get reasonable limits on fit
    Int_t x2binmax = h2->GetMaximumBin();
    Double_t x2binmaxX = x2min+x2binmax*(x2max-x2min)/x2totalbins;
    Double_t x2llim = x2binmaxX-obs_fwhm/2;
    Double_t x2ulim = x2binmaxX+obs_fwhm/2;
    if( x2llim < SF_hard_llim ){ //if the lower limit is below minimum, recalculate
      x2binmax = get_max_bin_on_range( h2, SF_hard_llim, x2max );
      x2binmaxX = x2min+x2binmax*(x2max-x2min)/x2totalbins;
      x2llim = SF_hard_llim;
      x2ulim = x2binmaxX + ( x2binmaxX - x2llim );
    }

    //Fit the pre-cal SF distribution with gaussian
    TF1 *fit1;
    h1->Fit("gaus","","",x1llim,x1ulim);
    fit1 = h1->GetFunction("gaus");
    Double_t fit1m = fit1->GetParameter(1);

    //Fit the pre-cal SF distribution with gaussian
    TF1 *fit2;
    h2->Fit("gaus","","",x2llim,x2ulim);
    fit2 = h2->GetFunction("gaus");
    Double_t fit2m = fit2->GetParameter(1);

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.03);
    legend->SetHeader(Form("HCal eCal Samp Frac SBS%d, Cal Set%d",config,set));
    legend->AddEntry(h1,Form("Pre-Cal, mean:%0.2f",fit1m),"l");
    legend->AddEntry(h2,Form("Post-Cal, mean:%0.2f",fit2m),"l");
    legend->Draw();

  }

  std::string sampFrac_datacomp_path = plotdir + "sampFrac_datacomp.png";
  c3->SaveAs(sampFrac_datacomp_path.c_str());

  ////////////////////////////
  // Get energy spectra canvas
  TCanvas *c4 = new TCanvas("c4","energy",1600/hcal_energy_Nsets,1200);

  // divide the canvas
  c4->Divide(1,hcal_energy_Nsets);

  //Generate adct multigraphs by calibration set
  for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
    c4->cd(set+1);

    //Get pre-ecal samp frac
    TFile *f1 = TFile::Open(hcal_energy_qr0_path.c_str());
    TH1D *h1= (TH1D*)f1->Get(Form("hE_set%d",set));    
    Double_t x1min = h1->GetXaxis()->GetXmin();
    Double_t x1max = h1->GetXaxis()->GetXmax();
    Int_t x1totalbins = h1->GetNbinsX();
    Int_t x1maxX = (E_hard_ulim-x1min)*x1totalbins/(x1max-x1min);
    h1->GetXaxis()->SetRange(0,x1maxX);
    h1->SetLineWidth(3);
    h1->SetLineColor(kRed);
    h1->Draw("hist");

    //Get reasonable limits on fit
    Int_t x1binmax = h1->GetMaximumBin();
    Double_t x1binmaxX = x1min+x1binmax*(x1max-x1min)/x1totalbins;
    Double_t x1llim = x1binmaxX-obs_fwhm/2;
    Double_t x1ulim = x1binmaxX+obs_fwhm/2;
    if( x1llim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
      x1binmax = get_max_bin_on_range( h1, E_hard_llim, x1max );
      x1binmaxX = x1min+x1binmax*(x1max-x1min)/x1totalbins;
      x1llim = E_hard_llim;
      x1ulim = x1binmaxX + ( x1binmaxX - x1llim );
    }

    //Get post-ecal energy
    TFile *f2 = TFile::Open(hcal_energy_qr1_path.c_str());
    TH1D *h2 = (TH1D*)f2->Get(Form("hE_set%d",set));
    Double_t x2min = h2->GetXaxis()->GetXmin();
    Double_t x2max = h2->GetXaxis()->GetXmax();
    Int_t x2totalbins = h2->GetNbinsX();
    h2->SetLineWidth(2);
    h2->SetLineColor(kBlack);
    h2->Draw("hist same");

    //Get reasonable limits on fit
    Int_t x2binmax = h2->GetMaximumBin();
    Double_t x2binmaxX = x2min+x2binmax*(x2max-x2min)/x2totalbins;
    Double_t x2llim = x2binmaxX-obs_fwhm/2;
    Double_t x2ulim = x2binmaxX+obs_fwhm/2;
    if( x2llim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
      x2binmax = get_max_bin_on_range( h2, E_hard_llim, x2max );
      x2binmaxX = x2min+x2binmax*(x2max-x2min)/x2totalbins;
      x2llim = E_hard_llim;
      x2ulim = x2binmaxX + ( x2binmaxX - x2llim );
    }

    //Fit the pre-cal E distribution with gaussian
    TF1 *fit1;
    h1->Fit("gaus","","",x1llim,x1ulim);
    fit1 = h1->GetFunction("gaus");
    Double_t fit1m = fit1->GetParameter(1);
    Double_t fit1s = fit1->GetParameter(2);
    Double_t fit1res = fit1s/fit1m*100; //Resolution is sig/peak in this context

    //Fit the pre-cal E distribution with gaussian
    TF1 *fit2;
    h2->Fit("gaus","","",x2llim,x2ulim);
    fit2 = h2->GetFunction("gaus");
    Double_t fit2m = fit2->GetParameter(1);
    Double_t fit2s = fit2->GetParameter(2);
    Double_t fit2res = fit2s/fit2m*100; //Resolution is FWHM/peak *100 (percent) where FWHM is 2.355*sigma

    //Add a legend
    auto legend = new TLegend(0.43,0.7,0.89,0.89);
    legend->SetTextSize(0.03);
    legend->SetHeader(Form("HCal eCal Energy Spectrum SBS%d, Cal Set%d",config,set));
    legend->AddEntry(h1,Form("Pre-Cal, mean:%0.2f, res:%0.2f",fit1m,fit1res),"l");
    legend->AddEntry(h2,Form("Post-Cal, mean:%0.2f, res:%0.2f",fit2m,fit2res),"l");
    legend->Draw();

  }

  std::string clusE_datacomp_path = plotdir + "clusE_datacomp.png";
  c4->SaveAs(clusE_datacomp_path.c_str());

  ///////////////////////////////////////////
  // Get monte carlo sampling fraction canvas
  TCanvas *c5 = new TCanvas("c5","MC sampling fraction",1600,1200);
  c5->cd();

  //get path to mc output files (no digitization, replay)
  std::string mc_path = outdir_path + Form("/hcal_calibrations/MC/hcalE_mc_class_%s_conf%d.root",experiment,config);

  //get path to mc output files (no digitization, replay)
  //std::string digmc_path = outdir_path + Form("/hcal_calibrations/MC/hcalE_mc_dig_%s_conf%d.root",experiment,config);

  //get path to mc output files (no digitization, replay)
  std::string digmc_path = outdir_path + Form("/hcal_calibrations/MC/hcalE_idx_mc_dig_%s_conf%d.root",experiment,config);

  //TFile *fSFmc = TFile::Open(mc_path.c_str());
  TFile *fSFmc = TFile::Open(digmc_path.c_str());
  //TH1D *hSFmc = (TH1D*)fSFmc->Get("hSampFrac_clus");
  TH1D *hSFmc = (TH1D*)fSFmc->Get("hSF_mc");
  //TH1D *hSFmc = (TH1D*)fSFmc->Get("hSF_mc_tridx");

  hSFmc->SetLineWidth(3);
  hSFmc->SetLineColor(kGreen);
  //hSFmc->SetFillStyle(3004);
  //hSFmc->SetFillColor(kGreen);
  hSFmc->Draw("hist");

  //get reasonable limits on fit
  Double_t SFmcmin = hSFmc->GetXaxis()->GetXmin();
  Double_t SFmcmax = hSFmc->GetXaxis()->GetXmax();
  Int_t SFmctotalbins = hSFmc->GetNbinsX();
  Int_t SFmcbinmax = hSFmc->GetMaximumBin();
  Double_t SFmcbinmaxX = SFmcmin+SFmcbinmax*(SFmcmax-SFmcmin)/SFmctotalbins;
  Double_t SFmcllim = SFmcbinmaxX-obs_fwhm/2;
  Double_t SFmculim = SFmcbinmaxX+obs_fwhm/2;
  if( SFmcllim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
    SFmcbinmax = get_max_bin_on_range( hSFmc, SF_hard_llim, SFmcmax );
    SFmcbinmaxX = SFmcmin+SFmcbinmax*(SFmcmax-SFmcmin)/SFmctotalbins;
    SFmcllim = SF_hard_llim;
    SFmculim = SFmcbinmaxX + ( SFmcbinmaxX - SFmcllim );
  }


  hSFmc->Fit("gaus","","",SFmcllim,SFmculim);
  TF1 *fitSFmc = hSFmc->GetFunction("gaus");
  Double_t fitSFmcp1 = fitSFmc->GetParameter(1)*100.;
  
  hSFmc->SetTitle(Form("MC Samp Frac SBS%d, peak %0.2f%%",config,fitSFmcp1));
  hSFmc->GetXaxis()->SetTitle("E_{clus}/KE");

  std::string MC_sampFrac_path = plotdir + "MC_sampFrac.png";
  c5->SaveAs(MC_sampFrac_path.c_str());

  fSFmc->Close();
  delete fSFmc;

  /////////////////////////////////////////
  // Get monte carlo energy spectrum canvas
  TCanvas *c6 = new TCanvas("c6","MC energy spectrum",1600,1200);
  c6->cd();

  //TFile *fEmc = TFile::Open(mc_path.c_str());
  TFile *fEmc = TFile::Open(digmc_path.c_str());
  //TH1D *hEmc = (TH1D*)fEmc->Get("hHCALe_clus");
  TH1D *hEmc = (TH1D*)fEmc->Get("hE_mc");
  //TH1D *hEmc = (TH1D*)fEmc->Get("hE_mc_tridx");

  hEmc->SetLineWidth(3);
  hEmc->SetLineColor(kGreen);
  //hEmc->SetFillStyle(3004);
  //hEmc->SetFillColor(kGreen);
  hEmc->Draw("hist");

  //get reasonable limits on fit
  Double_t Emcmin = hEmc->GetXaxis()->GetXmin();
  Double_t Emcmax = hEmc->GetXaxis()->GetXmax();
  Int_t Emctotalbins = hEmc->GetNbinsX();
  Int_t Emcbinmax = hEmc->GetMaximumBin();
  Double_t EmcbinmaxX = Emcmin+Emcbinmax*(Emcmax-Emcmin)/Emctotalbins;
  Double_t Emcllim = EmcbinmaxX-obs_fwhm/2;
  Double_t Emculim = EmcbinmaxX+obs_fwhm/2;
  if( Emcllim < E_hard_llim ){ //if the lower limit is below minimum, recalculate
    Emcbinmax = get_max_bin_on_range( hEmc, E_hard_llim, Emcmax );
    EmcbinmaxX = Emcmin+Emcbinmax*(Emcmax-Emcmin)/Emctotalbins;
    Emcllim = E_hard_llim;
    Emculim = EmcbinmaxX + ( EmcbinmaxX - Emcllim );
  }

  hEmc->Fit("gaus","0","",Emcllim,Emculim);
  TF1 *fitEmc = hEmc->GetFunction("gaus");
  Double_t fitEmcp1 = fitEmc->GetParameter(1);
  
  hEmc->SetTitle(Form("MC Primary Cluster Sampled Energy Spectrum SBS%d, peak %0.2f MeV",config,fitEmcp1*100));
  hEmc->GetXaxis()->SetTitle("GeV");

  std::string MC_energy_path = plotdir + "MC_energy.png";
  c6->SaveAs(MC_energy_path.c_str());

  fEmc->Close(); // Close the ROOT file
  delete fEmc; // Clean up the memory

  ///////////////////////////////////////////////////
  // Get dispersive direction uniformity check canvas
  TCanvas *c7 = new TCanvas("c7","Sampling Fraction by Row",1600/hcal_energy_Nsets,1200);
  
  // divide the canvas
  c7->Divide(1,hcal_energy_Nsets);
  c7->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);

  for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
    c7->cd(set+1);

    TFile *f1 = TFile::Open(hcal_energy_qr1_path.c_str());
    TH2D *hEvX = (TH2D*)f1->Get(Form("hSFvrow_set%d",set));
    vector<Double_t> Excell;
    vector<Double_t> Exmean;
    vector<Double_t> Exerr;
    
    slice_histo(hEvX, hcal::maxHCalRows, obs_fwhm, Excell, Exmean, Exerr );
    
    // Convert vectors to arrays
    Double_t* x = &Excell[0];  
    Double_t* y = &Exmean[0];
    Double_t* ey = &Exerr[0];

    hEvX->GetYaxis()->SetRangeUser(0.,0.2);
    hEvX->GetXaxis()->SetRangeUser(1.,23.);
    hEvX->GetYaxis()->SetTitle("E_{clus}/KE");
    hEvX->Draw("colz");
    
    //Make graphs with errors for reporting
    TGraphErrors *gEx = new TGraphErrors( hcal::maxHCalRows, x, y, 0, ey );
    
    gEx->SetTitle(Form("HCal SF vs X, SBS%d, set%d",config,set));
    gEx->SetMarkerStyle(33); // idx 20 Circles, idx 21 Boxes
    gEx->SetMarkerColor(kMagenta);
    gEx->SetMarkerSize(3);
    gEx->SetLineColor(kMagenta);
    gEx->SetLineWidth(3);
    for( Int_t r=0; r<hcal::maxHCalRows; r++ ){
      Double_t idx = r/hcal::maxHCalRows;
      if( y[r]==0 ) gEx->RemovePoint(idx);
    }
    gEx->Draw("same P");

    c7->Update();
  }

  std::string uni_X_path = plotdir + "uni_X.png";
  c7->SaveAs(uni_X_path.c_str());

  ///////////////////////////////////////////////////
  // Get transverse direction uniformity check canvas
  TCanvas *c8 = new TCanvas("c8","Sampling Fraction by Col",1600/hcal_energy_Nsets,1200);
  
  // divide the canvas
  c8->Divide(1,hcal_energy_Nsets);
  c8->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);

  for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
    c8->cd(set+1);

    TFile *f1 = TFile::Open(hcal_energy_qr1_path.c_str());
    TH2D *hEvY = (TH2D*)f1->Get(Form("hSFvcol_set%d",set));
    vector<Double_t> Eycell;
    vector<Double_t> Eymean;
    vector<Double_t> Eyerr;
    
    slice_histo(hEvY, hcal::maxHCalCols, obs_fwhm, Eycell, Eymean, Eyerr );
    
    // Convert vectors to arrays
    Double_t* x = &Eycell[0];  
    Double_t* y = &Eymean[0];
    Double_t* ey = &Eyerr[0];

    hEvY->GetYaxis()->SetRangeUser(0.,0.2);
    hEvY->GetXaxis()->SetRangeUser(1.,11.);
    hEvY->GetYaxis()->SetTitle("E_{clus}/KE");
    hEvY->Draw("colz");
    
    //Make graphs with errors for reporting
    TGraphErrors *gEy = new TGraphErrors( hcal::maxHCalCols, x, y, 0, ey );
    
    gEy->SetTitle(Form("HCal SF vs X, SBS%d, set%d",config,set));
    gEy->SetMarkerStyle(33); // idx 20 Circles, idx 21 Boxes
    gEy->SetMarkerColor(kMagenta);
    gEy->SetMarkerSize(3);
    gEy->SetLineColor(kMagenta);
    gEy->SetLineWidth(3);
    for( Int_t r=0; r<hcal::maxHCalCols; r++ ){
      Double_t idx = r/hcal::maxHCalCols;
      if( y[r]==0 ) gEy->RemovePoint(idx);
    }
    gEy->Draw("same P");

    c8->Update();
  }

  std::string uni_Y_path = plotdir + "uni_Y.png";
  c8->SaveAs(uni_Y_path.c_str());

  ///////////////////////////////////////////////
  // Get energy spectra MC/Data comparison canvas
  TCanvas *c9 = new TCanvas("c9","mc/data energy",1600/hcal_energy_Nsets,1200);

  // divide the canvas
  c9->Divide(1,hcal_energy_Nsets);

  //Generate adct multigraphs by calibration set
  for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
    c9->cd(set+1);
    
    //Get MC energy
    //TFile *fEmc = TFile::Open(mc_path.c_str());
    TFile *fEmc = TFile::Open(digmc_path.c_str());
    //TH1D *hEmc = (TH1D*)fEmc->Get("hHCALe_clus");
    TH1D *hEmc = (TH1D*)fEmc->Get("hE_mc");
    //TH1D *hEmc = (TH1D*)fEmc->Get("hE_mc_tridx");

    //get max value of histogram
    Double_t hEmc_max = hEmc->GetMaximum();

    Int_t hEmc_binMax = hEmc->GetMaximumBin();
    Double_t hEmc_binCenter = hEmc->GetBinCenter( hEmc_binMax );
    Double_t hEmc_fitLL = hEmc_binCenter - FWHM;
    Double_t hEmc_fitUL = hEmc_binCenter + FWHM;
    
    // Fit a Gaussian to this projection
    TF1 *hEmc_gfit = new TF1("hEmc_gfit", "gaus");
    hEmc->Fit( "hEmc_gfit", "QN", "", hEmc_fitLL, hEmc_fitUL ); // "Q" for quiet mode, "N" for no draw
    Double_t hEmc_fit_max = hEmc_gfit->GetMaximum(hEmc_fitLL,hEmc_fitUL);

    //Get data energy
    TFile *f2 = TFile::Open(hcal_energy_qr1_path.c_str());
    TH1D *hE = (TH1D*)f2->Get(Form("hE_set%d",set));

    //get max value of fit
    Int_t hE_binMax = hE->GetMaximumBin();
    Double_t hE_binCenter = hE->GetBinCenter( hE_binMax );
    Double_t hE_fitLL = hE_binCenter - FWHM;
    Double_t hE_fitUL = hE_binCenter + FWHM;
  
    // Fit a Gaussian to this projection
    TF1 *hE_gfit = new TF1( "hE_gfit", "gaus" );
    hE->Fit( "hE_gfit", "QN", "", hE_fitLL, hE_fitUL ); // "Q" for quiet mode, "N" for no draw
    Double_t hE_fit_max = hE_gfit->GetMaximum(hE_fitLL,hE_fitUL);

    //get scale factor, scale, and get max value
    Double_t scale_factor = hEmc_fit_max / hE_fit_max;
    hE->Scale(scale_factor);
    Double_t hE_max = hE->GetMaximum();

    Double_t yrangescale = hEmc_max*1.1;
    if( hEmc_max < hE_max )
      yrangescale = hE_max*1.1;

    //cout << "hE max: " << hE_max << " hEmc max: " << hEmc_max << " yrangescale " << yrangescale << endl;

    hEmc->SetLineWidth(3);
    hEmc->SetLineColor(kGreen);
    hEmc->SetFillStyle(3004);
    hEmc->SetFillColor(kGreen);
    hEmc->GetXaxis()->SetRangeUser(0,E_hard_ulim);
    hEmc->GetYaxis()->SetRangeUser(0,yrangescale);
    hEmc->GetXaxis()->SetTitle("GeV");
    hEmc->Draw("hist");

    hE->SetLineWidth(2);
    hE->SetLineColor(kBlack);
    hE->GetXaxis()->SetRangeUser(0,E_hard_ulim);
    hE->GetXaxis()->SetTitle("GeV");
    hE->Draw("hist same");

    //Add a legend
    auto legend = new TLegend(0.63,0.7,0.89,0.89);
    legend->SetTextSize(0.03);
    legend->SetHeader(Form("SBS%d, Cal Set%d",config,set));
    legend->AddEntry(hEmc,"MC","l");
    legend->AddEntry(hE,"Data (Scaled)","l");
    legend->Draw();

  }

  std::string energy_dmc_path = plotdir + "energy_dmc.png";
  c9->SaveAs(energy_dmc_path.c_str());

  //////////////////////////////////////////////////
  // Get sampling fraction MC/Data comparison canvas
  TCanvas *c10 = new TCanvas("c10","mc/data sampling fraction",1600/hcal_energy_Nsets,1200);

  // divide the canvas
  c10->Divide(1,hcal_energy_Nsets);

  //Generate adct multigraphs by calibration set
  for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
    c10->cd(set+1);
    
    //Get MC energy
    //TFile *fSFmc = TFile::Open(mc_path.c_str());
    TFile *fSFmc = TFile::Open(digmc_path.c_str());
    //TH1D *hSFmc = (TH1D*)fSFmc->Get("hSampFrac_clus");
    TH1D *hSFmc = (TH1D*)fSFmc->Get("hSF_mc");
    //TH1D *hSFmc = (TH1D*)fSFmc->Get("hSF_mc_tridx");

    //get max value of histogram
    Double_t hSFmc_max = hSFmc->GetMaximum();

    Int_t hSFmc_binMax = hSFmc->GetMaximumBin();
    Double_t hSFmc_binCenter = hSFmc->GetBinCenter( hSFmc_binMax );
    Double_t hSFmc_fitLL = hSFmc_binCenter - FWHM;
    Double_t hSFmc_fitUL = hSFmc_binCenter + FWHM;
  
    // Fit a Gaussian to this projection
    TF1 *hSFmc_gfit = new TF1("hSFmc_gfit", "gaus");
    hSFmc->Fit("hSFmc_gfit", "QN", "", hSFmc_fitLL, hSFmc_fitUL ); // "Q" for quiet mode, "N" for no draw
    Double_t hSFmc_fit_max = hSFmc_gfit->GetMaximum(hSFmc_fitLL,hSFmc_fitUL);

    //Get data energy
    TFile *f2 = TFile::Open(hcal_energy_qr1_path.c_str());
    TH1D *hSF = (TH1D*)f2->Get(Form("hSF_set%d",set));

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
    Double_t scale_factor = hSFmc_fit_max / hSF_fit_max;
    hSF->Scale(scale_factor);
    Double_t hSF_max = hSF->GetMaximum();

    Double_t yrangescale = hSFmc_max*1.1;
    if( hSFmc_max < hSF_max )
      yrangescale = hSF_max*1.1;

    hSFmc->SetLineWidth(3);
    hSFmc->SetLineColor(kGreen);
    hSFmc->SetFillStyle(3004);
    hSFmc->SetFillColor(kGreen);
    hSFmc->GetXaxis()->SetRangeUser(0,0.5);
    hSFmc->GetYaxis()->SetRangeUser(0,yrangescale);
    hSFmc->GetXaxis()->SetTitle("E_{clus}/KE");
    hSFmc->Draw("hist");

    hSF->SetLineWidth(2);
    hSF->SetLineColor(kBlack);
    hSF->GetXaxis()->SetRangeUser(0,0.5);
    hSF->GetXaxis()->SetTitle("E_{clus}/KE");
    hSF->Draw("hist same");

    //Add a legend
    auto legend = new TLegend(0.63,0.7,0.89,0.89);
    legend->SetTextSize(0.03);
    legend->SetHeader(Form("SBS%d, Cal Set%d",config,set));
    legend->AddEntry(hSFmc,"MC","l");
    legend->AddEntry(hSF,"Data (Scaled)","l");
    legend->Draw();

  }

  std::string sf_dmc_path = plotdir + "sf_dmc.png";
  c10->SaveAs(sf_dmc_path.c_str());

  /////////////////////////////////////////
  // Get TDC vs E before calibration canvas
  // TCanvas *c11 = new TCanvas("c11","TDC vs E Before Cal",1600/hcal_energy_Nsets,1200);

  // // divide the canvas
  // c11->Divide(1,hcal_energy_Nsets);

  // std::string hcal_tw_qr0_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/tdctw_class_%s_conf%d_qr0_pass%d.root",pass,experiment,config,pass);
  // std::string hcal_tw_qr1_path = outdir_path + Form("/hcal_calibrations/pass%d/timing/tdctw_class_%s_conf%d_qr1_pass%d.root",pass,experiment,config,pass);


  // //Generate histograms by calibration set
  // for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
  //   c11->cd(set+1);
  //   gStyle->SetOptStat(0);
  //   gStyle->SetOptFit(0);
  //   gStyle->SetPalette(53);

  //   TFile *f1 = TFile::Open(hcal_tw_qr0_path.c_str());
  //   TH2D *htvE = (TH2D*)f1->Get(Form("htdcvE_pblk_all_set%d",set));
  //   htvE->Draw("colz");

  // }

  // std::string tvE_before_path = plotdir + "tvE_before.png";
  // c11->SaveAs(tvE_before_path.c_str());

  // ////////////////////////////////////////
  // // Get TDC vs E after calibration canvas
  // TCanvas *c12 = new TCanvas("c12","TDC vs E After Cal",1600/hcal_energy_Nsets,1200);

  // // divide the canvas
  // c12->Divide(1,hcal_energy_Nsets);
  // gStyle->SetOptStat(0);
  // gStyle->SetOptFit(0);
  // gStyle->SetPalette(53);

  // //Generate histograms by calibration set
  // for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
  //   c12->cd(set+1);

  //   TFile *f1 = TFile::Open(hcal_tw_qr1_path.c_str());
  //   TH2D *htvE = (TH2D*)f1->Get(Form("htdcvE_pblk_all_set%d",set));
  //   htvE->Draw("colz");

  // }

  // std::string tvE_after_path = plotdir + "tvE_after.png";
  // c12->SaveAs(tvE_after_path.c_str());

  // ///////////////////////////////////////////////////////////////
  // // Get TDC vs TDC (timewalk corrected) after calibration canvas
  // TCanvas *c13 = new TCanvas("c13","TDC Correction Comparison",1600,1200/hcal_energy_Nsets);

  // // divide the canvas
  // c13->Divide(hcal_energy_Nsets,1);
  // gStyle->SetOptStat(0);
  // gStyle->SetOptFit(0);

  // Double_t approx_sigma = 5; //ns

  // //Generate histograms by calibration set
  // for( Int_t set=0; set<hcal_energy_Nsets; set++ ){
    
  //   c13->cd(set+1);

  //   TFile *f1 = TFile::Open(hcal_tw_qr1_path.c_str());
  //   TH1D *htdc = (TH1D*)f1->Get(Form("htdc_set%d",set));
  //   htdc->SetLineColor(kBlack);
  //   htdc->SetTitle(Form("HCal TDC Primary Block, Primary Cluster, Set %d",set));
  //   htdc->GetXaxis()->SetRangeUser(tdc_fit_llim,tdc_fit_ulim);
  //   htdc->Draw("hist");

  //   //get reasonable limits on fit
  //   Double_t tdcmin = htdc->GetXaxis()->GetXmin();
  //   Double_t tdcmax = htdc->GetXaxis()->GetXmax();
  //   Int_t tdctotalbins = htdc->GetNbinsX();
  //   Int_t tdcbinmax = htdc->GetMaximumBin();
  //   Double_t tdcbinmaxX = tdcmin+tdcbinmax*(tdcmax-tdcmin)/tdctotalbins;
  //   Double_t tdcllim = tdcbinmaxX-approx_sigma;
  //   Double_t tdculim = tdcbinmaxX+approx_sigma;
  //   htdc->Fit("gaus","0","",tdcllim,tdculim);
  //   TF1 *fittdc = htdc->GetFunction("gaus");
  //   Double_t fittdcp2 = fittdc->GetParameter(2);

  //   TH1D *htdc_tw = (TH1D*)f1->Get(Form("htdc_tw_set%d",set));
  //   htdc_tw->SetLineColor(kGreen);
  //   htdc_tw->SetTitle(Form("HCal TDC Primary Block, Primary Cluster, Set %d",set));
  //   htdc_tw->GetXaxis()->SetRangeUser(tdc_fit_llim,tdc_fit_ulim);
  //   htdc_tw->Draw("hist same");

  //   //get reasonable limits on fit
  //   Double_t tdc_twmin = htdc_tw->GetXaxis()->GetXmin();
  //   Double_t tdc_twmax = htdc_tw->GetXaxis()->GetXmax();
  //   Int_t tdc_twtotalbins = htdc_tw->GetNbinsX();
  //   Int_t tdc_twbinmax = htdc_tw->GetMaximumBin();
  //   Double_t tdc_twbinmaxX = tdc_twmin+tdc_twbinmax*(tdc_twmax-tdc_twmin)/tdc_twtotalbins;
  //   Double_t tdc_twllim = tdc_twbinmaxX-approx_sigma;
  //   Double_t tdc_twulim = tdc_twbinmaxX+approx_sigma;
  //   htdc_tw->Fit("gaus","0","",tdc_twllim,tdc_twulim);
  //   TF1 *fittdc_tw = htdc_tw->GetFunction("gaus");
  //   Double_t fittdc_twp2 = fittdc_tw->GetParameter(2);

  //   //Add a legend
  //   auto legend = new TLegend(0.43,0.7,0.89,0.89);
  //   legend->SetTextSize(0.03);
  //   legend->SetHeader(Form("HCal TDC Timewalk Correction SBS%d, Cal Set%d",config,set));
  //   legend->AddEntry(htdc,Form("Pre-Cal, sigma:%f",fittdcp2),"l");
  //   legend->AddEntry(htdc_tw,Form("Post-Cal, sigma:%f",fittdc_twp2),"l");
  //   legend->Draw();

  // }

  // std::string tdc_twcomp_path = plotdir + "tdc_twcomp.png";
  // c13->SaveAs(tdc_twcomp_path.c_str());

  /////////////////////////////////////////////////
  // Get NEV per channel per energy calibration set
  TCanvas *c14 = new TCanvas("c14","Cal Set Nev",1600/hcal_energy_Nsets,1200);

  std::vector<std::string> again_ts;
  ifstream again_const_file; again_const_file.open( adc_gain_path.c_str() );

  //Get all timestamps
  if(again_const_file.is_open()){
    std::string readline;

    while( getline( again_const_file, readline ) ){
      TString Tline = (TString)readline;
      
      if( Tline.BeginsWith("----") )
	again_ts.push_back(readline);
      
    }
  }

  // Bark if the number of timestamps doesn't match the number of calibration sets
  if( hcal_energy_Nsets != again_ts.size() )
    std::cout << "WARNING: Number of cal sets (" << hcal_energy_Nsets << ") not equal to number of ADC gain timestamps (" << again_ts.size() << ")" << std::endl;

  // divide the canvas
  c14->Divide(1,again_ts.size());
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TH2D *hNEV[hcal::gNstamp];

  //Generate histograms by calibration set
  for( Int_t set=0; set<again_ts.size(); set++ ){
    
    c14->cd(set+1);

    std::string nev_type = "#Number of events available for calibration";

    Double_t nev[hcal::maxHCalChan] = {0.};

    hNEV[set] = new TH2D(Form("NEV_set%d",set),Form("Number of events used for calibration per channel, set %d",set),12,0,12,24,0,24);

    util::readDB(adc_gain_path,again_ts[set],nev_type,nev);

    //Invert the histogram to match geometry facing away from target
    for( Int_t i=0; i<hcal::maxHCalChan; i++ ){
      Int_t row = (hcal::maxHCalRows) - i/hcal::maxHCalCols;
      Int_t col = (hcal::maxHCalCols) - i%hcal::maxHCalCols;

      // Int_t row = i/12;
      // Int_t col = i%12;


      // if( col==1 || col==hcal::maxHCalCols || row==1 || row==hcal::maxHCalRows )
      // 	hNEV[set]->SetBinContent(col,row,0);
      // else
      hNEV[set]->SetBinContent(col,row,nev[i]);

    }

    hNEV[set]->Draw("colz");

    // Create a TLatex object.
    TLatex l;
    l.SetTextSize(0.02);
    l.SetTextColor(kGreen);

    // Loop over the bins of the TH2D.
    for (int i = 1; i <= hNEV[set]->GetNbinsX(); ++i) {
        for (int j = 1; j <= hNEV[set]->GetNbinsY(); ++j) {
            // Get bin content.
            double binContent = hNEV[set]->GetBinContent(i, j);
            // Draw the bin content as text on the canvas.
            l.DrawLatex(hNEV[set]->GetXaxis()->GetBinCenter(i),
                        hNEV[set]->GetYaxis()->GetBinCenter(j),
                        Form("%d", int(binContent)));
        }
    }
    
    

  }

  std::string adcgain_nev_path = plotdir + "adcgain_nev.png";
  c14->SaveAs(adcgain_nev_path.c_str());
 
}
  
