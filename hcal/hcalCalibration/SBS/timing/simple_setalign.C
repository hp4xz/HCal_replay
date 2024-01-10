//SSeeds 12.12.23 - Rewrite of timing calibration script to run from diagnostic plots with tighter elastic cuts. Added functionality such that channels with insufficient statistics will be aligned to the row value.

#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

bool verb = true;

///////////////////////////////////////////////////////
//Manual overide/tune information
bool channel_override = false;
bool channel_override_tdc = false;
Int_t corrChan = 264; //Channel to correct
Int_t corrChan_tdc = -1; //Channel to correct
//Double_t corrChanOffset = -249.948; //sbs4 (-249.448, -249.948), sbs7 -239.234, sbs 11 (-271.841,-267.471,-271.841), sbs14 (-243.765,-243.765), sbs8 -253.499, sbs9 (-250.951)
Double_t corrChanOffset_tdc = 0.;
Double_t binXmax = 1e38;
Double_t binXmin = -1e38;
Double_t binXmax_tdc = 1e38;
Double_t binXmin_tdc = -1e38;
///////////////////////////////////////////////////////

//Main. If no arguments are passed to run_b and run_e, alignment over all runs will proceed.
void simple_setalign( const char *experiment = "gen", 
		      int config = 2, 
		      int pass = 1, 
		      int run_b = 0, 
		      int run_e = 0, 
		      int run_exclude_b = 0, 
		      int run_exclude_e = 0, 
		      bool best_clus = false, 
		      double adct_target = 0.,
		      double tdc_target = 0.,
		      double corrChanOffset = 0){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  std::string date = util::getDate();

  // Standardize pass 0/1, differences handled during histogram creation
  if( pass==0 )
    pass=1;

  // set style
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);

  std::string out_dir = gSystem->Getenv("OUT_DIR");
  std::string in_path = out_dir + Form("/hcal_calibrations/pass%d/diagnostics/timingdiag_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.root",pass,experiment,config,pass,run_b,run_e,run_exclude_b,run_exclude_e);

  std::ifstream file(in_path.c_str());
  if (!file.good()){
    cout << "ERROR: No file located at " << in_path << endl;
    return;
  }

  //Path to analysis file
  std::string out_path = out_dir + Form("/hcal_calibrations/pass%d/timing/setalign_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.root",pass,experiment,config,pass,run_b,run_e,run_exclude_b,run_exclude_e);

  //New Offset text paths
  std::string new_adctoffset_path = Form("parameters/simple_adctoffsets_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.txt",experiment,config,pass,run_b,run_e,run_exclude_b,run_exclude_e);
  std::string new_tdcoffset_path = Form("parameters/simple_tdcoffsets_%s_conf%d_pass%d_%d_to_%d_exclude_%d_to_%d.txt",experiment,config,pass,run_b,run_e,run_exclude_b,run_exclude_e);

  //Paths to database. These should point to database files which produced the data analyzed.
  std::string db_path;
  std::string exper = experiment;
  if( exper.compare("gmn")==0 ){
    if( pass==0 || pass==1 )
      db_path = hcal::gmn_p1_db_path;
    else if( pass==2 )
      db_path = hcal::gmn_p2_db_path;
    else{
      cout << "ERROR: Enter a valid pass for gmn." << endl;
      return;
    }
  }else if( exper.compare("gen")==0 ){
    db_path = hcal::gmn_p2_db_path;
  }else{
    cout << "ERROR: Enter a valid experiment (gmn or gen)" << endl;
    return;
  }

  std::string old_db_path = db_path + "/db_sbs.hcal.dat";
  std::string db_adctoffset_variable = "sbs.hcal.adc.timeoffset";
  std::string db_tdcoffset_variable = "sbs.hcal.tdc.offset";
  std::string db_tdccalib_variable = "sbs.hcal.tdc.calib";
  std::string db_tmax_variable = "sbs.hcal.tmax";

  // create output analysis file
  TFile *fout = new TFile( out_path.c_str(), "RECREATE" );

  // path to plots
  std::string plotdir = Form("../quality_plots/%s/conf%d/",experiment,config);

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s_p%d/",experiment,pass);

  int nruns = -1; //Analyze all available runs for this configuration

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
  std::string offset_timestamp = runs[0].adct_ts; //Assumes that each input file has only one timestamp in database
  vector<caldiag> cut;
  std::string targ = runs[0].target;
  Int_t mag = runs[0].sbsmag / 21; //convert to percent where max field is at 2100A
  util::ReadDiagnosticCutList(struct_dir,experiment,config,targ,mag,verb,cut); //needed only for min ev
  //Int_t min_epc = cut[0].min_ev;
  Int_t min_epc = 300; //TEST ONLY
  Double_t atime_fwhm = cut[0].atime_sig;

  //get vector of runnumbers
  vector<Double_t> runnumbers;
  for (Int_t irun=0; irun<nruns; ++irun){
    runnumbers.push_back(runs[irun].runnum);
  }

  // open file, get histograms
  TFile *f1 = TFile::Open(in_path.c_str());

  TH2D *hatimerun_ecut = (TH2D*)f1->Get("hatimerun_ecut");
  TH2D *hatimerow_ecut = (TH2D*)f1->Get("hatimerow_ecut");
  TH2D *hatimecol_ecut = (TH2D*)f1->Get("hatimecol_ecut");
  TH2D *hatimeid_ecut = (TH2D*)f1->Get("hatimeid_ecut");
  TH2D *hatimex_ecut = (TH2D*)f1->Get("hatimex_ecut");
  TH2D *hatimey_ecut = (TH2D*)f1->Get("hatimey_ecut");

  TH2D *htimerun_ecut = (TH2D*)f1->Get("htimerun_ecut");
  TH2D *htimerow_ecut = (TH2D*)f1->Get("htimerow_ecut");
  TH2D *htimecol_ecut = (TH2D*)f1->Get("htimecol_ecut");
  TH2D *htimeid_ecut = (TH2D*)f1->Get("htimeid_ecut");
  TH2D *htimex_ecut = (TH2D*)f1->Get("htimex_ecut");
  TH2D *htimey_ecut = (TH2D*)f1->Get("htimey_ecut");

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

  //Set up structures to hold calibration parameters
  calset adct_cal;
  calset tdc_cal;
  adct_cal.timestamp = offset_timestamp.c_str();
  tdc_cal.timestamp = offset_timestamp.c_str();
  util::readDB( old_db_path, offset_timestamp, db_adctoffset_variable, adct_cal.old_param );
  util::readDB( old_db_path, offset_timestamp, db_tdcoffset_variable, tdc_cal.old_param );
  util::readDB( old_db_path, offset_timestamp, db_tdccalib_variable, tdc_cal.tdc_calib );

  //Make histogram to store all adct slices
  TH1D *hadctAll = hatimeid->ProjectionY("hadctAll", 1, hatimeid->GetNbinsX(), "e");
  //Get fit parameters
  vector<Double_t> AllP = util::fitGaussianAndGetFineParams(hadctAll, atime_fwhm);

  //Make histogram to store all tdc slices
  TH1D *htdcAll = htimeid->ProjectionY("htdcAll", 1, htimeid->GetNbinsX(), "e");
  //Get fit parameters
  vector<Double_t> AllPtdc = util::fitGaussianAndGetFineParams(htdcAll, atime_fwhm);

  //Make canvas to write adct projection
  TCanvas *cAll = new TCanvas("cAll","ATime Diff All Channels",1600,1200);
  cAll->SetGrid();

  cAll->cd();
  hadctAll->Draw();

  cAll->Update();

  //Make canvas to write tdc projection
  TCanvas *cAlltdc = new TCanvas("cAlltdc","Time Diff All Channels",1600,1200);
  cAlltdc->SetGrid();

  cAlltdc->cd();
  htdcAll->Draw();

  cAlltdc->Update();

  /////////////////////
  //Course alignments

  // Output text file for new ADCt offsets
  ofstream writeADCtParFile;
  writeADCtParFile.open( new_adctoffset_path );

  ////////////////////
  // ADCtime vs ID

  /////////////////////////////////////////////////////////////////////////////////
  // Fits to electron arm elastic cuts only (very wide to get baseline statistics)

  TCanvas *c0 = new TCanvas("c0","ATime Diff vs Id (wide e-arm cuts only)",1600,1200);
  c0->Divide(1,2);
  c0->SetGrid();
    
  vector<Double_t> c0cell;
  vector<Double_t> c0mean;
  vector<Double_t> c0err;
    
  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *hadctWide = (TH2D*)hatimeid_ecut->Clone("hadctWide");

  Int_t c0binsX = hadctWide->GetNbinsX();

  TCanvas *c0v1 = new TCanvas("c0v1","ATime Diff fits HCal Top Half",1600,1200);
  TCanvas *c0v2 = new TCanvas("c0v2","ATime Diff fits HCal Bottom Half",1600,1200);

  if(verb){
    util::sliceHCalIDHisto(hadctWide, c0binsX, atime_fwhm, min_epc, c0v1, c0v2, c0cell, c0mean, c0err, binXmin, binXmax );
    c0v1->Write();
    c0v2->Write();
  } else
    util::sliceHisto(hadctWide, c0binsX, atime_fwhm, min_epc, c0cell, c0mean, c0err, binXmin, binXmax );

  c0->cd(1);

  // Convert vectors to arrays
  Double_t* c0x = &c0cell[0];  
  Double_t* c0y = &c0mean[0];
  Double_t* c0ey = &c0err[0];
  
  //hadctWide->GetYaxis()->SetRangeUser(-10.,10.);
  hadctWide->GetYaxis()->SetTitle("ns");
  hadctWide->GetXaxis()->SetTitle("Channel");
  hadctWide->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c0g = new TGraphErrors( c0cell.size(), c0x, c0y, 0, c0ey );
    
  c0g->SetTitle(Form("HCal adctime vs id (earm cuts only), SBS%d",config));
  c0g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c0g->SetMarkerColor(kBlack);
  c0g->SetMarkerSize(0.6);
  c0g->SetLineColor(kBlack);
  c0g->SetLineWidth(1);
  c0g->Draw("same P");

  c0->Update();
  c0->cd(2);

  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *hadctWideRow = (TH2D*)hatimerow_ecut->Clone("hadctWideRow");

  vector<Double_t> c0cell_row;
  vector<Double_t> c0mean_row;
  vector<Double_t> c0err_row;
    
  Int_t c0binsX_row = hadctWideRow->GetNbinsX();

  util::sliceHisto(hadctWideRow, c0binsX_row, atime_fwhm, min_epc, c0cell_row, c0mean_row, c0err_row, binXmin, binXmax );
    
  // Convert vectors to arrays
  Double_t* c0x_row = &c0cell_row[0];  
  Double_t* c0y_row = &c0mean_row[0];
  Double_t* c0ey_row = &c0err_row[0];
  
  //hadctWideRow->GetYaxis()->SetRangeUser(-10.,10.);
  hadctWideRow->GetYaxis()->SetTitle("ns");
  hadctWideRow->GetXaxis()->SetTitle("");
  hadctWideRow->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c0g_row = new TGraphErrors( c0cell_row.size(), c0x_row, c0y_row, 0, c0ey_row );
    
  c0g_row->SetTitle(Form("HCal adctime vs row, SBS%d",config));
  c0g_row->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c0g_row->SetMarkerColor(kBlack);
  c0g_row->SetMarkerSize(0.6);
  c0g_row->SetLineColor(kBlack);
  c0g_row->SetLineWidth(1);
  c0g_row->Draw("same P");

  c0->Update();
  c0->Write();
 
  /////////////////////////////////
  //Fits to tight elastic cut data

  TCanvas *c1 = new TCanvas("c1","ATime Diff vs ID Tight Cut",1600,1200);
  c1->Divide(1,2);
  c1->SetGrid();
    
  vector<Double_t> c1cell;
  vector<Double_t> c1mean;
  vector<Double_t> c1err;
    
  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *hadctAlign;
  if(best_clus)
    hadctAlign = (TH2D*)hatimeid_bc->Clone("hadctAlign");
  else
    hadctAlign = (TH2D*)hatimeid->Clone("hadctAlign");

  Int_t c1binsX = hadctAlign->GetNbinsX();

  TCanvas *c1v1 = new TCanvas("c0v1","ATime Diff fits All Cuts HCal Top Half",1600,1200);
  TCanvas *c1v2 = new TCanvas("c0v2","ATime Diff fits All Cuts HCal Bottom Half",1600,1200);
   
  if(verb){
    util::sliceHCalIDHisto(hadctAlign, c1binsX, atime_fwhm, min_epc, c1v1, c1v2, c1cell, c1mean, c1err, binXmin, binXmax );
    c1v1->Write();
    c1v2->Write();
  }else
    util::sliceHisto(hadctAlign, c1binsX, atime_fwhm, min_epc, c1cell, c1mean, c1err, binXmin, binXmax );
 
  c1->cd(1);

  // Convert vectors to arrays
  Double_t* c1x = &c1cell[0];  
  Double_t* c1y = &c1mean[0];
  Double_t* c1ey = &c1err[0];
  
  //hadctAlign->GetYaxis()->SetRangeUser(-10.,10.);
  hadctAlign->GetYaxis()->SetTitle("ns");
  hadctAlign->GetXaxis()->SetTitle("");
  hadctAlign->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c1g = new TGraphErrors( c1cell.size(), c1x, c1y, 0, c1ey );
    
  c1g->SetTitle(Form("HCal adctime vs id, SBS%d",config));
  c1g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c1g->SetMarkerColor(kBlack);
  c1g->SetMarkerSize(0.6);
  c1g->SetLineColor(kBlack);
  c1g->SetLineWidth(1);
  c1g->Draw("same P");

  c1->Update();
  c1->cd(2);

  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *hadctAlignRow;
  if(best_clus)
    hadctAlignRow = (TH2D*)hatimerow_bc->Clone("hadctAlignRow");
  else
    hadctAlignRow = (TH2D*)hatimerow->Clone("hadctAlignRow");

  vector<Double_t> c1cell_row;
  vector<Double_t> c1mean_row;
  vector<Double_t> c1err_row;
    
  Int_t c1binsX_row = hadctAlignRow->GetNbinsX();

  util::sliceHisto(hadctAlignRow, c1binsX_row, atime_fwhm, min_epc, c1cell_row, c1mean_row, c1err_row, binXmin, binXmax );
    
  // Convert vectors to arrays
  Double_t* c1x_row = &c1cell_row[0];  
  Double_t* c1y_row = &c1mean_row[0];
  Double_t* c1ey_row = &c1err_row[0];
  
  //hadctAlignRow->GetYaxis()->SetRangeUser(-10.,10.);
  hadctAlignRow->GetYaxis()->SetTitle("ns");
  hadctAlignRow->GetXaxis()->SetTitle("");
  hadctAlignRow->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c1g_row = new TGraphErrors( c1cell_row.size(), c1x_row, c1y_row, 0, c1ey_row );
    
  c1g_row->SetTitle(Form("HCal adctime vs row, SBS%d",config));
  c1g_row->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c1g_row->SetMarkerColor(kBlack);
  c1g_row->SetMarkerSize(0.6);
  c1g_row->SetLineColor(kBlack);
  c1g_row->SetLineWidth(1);
  c1g_row->Draw("same P");

  c1->Update();
  c1->Write();

  //Calculate new offsets

  //check to be sure the size of the tgraph x axis matches hcal channels
  if( hadctAlign->GetXaxis()->GetNbins() != hcal::maxHCalChan ){
    cout << "ERROR: adct histogram x axis size not equal to HCal channels." << endl;
    return;
  }

  //get an all fitted-channel average for last-resort on low statistics
  double cell_avg = AllP[1];

  if(verb)
    cout << "All ADCt Cell Mean Value: " << cell_avg << endl;

  int cell = 0;
  double mean_vals[hcal::maxHCalChan] = {0.};

  //Write to output text file
  writeADCtParFile << "#HCal ADCt offsets obtained " << date.c_str() << " for " << experiment << " config " << config << endl;
  writeADCtParFile << "#Offsets obtained from fits over ADCt distributions." << endl << endl;
  writeADCtParFile << db_adctoffset_variable << " =" << endl;

  //After suitable check on dimensions, loop over hcal rows and cols
  for( int r = 0; r<hcal::maxHCalRows; ++r){
    for( int c = 0; c<hcal::maxHCalCols; ++c){
      
      //mean_wide_val will serve as a default value if fine and course tight cuts don't do the job
      //if the cell cannot be found in the vector of x-axis tgraph points, no fit was performed
      //if no fit was performed, insufficient statistics exist for that cell

      //get the mean fitted/wide-cut value for this cell
      double mean_wide_val = -1000.;
      auto it_w = std::find(c0cell.begin(), c0cell.end(), (double)cell+0.5);
      if (it_w != c0cell.end()){
        int idx = std::distance(c0cell.begin(), it_w);
	mean_wide_val = c0mean[idx];
      }else{
	auto it_w_r = std::find(c0cell_row.begin(), c0cell_row.end(), (double)r+0.5);
	if (it_w_r != c0cell_row.end()){
	  int idx_r = std::distance(c0cell_row.begin(), it_w_r);
	  mean_wide_val = c0mean_row[idx_r];
	}
      }

      if( mean_wide_val == -1000. )
	mean_wide_val = cell_avg;

      //get the mean fitted/tight-cut value for this cell
      double mean_val = -1000.;
      auto it = std::find(c1cell.begin(), c1cell.end(), (double)cell+0.5);
      if (it != c1cell.end()){
        int idx = std::distance(c1cell.begin(), it);
	mean_val = c1mean[idx];
	if(verb)
	  cout << "Search Cell : Found Cell - " << (double)cell+0.5 << " : " << c1cell[idx] << endl;
      }else{
	auto it_r = std::find(c1cell_row.begin(), c1cell_row.end(), (double)r+0.5);
	if (it_r != c1cell_row.end()){
	  cout << "Insufficient statistics on cell " << cell << ", using fit to row." << endl;
	  int idx_r = std::distance(c1cell_row.begin(), it_r);
	  mean_val = c1mean_row[idx_r];
	  if(verb)
	    cout << "Search Row : Found Row - " << (double)r+0.5 << " : " << c1cell_row[idx_r] << endl;
	}
      }

      //default to wide value if no value can be obtained from id or row fits
      if( mean_val == -1000. ){
	cout << "WARNING: Cell " << cell << ", defaulting to fit over all data." << endl;
	mean_val = mean_wide_val;
      }

      //populate the adct offset array
      mean_vals[cell] = mean_val;

      if( channel_override && cell==corrChan )  //manual override option
	writeADCtParFile << corrChanOffset << " ";
      else
	writeADCtParFile << mean_vals[cell] + adct_cal.old_param[cell] - adct_target << " ";

      if(verb)
	cout << "Fitted ADCt value for cell " << cell << ": " << mean_vals[cell] << endl;

      //advance the id for the next loop
      cell++;

    } // endloop over columns 
    writeADCtParFile << endl;

  } // endloop over rows
  writeADCtParFile << endl << endl << endl;

  //Now write the new values in the format expected in the database for adc time
  cell = 0;
  writeADCtParFile << "(inverted) " << db_adctoffset_variable << " =" << endl;
  for( int r = 0; r<hcal::maxHCalRows; ++r){
    for( int c = 0; c<hcal::maxHCalCols; ++c){
      if( channel_override && cell==corrChan )  //manual override option
	writeADCtParFile << -corrChanOffset << " ";
      else
	writeADCtParFile << -(mean_vals[cell] - adct_cal.old_param[cell] - adct_target) << " ";
      cell++;
    } // endloop over columns 
    writeADCtParFile << endl;
  } // endloop over rows
  writeADCtParFile << endl << endl << endl;

  writeADCtParFile.close();

  ////////////////////
  // TDC time vs ID

  // Get tdc calibration constant for later use
  Double_t calib_const = tdc_cal.tdc_calib;

  // Output text file for new ADCt offsets
  ofstream writeTDCParFile;
  writeTDCParFile.open( new_tdcoffset_path );

   /////////////////////////////////////////////////////////////////////////////////
  // Fits to electron arm elastic cuts only (very wide to get baseline statistics)

  TCanvas *c6 = new TCanvas("c6","Time Diff vs Id (wide e-arm cuts only)",1600,1200);
  c6->Divide(1,2);
  c6->SetGrid();
    
  c6->cd(1);

  vector<Double_t> c6cell;
  vector<Double_t> c6mean;
  vector<Double_t> c6err;
    
  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *htdcWide = (TH2D*)htimeid_ecut->Clone("htdcWide");

  Int_t c6binsX = htdcWide->GetNbinsX();

  util::sliceHisto(htdcWide, c6binsX, atime_fwhm, min_epc, c6cell, c6mean, c6err, binXmin_tdc, binXmax_tdc );
    
  // Convert vectors to arrays
  Double_t* c6x = &c6cell[0];  
  Double_t* c6y = &c6mean[0];
  Double_t* c6ey = &c6err[0];
  
  //htdcWide->GetYaxis()->SetRangeUser(-10.,10.);
  htdcWide->GetYaxis()->SetTitle("ns");
  htdcWide->GetXaxis()->SetTitle("Channel");
  htdcWide->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c6g = new TGraphErrors( c6cell.size(), c6x, c6y, 0, c6ey );
    
  c6g->SetTitle(Form("HCal tdcime vs id (earm cuts only), SBS%d",config));
  c6g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c6g->SetMarkerColor(kBlack);
  c6g->SetMarkerSize(0.6);
  c6g->SetLineColor(kBlack);
  c6g->SetLineWidth(1);
  c6g->Draw("same P");

  c6->Update();
  c6->cd(2);

  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *htdcWideRow = (TH2D*)hatimerow_ecut->Clone("htdcWideRow");

  vector<Double_t> c6cell_row;
  vector<Double_t> c6mean_row;
  vector<Double_t> c6err_row;
    
  Int_t c6binsX_row = htdcWideRow->GetNbinsX();

  util::sliceHisto(htdcWideRow, c6binsX_row, atime_fwhm, min_epc, c6cell_row, c6mean_row, c6err_row, binXmin_tdc, binXmax_tdc );
    
  // Convert vectors to arrays
  Double_t* c6x_row = &c6cell_row[0];  
  Double_t* c6y_row = &c6mean_row[0];
  Double_t* c6ey_row = &c6err_row[0];
  
  //htdcWideRow->GetYaxis()->SetRangeUser(-10.,10.);
  htdcWideRow->GetYaxis()->SetTitle("ns");
  htdcWideRow->GetXaxis()->SetTitle("");
  htdcWideRow->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c6g_row = new TGraphErrors( c6cell_row.size(), c6x_row, c6y_row, 0, c6ey_row );
    
  c6g_row->SetTitle(Form("HCal tdcime vs row, SBS%d",config));
  c6g_row->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c6g_row->SetMarkerColor(kBlack);
  c6g_row->SetMarkerSize(0.6);
  c6g_row->SetLineColor(kBlack);
  c6g_row->SetLineWidth(1);
  c6g_row->Draw("same P");

  c6->Update();
  c6->Write();

  /////////////////////////////////////
  // TDC fits to both arm elastic cuts

  TCanvas *c2 = new TCanvas("c2","Time Diff vs Id",1600,1200);
  c2->Divide(1,2);
  c2->SetGrid();
    
  c2->cd(1);

  vector<Double_t> c2cell;
  vector<Double_t> c2mean;
  vector<Double_t> c2err;
    
  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *htdcAlign;
  if(best_clus)
    htdcAlign = (TH2D*)htimeid_bc->Clone("htdcAlign");
  else
    htdcAlign = (TH2D*)htimeid->Clone("htdcAlign");

  Int_t c2binsX = htdcAlign->GetNbinsX();

  util::sliceHisto(htdcAlign, c2binsX, atime_fwhm, min_epc, c2cell, c2mean, c2err );
    
  // Convert vectors to arrays
  Double_t* c2x = &c2cell[0];  
  Double_t* c2y = &c2mean[0];
  Double_t* c2ey = &c2err[0];
  
  //htdcAlign->GetYaxis()->SetRangeUser(-10.,10.);
  htdcAlign->GetYaxis()->SetTitle("ns");
  htdcAlign->GetXaxis()->SetTitle("");
  htdcAlign->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c2g = new TGraphErrors( c2cell.size(), c2x, c2y, 0, c2ey );
    
  c2g->SetTitle(Form("HCal tdc time vs id, SBS%d",config));
  c2g->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c2g->SetMarkerColor(kBlack);
  c2g->SetMarkerSize(0.6);
  c2g->SetLineColor(kBlack);
  c2g->SetLineWidth(1);
  c2g->Draw("same P");

  c2->Update();
  
  c2->cd(2);

  //Clone the histogram which corresponds to the cluster method chosen
  TH2D *htdcAlignRow;
  if(best_clus)
    htdcAlignRow = (TH2D*)htimerow_bc->Clone("htdcAlignRow");
  else
    htdcAlignRow = (TH2D*)htimerow->Clone("htdcAlignRow");

  vector<Double_t> c2cell_row;
  vector<Double_t> c2mean_row;
  vector<Double_t> c2err_row;
    
  Int_t c2binsX_row = htdcAlignRow->GetNbinsX();

  util::sliceHisto(htdcAlignRow, c2binsX_row, atime_fwhm, min_epc, c2cell_row, c2mean_row, c2err_row );
    
  // Convert vectors to arrays
  Double_t* c2x_row = &c2cell_row[0];  
  Double_t* c2y_row = &c2mean_row[0];
  Double_t* c2ey_row = &c2err_row[0];
  
  //htdcAlignRow->GetYaxis()->SetRangeUser(-10.,10.);
  htdcAlignRow->GetYaxis()->SetTitle("ns");
  htdcAlignRow->GetXaxis()->SetTitle("");
  htdcAlignRow->Draw("colz");
    
  //Make graphs with errors for reporting
  TGraphErrors *c2g_row = new TGraphErrors( c2cell_row.size(), c2x_row, c2y_row, 0, c2ey_row );
    
  c2g_row->SetTitle(Form("HCal adctime vs row, SBS%d",config));
  c2g_row->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
  c2g_row->SetMarkerColor(kBlack);
  c2g_row->SetMarkerSize(0.6);
  c2g_row->SetLineColor(kBlack);
  c2g_row->SetLineWidth(1);
  c2g_row->Draw("same P");

  c2->Update();
  c2->Write();  

  //Calculate new tdc offsets

  //check to be sure the size of the tgraph x axis matches hcal channels
  if( htdcAlign->GetXaxis()->GetNbins() != hcal::maxHCalChan ){
    cout << "ERROR: tdc histogram x axis size not equal to HCal channels." << endl;
    return;
  }

  //get an all fitted-channel average for last-resort on low statistics
  double cell_avg_tdc = AllPtdc[1];

  if(verb)
    cout << "All TDC Cell Mean Value: " << cell_avg_tdc << endl;

  int cell_tdc = 0;
  double mean_vals_tdc[hcal::maxHCalChan] = {0.};

  //Write to output text file
  writeTDCParFile << "#HCal TDC offsets obtained " << date.c_str() << " for " << experiment << " config " << config << endl;
  writeTDCParFile << "#Offsets obtained from fits over TDC distributions." << endl << endl;
  writeTDCParFile << db_tdcoffset_variable << " =" << endl;

  //After suitable check on dimensions, loop over hcal rows and cols
  for( int r = 0; r<hcal::maxHCalRows; ++r){
    for( int c = 0; c<hcal::maxHCalCols; ++c){
      
      //mean_wide_val will serve as a default value if fine and course tight cuts don't do the job
      //if the cell cannot be found in the vector of x-axis tgraph points, no fit was performed
      //if no fit was performed, insufficient statistics exist for that cell

      //get the mean fitted/wide-cut value for this cell
      double mean_wide_val = -1000.;
      auto it_w = std::find(c6cell.begin(), c6cell.end(), (double)cell_tdc+0.5);
      if (it_w != c6cell.end()){
        int idx = std::distance(c6cell.begin(), it_w);
	mean_wide_val = c6mean[idx];
      }else{
	auto it_w_r = std::find(c6cell_row.begin(), c6cell_row.end(), (double)r+0.5);
	if (it_w_r != c6cell_row.end()){
	  int idx_r = std::distance(c6cell_row.begin(), it_w_r);
	  mean_wide_val = c6mean_row[idx_r];
	}
      }

      if( mean_wide_val == -1000. )
	mean_wide_val = cell_avg_tdc;

      //get the mean fitted/tight-cut value for this cell
      double mean_val = -1000.;
      auto it = std::find(c2cell.begin(), c2cell.end(), (double)cell_tdc+0.5);
      if (it != c2cell.end()){
        int idx = std::distance(c2cell.begin(), it);
	mean_val = c2mean[idx];
	if(verb)
	  cout << "Search Cell : Found Cell - " << (double)cell_tdc+0.5 << " : " << c2cell[idx] << endl;
      }else{
	auto it_r = std::find(c2cell_row.begin(), c2cell_row.end(), (double)r+0.5);
	if (it_r != c2cell_row.end()){
	  cout << "Insufficient statistics on cell " << cell << ", using fit to row." << endl;
	  int idx_r = std::distance(c2cell_row.begin(), it_r);
	  mean_val = c2mean_row[idx_r];
	  if(verb)
	    cout << "Search Row : Found Row - " << (double)r+0.5 << " : " << c2cell_row[idx_r] << endl;
	}
      }

      //default to wide value if no value can be obtained from id or row fits
      if( mean_val == -1000. ){
	cout << "WARNING: Cell " << cell_tdc << ", defaulting to fit over all data." << endl;
	mean_val = mean_wide_val;
      }

      //populate the adct offset array
      mean_vals_tdc[cell_tdc] = mean_val;

      if( channel_override_tdc && cell_tdc==corrChan_tdc )  //manual override option
	writeTDCParFile << corrChanOffset << " ";
      else
	writeTDCParFile << mean_vals_tdc[cell_tdc]/calib_const + tdc_cal.old_param[cell_tdc] - tdc_target/calib_const << " ";

      if(verb)
	cout << "Fitted TDC value for cell " << cell_tdc << ": " << mean_vals_tdc[cell_tdc] << endl;

      //advance the id for the next loop
      cell_tdc++;

    } // endloop over columns 
    writeTDCParFile << endl;

  } // endloop over rows
  writeTDCParFile << endl << endl << endl;

  writeTDCParFile.close();


  fout->Write();
  st->Stop();

  cout << endl << "ADCt Analysis and fits written to " << out_path << endl << endl;
  cout << "New ADCt offsets written to " << new_adctoffset_path << endl << endl;
  cout << "New TDC offsets written to " << new_tdcoffset_path << endl << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
