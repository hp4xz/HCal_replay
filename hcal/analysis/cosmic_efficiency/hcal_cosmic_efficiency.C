// seeds 9.15.23 - script to determine hcal efficiency from cosmic data as a function of adc threshold. This script assumes a single threshold for all adcs.

#include <iostream>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TMath.h"

const int kNcell = 288; // Total number of HCal modules                      
const int kNrows = 24; // Total number of HCal rows                                         
const int kNcols = 12; // Total number of HCal columns                                 

const double PI = TMath::Pi();
const double Mp = 0.938272;
const double Mn = 0.939565;

//first, write this with a single set of thresholds in mind and produce some efficiency plots
const double thresh_amp = 0.15;
const double thresh_a = 0.05;
const int threshold_size = 10;

void hcal_cosmic_efficiency(const char *configfilename = "hcal_cosmic_efficiency.cfg"){

  // Define a clock to check macro processing time                                                                          
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  TChain *C = new TChain("T");

  //Declare configuration variables
  vector<double> a_thresh; //GeV                                                   
  vector<double> amp_thresh; //GeV                                                   

  // Reading config file                                                                       
  cout << "Opening the following files.." << endl;
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endthresh_a") ){
    if( !currentline.BeginsWith("#") ){
      double thresh = currentline.Atof();
      a_thresh.push_back(thresh);
    }
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endthresh_amp") ){
    if( !currentline.BeginsWith("#") ){
      double thresh = currentline.Atof();
      amp_thresh.push_back(thresh);
    }
  }

  if( threshold_size != a_thresh.size() ){
    cerr << "ERROR: static arrays that govern threshold must match const. Check threshold_size and try again" << endl;
    return;

  }
  
  if( a_thresh.size() != amp_thresh.size() ){
    cerr << "ERROR: amplitude and integrated adc threshold arrays must be equal for comparison" << endl;
    return;
  }
  
   // HCal params
  double hcala[kNcell], hcala_p[kNcell], hcalamp_p[kNcell], hcalamp[kNcell], hcalatime[kNcell], hcaltdc[kNcell];

  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.a_amp_p", 1 );
  C->SetBranchStatus( "sbs.hcal.a_amp", 1 );
  C->SetBranchStatus( "sbs.hcal.a_p", 1 );
  C->SetBranchStatus( "sbs.hcal.a", 1 );
  C->SetBranchStatus( "sbs.hcal.tdc", 1 );
  C->SetBranchStatus( "sbs.hcal.a_time", 1);

  C->SetBranchAddress( "sbs.hcal.a_amp_p", hcalamp_p );
  C->SetBranchAddress( "sbs.hcal.a_amp", hcalamp );
  C->SetBranchAddress( "sbs.hcal.a_p", hcala_p );
  C->SetBranchAddress( "sbs.hcal.a", hcala );
  C->SetBranchAddress( "sbs.hcal.tdc", hcaltdc );
  C->SetBranchAddress( "sbs.hcal.a_time", hcalatime );

  //Outfile                                                                                                                 
  TFile *fout = new TFile("cosmic_eff_out.root","RECREATE");

  //Declare some histograms
  TH2D *hamp_id = new TH2D("hamp_id","ADC amp vs ID",kNcell,0,kNcell,100,0,1);
  TH2D *ha_id = new TH2D("ha_id","int ADC vs ID",kNcell,0,kNcell,300,0,150);
  TH2D *hamp_p_id = new TH2D("hamp_p_id","ped sub ADC amp vs ID",kNcell,0,kNcell,100,0,1);
  TH2D *ha_p_id = new TH2D("h_p_id","ped sub int ADC vs ID",kNcell,0,kNcell,300,-0.5,1);
  TH2D *ha_time_id = new TH2D("ha_time_id","ADC time vs ID",kNcell,0,kNcell,160,-40,120);
  TH2D *htdc_id = new TH2D("htdc_id","TDC vs ID",kNcell,0,kNcell,100,0,1);
    
  //Arrays to calculate efficiencies per threshold
  int detected_a[threshold_size][kNcell] = {0};
  int expected_a[threshold_size][kNcell] = {0}; //will never have any expected on top/bottom row
  int detected_amp[threshold_size][kNcell] = {0};
  int expected_amp[threshold_size][kNcell] = {0};

  //all block variables
  int all_detected_a[threshold_size] = {0};
  int all_expected_a[threshold_size] = {0}; //will never have any expected on top/bottom row
  int all_detected_amp[threshold_size] = {0};
  int all_expected_amp[threshold_size] = {0};

  long Nevents = C->GetEntries();
  long nevent = 0;

  cout << "Processing events.." << endl;

  while( C->GetEntry( nevent++ ) ){

    cout << "Event " << nevent << "/" << Nevents << "\r";
    cout.flush();

    //setup bools for efficiency check over whole detector
    bool nomiss_a[threshold_size] = {true}; 
    bool nomiss_amp[threshold_size] = {true}; 
    bool expect_a[threshold_size] = {false}; 
    bool expect_amp[threshold_size] = {false}; 
    
    //now loop through all blocks and check if a block over threshold exists above and below
    for( int i=0; i<kNcell; ++i ){
      int row = i/kNcols;
      int col = i%kNcols;
      
      //prevent the above and below variables from being outside of assigned memory
      if( row == 0 || row == 23 )
	continue;

      double a_i = hcala[i];
      double amp_i = hcalamp[i];
      double a_p_i = hcala_p[i];
      double amp_p_i = hcalamp_p[i];
      double atime_i = hcalatime[i];
      double tdc_i = hcaltdc[i];

      double a_above = hcala_p[i-12];
      double a_below = hcala_p[i+12];
      double amp_above = hcalamp_p[i-12];
      double amp_below = hcalamp_p[i+12];

      //Fill histograms
      hamp_id->Fill(i,amp_i);
      //cout << "amp:" << amp_i << endl;
      ha_id->Fill(i,a_i);
      //cout << "a_i:" << a_i << endl;
      hamp_p_id->Fill(i,amp_p_i);
      //cout << "amp_p_i:" << amp_p_i << endl;
      ha_p_id->Fill(i,a_p_i);
      //cout << "a_p_i:" << a_p_i << endl;
      ha_time_id->Fill(i,atime_i);
      //cout << "atime_i:" << atime_i << endl;
      htdc_id->Fill(i,tdc_i);
      //cout << "tdc_i:" << tdc_i << endl;
      
      //now loop over all thresholds and fill arrays
      for( size_t t=0; t<a_thresh.size(); ++t ){
	//get current adc thresholds
	double a_threshold = a_thresh[t];
	double amp_threshold = amp_thresh[t];

	//switch bools if hits above and below imply a hit on the current block
	bool expect_hit_a = false;
	bool expect_hit_amp = false;
	
	//Fill expected arrays
	if( a_above>a_threshold && a_below>a_threshold ){
	  expected_a[t][i]++;
	  expect_a[t] = true;
	  expect_hit_a = true;
	}
      
	if( amp_above>amp_threshold && amp_below>amp_threshold ){
	  expected_amp[t][i]++;
	  expect_amp[t] = true;
	  expect_hit_amp = true;
	}

	//fill detected arrays with bools from earlier and current block adc
	if( expect_hit_a && a_i>a_threshold )
	  detected_a[t][i]++;
	else
	  nomiss_a[t] = false;
	if( expect_hit_amp && amp_i>amp_threshold )
	  detected_amp[t][i]++;
	else
	  nomiss_amp[t] = false;

      }//endloop over thresholds
      
    }//endloop over hcal blocks

    //have to loop again to add up all block variables
    for( size_t t=0; t<a_thresh.size(); ++t ){
      
      if( expect_a[t] )
	all_expected_a[t]++;
      
      if( expect_a[t] && nomiss_a[t] )
	all_detected_a[t]++;
      
      if( expect_amp[t] )
	all_expected_amp[t]++;
      
      if( expect_amp[t] && nomiss_amp[t] )
	all_detected_amp[t]++;
    }
      
  }//endloop over events

  //Vectors to hold all efficiency values
  vector<vector<double>> all_eff_a;
  vector<vector<double>> all_eff_amp;
  vector<double> allblock_eff_a;
  vector<double> allblock_eff_amp;
  
  //loop over thresholds, then blocks, get vector of vectors
  for( size_t t=0; t<a_thresh.size(); ++t ){
    // vector<double> efficiencies_a;
    // vector<double> efficiencies_amp;
    // for( int c=0; c<kNcell; ++c ){
    //   efficiencies_a.push_back(detected_a[t][c]/expected_a[t][c]);
    //   efficiencies_amp.push_back(detected_amp[t][c]/expected_amp[t][c]);
    // }
    // all_eff_a.push_back(efficiencies_a);
    // all_eff_a.push_back(efficiencies_amp);

    allblock_eff_a.push_back(all_detected_a[t]/all_expected_a[t]);
    cout << allblock_eff_a[t] << endl;
    allblock_eff_amp.push_back(all_detected_amp[t]/all_expected_amp[t]);
  }

  //get tgraph for all block efficiencies
  int N_ga = a_thresh.size();

  if( N_ga != allblock_eff_a.size() ){
    cout << "ERROR: efficiency and threshold not one-to-one" << endl;
    return;
  }
  
  double *allblock_thresholds_a = &a_thresh[0];
  double *allblock_efficiencies_a = &allblock_eff_a[0];

  cout << "vector sizes and overall " << a_thresh.size() << ":" << allblock_eff_a.size() << ":" << threshold_size << endl;
  for( int i=0; i<threshold_size; ++i ){
    cout << allblock_thresholds_a[i] << ":" << allblock_efficiencies_a[i] << endl;
  }

  TGraph *geff = new TGraph(N_ga,allblock_thresholds_a,allblock_efficiencies_a );
  
  TCanvas *c1 = new TCanvas("c1","hcal cosmic eff via adc amp",1600,700);
  geff->Draw("A");

  fout->Write();

}
