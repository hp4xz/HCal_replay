//sseeds 10.1.23 Script to create pdf of all cuts and experimental parameters used for energy and timing calibrations

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

const Int_t linecount = 22;
const Int_t tlinecount = 12;
const Double_t minEventPerCell = 100; //minimum elastic events detected in cell to be considered for calibration
const Double_t highDelta = 0.01; //max allowed discrepancy factor between energy deposited in block and expected
const Double_t timing_fit_ev_min = 50;

void plot_report( const char *experiment = "gmn", Int_t config = 7, Int_t pass = 0 ){
  
   
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  bool h2only = true; //only consider LH2 when calibrating HCal to constrain fermi smearing
  // outfile path
  std::string h2opt = "";
  if( h2only )
    h2opt = "_lh2only";

  // Get the date
  string date = util::getDate();

  //Set up path variables and output files
  string db_path = gSystem->Getenv("DB_DIR");
  std::string new_adcgain_path = Form("parameters/adcgaincoeff_%s%s_conf%d_pass%d.txt",experiment,h2opt.c_str(),config,pass);

  // Get information from .csv files
  std::string struct_dir = Form("../config/%s/",experiment); //unique to my environment for now
  Int_t nruns = -1; //Analyze all available runs for this configuration
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<calrun> runs; 
  util::ReadRunList(struct_dir,experiment,nruns,config,pass,verb,runs); //nruns=-1 modifies nruns to loop over all available
  
  //define structures for holding calibration data and indices
  calset gain_coeff[hcal::gNstamp];
  Int_t Ncal_set_size = 0;
  Int_t Ncal_set_idx = 0;

  //define structures for holding cut report data and indices	    
  reportset gain_report_set[hcal::gNmag];
  Int_t gain_report_set_size = 0;
  Int_t gain_report_set_idx = 0;

  //define structures for holding cut report data and indices	    
  reportset timing_report_set[hcal::gNmag];
  Int_t timing_report_set_size = 0;
  Int_t timing_report_set_idx = 0;

  //define array for modified sampling fraction
  Double_t sfracstar[hcal::gNmag];

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  //reporting indices
  Int_t target_change_index;
  Double_t config_sampling_fraction;
  Double_t config_e_sigma_ratio;
  bool first = true;
  Double_t test_adcg_coeff[hcal::maxHCalChan] = {1.};

  //Main loop over energy runs
  for (Int_t r=0; r<nruns; r++) {

    //Get run experimental parameters
    std::string current_timestamp = runs[r].adcg_ts;
    Int_t current_runnumber = runs[r].runnum;

    std::string current_target = runs[r].target;
    std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );
    Int_t mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A

    //Get run paths
    std::string rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-%s/pass%d/SBS%d/%s/rootfiles/",experiment,pass,config,targ_uppercase.c_str());
    std::string rootfile_path = rootfile_dir + Form("*%d*",current_runnumber);

    //Get target configuration
    SBStarget target_parameters(current_target);
    Int_t target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
    Double_t target_length = target_parameters.GetTargLength();
    Double_t target_rho = target_parameters.GetTargRho();
    Double_t cell_rho = target_parameters.GetCellRho();
    Double_t cell_diam = target_parameters.GetCellDiam();
    Double_t cell_dEdx = target_parameters.GetCelldEdx();
    Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
    Double_t target_dEdx = target_parameters.GetTargdEdx();
    Double_t M_avg = target_parameters.GetAvgMass();
        
    //check target and continue on deuterium of only calibrating with hydrogen
    if( h2only && target_index!=1 )
      continue;

    //Record old gain and offset parameters by tstamp from database (assumes one file for all timestamps)
    std::string old_db_path = db_path + "/db_sbs.hcal.dat";
    std::string new_adct_path = Form("../timing/parameters/adctoffsets_class_%s_conf%d_pass%d.txt",experiment,config,pass);
    std::string db_gain_variable = "sbs.hcal.adc.gain";
    std::string db_adct_variable = "sbs.hcal.adc.timeoffset";
    Double_t new_adct_offsets[hcal::maxHCalChan] = {0.};
    Double_t old_adct_offsets[hcal::maxHCalChan] = {0.};
    Double_t new_adcg_coeff[hcal::maxHCalChan] = {1.};
    util::readDB( old_db_path, runs[r].adct_ts, db_adct_variable, old_adct_offsets );
    //get new adct offsets
    std::string adct_active_timestamp;
    util::tsCompare( config_parameters.GetSBSTimestamp(), runs[r].adct_ts, adct_active_timestamp );
    util::readDB( new_adct_path, adct_active_timestamp, db_adct_variable, new_adct_offsets );

    //Look through available calibration sets by gain coeff timestamp for duplicate entries. Correct index where necessary.
    bool found_ts = false;
    for( Int_t cs=0; cs<=Ncal_set_size; cs++ ){
      if( gain_coeff[cs].timestamp == current_timestamp ){
	found_ts = true;
	Ncal_set_idx = cs;
      }
    }
    
    // If no, add new calibration set and correct the index.
    if( !found_ts ){
      Ncal_set_idx = Ncal_set_size;

      if( Ncal_set_size == hcal::gNstamp-1 ){
	cout << "Error: Allowable calibration sets exceed maximum of " << hcal::gNstamp << ". Check timestamps and adjust limits as necessary." << endl;
	return;
      }

      gain_coeff[Ncal_set_size].timestamp = current_timestamp.c_str();
      gain_coeff[Ncal_set_size].good_events = 0;

      //resize the matrices and vectors
      gain_coeff[Ncal_set_size].Ma.ResizeTo(hcal::maxHCalChan,hcal::maxHCalChan);
      gain_coeff[Ncal_set_size].Ma_oneblock.ResizeTo(hcal::maxHCalChan,hcal::maxHCalChan);
      gain_coeff[Ncal_set_size].ba.ResizeTo(hcal::maxHCalChan);
      gain_coeff[Ncal_set_size].ba_oneblock.ResizeTo(hcal::maxHCalChan);

      util::readDB( old_db_path, current_timestamp, db_gain_variable, gain_coeff[Ncal_set_size].old_param );      

      Ncal_set_size++;
      
    }

    //quick check on gain coefficients
    //Double_t test_adcg_coeff[hcal::maxHCalChan] = {1.};
    bool different = false;
    if( first )
      std::cout << "First calibration set, Ncal_set_idx: " << Ncal_set_idx << std::endl << std::endl;  
    for( int r=0; r<24; r++ ){
      for ( int c=0; c<12; c++){
	int i = r*12+c;
	if( first ){
	  test_adcg_coeff[i] = new_adcg_coeff[i];
	}else{
	  if( test_adcg_coeff[i] != new_adcg_coeff[i] )
	    different = true;
	}
      }
    }

    if( different ){
      std::cout << "ADC gain coeffients have changed on calibration set " << Ncal_set_idx << std::endl << std::endl;
      for( int r=0; r<24; r++ ){
	for ( int c=0; c<12; c++){
	  int i = r*12+c;
	  test_adcg_coeff[i] = new_adcg_coeff[i];
	}
      }
    }

    //Get available cuts for current config/target/field combination. Use first element (0) of cut
    vector<calcut> cut;
    util::ReadCutList(struct_dir,experiment,config,Ncal_set_idx,pass,current_target,mag,verb,cut);
    if( different || first )
      std::cout << cut[0];
    
    //Remove first calibration set trigger for future processing
    if( first )
      first = false;

    //Add globalcut and elastic cuts for reporting
    if( target_change_index != target_index ){
      std::cout << target_parameters;
      target_change_index = target_index;
    }

    //Look through available report sets by magnetic field setting and target for duplicate entries. Correct index where necessary.
    bool found_ri = false;
    for( Int_t rs=0; rs<=gain_report_set_size; rs++ ){
      if( gain_report_set[rs].mag == mag && gain_report_set[rs].targetidx == target_index ){
	found_ri = true;
	gain_report_set_idx = rs;
      }
    }

    if( !found_ri ){
      gain_report_set_idx = gain_report_set_size;

      gain_report_set[gain_report_set_idx].gcut = cut[0].gcut;
      gain_report_set[gain_report_set_idx].target = current_target;
      gain_report_set[gain_report_set_idx].mag = mag;
      gain_report_set[gain_report_set_idx].targetidx = target_index;
      gain_report_set[gain_report_set_idx].W2mean = cut[0].W2_mean;
      gain_report_set[gain_report_set_idx].W2sigma = cut[0].W2_sig;
      gain_report_set[gain_report_set_idx].dxmean_n = cut[0].dx0_n;
      gain_report_set[gain_report_set_idx].dxmean_p = cut[0].dx0_p;
      gain_report_set[gain_report_set_idx].dxsigma_n = cut[0].dx_sig_n;
      gain_report_set[gain_report_set_idx].dxsigma_p = cut[0].dx_sig_p;
      gain_report_set[gain_report_set_idx].dymean = cut[0].dy0;
      gain_report_set[gain_report_set_idx].dysigma = cut[0].dy_sig;
      gain_report_set[gain_report_set_idx].atimemean = cut[0].atime0;
      gain_report_set[gain_report_set_idx].atimesigma = cut[0].atime_sig;
      gain_report_set[gain_report_set_idx].minEv = minEventPerCell;
      gain_report_set[gain_report_set_idx].highdelta = highDelta;
    
      gain_report_set_size++;

      sfracstar[gain_report_set_idx] = cut[0].hcal_sf * cut[0].hcal_es;

    }


  }

  // Plot output directory
  std::string plotdir = "../quality_plots/reports/";

  //Add output report canvas
  TCanvas *c1[gain_report_set_size];
  Int_t gain_report_height = 1800;

  for( Int_t s=0; s<gain_report_set_size; s++ ){

    c1[s] = new TCanvas(Form("ecalgain_report%d",s), Form("Gain Configuration/Cut Information, Report %d",s), gain_report_height, 1000);
  
    // Set margin.
    c1[s]->SetLeftMargin(0.01);

    // Create a TText object.
    TText *t = new TText();

    // Set text align to the left (horizontal alignment = 1).
    t->SetTextAlign(11);
    t->SetTextSize(0.02);

    //make an array of strings
    std::string report[linecount] = {
      "General HCal Energy Calibration Info",
      Form("Experiment: %s, Configuration: %d, Pass: %d", experiment, config, pass),
      Form("Creation Date: %s", date.c_str() ),
      Form("Target: %s", gain_report_set[s].target.c_str() ),
      Form("SBS Field: %d%%", gain_report_set[s].mag ),
      "",
      "Elastic Cuts",
      Form("Global Elastic Cuts: %s", gain_report_set[s].gcut.c_str() ),
      Form("W2 mean (GeV): %.2f", gain_report_set[s].W2mean),
      Form("W2 sigma (GeV): %.2f", gain_report_set[s].W2sigma),
      Form("dx mean, proton (m): %.2f", gain_report_set[s].dxmean_p),
      Form("dx sigma, proton (m): %.2f", gain_report_set[s].dxsigma_p),
      Form("dy mean (m): %.2f", gain_report_set[s].dymean),
      Form("dy sigma (m): %.2f", gain_report_set[s].dysigma),
      Form("adc time mean (ns): %.2f", gain_report_set[s].atimemean),
      Form("adc time sigma (ns): %.2f", gain_report_set[s].atimesigma),
      "",
      "Other Cuts/Information",
      Form("Minimum Ev per Cell : %d", gain_report_set[s].minEv),
      Form("Minimum Energy Deposited in Cell (factor, vs expectation) : %.2f", gain_report_set[s].highdelta),
      Form("Sampling Fraction Target, Modified MC: %.4f", sfracstar[s]),
      "HCal Active Area (Projected Nucleon 1 row/col Within HCal Acceptance)"
    };
    // Loop to write the lines to the canvas.
    for( Int_t i = 0; i<linecount; i++ ) {
      // Vertical position adjusted according to line number.
      Double_t verticalPosition = 0.95 - i * 0.03;
      t->DrawTextNDC(0.1, verticalPosition, report[i].c_str());
    }


    std::string adcgain_report_path = plotdir + Form("adcgain_%s_report_sbs%d_set%d_pass%d.png",experiment,config,s,pass);
    c1[s]->SaveAs(adcgain_report_path.c_str());

  }

  TCanvas *c2 = new TCanvas("timing_report", "Timing Configuration/Cut Information, Report", gain_report_height, 500);
  
  

  // Set margin.
  c2->SetLeftMargin(0.01);

  // Create a TText object.
  TText *t = new TText();

  // Set text align to the left (horizontal alignment = 1).
  t->SetTextAlign(11);
  t->SetTextSize(0.04);

  //make an array of strings
  std::string report[tlinecount] = {
    "General HCal Timing Calibration Info",
    Form("Experiment: %s, Configuration: %d, Pass: %d", experiment, config, pass),
    Form("Creation Date: %s", date.c_str() ),
    "Target(s) Used: All Available",
    Form("Calibration Set: %s", config_parameters.GetSBSTimestamp().c_str()),
    "",
    "Elastic Cuts",
    Form("Global Elastic Cuts: %s", gain_report_set[0].gcut.c_str() ),
    "",
    "Other Cuts/Information",
    Form("Minimum Ev per Cell : %.1f", timing_fit_ev_min),
    "HCal Active Area (Projected Nucleon 1 row/col Within HCal Acceptance)"
  };
  // Loop to write the lines to the canvas.
  for( Int_t i = 0; i<tlinecount; i++ ) {
    // Vertical position adjusted according to line number.
    Double_t verticalPosition = 0.95 - i * 0.06;
    t->DrawTextNDC(0.1, verticalPosition, report[i].c_str());
  }


  std::string timing_report_path = plotdir + Form("timing_%s_report_sbs%d_pass%d.png",experiment,config,pass);
  c2->SaveAs(timing_report_path.c_str());


}
