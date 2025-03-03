//SSeeds 5.3.23 - Standalone script adapted from tcal.C at https://github.com/sebastianseeds/HCal_replay
//NOTE: requires $DB_DIR path set correctly in current environment. Script assumes hardware differences on tdcoffset timestamps and outputs alignment constants for all timestamps within the configuration provided.

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
#include "../include/hcal.h"

// TDC offset extraction constraints
const Int_t first_hcal_chan = 0;
const Int_t total_bins = 320;
const Int_t lower_lim = -140;
const Int_t upper_lim = 20;
const Int_t fit_event_min = 50;
const Double_t observed_tdc_sigma = 2.5; //rough estimate
const Double_t TDC_target = -75.; //Target value for peak tdc. Should be roughly centered in window.

// Set up data structure to hold configuration settings by field setting
typedef struct{
  vector<string> tstamp; //timestamp in database where {0,gain;1,tdcoffsets;2,tdccalib;3,adcoffsets}
  vector<string> run; //all runs
  string gcut; //globalcut
  Double_t SBSf;
  Double_t E_e;
  Double_t HCal_d;
  Double_t HCal_th;
  Double_t BB_th;
  Double_t W2_mean;
  Double_t W2_sig;
  Double_t dx0_n;
  Double_t dx0_p;
  Double_t dy0;
  Double_t dx_sig_n;
  Double_t dx_sig_p;
  Double_t dy_sig;
  Double_t atime0;
  Double_t atime_sig;
  Int_t useAlshield;
  vector<Double_t> otdcoff;
  Double_t tcalib;
} PARA;

// Use nested structure to index parameters by target
typedef struct{
  string targ;
  vector<string> config;
  PARA para[hcal::gMaxfset];
} TARPARA;

// Create structure to define calibration sets for categorical analysis
typedef struct{
  vector<Int_t> target_index;
  vector<Int_t> field_index;
  vector<Double_t> tcalib;
  vector<string> timestamp;
} OFFSET;

//Main <experiment> <configuration> <quasi-replay>; qreplay should only be performed after new offsets obtained
void tdc_align( const char *exp = "gmn", Int_t config=4, bool qreplay = true ){

  // Setup qreplay int index for reporting
  Int_t qreplayIdx = 0;
  if( qreplay )
    qreplayIdx = 1;

  // Get index for timestamp type
  Int_t calidx;
  Int_t offidx;
  for( Int_t s=0; s<hcal::gNstamp; s++ ){
    if( hcal::gStamp[s].compare("tdcoff")==0 )
      offidx=s;
    if( hcal::gStamp[s].compare("tdccalib")==0 )
      calidx=s;
  }

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = hcal::getDate();

  // Set up output analysis file
  TFile *fout = new TFile( Form("outfiles/tdcalign_%s_conf%d_qr%d.root",exp,config,qreplayIdx), "RECREATE" );

  // Declare path and open log file for reporting
  ofstream report;
  string logPath = Form("logs/log_tdcalign_%s_conf%d_qreplay%d.log",exp,config,qreplayIdx);

  report.open( logPath );
  report << "#Error and performance report from " << exp << " config " << config << " obtained " << date.c_str() << endl << endl;

  // Read in meta-config for multiple field settings within a given SBS config
  string metaconfigname = Form("../config/%s/conf%d/%s_conf%d_meta.cfg",exp,config,exp,config);
  ifstream metaconfigfile(metaconfigname);

  // Create array of structures to hold all parameters by target
  TARPARA allp[hcal::gNtar];
  // Create structure to define calibration sets
  OFFSET oset;

  TString mcline;
  for( Int_t t=0; t<hcal::gNtar; t++ ){

    //targ_config[t].targ = hcal::gTar[t];
    allp[t].targ = hcal::gTar[t];
    string end = "end" + hcal::gTar[t];

    while( mcline.ReadLine( metaconfigfile ) && !mcline.BeginsWith(end) ){
      if( !mcline.BeginsWith("#") ){
	string configpath = Form("../config/%s/conf%d/",exp,config);

	configpath += (string)mcline;

	//targ_config[t].config.push_back(configpath);
	allp[t].config.push_back(configpath);
	report << "Recording " << hcal::gTar[t] << " config file: " << configpath << ".." << endl;
      }    
    }
  } 

  cout << endl << "Configuration " << config << " " << exp << " TDC parameters loaded for " << endl;

  // Most recent tdc calibration constant for combined offsets in ns if applicable
  Double_t tdccalib_master;

  // Loop over potential targets
  for( Int_t t=0; t<hcal::gNtar; t++ ){
    string target = hcal::gTar[t];
    Int_t tarfset = allp[t].config.size();

    report << endl << endl << "Total field settings for target " << target << ": " << tarfset << endl;

    // Loop over all available field settings
    for( Int_t f=0; f<tarfset; f++ ){

      report << endl << endl << "Reading in " << exp << " config " << config << " " << target << " field set " << f << " analysis parameters.." << endl;

      // Record timestamps
      ifstream configfile(allp[t].config[f]);
      TString currentline;
      Int_t stampidx=0;
      while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endstamp") ){
	if( !currentline.BeginsWith("#") ){
	  allp[t].para[f].tstamp.push_back((string)currentline);
	  report << "Recording " << hcal::gStamp[stampidx] << " timestamp: " << currentline << ".." << endl;
	  stampidx++;
	}    
      }

      // Record data files
      while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
	if( !currentline.BeginsWith("#") ){
	  allp[t].para[f].run.push_back((string)currentline);
	  report << "Recording file: " << currentline << ".." << endl;
	}    
      }

      // Record globalcuts
      while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
	if( !currentline.BeginsWith("#") ){
	  allp[t].para[f].gcut += (string)currentline;
	  report << "Recording globalcut: " << currentline << ".." << endl;
	}    
      }
      
      // Record physics cuts and parameters
      while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");
	Int_t ntokens = tokens->GetEntries();
	if( ntokens>1 ){
	  TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
	  if( skey == "SBSf" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].SBSf = sval.Atof();
	    report << "Recording beam energy: " << allp[t].para[f].SBSf << endl;
	  }
	  if( skey == "E_e" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].E_e = sval.Atof();
	    report << "Loading beam energy: " << allp[t].para[f].E_e << endl;
	  }
	  if( skey == "HCal_d" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].HCal_d = sval.Atof();
	    report << "Loading HCal distance: " << allp[t].para[f].HCal_d << endl;
	  }
	  if( skey == "HCal_th" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].HCal_th = sval.Atof() * TMath::DegToRad();	
	    report << "Loading HCal angle: " << allp[t].para[f].HCal_th << endl;
	  }
	  if( skey == "BB_th" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].BB_th = sval.Atof() * TMath::DegToRad();	
	    report << "Loading BBCal angle: " << allp[t].para[f].BB_th << endl;
	  }
	  if( skey == "W2_mean" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].W2_mean = sval.Atof();
	    report << "Loading W2 mean cut: " << allp[t].para[f].W2_mean << endl;
	  }
	  if( skey == "W2_sig" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].W2_sig = sval.Atof();
	    report << "Loading W2 sigma cut: " << allp[t].para[f].W2_sig << endl;
	  }
	  if( skey == "dx0_n" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dx0_n = sval.Atof();
	    report << "Loading x position of neutron spot: " << allp[t].para[f].dx0_n << endl;
	  }
	  if( skey == "dx0_p" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dx0_p = sval.Atof();
	    report << "Loading y position of proton spot: " << allp[t].para[f].dx0_p << endl;
	  }
	  if( skey == "dy0" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dy0 = sval.Atof();
	    report << "Loading y position of both hadron spots: " << allp[t].para[f].dy0 << endl;
	  }
	  if( skey == "dx_sig_n" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dx_sig_n = sval.Atof();
	    report << "Loading x sigma of neutron spot: " << allp[t].para[f].dx_sig_n << endl;
	  }
	  if( skey == "dx_sig_p" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dx_sig_p = sval.Atof();
	    report << "Loading x sigma of proton spot: " << allp[t].para[f].dx_sig_p << endl;
	  }
	  if( skey == "dy_sig" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].dy_sig = sval.Atof();
	    report << "Loading y sigma of both hadron spots: " << allp[t].para[f].dy_sig << endl;
	  }
	  if( skey == "atime0" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].atime0 = sval.Atof();
	    report << "Loading ADC time mean: " << allp[t].para[f].atime0 << endl;
	  }
	  if( skey == "atime_sig" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].atime_sig = sval.Atof();
	    report << "Loading ADC time sigma: " << allp[t].para[f].atime_sig << endl;
	  }
	  if( skey == "useAlshield" ){
	    TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	    allp[t].para[f].useAlshield = sval.Atoi();
	    report << "Loading Aluminum absorber option: " << allp[t].para[f].useAlshield << endl;
	  }
	}
	delete tokens;
      } 

      //Record old tdc offsets by tstamp from database
      string DBpath = gSystem->Getenv("DB_DIR");
      TString tdcOffsetPath = DBpath + "/db_sbs.hcal.dat";
      ifstream tdcOffsetFile( tdcOffsetPath );
      report << endl << "Recording previous offsets from database file " << tdcOffsetPath << " with timestamp " << allp[t].para[f].tstamp[offidx] << ".." << endl;
      if( !tdcOffsetFile ){
	cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
	return 0;
      }

      Int_t n0=0;
      Double_t d0=0;
      string line;
      bool read_offset = false;
      bool found_tstamp = false;
      if( allp[t].para[f].tstamp[offidx].compare("none")==0 )
	found_tstamp = true;
      
      while( getline( tdcOffsetFile, line ) ){
	
	if( n0==hcal::maxHCalChan ) break;
	
	TString Tline = (TString)line;
	
	if( Tline.BeginsWith(allp[t].para[f].tstamp[offidx]) && !found_tstamp ){
	  found_tstamp = true;
	  continue;
	}
	
	if( Tline.BeginsWith("sbs.hcal.tdc.offset") && found_tstamp && !read_offset ){
	  read_offset = true;
	  continue;
	}

	if( found_tstamp && read_offset ){
	  istringstream iss( line );
	  while( iss >> d0 ){
	    allp[t].para[f].otdcoff.push_back( d0 );
	    n0++;
	  }
	}
      }
	
      for( Int_t r=0; r<hcal::maxHCalRows; r++){
	for( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  report << allp[t].para[f].otdcoff[i] << " ";
	}
	report << endl;
      }
	
      //Record tdc calibration constant from database
      TString tdcCalibPath = DBpath + "/db_sbs.hcal.dat";
      ifstream tdcCalibFile( tdcCalibPath );
      report << endl << "Recording tdc calibration constant from database file " << tdcCalibPath << " with timestamp " << allp[t].para[f].tstamp[calidx] << ".." << endl;
      if( !tdcCalibFile ){
	cerr << endl << "ERROR: No input constant file present -> path to db_sbs.hcal.dat expected." << endl;
	return 0;
      }

      n0=0;
      d0=0;
      string cline;
      read_offset = false;
      found_tstamp = false;
      if( allp[t].para[f].tstamp[calidx].compare("none")==0 )
	found_tstamp = true;
      
      while( getline( tdcCalibFile, cline ) ){
	
	TString Tline = (TString)cline;

	if( f==0 && t==0 && Tline.BeginsWith("sbs.hcal.tdc.calib") ){
	  string temp;

	  bool val = false;
	  for( Int_t l=0; l<cline.length(); l++){
	    if( val==true )
	      temp+=cline[l];
	    if( cline[l] == '=' ) 
	      val=true;
	  }
	  tdccalib_master = stof( temp );
	  report << "Recorded tdccalib_master: " << tdccalib_master << endl;
	}
	
	if( Tline.BeginsWith(allp[t].para[f].tstamp[calidx]) && !found_tstamp ){
	  found_tstamp = true;
	  continue;
	}

	if( found_tstamp && n0==0 && Tline.BeginsWith("sbs.hcal.tdc.calib") ){
	  string temp;

	  bool val = false;
	  for( Int_t l=0; l<cline.length(); l++){
	    if( val==true )
	      temp+=cline[l];
	    if( cline[l] == '=' ) 
	      val=true;
	  }
	  allp[t].para[f].tcalib = stof( temp );
	  report << "Recorded allp[" << t << "].para[" << f << "].tcalib: " << allp[t].para[f].tcalib << endl;
	  n0++;
	}
      }

      cout << "  TDC calibration constant " << allp[t].para[f].tcalib;

      //Check for changes in calibration set from tdc offset timestamp
      if( std::find(oset.timestamp.begin(),oset.timestamp.end(),allp[t].para[f].tstamp[offidx]) == oset.timestamp.end() ){
	oset.timestamp.push_back( allp[t].para[f].tstamp[offidx] );
	oset.target_index.push_back( t );
	oset.field_index.push_back( f );
	oset.tcalib.push_back( allp[t].para[f].tcalib );
	report << "FOR t " << t << " f " << f << " ADDING timestamp " << oset.timestamp.back() << " for change in tdc offset configuration" << endl;
      }

      //Check for changes in calibration set from tdc calibration constant timestamp
      if( std::find(oset.timestamp.begin(),oset.timestamp.end(),allp[t].para[f].tstamp[calidx]) == oset.timestamp.end() ){
      	oset.timestamp.push_back( allp[t].para[f].tstamp[calidx] );
      	oset.target_index.push_back( t );
      	oset.field_index.push_back( f );
      	oset.tcalib.push_back( allp[t].para[f].tcalib );
      	report << "FOR t " << t << " f " << f << " ADDING timestamp " << " for change in tdc calibration constant configuration" <<  oset.timestamp.back();
      }

      cout << ", target " << hcal::gTar[t] << ", field set " << allp[t].para[f].SBSf << "%.. " << endl;
      usleep( hcal::second );

    }//endloop over field settings
  }//endloop over targets

  Int_t tssize = oset.timestamp.size(); //Number of calibration points

  cout << endl << "Total number of calibration sets: " << oset.timestamp.size() << "." << endl << endl;
  report << "Total number of calibration sets: " << oset.timestamp.size() << "." << endl;

  cout << "Master TDC calibration constant: " << tdccalib_master << endl;

  // Declare paths for new tdc offsets
  string writeParPath = Form("parameters/tdcoffsets_%s_conf%d.txt",exp,config);
  string readParPath = Form("parameters/tdcoffsets_%s_conf%d.txt",exp,config);
  if( qreplay )
    writeParPath = "/dev/null"; //Safety to prevent overwriting constants on quasi-replay

  Double_t ntdcoff[tssize][hcal::maxHCalChan];

  // Recording fit parameters if quasi-replay selected
  if( qreplay ){
    for( Int_t set=0; set<tssize; set++ ){

      ifstream readParFile( readParPath );
      if( !readParFile ){
	cerr << endl << "ERROR: No input constant file present -> parameters/<tdcoffsetfile> expected." << endl;
	return 0;
      }

      report << "Loading newly calculated tdc offset parameters.." << endl;

      // Record new tdc offsets
      Int_t n0=0;
      Double_t d0=0;
      string line;
      bool read_offset = false;
      bool found_tstamp = false;
      if( oset.timestamp[set].compare("none")==0 )
	found_tstamp = true;

      while( getline( readParFile, line ) ){
      
	if( n0==hcal::maxHCalChan ) break;
      
	TString Tline = (TString)line;
  
	if( Tline.BeginsWith(oset.timestamp[set]) && !found_tstamp ){
	  found_tstamp = true;
	  continue;
	}

	if( Tline.BeginsWith("sbs.hcal.tdc.offset") && !read_offset ){
	  read_offset = true;
	  continue;
	}
      
	if( found_tstamp && read_offset ){
	  istringstream iss( line );
	  while( iss >> d0 ){
	    ntdcoff[set][n0]=d0;
	    n0++;
	  }
	}
      }
    
      cout << endl << endl << "New TDC offsets from file: " << endl;
      report << endl << endl << "New TDC offsets from file: " << endl;
    
      for( Int_t r=0; r<hcal::maxHCalRows; r++){
	for( Int_t c=0; c<hcal::maxHCalCols; c++){
	  Int_t i = r*hcal::maxHCalCols+c;
	  cout << ntdcoff[set][i] << " ";
	  report << ntdcoff[set][i] << " ";
	}
	cout << endl;
	report << endl;
      }

      cout << "New TDC offsets loaded for quasi-replay.." << endl;

    }//endloop over sets
  }//endif qreplay

  // Re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // Create simple output tree
  TTree *P = new TTree("P","Analysis Tree");

  // Timing
  Double_t pblkid_out; //hcal primary cluster, primary block id
  Double_t tdc_out; //hcal primary cluster tdc time, tree
  Double_t atime_out; //hcal primary cluster adc time, tree
  Double_t hodotmean_out; //hodoscope primary cluster mean tdc time

  // Physics
  Double_t dx_out; //hcal primary cluster dx
  Double_t dy_out; //hcal primary cluster dy
  Double_t W2_out; //W2
  Double_t Q2_out; //Q2
  Double_t hcale_out; //hcal primary cluster energy
  Int_t mag_out; //sbs magnetic field strength (percent)
  Int_t run_out; //run number
  Int_t tar_out; //target, LH2 or LD2

  // Set output tree branches
  P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
  P->Branch( "tdc", &tdc_out, "tdc/D" );
  P->Branch( "atime", &atime_out, "atime/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "mag", &mag_out, "mag_out/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "tar", &tar_out, "tar_out/I" );

  cout << "Analysis output tree \"P\" created.." << endl;
  report << "Analysis output tree \"P\" created.." << endl;

  // Declare necessary histograms for fits, correct for trigger jitter to get more accurate mean
  TH2D *htp_hodocorr_ID[NTSTAMP];

  for( Int_t set=0; set<tssize; set++ ){
    htp_hodocorr_ID[set] = new TH2D(Form("htp_hodocorr_ID_set%d",set),
				    Form("TDC Primary Block - TDC hodo, set: %s;Channel;TDC_{HCAL}-TDC_{HODO} (ns)",oset.timestamp[set].c_str()),
				    hcal::maxHCalChan,
				    first_hcal_chan,
				    hcal::maxHCalChan,
				    total_bins,
				    lower_lim,
				    upper_lim);

  }

  //////////////////////////
  // Begin looping over data

  // Loop over potential targets
  for( Int_t t=0; t<hcal::gNtar; t++ ){
    
    string target = hcal::gTar[t];
    //Int_t tarfset = targ_config[t].config.size();
    Int_t tarfset = allp[t].config.size();

    // Loop over all available field settings
    for( Int_t f=0; f<tarfset; f++ ){
      
      Double_t mag = allp[t].para[f].SBSf;

      cout << "Looping over data in " << exp << " config " << config << ", target " << target << ", field " << mag << ".."  << endl;
      report << "Looping over data in " << exp << " config " << config << ", target " << target << ", field " << mag << ".."  << endl;
      
      Int_t nruns = allp[t].para[f].run.size();

      for( Int_t r=0; r<nruns; r++ ){
	
	//Add the run(s) from config structure
	C = new TChain("T");
	C->Add(allp[t].para[f].run[r].c_str());
	
	//Get run from string
	string runpathstr = allp[t].para[f].run[r];
	string runstr = "";
	string segstr = "";
	Int_t skipN = 0;
	Int_t rnct = 5; //run number digits

	for( Int_t l=0; l<runpathstr.length(); l++ ){
	  if( runpathstr[l]=='_' && skipN==1 ){ //read 5 digits from the second underscore
	    for( Int_t d=1; d<=rnct; d++ ){
	      runstr += runpathstr[l+d];
	    }
	  }
	  if( runpathstr[l]=='_' && skipN==4 ){ //read all remaining digits from last underscore
	    Int_t b = l;
	    while( runpathstr[b+1]!='.' ){
	      segstr += runpathstr[b+1];
	      b++;
	    }
	  }
	  if( runpathstr[l]=='_' ){ //skip first underscore
	    skipN++;
	  }
	}
 
	Int_t run = stoi( runstr );
	Int_t seg = stoi( segstr );

	report << "Analyzing run " << run << " segment " << seg << ".." << endl;

	// Establish set index to build tdc vs ID plots and quasi-replay with new offsets
	Int_t setidx=-1;

	for( Int_t set=0; set<tssize; set++ ){
	  report << "Set number " << set << " oset: " << oset.timestamp[set] << ", allp: " << allp[t].para[f].tstamp[offidx] << endl;
	  if( oset.timestamp[set].compare(allp[t].para[f].tstamp[offidx]) == 0 )
	    setidx = set;
	}
	if( setidx==-1 ){
	  cout << "ERROR: run DB timestamp not in calibration set on set index " << setidx << ", with parameter timestamp: " << allp[t].para[f].tstamp[offidx] << ".." << endl;
	  return;
	}

	// Get TDC calibration constant for this run
	Double_t tdcCalib = allp[t].para[f].tcalib;

	// Setting up ROOT tree branch addresses
	C->SetBranchStatus("*",0);    
	Double_t BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
	Double_t BBtr_vz[hcal::maxTracks];
	Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;	
	Double_t HCALx, HCALy, HCALe;
	Double_t crow, ccol, nblk;
	Double_t cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan];
	Double_t HODOtmean;
	// Speed up processing by switching on only useful branches
	C->SetBranchStatus( "*", 0 );
	C->SetBranchStatus( "sbs.hcal.x", 1 );
	C->SetBranchStatus( "sbs.hcal.y", 1 );
	C->SetBranchStatus( "sbs.hcal.e", 1 );
	C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
	C->SetBranchStatus( "sbs.hcal.colblk", 1 );
	C->SetBranchStatus( "sbs.hcal.nblk", 1 );
	C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
	C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
	C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
	C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
	C->SetBranchStatus( "bb.tr.n", 1 );
	C->SetBranchStatus( "bb.tr.px", 1 );
	C->SetBranchStatus( "bb.tr.py", 1 );
	C->SetBranchStatus( "bb.tr.pz", 1 );    
	C->SetBranchStatus( "bb.tr.vz", 1 );
	C->SetBranchStatus( "bb.tr.p", 1 );
	C->SetBranchStatus( "bb.ps.e", 1 );
	C->SetBranchStatus( "bb.ps.x", 1 );
	C->SetBranchStatus( "bb.ps.y", 1 );
	C->SetBranchStatus( "bb.sh.e", 1 );
	C->SetBranchStatus( "bb.sh.x", 1 );
	C->SetBranchStatus( "bb.sh.y", 1 );
	C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
	C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
	C->SetBranchStatus( "bb.gem.track.nhits", 1 );
	C->SetBranchStatus( "bb.etot_over_p", 1 );

	// Linking memory
	C->SetBranchAddress( "sbs.hcal.x", &HCALx );
	C->SetBranchAddress( "sbs.hcal.y", &HCALy );
	C->SetBranchAddress( "sbs.hcal.e", &HCALe );
	C->SetBranchAddress( "sbs.hcal.rowblk", &crow );
	C->SetBranchAddress( "sbs.hcal.colblk", &ccol );
	C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total2 number of blocks in highest E clus
	C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
	C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
	C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime ); // Array of block TDC times
	C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times
	C->SetBranchAddress( "bb.tr.n", &BBtr_n );
	C->SetBranchAddress( "bb.tr.px", BBtr_px );
	C->SetBranchAddress( "bb.tr.py", BBtr_py );
	C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
	C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
	C->SetBranchAddress( "bb.tr.p", BBtr_p );
	C->SetBranchAddress( "bb.ps.e", &BBps_e );
	C->SetBranchAddress( "bb.ps.x", &BBps_x );
	C->SetBranchAddress( "bb.ps.y", &BBps_y );
	C->SetBranchAddress( "bb.sh.e", &BBsh_e );
	C->SetBranchAddress( "bb.sh.x", &BBsh_x );
	C->SetBranchAddress( "bb.sh.y", &BBsh_y ); 
	C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );

	TCut GCut = allp[t].para[f].gcut.c_str();
	TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

	// Set up hcal coordinate system with hcal angle wrt exit beamline
	vector<TVector3> hcalaxes; hcal::sethcalaxes( allp[t].para[f].HCal_th, hcalaxes );
	TVector3 hcalorigin = allp[t].para[f].HCal_d*hcalaxes[2] + hcal::HCalvoff*hcalaxes[0];
	Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
	Double_t Eloss_outgoing = hcal::celldiameter/2.0/sin(allp[t].para[f].BB_th) * hcal::tarrho[t] * hcal::dEdx[t];

	long nevent = 0, nevents = C->GetEntries(); 
	Int_t treenum = 0, currenttreenum = 0;

	while (C->GetEntry(nevent++)) {

	  cout << "Analyzing run " << run << ", segment " << seg << ": " << nevent << "/" << nevents << " \r";
	  cout.flush();

	  ///////
	  //Single-loop globalcut method. Save pass/fail for output tree.
	  currenttreenum = C->GetTreeNumber();
	  if( nevent == 1 || currenttreenum != treenum ){
	    treenum = currenttreenum; 
	    GlobalCut->UpdateFormulaLeaves();
	    report << "Updating formula leaves and switching segment at event: " << nevent << endl;
	  }
	  bool failedglobal = GlobalCut->EvalInstance(0) == 0;

	  if( failedglobal ) continue;
	  
	  ///////
	  //Physics calculations
	  //correct beam energy with vertex information, primary track
	  Double_t ebeam_c = allp[t].para[f].E_e - ( (BBtr_vz[0]+hcal::l_tgt[t]/2.0) * hcal::tarrho[t] * hcal::dEdx[t] + hcal::uwallthick[t] * hcal::crho[t] * hcal::cdEdx[t] );

	  TVector3 vertex( 0., 0., BBtr_vz[0] );

	  //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
	  Double_t precon = BBtr_p[0] + Eloss_outgoing;

	  //set up four-momenta with some empty for various calculation methods
	  TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
	  //TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' plvect
	  TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
	  TLorentzVector ptarg; //target momentum
	  ptarg.SetPxPyPzE( 0., 0., 0., hcal::M_t[t] );
	  TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
	  TVector3 qv = q.Vect();
	  TLorentzVector pN; //N' momentum
      
	  //simple calculations for e' and N'
	  Double_t etheta = acos( pe.Pz() / pe.E() );
	  Double_t ephi = atan2( pe.Py(), pe.Px() );
	  Double_t pcent = ebeam_c/( 1. + ( ebeam_c/hcal::M_t[t] )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
	  Double_t phNexp = ephi + hcal::PI;
	  Double_t Q2, W2;

	  //e' p reconstruction with track angles (not momentum)
	  Double_t nu = pbeam.E() - pcent;
	  Double_t pNexp = sqrt( pow(nu, 2.) + 2. * hcal::M_t[t] * nu );
	  Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
	  W2 = pow( hcal::M_t[t], 2.0 ) + 2.0*hcal::M_t[t] * (ebeam_c-pe.E()) - Q2;

	  //Calculate h-arm quantities
	  vector<Double_t> xyhcalexp; hcal::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	  TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
	  Double_t dx = HCALx - xyhcalexp[0];
	  Double_t dy = HCALy - xyhcalexp[1];

	  ///////
	  //BigBite/SBS Acceptance matching.
	  bool failedaccmatch = 
	    xyhcalexp[1] > hcal::posHCalYf ||
	    xyhcalexp[1] < hcal::posHCalYi ||
	    xyhcalexp[0] > hcal::posHCalXf ||
	    xyhcalexp[0] < hcal::posHCalXi;
	  
	  if( failedaccmatch ) continue;

	  Double_t pblkid = (double)cblkid[0]-1;
	  Double_t hcaltdc = cblktime[0];
	  if( qreplay ) 
	    hcaltdc = hcaltdc + allp[t].para[f].otdcoff[(Int_t)pblkid]*tdcCalib - ntdcoff[setidx][(Int_t)pblkid]*tdcCalib;

	  Double_t tdc_tc = hcaltdc-HODOtmean; //Primary cluster, primary block tdc with trigger correction (ns)

	  htp_hodocorr_ID[setidx]->Fill( pblkid, tdc_tc );

	  pblkid_out = pblkid;
	  tdc_out = cblktime[0];
	  atime_out = cblkatime[0];
	  hcale_out = HCALe;
	  dx_out = dx;
	  dy_out = dy;
	  W2_out = W2;
	  Q2_out = Q2;
	  mag_out = mag;
	  run_out = run;
	  tar_out = t;
	  hodotmean_out = HODOtmean;

	  P->Fill();

	}//end loop over event

	// getting ready for the next run
	C->Reset();

      }//endloop over runs
      
    }//endloop over fields
  
  }//endloop over targets
  
  report << "Proceeding to TH2D slices and fits.." << endl;

  // Declare canvas and graph for plots
  TCanvas *TDC_top[tssize];
  TCanvas *TDC_bot[tssize];

  // Set up arrays for tgraph
  //No error on cell location
  Double_t cellerr[hcal::maxHCalChan] = {0.};

  //Make arrays for tdc tgraphs
  Double_t tcell[tssize][hcal::maxHCalChan];
  Double_t tcval[tssize][hcal::maxHCalChan];
  Double_t tcvalw[tssize][hcal::maxHCalChan];
  Double_t tcerr[tssize][hcal::maxHCalChan];
  TH1D *tcellslice[tssize][hcal::maxHCalChan];

  //Get averages for empty channels
  Double_t tcval_avg[tssize];
  Int_t tcval_Ng[tssize];

  for( Int_t set=0; set<tssize; set++ ){
    
    Int_t ti = oset.target_index[set]; //target index
    Int_t fi = oset.field_index[set]; //field index
    Double_t tc = oset.tcalib[set]; //tdc calibration constant
    string ts = oset.timestamp[set]; //timestamp for this index
    tcval_avg[set] = 0.;
    tcval_Ng[set] = 0;

    report << "Looping over timestamp " << ts << ", set " << set << ", target index " << ti << ", field index " << fi << ".." << endl;

    //Set up canvases for quick inspections
    TDC_top[set] = new TCanvas(Form("TDC_top, timestamp: %s",ts.c_str()),"TDC_top",1600,1200);
    TDC_bot[set] = new TCanvas(Form("TDC_bot, timestamp: %s",ts.c_str()),"TDC_bot",1600,1200);
    TDC_top[set]->Divide(12,12);
    TDC_bot[set]->Divide(12,12);
    gStyle->SetOptStat(0);
    
    for(Int_t c=0; c<hcal::maxHCalChan; c++){
      //initialize the graph arrays
      tcvalw[set][c] = 0.;
      tcerr[set][c] = 0.;

      Int_t half_chan = hcal::maxHCalChan/2;

      //Index through the canvas
      TDC_top[set]->cd(c+1);
      if( c>=half_chan )
	TDC_bot[set]->cd(c+1-half_chan);
      

      //Get slices from htDiff_ID and fit for mean vals
      Double_t tfitl = 0.;
      Double_t tfith = 0.;
      tcell[set][c] = c;
      tcellslice[set][c] = htp_hodocorr_ID[set]->ProjectionY(Form("tcellslice_%d_%d",set,c+1),c+1,c+1);
      tcval[set][c] = allp[ti].para[fi].otdcoff[c]; 
      
      tcvalw[set][c] = 0; 

      Int_t sliceN = tcellslice[set][c]->GetEntries();
      Int_t binmax = tcellslice[set][c]->GetMaximumBin();
      Double_t binmaxX = lower_lim+binmax*(upper_lim-lower_lim)/total_bins;
      Double_t binmaxY = tcellslice[set][c]->GetBinContent(binmax);
      tfitl = binmaxX - 4*observed_tdc_sigma;
      tfith = binmaxX + 4*observed_tdc_sigma;

      if( sliceN<fit_event_min || binmaxX<=lower_lim || binmaxX>=upper_lim ){
	tcellslice[set][c]->Draw();
	continue;
      }

      TF1 *sgausfit = new TF1("sgausfit",hcal::g_sgfit_bg,tfitl,tfith,5);
      sgausfit->SetLineWidth(4);
      sgausfit->SetParameter(0,binmaxY);
      sgausfit->SetParameter(1,binmaxX);
      sgausfit->SetParLimits(1,tfitl,tfith);
      sgausfit->SetParameter(2,observed_tdc_sigma);
      sgausfit->SetParLimits(2,1.,3*observed_tdc_sigma);
      sgausfit->SetParameter(3,1.1);
      sgausfit->SetParLimits(3,0.3,8);
      sgausfit->SetParameter(4,25);
      sgausfit->SetParLimits(4,5,30);
      tcellslice[set][c]->Fit("sgausfit","RBMQ");

      tcellslice[set][c]->Draw();
      
      tcval[set][c] = sgausfit->GetParameter(1);
      tcvalw[set][c] = sgausfit->GetParameter(1);
      tcerr[set][c] = sgausfit->GetParameter(2);
      tcellslice[set][c]->SetTitle(Form("Set:%d N:%d MaxX:%f Mean:%f Sigma:%f Alpha:%f",set,sliceN,binmaxX,tcval[set][c],tcerr[set][c],sgausfit->GetParameter(3)));    
      tcellslice[set][c]->Write();

      tcval_avg[set] += tcval[set][c];
      tcval_Ng[set]++;

    }    
    TDC_top[set]->Write();
    TDC_bot[set]->Write();

    tcval_avg[set] /= tcval_Ng[set];

  }//endloop over slice fits

  TCanvas *TDC_graph[tssize];
  TGraphErrors *gtdc_c[tssize];

  // Output text file for new TDC offsets
  ofstream writeParFile;
  writeParFile.open( writeParPath );

  // Get tgraph and record tdc offsets
  for( Int_t set=0; set<tssize; set++ ){

    Int_t ti = oset.target_index[set]; //target index
    Int_t fi = oset.field_index[set]; //field index
    Double_t tc = oset.tcalib[set]; //tdc calibration constant		
    string ts = oset.timestamp[set]; //timestamp for this index

    TDC_graph[set] = new TCanvas(Form("TDC_graph, timestamp: %s",ts.c_str()),"TDC vs ID, means",1600,1200);
    TDC_graph[set]->cd();

    //Make graphs with errors for reporting. All failed fits are zero here
    gtdc_c[set] = new TGraphErrors( hcal::maxHCalChan, tcell[set], tcvalw[set], cellerr, tcerr[set] );
    gtdc_c[set]->GetXaxis()->SetLimits(-10,290);  
    gtdc_c[set]->GetYaxis()->SetLimits(lower_lim,upper_lim);
    gtdc_c[set]->SetTitle(Form("TDC_{hcal}-TDCmean_{hodo} vs Cell, timestamp: %s, qreplay: %d",ts.c_str(),qreplayIdx ) );
    gtdc_c[set]->GetXaxis()->SetTitle("Cell");
    gtdc_c[set]->GetYaxis()->SetTitle("TDC_{HCAL}-TDC_{MEAN,HODO}");
    gtdc_c[set]->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
    gtdc_c[set]->Draw();
    gtdc_c[set]->Write(Form("gtdc_c_s%d",set));
        
    //Write to output text file
    if( qreplay ) continue; //Don't write if quasi-replay
    writeParFile << "#HCal TDC offsets obtained " << date.c_str() << " for db tdc offset timestamp " << ts << endl;
    writeParFile << "#Offsets obtained from fits over TDC distributions." << endl;
    if( oset.timestamp[set].compare("none")!=0 )
      writeParFile << oset.timestamp[set] << endl;
    writeParFile << "sbs.hcal.tdc.offset =" << endl;
    report << endl << endl << "#HCal TDC offsets obtained " << date.c_str() << " for db tdc offset timestamp " << ts << endl;
    report << "#Offsets obtained from fits over TDC distributions = " << endl;
  
    Int_t cell = 0;

    for( Int_t r = 0; r<hcal::maxHCalRows; r++){
      for( Int_t c = 0; c<hcal::maxHCalCols; c++){
	
	if( allp[ti].para[fi].otdcoff[cell]>0 && ( tcval[set][cell]==0 || abs(tcerr[set][cell])<0.05 ) ){ //This will only work if most channels are closely aligned already and oldtdcoffsets aren't believable
	  writeParFile << tcval_avg[set]/tc + allp[ti].para[fi].otdcoff[cell] - TDC_target/tc << " ";
	  report << tcval_avg[set]/tc + allp[ti].para[fi].otdcoff[cell] - TDC_target/tc << " with adjustment.. ";
	}else{
	  writeParFile << tcval[set][cell]/tc + allp[ti].para[fi].otdcoff[cell] - TDC_target/tc << " ";
	  report << tcval[set][cell]/tc + allp[ti].para[fi].otdcoff[cell] - TDC_target/tc << " ";
	}
	cell++;
      } // endloop over columns 
      writeParFile << endl;
      report << endl;
    } // endloop over rows
    writeParFile << endl << endl;
    report << endl << endl;
    
  } // endloop over calibration sets

  writeParFile.close();
  fout->Write();
  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
  report << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  report.close();

  cout << endl << endl << "Detailed processing report can be found here: " << logPath << "." << endl;

}
