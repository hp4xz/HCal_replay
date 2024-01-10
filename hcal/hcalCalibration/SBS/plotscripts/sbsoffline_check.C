//sseeds 10.1.23 Script to analyze the replay of a single data file with updated database parameters. Extracts hcal energy, sampling fraction, adc time per channel, tdc time per channel. Compares plots to qreplay result.

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../include/sbs.h"

const Double_t minEventPerCell = 100; //minimum elastic events detected in cell to be considered for calibration
const Double_t highDelta = 0.01; //max allowed discrepancy factor between energy deposited in block and expected
const Int_t hcal_first_chan = 0; //first channel in hcal index
const Int_t bins_magnet_percent = 21; //SBS field at 2100A, so divide by this factor to get percent
const Double_t start_magnet_percent = 0.;
const Double_t end_magnet_percent = 105.;
const Int_t total_tdc_bins = 320;
const Int_t lower_tdc_lim = -120;
const Int_t upper_tdc_lim = 40;
const Int_t total_adct_bins = 320;
const Int_t lower_adct_lim = -60;
const Int_t upper_adct_lim = 100;
const Int_t bins_dxdy = 250;
const Double_t start_dx = -4.;
const Double_t end_dx = 3.;
const Double_t start_dy = -1.25;
const Double_t end_dy = 1.25;
const Double_t bins_SFE = 400;
const Double_t start_SFE = 0.;
const Double_t end_SFE = 1.;
const Double_t bins_W2 = 280;
const Double_t start_W2 = 0.;
const Double_t end_W2 = 1.4;
const Int_t linecount = 25;
const Int_t atimeNsig = 6;
const Double_t FWHM = 0.03;
const Double_t E_hard_ulim = 0.95; //all GMn kine other than 11/7
const Double_t SF_hard_ulim = 0.39;
const Double_t E_approx_FWHM = 0.3;
const Double_t SFvrowcol_yrange = 0.12;

void overlayWithGaussianFits(TH2D* hist, TCanvas* canvas, Double_t ymin, Double_t ymax, string axis) {
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
  hist->GetYaxis()->SetRangeUser(ymin,ymax);
  hist->GetYaxis()->SetTitle(axis.c_str());
  hist->GetXaxis()->SetTitle("Channel");
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



void sbsoffline_check( const char *experiment = "gmn", Int_t runN = 11500, Int_t config = 4, Int_t pass = 0 ){  
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
    
  // Get the date
  string date = util::getDate();

  // Get outdir path and declare outfile
  std::string plotdir = Form("../quality_plots/%s/conf%d%s/",experiment,config,"");
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string sbsofflinecheck_path = outdir_path + Form("/hcal_calibrations/qreplay/sbsoffline_check_run%d_pass%d.root",runN,pass);

  // Get path to MC for comparisons
  std::string digmc_path = outdir_path + Form("/hcal_calibrations/MC/hcalE_idx_mc_dig_%s_conf%d.root",experiment,config);
  std::string energy_dmc_path = plotdir + "sbsofflinecheck_energy_dmc.png";
  std::string sf_dmc_path = plotdir + "sbsofflinecheck_sf_dmc.png";
  std::string adctID_path = plotdir + "sbsofflinecheck_adctvID.png";
  std::string tdcID_path = plotdir + "sbsofflinecheck_tdcvID.png";

  //Set up path variables and output files
  TFile *fout = new TFile( sbsofflinecheck_path.c_str(), "RECREATE" );

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
  reportset report_set[hcal::gNmag];
  Int_t report_set_size = 0;
  Int_t report_set_idx = 0;

  //Get experimental configuration parameters
  SBSconfig config_parameters(experiment,config);    
  std::cout << config_parameters;
  Double_t ebeam = config_parameters.GetEbeam();
  Double_t hcaltheta_rad = config_parameters.GetHCALtheta_rad();
  Double_t hcaldist = config_parameters.GetHCALdist();
  Double_t sbsdist = config_parameters.GetSBSdist();
  Double_t bbtheta_rad = config_parameters.GetBBtheta_rad(); //in radians
  std::string sbs_timestamp = config_parameters.GetSBSTimestamp();

  std::cout << "Setting up histograms.." << std::endl;

  TH1D *hW2 = new TH1D("hW2",
		       "W2; GeV",
		       bins_W2,
		       start_W2,
		       end_W2);

  TH1D *hE = new TH1D("hE",
		       "HCal Primary Cluster Energy; GeV",
		       bins_SFE,
		       start_SFE,
		       end_SFE);

  TH1D *hSF = new TH1D("hSF",
		       "HCal Primary Cluster Sampling Fraction",
		       bins_SFE,
		       start_SFE,
		       end_SFE);

  TH2D *htp_hodocorr_ID = new TH2D("htp_hodocorr_ID",
				   "TDC (hodo corrected) vs ID",
				   hcal::maxHCalChan,
				   hcal_first_chan,
				   hcal::maxHCalChan,
				   total_tdc_bins,
				   lower_tdc_lim,
				   upper_tdc_lim);

  TH2D *hap_hodocorr_ID = new TH2D("hap_hodocorr_ID",
				   "ADCt (hodo corrected) vs ID",
				   hcal::maxHCalChan,
				   hcal_first_chan,
				   hcal::maxHCalChan,
				   total_adct_bins,
				   lower_adct_lim,
				   upper_adct_lim);

  TH1D *hdx = new TH1D("hdx",
		       "dx; m",
		       bins_dxdy,
		       start_dx,
		       end_dx);


  TH1D *hdy = new TH1D("hdy",
		       "dy; m",
		       bins_dxdy,
		       start_dy,
		       end_dy);

  //Find some relevant parameters from the runlist
  Int_t mag;
  std::string current_target;
  std::string current_timestamp;
  
  for (Int_t r=0; r<nruns; r++) {
    
    Int_t current_runnumber = runs[r].runnum;

    if( current_runnumber == runN ){
      mag = runs[r].sbsmag / 21; //convert to percent where max field is at 2100A
      current_timestamp = runs[r].adcg_ts;
      current_target = runs[r].target;
      break;
    }

    if(r==nruns-1){
      cout << "Did not find run number in run list, aborting.." << endl;
      return;
    }

  }
  
  //if necessary, convert target string to uppercase
  std::string targ_uppercase = current_target; transform(targ_uppercase.begin(), targ_uppercase.end(), targ_uppercase.begin(), ::toupper );

  //Get run path - hardcoded to my environment, should contain a test replay of a single run, min 100k events
  std::string rootfile_dir = "/volatile/halla/sbs/seeds/";
  //std::string rootfile_path = rootfile_dir + Form("*%d*",runN);
  //std::string rootfile_path = rootfile_dir + Form("e1209019_replayed_13453_stream0_seg0_0_firstevent0_nevent100005.root");

  //std::string rootfile_path = rootfile_dir + "e1209019_fullreplay_11587_stream0_seg0_0.root";
  //std::string rootfile_path = rootfile_dir + "e1209019_fullreplay_13379_stream0_seg0_60.root";

  //std::string rootfile_path = rootfile_dir + Form("e1209019_fullreplay_%d_stream0_seg0_0.root",runN);
  std::string rootfile_path = rootfile_dir + Form("rootfiles/e1209019_fullreplay_%d*.root",runN);


  cout << "Loaded file at " << rootfile_path << endl;

  //Get target configuration
  SBStarget target_parameters(current_target);
  std::cout << target_parameters;
  Int_t target_index = target_parameters.GetTargIndex();  //Target index (1:lh2,2:ld2,3:he3)
  Double_t target_length = target_parameters.GetTargLength();
  Double_t target_rho = target_parameters.GetTargRho();
  Double_t cell_rho = target_parameters.GetCellRho();
  Double_t cell_diam = target_parameters.GetCellDiam();
  Double_t cell_dEdx = target_parameters.GetCelldEdx();
  Double_t upstream_wthick = target_parameters.GetUpstreamWallThick();
  Double_t target_dEdx = target_parameters.GetTargdEdx();
  Double_t M_avg = target_parameters.GetAvgMass();
  
  //Get the relevant index. Note: will need to be updated if more than two calibration sets exist on a single configuration
  Ncal_set_idx = 0;
  std::string active_timestamp;
  util::tsCompare( current_timestamp, config_parameters.GetSBSTimestamp(), active_timestamp );
  //if the current timestamp ever becomes the active timestamp (ts) it is because the current ts from the database is newer than the time when the configuration started. This indicates that a new set of cuts needs to be loaded.
  if( current_timestamp.compare( active_timestamp )==0 )
    Ncal_set_idx = 1;

  //Get available cuts for current config/target/field combination. Use first element (0) of cut
  vector<calcut> cut;
  util::ReadCutList(struct_dir,experiment,config,Ncal_set_idx,pass,current_target,mag,verb,cut);
  std::cout << cut[0];

  // Setting up chain and branch addresses
  TChain *C = new TChain("T");

  C->Add(rootfile_path.c_str());
  
  C->SetBranchStatus("*",0);    
  Double_t BBtr_p[hcal::maxTracks], BBtr_px[hcal::maxTracks], BBtr_py[hcal::maxTracks], BBtr_pz[hcal::maxTracks];
  Double_t BBtr_vz[hcal::maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;	
  Double_t HCALx, HCALy, HCALe;
  Double_t pblkrow, pblkcol, nblk, nclus;
  Int_t Ncid;
  Double_t cblkid[hcal::maxHCalChan], cblke[hcal::maxHCalChan], cblkatime[hcal::maxHCalChan], cblktime[hcal::maxHCalChan], cblkagain[hcal::maxHCalChan];
  Double_t cid[hcal::maxHCalClus], crow[hcal::maxHCalClus], ccol[hcal::maxHCalClus], ce[hcal::maxHCalClus], cx[hcal::maxHCalClus], cy[hcal::maxHCalClus], catime[hcal::maxHCalClus];
  Double_t HODOtmean;

  // Speed up processing by switching on only useful branches
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
  C->SetBranchStatus( "sbs.hcal.colblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  //C->SetBranchStatus( "sbs.hcal.clus_blk.again", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.x", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.y", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.atime", 1 );
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
  C->SetBranchStatus( "Ndata.sbs.hcal.clus.id", 1 ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

  // Linking memory
  C->SetBranchAddress( "sbs.hcal.x", &HCALx );
  C->SetBranchAddress( "sbs.hcal.y", &HCALy );
  C->SetBranchAddress( "sbs.hcal.e", &HCALe );
  C->SetBranchAddress( "sbs.hcal.rowblk", &pblkrow );
  C->SetBranchAddress( "sbs.hcal.colblk", &pblkcol );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk ); // Total number of blocks in highest E clus
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus ); // Total number of clusters
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid ); // kNcell-1 index for each block
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke ); // Array of block energies
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", cblktime ); // Array of block TDC times
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cblkatime ); // Array of block ADC times
  //C->SetBranchAddress( "sbs.hcal.clus_blk.again", cblkagain ); // Array of block ADC gain coeff
  C->SetBranchAddress( "sbs.hcal.clus.id", cid );
  C->SetBranchAddress( "sbs.hcal.clus.row", crow );
  C->SetBranchAddress( "sbs.hcal.clus.col", ccol );
  C->SetBranchAddress( "sbs.hcal.clus.e", ce );
  C->SetBranchAddress( "sbs.hcal.clus.x", cx );
  C->SetBranchAddress( "sbs.hcal.clus.y", cy );
  C->SetBranchAddress( "sbs.hcal.clus.atime", catime );
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
  C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Ncid ); //Odd maxing out at 10 clusters on all cluster Ndata objects, so this is needed in addition to sbs.hcal.nclus

  //globalcut enables
  C->SetBranchStatus( "bb.tr.tg_th", 1 );
  C->SetBranchStatus( "bb.tr.tg_ph", 1 );

  //Use TTreeFormula to avoid looping over data an additional time
  TCut GCut = cut[0].gcut.c_str();

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

  // Set up hcal coordinate system with hcal angle wrt exit beamline
  vector<TVector3> hcalaxes; util::sethcalaxes( hcaltheta_rad, hcalaxes );
  TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal::HCalvoff_p2*hcalaxes[0]; //use pass 2 corrected vertical offset
  Double_t BdL = hcal::maxSBSfield * hcal::sbsdipolegap * (mag/100); //scale crudely by magnetic field
  Double_t Eloss_outgoing = cell_diam/2.0/sin(bbtheta_rad) * target_rho * target_dEdx;

  long nevent = 0, nevents = C->GetEntries(); 
  Int_t treenum = 0, currenttreenum = 0;

  Double_t config_sampling_fraction;
  Double_t config_e_sigma_ratio;

  //Main loop over events in run
  while (C->GetEntry(nevent++)) {

    cout << "Analyzing run " << runN << ": " <<  nevent << "/" << nevents << " \r";
    cout.flush();

    ///////
    //Single-loop elastic globalcut method. Save pass/fail for output tree.
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum; 
      GlobalCut->UpdateFormulaLeaves();
    }
    bool failedglobal = GlobalCut->EvalInstance(0) == 0;
	  
    if( failedglobal ) 
      continue;

    ///////
    //HCal Active Area Cut
    bool failedactivearea = 
      pblkrow==0 || 
      pblkrow==23 || 
      pblkcol==0 || 
      pblkcol==11;

    if( failedactivearea ) 
      continue; //All events with primary cluster element on edge blocks cut

    ///////
    //HCal primary cluster coincidence time cut (using adctime while hcal tdc suspect, new offsets)
    // Int_t pblkid = cblkid[0]-1; //define primary block, primary cluster ID

    // Double_t natime = cblkatime[0]+old_adct_offsets[pblkid]-new_adct_offsets[pblkid]; //new atime
    // Double_t atime0 = cut[0].atime0; //observed elastic peak in adc time
    // Double_t atimesig = cut[0].atime_sig; //observed width of elastic peak in adc time

    // bool failedcoin = abs(natime-atime0)>atimeNsig*atimesig;

    //remove cuts on timing post alignment
    //if( failedcoin ) 
      //continue; //All events where adctime outside of reasonable window cut


    ///////
    //Physics calculations
    //correct beam energy with vertex information, primary track
    Double_t ebeam_c = ebeam - ( (BBtr_vz[0]+target_length/2.0) * target_rho * target_dEdx + upstream_wthick * cell_rho * cell_dEdx );

    TVector3 vertex( 0., 0., BBtr_vz[0] );

    //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
    Double_t precon = BBtr_p[0] + Eloss_outgoing;

    //set up four-momenta with some empty for various calculation methods
    TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
    TLorentzVector pe( precon*BBtr_px[0]/BBtr_p[0], precon*BBtr_py[0]/BBtr_p[0], precon*BBtr_pz[0]/BBtr_p[0], precon ); //e' recon plvect
    TLorentzVector ptarg; //target momentum
    ptarg.SetPxPyPzE( 0., 0., 0., M_avg );
    TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
    TVector3 qv = q.Vect();
    TLorentzVector pN; //N' momentum
      
    //simple calculations for e' and N'
    Double_t etheta = acos( pe.Pz() / pe.E() );
    Double_t ephi = atan2( pe.Py(), pe.Px() );
    Double_t pcent = ebeam_c/( 1. + ( ebeam_c/M_avg )*( 1.0 - cos(etheta) ) ); //e' p reconstructed by angles
    Double_t phNexp = ephi + hcal::PI;
    Double_t Q2, W2;

    //e' p reconstruction with track angles (not momentum)
    Double_t nu = pbeam.E() - pcent;
    Double_t pNexp = sqrt( pow(nu, 2.) + 2. * M_avg * nu );
    Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
    TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );
    pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
    Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
    W2 = pow( M_avg, 2.0 ) + 2.0*M_avg * (ebeam_c-pe.E()) - Q2;

    //Calculate h-arm quantities
    vector<Double_t> xyhcalexp; util::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
    TVector3 hcalpos = hcalorigin + HCALx*hcalaxes[0] + HCALy*hcalaxes[1]; //from primary blk
    Double_t KE_p = nu; //For elastics total predicted by earm
    Double_t SFrac = HCALe/KE_p; //Measured
    Double_t dx = HCALx - xyhcalexp[0];
    Double_t dy = HCALy - xyhcalexp[1];
    TVector3 neutdir = (hcalpos - vertex).Unit();
    Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + hcal::sbsdipolegap/2.0) );
    TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
    Double_t thetapq_p = acos( protdir.Dot( qv.Unit() ) );
    Double_t thetapq_n = acos( neutdir.Dot( qv.Unit() ) );

    //Target Energy in HCal for calibrations
    Double_t hcal_samp_frac = cut[0].hcal_sf; config_sampling_fraction = hcal_samp_frac;
    Double_t hcal_esratio = cut[0].hcal_es; config_e_sigma_ratio = hcal_esratio;
    Double_t KE_exp = KE_p*hcal_samp_frac/hcal_esratio; //Expected E in HCal, equivalent to KE*modified_samp_frac

    ///////
    //BigBite/SBS Acceptance matching. Redundant here with dxdy cuts. Removed.
    bool failedaccmatch = 
      xyhcalexp[1] > hcal::posHCalYf ||
      xyhcalexp[1] < hcal::posHCalYi ||
      xyhcalexp[0] > hcal::posHCalXf ||
      xyhcalexp[0] < hcal::posHCalXi;


    //Fill timing histograms before more stringent cuts
    Double_t pblkid = (double)cblkid[0]-1;
    Double_t hcaltdc = cblktime[0];
    Double_t tdc_tc = hcaltdc-HODOtmean; //Primary cluster, primary block tdc with trigger correction (ns)

    htp_hodocorr_ID->Fill( pblkid, tdc_tc );

    Double_t hcalatime = cblkatime[0];
    Double_t adct_tc = hcalatime-HODOtmean; //Primary cluster, primary block tdc with trigger correction (ns)

    hap_hodocorr_ID->Fill( pblkid, adct_tc );

    ////////////////////////////
    //Primary W2 cut on elastics
    Double_t W2_mean = cut[0].W2_mean;
    Double_t W2_sig = cut[0].W2_sig;
    

    bool failedW2 = fabs(W2-W2_mean)>W2_sig;
    hW2->Fill(W2);

    if( failedW2 ) 
      continue; //Observed mean W2 cut on elastic peak

    hdx->Fill(dx);
    hdy->Fill(dy);

    //proceed to energy with final cuts

    ///////
    //dy elastic cut
    Double_t dy0 = cut[0].dy0;
    Double_t dysig = cut[0].dy_sig;
    bool faileddy = abs(dy-dy0)>3*dysig;

    if( faileddy ) 
      continue;

    ///////
    //PID
    Int_t pid = -1;
    //util::checkPID(current_target,cut[0].dx0_p,cut[0].dx0_n,cut[0].dx_sig_p,cut[0].dx_sig_n,dx,pid);
    //util::checkPID(current_target,0.0,0.07,cut[0].dx_sig_p,cut[0].dx_sig_n,dx,pid);
    util::checkPID(current_target,-0.725,0.098,cut[0].dx_sig_p,cut[0].dx_sig_n,dx,pid);

    ///////
    //dx elastic cut
    if( pid==-1 ) 
      continue;

    Double_t clusE = 0.0;
    for( Int_t blk = 0; blk<nblk; blk++ ){
      Int_t blkid = Int_t(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
      Double_t blke = cblke[blk];
      clusE += blke;
    }

    if( clusE != HCALe )
      cout << endl << "ERROR: primary cluster energy as sum of blocks does not match tree value for total cluster energy" << endl << endl;
    
    Double_t sampling_fraction = clusE/KE_p; //This must be KE_p NOT KE_exp

    hE->Fill( clusE );
    hSF->Fill( sampling_fraction );

  }

  ////
  //Get E and SF comparison plots
  TCanvas *c1 = new TCanvas("c1","mc/data energy",1600,1200);
  c1->cd();
  gStyle->SetOptStat(0);

  TH1D *hE_clone = (TH1D*)hE->Clone("hE_clone");
  Int_t hE_clone_entries = hE_clone->GetEntries();

  //get max value of fit to peak for hE_clone
  Int_t hE_clone_binMax = hE_clone->GetMaximumBin();
  Double_t hE_clone_binCenter = hE_clone->GetBinCenter( hE_clone_binMax );
  Double_t hE_clone_fitLL = hE_clone_binCenter - FWHM;
  Double_t hE_clone_fitUL = hE_clone_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hE_clone_gfit = new TF1("hE_clone_gfit", "gaus");
  hE_clone->Fit("hE_clone_gfit", "QN", "", hE_clone_fitLL, hE_clone_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hE_clone_fit_max = hE_clone_gfit->GetMaximum(hE_clone_fitLL,hE_clone_fitUL);

  //cout << "hE_clone fit max " << hE_clone_fit_max << endl;

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
  Double_t escale_factor = hEmc_fit_max / hE_clone_fit_max;
  hE_clone->Scale(escale_factor);
  Double_t escale_nocut_factor = hEmc_fit_max / hEmc_nocut_fit_max;
  hEmc_nocut->Scale(escale_nocut_factor);
  Double_t hE_clone_max = hE_clone->GetMaximum();

  //cout << hEmc_max << ":" << hE_max << endl;

  Double_t eyrangescale = hEmc_max*1.1;
  if( hEmc_max < hE_clone_max )
    eyrangescale = hE_clone_max*1.1;

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
    
  hE_clone->SetLineWidth(2);
  hE_clone->SetLineColor(kBlack);
  hE_clone->GetXaxis()->SetRangeUser(0,E_hard_ulim);
  hE_clone->GetXaxis()->SetTitle("GeV");
  hE_clone->Draw("hist same");
  
  //Add a legend
  auto elegend = new TLegend(0.33,0.7,0.89,0.89);
  elegend->SetTextSize(0.03);
  elegend->SetHeader(Form("SBS%d SBS-offline Check, E",config));
  elegend->AddEntry(hEmc,Form("MC, %d entries",hEmc_entries),"l");
  elegend->AddEntry(hEmc_nocut,Form("MC no cuts (scale: %0.2f), %d entries",escale_nocut_factor,hEmc_nocut_entries),"l");
  elegend->AddEntry(hE_clone,Form("Run %d (scale: %0.2f), %d entries",runN,escale_factor,hE_clone_entries),"l");
  elegend->Draw();

  c1->SaveAs(energy_dmc_path.c_str());
  c1->Write();

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
  TH1D *hSF_clone = (TH1D*)hSF->Clone("hE_clone");
  Int_t hSF_clone_entries = hSF_clone->GetEntries();

  //TH1D *hSF = (TH1D*)qr_file->Get("hSF_new");
  //Int_t hSF_entries = hSF->GetEntries();

  //get max value of fit
  Int_t hSF_clone_binMax = hSF->GetMaximumBin();
  Double_t hSF_clone_binCenter = hSF->GetBinCenter( hSF_clone_binMax );
  Double_t hSF_clone_fitLL = hSF_clone_binCenter - FWHM;
  Double_t hSF_clone_fitUL = hSF_clone_binCenter + FWHM;
  
  // Fit a Gaussian to this projection
  TF1 *hSF_clone_gfit = new TF1("hSF_clone_gfit", "gaus");
  hSF->Fit("hSF_clone_gfit", "QN", "", hSF_clone_fitLL, hSF_clone_fitUL ); // "Q" for quiet mode, "N" for no draw
  Double_t hSF_clone_fit_max = hSF_clone_gfit->GetMaximum(hSF_clone_fitLL,hSF_clone_fitUL);

  //get scale factor
  Double_t sfscale_factor = hSFmc_fit_max / hSF_clone_fit_max;
  hSF->Scale(sfscale_factor);
  Double_t sfscale_nocut_factor = hSFmc_fit_max / hSFmc_nocut_fit_max;
  hSFmc_nocut->Scale(sfscale_nocut_factor);
  Double_t hSF_clone_max = hSF->GetMaximum();

  Double_t sfyrangescale = hSFmc_max*1.1;
  if( hSFmc_max < hSF_clone_max )
    sfyrangescale = hSF_clone_max*1.1;

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
  sflegend->SetHeader(Form("SBS%d SBS-offline Check, SF",config));
  sflegend->AddEntry(hSFmc,Form("MC, %d entries",hSFmc_entries),"l");
  sflegend->AddEntry(hSFmc_nocut,Form("MC no cuts (scale: %0.2f), %d entries",sfscale_nocut_factor,hSFmc_nocut_entries),"l");
  sflegend->AddEntry(hSF,Form("Run %d (scale: %0.2f), %d entries",runN,sfscale_factor,hSF_clone_entries),"l");
  //sflegend->AddEntry(hSFmc,"MC","l");
  //sflegend->AddEntry(hSF,"Data (Scaled)","l");
  sflegend->Draw();

  c2->SaveAs(sf_dmc_path.c_str());
  c2->Write();

  //ADC time alignments
  TCanvas *c3 = new TCanvas("c3","adct alignments",1600,1200);
  c3->cd();
  gStyle->SetOptStat(0);

  std::string adc_marker = "HCal adct - hodo tmean (ns)";
  overlayWithGaussianFits(hap_hodocorr_ID,c3,-20,20,adc_marker);

  c3->SaveAs(adctID_path.c_str());
  c3->Write();

  //TDC alignments
  TCanvas *c4 = new TCanvas("c4","tdc alignments",1600,1200);
  c4->cd();
  gStyle->SetOptStat(0);

  std::string tdc_marker = "HCal tdc - hodo tmean (ns)";
  overlayWithGaussianFits(htp_hodocorr_ID,c4,-20,20,tdc_marker);

  c4->SaveAs(tdcID_path.c_str());
  c4->Write();

  fout->Write();

  cout << "Analysis complete. Outfile located at " << sbsofflinecheck_path << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    
 

}
