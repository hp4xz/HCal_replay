//seeds Oct 10, waveform check with cuts on adct difference between pblk and all secondary blks in pclus. Assumes min sample is zero (left side of adc window)

#include <iostream>
#include <sstream>
#include <fstream>
#include <TMath.h>
#include <cstdio>
#include <cstdlib>
#include "TLatex.h"

const int maxGB = 1000;
const int maxPB = 1000;
const int maxSamps = 20000;
const int maxSHisto = 20;
const int maxSampBins = 40;

const double hcalheight = -0.2897;
const int hcalrows = 24;
const int hcalcols = 12;
const int maxTracks = 16; 
const double PI = TMath::Pi();
const double M_e = 0.0051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const int total_bins_short = 480;
const int lower_lim_short = -120;
const int upper_lim_short = 120;

const double thresh = 0.10;

int CountLocalMaxima(TH1D* hist) {
    if (!hist) {
        // Handle the case where the input histogram is nullptr
        return 0;
    }

    int numMaxima = 0;
    int numBins = hist->GetNbinsX();

    // Loop over the bins in the histogram
    for (int bin = 2; bin < numBins; bin++) {
        double binContent = hist->GetBinContent(bin);
        double prevBinContent = hist->GetBinContent(bin - 1);
        double nextBinContent = hist->GetBinContent(bin + 1);

        // Check if the current bin has a higher content than its neighbors
        if (binContent > prevBinContent && binContent > nextBinContent) {
            numMaxima++;
        }
    }

    return numMaxima;
}

int CountLocalMaximaAboveThreshold(TH1D* hist) {
    if (!hist) {
        // Handle the case where the input histogram is nullptr
        return 0;
    }

    int numMaxima = 0;
    int numBins = hist->GetNbinsX();
    double maxBinContent = hist->GetMaximum() * thresh; // 10 percent threshold

    // Loop over the bins in the histogram
    for (int bin = 2; bin < numBins; bin++) {
        double binContent = hist->GetBinContent(bin);
        double prevBinContent = hist->GetBinContent(bin - 1);
        double nextBinContent = hist->GetBinContent(bin + 1);

        // Check if the current bin has a higher content than its neighbors
        if (binContent > prevBinContent && binContent > nextBinContent && binContent >= maxBinContent) {
            numMaxima++;
        }
    }

    return numMaxima;
}

bool adjacent(int seed, int block){

  int seed_row = seed/hcalcols;
  int seed_col = seed%hcalcols;
  int block_row = block/hcalcols;
  int block_col = block%hcalcols;

  bool rowclose = abs(block_row-seed_row)<=1;
  bool colclose = abs(block_col-seed_col)<=1;

  return rowclose && colclose;

}

//Main. Configured for GMn. Args->sbs=kinematic, tdiff_cut=cut in ns around pblk-sblk peak, select_bidx=selected tdiff sblk index to populate example histograms (-1 allows all), rootfilename=replayed root file with hcal samples
void adct_diff_waveform_check(int sbs = 11, double tdiff_cut = 9., int select_bidx = -1, string rootfilename = "hcal_gmn_replayed_12916_stream0_seg0_0_firstevent0_nevent6000.root"){//main
 
  //hcal_gmn_replayed_12916_stream0_seg0_0_firstevent0_nevent6000.root emin=0.005, tmax=1000.

  TChain *C = new TChain("T");

  string rootfilepath = "/lustre19/expphy/volatile/halla/sbs/seeds/" + rootfilename;

  string outputfilename = Form("root_out/tdiff_waveforms_sbs%d.root",sbs);

  C->Add(rootfilepath.c_str());

  cout << "Added rootfile " << rootfilepath << " to chain." << endl;

  double E_e;
  double BB_d;
  double BB_th;
  double HCal_d;
  double HCal_th;
  double W2_mean;
  double W2_sig;
  double dy_mean;
  double dy_sig;

  //SBS 4
  if(sbs==4){
    E_e = 3.7393; //Beam energy in GeV
    BB_d = 1.7988; //Distance to bigbite from target in m
    BB_th = TMath::DegToRad() * 36.0; //Angle wrt beamline for bigbite arm in radians
    HCal_d =11.0; //Distance to HCal from target in m
    HCal_th = TMath::DegToRad() * 31.9; //Angle wrt beamline for HCal in radians
    W2_mean = 1.00; //W^2 location of the mean of the elastic peak
    W2_sig = 0.24; //W^2 width of the elastic peak
    dy_mean = -1.687; //Elastic peak in dy mean
    dy_sig = 0.69; //Elastic peak in dy sigma
  }

  // //SBS 7
  if(sbs==7){
    E_e = 7.9072;
    BB_d = 1.850;
    BB_th = 40.0*TMath::DegToRad();
    HCal_d =14.0;
    HCal_th = 16.1*TMath::DegToRad();
    W2_mean = 0.98;
    W2_sig = 0.35;
    dy_mean = -0.0276561;
    dy_sig = 0.283612;
  }

  //SBS 11
  if(sbs==11){

    E_e = 9.8594;
    BB_d = 1.550;
    BB_th = 42.0*TMath::DegToRad();
    HCal_d =14.5;
    HCal_th = 13.3*TMath::DegToRad();
    W2_mean = 0.98;
    W2_sig = 0.35;
    dy_mean = 0.015;
    dy_sig = 0.083;
  }

  //SBS 14
  if(sbs==14){

    E_e = 5.9649;
    BB_d = 1.850;
    BB_th = 46.5*TMath::DegToRad();
    HCal_d =14.0;
    HCal_th = 17.3*TMath::DegToRad();
    W2_mean = 0.98;
    W2_sig = 0.35;
    dy_mean = 0.001;
    dy_sig = 0.045;
  }
  
  //SBS 8
  if(sbs==8){

    E_e = 5.965;
    BB_d = 1.97473;
    BB_th = TMath::DegToRad() * 26.5;
    HCal_d =11.0;
    HCal_th = TMath::DegToRad() * 29.4;
    W2_mean = 1.00;
    W2_sig = 0.24;
    dy_mean = -1.687;
    dy_sig = 0.69;
  }

  //SBS 9
  if(sbs==9){

    E_e = 4.013;
    BB_d = 1.550;
    BB_th = 49.0*TMath::DegToRad();
    HCal_d =11.0;
    HCal_th = 22.0*TMath::DegToRad();
    W2_mean = 0.98;
    W2_sig = 0.35;
    dy_mean = -0.0276561;
    dy_sig = 0.283612;
  }

  //reasonable globalcut (wide) for all kinematics
  //TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&sbs.hcal.nclus>0&&bb.sh.nclus>0";

  if(sbs==-1){
    cout << "ERROR: No valid sbs runnumber given." << endl;
    return;
  }

  //Manage TTree
  double HCALx, HCALy, HCALe, ekineW2;
  double hcalGB_id[maxGB], hcalGB_cid[maxGB];
  int NhcalGB_id,Ndata;

  double hcalB_id[maxPB], cblkatime[maxPB], cblktime[maxPB];
  double hcal_samps[maxSamps], hcal_sampsidx[maxPB], hcal_nsamps[maxPB], hcal_adc[maxPB], hcal_tdc[maxPB], hcal_row[maxPB], hcal_col[maxPB]; 
  //double hcal_sampsid[maxPB];
  double nblk, index, nclus;

  C->SetBranchStatus("*",0);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);

  C->SetBranchStatus("sbs.hcal.nblk",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.id",1);
  C->SetBranchStatus("sbs.hcal.index",1);
  C->SetBranchStatus("sbs.hcal.nclus",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.atime", 1);
  C->SetBranchStatus("sbs.hcal.clus_blk.tdctime", 1);
  C->SetBranchStatus("sbs.hcal.samps", 1);
  C->SetBranchStatus("sbs.hcal.samps_idx", 1);
  C->SetBranchStatus("sbs.hcal.nsamps", 1);
  //C->SetBranchStatus("sbs.hcal.samps_elemID", 1);
  C->SetBranchStatus("sbs.hcal.a", 1);
  C->SetBranchStatus("sbs.hcal.tdc", 1);
  C->SetBranchStatus("sbs.hcal.adcrow", 1);
  C->SetBranchStatus("sbs.hcal.adccol", 1);
  C->SetBranchStatus("Ndata.sbs.hcal.adcrow",1);
  //C->SetBranchStatus("Ndata.sbs.hcal.samps_elemID", 1);

  C->SetBranchAddress("sbs.hcal.x", &HCALx);
  C->SetBranchAddress("sbs.hcal.y", &HCALy);
  C->SetBranchAddress("sbs.hcal.e", &HCALe);

  C->SetBranchAddress("sbs.hcal.clus_blk.id", hcalB_id);
  C->SetBranchAddress("sbs.hcal.nblk", &nblk);
  C->SetBranchAddress("sbs.hcal.index", &index);
  C->SetBranchAddress("sbs.hcal.nclus", &nclus);
  C->SetBranchAddress("sbs.hcal.clus_blk.atime", cblkatime); // Array of block ADC times
  C->SetBranchAddress("sbs.hcal.clus_blk.tdctime", cblktime); // Array of block ADC times
  C->SetBranchAddress("sbs.hcal.samps",hcal_samps);
  C->SetBranchAddress("sbs.hcal.samps_idx",hcal_sampsidx);
  C->SetBranchAddress("sbs.hcal.nsamps",hcal_nsamps);
  //C->SetBranchAddress("sbs.hcal.samps_elemID",hcal_sampsid);
  C->SetBranchAddress("sbs.hcal.a",hcal_adc);
  C->SetBranchAddress("sbs.hcal.tdc",hcal_tdc);
  C->SetBranchAddress("sbs.hcal.adcrow",hcal_row);
  C->SetBranchAddress("sbs.hcal.adccol",hcal_col);
  C->SetBranchAddress("Ndata.sbs.hcal.adcrow",&Ndata);
  //C->SetBranchAddress("Ndata.sbs.hcal.samps_elemID",&Ndata);

  TFile *fout = new TFile(outputfilename.c_str(),"RECREATE");

  TH1D *ht_all_diff = new TH1D("ht_all_diff",
			       "adct all channels/runs all cluster block difference",
			       total_bins_short,
			       lower_lim_short,
			       upper_lim_short);

  TH1D *ht_all_tdiff = new TH1D("ht_all_tdiff",
			       "tdc all channels/runs all cluster block difference",
			       total_bins_short,
			       lower_lim_short,
			       upper_lim_short);

  TH1D *hnmaxima_pblk = new TH1D("hnmaxima_pblk",
			       "N local maxima in pblk waveform, tdiff cut blk",
			       20,
			       0,
			       20);

  TH1D *hnmaxima_gblk = new TH1D("hnmaxima_gblk",
			       "N local maxima in gblk waveform, tdiff cut blk",
			       20,
			       0,
			       20);

  TH1D *hnmaxima_bblk = new TH1D("hnmaxima_bblk",
			       "N local maxima in bblk waveform, tdiff anticut blk",
			       20,
			       0,
			       20);

  //samps histograms
  TH1D *hWaveform_pblk[maxSHisto+10];
  TH1D *hWaveform_gblk[maxSHisto+10];
  TH1D *hWaveform_bblk[maxSHisto+10];

  for( int i=0; i<maxSHisto; ++i ){
    hWaveform_pblk[i] = new TH1D(Form("hWaveform_pblk_ex%d",i),
				Form("hcal adct pblk waveform, example %d",i),
				maxSampBins,
				0,
				maxSampBins);

    hWaveform_gblk[i] = new TH1D(Form("hWaveform_gblk_ex%d",i),
				Form("hcal adct sblk waveform, tdiff cut, example %d",i),
				maxSampBins,
				0,
				maxSampBins);

    hWaveform_bblk[i] = new TH1D(Form("hWaveform_bblk_ex%d",i),
				Form("hcal adct sblk waveform, tdiff anticut, example %d",i),
				maxSampBins,
				0,
				maxSampBins);
  }

  Long64_t Nevents = C->GetEntries();
  int pblk_counter = 0;
  int gblk_counter = 0;
  int bblk_counter = 0;
  bool fullhistos = false;

  cout << Nevents << " in tree." << endl;

  for( Long64_t nevent = 1; nevent <Nevents; nevent++){

    // cout << " Entry = " << nevent << "/" << Nevents << ", histos full: " << fullhistos << "\r";
    // cout.flush();
    
    C->GetEntry(nevent); 

    //coin
    double pblkatime = cblkatime[0];
    double pblktime = cblktime[0];
    double pblkid = hcalB_id[0]-1;
    int pblkrow = pblkid/hcalcols;
    int pblkcol = (int)pblkid%hcalcols;

    bool hcalaa = pblkrow!=0 && pblkrow!=23 && pblkcol!=0 && pblkcol!=11;

    if( !hcalaa )
      continue;

    //Quality checks on tmax, start on block index 1 (second block) to avoid comparison between pblk and pblk
    vector<int> goodBlockID;
    vector<int> badBlockID;
    for( int b=1; b<nblk; b++ ){

      if( b!=select_bidx && select_bidx!=-1 )
	continue;

      double blkatime = cblkatime[b];
      double blktime = cblktime[b];
      double blkid = hcalB_id[b]-1;
      int blkrow = blkid/hcalcols;
      int blkcol = (int)blkid%hcalcols;

      double diff = pblkatime - blkatime;
      double tdiff = pblktime - blktime;
      ht_all_diff->Fill(diff);
      if(tdiff!=0)
	ht_all_tdiff->Fill(tdiff);

      bool failedtdiff = abs(diff)>tdiff_cut;

      if(failedtdiff)
	badBlockID.push_back((int)blkid);
      else
	goodBlockID.push_back((int)blkid);

      bool close = adjacent((int)pblkid,(int)blkid);
      
      int close_idx=-1;
      if(close){
	if( pblkrow-blkrow == 1 )
	  close_idx = blkid-pblkid+14;
	else if( pblkrow-blkrow == 0 )
	  close_idx = blkid-pblkid+5;
	else if( pblkrow-blkrow ==-1 )
	  close_idx = pblkid-blkid+20;
      }

    }//endloop over cluster blocks

    //Check if enough example waveforms have been saved
    if( pblk_counter>=maxSHisto && gblk_counter>=maxSHisto && bblk_counter>=maxSHisto ){
      fullhistos=true;
      continue;
    }

    //Check if no block IDs pass selection criteria
    if(goodBlockID.empty()&&badBlockID.empty())
      continue;

    //Look through samples for waveforms corresponding to pblk and good/bad cluster blocks
    int r,c,idx,n;
    double adc, tdc;
    vector<string> pblk_title;
    vector<string> gblk_title;
    vector<string> bblk_title;
    for(int m = 0; m < Ndata; m++) {
      r = hcal_row[m];
      c = hcal_col[m];

      if(r>= hcalrows || c >= hcalcols || r<0 || c<0 ){
	cerr << "WARNING: indexing failure. Check hcal row/col indices" << endl;
	continue;
      }

      bool pblk_match = r==pblkrow && c==pblkcol;

      bool gblk_match=false;
      for( size_t g=0; g<goodBlockID.size(); ++g ){
	int blkr = goodBlockID[g]/hcalcols;
	int blkc = goodBlockID[g]%hcalcols;
	if( blkr==r && blkc==c ){
	  gblk_match=true;
	  break;
	}
      }

      bool bblk_match=false;
      for( size_t g=0; g<badBlockID.size(); ++g ){
	int blkr = badBlockID[g]/hcalcols;
	int blkc = badBlockID[g]%hcalcols;
	if( blkr==r && blkc==c ){
	  bblk_match=true;
	  break;
	}
      }

      //Don't bother to proceed if no matches exist
      // if( !gblk_match && !bblk_match && !pblk_match )
      // 	continue;

      idx = hcal_sampsidx[m];
      n = hcal_nsamps[m];
      adc = hcal_adc[m];
      tdc = hcal_tdc[m];

      if(tdc>1200 || tdc<-800)
	tdc=-1000;

      TH1D *h1 = new TH1D("h1","h1",maxSampBins,0,maxSampBins);

      for(int s = 0; s < maxSampBins && s < n; s++)
	h1->SetBinContent(s+1,hcal_samps[idx+s]);

      int lmaxima = CountLocalMaximaAboveThreshold(h1);
      //cout << lmaxima << endl;

      if(pblk_match)
	hnmaxima_pblk->Fill(lmaxima);
      if(gblk_match)
	hnmaxima_gblk->Fill(lmaxima);
      if(bblk_match)
	hnmaxima_bblk->Fill(lmaxima);

      h1->Delete();

      if(pblk_match && pblk_counter<maxSHisto-1){
	pblk_title.push_back(Form("ev:%lld r-c:%d-%d (ADC=%0.2f,TDC=%0.2f) lmax:%d",nevent,r,c,adc,tdc,lmaxima));
	pblk_counter++;
	for(int s = 0; s < maxSampBins && s < n; s++)
	  hWaveform_pblk[pblk_counter]->SetBinContent(s+1,hcal_samps[idx+s]);
      }
      if(gblk_match && gblk_counter<maxSHisto-1){
	gblk_title.push_back(Form("ev:%lld r-c:%d-%d (ADC=%0.2f,TDC=%0.2f) lmax:%d",nevent,r,c,adc,tdc,lmaxima));
	gblk_counter++;
	for(int s = 0; s < maxSampBins && s < n; s++)
	  hWaveform_gblk[gblk_counter]->SetBinContent(s+1,hcal_samps[idx+s]);
      }
      if(bblk_match && bblk_counter<maxSHisto-1){
	bblk_title.push_back(Form("ev:%lld r-c:%d-%d (ADC=%0.2f,TDC=%0.2f) lmax:%d",nevent,r,c,adc,tdc,lmaxima));
	bblk_counter++;

	//cout << bblk_counter << endl;

	for(int s = 0; s < maxSampBins && s < n; s++){

	  //cout << s+1 << " " << hcal_samps[idx+s] << endl;

	  hWaveform_bblk[bblk_counter]->SetBinContent(s+1,hcal_samps[idx+s]);


	}
      }


      // //Build relevant histograms
      // for(int s = 0; s < maxSampBins && s < n; s++){
	
      // 	if(pblk_match){
      // 	  hWaveform_pblk[pblk_counter]->SetBinContent(s+1,hcal_samps[idx+s]);
      // 	  pblk_title.push_back(Form("ev:%lld r-c:%d-%d (ADC=%0.2f,TDC=%0.2f)",nevent,r,c,adc,tdc));
      // 	  pblk_counter++;
      // 	}
      // 	if(gblk_match){
      // 	  hWaveform_gblk[gblk_counter]->SetBinContent(s+1,hcal_samps[idx+s]);
      // 	  gblk_title.push_back(Form("ev:%lld r-c:%d-%d (ADC=%0.2f,TDC=%0.2f)",nevent,r,c,adc,tdc));
      // 	  gblk_counter++;
      // 	}
      // 	if(bblk_match){
      // 	  hWaveform_bblk[bblk_counter]->SetBinContent(s+1,hcal_samps[idx+s]);
      // 	  bblk_title.push_back(Form("ev:%lld r-c:%d-%d (ADC=%0.2f,TDC=%0.2f)",nevent,r,c,adc,tdc));
      // 	  bblk_counter++;
      // 	}

      // }

    }//endloop over samples


    //Set titles for waveform histograms with vectors built before
    int phist_last_idx = pblk_counter;
    for( size_t b=0; b<pblk_title.size(); ++b ){
      int phist_idx_mod = pblk_title.size()-b;

      //cout << phist_last_idx-phist_idx_mod << " " << pblk_title[b] << endl;

      hWaveform_pblk[phist_last_idx-phist_idx_mod]->SetTitle(pblk_title[b].c_str());
    }
    

    int ghist_last_idx = gblk_counter;
    for( size_t b=0; b<gblk_title.size(); ++b ){
      int ghist_idx_mod = gblk_title.size()-b;

      //cout << ghist_last_idx-ghist_idx_mod << " " << gblk_title[b] << endl;

      hWaveform_gblk[ghist_last_idx-ghist_idx_mod]->SetTitle(gblk_title[b].c_str());
    }
    

    int bhist_last_idx = bblk_counter;
    for( size_t b=0; b<bblk_title.size(); ++b ){
      int bhist_idx_mod = bblk_title.size()-b;

      //cout << bhist_last_idx-bhist_idx_mod << " " << bblk_title[b] << endl;

      hWaveform_bblk[bhist_last_idx-bhist_idx_mod]->SetTitle(bblk_title[b].c_str());
    }
    
  }//endloop over events
  
  fout->Write();

  cout << endl << "Analysis complete. Output written to " << outputfilename << endl;

}// end main







