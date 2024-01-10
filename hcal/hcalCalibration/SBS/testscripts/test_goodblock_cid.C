//seeds 10.1.23 Test output of new cluser index on goodblocks in calorimeter class, SBS-offline

#include <iostream>
#include <sstream>
#include <fstream>
#include <TMath.h>
#include <cstdio>
#include <cstdlib>
#include "TLatex.h"

const int maxGB = 1000;
const int maxPB = 100;

void test_goodblock_cid(int runnum=12049,int events=15000){//main
  
  TChain *C = new TChain("T");
  
  string rootfilepath = Form("/lustre19/expphy/volatile/halla/sbs/seeds/e1209019_replayed_%d_stream0_seg0_0_firstevent0_nevent%d.root",runnum,events);

  C->Add(rootfilepath.c_str());

  Long64_t Nevents = C->GetEntries();
	 
  cout << "Total entries from " << rootfilepath << ": " << Nevents << endl;
	 
  double hcalGB_id[maxGB], hcalGB_cid[maxGB];
  int NhcalGB_id;

  double hcalB_id[maxPB];
  double nblk, index, nclus;

  C ->SetBranchStatus("*",0);

  C->SetBranchStatus("sbs.hcal.goodblock.id",1);
  C->SetBranchStatus("sbs.hcal.goodblock.cid",1);
  C->SetBranchStatus("Ndata.sbs.hcal.goodblock.id",1);
  C->SetBranchStatus("sbs.hcal.nblk",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.id",1);
  C->SetBranchStatus("sbs.hcal.index",1);
  C->SetBranchStatus("sbs.hcal.nclus",1);

  C->SetBranchAddress("sbs.hcal.goodblock.id", hcalGB_id);
  C->SetBranchAddress("sbs.hcal.goodblock.cid", hcalGB_cid);
  C->SetBranchAddress("Ndata.sbs.hcal.goodblock.id", &NhcalGB_id);
  C->SetBranchAddress("sbs.hcal.clus_blk.id", hcalB_id);
  C->SetBranchAddress("sbs.hcal.nblk", &nblk);
  C->SetBranchAddress("sbs.hcal.index", &index);
  C->SetBranchAddress("sbs.hcal.nclus", &nclus);

  //histograms
  TH1D *hblkid_cb = new TH1D("hblkid_cb","Block ID primary cluster (all blocks), clus_blk branch; channel",310,-10,300);
  TH1D *hblkid_gb = new TH1D("hblkid_gb","Block ID cluster id primary (all blocks), goodblock branch cut on sbs.hcal.index; channel",300,-10,290);

  for( Long64_t nevent = 1; nevent <Nevents; nevent++){

    // cout << " Entry = " << nevent << "/" << Nevents << "\r";
    // cout.flush();
    
    C->GetEntry(nevent); 
    
    if(nclus<1)
      continue;

    cout << "Nclus>0 event " << nevent << " best block index " << index << endl;

    //write out the blocks in the primary cluster
    int NB = nblk;
    cout << "Number of blocks in primary cluster: " << NB << endl;
    cout << "All blocks in cluster:" << endl; 
    for(int b=0; b<NB; ++b){
      double bid = hcalB_id[b];
      hblkid_cb->Fill(bid);
      cout << "  block:id = " << b << ":" << bid << endl;
    }
    
    //write out the goodblocks with cluster id of zero
    int NGB = NhcalGB_id;
    int counter = 0;
    cout << "Number of goodblocks: " << NGB << endl;
    // cout << "All blocks whose cluster id is zero:" << endl;
    for(int b=0; b<NGB; ++b){

      //cout << " GBclusid:GBid " << hcalGB_cid[b] << ":" << hcalGB_id[b] << endl;

      double gbid = hcalGB_id[b];
      if(hcalGB_cid[b]==index){
    	cout << "  goodblock:id = " << b << ":" << gbid << endl;
	hblkid_gb->Fill(gbid);
	counter++;
      }
    }
    cout << "Number of goodblocks with cluster id = best block index: " << counter << endl;

    cout << endl << endl;

    if( counter != NB ){
      cout << "ERROR: Goodblock primary index total blocks not equal to clus_blk total blocks" << endl;
      return;
    }
  }

  //Make canvas for direct comparison on highest E cluster
  TCanvas *c1 = new TCanvas("c1","goodblock Primary Index vs clus_blk IDs",1200,500);
  c1->cd();

  hblkid_cb->SetTitle("Block ID, primary cluster (all cluster blocks)");
  hblkid_cb->SetLineColor(kGreen);
  hblkid_cb->SetLineWidth(2);
  //hblkid_cb->SetLineStyle(3);
  hblkid_cb->SetFillColor(kGreen);
  hblkid_cb->SetFillStyle(3005);
  hblkid_cb->Draw();

  hblkid_gb->SetLineColor(kBlue);
  hblkid_gb->SetLineWidth(2);
  hblkid_gb->SetLineStyle(2);
  hblkid_gb->SetFillColor(kBlue);
  hblkid_gb->SetFillStyle(3004);
  hblkid_gb->Draw("same");

  //Add a legend to the canvas
  auto leg = new TLegend(0.1,0.6,0.5,0.9);
  leg->AddEntry( hblkid_cb, "clus_blk, all blocks", "l");
  leg->AddEntry( hblkid_gb, "goodblock, cid==sbs.hcal.index, all blocks", "l");
  leg->Draw();

  c1->Write();

  cout << endl << endl << "Diagnostic completed successfully." << endl;

}// end main







