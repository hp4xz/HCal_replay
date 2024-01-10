//seeds 10.1.23 Test tmax db_sbs.hcal.dat output for comparison with MC elastic, hcal nblk

#include <iostream>
#include <sstream>
#include <fstream>
#include <TMath.h>
#include <cstdio>
#include <cstdlib>
#include "TLatex.h"
#include "../include/sbs.h"

const int maxGB = 1000;
const int maxPB = 100;

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

bool adjacent(int seed, int block){

  int seed_row = seed/hcalcols;
  int seed_col = seed%hcalcols;
  int block_row = block/hcalcols;
  int block_col = block%hcalcols;

  bool rowclose = abs(block_row-seed_row)<=1;
  bool colclose = abs(block_col-seed_col)<=1;

  return rowclose && colclose;

}

void test_tmax(int runnum=12674,int events=50002){//main
 
  TChain *C = new TChain("T");

  //string rootfilepath = Form("/lustre19/expphy/volatile/halla/sbs/seeds/e1209019_replayed_%d_stream0_seg0_0_firstevent0_nevent%d.root",runnum,events);

  //string rootfilepath = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS11/LH2/rootfiles/e1209019_fullreplay_12380*.root";

  //string rootfilepath = "/lustre19/expphy/volatile/halla/sbs/seeds/e1209019_replayed_12916_stream0_seg0_0_firstevent0_nevent10013.root";

  string rootfilepath = Form("/lustre19/expphy/volatile/halla/sbs/seeds/vary_emin/e1209019_replayed_%d_stream0_seg0_0_firstevent0_nevent%d.root",runnum,events);

  //string outputfilename = Form("root_out/tmaxcompOUT_run%d_Nev%d.root",runnum,events);
  string outputfilename = Form("root_out/tmaxcompOUT_control_run%d_Nev%d.root",runnum,events);

  C->Add(rootfilepath.c_str());

  cout << "Added rootfile " << rootfilepath << " to chain." << endl;

  int sbs = -1;
  if(runnum==11587)
    sbs=4;
  if(runnum==12049)
    sbs=7;
  if(runnum==12369||runnum==12674||runnum==12916||runnum==13057)
    sbs=11;
  if(runnum==13242||runnum==13314)
    sbs=14;
  if(runnum==13486)
    sbs=8;
  if(runnum==13662||runnum==13682||runnum==13793)
    sbs=9;


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
  //TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&sbs.hcal.e>0.025&&bb.ps.e+bb.sh.e>1.7";
  //TCut globalcut = "bb.ps.e>0.2&&sbs.hcal.e>0.040";
  //TCut globalcut = "sbs.hcal.e>0.040";


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
  //TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&sbs.hcal.e>0.13&&bb.ps.e+bb.sh.e>3.2";
  

  //TCut globalcut = "bb.tr.n==1";

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
  TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&sbs.hcal.nclus>0&&bb.sh.nclus>0";

  if(sbs==-1){
    cout << "ERROR: No valid sbs runnumber given." << endl;
    return;
  }
  //TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.01&&bb.ps.e+bb.sh.e>1.7";

  //TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.025&&bb.ps.e+bb.sh.e>1.5";

  //TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&sbs.hcal.e>0.130&&bb.ps.e+bb.sh.e>1.5";

  //TCut globalcut = "bb.tr.n==1";
  
  // TEventList *elist = new TEventList("elist","Elastic Event List");
  // C->Draw(">>elist",globalcut);

  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vz[maxTracks];
  //double BBtgt_x[maxTracks], BBtgt_y[maxTracks], BBtgt_th[maxTracks], BBtgt_ph[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, BBsh_atime;
  double HCALx, HCALy, HCALe, ekineW2;
  double hcalGB_id[maxGB], hcalGB_cid[maxGB];
  int NhcalGB_id;

  double hcalB_id[maxPB], cblkatime[maxPB], cblktime[maxPB], cblke[maxPB];
  double nblk, index, nclus;

  C ->SetBranchStatus("*",0);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);
  C->SetBranchStatus("bb.sh.atimeblk",1);
  C->SetBranchStatus("bb.sh.nclus",1);
  C->SetBranchStatus("bb.gem.track.nhits",1);
  C->SetBranchStatus("bb.etot_over_p",1);
  //C->SetBranchStatus("e.kine.W2",1);
  // C->SetBranchStatus("sbs.hcal.goodblock.id",1);
  // C->SetBranchStatus("sbs.hcal.goodblock.cid",1);
  // C->SetBranchStatus("Ndata.sbs.hcal.goodblock.id",1);
  C->SetBranchStatus("sbs.hcal.nblk",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.id",1);
  C->SetBranchStatus("sbs.hcal.index",1);
  C->SetBranchStatus("sbs.hcal.nclus",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.atime", 1);
  C->SetBranchStatus("sbs.hcal.clus_blk.tdctime", 1);
  C->SetBranchStatus("sbs.hcal.clus_blk.e", 1);

  C->SetBranchStatus("bb.tr.n", 1);

  C->SetBranchAddress("sbs.hcal.x", &HCALx);
  C->SetBranchAddress("sbs.hcal.y", &HCALy);
  C->SetBranchAddress("sbs.hcal.e", &HCALe);
  C->SetBranchAddress("bb.tr.n", &BBtr_n);
  C->SetBranchAddress("bb.tr.px", BBtr_px);
  C->SetBranchAddress("bb.tr.py", BBtr_py);
  C->SetBranchAddress("bb.tr.pz", BBtr_pz);
  C->SetBranchAddress("bb.tr.p", BBtr_p);
  C->SetBranchAddress("bb.tr.vz", BBtr_vz);
  C->SetBranchAddress("bb.ps.e", &BBps_e);
  C->SetBranchAddress("bb.ps.x", &BBps_x);
  C->SetBranchAddress("bb.ps.y", &BBps_y);
  C->SetBranchAddress("bb.sh.e", &BBsh_e);
  C->SetBranchAddress("bb.sh.x", &BBsh_x);
  C->SetBranchAddress("bb.sh.y", &BBsh_y);
  C->SetBranchAddress("bb.sh.atimeblk", &BBsh_atime);
  // C->SetBranchAddress("sbs.hcal.goodblock.id", hcalGB_id);
  // C->SetBranchAddress("sbs.hcal.goodblock.cid", hcalGB_cid);
  // C->SetBranchAddress("Ndata.sbs.hcal.goodblock.id", &NhcalGB_id);
  C->SetBranchAddress("sbs.hcal.clus_blk.id", hcalB_id);
  C->SetBranchAddress("sbs.hcal.nblk", &nblk);
  C->SetBranchAddress("sbs.hcal.index", &index);
  C->SetBranchAddress("sbs.hcal.nclus", &nclus);
  C->SetBranchAddress("sbs.hcal.clus_blk.atime", cblkatime); // Array of block ADC times
  C->SetBranchAddress("sbs.hcal.clus_blk.tdctime", cblktime); // Array of block ADC times
  C->SetBranchAddress("sbs.hcal.clus_blk.e", cblke); // Array of block ADC times

  TFile *fout = new TFile(outputfilename.c_str(),"RECREATE");

  TH1D *hW2 = new TH1D("hW2", " ;GeV2  ", 100,0,5);
  TH1D *hQ2 = new TH1D("hQ2", " ;GeV2  ", 150,0,15);
  TH1D *hblkid = new TH1D("hblkid", " ;ID  ", 310,-10,300);
  TH1D *hcloseid = new TH1D("closid", " ;ID  ", 35,-5,30);

  TH2D *hEvblkid = new TH2D("hEvblkid",
			    "Block E vs cluster block id",
			    20,
			    0,
			    20,
			    1000.,
			    0.,
			    2.0);
  
  TH2D *hEvblkid_elascut = new TH2D("hEvblkid_elascut",
			       "Block E vs cluster block id elastic cut",
			       20,
			       0,
			       20,
			       1000.,
			       0.,
			       2.0);

  TH2D *ht_alldiff_v_closeid = new TH2D("ht_alldiff_v_closeid",
			       "adct all channels/runs all cluster block difference vs adjacent block id",
			       11,
			       -1,
			       10,
			       total_bins_short,
			       lower_lim_short,
			       upper_lim_short);

  TH2D *ht_alldiff_v_eid = new TH2D("ht_alldiff_v_eid",
			       "adct all channels/runs all cluster block difference vs descending block e id",
			       11,
			       -1,
			       10,
			       total_bins_short,
			       lower_lim_short,
			       upper_lim_short);

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

  TH1D *ht_all_diff_border = new TH1D("ht_all_border",
			       "adct all channels/runs all cluster block difference, bordering",
			       total_bins_short,
			       lower_lim_short,
			       upper_lim_short);

  TH1D *ht_all_diff_separated = new TH1D("ht_all_separated",
			       "adct all channels/runs all cluster block difference, separated",
			       total_bins_short,
			       lower_lim_short,
			       upper_lim_short);

  TH1D *ht_all_diff_close = new TH1D("ht_all_diff_close",
				     "adct all channels/runs all cluster block difference, adjacent",
				     total_bins_short,
				     lower_lim_short,
				     upper_lim_short);

  TH1D *ht_all_diff_far = new TH1D("ht_all_diff_far",
				   "adct all channels/runs all cluster block difference, not adjacent",
				   total_bins_short,
				   lower_lim_short,
				   upper_lim_short);

  TH1D *ht_all_diff_elascut = new TH1D("ht_all_diff_elascut",
				       "adct all channels/runs all cluster block difference wide elastic cut",
				       total_bins_short,
				       lower_lim_short,
				       upper_lim_short);

  TH1D *ht_all_diff_coincut = new TH1D("ht_all_diff_coincut",
				       "adct all channels/runs all cluster block difference BBCal/HCal coin cut",
				       total_bins_short,
				       lower_lim_short,
				       upper_lim_short);

  TH1D *ht_all_diff_allcut = new TH1D("ht_all_diff_allcut",
				       "adct all channels/runs all cluster block difference BBCal/HCal all cut",
				       total_bins_short,
				       lower_lim_short,
				       upper_lim_short);
  
  TH1D *hnblk = new TH1D("hnblk", " ;N", 20,0,20);
  TH1D *hnblk_elascut = new TH1D("hnblk_elascut", " ;N", 20,0,20);
  TH1D *hnblk_coincut = new TH1D("hnblk_coincut", " ;N", 20,0,20);
  TH1D *hnblk_allcut = new TH1D("hnblk_allcut", " ;N", 20,0,20);

  TH1D *hcoin = new TH1D("hcoin", "BBCal SH atimeblk - HCal atimeblk ;ns", total_bins_short,
			 lower_lim_short,
			 upper_lim_short);

  //set up the global cut formula
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  Long64_t Nevents = C->GetEntries();
  
  //ttree formula markers
  int treenum = 0, currenttreenum = 0;

  cout << Nevents << " in tree." << endl;

  for( Long64_t nevent = 1; nevent <Nevents; nevent++){

    cout << " Entry = " << nevent << "/" << Nevents << "\r";
    cout.flush();
    
    C->GetEntry(nevent); 
    
    //Single-loop globalcut method. Save pass/fail for output tree.
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum; 
      GlobalCut->UpdateFormulaLeaves();
    }
    bool failedglobal = GlobalCut->EvalInstance(0) == 0;

    //elastic kinematics
    double etheta = acos( BBtr_pz[0]/BBtr_p[0]);
    double ephi = atan2(BBtr_py[0],BBtr_px[0]);
    
    TVector3 vertex(0,0,BBtr_vz[0]);
    TLorentzVector Pbeam(0,0,E_e,E_e);
    TLorentzVector kprime(BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0]);
    TLorentzVector Ptarg(0, 0, 0, M_p);

    TLorentzVector q = Pbeam - kprime;
    TLorentzVector PgammaN = Ptarg + q; //should go through and write this out. Momentum of virtual photon

    double pel = E_e/ (1. +E_e/M_p*(1.-cos(etheta)));//momentum of elastically scattered electron 
    double nu = E_e -BBtr_p[0]; //kinetic energy of the elasticlly scattered electron 
    double pp = sqrt(pow(nu,2)+2 *M_p*nu); 
    double phinucleon = ephi + PI; //coplanar 
    double thetanucleon = acos((E_e - BBtr_pz[0])/pp);

    TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );
    
    TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
    TVector3 HCAL_xaxis(0,-1,0);
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();

    TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

    double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis) / (pNhat.Dot( HCAL_zaxis ) );
    
    TVector3 HCAL_intersect = vertex + sintersect * pNhat; 

    double yexpect_HCAL = ( HCAL_intersect - HCAL_origin ).Dot( HCAL_yaxis );
    double xexpect_HCAL = ( HCAL_intersect - HCAL_origin ).Dot( HCAL_xaxis );

    double Q2 = 2*E_e *BBtr_p[0]*( 1-cos(etheta));
    double W2 = pow( M_p,2 )+2*M_p*nu-Q2;
  
    hW2->Fill(W2);
    hQ2->Fill(Q2);
    
    double dx = HCALx - xexpect_HCAL;
    double dy = HCALy - yexpect_HCAL;

    //coin
    double pblkatime = cblkatime[0];
    double pblke = cblke[0];
    double pblktime = cblktime[0];
    double pblkid = hcalB_id[0]-1;
    int pblkrow = pblkid/hcalcols;
    int pblkcol = (int)pblkid%hcalcols;
    double cointime = BBsh_atime - pblkatime;
    bool failedcoin = cointime<-10 || cointime>5;
    bool failedW2 = W2>1.2;

    hcoin->Fill(cointime);

    bool hcalaa = pblkrow!=0 && pblkrow!=23 && pblkcol!=0 && pblkcol!=11;

    if( !hcalaa )
      continue;

    hEvblkid->Fill(0.,pblke);
    if( !failedglobal && !failedW2 )
      hEvblkid_elascut->Fill(0.,pblke);

    //BigBite/SBS Acceptance matching. Only valid for neutrons or zero field
    // bool failedaccmatch = 
    //   yexpect_HCAL > hcal::posHCalYf ||
    //   yexpect_HCAL < hcal::posHCalYi ||
    //   xexpect_HCAL > hcal::posHCalXf ||
    //   xexpect_HCAL < hcal::posHCalXi;

    //Quality checks on tmax
    for( int b=1; b<nblk; b++ ){

      double blkatime = cblkatime[b];
      double blktime = cblktime[b];
      double blkid = hcalB_id[b]-1;
      double blke = cblke[b];
      int blkrow = blkid/hcalcols;
      int blkcol = (int)blkid%hcalcols;

      hblkid->Fill(blkid);
      hEvblkid->Fill(b,blke);
      if( !failedglobal && !failedW2 )
	hEvblkid_elascut->Fill(b,pblke);

      double diff = pblkatime - blkatime;
      double tdiff = pblktime - blktime;

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

      if( close_idx%2==0 )
	ht_all_diff_border->Fill(diff);
      else
	ht_all_diff_separated->Fill(diff);

      hcloseid->Fill(close_idx);

      ht_alldiff_v_closeid->Fill(close_idx,diff);

      //Fill atime diff histos
      ht_all_diff->Fill(diff);
      if(tdiff!=0)
	ht_all_tdiff->Fill(tdiff);
      ht_alldiff_v_eid->Fill(b,diff);

      if(close)
	ht_all_diff_close->Fill(diff);
      else
	ht_all_diff_far->Fill(diff);

      if( !failedglobal && !failedW2 )
	ht_all_diff_elascut->Fill(diff);

      if( !failedcoin )
	ht_all_diff_coincut->Fill(diff);

    }

    //Fill nblk histos
    hnblk->Fill(nblk);
    if( !failedglobal && !failedW2 && !failedcoin )
      hnblk_allcut->Fill(nblk);

    if( !failedglobal && !failedW2 )
      hnblk_elascut->Fill(nblk);

    if( !failedcoin )
      hnblk_coincut->Fill(nblk);

  }

  //cout << "Total elastics = " << elastic;

  fout->Write();

  cout << endl << "Analysis complete. Output written to " << outputfilename << endl;

}// end main







