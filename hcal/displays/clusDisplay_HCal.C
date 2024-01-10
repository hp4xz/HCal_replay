#include <TH2.h>
#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "TVectorD.h"
#include "TVector3.h"
#include "TMath.h"
#include "hcal.h"

//SBS8 contants
const Double_t E_e = 5.965; 
const Double_t HCal_th = 29.4 * TMath::DegToRad();
const Double_t HCal_d = 11.;

const Double_t HCal_pRes = 0.07; //HCal position resolution in cm
const Int_t kNrows = 24;
const Int_t kNcols = 12;
const Int_t kNrows_clus = 5;
const Int_t kNcols_clus = 5;
const Int_t Tabs = 4;
const Double_t PI = TMath::Pi();
const Double_t HCal_Xi = -2.090; // Distance from beam center to top of HCal in analysis coordinates (m)
const Double_t HCal_Xf = 1.360; // Distance from beam center to bottom of HCal in analysis coordinates (m)
const Double_t HCal_Yi = -0.825; // Distance from beam center to opposite-beam side of HCal in analysis coordinates (m)
const Double_t HCal_Yf = 0.825; // Distance from beam center to beam side of HCal in analysis coordinates (m)
const Double_t HCal_cd = 0.15; // Cell width (X and Y) (m)
const Double_t hcalheight = 0.365; // Height of HCal above beamline
const Double_t M_p = 0.938272; //Mass proton (GeV)
const Double_t M_n = 0.939565; //Mass neutron (GeV)
Int_t Nrows;
Int_t Ncols;

const Int_t DISP_MIN_SAMPLE = 0;
const Int_t DISP_MAX_SAMPLE = 40;

const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);

const Int_t minADC = 0;
const Int_t maxADC = 4000;
const Int_t kCanvSize = 100;

std::string user_input;

Int_t gCurrentEntry = -1;

TChain *T = 0;
Int_t foundModules = 0;
TCanvas *canvas = 0;

TCanvas *subCanv[4];

void clicked_displayEntryButton();
void clicked_displayNextButton();

namespace hcalgui {
  TGMainFrame *main = 0;
  TGHorizontalFrame *frame1 = 0;
  TGTab *fTab;
  TGLayoutHints *fL3;
  TGCompositeFrame *tf;
  TGTextButton *exitButton;
  TGTextButton *displayEntryButton;
  TGTextButton *displayNextButton;
  TGNumberEntry *entryInput;
  TGLabel *ledLabel;

  TRootEmbeddedCanvas *canv[4];

  TGCompositeFrame* AddTabSub(Int_t sub) {
    if( sub==0 ) tf = fTab->AddTab(Form("PClus Disp"));
    if( sub==1 ) tf = fTab->AddTab(Form("PClus WForm"));
    if( sub==2 ) tf = fTab->AddTab(Form("PClus Err"));
    if( sub==3 ) tf = fTab->AddTab(Form("Other Clus"));
    //tf = fTab->AddTab(Form("HCAL %d",sub+1));

    TGCompositeFrame *fF5 = new TGCompositeFrame(tf, (12+1)*kCanvSize,(6+1)*kCanvSize , kHorizontalFrame);
    TGLayoutHints *fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
        kLHintsExpandY, 5, 5, 5, 5);
    TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(Form("hcalSubCanv%d",sub), fF5, 6*kCanvSize,8*kCanvSize);
    canv[sub] = fEc1;
    fF5->AddFrame(fEc1,fL4);
    tf->AddFrame(fF5,fL4);
    return tf;
  }

  void SetupHMGUI() {
    if(!main) {
      main = new TGMainFrame(gClient->GetRoot(), 1000, 900);
      frame1 = new TGHorizontalFrame(main, 150, 20, kFixedWidth);
      ledLabel = new TGLabel(frame1,"LED Bit:    , NClus:      ");
      displayEntryButton = new TGTextButton(frame1,"&Display Entry","clicked_displayEntryButton()");
      entryInput = new TGNumberEntry(frame1,0,5,-1,TGNumberFormat::kNESInteger);
      displayNextButton = new TGTextButton(frame1,"&Next Entry","clicked_displayNextButton()");
      exitButton = new TGTextButton(frame1, "&Exit", 
          "gApplication->Terminate(0)");
      TGLayoutHints *frame1LH = new TGLayoutHints(kLHintsTop|kLHintsLeft|
          kLHintsExpandX,2,2,2,2);
      frame1->AddFrame(ledLabel,frame1LH);
      frame1->AddFrame(displayEntryButton,frame1LH);
      frame1->AddFrame(entryInput,frame1LH);
      frame1->AddFrame(displayNextButton,frame1LH);
      frame1->AddFrame(exitButton,frame1LH);
      frame1->Resize(800, displayNextButton->GetDefaultHeight());
      main->AddFrame(frame1, new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1));

      // Create the tab widget
      fTab = new TGTab(main, 300, 300);
      fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

      // Create Tabs
      for(Int_t i = 0; i < 4; i++) {
        tf = AddTabSub(i);
      }
      main->AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1));
      main->MapSubwindows();
      main->Resize();   // resize to default size
      main->MapWindow();

      subCanv[0] = canv[0]->GetCanvas();
      subCanv[0]->Divide(2,0,0.001,0.001);

      for( Int_t i=1; i<3; i++ ){
	subCanv[i] = canv[i]->GetCanvas();
	subCanv[i]->Divide(kNcols_clus,kNrows_clus,0.001,0.001);
      }

      subCanv[3] = canv[3]->GetCanvas();
      subCanv[3]->Divide(2,0,0.001,0.001);
      
    }
  }
}


Double_t nhit = 0;
TH1F *histos[kNrows][kNcols];
TH1F *nullBlock;
TH2D *heatMapHisto;
TH2D *clustHisto;
TH2D *clust2Histo;
TH2D *clust3Histo;
Double_t proj_check[2];

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
		     TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

bool is_number(const std::string& mystring)
{
  std::string::const_iterator it = mystring.begin();
  while (it != mystring.end() && std::isdigit(*it)) ++it;
  return !mystring.empty() && it == mystring.end();
}

void displayEvent(Int_t entry = -1)
{
    
    // Check event increment and increment
    if( entry == -1 ) {
      gCurrentEntry++;
    } else {
      gCurrentEntry = entry;
    }

    if( gCurrentEntry < 0 ) {
      gCurrentEntry = 0;
    }

    // Get the event from the TTree
    T->GetEntry( gCurrentEntry );
    std::cout << "Displaying event " << gCurrentEntry << std::endl;
    hcalgui::ledLabel->SetText(TString::Format("LED Bit: %02d, NClus: %5d",Int_t(hcalt::ledbit),Int_t(hcalt::nclus)));
  
    int r,c,idx,n,sub;
    int nblk = hcalt::nblk;
    int cblkid[nblk];
    int cid = hcalt::cid[0];

    // Clear old Displays
    for( Int_t b=1; b<=kNcols_clus*kNrows_clus; b++){
      subCanv[1]->cd(b);
      gStyle->SetOptStat(0);
      nullBlock->Draw();
      //gPad->DrawFrame();
      gPad->SetFillColor(1);
      gPad->Update();

      subCanv[2]->cd(b);
      gStyle->SetOptStat(0);
      nullBlock->Draw();
      //gPad->DrawFrame();
      gPad->SetFillColor(1);
      gPad->Update();
    }

    // Clear old histograms, just in case modules are not in the tree
    for(r = 0; r < kNrows; r++) {
      for(c = 0; c < kNcols; c++) {
	histos[r][c]->Reset("ICES M");
      }
    }

    // Clear old histogram  
    heatMapHisto->Reset("ICES");
    clustHisto->Reset("ICES");
    clust2Histo->Reset("ICES");
    clust3Histo->Reset("ICES");

    Float_t peak[kNrows][kNcols];
    Double_t adc[kNrows][kNcols];
    Double_t adc_p[kNrows][kNcols];
    Double_t tdc[kNrows][kNcols];
    Double_t tdcmult[kNrows][kNcols];
    Double_t pBlockPeak;
    Int_t errIdx = 1;
    for(r  = 0; r < kNrows; r++) {
      for(c  = 0; c < kNcols; c++) {
	peak[r][c] = 0.0;
	adc[r][c] = 0.0;
	adc_p[r][c] = 0.0;
	tdc[r][c] = 0.0;
	tdcmult[r][c] = 0.0;
      }
    }

    for(Int_t m = 0; m < hcalt::ndata; m++) {
      r = hcalt::row[m];
      c = hcalt::col[m];
      Int_t chan = kNcols*r+c;

      if(r < 0 || c < 0) {
	std::cerr << "Why is row negative? Or col?" << std::endl;
	continue;
      }
      if(r>= kNrows || c >= kNcols)
	continue;
      idx = hcalt::samps_idx[m];
      n = hcalt::nsamps[m];
      adc[r][c] = hcalt::a[m];
      adc_p[r][c] = hcalt::a_p[m];

      heatMapHisto->SetBinContent( c+1, r+1, adc_p[r][c] );
      
      tdc[r][c] = hcalt::tdc[m];
      tdcmult[r][c] = hcalt::tdc_mult[m];
      bool saturated = hcalt::a_amp_p[m]>3500;
      bool mansat = false;
      bool negped = adc_p[r][c]<-5;
      bool skip = true;
      for( int b = 0; b < nblk; b++ ){
	if( chan==(hcalt::cblkid[b]-1) ) skip=false;
	//cout << chan << " " << hcalt::cblkid[b] << endl;
      }
      if( !saturated && !negped && skip==true ) continue;
      
      bool displayed = false;
      for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
	displayed = true;
	histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
	if(peak[r][c]<hcalt::samps[idx+s])
	  peak[r][c]=hcalt::samps[idx+s];
	if(peak[r][c]>2) {
	  mansat = true;
	}
      }
      histos[r][c]->SetLineColor(kBlue+1);
      if( saturated==true || negped==true ){
	cout << "Bad block found at r:" << r+1 << " c:" << c+1;
	if( saturated==true ) cout << " -saturated";
	if( negped==true ) cout << " -negative pedestal";
	cout << endl;
	subCanv[2]->cd(errIdx);
	subCanv[2]->SetGrid();
	histos[r][c]->SetTitle(TString::Format("%d-%d (TDC=%g,TMu=%g)",r+1,c+1,tdc[r][c],tdcmult[r][c]));
        histos[r][c]->SetLineColor(kRed+1);
	histos[r][c]->Draw();
	gPad->SetFillColor(18);
	gPad->Update();
	errIdx++;
      }
      if(!displayed) {
	std::cerr << "Skipping empty module: " << m << std::endl;
	for(Int_t s = 0;  s < DISP_FADC_SAMPLES; s++) {
	  histos[r][c]->SetBinContent(s+1,-404);
	}
      }
    }

    // Declare signal amplitude, adc, and tdc arrays for this event
    Int_t zoomCenter = (int)kNrows_clus*kNcols_clus/2+1;
    Int_t clsOffset = zoomCenter - hcalt::cblkid[0];
    Int_t cblkrow_prim = (int)((hcalt::cblkid[0]-1)/kNcols);
    Int_t cblkcol_prim = ((int)hcalt::cblkid[0]-1)%kNcols;
    

    cout << "Number of blocks in primary cluster: " << hcalt::nblk << endl;
    for( int b = 0; b < nblk; b++ ){
      Int_t r = (int)((hcalt::cblkid[b]-1)/kNcols);
      Int_t c = ((int)hcalt::cblkid[b]-1)%kNcols;
      Int_t idx = (int)(hcalt::cblkid[b]-1);
      cout << "                ..at r:" << r << " c:" << c << " with ADCamp_p=" << peak[r][c]*1000 << "(mV) E=" << hcalt::ceblk[b]*1000 << "(MeV) TDC=" << hcalt::tdc[idx] << "(ns) TDCclusblk=" << hcalt::cblktdc[b] << endl;
    }
    // Check the cluster and set array values
    for( int el = 0; el < kNrows*kNcols; el++ ){

      r = el/kNcols;
      c = el%kNcols;
      
      //Populate no block background
      clustHisto->SetBinContent( c+1, r+1, 0.001 );
      clust2Histo->SetBinContent( c+1, r+1, 0.001 );
      clust3Histo->SetBinContent( c+1, r+1, 0.001 );
      
      if( nblk==0 ) clustHisto->SetBinContent( c+1, r+1, -1 );
     
      for( int b = 0; b < nblk; b++ ){
	Int_t cblkid = hcalt::cblkid[b];
	Int_t cblkrow = (int)((hcalt::cblkid[b]-1)/kNcols);
	Int_t cblkcol = ((int)hcalt::cblkid[b]-1)%kNcols;
	Double_t ceblk = hcalt::ceblk[b];
	if( b==0 ){
	  subCanv[1]->cd(cblkid+clsOffset);
	  subCanv[1]->SetGrid();
	}else if( el == cblkid-1 ){
	  Int_t b2relrowcorr = (cblkrow_prim-cblkrow)*kNcols_clus;
	  Int_t b2relcolcorr = cblkcol_prim-cblkcol;
	  subCanv[1]->cd(zoomCenter+b2relrowcorr-b2relcolcorr);
	  subCanv[1]->SetGrid();
	}
	if( el == cblkid-1 ) {
	  clustHisto->SetBinContent( c+1, r+1, ceblk*1000 );
	  histos[r][c]->SetTitle(TString::Format("%d-%d (TDC=%g,TMu=%g)",r+1,c+1,tdc[r][c],tdcmult[r][c]));
	  histos[r][c]->Draw();
	  if( tdcmult[r][c]>1 ){
	    cout << "High TDC Multiplicity: " << tdcmult[r][c] << endl;
	    gPad->SetFillColor(2);
	  }else{
	    gPad->SetFillColor(30);
	  }	    
	}
	gPad->Update();
      }
	
    }

    Int_t cl2r = (int)((hcalt::cid[1])/kNcols);
    Int_t cl3r = (int)((hcalt::cid[2])/kNcols);
    Int_t cl2c = ((int)hcalt::cid[1])%kNcols;
    Int_t cl3c = ((int)hcalt::cid[2])%kNcols;
    
    //cout << hcalt::cid[1] << "=" << cl2r << ":" << cl2c << " " << hcalt::cid[2] << "=" << cl3r << ":" << cl3c << endl;

    if( hcalt::nclus == 3 ){
      clust2Histo->SetBinContent(cl2c,cl2r,10.0);
      clust3Histo->SetBinContent(cl3c,cl3r,10.0);
    }else if( hcalt::nclus == 2 ){
      clust2Histo->SetBinContent(cl2c,cl2r,10.0);
    }
    /*
    cout << hcalt::cx[0] << ":" << hcalt::cblkx[0] << endl;
    cout << hcalt::cy[0] << ":" << hcalt::cblky[0] << endl;
    cout << hcalt::x[0]*(1/HCal_cd)+(-HCal_Xi/HCal_cd)+0.5 << endl;
    cout << hcalt::y[0]*(1/HCal_cd)+(HCal_Yf/HCal_cd)+0.5 << endl;
    */
    //Calculate the projected hadron position from e' reconstruction
    Double_t ephi = atan2( hcalt::BBtr_py[0], hcalt::BBtr_px[0] );
    Double_t nu = E_e - hcalt::BBtr_p[0];
    Double_t pp = sqrt( pow(nu,2)+2.*M_p*nu );
    TVector3 vertex( 0, 0, hcalt::BBtr_vz[0] ); // z location of vertex in hall coordinates
    Double_t phinucleon = ephi + PI; //assume coplanarity
    Double_t thetanucleon = acos( (E_e - hcalt::BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics
    TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );
    TVector3 HCAL_zaxis( sin(-HCal_th), 0, cos(-HCal_th) );
    TVector3 HCAL_xaxis( 0, -1, 0 );
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
    TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;
    Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );
    TVector3 HCAL_intersect = vertex + sintersect * pNhat;
    Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
    Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
    
    //cout << "Raw projected x:y " << xexpect_HCAL << ":" << yexpect_HCAL << endl;
    //cout << "Raw cluster x:y " << hcalt::x[0] << ":" << hcalt::y[0] << endl;
    
    Double_t pRes_conv = HCal_pRes/HCal_cd;

    Double_t exposX = xexpect_HCAL*(1/HCal_cd)+(-HCal_Xi/HCal_cd)+0.5;
    Double_t exposY = yexpect_HCAL*(1/HCal_cd)+(HCal_Yf/HCal_cd)+0.5;
    TEllipse *expos;
    if( xexpect_HCAL>HCal_Xi &&
	xexpect_HCAL<HCal_Xf &&
	yexpect_HCAL>HCal_Yi &&
	yexpect_HCAL<HCal_Yf ){
      expos = new TEllipse(exposY,exposX,pRes_conv,pRes_conv);
      //cout << "Projected ON" << endl;
    }else{
      expos = new TEllipse(0,0,0.,0.);
      cout << "Projected OFF HCal" << endl;
    }
    if( proj_check[0]==xexpect_HCAL || proj_check[1]==yexpect_HCAL ){
      expos = new TEllipse(0,0,0.,0.);
      cout << "No e' track information." << endl;
      //cout << "              ..with px=" << hcalt::BBtr_px[0] << endl;
      //cout << "              ..with py=" << hcalt::BBtr_py[0] << endl;
      //cout << "              ..with pz=" << hcalt::BBtr_pz[0] << endl;
    } 
    proj_check[0]=xexpect_HCAL;
    proj_check[1]=yexpect_HCAL;

    expos->SetFillStyle(0);
    expos->SetLineWidth(2);
    expos->SetLineColor(kGreen);

    //Get cluster position and size info from tree
    Double_t clposX = hcalt::x[0]*(1/HCal_cd)+(-HCal_Xi/HCal_cd)+0.5;
    Double_t clposY = hcalt::y[0]*(1/HCal_cd)+(HCal_Yf/HCal_cd)+0.5;
    TEllipse *clrad = new TEllipse(clposY,clposX,0.3*(1/HCal_cd),0.3*(1/HCal_cd));
    clrad->SetFillStyle(0);
    clrad->SetLineWidth(2);
    clrad->SetLineColor(kWhite);
    TEllipse *clpos = new TEllipse(clposY,clposX,0.05,0.05);
    clpos->SetFillStyle(1001);
    clpos->SetFillColor(0);
    clpos->SetLineWidth(3);
    clpos->SetLineColor(kRed);

    gStyle->SetPalette(53);
    //cout << "cluster x:y " << clposX << ":" << clposY << endl;
    //cout << "projected x:y " << exposX << ":" << exposY << endl;

    subCanv[0]->cd(1); //Heatmap 2d histogram
    subCanv[0]->SetGridx();
    subCanv[0]->SetGridy();
    heatMapHisto->SetStats(0);
    heatMapHisto->GetYaxis()->SetNdivisions(kNrows);
    heatMapHisto->GetXaxis()->SetNdivisions(kNcols);
    heatMapHisto->SetTitle("HCal Pedestal Subtracted ADC (sRAU)");
    heatMapHisto->Draw("colz");
    gPad->Update();

    subCanv[0]->cd(2); //Cluster 2d histogram
    subCanv[0]->SetGridx();
    subCanv[0]->SetGridy();
    clustHisto->SetStats(0);
    clustHisto->GetYaxis()->SetNdivisions(kNrows);
    clustHisto->GetXaxis()->SetNdivisions(kNcols);
    clustHisto->SetTitle("HCal Cluster Elements");
    clustHisto->Draw("colz");
    clpos->Draw("same");
    clrad->Draw("same");
    expos->Draw("same");
    gPad->Update();

    subCanv[3]->cd(1); //Cluster 2d histogram
    subCanv[3]->SetGrid();
    clust2Histo->SetStats(0);
    clust2Histo->SetTitle("HCal Secondary Cluster Elements");
    clust2Histo->Draw("colz");
    gPad->Update();

    subCanv[3]->cd(2); //Cluster 2d histogram
    subCanv[3]->SetGrid();
    clust3Histo->SetStats(0);
    clust3Histo->SetTitle("HCal Tertiary Cluster Elements");
    clust3Histo->Draw("colz");
    gPad->Update();
    
}




void clicked_displayNextButton()
{
  //if(gCurrentEntry>gMaxEntries);
  hcalgui::entryInput->SetIntNumber(++gCurrentEntry);
  displayEvent(gCurrentEntry);
}

void clicked_displayEntryButton()
{
  gCurrentEntry = hcalgui::entryInput->GetIntNumber();
  displayEvent(gCurrentEntry);
}

//Hardcode set for SBS8
Int_t clusDisplay_HCal(Int_t run = 13559, Int_t event = -1){
  
  // Initialize function with user commands
  /*
  cout << "Enter run number for analysis." << endl;
  cin >> run;
  */

  Nrows=1;
  Ncols=1;
  heatMapHisto = new TH2D( "Hits by Module", "", kNcols, 0., kNcols, kNrows, 0., kNrows );
  clustHisto = new TH2D( "Cluster Elements", "", kNcols, 0., kNcols, kNrows, 0., kNrows );
  clust2Histo = new TH2D( "Secondary Cluster Elements", "", kNcols, 0., kNcols, kNrows, 0., kNrows );
  clust3Histo = new TH2D( "Tertiary Cluster Elements", "", kNcols, 0., kNcols, kNrows, 0., kNrows );


  hcalgui::SetupHMGUI();
  gStyle->SetLabelSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.08);  

  if(!T) { 
    T = new TChain("T");
    //T->Add(TString::Format("%s/e1209019_%d*.root",getenv("DATA_DIR"),run));
    T->Add(TString::Format("/volatile/halla/sbs/seeds/rootfiles/hcal_general_%d*",run));
    //T->Add(TString::Format("/volatile/halla/sbs/seeds/rootfiles/hcal_gmn_fullreplay_%d*",run));
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.a_p",hcalt::a_p);
    T->SetBranchAddress("sbs.hcal.a_amp",hcalt::a_amp);
    T->SetBranchAddress("sbs.hcal.a_amp_p",hcalt::a_amp_p);
    T->SetBranchAddress("sbs.hcal.ped",hcalt::ped);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.adcrow",hcalt::row);
    T->SetBranchAddress("sbs.hcal.adccol",hcalt::col);
    T->SetBranchAddress("sbs.hcal.tdc_mult",hcalt::tdc_mult);
    T->SetBranchAddress("sbs.hcal.idblk",hcalt::idblk);
    T->SetBranchAddress("sbs.hcal.x",hcalt::x);
    T->SetBranchAddress("sbs.hcal.y",hcalt::y);

    // Reference TDC branches
    //T->SetBranchAddress("sbs.hcal.Ref.tdc",hcalt::tdc);
    //T->SetBranchAddress("sbs.hcal.Ref.tdc",hcalt::tdc);
    //T->SetBranchAddress("sbs.hcal.Ref.tdc",hcalt::tdc);
    //T->SetBranchAddress("sbs.hcal.Ref.tdcelemID",hcalt::tdc);
    //T->SetBranchAddress("sbs.hcal.Ref.tdc_mult",hcalt::tdc);


    // Add track branches
    T->SetBranchAddress("bb.tr.p",hcalt::BBtr_p);
    T->SetBranchAddress("bb.tr.px",hcalt::BBtr_px);
    T->SetBranchAddress("bb.tr.py",hcalt::BBtr_py);
    T->SetBranchAddress("bb.tr.pz",hcalt::BBtr_pz);
    T->SetBranchAddress("bb.tr.vz",hcalt::BBtr_vz);

    // Add clustering branches
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid);
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk);
    T->SetBranchAddress("sbs.hcal.clus.nblk",hcalt::cnblk);
    T->SetBranchAddress("sbs.hcal.nclus",&hcalt::nclus);
    T->SetBranchAddress("sbs.hcal.nblk",&hcalt::nblk);
    T->SetBranchAddress("sbs.hcal.clus_blk.id",hcalt::cblkid);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk);
    T->SetBranchAddress("sbs.hcal.clus.x",hcalt::cx);
    T->SetBranchAddress("sbs.hcal.clus.y",hcalt::cy);
    T->SetBranchAddress("sbs.hcal.clus_blk.x",hcalt::cblkx);
    T->SetBranchAddress("sbs.hcal.clus_blk.y",hcalt::cblky);
    T->SetBranchAddress("sbs.hcal.clus_blk.tdctime",hcalt::cblktdc);
    //T->SetBranchAddress("sbs.hcal.tdctimeblk",hcalt::);
    //T->SetBranchAddress("sbs.hcal.atimeblk",hcalt::);

    T->SetBranchStatus("Ndata.sbs.hcal.adcrow",1);
    T->SetBranchStatus("Ndata.sbs.hcal.Ref.tdcelemID",1);
    T->SetBranchAddress("Ndata.sbs.hcal.adcrow",&hcalt::ndata);
    std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
    for(Int_t r = 0; r < kNrows; r++) {
      for(Int_t c = 0; c < kNcols; c++) {
        histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
      }
    }
    nullBlock = new TH1F("nullBlock","nullblock",DISP_FADC_SAMPLES,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  }

  if(T->GetEntries()<=0){
    cout << "Error: no data for this run available." << endl;
    //return 0;
  }
  return 0;
  /*
  gCurrentEntry = event;
  while( user_input != "q" ) {
    if(is_number(user_input)) {
      gCurrentEntry = std::stoi(user_input);
    } else {
      gCurrentEntry++;
    }
    displayEvent(gCurrentEntry);
    //std::cout << "Display options: <enter> == next event, or q to stop." << std::endl;
    getline(std::cin,user_input);
  }
  
  return 0;
  */
}

