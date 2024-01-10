#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TString.h"

namespace util {

  //Console
  string getDate();

  //General
  void readParam(std::string const_path,vector<Double_t> &param); //reads parameters from txt file
  Int_t countSets(std::string const_path,std::string type); //Counts the number of calibration sets for plotting
  void readDB(std::string const_path,std::string timestamp,std::string type,Double_t param[]); 
  void readDB(std::string const_path,std::string timestamp,std::string type,Double_t &param); //overload 1
  void readDB(std::string const_path,std::string timestamp,std::string type,vector<Double_t> &param); //overload 2
  void tsCompare(std::string tsA,std::string tsB,std::string &tsLate);

  //Functions to read csv files
  void ReadRunList(std::string runsheet_dir,    // Dir. path containing CSV files with run info
		   std::string experiment,      // SBS experiment {gmn, gen, etc.}
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   Int_t replay_pass,           // replay pass
		   Int_t verbose,               // verbosity
		   vector<calrun> &crun);       // Output: Vector of crun structs

  void ReadCutList(std::string cutsheet_dir,    // Dir. path containing CSV files with cut info
		   std::string experiment,      // experiment {gmn,gen,genrp,gep,etc.}
		   Int_t sbsconfig,             // SBS configuration
		   Int_t replay_pass,           // replay pass
		   std::string target,          // target
		   Int_t field,                 // sbs field in percent
		   Int_t verbose,               // verbosity
		   vector<calcut> &ccut);       // Output: Vector of ccut structs

  void ReadDiagnosticCutList(std::string cutsheet_dir,    // Dir. path containing CSV files with cut info
			     std::string experiment,      // experiment {gmn,gen,genrp,gep,etc.}
			     Int_t sbsconfig,             // SBS configuration
			     std::string target,          // target
			     Int_t field,                 // sbs field in percent
			     Int_t verbose,               // verbosity
			     vector<caldiag> &cdcut);     // Output: Vector of ccut structs
  
  void diagnosticReport(const std::vector<caldiag>& data, // Vector of many caldiag objects by target/sbs-field
			const suppset& supplement,        // Set of expermental objects to augment output
			std::string report_path);         // Path to report file holding canvases

  //Geometry/physics
  void sethcalaxes( Double_t sbstheta_rad,                 //SBS arm theta in radians
		    vector<TVector3> &hcal_axes );         //Pass 3-vector (unit) for modification in HCal coordinates

  void getxyhcalexpect( TVector3 vertex,                   //Unit 3-vector from the z-component of the vertex
			TVector3 pNhat,                    //Unit direction of momentum transfer to recoil nucleon
			TVector3 hcal_origin,              //Location of hcal origin in hall coordinates
			vector<TVector3> hcal_axes,        //Define the hcal axes (unit)
			vector<Double_t> &xyhcalexpect );  //Expected location of elastic recoil nucleon on hcal

  void checkPID( std::string target,                       //Target can sometimes disambiguate p/n peaks
		 Double_t dx0p,                            //Distribution dx location of proton peak
		 Double_t dx0n,                            //Distribution dx location of neutron peak
		 Double_t dxsigp,                          //Distribution dx proton sigma
		 Double_t dxsign,                          //Distribution dx neutron sigma
		 Double_t dx,                              //Event dx
		 Int_t &pid );                             //Modified pid(-1:neither,0:ambiguous,1:proton,2:neutron)
  
  bool Nspotcheck(Double_t dy,                             //data dy
		  Double_t dx,                             //data dx
		  Double_t dy_mean,                        //mean location of nucleon dy
		  Double_t dx_mean,                        //mean location of nucleon dx
		  Double_t dy_sigma,                       //dy sigma of nucleon
		  Double_t dx_sigma,                       //dx sigma of nucleon
		  Double_t rotationAngle=0);               //rotation angle of ellipse (if ever necessary)

  Double_t assignScore( Double_t val_1,                    //Cluster energy
			Double_t val_2,                    //Cluster coin time
			Double_t max1,                     //Max cluster energy for event
			const vector<Double_t> &gfit );     //coin time gaussian profile 

  //Histogram utility
  Int_t get_max_bin_on_range( TH1D *h,                     //1D histogram to find max bin
			      Double_t xlow,               //lower limit on find range
			      Double_t xhigh );            //upper limit on find range

  vector<Double_t> fitGaussianAndGetFineParams(TH1D* hist,           //1D histogram to fit
					       Double_t sig,         //guess at std dev for course fit
					       Double_t low,         //min limit on fit range  
					       Double_t high);       //max limit on fit range

  void sliceHisto( TH2D *h2,                               //2D histogram to slice
		   Int_t Nslices,                          //number of x bins for slicing
		   Double_t fwhm,                          //guess at std dev of gaus fit for preliminary fits
		   Int_t min_ev,                           //minimum events to fit the x-cell
		   vector<Double_t> &cell,                 //vector to store bin x values
		   vector<Double_t> &mean,                 //vector to store mean values of fit
		   vector<Double_t> &err,                  //vector to store std dev values of fit
		   Double_t meanXmin,                      //hard min limit on finding mean value (if known)
		   Double_t meanXmax);                     //hard max limit on finding mean value (if known)

  void sliceHCalIDHisto( TH2D *h2,                               //2D histogram to slice
			 Int_t Nslices,                          //number of x bins for slicing
			 Double_t fwhm,                          //guess at std dev of gaus fit for preliminary fits
			 Int_t min_ev,                           //minimum events to fit the x-cell
			 TCanvas *c1,                            //Canvas to hold fits for top half cells
			 TCanvas *c2,                           //Canvas to hold fits for bottom half cells
			 vector<Double_t> &cell,                 //vector to store bin x values
			 vector<Double_t> &mean,                 //vector to store mean values of fit
			 vector<Double_t> &err,                  //vector to store std dev values of fit
			 Double_t meanXmin,                      //hard min limit on finding mean value (if known)
			 Double_t meanXmax);                           //Canvas to hold fits for bottom half cells

  void sliceHistoFine( TH2D *h2,                           //2D histogram to slice
		       Int_t Nslices,                      //number of x bins for slicing
		       Double_t fwhm,                      //guess at std dev of gaus fit for preliminary fits
		       Int_t min_ev,                       //minimum events to fit the x-cell
		       vector<Double_t> &cell,             //vector to store bin x values
		       vector<Double_t> &mean,             //vector to store mean values of fit
		       vector<Double_t> &err );            //vector to store std dev values of fit

  //Fit functions
  Double_t g_gfit(Double_t *x, Double_t *par);
  Double_t g_gfit_bg(Double_t *x, Double_t *par);
  Double_t g_sgfit(Double_t *x, Double_t *par);
  Double_t g_sgfit_bg(Double_t *x, Double_t *par);
  Double_t g_expofit(Double_t *x, Double_t *par);
  Double_t g_twfit(Double_t *x, Double_t *par);
  Double_t g_scexpofit(Double_t *x, Double_t *par);
  Double_t g_lfit(Double_t *x, Double_t *par);
  Double_t g_p3fit(Double_t *x, Double_t *par);
  Double_t g_p4fit(Double_t *x, Double_t *par);
  Double_t g_p6fit(Double_t *x, Double_t *par);

}

#endif
