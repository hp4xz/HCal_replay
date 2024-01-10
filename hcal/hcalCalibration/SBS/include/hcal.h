#ifndef HCAL_H
#define HCAL_H

#include <cmath>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"

const Int_t NTARGET = 3;
const Int_t NSTAMP = 4;
const Int_t NTSTAMP = 10;
const Int_t GMNNCONFIG = 6;

namespace hcal{
  //Paths
  //GMn, pass 1 database
  std::string gmn_p1_db_path = "/work/halla/sbs/seeds/alt_sbsreplay/SBS-replay/DB";
  //GMn, pass 2 database
  std::string gmn_p2_db_path = "/w/halla-scshelf2102/sbs/seeds/sbsoffline/SBS-replay/DB/";

  //HCal Detector
  const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event
  const Int_t maxClus = 35;  //larger than max per tree
  const Int_t maxBlk = 25;  //larger than max per geometry
  //HCal - Note that measurement of vertical length is 381.6cm (imprecise, by hand), indicating that the MC figures are correct
  const Int_t maxHCalChan = 288; // Total HCal channels
  const Int_t maxHCalRows = 24; // Total HCal rows
  const Int_t maxHCalCols = 12; // Total HCal cols
  const Int_t maxHCalClus = 100; // Maximum HCal clusters
  const Int_t maxHCalBlk = 25; // Total HCal cols
  const Int_t maxClusters = 10; // Total HCal clusters with information saved
  const Double_t HCalvoff = -0.2897; // Height of HCal above beamline in m
  const Double_t HCalvoff_p2 = 0.0; // Post pass 2, DB corrected, no offset necessary.
  const Double_t HCalblk_l_h = 0.15494; // Horizontal length of all HCAL blocks in m from MC database
  const Double_t HCalblk_l_v = 0.15875; // Vertical length of all HCAL blocks in m from MC database
  /* const Double_t posHCalXi = -2.3531; // Distance from beam center to top of HCal in m from MC database */
  /* const Double_t posHCalXf = 1.45309; // Distance from beam center to bottom of HCal in m from MC database */
  const Double_t posHCalXi = -2.655; // Distance from beam center to top of HCal in m from MC database
  const Double_t posHCalXf = 1.155; // Distance from beam center to bottom of HCal in m from MC database
  /* const Double_t posHCalYi = -0.93155; // Distance from beam center to opposite-beam side of HCal in m from MC database */
  /* const Double_t posHCalYf = 0.93155; // Distance from beam center to beam side of HCal in m from MC database */
  const Double_t posHCalYi = -0.92964; // Distance from beam center to opposite-beam side of HCal in m from MC database
  const Double_t posHCalYf = 0.92964; // Distance from beam center to beam side of HCal in m from MC database
  const Double_t HCalSampFrac = 0.077;  //Re-evaluated with MC GEn settings using second to outermost shower column for kin2
  //Trigger TDC and timestamps
  const Int_t maxTDCTrigChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
  const Double_t tdiffwidecut = 50; // Set at 50ns from nominal 510ns (passed by user) from GMn
  const Int_t gNstamp = NTSTAMP;
  const string gStamp[NSTAMP] = { "gain", "tdcoff", "tdccalib", "adctoff" }; //Ordered list timestamp types

  //BBCal
  const Int_t maxBBCalShChan = 189; // Total BBCal Shower Channels
  const Int_t maxBBCalShRows = 27;
  const Int_t maxBBCalShCols = 7;
  const Int_t maxBBCalPSChan = 52; // Total BBCal Preshower Channels
  const Int_t maxBBCalPSRows = 26;
  const Int_t maxBBCalPSCols = 2;

  //Beamline and target
  const Int_t chargeConvert = 3318; // See D.Flay Doc DB sbs.jlab.org/DocDB/0001/000164/002/dflay_bcm-ana-update_02-21-22.pdf p.8
  const Int_t clockActual = 103700; // Needed to convert the 104kHz clock to the actual counting rate
  const Int_t gMaxfset = 5; //Reasonable maximum on total SBS field settings per configuration
  const Int_t gNtar = 3; //Total number of possible targets in SBS (LH2, LD2, Helium-3)
  const string gTar[NTARGET] = { "lh2", "ld2", "he3" }; //Ordered list of available targets
  
  //SBS Magnet
  const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
  const Double_t maxSBSfield = 1.26; //Tesla
  const Double_t SBSfield = 1.0; //fraction of max field. TODO: should be variable per run
  const Double_t SBSdist = 2.25; //m
  const Double_t sbsdipolegap = 1.22;  //m
  const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;
  const Int_t gNmag = 20; //Field setting maximum * max targets (reasonable limit)

  //GEMs
  const Double_t GEMpitch = 10*TMath::DegToRad();

  ///////////////
  ///Physics/Math
  const Double_t PI = TMath::Pi();
  const Double_t M_e = 0.00051;
  const Double_t M_p = 0.938272;
  const Double_t M_n = 0.939565;
  const Int_t protonMCidx = 2212; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
  const Int_t neutronMCidx = 2112; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
  const Int_t electronMCidx = 11; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
  const Int_t photonMCidx = 22; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
  const UInt_t second = 1000000; //For conversion to seconds used by reporting time delays (us)

  // Following quantities vary with experiment/configuration or target
  Double_t     ebeam(std::string experiment,Int_t config);        //GeV
  Double_t     bbtheta(std::string experiment,Int_t config);      //deg
  Double_t     bbdist(std::string experiment,Int_t config);       //m
  Double_t     sbstheta(std::string experiment,Int_t config);     //deg
  Double_t     sbsdist(std::string experiment,Int_t config);      //m
  Double_t     hcaltheta(std::string experiment,Int_t config);    //deg
  Double_t     hcaldist(std::string experiment,Int_t config);     //m
  std::string  sbsts(std::string experiment,Int_t config);        //date and time immediately prior to config

  Int_t        tidx(std::string target);                   //index (1:lh2,2:ld2,3:he3)
  Double_t     ltgt(std::string target);                   //m
  Double_t     tarrho(std::string target);                 //g/cc
  Double_t     crho(std::string target);                   //g/cc
  Double_t     cdiam(std::string target);                  //m
  Double_t     cdEdx(std::string target);                  //GeV*m/g
  Double_t     cthick(std::string target);                 //m
  Double_t     uwallthick(std::string target);             //cm
  Double_t     dwallthick(std::string target);             //cm
  Double_t     dEdx(std::string target);                   //GeV*m/g
  Double_t     M_t(std::string target);                    //GeV

}

// a class for SBS config
class SBSconfig {
 public:

  Int_t        GetSBSconf()       const { return fSBSconf; }
  Double_t     GetEbeam()         const { return fEbeam; }
  Double_t     GetBBtheta()       const { return fBBtheta; }
  Double_t     GetBBtheta_rad()   const { return fBBtheta_rad; }
  Double_t     GetBBdist()        const { return fBBdist; }
  Double_t     GetSBStheta()      const { return fSBStheta; }
  Double_t     GetSBStheta_rad()  const { return fSBStheta_rad; }
  Double_t     GetSBSdist()       const { return fSBSdist; }
  Double_t     GetHCALtheta()     const { return fHCALtheta; }
  Double_t     GetHCALtheta_rad() const { return fHCALtheta_rad; }
  Double_t     GetHCALdist()      const { return fHCALdist; }
  std::string  GetSBSTimestamp()  const { return fSBSts; }
  
  // constructor
  SBSconfig(std::string experiment, Int_t config) {
    fSBSconf       = config;
    fEbeam         = hcal::ebeam(experiment,config);
    fBBtheta       = hcal::bbtheta(experiment,config);
    fBBtheta_rad   = hcal::bbtheta(experiment,config)*TMath::DegToRad();
    fBBdist        = hcal::bbdist(experiment,config);
    fSBStheta      = hcal::sbstheta(experiment,config);
    fSBStheta_rad  = hcal::sbstheta(experiment,config)*TMath::DegToRad();
    fSBSdist       = hcal::sbsdist(experiment,config);
    fHCALtheta     = hcal::hcaltheta(experiment,config);
    fHCALtheta_rad = hcal::hcaltheta(experiment,config)*TMath::DegToRad();
    fHCALdist      = hcal::hcaldist(experiment,config);
    fSBSts         = hcal::sbsts(experiment,config);
  }

  // define an ostream operator to print to screen conveniently
  friend ostream& operator <<(ostream &out, const SBSconfig& sbsconf) {
    out  << " -------------------------- "                                             << std::endl
	 << Form(" SBS Config: %d, "                   , sbsconf.fSBSconf)             << std::endl
    	 << Form(" Beam energy: %0.4f (GeV),"          , sbsconf.fEbeam)               << std::endl
    	 << Form(" BigBite angle: %0.1f (deg),"        , sbsconf.fBBtheta)             << std::endl
      	 << Form(" BigBite distance: %0.5f (m),"       , sbsconf.fBBdist)              << std::endl
    	 << Form(" Super BigBite angle: %0.1f (deg),"  , sbsconf.fSBStheta)            << std::endl
      	 << Form(" Super BigBite distance: %0.2f (m)," , sbsconf.fSBSdist)             << std::endl
	 << Form(" HCAL angle: %0.1f (deg),"           , sbsconf.fHCALtheta)           << std::endl
    	 << Form(" HCAL distance: %0.1f (m)"           , sbsconf.fHCALdist)            << std::endl
    	 << Form(" SBS config timestamp: %s"           , sbsconf.fSBSts.c_str())       << std::endl
	 << " -------------------------- "                        << std::endl         << std::endl;
    return out;
  }

 private:
  Int_t        fSBSconf;             // SBS configuration number
  Double_t     fEbeam;               // beam energy (better to get this from tree) (GeV)
  Double_t     fBBtheta;             // BigBite magnet angle (deg)
  Double_t     fBBtheta_rad;         // BigBite magnet angle (rad)
  Double_t     fBBdist;              // BigBite magnet distance from target (m)
  Double_t     fSBStheta;            // Super BigBite magnet angle (deg)
  Double_t     fSBStheta_rad;        // Super BigBite magnet angle (rad)
  Double_t     fSBSdist;             // Super BigBite magnet distance from target (m)
  Double_t     fHCALtheta;           // HCAL angle (deg)
  Double_t     fHCALtheta_rad;       // HCAL angle (rad)
  Double_t     fHCALdist;            // HCAL distance from target (m)
  std::string  fSBSts;               // Date and time immediately prior to config in timestamp format
};

// a class for SBS target
class SBStarget {
 public:

  std::string  GetTarget()              const { return ftarget; }
  Int_t        GetTargIndex()           const { return ftidx; }
  Double_t     GetTargLength()          const { return fltgt; }
  Double_t     GetTargRho()             const { return ftarrho; }
  Double_t     GetCellRho()             const { return fcrho; }
  Double_t     GetCellDiam()            const { return fcdiam; }
  Double_t     GetCelldEdx()            const { return fcdEdx; }
  Double_t     GetCellThick()           const { return fcthick; }
  Double_t     GetUpstreamWallThick()   const { return fuwallthick; }
  Double_t     GetDownstreamWallThick() const { return fdwallthick; }
  Double_t     GetTargdEdx()            const { return fdEdx; }
  Double_t     GetAvgMass()             const { return fM_t; }
  
  // constructor
  SBStarget(std::string target) {
    ftarget        = target;
    ftidx          = hcal::tidx(target);
    fltgt          = hcal::ltgt(target);
    ftarrho        = hcal::tarrho(target);
    fcrho          = hcal::crho(target);
    fcdiam         = hcal::cdiam(target);
    fcdEdx         = hcal::cdEdx(target);
    fcthick        = hcal::cthick(target);
    fuwallthick    = hcal::uwallthick(target);
    fdwallthick    = hcal::dwallthick(target);
    fdEdx          = hcal::dEdx(target);
    fM_t           = hcal::M_t(target);
  }

  // define an ostream operator to print to screen conveniently
  friend ostream& operator <<(ostream &out, const SBStarget& sbstarg) {
    out  << " -------------------------- "                                             << std::endl
	 << Form(" Target: %s, "                       , sbstarg.ftarget.c_str())      << std::endl
	 << Form(" Target index: %d, "                 , sbstarg.ftidx)                << std::endl
	 << Form(" Length of Target: %0.2f (m), "      , sbstarg.fltgt)                << std::endl
    	 << Form(" Target Rho: %0.3f (g/cc),"          , sbstarg.ftarrho)              << std::endl
    	 << Form(" Cell Rho: %0.2f (m),"               , sbstarg.fcrho)                << std::endl
    	 << Form(" Cell Diameter: %0.3f (m),"          , sbstarg.fcdiam)               << std::endl
      	 << Form(" Cell dEdx: %0.4f (GeV*m/g),"        , sbstarg.fcdEdx)               << std::endl
    	 << Form(" Cell Thickness: %0.3f (m),"         , sbstarg.fcthick)              << std::endl
      	 << Form(" Upstream Wall Thickness: %0.4f (m),", sbstarg.fuwallthick)          << std::endl
	 << Form(" Downst Wall Thickness: %0.4f (m),"  , sbstarg.fdwallthick)          << std::endl
    	 << Form(" Target dEdx: %0.4f (GeV*m/g)"       , sbstarg.fdEdx)                << std::endl
    	 << Form(" Average Target Mass: %0.3f (GeV)"   , sbstarg.fM_t)                 << std::endl
	 << " -------------------------- "                        << std::endl         << std::endl;
    return out;
  }

 private:
  std::string  ftarget;                // Target
  Int_t        ftidx;                  // Target index (1:lh2,2:ld2,3:he3)
  Double_t     fltgt;                  // Length of the target (m)
  Double_t     ftarrho;                // Target density (g/cc), He3 is a guess
  Double_t     fcrho;                  // Density of {Al,Al,glass}
  Double_t     fcdiam;                 // Cell diameter (m), He3 is a guess
  Double_t     fcdEdx;                 // Cell wall collisional stopping power (GeV*m/g)
  Double_t     fcthick;                // Target cell thickness (m), He3 is a guess
  Double_t     fuwallthick;            // Upstream wall thickness (cm)
  Double_t     fdwallthick;            // Downstream wall thickness (cm)
  Double_t     fdEdx;                  // Target collisional stopping power (GeV*m/g)
  Double_t     fM_t;                   // Average Target Mass (GeV)
};

#endif // HCAL_H
