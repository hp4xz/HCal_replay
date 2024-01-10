#include "../include/hcal.h"

namespace hcal {

  //Beam energy in GeV
  Double_t ebeam(std::string experiment,Int_t config) {
    if( experiment.compare("gmn")==0 ){

      if(config==1)
	return 1.916;
      else if(config==4)
	//return 3.7278;
	return 3.7393;
      else if(config==7)
	//return 7.906;
	return 7.9308;
      else if(config==11)
	//return 9.91;
	return 9.889;
      else if(config==14)
	//return 5.965;
	return 5.9827;
      else if(config==8)
	//return 5.965;
	return 5.9826;
      else if(config==9)
	//return 4.013;
	return 4.0268;
      else if(config==4363)
	return 6.373;

    }else if( experiment.compare("gen")==0 ){ //This section to be finished as new experiments are performed
      
      if(config==1)
	return 4.291;
      if(config==2)
	return 4.291;
      if(config==3)
	return 6.373;
      if(config==4||config==45)
	return 8.4;

    }else{
      std::cerr << "Error: enter a valid SBS configuration." << std::endl;
      return -1;
    }
    return -1;
  }

  //Angle of the BigBite (electron) arm wrt exit beamline in degrees
  Double_t bbtheta(std::string experiment,Int_t config){
    if( experiment.compare("gmn")==0 ){

      if(config==1)
	return 51.0;
      else if(config==4)
	return 36.0;
      else if(config==7)
	return 40.0;
      else if(config==11)
	return 42.0;
      else if(config==14)
	return 46.5;
      else if(config==8)
	return 26.5;
      else if(config==9)
	return 49.0;
      else if(config==4363)
	return 36.5;

    }else if( experiment.compare("gen")==0 ){ //This section to be finished as new experiments are performed
      
      if(config==2||config==1)
	return 29.5;
      else if(config==3)
	return 36.5;
      else if(config==4||config==45)
	return 35;

    }else{
      std::cerr << "Error: enter a valid SBS configuration." << std::endl;
      return -1;
    }
    return -1;
  }

  //Distance from the target to the BigBite magnet in m
  Double_t bbdist(std::string experiment,Int_t config){
    if( experiment.compare("gmn")==0 ){

      if(config==1)
	return 1.8518;
      else if(config==4)
	return 1.7988;
      else if(config==7)
	return 1.84896;
      else if(config==11)
	return 1.55146;
      else if(config==14)
	return 1.84787;
      else if(config==8)
	return 1.97473;
      else if(config==9)
	return 1.550;
      else if(config==4363)
	return 1.63;

    }else if( experiment.compare("gen")==0 ){ //This section to be finished as new experiments are performed
      
      if(config==2||config==3||config==4||config==1||config==45)
	return 1.63; 

    }else{
      std::cerr << "Error: enter a valid SBS configuration." << std::endl;
      return -1;
    }
    return -1;
  }

  //Angle of the SBS magnet (in hadron arm) wrt exit beamline in degrees
  Double_t sbstheta(std::string experiment,Int_t config){
    if( experiment.compare("gmn")==0 ){

      if(config==1)
	return 33.5;
      else if(config==4)
	return 31.9;
      else if(config==7)
	return 16.1;
      else if(config==11)
	return 13.3;
      else if(config==14)
	return 17.3;
      else if(config==8)
	return 29.9;
      else if(config==9)
	return 22.5;
      else if(config==4363)
	return 22.1;

    }else if( experiment.compare("gen")==0 ){ //This section to be finished as new experiments are performed
      
      if(config==2||config==1)
	return 34.7;
      else if(config==3)
	return 22.1;
      else if(config==4||config==45)
	return 18;

    }else{
      std::cerr << "Error: enter a valid SBS configuration." << std::endl;
      return -1;
    }
    return -1;
  }

  //Distance from target to the SBS magnet in m
  Double_t sbsdist(std::string experiment,Int_t config){
    if( experiment.compare("gmn")==0 ){
      
      if(config==1||config==4||config==7||config==11
	 ||config==14||config==8||config==9)
	return 2.25;
      else if(config==4363)
	return 2.8;

    }else if( experiment.compare("gen")==0 ){ //This section to be finished as new experiments are performed
      
      if(config==45||config==1||config==2||config==3||config==4)
	return 2.8;

    }else{
      std::cerr << "Error: enter a valid SBS configuration." << std::endl;
      return -1;
    }
    return -1;
  }

  //Distance from target to the hadron calorimeter in m
  Double_t hcaldist(std::string experiment,Int_t config){
    if( experiment.compare("gmn")==0 ){
      
      if(config==1)
	return 13.5;
      else if(config==4||config==8||config==9)
	return 11.0;
      else if(config==7||config==14)
	return 14.0;
      else if(config==11)
	return 14.5;
      else if(config==4363)
	return 17.0;

    }else if( experiment.compare("gen")==0 ){ //This section to be finished as new experiments are performed
     
      if(config==45||config==1||config==2||config==3||config==4)
	return 17;

    }else{
      std::cerr << "Error: enter a valid SBS configuration." << std::endl;
      return -1;
    }
    return -1;
  }

  //Angle hadron calorimeter makes wrt exit beamline in degrees
  Double_t hcaltheta(std::string experiment,Int_t config){
    if( experiment.compare("gmn")==0 ){

      if(config==1)
	return 33.5;
      else if(config==4)
	return 31.9;
      else if(config==7)
	return 16.1;
      else if(config==11)
	return 13.3;
      else if(config==14)
	return 17.3;
      else if(config==8)
	return 29.4;
      else if(config==9)
	return 22.0;
      else if(config==4363)
	return 21.6;

    }else if( experiment.compare("gen")==0 ){ //This section to be finished as new experiments are performed

      if(config==1||config==2)
	return 34.7;
      else if(config==3)
	return 21.6;
      else if(config==45||config==4)
	return 18;

    }else{
      std::cerr << "Error: enter a valid SBS configuration." << std::endl;
      return -1;
    }
    return -1;
  }

  //Timestamp for each configuration
  std::string sbsts(std::string experiment,Int_t config){
    if( experiment.compare("gmn")==0 ){

      if(config==1)
	return "none";
      else if(config==4)
	return "--------[ 2021-10-18 05:30:00 ]";
      else if(config==7)
	return "--------[ 2021-11-08 17:00:00 ]";
      else if(config==11)
	return "--------[ 2021-11-17 08:00:00 ]";
      else if(config==14)
	return "--------[ 2022-01-10 22:30:00 ]";
      else if(config==8)
	return "--------[ 2022-01-20 07:30:00 ]";
      else if(config==9)
	return "--------[ 2022-02-01 17:42:00 ]";
      else if(config==4363)
	return "-------[ 2021-10-10 00:00:00 ]";

    }else if( experiment.compare("gen")==0 ){ //This section to be finished as new experiments are performed
     
      if(config==1||config==2)
	return "--------[ 2022-10-17 00:00:00 ]";
      else if(config==3)
	return "--------[ 2022-11-11 00:00:00 ]";
      else if(config==4)
	return "--------[ 2023-01-15 00:00:00 ]";
      else if(config==45)
	return "--------[ 2023-09-01 00:00:00 ]";

    }else{
      std::cerr << "Error: enter a valid SBS configuration." << std::endl;
      return "";
    }
    return "";
  }

  //target index for analysis tree
  Int_t tidx(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 1;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 2;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 3;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //length of the target in meters
  Double_t ltgt(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 0.15;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 0.15;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 0.30;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //target density in g/cc
  Double_t tarrho(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 0.0723;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 0.169;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 0.26;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //cell density in g/cc
  Double_t crho(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 2.71;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 2.71;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 2.71;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //cell diameter in cm
  Double_t cdiam(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 4.064;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 4.064;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 4.064;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //collisional stopping power of cell wall in GeV*m/g
  Double_t cdEdx(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 0.0021;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 0.0021;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 0.0021;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //cell wall thickness in m
  Double_t cthick(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 0.02;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 0.02;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 0.5;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //Upstream cell wall thickness in m
  Double_t uwallthick(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 0.0145;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 0.0145;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 0.5;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //Downstream cell wall thickness in m
  Double_t dwallthick(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 0.015;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 0.015;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 0.5;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //target collisional stopping power in GeV*m/g
  Double_t dEdx(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 0.00574;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return 0.01148;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return 0.01722;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

  //average target mass
  Double_t M_t(std::string tar){
    if( tar.compare("lh2")==0 || tar.compare("LH2")==0 ){
      return 0.938272;
    }else if( tar.compare("ld2")==0 || tar.compare("LD2")==0 ){
      return (0.938272+0.939565)/2;
    }else if( tar.compare("he3")==0 || tar.compare("He3")==0 || tar.compare("HE3")==0 ){
      return (0.938272+0.938272+0.939565)/3;
    }else{
      std::cerr << "Error: enter a valid target." << std::endl;
      return -1;
    } 
    return -1;
  }

} //::hcal
