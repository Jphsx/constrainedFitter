#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "gear/BField.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "MassConstraintFitter.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;
typedef CLHEP::Hep3Vector Vector3D ;
//
// MarlinKinfit stuff
#include "LeptonFitObject.h"
#include "JetFitObject.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"

// Marlin stuff
#include <marlin/Global.h>
// the event display

// ROOT stuff
#include "TMath.h"
#include "TMatrixD.h"

#include <cstdlib>
#include <cmath>

using namespace lcio;

MassConstraintFitter aMassConstraintFitter;

//////////////////Adapted from GammaGammaCandidateFinder.cc//////////////////////////////
MassConstraintFitter::MassConstraintFitter() : marlin::Processor("MassConstraintFitter") {

  registerProcessorParameter( "Printing" , 
			      "Print certain messages"  ,
			      _printing,
			       (int)1 ) ;

  registerProcessorParameter("RootFile" ,
			     "Name of the output root file" ,
			     m_rootFile,
			     std::string("MassConstraintFitAnalysis.root"));

  registerProcessorParameter("parentPdg",
			     "PDG code of parent particle",
			     _parentPdg,
			     (int) 221);
  
  registerProcessorParameter("parentMass",
			     "Parent particle mass [GeV]",
			     _parentMass,
			     (double) 0.547862);
  
  registerProcessorParameter("parentCharge",
			     "Parent particle charge",
			     _parentCharge,
			     (double) 0.0);

  registerProcessorParameter("nDaughters",
			    "Number of daughter decay particles",
			     _nDaughters,
			     (int) 3);

  registerProcessorParameter("nCharged",
			    "Number of daughter decay particles with charge !=0",
			     _nCharged,
			     (int) 2);

  registerProcessorParameter("nNeutral",
			     "Number of daughter decay particles with charge ==0",
			     _nNeutral,
			     (int) 1);

  std::vector<int> daughterChargedPdgs;
  daughterChargedPdgs.push_back(211);
  registerProcessorParameter("daughterChargedPdgs",
			     "PDG codes of daughter decay particles with charge != 0",
			     _daughterChargedPdgs,
			      daughterChargedPdgs);

  std::vector<int> daughterNeutralPdgs;
  daughterNeutralPdgs.push_back(22);
  registerProcessorParameter("daughterNeutralPdgs",
			     "PDG codes of daughter decay particles with charge == 0",
			     _daughterNeutralPdgs,
			      daughterNeutralPdgs);
			     
  std::vector<float> daughterChargedMass;
  daughterChargedMass.push_back(0.13957018);
  registerProcessorParameter("daughterChargedMass",
			     "Mass [GeV] of daughter decay particles with charge != 0",
			     _daughterChargedMass,
			      daughterChargedMass); 

  std::vector<float> daughterNeutralMass;
  daughterNeutralMass.push_back("0.0");
  registerProcessorParameter("daughterNeutralMass",
			     "Mass [GeV] of daugher decay particles with charge == 0",
			     _daughterNeutralMass,
			      daughterNeutralMass);

  std::string inputParticleCollectionName = "PandoraPFOs";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName" , 
			     "Input Particle Collection Name "  ,
			     _inputParticleCollectionName,
			      inputParticleCollectionName);

  std::string inputTrackCollectionName = "MarlinTrkTracks";
  registerInputCollection( LCIO::TRACK,
				"InputTrackCollectionName" ,
				"Input Track Collection Name " ,
				_inputTrackCollectionName,
				inputTrackCollectionName);
  
  registerInputCollection( LCIO::MCPARTICLE,
				"MCParticleCollection" ,
				"Name of the MCParticle input collection" ,
				_mcParticleCollectionName,
				std::string("MCDecayParticles"));

 
  std::string outputParticleCollectionName = "MassConstraintCandidates";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "OutputParticleCollectionName" , 
			     "Output Particle Collection Name "  ,
			     _outputParticleCollectionName,
			     outputParticleCollectionName);

  std::string outputTrackCollectionName = "MassConstraintCandidates";
  registerOutputCollection( LCIO::TRACK,
			    "OutputTrackCollectionName",
			    "Output Particle Collection Name" ,
			    _outputParticleCollectionName,
			    outputTrackCollectionName);

  registerProcessorParameter( "FitProbabilityCut" , 
			      "Minimum fit probability"  ,
			      _fitProbabilityCut,
			      (double)0.001);  

  registerProcessorParameter( "fitter" ,
                              "0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
                              _ifitter,
                              (int)0);

  registerProcessorParameter("FitAnalysis" ,
			      "0 = No TTree , 1 = make TTree,  generates a TTree of fit and measured particle parameters along with pull distributions",
			      _fitAnalysis,
			      (int)1);
  
  registerProcessorParameter("GeneratorAnalysis",
			     "0 = No TTree , 1 = make TTree, includes generator level information in the TTree if available",
			     _genAnalysis,
			     (int)1);

  
  return;

}

//===================================================================================

void MassConstraintFitter::init() {
  if(_printing>1)printParameters(); 

  evtNo=0;

  if(_fitAnalysis)
	rootFile = new TFile(m_rootFile.c_str(),"RECREATE");
 
////////////////////////////////////////////////////////
///fit and measured analysis tree init
//_fitAnalysis includes the general analysis from the reconstructed particles
//it produces a ttree that contains the reconstructed (measured)  particles, fit particles, and pull distributions between
//the fit and measured particles
////////////////////////////////////////////////////////
//initialize tree and variables 
  if(_fitAnalysis){
	tree = new TTree("tree", "tree");
  	tree->SetDirectory(rootFile);

  	//Some useful plots
  	tree->Branch("RecoEnergy", &RecoEnergy);
  	tree->Branch("FitEnergy", &FitEnergy);
        tree->Branch("RecoMass", &RecoMass);
  	tree->Branch("FitProbability", &FitProbability );
	tree->Branch("Chisq", &Chisq);
	tree->Branch("evtNo", &evtNo);

	tree->Branch("measNeutral.", "vector<TLorentzVector>", &measNeutral);
        tree->Branch("measCharged.", "vector<TLorentzVector>", &measCharged);
	tree->Branch("measTrack.","vector<vector<double> >", &measTrack);
        tree->Branch("fitNeutral.", "vector<TLorentzVector>", &fitNeutral);
        tree->Branch("fitCharged.", "vector<TLorentzVector>", &fitCharged);
	tree->Branch("fitTrack.", "vector<vector<double> >",&fitTrack);

	tree->Branch("measNeutral_err.", "vector<vector<double> >", &measNeutral_err);
	tree->Branch("measCharged_err.", "vector<vector<double> >", &measCharged_err);
	tree->Branch("measTrack_err.","vector<vector<double> >",&measTrack_err);
	tree->Branch("fitNeutral_err.", "vector<vector<double> >", &fitNeutral_err);
	tree->Branch("fitCharged_err.", "vector<vector<double> >", &fitCharged_err);
	tree->Branch("fitTrack_err.", "vector<vector<double> >",&fitTrack_err);

	tree->Branch("measNeutralPdg.", "vector<double>", &measNeutralPdg);
	tree->Branch("measChargedPdg.", "vector<double>", &measChargedPdg);

	tree->Branch("fitmeas_NeutralPulls.", "vector<vector<double> >", &fitmeas_NeutralPulls);
	tree->Branch("fitmeas_ChargedPulls.", "vector<vector<double> >", &fitmeas_ChargedPulls);
	
	tree->Branch("measParent.", "TLorentzVector", &measParent);
	tree->Branch("fitParent.", "TLorentzVector", &fitParent);
	
	tree->Branch("measParent_err.", "vector<double>", &measParent_err);
	tree->Branch("fitParent_err.", "vectore<double>", &fitParent_err);

	tree->Branch("fitmeas_NeutralPulls.", "vector<vector<double> >",&fitmeas_NeutralPulls);
	tree->Branch("fitmeas_ChargedPulls.", "vector<vector<double> >",&fitmeas_ChargedPulls);
   }
	//////////////mc stuff
	if(_genAnalysis){
	tree->Branch("genNeutral.", "TLorentzVector", &genNeutral);
	tree->Branch("genCharged.", "TLorentzVector", &genCharged);
	
	tree->Branch("measgen_NeutralPulls.", "vector<vector<double> >", &measgen_NeutralPulls);
	tree->Branch("measgen_ChargedPulls.", "vector<vector<double> >", &measgen_ChargedPulls);
	tree->Branch("fitgen_NeutralPulls.", "vector<vector<double> >" , &fitgen_NeutralPulls);
	tree->Branch("fitgen_ChargedPulls.", "vector<vector<double> >" , &fitgen_ChargedPulls);

	tree->Branch("genParent.", "TLorentzVector", &genParent);
	
	tree->Branch("genNeutralPdg.", "vector<double>", &genNeuralPdg);
	tree->Branch("genChargedPdg.", "vector<double>", &genChargedPdg);
	
    }


   ////////////////////////////////////////////////////////////////////////
  //generator tree init
  //the generator analysis houses the measured particles, the monte carlo particles
  //and pull distributions between the measured/fit particles and monte carlo particles
  //only properly matched particles are added to the tree, match is governed
  //by cuts on the fractional error between the measured and mc particle
  /////////////////////////////////////////////////////////////////////
   //initialize generator tree 
   
    
 
  return;
}

//===================================================================================

void MassConstraintFitter::processRunHeader( LCRunHeader* run) { 
  return;
}

//===================================================================================

void MassConstraintFitter::processEvent( LCEvent * evt ) { 

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  // Access PFO collection
  bool found = this->FindPFOs(evt);
  bool found2 = this->FindTracks(evt);
  bool found3 = 0;

 if(_genAnalysis){
   found3 = this->FindMCParticles(evt);
 }

  if((found && found2) || (found && found && _genAnalysis && found3)){
    if(_printing>1)std::cout << "Analysis of resonance with PDG code: " << _parentPdg << std::endl; 
    this->FindMassConstraintCandidates(recparcol);
  }} 
  
  
  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  return;
  
}
//===================================================================================
void MassConstraintFitter::end(){
 if(_fitAnalysis){
 	rootFile->cd();
 	tree->Write();
        genTree->Write(); 
 }
  return;
}
//===================================================================================
bool MassConstraintFitter::FindPFOs( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _pfovec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){
 // std::cout<<"collection item: "<< *name << " with " << _inputParticleCollectionName<<std::endl;    
    if(*name==_inputParticleCollectionName){ 
//	std::cout<< "found matching collection name for photon" <<std::endl;
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
//	std::cout<<" number of elements "<<nelem<<std::endl;
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	_pfovec.push_back(recoPart);
      }
    }
  }

  if(_printing>1)std::cout << "Find PFOs : " << tf << std::endl; 

  return tf;
}
//===================================================================================
bool MassConstraintFitter::FindTracks( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _trackvec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){    
    if(*name==_inputTrackCollectionName){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	Track* track = dynamic_cast<Track*>(col->getElementAt(i));
	_trackvec.push_back(track);
      }
    }
  }

  if(_printing>1)std::cout << "FindTracks : " << tf << std::endl; 

  return tf;
}
//===================================================================================  
bool MassConstraintFitter::FindMCParticles( LCEvent* evt ){
	  bool tf = false;

  // clear old vector
  _mcpartvec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){    
    if(*name==_mcParticleCollectionName){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	MCParticle* mcPart = dynamic_cast<MCParticle*>(col->getElementAt(i));
	_mcpartvec.push_back(mcPart);
      }
    }
  }

  if(_printing>1)std::cout << "Find MCParticles : " << tf << std::endl; 

  return tf;
}

int MassConstraintFitter::getCorrespondingMCParticleIndex(TLorentzVector rec){
	 if(_mcpartvec.size() > 1) return -1 ;
        if(_mcpartvec.size() == 0) return -1;
        int closest_match_index=-1;
        double theta_residual=-1;
        double phi_residual=-1;
        double e_residual=0;
        double tempresidual1=-1;
        double tempresidual2=-1;
        double tempresidual3=0;//was going to be residual for energy but dont include because no distinguishing between neutral and charged particles
       TLorentzVector mc;
      
        for(int i=0; i<_mcpartvec.size(); i++){
                //compare angles
                mc.SetPxPyPzE(_mcpartvec[i]->getMomentum()[0],_mcpartvec[i]->getMomentum()[1],_mcpartvec[i]->getMomentum()[2],_mcpartvec[i]->getEnergy());
		if((mc.M()-rec.M()) < 0.001){//mass check, look at difference incase of rounding. Is it within 5%? dont divide because photons
                	if(closest_match_index==-1){
                        	closest_match_index = 0;
                        	theta_residual = abs(rec.Theta() - mc.Theta());
                        	phi_residual = abs(getPhi(rec) - getPhi(mc));
                    	//    e_residual = abs(rec.E() - mc.E());
                	}
                	tempresidual1 = abs(rec.Theta() - mc.Theta());
                	tempresidual2 = abs(getPhi(rec) - getPhi(mc));
             	//   tempresidual3 = abs(rec.E() - mc.E());
                	if((tempresidual1+tempresidual2+tempresidual3) < (theta_residual+phi_residual+e_residual) ){
                	        closest_match_index=i;
                	        theta_residual = tempresidual1;
                	        phi_residual = tempresidual2;
                	    //    e_residual = tempresidual3;
                	}
		}

        }
	if(closest_match_index != -1){
		mc.SetPxPyPzE(_mcpartvec[closest_match_index]->getMomentum()[0],_mcpartvec[closest_match_index]->getMomentum()[1],_mcpartvec[closest_match_index]->getMomentum()[2],_mcpartvec[closest_match_index]->getEnergy());
        	if(_printing>3){
			std::cout<<"MC Match: "<<std::endl;
			std::cout<<"Reco (E,theta,phi) "<<rec.E()<<", "<<rec.Theta()<<", "<<rec.Phi()<<" "<<std::endl;
			std::cout<<"MC   (E,theta,phi) "<<mc.E()<<", "<<mc.Theta()<<", "<<mc.Phi()<<" " <<std::endl;
		}
	}
	else{
		if(_printing>3){
		std::cout<<"Particle not matched "<<std::endl;
		std::cout<<"Reco (E,theta,phi)"<<rec.E()<<", "<<rec.Theta()<<", "<<rec.Phi()<<" "<<std::endl;
		}
	}
        return closest_match_index;

}

//===================================================================================
//returns error array with parameterization {dE, dTheta, dPhi}
//error for photon energy is sigma/E = 0.18/sqrt(E)
//error on photon angle is sigma/E = 0.001/sqrt(E)
//i need to get errors from the recoparticle error matrix probably not make my own?
std::vector<double> MassConstraintFitter::getNeutralErrors(TLorentzVector pneutral, ReconstructedParticle* pNeutral ){

	std::vector<double> errors;/// = new double[3];
	if(pNeutral->getType() == 22){
	 	errors.push_back(0.16*std::sqrt(pneutral.E()) );
		errors.push_back(0.001/std::sqrt(pneutral.E()) );
		errors.push_back(0.001/std::sqrt(pneutral.E()) );
	}
	else{
		errors.push_back(0.55*std::sqrt(pneutral.E()) );
		errors.push_back(0.001/std::sqrt(pneutral.E()) );
		errors.push_back(0.001/std::sqrt(pneutral.E()) );
	}
	return errors;

}
//===================================================================================
//returns error array with parameterization {dk, dTheta, dPhi} with the standard parameterization used in the code
//errors come from track record with (d0,z0,omega,tanLambda,theta, lamda)
std::vector<double> MassConstraintFitter::getChargedParticleErrors(TLorentzVector pcharged, Track* ptrk){
	std::vector<double> errors; 

	errors.push_back(std::sqrt(ptrk->getCovMatrix()[0];//d0
	errors.push_back(std::sqrt(ptrk->getCovMatrix()[5];//phi
	errors.push_back(std::sqrt(ptrk->getCovMatrix()[9];//ome
	errors.push_back(std::sqrt(ptrk->getCovMatrix()[12];//z0
	errors.push_back(std::sqrt(ptrk->getCovMatrix()[14];//tan

	return errors;
}
//===================================================================================
void MassConstraintFitter::PrintCov(FloatVec cov, int dim){
	int index = 0;
	for(int i=0; i<dim; i++){
		std::cout<< " | " ;
		for(int j=0; j<dim; j++){
				std::cout<<cov[index]<<" ";
				index++;
		}
		std::cout<<" | " << std::endl;
	}
	std::cout<< std::endl;
}
//===================================================================================
void MassConstraintFitter::PrintCov(double* cov, int dim){
	int index = 0;
	for(int i=0; i<dim; i++){
		std::cout<< " | " ;
		for(int j=0; j<dim; j++){
				std::cout<<cov[index]<<" ";
				index++;
		}
		std::cout<<" | " << std::endl;
	}
	std::cout<<std::endl;
}
//===================================================================================
//take in the fit covariance matrix and convert it to global fit error vectors for readability
//*updated now takes in cov, and based on #of particles/neutrals/charged builds the error data structures for fit particles
//it goes charged and then neutrals
void MassConstraintFitter::setFitErrors(double* cov, int dim){
//using the big fit cov matrix, extract the diagonall and put the fit errors onto global vectors
	//diagonals are every 10 indices so add 10
	std::vector<double> temp_err;
	int chargeCounter=0;
	for(int i=0; i<dim*dim; i = i+10){
		temp_err.push_back(cov[i]);
		if( temp_err.size() == 3){
			if(chargeCounter < chargedFO.size()){
				fitCharged_err.push_back(temp_err);
				chargeCounter++;
				temp_err.clear();
			}
			else{
				fitNeutral_err.push_back(temp_err);
				temp_err.clear();
			}
		}
	}
 

}
//===================================================================================
//manually constructs the  covariance matrix for parent particle, meant to be transformed into 4x4 4vector form
//***********This is to be used after the global measured errors vectors are  populated**************************
double* MassConstraintFitter::ConstructParentMeasCovMatrix(){
//concat vectors	
	std::vector<double> params;  //THIS may be missing START/END param

	for(int i=0; i< measCharged_err.size(); i++){
		params.insert(params.end(), measCharged_err.at(i).begin(), measCharged_err.at(i).end());
	}
	for(int i=0; i< measNeutral_err.size(); i++){
		params.insert(params.end(), measNeutral_err.at(i).begin(), measNeutral_err.at(i).end());
	}
	int dim = params.size();

	double* meascov = new double[dim*dim];
	int index = 0;
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			if(i==j){
				meascov[index]= params[i]*params[j];
			}
			else{
				meascov[index] = 0;
			}
			index ++;
		}
	}
//	std::cout<<"DEBUGGING PRINT ENTIRE MATRIX "<< std::endl;
//	for(int i= 0; i<dim*dim; i++){
//		std::cout<< meascov[i] << " ";
//		if(i+1 % 9 == 0) std::cout<<std::endl; 
//	}

return meascov;	
}
//===================================================================================
//take in the fit and measured covariance matrices (4x4) 4 vector form, and push diagonals onto vectors to be stored in tree
//convention is [dPx, dPy, dPz, dE]
void MassConstraintFitter::setParentErrors(FloatVec meascov, FloatVec fitcov){

	measparent_err.push_back( std::sqrt(meascov[0]) );
	measparent_err.push_back( std::sqrt(meascov[5]) );
	measparent_err.push_back( std::sqrt(meascov[10]) );
	measparent_err.push_back( std::sqrt(meascov[15]) );
	
	fitparent_err.push_back( std::sqrt(fitcov[0]) );
	fitparent_err.push_back( std::sqrt(fitcov[5]) );
	fitparent_err.push_back( std::sqrt(fitcov[10]) );
	fitparent_err.push_back( std::sqrt(fitcov[15]) );
}
//input 6 parameter track stuff
std::vector<double> MassConstraintFitter::ConstructChargedSubMatrix(std::vector<double> p, TLorentzVector ptlv){
	std::vector<double> submat;
	//submatrix jacobian is this (derivatives of parent wrt to charged parameter


	/*  
		dpx/dd0 dpx/dphi dpx/dome dpx/dz0 dpx/dtanl dpx/dstart
		dpy/dd0 dpy/dphi dpy/dome dpy/dz0 dpy/dtanl dpy/dstart
		dpz/dd0 dpz/dphi dpz/dome dpz/dz0 dpz/dtanl dpz/dstart
		de/dd0 de/dphi de/dome de/dz0 de/dtanl de/dstart

	*/
	double lambda = atan(p.at(4));
	double phi = p.at(1);
	double omega = p.at(2);
	double P = ptlv.P();

	const double c = 2.99792458e8; // m*s^-1
  	const double B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();;          
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = B*c*mm2m*eV2GeV;
	
	

	submat.push_back(0); //dpx/dd0
	submat.push_back(P*cos(lambda)*cos(phi); //dpx/dphi ...
 	submat.push_back(-eB*cos(phi)/(omega*omega)); //dpx/dome
	submat.push_back(0); //dpx/dz0
	submat.push_back(-P*cos(phi)*tan(lambda)/pow(1+(tan(lambda)*tan(lambda)),3.0/2.0); //dpx/dtanl
	submat.push_back(0);// dpx/dstart
	submat.push_back(0); //dpy/dd0
	submat.push_back(P*cos(lambda)*cos(phi)); //dpy/dphi
	submat.push_back(-eB*sin(phi)/(omega*omega)); //dpy/dome
	submat.push_back(0);//dpy/dz0
	submat.push_back(-P*sin(phi)*tan(lambda)/pow(1+(tan(lambda)*tan(lambda)),3.0/2.0); //dpy/dtanl
	submat.push_back(0);//dpy/dstart
	submat.push_back(0);//dpz/dd0'
	submat.push_back(0);//dpz/dphi
	submat.push_back(-eB*tan(lambda)/(omega*omega)); //dpz/dome
	submat.push_back(0);//dpz/dz0
	submat.push_back(1/pow(1+(tan(lambda)*tan(lambda)),3.0/2.0);//dpz/dtanl
	submat.push_back(0);//dpz/dstart
	submat.push_back(0);//de/dd0
	submat.push_back(0);//de/dphi
	submat.push_back(0);//de/dome ?? is this right ?? probably not
	submat.push_back(0);//de/dz0
	submat.push_back(0);//de/dtanl
	submat.push_back(0);//de/dstart
	
	

	return submat;
}
std::vector<double> MassConstraintFitter::ConstructNeutralSubMatrix(TLorentzVector p){
	std::vector<double> submat;
	//submatrix jacobian is this (derivatives of parent(Px,Py,Pz,E) wrt to neutral parameter)
	/*
	  dPx/dE dPx/dtheta dPx/dphi
	  dPy/dE dPy/dtheta dPy/dphi
	  dPz/dE dPz/dtheta dPz/dphi
	  dEp/dEc dEp/dtheta dEp/dphi
	*/ 
	submat.push_back(p.Px()/p.P());
	submat.push_back(p.Pz()*p.Px()/p.Perp());
	submat.push_back(-p.Py());
	submat.push_back(p.Py()/p.P());
	submat.push_back(p.Pz()*p.Py()/p.Perp());
	submat.push_back(p.Px());
	submat.push_back(p.Pz()/p.P());
	submat.push_back(-p.Perp());
	submat.push_back(0.0);
	submat.push_back(1.0);
	submat.push_back(0.0);
	submat.push_back(0.0);

	return submat;
	
}
double* MassConstraintFitter::ConcatSubMatrices(std::vector<std::vector<double> > matrices ){
	//this is the correctly arranged jacobian made in a vector and then copied to the desired type double*
	std::vector<double> jacobian;
	std::vector<std::vector<double>::iterator> its;
	for(int i=0; i<matrices.size(); i++){
		std::vector<double>::iterator it = matrices.at(i).begin();
		its.push_back(it);
	}
	
	while(its.at(_nCharged+_nNeutral-1) != matrices.at(_nCharged+_n_Neutral-1).end()){

		for(int i=0; i<_nCharged; i++){
			if(its.at(i) <= matrices.at(i).end()){
				jacobian.insert(jacobian.end(), its.at(i), its.at(i)+6);
				its.at(i) = its.at(i) + 6;
			}
		}
		for(int i=0; i<_nNeutral; i++){
			if(its.at(i) <= matices.at(i).end()){
				jacobian.insert(jacobian.end(), its.at(i), its.at(i)+3);
				its.at(i) = its.at(i) + 3;
			}
		}

	}

	
	//copy to double* because thats efficient ;)
	double* jacobian2 = new double[jacobian.size()];
	for(int i=0; i<jacobian.size(); i++){
		jacobian2[i] = jacobian.at(i);
	}

	return jacobian2;
	
}
//===================================================================================
//FloatVec MassConstraintFitter::ConstructCovMatrix(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, double*  cov){
FloatVec MassConstraintFitter::ConstructCovMatrix(std::vector<TLorentzVector> charged, std::vector<TLorentzVector> neutral, double* cov){

      
            FloatVec covP(16,0.0);
      //  (px,py,pz,E)

      //  Follow for now implementation in FourMomentumCovMat based on TMatrixD

      // Note some cases of cov_dim != 6 have been seen ...
             const int nrows  = (6*_nCharged+3*_nNeutral);  
             const int ncolumns  =4; 


	std::vector<std::vector<double> > matrices;
	for(int i =0; i< neutral.size(); i++){
		matrices.push_back(ConsructChargedSubMatrix(neutral.at(i));
	}
	for(int i= 0; i< charged.size(); i++){	
		matrices.push_back(ConstructNeutralSubMatrix(charged.at(i));
	}

	double* jacobian = ConcatSubMatrices(matrices);


      //  Now set up to calculate the new covariance matrix, namely  V' = D^T V D
             TMatrixD Dmatrix(nrows,ncolumns, jacobian_by_columns, "F");
             TMatrixD Vmatrix(nrows,nrows, cov, "F");                      // May need to have the filling array explicitly dimensioned ??
             TMatrixD Covmatrix(ncolumns,ncolumns);                        // Hopefully this is good enough for initialization ?

             Covmatrix.Mult( TMatrixD( Dmatrix, TMatrixD::kTransposeMult, Vmatrix) ,Dmatrix);   // Is this what we want ?
	   int index=0;
             for(int ii=0;ii<4;ii++){
                 for(int jj=0;jj<4;jj++){
        //            if(_printing>3)std::cout << "4-vector cov. " << ii << " " << jj << " " << Covmatrix(ii,jj) << " Corr: " << Covmatrix(ii,jj)/std::sqrt(Covmatrix(ii,ii)*Covmatrix(jj,jj)) << std::endl;
        		covP[index] = Covmatrix(ii,jj);
			index++;
                 }                  

	return covP;
} 
//===================================================================================
//Current implementation only allows one type of fitter because using the base class is not working correctly
//BaseFitter* MassConstraintFitter::setUpFit(TLorentzVector gamma, TLorentzVector p1, TLorentzVector p2, Track* p1Track, Track* p2Track ){
//OPALFitterGSL* MassConstraintFitter::setUpFit(TLorentzVector gamma, TLorentzVector p1, TLorentzVector p2, Track* p1Track, Track* p2Track){
BaseFitter* MassConsraintFitter::setUpFit(vector<int> neutralIndices, vector<int> chargedIndices,
					vector<TLorentzVector> pneutral, vector<TLorentzVector> ptrack,
					vector<ReconstructedParticle*> pNeutralVec, vector<Track*> pTrackVec){                                       

	//set up mass constraint
	MassConstraint mc( (double)_parentMass );

	//construct fit objects for gamma, and 2 charged particles
	std::vector<std::vector<double> > neutralErrors;
	std::vector<std::vector<double> > chargedErrors;

	
//Working here today? implement these get error functions
	for(unsigned int i =0; i< neutralIndices.size(); i++){
		neutralErrors.push_back(getNeutralErrors(pNeutralVec.at( neutralIndices.at(i)) ) );//should i pass recopart or tlv here?
	}
	for(unsigned int i=0; i < chargedIndices.size(); i++){
		chargedErrors.push_back(getChargedErrors(ptrack.at( chargedIndices.at(i)), pTrackVec.at( chargedIndices.at(i)) ) );
	}


	//set up JetFitObjects
	for(unsigned int i=0; i < neutralIndices.size(); i++){
		neutralJets.push_back( new JetFitObject(pneutral.at(neutralIndices.at(i)).E(),pneutral.at(neutralIndices.at(i)).Theta(),pneutral.at(neutralIndices.at(i)).Phi(), neutralErrors.at(i)[0], neutralErrors.at(i)[1], neutralErrors.at(i)[2]));
	}					

	for(unsigned int i =0; i < chargedIndices.size(); i++){
		//chargedFO.push_back( new LeptonFitObject( 1/ptrack.at(chargedIndices.at(i)).Perp(), ptrack.at(chargedIndices.at(i)).Theta(), ptrack.at(chargedIndices.at(i)).Phi(), chargedErrors.at(i)[0], chargedErrors.at(i)[1], chargedErrors.at(i)[2]), ptrack.at(chargedIndices.at(i)).M());
		trackFO.push_back( new TrackParticleFitObject(pTrackVec.at(chargedIndices.at(i)),ptrack.at(chargedIndices.at(i)).M()  ));
	} 	


	if(_printing>4)std::cout <<" Fitting Measured Quantities:" <<std::endl;						
	
	if(_printing>4)std::cout <<" Neutral Particles "<<std::endl;
	for(unsigned int i =0; i< neutralIndices.size(); i++){
		if(_printing>4)std::cout << "neutral "<< i <<" E, theta, phi, M, Pdg: " <<  pneutral.at(neutralIndices.at(i)).E()<< " "<< pneutral.at(neutralIndices.at(i)).Theta()<< " "<< pneutral.at(neutralIndices.at(i)).Phi()<< " "<< pneutral.at(neutralIndices.at(i)).M()<< " "<< pNeutralVec.at(neutralIndices.at(i))->getType() << std::endl;
		if(_printing>4)std::cout << "errors "<< i <<" dE, dTheta, dPhi: " << neutralErrors.at(i)[0] << " "<< neutralErrors.at(i)[1] << " "<<neutralErrors.at(i)[2]<< std::endl;
	}	
	if(_prining>4)std::cout <<" Charged Particles "<<std::endl;
	for(unsigned int i=0; i< chargedIndices.size(); i++){
		if(_printing>4)std::cout <<"charged "<< i <<" d0, phi, omega, z0, tanL, M: "<< pTrackVec.at(i)->getD0()<<" "<<pTrackVec.at(i)->getPhi()<<" "<<pTrackVec.at(i)->getOmega()<<" "<<pTrackVec.at(i)->getZ0()<<" "<<pTrackVec.at(i)->getTanLambda()<<" "<<ptrack.at(chargedIndices.at(i)).M() << std::endl;
		if(_printing>4)std::cout<<"errors "<< i << "dd0 dPhi ddome dz0 dtanL ";
		for(unsigned int j=0; j<chargedErrors.at(i).size(); j++){
			if(_printing>4)std::cout<< chargedErrors.at(j) << " ";
		}
		if(_printing>4) std::cout<<std::endl;
	}


	//clean up memory with errors that are not going to be used again in this scope
	neutralErrors.clear();
	chargedErrors.clear();
					
	//add fitobjects to mass constraint
	for(unsigned int i= 0; i<chargedIndices.size(); i++){
		mc.addtoFOList(*chargedFO[i]);
	}
	for(unsigned int i =0; i<neutralIndices.size(); i++){
		mc.addtoFOList(*neutralJets[i]);
	}

	//declare fitter
	//include options for other fitters
	OPALFitterGSL * fitter = new OPALFitterGSL();
	///////////////////	

	//add constraints and fit objects to fitter
	for(unsigned int i=0; i< chargedIndices.size(); i++){
		fiter->addFitObject( chargedFO[i] );
	}
	for(unsigned int i=0; i< neutralIndices.size(); i++){
		fitter->addFitObject( neutralJets[i] );
	}
	//add the things, what order should i use? charged then neutral?
	fitter->addConstraint(mc);
	
	//do the fit (this also returns a fit probability)
	fitter->fit();

	//TODO:: look up the function that returns fit prob, dont fit twice if printing	
//	if(_printing>4)std::cout <<" Fit Probability: "<< fitter->fit()<<std::endl;

	
	return fitter;

}
void MassConstraintFitter::generateSubsets(std::vector<int> v, int k, int start, int currLen, std::vector<bool> used, std::vector<vector<int> >& combinations){
	if(currLen == k) {
		std::vector<int> indices;
		for(unsigned int i=0; i<v.size(); i++){
			if(used.at(i) == true) {
				indices.push_back(i);
			}
		}
		combinations.push_back(indices);
		return;
	}
	if (start == v.size()){
		return;
	}
	used[start] = true;
	generateSubset(v,k,start+1,currLen+1,used, combinations);

	used[start] = false;
	generateSubset(v,k,start+1,currLen,used, combinations);
}
//n choose k
std::vector<std::vector<int> > MassConstraintFitter::generateIndicesCombinations(int vectorsize, int nparticles){
	int n = vectorsize;
	int k = nparticles;
		std::vector<int> master;
		std::vector<bool> used(n);
		std::vector<vector<int> > combinations;
		for(int i=0; i<n; i++){
			master.push_back(i);
			used.at(i)=false;
		}
	generateSubsets(master,k,0,0,used,combinations);
	return combinations;
}
//clone the old RP and update the momentum/Energy from the fit
//these construction methods need the fit errrors too
ReconstructedParticleImpl* MassConstraintFitter::constructFitParticle(TLorentzVector fitp, ReconstructedParticle* measrp){
	ReconstructedParticleImpl* fitp = new ReconstructedParticleImpl();
	
	fitrp = measrp->clone();
	fitrp->setEnergy(fitp.E());
	double* mom = new double[3];
	mom[0]=fitp.Px();
	mom[1]=fitp.Py();
	mom[2]=fitp.Pz();
	fitrp->setMomentum(mom);

	return fitrp;
	
}
//TODO ADD NEW FIT COV MATRIX
Track* MassConstraintFitter::constructTrack(TLorentzVector fitp, Track* meast){
	Track* fitt = new Track();
	fitt = meast.clone();
	const double c = 2.99792458e8; // m*s^-1
  	const double B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();;          
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = B*c*mm2m*eV2GeV;
	//update 
	double sinlambda = fitp.Pz()/fitp.P();
	double coslambda = fitp.Px()/(fitp.P()*cos(fitp.Phi()));
	double tanlambda = sinlambda/coslambda;
	double omega = eB/(fitp.P()*coslambda);
	if(t->getOmega() < 0){
		omega = -omega;
	}
	fitt->setPhi(fitp.Phi());
	fitt->setTanLambda(tanlambda);
	fitt->setOmega(omega);

	return fitt;
}
std::vector<double> MassConstraintFitter::buildTrackVector(Track* t){

	std::vector<double> tvect;
	tvect.push_back(t->getD0());
	tvect.push_back(t->getPhi());
	tvect.push_back(t->getOmega());
	tvect.push_back(t->getZ0());
	tvect.push_back(t->getTanLambda());
	tvect.push_back(0.0);//sStart not measured
	return tvect;
}
std::vector<double> MassConstraintFitter::buildFitTrackVector(TrackParticleFitObject* tfo){
	std::vector<double> tvect;
	// 6 parameters d0 phi omega z0 tanL sStart||sEnd

	for(int i=0; i<6; i++){
		tvect.push_back(tfo->getParam(i));
	}
	return tvect;
}
//===================================================================================
void MassConstraintFitter::FindMassConstraintCandidates(LCCollectionVec * recparcol) {

  if(_printing>1)std::cout << "FindMassConstraintCandidates : (nPFOs = " << _pfovec.size() << " )" << std::endl; 
  if(_printing>1)std::cout << "FindMassConstraintCandidates : evtNo = " << evtNo << std::endl;
  // Look for candidates
 
  std::vector<TLorentzVector>pneutral;
  std::vector<TLorentzVector>ptrack;

  std::vector<ReconstructedParticle*>pNeutralVec;
  std::vector<Track*>pTrackVec;


//neutral loop, load collection photons into local variables
  for(unsigned int i=0;i<_pfovec.size();i++){

	if(_pfovec[i]->getCharge() == 0 ){
		 if(_printing>1)std::cout << "FindMassConstraintCandidates : " << _pfovec[i]->getType() << std::endl;

	         TLorentzVector tempNeutral(_pfovec[i]->getMomentum()[0], _pfovec[i]->getMomentum()[1], _pfovec[i]->getMomentum()[2], _pfovec[i]->getEnergy() );
		
		if(_printing>1)std::cout << "Neutral candidate (px,py,pz,E) "<<tempNeutral.Px()<<" "<<tempNeutral.Py()<<" "<<tempNeutral.Pz()<<" "<<pGamma.E()<<std::endl;
		if(_printing>1)std::cout << "Neutral candidate (E, theta, phi) "<< tempNeutral.E()<<" "<<tempNeutral.Theta()<<" "<<tempNeutral.Phi()<<std::endl;

		/*if(_usePhotonCalibration){//check for ->getType()==22
			calibratePhoton(pGamma);
			if(_printing>1)std::cout <<"g candidate calibration (px,py,pz,E) " <<pGamma.Px()<<" "<<pGamma.Py()<<" "<<pGamma.Pz()<<" "<<pGamma.E()<<std::endl;
			if(_printing>1)std::cout <<"g candidate calibration (E,theta,phi) "<<pGamma.E()<<" "<<pgamma.Theta()<<" "<<pgamma.Phi()<<std::endl;
		}*/

		pneutral.push_back(tempNeutral);
		pNeutralVec.push_back(_pfovec[i]);
	}
  }//end pfovec loop
//start track loop



  const double c = 2.99792458e8; // m*s^-1
  const double B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();;          
  const double mm2m = 1e-3;
  const double eV2GeV = 1e-9;
  const double eB = B*c*mm2m*eV2GeV;
 //
 for(unsigned int i=0; i<_trackvec.size();i++){
	for(unsigned int j=0; j<_daughterChargedMass.size(); j++{
	
		double cosLambda = 1 / std::sqrt(1 + _trackvec[i]->getTanLambda()*_trackvec[i]->getTanLambda() );
		double P = (eB/fabs(_trackvec[i]->getOmega()))/cosLambda;
		double sinLambda = _trackvec[i]->getTanLambda()*cosLambda;
		double cosPhi = cos(_trackvec[i]->getPhi());
		double sinPhi = sin(_trackvec[i]->getPhi());
		double px = P*cosLambda*cosPhi;
		double py = P*cosLambda*sinPhi;
		double pz = P*sinLambda;
		double E = std::sqrt( _daughterChargedMass[j]*_daughterChargedMass[j] + px*px + py*py + pz*pz );
		

			
			pTrackVec.push_back(_trackvec[i]);

			TLorentzVector trk(px,py,pz,E);

			if(_printing>1)std::cout << "-+ candidate "<< i <<" (px,py,pz,E,M)  "<< px << " " << py << " " << pz << " " << E << " " << M << std::endl;
			if(_printing>1)std::cout << "-+ candidate "<< i << " (d0,phi,k,z0,tanLambda) " <<_trackvec[i]->getD0()<<" "<<_trackvec[i]->getPhi()<< " "<<_trackvec[i]->getOmega()<< " "<< _trackvec[i]->getZ0()<< " "<<_trackvec[i]->getTanLambda()<<std::endl;
			
			ptrack.push_back(trk);		

	
	}
	
}
  	
  
  if(_printing>1)std::cout << "FindMassConstraintCandidates : (nNeutral = " << pNeutralVec.size() << " " << pneutral.size() << " )" << std::endl;   
  if(_printing>1)std::cout << "FindMassConstraintCandidates : (n+- = "<< pTrackVec.size() << " " << ptrack.size() << " )" <<std::endl;
  // loop over pairs of candidate photons and keep all combinations consistent with the given parent mass assumption

   //size check
    //if( pgamma.size() >= 1 && pplus.size() >= 1 && pminus.size() >= 1){
    if(pneutral.size() >= _nNeutral && _trackvec.size() >=  _nCharged){  //start back up here need to make ijk arrays with indices for both charged and neutral candidates
    

//store the best fit prob i,j,k then after iterating through all permutations recompute fit and compute analysis vars
	double fitprobmax = -1;
	//double candidate_ijk[3] = {-1,-1,-1};
        vector<vector<int> > neutralCandidateIndices = generateIndicesCombinations(pneutral.size(), _nNeutral);
	vector<vector<int> > chargedCandidateIndices = generateIndicesCombinations(ptrack.size(), _nCharged);
	double fitprob = -1;

	vector<int> bestCandidatesNeutral;
	vector<int> bestCandidatesCharged;
	
	//iterate and fit all particle combinations, save the best one (highest fit probability)
	double chargeSum= 0.0;
	for(unsigned int i=0; i< neutralCandidateIndices.size(); i++){
		for(unsigned int j=0; j< chargedCandidateIndices.size(); i++){
			//iterate over charged particles and find total charge
			for( unsigned int k=0; k<chargeCandidateIndices.at(j).size(); k++){
				 (pTrackVec[chargeCandidateIndices.at(j).at(k)].getOmega() < 0 ) ? chargeSum+=-1.0 : chargeSum+=1.0) ;
			}
			//if charge is consistent with parent, try fitting
			//TODO check and see if neutral types are consistent?
			if(chargeSum == _parentCharge){
				OPALFitterGSL* fitter = setUpFit(neutralCandidateIndices.at(i), chargedCandidateIndices.at(j), pneutral, ptrack, pNeutralVec, pTrackVec);
			}
			if(fitprobmax == -1 && (fitprob > _fitPobabilityCut){
				fitprobmax = fitprob;
				bestCandidatesNeutral = neutralCandidateIndices.at(i);
				bestCandidatesCharged = chargedCandidateIndices.at(j);
			}
			if( (fitprob > fitprobmax ) && (fitprob > _fitProbabilityCut)){
				fitprobmax = fitprob;
				bestCandidatesNeutral = neutralCandidateIndices.at(i);
				bestCandidatesCharged = chargedCandidateIndices.at(j);
			}
			neutralJets.clear();
			chargedFO.clear();
			delete fitter;
		}
	}
     
	std::cout<<"loop completed"<<std::endl;
	//recompute params for best fit probability
//if nothing makes it for this event just return
	if(bestCandidatesNeutral.size()==0 && bestCandidatesCharged.size()==0){
		evtNo++; 
		return;
	}  

	//now refit with best candidates
	OPALFitterGSL* fitter = setUpFit(neutralIndices, chargedIndices, pneutral, ptrack, pNeutralVec, pTrackVec);
	
					
	int cov_dim;
	double * cov = fitter->getGlobalCovarianceMatrix(cov_dim);  //3*(r+k) === 3*(_nNeutral + _nCharged)
	//construct measured cov matrix
//	double * mcov = ConstructParentMeasCovMatrix();
	//if cov matrix dimension is not correct return, fit did not converge
	if(cov_dim == 0){
		evtNo++;
		return;
	}
	
	//build up fit data structures TLVs first
	
	//neutral
	TLorentzVector temp;
	for(int i=0; i<neutralJets.size() ; i++){
		temp.SetPxPyPzE(neutralJets.at(i)->getPx(),neutralJets.at(i)->getPy(), neutralJets.at(i)->getPz(), neutralJets.at(i)->getE());
		fitNeutral.push_back(temp);
	}
	for(int i=0; i<chargedFO.size() ; i++){
		temp.SetPxPyPzE(trackFO.at(i)->getPx(), TrackFO.at(i)->getPy(), TrackFO.at(i)->getPz(), TrackFO.at(i)->getE());
		
		fitCharged.push_back(temp);
		fitTrack.push_back(buildFitTrackVector(TrackFO.at(i));
		
	}
	for(int i=0; i<fitNeutral.size(); i++){
		fitParent += fitNeutral.at(i);	
	}
	for(int i=0; i<fitCharged.size(); i++){
		fitParent += fitCharged.at(i);
	}



	
	//sets an array of charged and neutral errors arrays from cov
	setFitErrors(cov, cov_dim);
	std::vector<ReconstructedParticle*> recoNeutraltemp; // be safe and clear these at the end
	//std::vector<Track*> recoTracktemp;
	//put candidates into local LorentzVector for ease of coding and readability
	//do TLVs first then make new RecoParts to get generate errors
	for(int i=0; i<bestCandidatesCharged.size(); i++){
		measCharged.push_back( ptrack.at(bestCandidatesCharged[i]));
		measTrack.push_back( pTrackVec.at(bestCandidatesCharged[i]));
		//recoTracktemp.push_back(pTrackVec.at(bestCandidatesCharged[i]));
	}
	for(int i=0; i<bestCandidatesNeutral.size(); i++){
		measNeutral.push_back( pneutral.at(bestCandidatesNeutral[i]));
		recoNeutraltemp.push_back( pNeutralVec.at(bestCandidatesNeutral[i]));
	}
	//make reco particles


	//also have error arrays ready locally for readability
	//Photon Errors go {dE,dTheta,dPhi}
	//ChargedParticle errors go {dk, dTheta, dPhi}
	for(int i=0; i<measCharged.size(); i++){
		measCharged_err.push_back( getChargedParticleErrors(measCharged.at(i), recoTracktemp.at(i)) );
	}
	for(int i=0; i<measNeutral.size(); i++){
		measNeutral_err.push_back( getNeutralErrors( measNeutral.at(i), recoNeutraltemp.at(i)) );
	}
	
	//construct measured cov matrix must be called after measured global vectors are populated
	double * mcov = ConstructParentMeasCovMatrix();

	//populate measured  parent errors from covariance matrix
	FloatVec fitParentCov =  ConstructCovMatrix(fitTrack,fitNeutral,cov);  
	FloatVec measParentCov = ConstructCovMatrix(measCharged,measNeutral,mcov);
	setParentErrors(measParentCov, fitParentCov);
	
	


///im to this point
	if(_printing>3){
		std::cout <<"Optimal Candidate Fit Quantities-------------"<<std::endl;
		std::cout << "Constrained fit results  RC: " << fitter->getError() << std::endl;
		std::cout << "Measured mass = " << measParent.M() << std::endl;
		std::cout << "No. of iterations " << fitter->getIterations() << std::endl;
		std::cout << "Fit probability = " << fitprob << std::endl;
		std::cout << "Covariance matrix dimension " << cov_dim << std::endl;
		std::cout << "Candidates Post Fit: "<< std::endl;
		
		for(int i=0; i< fitCharged.size(); i++){
			std::cout << "Track "<< i <<" (d0,phi,omega,z0,tanL,M): "<<std::endl;
			std::cout << " Measured: ";//<< 1/measCharged.at(i).Perp() <<" "<< measCharged.at(i).Theta() <<" "<<measCharged.at(i).Phi()<<" "<<measCharged.at(i).M()<<" ";
			for(int j=0; j<measTrack.size(); j++){
				std::cout<< measTrack.at(i) << " ";
			}
			std::cout << measCharged.at(i).M() <<std::endl;

			std::cout << " Fit: ";// << 1/fitCharged.at(i).Perp() <<" "<<fitCharged.at(i).Theta() <<" "<<fitCharged.at(i).Phi()<<" "<<fitChaged.at(i).M()<<std::endl;
			for(int j=0; j<fitTrack.size(); j++){
				std::cout<< fitTrack.at(i) << " ";
			}
			std::cout << fitCharged.at(i).M() <<std::endl;

			std::cout << "Track "<< i <<" Errors (dd0 dphi domega dz0 dtanL): "<<std::endl;
			std::cout << "Measured: ";
			for(int j=0; j<measCharged_err.at(i).size(); j++){
				std::cout<< measCharged_err.at(i).at(j);
			}
			std::cout <<" Fit: ";
			for(int j=0; j<fitCharged_err.at(i).size(); j++){
				std::cout<< fitCharged_err.at(i).at(j);
			}
			std::cout<<std::endl;

		}

		for(int i=0; i< fitNeutral.size(); i++){
			std::cout << "Neutral "<< i <<" (E,theta,phi,M): "<<std::endl;
			std::cout << " Measured: "<< measNeutral.at(i).E() <<" "<< measNeutral.at(i).Theta() <<" "<<measNeutral.at(i).Phi()<<" "<<measNeural.at(i).M()<<" ";
			std::cout << " Fit: " << fitNeutral.at(i).E() <<" "<<fitNeutral.at(i).Theta() <<" "<<fitNeutral.at(i).Phi()<<" "<<fitNeural.at(i).M()<<std::endl;
			std::cout << "Neutral "<< i <<" Errors (dE, dtheta, dphi): "<<std::endl;
			std::cout << "Measured: ";
			for(int j=0; j<measNeutral_err.at(i).size(); j++){
				std::cout<< measNeutral_err.at(i).at(j);
			}
			std::cout <<" Fit: ";
			for(int j=0; j<fitNeutral_err.at(i).size(); j++){
				std::cout<< fitNeutral_err.at(i).at(j);
			}
			std::cout<<std::endl;

		}
		
		 
		std::cout << "Parent Measured Energy "<< measParent.E() <<" Parent Fit Energy "<< fitParent.E()<<std::endl;
		std::cout <<"----------------------"<<std::endl;

		std::cout<<"Measured Covariance Matrix: "<<std::endl;
		PrintCov(mcov,(6*_nCharged+3*_nNeutral));
		std::cout<<"Fit Covariance Matrix: "<<std::endl;
		PrintCov(cov,6*_nCharged+3*_nNeutral));
		std::cout<<"Measured Parent 4vector 4x4 Covariance Matrix: "<<std::endl;
		PrintCov(measParentCov,4);
		std::cout<<"Fit Parent 4vector 4x4 Covariance Matrix: "<<std::endl;
		PrintCov(fitParentCov,4);					     
	}//end optimal fit printing
				//////////////IM up to this point, add pull stuff				
			              
	if(cov_dim > 0){//condition required because there is not always convergence in the fitting	
		if(fitprob> _fitProbabilityCut){
			//TODO:: add qualifying fits to collection that are less than optimal
			

			if(_parentCharge == 0){
				//parent is neutral store as a reco
				ReconstrucedParticleImpl * recoPart = new ReconstructedParticleImpl();
				//adding old particles to fit particle because lazy and dont want to make new RecoParts*
				for(int i=0; i<fitNeutral.size(); i++){
					//NOTE: these constiuent particles dont have updated covariance matrices only P,E
					recoPart->addParticle(constructFitParticle( fitNeutral.at(i),recoNeutraltemp.at(i)));
				}
				for(int i=0; i<recoTracktemp.size(); i++){
					recoPart->addTrack(constructFitTrack( fitCharged.at(i),rectoTracktemp.at(i)));
				}
				recoPart->setEnergy(fitParent.E());
				double mom[3] = {fitParent.Px(), fitParent.Py(), fitParent.Pz() };
				recoPart->setMomentum(mom);
				recoPart->setMass( _parentMass );
				recoPart->setCharge( _parentCharge );
				recoPart->setCovMatrix( fitParentCov );
				recoPart->setType( _parentPdg );
				recoPart->setGoodnessOfPID((float)fitprob);
				if(_printing>1)std::cout << "candidate passes fit probability cut - KEEP" << std::endl;
				recparcol->addElement( recoPart );
				 
			}
			//if the parent is charged store a track, jacobians have not been derived for this yet
			if(_parentCharge != 0){
				//parent is charged store as a track
				TrackImpl* recoTrack = TrackImpl();
				//recoTrack = constructFitTrack( fitCharged
				
			}
			
			//add stuff to regular ttree
			if(_fitAnalysis){

 				RecoEnergy=parentmeas.E();
				FitEnergy=parentFit.E();
				RecoMass=parentmeas.M();
				FitProbability=fitprob; 
				/*P1_k_FitMeas_pull= (1/p1Fit.Perp() -1/p1meas.Perp())/ std::sqrt( p1meas_err[0]*p1meas_err[0] - cov[0] )  ;
				P1_Theta_FitMeas_pull= (p1Fit.Theta() - p1meas.Theta())/ std::sqrt(p1meas_err[1]*p1meas_err[1] - cov[10]) ;	
				P1_Phi_FitMeas_pull= getPhiResidual(p1Fit.Phi(), p1meas.Phi())/ std::sqrt( p1meas_err[2]*p1meas_err[2] - cov[20]);
		
				P2_k_FitMeas_pull= (1/p2Fit.Perp() - 1/p2meas.Perp())/ std::sqrt( p2meas_err[0]*p2meas_err[0] - cov[30]) ;
				P2_Theta_FitMeas_pull= (p2Fit.Theta() - p2meas.Theta())/ std::sqrt( p2meas_err[1]*p2meas_err[1] - cov[40] );
				P2_Phi_FitMeas_pull= getPhiResidual(p2Fit.Phi(), p2meas.Phi())/ std::sqrt( p2meas_err[2]*p2meas_err[2] - cov[50]);/home/bears/work/research/constrainedFitter/MassConstraintFitter.h
				
		
				Gamma_E_FitMeas_pull= (gammaFit.E() - gammameas.E())/ std::sqrt( gammameas_err[0]*gammameas_err[0] - cov[60] );
				Gamma_Theta_FitMeas_pull= (gammaFit.Theta() - gammameas.Theta())/ std::sqrt( gammameas_err[1]*gammameas_err[1]- cov[70]);
				Gamma_Phi_FitMeas_pull=getPhiResidual(gammaFit.Phi(), gammameas.Phi())/std::sqrt( gammameas_err[2]*gammameas_err[2] - cov[80]) ;*/
			
				//charged pulls fit and measured
				for(int i=0; i<measTrack.size(); i++){
					for(int j=0; j<measTrack.at(i).size(); j++){
						fitmeas_ChargedPulls.push_back( measTrack.at(j)-fitTrack.at(j)/ ( std::sqrt( measCharged_err.at(j)*measCharged_err.at(j) - fitCharged_err.at(j)*fitCharged_err.at(j)) ) );
					}
				}
				//neutral pulls fit andmeasured
				for(int i=0; measNeutral.size(); i++){
					for(int j=0; j<measNeutral.at(i).size();j++){
						fitmeas_NeutralPulls.push_Back( measNeutral.at(j)-fitNeutral.at(j)/ ( std::sqrt( measNeutral_err.at(j)*measNeutral_err.at(j) - fitNeutral_err.at(j)*fitNeutral_err.at(j)) ) );
					}
				}
				
		//		if(Gamma_Phi_FitMeas_pull > 100){
		//			 std::cout<<"crazy!!"<<std::endl;
		//			std::cout<<gammaFit.Phi()<<" "<<gammameas.Phi()<<" "<<getPhiResidual(gammaFit.Phi(), gammameas.Phi()) <<" "<<gammameas_err[2]<<" "<<cov[20]<<std::endl;
		//		}

				fitmeas_ParentPulls.push_back(( fitparent.Px() - measparent.Px())/ std::sqrt(measparent_err[0]*measparent_err[0]-fitparent_err[0]*fitparent_err[0])) ;
				fitmeas_ParentPulls.push_back(( fitparent.Py() - measparent.Py())/ std::sqrt(measparent_err[1]*measparent_err[1]-fitparent_err[1]*fitparent_err[1])) ;
				fitmeas_ParentPulls.push_back(( fitparent.Pz() - measparent.Pz())/ std::sqrt(measparent_err[2]*measparent_err[2]-fitparent_err[2]*fitparent_err[2])) ;
				fitmeas_ParentPulls.push_back(( fitparent.E() - measparent.E())/ std::sqrt(measparent_err[3]*measparent_err[3]-fitparent_err[3]*fitparent_err[3])) ;

				Chisq = fitter->getChi2();
				//tree->Fill();
			}//end fit analysis option
	
			//add stuff to generator ttree
			if(_genAnalysis){

				
			//
			//	std::cout<<"test parent meas vector after tree fill "<<std::endl;
			//	for(int i= 0; i<parentmeas_err.size(); i++){
			//	std::cout<<parentmeas_err[i]<<" ";
			//	}
			//	std::cout<<std::endl;
			//get montecarlo 4 vectors

			//note particle matching is done with mcparticle so only 4vector parameterizations are currently implemented
				//int mcp1index = getCorrespondingMCParticleIndex(p1meas, 1, _daughterPID);
				//int mcp2index = getCorrespondingMCParticleIndex(p2meas, -1, -_daughterPID);
				//int mcpgindex = getCorrespondingMCParticleIndex(gammameas, 0,_neutralPID);
			//these pulls can be made and valid only when each pgetCorrespondingMCParticleIndexarticle has a MC match
				std::vector<int> neutralMCPIndices;
				std::vector<int> chargedMCPIndices;

				for(int i=0; i<fitCharged.size(); i++){
					chargedMCPIndices.push_back(getCorrespondingMCParticleIndex(fitCharged.at(i)) );
				}
				for(int i=0; i<fitNeutral.size(); i++){
					neutralMCPIndices.push_back(getCorrespondingMCParticleIndex(fitNeutral.at(i)) );
				}
				
				
				/*if(mcp1index != -1 && mcp2index != -1 && mcpgindex !=-1){
					 mcp1.SetPxPyPzE(_mcpartvec[mcp1index]->getMomentum()[0],
						_mcpartvec[mcp1index]->getMomentum()[1],
						_mcpartvec[mcp1index]->getMomentum()[2],
						_mcpartvec[mcp1index]->getEnergy());
					 mcp2.SetPxPyPzE(_mcpartvec[mcp2index]->getMomentum()[0],
						_mcpartvec[mcp2index]->getMomentum()[1],
						_mcpartvec[mcp2index]->getMomentum()[2],
						_mcpartvec[mcp2index]->getEnergy()); 
					mcgamma.SetPxPyPzE(_mcpartvec[mcpgindex]->getMomentum()[0],
						_mcpartvec[mcpgindex]->getMomentum()[1],
						_mcpartvec[mcpgindex]->getMomentum()[2],
						_mcpartvec[mcpgindex]->getEnergy());
				*/
			//build up mcp 4vetors
				for(int i=0; i< chargedMCPIndices.size(); i++){
					TLorentzvector mcp;
					if( chargedMCPIndices.at(i) != -1){
						
						mcp.SetPxPyPzE(_mcpartvec[chargedMCPIndices.at(i)]->getMomentum()[0],
						_mcpartvec[chargedMCPIndices.at(i)]->getMomentum()[1],
						_mcpartvec[chargedMCPIndices.at(i)]->getMomentum()[2],
						_mcpartvec[chargedMCPIndices.at(i)]->getEnergy());
					}
					else{
						mcp.SetPxPyPzE(0,0,0,0);
					}
					mcCharged.push_back(mcp);

				}
				for(int i=0; i< neutralMCPIndices.size(); i++){
					TLorentzVector mcp;
					if( neutralMCPIndices.at(i) != -1){
						mcp.SetPxPyPzE(_mcpartvec[neutralMCPIndices.at(i)]->getMomentum()[0],
						_mcpartvec[neutralMCPIndices.at(i)]->getMomentum()[1],
						_mcpartvec[neutralMCPIndices.at(i)]->getMomentum()[2],
						_mcpartvec[neutralMCPIndices.at(i)]->getEnergy());
					}
					else{
						mcp.SetPxPyPzE(0,0,0,0);
					}
					mcNeutral.push_back(mcp);
					
				}
				//all MCP pulls are going to be 3 parameters E theta phi

				//measgen charged
				//vector<double> temp;
				//for(int i=0; i<measNeutral.size(); i++){
				//	temp.pushback(measNeutral );
	
				//}

				

				/*	P1_k_MeasGen_pull = (1/p1meas.Perp() -  1/mcp1.Perp())/ p1meas_err[0]; 
  					P1_Theta_MeasGen_pull = ( p1meas.Theta() - mcp1.Theta())/ p1meas_err[1];
		  			P1_Phi_MeasGen_pull = getPhiResidual(p1meas.Phi(), mcp1.Phi())/p1meas_err[2] ;
	
  					P2_k_MeasGen_pull = (1/p2meas.Perp() - 1/mcp2.Perp())/ p2meas_err[0];
		  			P2_Theta_MeasGen_pull = ( p2meas.Theta() - mcp2.Theta())/ p2meas_err[1];
  					P2_Phi_MeasGen_pull = getPhiResidual(p2meas.Phi(), mcp2.Phi()) / p2meas_err[2];

					//gama gen err
					gammagen_err = getPhotonErrors(mcgamma);

		  			Gamma_E_MeasGen_pull = ( gammameas.E() - mcgamma.E() )/ gammagen_err[0];
  					Gamma_Theta_MeasGen_pull =  ( gammameas.Theta() - mcgamma.Theta())/ gammagen_err[1];
					Gamma_Phi_MeasGen_pull = getPhiResidual(gammameas.Phi(), mcgamma.Phi()) / gammagen_err[2];
			
                       
			
					P1_k_FitGen_pull = (1/p1Fit.Perp() - 1/mcp1.Perp())/ std::sqrt(cov[0]);
		 	 		P1_Theta_FitGen_pull = (p1Fit.Theta() - mcp1.Theta())/ std::sqrt(cov[10]);
					P1_Phi_FitGen_pull = getPhiResidual(p1Fit.Phi(), mcp1.Phi()) / std::sqrt(cov[20]);			

		  			P2_k_FitGen_pull = (1/p2Fit.Perp() - 1/mcp2.Perp())/ std::sqrt(cov[30]);
  					P2_Theta_FitGen_pull = (p2Fit.Theta() - mcp2.Theta())/ std::sqrt(cov[40]);
					P2_Phi_FitGen_pull = getPhiResidual(p2Fit.Phi(), mcp2.Phi()) / std::sqrt(cov[50]);

  					Gamma_E_FitGen_pull = (gammaFit.E() - mcgamma.E())/ std::sqrt(cov[60]);
		  			Gamma_Theta_FitGen_pull = (gammaFit.Theta() - mcgamma.Theta())/ std::sqrt(cov[70]);
					Gamma_Phi_FitGen_pull = getPhiResidual(gammaFit.Phi(), mcgamma.Phi()) / std::sqrt(cov[80]);	*/		


/*					if(fabs(Gamma_Phi_FitGen_pull)> 100 ||fabs( P2_Phi_FitGen_pull) > 100 || fabs(P1_Phi_FitGen_pull) >100){
						std::cout<<"crazy!! phi from fit gen"<<std::endl;
						std::cout<<"(P1,P2,Gamma) "<<P1_Phi_FitGen_pull<<" "<<P2_Phi_FitGen_pull<<" "<<Gamma_Phi_FitGen_pull<<std::endl;
						std::cout<<"residual (P1,P2,Gamma)"<<getPhiResidual(p1Fit.Phi(), mcp1.Phi())<<" "<<getPhiResidual(p2Fit.Phi(),mcp2.Phi())<<" "<<getPhiResidual(gammaFit.Phi(), mcgamma.Phi())<<std::endl;
						std::cout<<"errors (P1,P2,Gamma)"<<std::sqrt(cov[50])<<" "<<std::sqrt(cov[80])<<" "<<std::sqrt(cov[20])<<std::endl;
					
					}
					if(fabs(Gamma_Phi_MeasGen_pull)> 100 ||fabs(P2_Phi_MeasGen_pull)> 100 || fabs(P1_Phi_MeasGen_pull) >100){
						std::cout<<"crazzy!! phi from meas gen"<<std::endl;
					}*/
					//genTree->Fill();	
                        	}
 				else{
					std::cout<<"Could not locate all corresponding MC particles"<<std::endl;
				}
			}//end genfit analysis option
                }//end prob cut
	}// end cov dim check
     
	tree->Fill();
	//memory management
	pneutral.clear();
	ptrack.clear();
	pNeutralVec.clear();
	pTrackVec.clear();

	measNeutral.clear();
	measCharged.clear();
	mesaTrack.clear();
	fitNeutral.clear();
	fitCharged.clear();
	fitTrack.clear();
	
	measCharged_err.clear();
	measNeutral_err.clear();
	parentmeas_err.clear();

	fitCharged_err.clear();
	fitNeutral_err.clear();
	
	measParent_err.clear();
	fitParent_err.clear();

	//gammagen_err.clear();
		neutralJets.clear();
		TrackFO.clear();
	delete fitter;
    }//end pvector size check
	//track events for each call
	evtNo++;
	
	return;
}


