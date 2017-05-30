#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "gear/BField.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "DiTrackGammaCandidateFinder.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;
typedef CLHEP::Hep3Vector Vector3D ;

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

DiTrackGammaCandidateFinder aDiTrackGammaCandidateFinder;

//////////////////Adapted from GammaGammaCandidateFinder.cc//////////////////////////////
DiTrackGammaCandidateFinder::DiTrackGammaCandidateFinder() : marlin::Processor("DiTrackGammaCandidateFinder") {

  registerProcessorParameter( "Printing" , 
			      "Print certain messages"  ,
			      _printing,
			       (int)1 ) ;

  registerProcessorParameter("RootFile" ,
			     "Name of the output root file" ,
			     m_rootFile,
			     std::string("DiTrackGammaFitAnalysis.root"));

  
  registerProcessorParameter( "DiTrackGammaResonanceName" , 
			      "Particle decaying to x+ x- Gamma"  ,
			      _resonanceName,
			    std::string("diTrackGammaResonanceName")) ;
  
  registerProcessorParameter( "DiTrackGammaResonancePID" ,
			      "PID of Particle decaying to x+ x- Gamma" ,
			      _resonancePID,
			      (int)211) ;

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

 // collections that are from fast monte carlo benchmarking simulation
/* registerProcessorParameter("BenchmarkSimulation",
		 	    "Option to use Filtered MCParticles and Fast Reconstucted Particles",
			    _BenchmarkOption,
			    (int)0); 
 
  std::string inputPhotonCollectionName="FastReconstructedPhotons";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "FastPhotonCollection",
			   "Fast Smeared Photons from MCParticle for benchmark",
			   _inputPhotonCollectionName,
			   inputPhotonCollectionName);
  std::string inputChargedParticleCollectionName="FastReconstructedChargedParticles";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "FastChargedParticleCollection",
			   "Fast Smeared Charged Particles from MCParticle for benchmark",
			   _inputChargedParticleCollectionName,
			   inputChargedParticleCollectionName); 
/////////// */
  std::string outputParticleCollectionName = "DiTrackGammaCandidates";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "OutputParticleCollectionName" , 
			     "Output Particle Collection Name "  ,
			     _outputParticleCollectionName,
			     outputParticleCollectionName);

  registerProcessorParameter( "GammaMomentumCut" , 
			      "Minimum Momentum Cut for Photon (GeV)"  ,
			      _gammaMomentumCut,
			       (float)0.5) ;   

  registerProcessorParameter( "MinPtCut" ,
			      "Minimum Pt Charged Particle (GeV)" ,
			      _minPtCut,
			      (float)0.8) ;

  registerProcessorParameter( "DiTrackGammaResonanceMass" , 
			      "Nominal Mass of Particle Decaying to x+ x- Gamma (GeV)"  ,
			      _resonanceMass,
			       (float)0.547) ;
  
  registerProcessorParameter( "DaughterDecayMass" ,
			      "Daughter Charged Particle Mass for Mass Constraint",
			      _daughterDecayMass,
				(float)0.135) ;
  
  registerProcessorParameter("DaughterDecayPID" ,
			     "PID of charged particle pair" ,
			     _daughterPID,
			     (int)211);
 
  registerProcessorParameter("NeutralPID",
			     "PID of neutral particle",
			     _neutralPID,
			     (int)22);

  registerProcessorParameter( "MaxDeltaM" ,
				"Maximum difference between candidate mass and Parent Resonance mass (GeV)"  ,
				_dmcut,
				(float)0.07) ;

  registerProcessorParameter( "FitProbabilityCut" , 
			      "Minimum fit probability"  ,
			      _fitProbabilityCut,
			      (double)0.001);  

  registerProcessorParameter( "fitter" ,
                              "0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
                              _ifitter,
                              (int)0);

  registerProcessorParameter("FitAnalysis" ,
			      "0 = No TTree , 1 = make TTree",
			      _fitAnalysis,
			      (int)1);
  
  registerProcessorParameter("GeneratorAnalysis",
			     "0 = No TTree , 1 = make TTree",
			     _genAnalysis,
			     (int)1);

  registerProcessorParameter("MCDeltaCharged",
			     "Fractional Deviation allowed for Charged Reco Particles matched to a MC particle",
			     _mcDeltaCharged,
			    (float) 0.01);

  registerProcessorParameter("MCDeltaGamma",
			     "Fractional Deviation allowed for Gamma Reco Particles matched to a MC particle",
			     _mcDeltaGamma,
			     (float) 0.6); 			

  registerProcessorParameter("MCVertex",
			     "Allowed MC distance from interaction point",
			     _mcVertex,
			    (float) 0.0);

  registerProcessorParameter("usePhotonCalibration",
			     "use pandora reconstruction photon energy adjustments",
			     _usePhotonCalibration,
			    (int) 1);
  registerProcessorParameter("useTrackCalibration",
			     "use track error recalibration",
			     _useTrackCalibration,
			    (int) 1);
  return;

}

//===================================================================================

void DiTrackGammaCandidateFinder::init() {
  if(_printing>1)printParameters(); 

  evtNo=0;

  if(_fitAnalysis || _genAnalysis)
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

	//measured particle TLorentzVectors and errors (E,theta,phi) or (k,theta,phi)
        tree->Branch("gammameas.",&gammameas);
	tree->Branch("p1meas.",&p1meas);
	tree->Branch("p2meas.",&p2meas);
	tree->Branch("gammameas_err", "vector<double>", &gammameas_err);
	tree->Branch("p1meas_err", "vector<double>", &p1meas_err);
	tree->Branch("p2meas_err", "vector<double>", &p2meas_err);

	//fit particle TLorentzVectors and errors (E,theta,phi) or (k,theta,phi)
	tree->Branch("gammaFit.", &gammaFit);
	tree->Branch("p1Fit.", &p1Fit);
	tree->Branch("p2Fit.", &p2Fit);
	tree->Branch("gammafit_err", "vector<double>", &gammafit_err);
	tree->Branch("p1fit_err", "vector<double>", &p1fit_err);
	tree->Branch("p2fit_err", "vector<double>",&p2fit_err);

	//parent particle
	tree->Branch("parentmeas.",&parentmeas);
	tree->Branch("parentFit.", &parentFit);
	tree->Branch("parentmeas_err", "vector<double>", &parentmeas_err);
	tree->Branch("parentfit_err", "vector<double>", &parentfit_err);

	//pull distributions between fit and measured particles
	// fit-meas / sqrt( var(meas)-var(fit) )
  	tree->Branch("P1_k_FitMeas_pull", &P1_k_FitMeas_pull);
  	tree->Branch("P1_Theta_FitMeas_pull", &P1_Theta_FitMeas_pull);
  	tree->Branch("P1_Phi_FitMeas_pull", &P1_Phi_FitMeas_pull);
  	tree->Branch("P2_k_FitMeas_pull", &P2_k_FitMeas_pull);
  	tree->Branch("P2_Theta_FitMeas_pull", &P2_Theta_FitMeas_pull);
  	tree->Branch("P2_Phi_FitMeas_pull", &P2_Phi_FitMeas_pull);
  	tree->Branch("Gamma_E_FitMeas_pull", &Gamma_E_FitMeas_pull);
  	tree->Branch("Gamma_Theta_FitMeas_pull", &Gamma_Theta_FitMeas_pull);
  	tree->Branch("Gamma_Phi_FitMeas_pull", &Gamma_Phi_FitMeas_pull);
 
	//pull distributions for fit and measured parent particlei
	//uses 4vector parameterization
	//fit-meas/ sqrt(var(meas)-var(fit))
	tree->Branch("Parent_Px_FitMeas_pull", &Parent_Px_FitMeas_pull);
	tree->Branch("Parent_Py_FitMeas_pull", &Parent_Py_FitMeas_pull); 
	tree->Branch("Parent_Pz_FitMeas_pull", &Parent_Pz_FitMeas_pull);
	tree->Branch("Parent_E_FitMeas_pull" , &Parent_E_FitMeas_pull);
    }


   ////////////////////////////////////////////////////////////////////////
  //generator tree init
  //the generator analysis houses the measured particles, the monte carlo particles
  //and pull distributions between the measured/fit particles and monte carlo particles
  //only properly matched particles are added to the tree, match is governed
  //by cuts on the fractional error between the measured and mc particle
  /////////////////////////////////////////////////////////////////////
   //initialize generator tree 
   if(_genAnalysis){
	genTree = new TTree("gentree","gentree");
	genTree->SetDirectory(rootFile);
 
	 //pulls for meas-gen/sigma(meas)
  	genTree->Branch("P1_k_MeasGen_pull",&P1_k_MeasGen_pull);
  	genTree->Branch("P1_Theta_MeasGen_pull",&P1_Theta_MeasGen_pull);
  	genTree->Branch("P1_Phi_MeasGen_pull",&P1_Phi_MeasGen_pull);
  	genTree->Branch("P2_k_MeasGen_pull",&P2_k_MeasGen_pull);
  	genTree->Branch("P2_Theta_MeasGen_pull",&P2_Theta_MeasGen_pull);
  	genTree->Branch("P2_Phi_MeasGen_pull",&P2_Phi_MeasGen_pull);
  	genTree->Branch("Gamma_E_MeasGen_pull",&Gamma_E_MeasGen_pull);
  	genTree->Branch("Gamma_Theta_MeasGen_pull",&Gamma_Theta_MeasGen_pull);
  	genTree->Branch("Gamma_Phi_MeasGen_pull",&Gamma_Phi_MeasGen_pull);

	  //pulls for fit-gen/sigma(fit)
  	genTree->Branch("P1_k_FitGen_pull",&P1_k_FitGen_pull);
  	genTree->Branch("P1_Theta_FitGen_pull",&P1_Theta_FitGen_pull);
  	genTree->Branch("P1_Phi_FitGen_pull",&P1_Phi_FitGen_pull);
  	genTree->Branch("P2_k_FitGen_pull",&P2_k_FitGen_pull);
  	genTree->Branch("P2_Theta_FitGen_pull",&P2_Theta_FitGen_pull);
 	genTree->Branch("P2_Phi_FitGen_pull",&P2_Phi_FitGen_pull);
  	genTree->Branch("Gamma_E_FitGen_pull",&Gamma_E_FitGen_pull);
  	genTree->Branch("Gamma_Theta_FitGen_pull",&Gamma_Theta_FitGen_pull);
  	genTree->Branch("Gamma_Phi_FitGen_pull",&Gamma_Phi_FitGen_pull);

	//measured 4 vector stuff for easy comparison
	genTree->Branch("gammameas.", &gammameas);
	genTree->Branch("p1meas.", &p1meas);
	genTree->Branch("p2meas.", &p2meas);
	genTree->Branch("gammameas_err", "vector<double>", &gammameas_err);
	genTree->Branch("p1meas_err", "vector<double>", &p1meas_err);
	genTree->Branch("p2meas_err", "vector<double>",&p2meas_err);


	//parent particle
	genTree->Branch("parentmeas.",&parentmeas);
	genTree->Branch("parentFit.", &parentFit);
	genTree->Branch("parentmeas_err", "vector<double>", &parentmeas_err);
	genTree->Branch("parentfit_err", "vector<double>", &parentfit_err);

	
	//mc 4 vector stuff
	genTree->Branch("mcp1.", &mcp1);
	genTree->Branch("mcp2.", &mcp2);
	genTree->Branch("mcgamma.", &mcgamma);
	genTree->Branch("gammagen_err", "vector<double>", &gammagen_err);
	genTree->Branch("evtNo", &evtNo);

	//standard useful plots again
        genTree->Branch("RecoEnergy", &RecoEnergy);
	genTree->Branch("RecoMass", &RecoMass);
    } 
 
  return;
}

//===================================================================================

void DiTrackGammaCandidateFinder::processRunHeader( LCRunHeader* run) { 
  return;
}

//===================================================================================

void DiTrackGammaCandidateFinder::processEvent( LCEvent * evt ) { 

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
std::cout<<"ENTERED PROCESS EVENT"<<std::endl;
  // Access PFO collection
  bool found = this->FindPFOs(evt);
  bool found2 = this->FindTracks(evt);
  bool found3 = 0;
//std::cout<<"CONDITONAL FOR MCPARTICLE"<<std::endl;
 if(_genAnalysis){
   found3 = this->FindMCParticles(evt);
 }
// std::cout<<"conditional befor DITRACK CALL"<<std::endl;
  if((found && found2) || (found && found && _genAnalysis && found3)){
    if(_printing>1)std::cout << "Analysis : " << _resonanceName << std::endl; 
    this->FindDiTrackGammaCandidates(recparcol);
  }
  
  
  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  return;
  
}
//===================================================================================
void DiTrackGammaCandidateFinder::end(){
 if(_fitAnalysis){
 	rootFile->cd();
 	tree->Write();
        genTree->Write(); 
 }
  return;
}
//===================================================================================
bool DiTrackGammaCandidateFinder::FindPFOs( LCEvent* evt ) {

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
bool DiTrackGammaCandidateFinder::FindTracks( LCEvent* evt ) {

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
bool DiTrackGammaCandidateFinder::FindMCParticles( LCEvent* evt ){
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
//===================================================================================
//this function adjusts phi from [-pi,pi] -> [0,2pi] and returns the difference between the two arguments
double DiTrackGammaCandidateFinder::getPhiResidual(double phi1, double phi2){
	//since phi goes -pi,pi a method for finding this quantity is useful
	//always have input phi1 - phi2
	double residual= 0.0;
	double phi_1=0.0;
	double phi_2=0.0;

	phi1 < 0 ? phi_1 = 2*M_PI+phi1 : phi_1=phi1;
	phi2 < 0 ? phi_2 = 2*M_PI+phi2 : phi_2=phi2;

	residual = phi_1 - phi_2;
	if(fabs(residual) > M_PI && residual > 0.0) residual =residual - 2*M_PI;
	if(fabs(residual) > M_PI && residual < 0.0) residual = 2*M_PI + residual;
	return residual;
}
//===================================================================================
//this function expects a global MCParticle vector of all the MCParticles from the current event, the index returned corresponds to that global MCparticle vector
//in the future a better implementation might return  MCParticle* instead, not sure which would be more versatile
int DiTrackGammaCandidateFinder::getCorrespondingMCParticleIndex(TLorentzVector vReco, int recoCharge, int recoPdg){
//note so pdgs may be mixed up because +pdgcode implied +charge, which isnt true for elctron, so disabled pdg matching and just used E,theta,phi
//it would be better in the future to look at signed curvature as a matching condition	
	//each vector index holds a score which is the sum of residuals of the MCparticle @ index i with vReco
	//returns the lowest score, which is the most closely matching index of score and of MCparticle
	//generate score based on (E,theta,phi), (k,theta,phi)
	
	//when searcing for the particle, store the best matching particle, its score, 
	//and index
	//also (temporary) look at distance from interaction point, all events considered
	//are created at interaction point
	TLorentzVector lowestMCvec;
	double pScore= 0.0;
	double lowestScore = -1;
	double lowestIndex = -1;
	double lowestMCdistance = -1;
	double MCDistanceFromIP = -1;

	double _mcDelta;
        //mc delta is a % (fractional error) based on the sum of fractional errors between the parameters (E,theta,phi)
        //the % error is different for charged and neutral particles
	if(recoPdg == 22){
		_mcDelta = _mcDeltaGamma;
	}
	else{
		_mcDelta = _mcDeltaCharged;
	}
	//iterate through the list of mc particles for this event, calculate scores, and find the best match
	for(unsigned int i=0; i<_mcpartvec.size(); i++){ 
		//match charge and pdg
	//	if((int)_mcpartvec[i]->getCharge()==recoCharge && _mcpartvec[i]->getPDG() == recoPdg ){
		if(_mcpartvec[i]->getPDG() == recoPdg ){
			MCDistanceFromIP = std::sqrt(_mcpartvec[i]->getVertex()[0]*_mcpartvec[i]->getVertex()[0] + _mcpartvec[i]->getVertex()[1]*_mcpartvec[i]->getVertex()[1] + _mcpartvec[i]->getVertex()[2]*_mcpartvec[i]->getVertex()[2]  );
			//check vertex cut
			if(MCDistanceFromIP <= _mcVertex){
				//normalize all scores
				TLorentzVector mcp(_mcpartvec[i]->getMomentum()[0], _mcpartvec[i]->getMomentum()[1], _mcpartvec[i]->getMomentum()[2], _mcpartvec[i]->getEnergy());
				//add energy pointes to score
				//try not considering E for score with photons because of sqrt(E) error scaling
				if(recoPdg != 22){
				pScore += fabs(mcp.E()-vReco.E())/mcp.E();
				}
				//add theta points to score
				pScore += fabs(mcp.Theta()-vReco.Theta())/mcp.Theta();

				//add phi points to score
				if(mcp.Phi()>=0){ //phi goes from -pi to pi so extra helper functions are added for phi residuals and 2pi is added where necessary
					pScore += fabs( getPhiResidual(mcp.Phi(), vReco.Phi())) / mcp.Phi();
				} 
				else{
					pScore += fabs( getPhiResidual(mcp.Phi(), vReco.Phi()))/ (2*M_PI + mcp.Phi());
				}
				//if this particle is the best match and passes the cuts update it as the closest match
				if( (lowestScore == -1 && (pScore < _mcDelta)) || ((pScore < lowestScore) && lowestScore <_mcDelta) ){
					lowestMCvec = mcp;
					lowestScore = pScore;
					lowestIndex = i;
					lowestMCdistance = MCDistanceFromIP;
				}
			}		
		}
	}
	//if a particle was matched print relevant information about the matching
	if(lowestIndex != -1){
        	if(_printing>3){
			std::cout<<"MC Match: "<<std::endl;
			std::cout<<"Reco (E,theta,phi, Charge) "<<vReco.E()<<", "<<vReco.Theta()<<", "<<vReco.Phi()<<" "<<recoCharge<<std::endl;
			std::cout<<"MC   (E,theta,phi, Charge) "<<lowestMCvec.E()<<", "<<lowestMCvec.Theta()<<", "<<lowestMCvec.Phi()<<" "<<_mcpartvec[lowestIndex]->getCharge() <<std::endl;
			std::cout<<"MC Distance from IP "<<lowestMCdistance<<std::endl;
			std::cout<<"Match Fractional Error "<< lowestScore <<std::endl;
		}
	}
	else{
		if(_printing>3){
		std::cout<<"Particle not matched "<<std::endl;
		std::cout<<"Reco (E,theta,phi) Charge "<<vReco.E()<<", "<<vReco.Theta()<<", "<<vReco.Phi()<<" "<<recoCharge<<std::endl;
		}
	}
	return lowestIndex;
	
}
//adjusts photon energy
void DiTrackGammaCandidateFinder::calibratePhoton(TLorentzVector& v){
	//reconstruction sets photon energy too high, reduce by 10%
	double cFactor =0.91;
	v.SetPxPyPzE(cFactor*v.Px(),cFactor*v.Py(),cFactor*v.Pz(),cFactor*v.E());
}
//adjusts the curvature errors
void DiTrackGammaCandidateFinder::calibratePionError(std::vector<double>& errors){
	//track curvature errors are 30% too low
	double cFactor = 0.285;
	errors[0] = errors[0]+cFactor*errors[0];	
}
//===================================================================================
//returns error array with parameterization {dE, dTheta, dPhi}
//error for photon energy is sigma/E = 0.18/sqrt(E)
//error on photon angle is sigma/E = 0.001/sqrt(E)
//i need to get errors from the recoparticle error matrix probably not make my own?
std::vector<double> DiTrackGammaCandidateFinder::getPhotonErrors(TLorentzVector pgamma){
	std::vector<double> errors;/// = new double[3];
	 errors.push_back(0.18*std::sqrt(pgamma.E()) );
	// errors.push_back( 0.001/std::sqrt(pgamma.E()) );
	// errors.push_back( 0.001/std::sqrt(pgamma.E()) );
	errors.push_back(0.001/std::sqrt(pgamma.E()) );
	errors.push_back(0.0013/std::sqrt(pgamma.E()) );
	return errors;

}
//===================================================================================
//returns error array with parameterization {dk, dTheta, dPhi} with the standard parameterization used in the code
//errors come from track record with (d0,z0,omega,tanLambda,theta, lamda)
std::vector<double> DiTrackGammaCandidateFinder::getChargedParticleErrors(TLorentzVector pcharged, Track* ptrk){
	std::vector<double> errors; // = new double[3];
		/*dk = 1/pt*dOmega/Omega */
	errors.push_back( fabs( 1/pcharged.Perp()*(std::sqrt(ptrk->getCovMatrix()[5])/ptrk->getOmega()) ) );
		//dk is dOmega going to change 1/pt = k to k = 1/R
	//	errors.push_back(ptrk->getCovMatrix()[5]);
		/*dTheta=cos^2(lamda)*dtanLamda*/
	errors.push_back( std::sqrt(ptrk->getCovMatrix()[14])/(ptrk->getTanLambda()*ptrk->getTanLambda() + 1) );
	errors.push_back( std::sqrt(ptrk->getCovMatrix()[2]) );

	//calibration TODO::add calibration parameter option////////////////////////////////////////////////////////
	if(_daughterPID == 211 && _useTrackCalibration){
		calibratePionError(errors);
	}
	return errors;
}
//===================================================================================
void DiTrackGammaCandidateFinder::PrintCov(FloatVec cov, int dim){
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
void DiTrackGammaCandidateFinder::PrintCov(double* cov, int dim){
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
void DiTrackGammaCandidateFinder::setFitErrors(double* cov){

	gammafit_err.push_back(  std::sqrt(cov[60]) );
	gammafit_err.push_back(  std::sqrt(cov[70]) );
	gammafit_err.push_back(  std::sqrt(cov[80]) );
	
	p1fit_err.push_back(  std::sqrt(cov[0]) );
	p1fit_err.push_back(  std::sqrt(cov[10]) );
	p1fit_err.push_back(  std::sqrt(cov[20]) );
	
	p2fit_err.push_back(  std::sqrt(cov[30]) );
	p2fit_err.push_back(  std::sqrt(cov[40]) );
	p2fit_err.push_back(  std::sqrt(cov[50]) );
 

}
//===================================================================================
//manually constructs the 9x9 covariance matrix for parent particle, meant to be transformed into 4x4 4vector form
//***********This is to be used after the global measured errors vectors are  populated**************************
double* DiTrackGammaCandidateFinder::ConstructParentMeasCovMatrix(){
//concat vectors	
	std::vector<double> params;
	params.insert(params.end(), p1meas_err.begin(), p1meas_err.end());
	params.insert(params.end(), p2meas_err.begin(), p2meas_err.end());
	params.insert(params.end(), gammameas_err.begin(), gammameas_err.end());
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
void DiTrackGammaCandidateFinder::setParentErrors(FloatVec meascov, FloatVec fitcov){

	parentmeas_err.push_back( std::sqrt(meascov[0]) );
	parentmeas_err.push_back( std::sqrt(meascov[5]) );
	parentmeas_err.push_back( std::sqrt(meascov[10]) );
	parentmeas_err.push_back( std::sqrt(meascov[15]) );
	
	parentfit_err.push_back( std::sqrt(fitcov[0]) );
	parentfit_err.push_back( std::sqrt(fitcov[5]) );
	parentfit_err.push_back( std::sqrt(fitcov[10]) );
	parentfit_err.push_back( std::sqrt(fitcov[15]) );
}
//===================================================================================
FloatVec DiTrackGammaCandidateFinder::ConstructCovMatrix(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, double*  cov){
        
          double Energy; //parent energy
          double Mom[3]; // parent momentum
      // The 4-vector of the fitted gamma-gamma is the sum of the two fitted 4-vectors.

      // Get the individual fitted photon 4-vectors - needed to calculate the four-momentum covariance matrix of the gamma-gamma system
          double Mom1[3],Mom2[3], Mom3[3]; //momentum of +,-,gamma
          double E1,E2,E3; // energy of +,-,gamma
          double pT1,pT2, pT3; // transverse momentum of +,-,gamma
          double P1,P2,P3;

          E1 = p1.E();
          Mom1[0] = p1.Px();
          Mom1[1] = p1.Py();
          Mom1[2] = p1.Pz();
          pT1 = sqrt(Mom1[0]*Mom1[0] + Mom1[1]*Mom1[1]);
	  P1 = sqrt( Mom1[0]*Mom1[0] + Mom1[1]*Mom1[1] + Mom1[2]*Mom1[2]);

          E2 = p2.E();
          Mom2[0] = p2.Px();
          Mom2[1] = p2.Py();
          Mom2[2] = p2.Pz();
          pT2 = sqrt(Mom2[0]*Mom2[0] + Mom2[1]*Mom2[1]);
	  P2 = sqrt(Mom2[0]*Mom2[0] + Mom2[1]*Mom2[1] + Mom2[2]*Mom2[2]);

	  E3 = p3.E();
	  Mom3[0] = p3.Px();
	  Mom3[1] = p3.Py();
	  Mom3[2] = p3.Pz();
	  pT3 = sqrt(Mom3[0]*Mom3[0] + Mom3[1]*Mom3[1]);
          P3 = sqrt(Mom3[0]*Mom3[0] + Mom3[1]*Mom3[1] + Mom3[2]*Mom3[2]);

          Energy = E1 + E2 + E3;
          Mom[0] = Mom1[0] + Mom2[0] + Mom3[0];
          Mom[1] = Mom1[1] + Mom2[1] + Mom3[1];
          Mom[2] = Mom1[2] + Mom2[2] + Mom3[2];

       //   FloatVec covP(10, 0.0);
            FloatVec covP(16,0.0);
      //  (px,py,pz,E)

      //  Follow for now implementation in FourMomentumCovMat based on TMatrixD

         // if(cov_dim == 9 ){
         if(true){
      // Note some cases of cov_dim != 6 have been seen ...
             const int nrows      = 9; // n rows jacobian
             const int ncolumns   = 4; // n columns jacobian //these look backwards

      //  Transformation matrix D (6*4) between (E1, theta1, phi1, E2, theta2, phi2) and (px, py, pz, E)
            // double jacobian_by_columns[nrows*ncolumns] = {
            //    Mom1[0]/E1, Mom1[0]*Mom1[2]/pT1, -Mom1[1], Mom2[0]/E2, Mom2[0]*Mom2[2]/pT2, -Mom2[1],
            //    Mom1[1]/E1, Mom1[1]*Mom1[2]/pT1,  Mom1[0], Mom2[1]/E2, Mom2[1]*Mom2[2]/pT2,  Mom2[0],
             //   Mom1[2]/E1,        -pT1        ,    0.0  , Mom2[2]/E2,        -pT2        ,     0.0 ,
              //     1.0    ,         0.0        ,    0.0  ,    1.0    ,         0.0        ,     0.0   };
	 //Transformation matrix D (9*4) between (E1,theta1,phi1,E2,theta2,phi2,E3,theta3,phi3) and (px,py,pz,E)
	 double jacobian_by_columns[nrows*ncolumns] = {
-pT1*Mom1[0], Mom1[0]*Mom1[2]/pT1, -Mom1[1], -pT2*Mom2[0], Mom2[0]*Mom2[2]/pT2, -Mom2[1], Mom3[0]/E3, Mom3[0]*Mom3[2]/pT3, -Mom3[1],
-pT1*Mom1[1], Mom1[1]*Mom1[2]/pT1,  Mom1[0], -pT2*Mom2[1], Mom2[1]*Mom2[2]/pT2,  Mom2[0], Mom3[1]/E3, Mom3[1]*Mom3[2]/pT3,  Mom3[0],
-pT1*P1     ,        -pT1,              0,   -pT2*P2     ,         -pT2,            0,    Mom3[2]/E3,        -pT3,            0,
-pT1*P1*P1/E1,         0,               0,   -pT2*P2*P2/E2,          0,             0,         1,              0,             0   };


//test:: try particle order 312 instead of 123 (photon first in jacobian)
//	double jacobian_by_columns[nrows*ncolumns] = {
// Mom3[0]/E3, Mom3[0]*Mom3[2]/pT3, -Mom3[1], -pT1*Mom1[0], Mom1[0]*Mom1[2]/pT1, -Mom1[1], -pT2*Mom2[0], Mom2[0]*Mom2[2]/pT2, -Mom2[1],
// Mom3[1]/E3, Mom3[1]*Mom3[2]/pT3,  Mom3[0], -pT1*Mom1[1], Mom1[1]*Mom1[2]/pT1,  Mom1[0], -pT2*Mom2[1], Mom2[1]*Mom2[2]/pT2,  Mom2[0],
// Mom3[2]/E3,        -pT3,            0, -pT1*P1     ,        -pT1,              0,   -pT2*P2     ,         -pT2,            0, 
 //    1,              0,             0,    -pT1*P1*P1/E1,         0,               0,   -pT2*P2*P2/E2,          0,             0 };

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
             }

        /*     covP[0] = Covmatrix(0,0); // px-px
             covP[1] = Covmatrix(0,1); // px-py
             covP[2] = Covmatrix(1,1); // py-py
             covP[3] = Covmatrix(0,2); // px-pz
             covP[4] = Covmatrix(1,2); // py-pz
             covP[5] = Covmatrix(2,2); // pz-pz
             covP[6] = Covmatrix(0,3); // px-E
             covP[7] = Covmatrix(1,3); // py-E
             covP[8] = Covmatrix(2,3); // pz-E
             covP[9] = Covmatrix(3,3); // E-E
	*/
          }
          else{
              // if(_printing>1)std::cout << " Fit has no valid covariance information (cov_dim = " << cov_dim << " )" << std::endl;              
          }

	return covP;
} 
//===================================================================================
//Current implementation only allows one type of fitter because using the base class is not working correctly
//BaseFitter* DiTrackGammaCandidateFinder::setUpFit(TLorentzVector gamma, TLorentzVector p1, TLorentzVector p2, Track* p1Track, Track* p2Track ){
OPALFitterGSL* DiTrackGammaCandidateFinder::setUpFit(TLorentzVector gamma, TLorentzVector p1, TLorentzVector p2, Track* p1Track, Track* p2Track){

	//set up mass constraint
	MassConstraint mc( (double)_resonanceMass );

	//construct fit objects for gamma, and 2 charged particles
	std::vector<double> PhotonErrors = getPhotonErrors(gamma);
//	JetFitObject*	
	 gammaJet = new JetFitObject(gamma.E(), gamma.Theta(), gamma.Phi(), PhotonErrors[0], PhotonErrors[1], PhotonErrors[2]);
	gammaJet->setName("Photon");					

	std::vector<double> P1Errors = getChargedParticleErrors(p1,p1Track);
//	LeptonFitObject* 
	part1 = new LeptonFitObject(1/p1.Perp(),p1.Theta(),p1.Phi(),P1Errors[0],P1Errors[1],P1Errors[2],_daughterDecayMass);
	part1->setName("Particle+");
					
	std::vector<double> P2Errors = getChargedParticleErrors(p2,p2Track);
//	LeptonFitObject* 
	part2 = new LeptonFitObject(1/p2.Perp(),p2.Theta(),p2.Phi(),P2Errors[0],P2Errors[1],P2Errors[2],_daughterDecayMass);
	part2->setName("Particle-");
	if(_printing>4)std::cout <<" Measured Quantities" <<std::endl;						
	if(_printing>4)std::cout <<" Track1 particle (k,theta,phi): " << 1/p1.Perp() << " " << p1.Theta() << " " << p1.Phi() << std::endl;
        if(_printing>4)std::cout <<" Errors (dk, dtheta, dphi): " << P1Errors[0]<< " " << P1Errors[1] << " " <<P1Errors[2] <<std::endl;				
	if(_printing>4)std::cout <<" Track2 particle (k,theta,phi): " << 1/p2.Perp() << " " << p2.Theta() << " " << p2.Phi() << std::endl;
        if(_printing>4)std::cout <<" Errors (dk,dtheta,dphi): " << P2Errors[0]<< " " << P2Errors[1] << " " << P2Errors[2] << std::endl;
	if(_printing>4)std::cout <<" Photon (E,theta,phi): " << gamma.E() << " " << gamma.Theta() << " " << gamma.Phi() << std::endl;
	if(_printing>4)std::cout <<" Errors (dE,dtheta,dphi): " << PhotonErrors[0] << " " << PhotonErrors[1] << " " << PhotonErrors[2] << std::endl;

	//clean up memory with errors that are not going to be used again in this scope
	 PhotonErrors.clear();
	 P1Errors.clear();
	 P2Errors.clear();
					
	//add fitobjects to mass constraint
//	mc.addToFOList(*gammaJet);
        mc.addToFOList(*part1);
	mc.addToFOList(*part2);
	mc.addToFOList(*gammaJet);

	//declare fitter
	OPALFitterGSL * fitter = new OPALFitterGSL();
	
	//add constraints and fit objects to fitter
//	fitter->addFitObject(gammaJet);
	fitter->addFitObject(part1);
	fitter->addFitObject(part2);
	fitter->addFitObject(gammaJet);
	fitter->addConstraint(mc);
	
	//do the fit (this also returns a fit probability)
	fitter->fit();

	//TODO:: look up the function that returns fit prob, dont fit twice if printing	
//	if(_printing>4)std::cout <<" Fit Probability: "<< fitter->fit()<<std::endl;

	
	return fitter;

}
//===================================================================================
void DiTrackGammaCandidateFinder::FindDiTrackGammaCandidates(LCCollectionVec * recparcol) {

  if(_printing>1)std::cout << "FindDiTrackGammaCandidates : (nPFOs = " << _pfovec.size() << " )" << std::endl; 
  if(_printing>1)std::cout << "FindDiTrackGammaCandidates : evtNo = " << evtNo << std::endl;
  // Look for candidates
  std::vector<TLorentzVector>pplus;
  std::vector<TLorentzVector>pminus;
  std::vector<TLorentzVector>pgamma;

  std::vector<ReconstructedParticle*>pGammaPfoVec;
  std::vector<Track*>pPlusTrackVec;
  std::vector<Track*>pMinusTrackVec;
//photon loop, load collection photons into local variables
  for(unsigned int i=0;i<_pfovec.size();i++){
//get type should be changed to getPDG? will this mess up other collections?
	if(_pfovec[i]->getType() == 22 && _pfovec[i]->getEnergy()>_gammaMomentumCut){//gamma
		 if(_printing>1)std::cout << "FindDiTrackGammaCandidates : Photon " << _pfovec[i]->getType() << std::endl;
		TLorentzVector pGamma(_pfovec[i]->getMomentum()[0], _pfovec[i]->getMomentum()[1], _pfovec[i]->getMomentum()[2], _pfovec[i]->getEnergy() );
		if(_printing>1)std::cout << "g candidate (px,py,pz,E) "<<pGamma.Px()<<" "<<pGamma.Py()<<" "<<pGamma.Pz()<<" "<<pGamma.E()<<std::endl;
	//	if(_printing>1)std::cout << "g candidate (E, theta, phi) "<< pGamma.E()<<" "<<pgamma.Theta()<<" "<<pgamma.Phi()<<std::endl;
		if(_usePhotonCalibration){//placeholder for photon calibration parameter
			calibratePhoton(pGamma);
			if(_printing>1)std::cout <<"g candidate calibration (px,py,pz,E) " <<pGamma.Px()<<" "<<pGamma.Py()<<" "<<pGamma.Pz()<<" "<<pGamma.E()<<std::endl;
		//	if(_printing>1)std::cout <<"g candidate calibration (E,theta,phi) "<<pGamma.E()<<" "<<pgamma.Theta()<<" "<<pgamma.Phi()<<std::endl;
		}
		pgamma.push_back(pGamma);
		pGammaPfoVec.push_back(_pfovec[i]);
	}
  }//end pfovec loop
//start track loop



  const double c = 2.99792458e8; // m*s^-1
  const double B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();;          
  const double mm2m = 1e-3;
  const double eV2GeV = 1e-9;
  const double eB = B*c*mm2m*eV2GeV;
  //const double pionMass = 0.13497;

 for(unsigned int i=0; i<_trackvec.size();i++){

		
		//TODO: add momentum cuts
		if(eB/fabs(_trackvec[i]->getOmega()) < _minPtCut) continue;
		//there must be an error in here somewhere in these calculations is it E or the Px,py,pz?	
		double cosLambda = 1 / std::sqrt(1 + _trackvec[i]->getTanLambda()*_trackvec[i]->getTanLambda() );
		double P = (eB/fabs(_trackvec[i]->getOmega()))/cosLambda;
		double sinLambda = _trackvec[i]->getTanLambda()*cosLambda;
		double cosPhi = cos(_trackvec[i]->getPhi());
		double sinPhi = sin(_trackvec[i]->getPhi());
		double px = P*cosLambda*cosPhi;
		double py = P*cosLambda*sinPhi;
		double pz = P*sinLambda;
		double E = std::sqrt( _daughterDecayMass*_daughterDecayMass + px*px + py*py + pz*pz );
		

		if(_trackvec[i]->getOmega() > 0.0){
			pPlusTrackVec.push_back(_trackvec[i]);
			TLorentzVector Plus(px,py,pz,E);
			if(_printing>1)std::cout << "+ candidate "<< i <<" (px,py,pz,E)  "<< px << " " << py << " " << pz << " " << E << std::endl;
			if(_printing>1)std::cout << "+ candidate "<< i << " (d0,phi,k,z0,tanLambda) " <<_trackvec[i]->getD0()<<" "<<_trackvec[i]->getPhi()<< " "<<_trackvec[i]->getOmega()<< " "<< _trackvec[i]->getZ0()<< " "<<_trackvec[i]->getTanLambda()<<std::endl;
			pplus.push_back(Plus);		
		}
		if(_trackvec[i]->getOmega() < 0.0){
			pMinusTrackVec.push_back(_trackvec[i]);
			TLorentzVector Minus(px,py,pz,E);
			if(_printing>1)std::cout << "- candidate "<< i << " (px,py,pz,E) "<< px << " " << py << " " << pz << " " << E << std::endl;
			if(_printing>1)std::cout << "- candidate "<< i << " (d0,phi,k,z0,tanLambda) " <<_trackvec[i]->getD0()<<" "<<_trackvec[i]->getPhi()<< " "<<_trackvec[i]->getOmega()<< " "<< _trackvec[i]->getZ0()<< " "<<_trackvec[i]->getTanLambda()<<std::endl;

			pminus.push_back(Minus);
		}
	
}
  	
  
  if(_printing>1)std::cout << "FindDiTrackGammaCandidates : (nphotons = " << pGammaPfoVec.size() << " " << pgamma.size() << " )" << std::endl;  
  if(_printing>1)std::cout << "FindDiTrackGammaCandidates : (npi+ = " << pPlusTrackVec.size() << " " << pplus.size() << " )" << std::endl; 
  if(_printing>1)std::cout << "FindDiTrackGammaCandidates : (npi- = " << pMinusTrackVec.size() << " " << pminus.size() << " )" << std::endl; 
  // loop over pairs of candidate photons and keep all combinations consistent with the given parent mass assumption

    if( pgamma.size() >= 1 && pplus.size() >= 1 && pminus.size() >= 1){

    

//store the best fit prob i,j,k then after iterating through all permutations recompute fit and compute analysis vars
	double fitprobmax = -1;
	double candidate_ijk[3] = {-1,-1,-1};
	double fitprob = -1;

	//iterate and fit all particle combinations, save the best one
	for(unsigned int i=0; i<pgamma.size();i++){
		for(unsigned int j=0; j<pplus.size(); j++){
			for(unsigned int k=0; k<pminus.size(); k++){

				TLorentzVector parent = pgamma[i] + pplus[j] + pminus[k];
				
				 if(_printing>2)std::cout << "Testing for " << _resonanceName << " " << i << " " << j << " " << k << "  M = " << parent.M() << std::endl; 
				if( fabs(parent.M() - _resonanceMass) < _dmcut ){ //particle combination passes, perform mass constrained fits
			
					
					

					OPALFitterGSL* fitter = setUpFit(pgamma[i], pplus[j], pminus[k], pPlusTrackVec[j], pMinusTrackVec[k]);
			  
					fitprob = (double) fitter->getProbability();
					std::cout<<"fit prob "<<fitprob<<std::endl;     	
					if(fitprobmax == -1 && fitprob > _fitProbabilityCut){
						fitprobmax = fitprob;
						candidate_ijk[0] = i;
						candidate_ijk[1] = j;
						candidate_ijk[2] = k;
					}
					if(fitprob > fitprobmax && fitprob > _fitProbabilityCut){
						fitprobmax = fitprob;
						candidate_ijk[0] = i;
						candidate_ijk[1] = j;
						candidate_ijk[2] = k;
					}
				// memory management, since fit objects are stored globally delete them after each fit, only store best candidate indices and refit later
					delete gammaJet;
					delete part1;
					delete part2;
					delete fitter;	
				}
			}
		}
	}
     
	std::cout<<"loop completed"<<std::endl;
	//recompute params for best fit probability
//if nothing makes it for this event just return
	if(candidate_ijk[0]==-1 || candidate_ijk[1]==-1 || candidate_ijk[2]==-1 ){
		evtNo++; 
		return;
	}
     
	

	OPALFitterGSL* fitter = setUpFit(pgamma[candidate_ijk[0]], pplus[candidate_ijk[1]], pminus[candidate_ijk[2]], pPlusTrackVec[candidate_ijk[1]], pMinusTrackVec[candidate_ijk[2]]);
	
					
	int cov_dim;
	double * cov = fitter->getGlobalCovarianceMatrix(cov_dim);  //9x9
	//construct measured cov matrix
//	double * mcov = ConstructParentMeasCovMatrix();
	//if cov matrix dimension is not correct return, fit did not converge
	if(cov_dim == 0){
		evtNo++;
		return;
	}
		
	gammaFit.SetPxPyPzE(gammaJet->getPx(),gammaJet->getPy(),gammaJet->getPz(),gammaJet->getE());
	p1Fit.SetPxPyPzE(part1->getPx(), part1->getPy(), part1->getPz(), part1->getE());
	p2Fit.SetPxPyPzE(part2->getPx(), part2->getPy(), part2->getPz(), part2->getE());

        parentFit = gammaFit + p1Fit + p2Fit;

	//sets gammafit_err
	setFitErrors(cov);

	//put candidates into local LorentzVector for ease of coding and readability
	p1meas = pplus[candidate_ijk[1]];
	p2meas = pminus[candidate_ijk[2]];
	gammameas = pgamma[candidate_ijk[0]];
        parentmeas = p1meas + p2meas + gammameas;
	//also have error arrays ready locally for readability
	//Photon Errors go {dE,dTheta,dPhi}
	//ChargedParticle errors go {dk, dTheta, dPhi}
	p1meas_err = getChargedParticleErrors( p1meas, pPlusTrackVec[candidate_ijk[1]]);
	p2meas_err = getChargedParticleErrors( p2meas, pMinusTrackVec[candidate_ijk[2]]);
	gammameas_err = getPhotonErrors( gammameas );
	
	//construct measured cov matrix must be called after measured global vectors are populated
	double * mcov = ConstructParentMeasCovMatrix();

	//populate measured  parent errors from covariance matrix
	FloatVec fitParentCov =  ConstructCovMatrix(p1Fit,p2Fit,gammaFit,cov);  
	FloatVec measParentCov = ConstructCovMatrix(p1meas,p2meas,gammameas,mcov);
	setParentErrors(measParentCov, fitParentCov);
	

	if(_printing>3){
		std::cout <<"Optimal Candidate Fit Quantities-------------"<<std::endl;
		std::cout << "Constrained fit results  RC: " << fitter->getError() << std::endl;
		std::cout << "Measured mass = " << parentmeas.M() << std::endl;
		std::cout << "No. of iterations " << fitter->getIterations() << std::endl;
		std::cout << "Fit probability = " << fitprob << std::endl;
		std::cout << "Covariance matrix dimension " << cov_dim << std::endl;
		std::cout << "Candidates Post Fit: "<< std::endl;
		std::cout << "track1 particle (k,theta,phi): "<< 1/p1Fit.Perp() <<" " << p1Fit.Theta() << " "<< p1Fit.Phi()<<std::endl;
		std::cout << "track2 particle (k,theta,phi): "<< 1/p2Fit.Perp() <<" " << p2Fit.Theta() << " "<< p2Fit.Phi()<<std::endl;
		std::cout << "gamma (E,theta, phi): "<< gammaFit.E()<<" " <<gammaFit.Theta()<< " "<< " " <<gammaFit.Phi()<<std::endl;
		std::cout << "Gamma residual (fit-meas) (dE,dTheta,dPhi) " << gammaFit.E()-gammameas.E() << " " << gammaFit.Theta()-gammameas.Theta() << " " << gammaFit.Phi()-gammameas.Phi() << std::endl;
		std::cout << "DEBUG VALS (fit meas) " << gammaFit.E() <<" "<< gammameas.E() << " | " << gammaFit.Theta() << " " << gammameas.Theta()<< " | "<< gammaFit.Phi()<<" "<<gammameas.Phi()<<std::endl;
		std::cout << "P1 residual (fit-meas) (dk,dTheta,dPhi) " << 1/p1Fit.Perp()-1/p1meas.Perp() << " " << p1Fit.Theta()-p1meas.Theta() << " " << p1Fit.Phi()-p1meas.Phi() << std::endl;
		std::cout << "DEBUG VALS (fit meas) " << 1/p1Fit.Perp() << " " << 1/p1meas.Perp() << " | " << p1Fit.Theta() << " "<<p1meas.Theta()<<" | "<<p1Fit.Phi()<<" "<<p1meas.Phi()<<std::endl;
		std::cout << "P2 residual (fit-meas) (dk,dTheta,dPhi) " << 1/p2Fit.Perp()-1/p2meas.Perp() << " " << p2Fit.Theta()-p2meas.Theta() << " " << p2Fit.Phi()-p2meas.Phi() << std::endl;
		std::cout << "DEBUG VALS (fit meas) " << 1/p2Fit.Perp() <<" "<< 1/p2meas.Perp() <<" | "<<p2Fit.Theta() <<" "<<p2meas.Theta()<<" | "<<p2Fit.Phi()<<" "<<p2meas.Phi() << std::endl;	
		std::cout << "Parent residual (particle fit sum - particle meas sum) (dE,dTheta,dPhi) " << parentFit.E()-parentmeas.E() << " " << parentFit.Theta()-parentmeas.Theta() << " " << getPhiResidual(parentFit.Phi(),parentmeas.Phi()) << std::endl;
		std::cout << "Parent Measured Energy "<< parentmeas.E() <<" Parent Fit Energy "<< parentFit.E()<<std::endl;
		std::cout <<"----------------------"<<std::endl;

		//fit cov matrix is in the order in which fit objects are added, in a column-wise array. Set up does Photon,particle1+, particle2-. So, the diagonal for variances of the 9 parameters are as follows:
		//Particle:  Photon		P1		P2
		//Parameter: dE,dTheta,dPhi	k,dTheta,dPhi	k,dTheta,dPhi
		//Indices:   60,70,80		0,10,20 	30,40,50
//		for ( int i=0; i<cov_dim*cov_dim; i=( i+cov_dim+1) ){
//			std::cout << "Covariance matrix element sqrt(diagonal) " << i << " fit:  " << std::sqrt(cov[i]) << " measured: "<< std::sqrt(mcov[i]) <<  std::endl;                   
//		}
		std::cout<<"Measured Covariance Matrix 9x9: "<<std::endl;
		PrintCov(mcov,9);
		std::cout<<"Fit Covariance Matrix 9x9: "<<std::endl;
		PrintCov(cov,9);
		std::cout<<"Measured Parent 4vector 4x4 Covariance Matrix: "<<std::endl;
		PrintCov(measParentCov,4);
		std::cout<<"Fit Parent 4vector 4x4 Covariance Matrix: "<<std::endl;
		PrintCov(fitParentCov,4);					     
	}//end optimal fit printing
								
			              
	if(cov_dim > 0){//condition required because there is not always convergence in the fitting	
		if(fitprob> _fitProbabilityCut){
			//TODO:: add qualifying fits to collection that are less than optimal
			ReconstructedParticleImpl * recoPart = new ReconstructedParticleImpl();
			recoPart->addParticle(pGammaPfoVec[candidate_ijk[0]]);
			recoPart->addTrack(pPlusTrackVec[candidate_ijk[1]]);
			recoPart->addTrack(pMinusTrackVec[candidate_ijk[2]]);

			recoPart->setEnergy(parentFit.E());
				
			double mom[3] = {parentFit.Px(),parentFit.Py(),parentFit.Pz() };
			recoPart->setMomentum(mom);
			recoPart->setMass( _resonanceMass );
			recoPart->setCharge(0.0);
//			RETURN transformed cov matrix here
			recoPart->setCovMatrix( fitParentCov );// this call needs fancy cov matrix
			recoPart->setGoodnessOfPID((float)fitprob);

			recoPart->setType( _resonancePID ); // make some selection for pid stuff, or add parameter to xm
			
			//add best fit candidate to output collections
			if(_printing>1)std::cout << "candidate passes fit probability cut - KEEP" << std::endl;
			recparcol->addElement( recoPart );

			
			//add stuff to regular ttree
			if(_fitAnalysis){

 				RecoEnergy=parentmeas.E();
				FitEnergy=parentFit.E();
				RecoMass=parentmeas.M();
				FitProbability=fitprob; 
				P1_k_FitMeas_pull= (1/p1Fit.Perp() -1/p1meas.Perp())/ std::sqrt( p1meas_err[0]*p1meas_err[0] - cov[0] )  ;
				P1_Theta_FitMeas_pull= (p1Fit.Theta() - p1meas.Theta())/ std::sqrt(p1meas_err[1]*p1meas_err[1] - cov[10]) ;	
				P1_Phi_FitMeas_pull= getPhiResidual(p1Fit.Phi(), p1meas.Phi())/ std::sqrt( p1meas_err[2]*p1meas_err[2] - cov[20]);
		
				P2_k_FitMeas_pull= (1/p2Fit.Perp() - 1/p2meas.Perp())/ std::sqrt( p2meas_err[0]*p2meas_err[0] - cov[30]) ;
				P2_Theta_FitMeas_pull= (p2Fit.Theta() - p2meas.Theta())/ std::sqrt( p2meas_err[1]*p2meas_err[1] - cov[40] );
				P2_Phi_FitMeas_pull= getPhiResidual(p2Fit.Phi(), p2meas.Phi())/ std::sqrt( p2meas_err[2]*p2meas_err[2] - cov[50]);
				
		
				Gamma_E_FitMeas_pull= (gammaFit.E() - gammameas.E())/ std::sqrt( gammameas_err[0]*gammameas_err[0] - cov[60] );
				Gamma_Theta_FitMeas_pull= (gammaFit.Theta() - gammameas.Theta())/ std::sqrt( gammameas_err[1]*gammameas_err[1]- cov[70]);
				Gamma_Phi_FitMeas_pull=getPhiResidual(gammaFit.Phi(), gammameas.Phi())/std::sqrt( gammameas_err[2]*gammameas_err[2] - cov[80]) ;			
		//		if(Gamma_Phi_FitMeas_pull > 100){
		//			 std::cout<<"crazy!!"<<std::endl;
		//			std::cout<<gammaFit.Phi()<<" "<<gammameas.Phi()<<" "<<getPhiResidual(gammaFit.Phi(), gammameas.Phi()) <<" "<<gammameas_err[2]<<" "<<cov[20]<<std::endl;
		//		}

				Parent_Px_FitMeas_pull=( parentFit.Px() - parentmeas.Px())/ std::sqrt(parentmeas_err[0]-parentfit_err[0]) ;
				Parent_Py_FitMeas_pull=( parentFit.Py() - parentmeas.Py())/ std::sqrt(parentmeas_err[1]-parentfit_err[1]) ;
				Parent_Pz_FitMeas_pull=( parentFit.Pz() - parentmeas.Pz())/ std::sqrt(parentmeas_err[2]-parentfit_err[2]) ;
				Parent_E_FitMeas_pull=( parentFit.E() - parentmeas.E())/ std::sqrt(parentmeas_err[3]-parentfit_err[3]) ;

				Chisq = fitter->getChi2();
				tree->Fill();
			}//end fit analysis option
	
			//add stuff to generator ttree
			if(_genAnalysis){

				
			//
				std::cout<<"test parent meas vector after tree fill "<<std::endl;
				for(int i= 0; i<parentmeas_err.size(); i++){
				std::cout<<parentmeas_err[i]<<" ";
				}
				std::cout<<std::endl;
			//get montecarlo 4 vectors
				int mcp1index = getCorrespondingMCParticleIndex(p1meas, 1, _daughterPID);
				int mcp2index = getCorrespondingMCParticleIndex(p2meas, -1, -_daughterPID);
				int mcpgindex = getCorrespondingMCParticleIndex(gammameas, 0,_neutralPID);
			//these pulls can be made and valid only when each particle has a MC match
				if(mcp1index != -1 && mcp2index != -1 && mcpgindex !=-1){
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

					P1_k_MeasGen_pull = (1/p1meas.Perp() -  1/mcp1.Perp())/ p1meas_err[0]; 
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
					Gamma_Phi_FitGen_pull = getPhiResidual(gammaFit.Phi(), mcgamma.Phi()) / std::sqrt(cov[80]);			


/*					if(fabs(Gamma_Phi_FitGen_pull)> 100 ||fabs( P2_Phi_FitGen_pull) > 100 || fabs(P1_Phi_FitGen_pull) >100){
						std::cout<<"crazy!! phi from fit gen"<<std::endl;
						std::cout<<"(P1,P2,Gamma) "<<P1_Phi_FitGen_pull<<" "<<P2_Phi_FitGen_pull<<" "<<Gamma_Phi_FitGen_pull<<std::endl;
						std::cout<<"residual (P1,P2,Gamma)"<<getPhiResidual(p1Fit.Phi(), mcp1.Phi())<<" "<<getPhiResidual(p2Fit.Phi(),mcp2.Phi())<<" "<<getPhiResidual(gammaFit.Phi(), mcgamma.Phi())<<std::endl;
						std::cout<<"errors (P1,P2,Gamma)"<<std::sqrt(cov[50])<<" "<<std::sqrt(cov[80])<<" "<<std::sqrt(cov[20])<<std::endl;
					
					}
					if(fabs(Gamma_Phi_MeasGen_pull)> 100 ||fabs(P2_Phi_MeasGen_pull)> 100 || fabs(P1_Phi_MeasGen_pull) >100){
						std::cout<<"crazzy!! phi from meas gen"<<std::endl;
					}*/
					genTree->Fill();	
                        	}
 				else{
					std::cout<<"Could not locate all corresponding MC particles"<<std::endl;
				}
			}//end genfit analysis option
                }//end prob cut
	}// end cov dim check
     
	//memory management
	pgamma.clear();
	pplus.clear();
	pminus.clear();
	pGammaPfoVec.clear();
	pPlusTrackVec.clear();
	pMinusTrackVec.clear();
	
	p1meas_err.clear();
	p2meas_err.clear();
	gammameas_err.clear();
	parentmeas_err.clear();

	p1fit_err.clear();
	p2fit_err.clear();
	gammafit_err.clear();
	parentfit_err.clear();

	gammagen_err.clear();
	delete gammaJet;
  	delete part1;
  	delete part2;
	delete fitter;
    }//end pvector size check
	//track events for each call
	evtNo++;
	
	return;
}

/*
      // Make the reconstructed particle using the result of the constrained fit
      // likely need some minimal fit probability cut - configurable for each gamma-gamma hypothesis
          ReconstructedParticleImpl * recoPart = new ReconstructedParticleImpl();
          recoPart->addParticle(pGammai);
          recoPart->addParticle(pGammaj);
          double Energy;
          double Mom[3];
      // The 4-vector of the fitted gamma-gamma is the sum of the two fitted 4-vectors.

      // Get the individual fitted photon 4-vectors - needed to calculate the four-momentum covariance matrix of the gamma-gamma system
          double Mom1[3],Mom2[3];
          double E1,E2;
          double pT1,pT2;

          E1 = j1.getE();
          Mom1[0] = j1.getPx();
          Mom1[1] = j1.getPy();
          Mom1[2] = j1.getPz();
          pT1 = sqrt(Mom1[0]*Mom1[0] + Mom1[1]*Mom1[1]);

          E2 = j2.getE();
          Mom2[0] = j2.getPx();
          Mom2[1] = j2.getPy();
          Mom2[2] = j2.getPz();
          pT2 = sqrt(Mom2[0]*Mom2[0] + Mom2[1]*Mom2[1]);

          Energy = E1 + E2;
          Mom[0] = Mom1[0] + Mom2[0];
          Mom[1] = Mom1[1] + Mom2[1];
          Mom[2] = Mom1[2] + Mom2[2];

          FloatVec covP(10, 0.0);
      //  (px,py,pz,E)

      //  Follow for now implementation in FourMomentumCovMat based on TMatrixD

          if(cov_dim == 6 ){
      // Note some cases of cov_dim != 6 have been seen ...
             const int nrows      = 6; // n rows jacobian
             const int ncolumns   = 4; // n columns jacobian //these look backwards

      //  Transformation matrix D (6*4) between (E1, theta1, phi1, E2, theta2, phi2) and (px, py, pz, E)
             double jacobian_by_columns[nrows*ncolumns] = {
                Mom1[0]/E1, Mom1[0]*Mom1[2]/pT1, -Mom1[1], Mom2[0]/E2, Mom2[0]*Mom2[2]/pT2, -Mom2[1],
                Mom1[1]/E1, Mom1[1]*Mom1[2]/pT1,  Mom1[0], Mom2[1]/E2, Mom2[1]*Mom2[2]/pT2,  Mom2[0],
                Mom1[2]/E1,        -pT1        ,    0.0  , Mom2[2]/E2,        -pT2        ,     0.0 ,
                   1.0    ,         0.0        ,    0.0  ,    1.0    ,         0.0        ,     0.0   };

      //  Now set up to calculate the new covariance matrix, namely  V' = D^T V D
             TMatrixD Dmatrix(nrows,ncolumns, jacobian_by_columns, "F");
             TMatrixD Vmatrix(nrows,nrows, cov, "F");                      // May need to have the filling array explicitly dimensioned ??
             TMatrixD Covmatrix(ncolumns,ncolumns);                        // Hopefully this is good enough for initialization ?

             Covmatrix.Mult( TMatrixD( Dmatrix, TMatrixD::kTransposeMult, Vmatrix) ,Dmatrix);   // Is this what we want ?

             for(int ii=0;ii<4;ii++){
                 for(int jj=0;jj<4;jj++){
                    if(_printing>3)std::cout << "4-vector cov. " << ii << " " << jj << " " << Covmatrix(ii,jj) << std::endl;
                 }                  
             }

             covP[0] = Covmatrix(0,0); // px-px
             covP[1] = Covmatrix(0,1); // px-py
             covP[2] = Covmatrix(1,1); // py-py
             covP[3] = Covmatrix(0,2); // px-pz
             covP[4] = Covmatrix(1,2); // py-pz
             covP[5] = Covmatrix(2,2); // pz-pz
             covP[6] = Covmatrix(0,3); // px-E
             covP[7] = Covmatrix(1,3); // py-E
             covP[8] = Covmatrix(2,3); // pz-E
             covP[9] = Covmatrix(3,3); // E-E
          }
          else{
               if(_printing>1)std::cout << " Fit has no valid covariance information (cov_dim = " << cov_dim << " )" << std::endl;              
          }

          recoPart->setEnergy( Energy );
          recoPart->setMomentum( Mom );
          recoPart->setMass( mass_constraint );
          recoPart->setCharge( 0.0 );
          recoPart->setCovMatrix( covP );

      // For now - store the fit probability in the goodnessofPiD variable
          float goodnessOfPID = float(fit_probability);
//          if( mgg  < _ggResonanceMass)goodnessOfPID = - goodnessOfPID;
          recoPart->setGoodnessOfPID( goodnessOfPID );

      // Default choice is Pi0
          recoPart->setType( 111 );
          if(_ggResonanceName=="Eta")recoPart->setType( 221 );
          if(_ggResonanceName=="EtaPrime")recoPart->setType( 331 );

          if(_printing>1)std::cout << "Fitting " << _ggResonanceName << " gg candidate  M= " 
                                   << mgg << " " << i << " " << j << " " << pgi.e() << " " << pgj.e() << 
                                   " Fitted energies: " << j1.getE() << " " << j2.getE() 
                                   << " Fit probability = " << fit_probability << std::endl; 

      // Add it to the collection if it exceeds the required miminum fit probability
          if(fit_probability > _fitProbabilityCut){
             if(_printing>1)std::cout << "GG candidate passes fit probability cut - KEEP" << std::endl;
             recparcol->addElement( recoPart );
          }
        }
      }
    }
  }

  return;
} */
