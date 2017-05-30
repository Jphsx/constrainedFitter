#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "TFile.h"
#include "TTree.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TLorentzVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;

#include "LeptonFitObject.h"
#include "JetFitObject.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"


using namespace lcio ;


/** DiTrackGammaCandidateFinder:<br>
 *
 * (modelled after GammaGammaCandidateFinder processor)
 * 
 * @author Justin Anguiano, University of Kansas
 * Code follows convention X -> p1+ p2- gamma
 */

class DiTrackGammaCandidateFinder : public marlin::Processor {
  
 public:
  
 virtual marlin::Processor*  newProcessor() { return new DiTrackGammaCandidateFinder ; }
  
  DiTrackGammaCandidateFinder() ;

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the proscessor, e.g. book histograms.
   */
  virtual void init() ;
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 


  /** Called after data processing for clean up.
   */
  virtual void end() ;

  bool FindPFOs( LCEvent* evt );
  bool FindTracks( LCEvent* evt );
  bool FindMCParticles( LCEvent* evt);
  
  int getCorrespondingMCParticleIndex(TLorentzVector vReco, int recoCharge, int recoPdg);
  void FindDiTrackGammaCandidates( LCCollectionVec* recparcol);
// BaseFitter* setUpFit(TLorentzVector gamma, TLorentzVector p1, TLorentzVector p2, Track* p1Track, Track* p2Track );
  OPALFitterGSL* setUpFit(TLorentzVector gamma, TLorentzVector p1, TLorentzVector p2, Track* p1Track, Track* p2Track);
  std::vector<double> getChargedParticleErrors(TLorentzVector pcharged, Track* ptrk);
  void calibratePhoton(TLorentzVector& v);
  void calibratePionError(std::vector<double>& errors);
  std::vector<double> getPhotonErrors(TLorentzVector pgamma);
  double getPhiResidual(double phi1, double phi2);
  void PrintCov(FloatVec cov, int dim);
  void PrintCov(double* cov, int dim);
  void setFitErrors(double* cov);
  double* ConstructParentMeasCovMatrix();
  void setParentErrors(FloatVec meascov, FloatVec fitcov);
//transforms 9x9 to 4x4
  FloatVec ConstructCovMatrix(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, double* cov);

//TTree stuff for fit analysis
  TFile *rootFile;
  //normal analysis tree for comparing fits to measured
  TTree *tree;

  //tree that includes generator level information
  TTree *genTree; 
  
 ///////////////////////////////////////////////////////
 // Vars for regular tree
 /////////////////////////////////////////////////////
  //TTree vars for fit analysis
  //money plots
  double RecoEnergy;
  double RecoMass;
  double FitEnergy;
  double FitProbability;
  double Chisq; 
  int evtNo;

//add variables
  vector<TLorentzVector> measNeutral;
  vector<TLorentzVector> measCharged;
  vector<TLorentzVector> fitNeutral;
  vector<TLorentzVector> fitCharged;

  vector<vector<double> > measNeutral_err;
  vector<vector<double> > measCharged_err;
  vector<vector<double> > fitNeutral_err;
  vector<vector<double> > fitCharged_err;

  vector<double> measNeutralPdg;
  vector<double> measChargedPdg;

	//vector<double>: (E|k, Theta, Phi); pull order
  vector<vector<double> > fitmeas_NeutralPulls;
  vector<vector<double> > fitmeas_ChargedPulls;

  TLorentzVector measParent;
  TLorentzVector fitParent;

  vector<double> measParent_err;
  vector<double> fitParent_err;

		//(Px, Py, Pz, E) pull order
  vector<double> fitmeas_ParentPulls;

  //monte carlo variables

  vector<TLorentzVector> genNeutral;
  vector<TLorentzVector> genCharged;

	//vector<double>? (E|k, Theta, Phi); pull order
  vector<vector<double> > measgen_NeutralPulls;
  vector<vector<double> > measgen_ChargedPulls;
  vector<vector<double> > fitgen_NeutralPulls;
  vector<vector<double> > fitgen_ChargedPulls;

  TLorentzVector genParent;

  vector<double> genNeutralPdg;
  vector<double> genChargedPdg;


  
private:



  //vectors for particle input information
	int _parentPdg;
	double _parentMass;
	double _parentCharge;
	int _nDaughters;
	int _nCharged;
	int _nNeutral;
	vector<double> _daughterChargedPdgs;
	vector<double> _daughterNeutralPdgs;
	vector<double> _daugherChargedMass;
	vector<double> _daughterNeutralMass;


  //global fitobject pointers
  //return Fit objects in the fitter is not retaining the derived class and chopping to base class
  JetFitObject* gammaJet;
  LeptonFitObject* part1;
  LeptonFitObject* part2; 
//  OPALFitterGSL* fitter;
    BaseFitter* ftest;

  std::vector<ReconstructedParticle*>_pfovec;
  std::vector<Track*> _trackvec;
  std::vector<MCParticle*> _mcpartvec;
  int   _printing;
  std::string _inputTrackCollectionName;
  std::string _inputParticleCollectionName;
  std::string _mcParticleCollectionName;
  std::string _outputParticleCollectionName;
  double _fitProbabilityCut;
 
  int _ifitter;
  int _fitAnalysis;
  int _genAnalysis;
 
  std::string m_rootFile;

protected:

} ;
