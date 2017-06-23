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

//#include "LeptonFitObject.h"
#include "TrackParticleFitObject.h"
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
  

  int getCorrespondingMCParticleIndex(TLorentzVector rec);
  void FindDiTrackGammaCandidates( LCCollectionVec* recparcol);

   OPALFitterGSL* setUPFit(vector<int> neutralIndices, vector<int> chargedIndices, vector<TLorentzVector> pneutral, vector<TLorentzVector> ptrack, vector<ReconstructedParticle*> pNeutralVec, vector<Track*> pTrackVec);
  std::vector<double> getChargedParticleErrors(TLorentzVector pcharged, Track* ptrk);
  std::vector<double> getNeutralErrors(TLorentzVector pneutral, ReconstructedParticle* pNeutral);

  void PrintCov(FloatVec cov, int dim);
  void PrintCov(double* cov, int dim);
  void setFitErrors(double* cov);
  double* ConstructParentMeasCovMatrix();
  std::vector<double> ConstructChargedSubMatrix(std::vector<double> p, TLorentzVector ptlv);
  std::vector<double> ConstructNeutralSubMatrix(TLorenztVector p);
  double* ConcatSubMatrices(std::vector<std::vector<double> > matrices);
  void setParentErrors(FloatVec meascov, FloatVec fitcov);
  void generateSubsets(std::vector<int> v, int k, int start, int currLen, std::vector<bool> used, std::vector<vector<int> >& combinations);
  
	std::vector<std::vector<int> > generateIndicesCombinations(int vectorsize, int nparticles);

	ReconstructedParticleImpl* constructFitParticle(TLorentzVector fitp, ReconstructedParticle* measrp);

  Track* constructTrack(TLorentzVector fitp, Track* meast);
  std::vector<double> buildTrackVector(Track* t);
  std::vector<double> buildFitTrackVector(TrackParticleFitObject* tfo);

  FloatVec ConstructCovMatrix(std::vector<TLorentzVector> charged, std::vector<TLorentzVector> neutral, double* cov);

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
  vector<vector<double> > measTrack;
  vector<TLorentzVector> fitNeutral;
  vector<TLorentzVector> fitCharged;
  vector<vector<double> > fitTrack;

  vector<vector<double> > measNeutral_err;
  vector<vector<double> > measCharged_err;
  vector<vector<double> > fitNeutral_err;
  vector<vector<double> > fitCharged_err;

  vector<double> measNeutralPdg;
  vector<double> measChargedPdg;

	//vector<double>: (E || k, Theta, Phi); pull order
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

	//vector<double>? (E , Theta, Phi); pull order even for charged particles
  vector<vector<double> > measgen_NeutralPulls;  
  vector<vector<double> > measgen_ChargedPulls;
  vector<vector<double> > fitgen_NeutralPulls;
  vector<vector<double> > fitgen_ChargedPulls;

  TLorentzVector genParent;

  vector<double> genNeutralPdg;
  vector<double> genChargedPdg;

  vector<vector<double> > mcTrack;//not implemented yet
  vector<TLorentzVector> mcCharged;
  vector<TLorentzVector> mcNeutral;


  
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
 
// JetFitObject* gammaJet;
  std::vector<JetFitObject*> neutralJets;

  //LeptonFitObject* part1;
 // LeptonFitObject* part2;
//  std::vector<LeptonFitObject*> chargedFO; 
    std::vector<TrackParticleFitObject*> TrackFO;
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
  std::string _outputTrackCollectionName;
  double _fitProbabilityCut;
 
  int _ifitter;
  int _fitAnalysis;
  int _genAnalysis;
 
  std::string m_rootFile;

protected:

} ;
