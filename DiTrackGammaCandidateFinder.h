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

 //four vector quantities for measured and fit
  TLorentzVector gammameas;
  TLorentzVector p1meas;
  TLorentzVector p2meas;
  TLorentzVector parentmeas;

  TLorentzVector gammaFit;
  TLorentzVector p1Fit;
  TLorentzVector p2Fit;
  TLorentzVector parentFit;

  //errors on parameters
  //[ dE, dTheta, dPhi ] and [ dK, dTheta, dPhi ]
  std::vector<double> gammameas_err;
  std::vector<double> p1meas_err;
  std::vector<double> p2meas_err;

  std::vector<double> gammafit_err;
  std::vector<double> p1fit_err;
  std::vector<double> p2fit_err;

  //error on parent
  //[dPx, dPy, dPz, dE]
  std::vector<double> parentmeas_err;
  std::vector<double> parentfit_err;

 //photon errors based on generator vales
 std::vector<double> gammagen_err;

  double FitProbability;
  double Chisq; 
  int evtNo;

  //pulls for fit-meas/sqrt( var(meas)-var(fit))
  double P1_k_FitMeas_pull;
  double P1_Theta_FitMeas_pull;
  double P1_Phi_FitMeas_pull;
  double P2_k_FitMeas_pull;
  double P2_Theta_FitMeas_pull;
  double P2_Phi_FitMeas_pull;
  double Gamma_E_FitMeas_pull;
  double Gamma_Theta_FitMeas_pull;
  double Gamma_Phi_FitMeas_pull;

 //pulls for fit-emas/sqrt( var(meas)-var(fit)) parent particle
  double Parent_Px_FitMeas_pull;
  double Parent_Py_FitMeas_pull;
  double Parent_Pz_FitMeas_pull;
  double Parent_E_FitMeas_pull;

 ////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////
  // Vars for gen tree
  /////////////////////////////////////////////////////////
  //pulls for meas-gen/sigma(meas)
  double P1_k_MeasGen_pull;
  double P1_Theta_MeasGen_pull;
  double P1_Phi_MeasGen_pull;
  double P2_k_MeasGen_pull;
  double P2_Theta_MeasGen_pull;
  double P2_Phi_MeasGen_pull;
  double Gamma_E_MeasGen_pull;
  double Gamma_Theta_MeasGen_pull;
  double Gamma_Phi_MeasGen_pull;

  //pulls for fit-gen/sigma(fit)
  double P1_k_FitGen_pull;
  double P1_Theta_FitGen_pull;
  double P1_Phi_FitGen_pull;
  double P2_k_FitGen_pull;
  double P2_Theta_FitGen_pull;
  double P2_Phi_FitGen_pull;
  double Gamma_E_FitGen_pull;
  double Gamma_Theta_FitGen_pull;
  double Gamma_Phi_FitGen_pull;

  //generator money plots
  double mcESum;
  double mcM;
  //generator four vector quantities
  TLorentzVector mcp1;
  TLorentzVector mcp2;
  TLorentzVector mcgamma;   
  ////////////////////////////////////////////////////////

  
private:

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
  std::string _resonanceName;
  float _gammaMomentumCut;
  float _resonanceMass;
  float _dmcut;
  float _daughterDecayMass;
  double _fitProbabilityCut;
  float _minPtCut;
  int _ifitter;
  int _fitAnalysis;
  int _genAnalysis;
  int _resonancePID;
  int _daughterPID;
  int _neutralPID;
  float _mcDeltaCharged;
  float _mcDeltaGamma;
  float _mcVertex;
  int _usePhotonCalibration;
  int _useTrackCalibration;
  std::string m_rootFile;

protected:

} ;
