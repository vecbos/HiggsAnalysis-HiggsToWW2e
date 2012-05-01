// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    HiggsAnalysis/HiggsToWW2e
// Class:      HWWTreeDumper
// 
//-----------------------------------------------------------------------

#ifndef HWWTreeDumper_h
#define HWWTreeDumper_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsHLTObjectFiller.h"

#include <TTree.h>
#include <TFile.h>


#include "HiggsAnalysis/HiggsToGammaGamma/interface/EGEnergyCorrector.h"
#include "HiggsAnalysis/HiggsToGammaGamma/interface/PhotonFix.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "CondFormats/EcalObjects/interface/EcalFunctionParameters.h" 

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

class HWWTreeDumper : public edm::EDAnalyzer {
public:
  explicit HWWTreeDumper(const edm::ParameterSet&);
  ~HWWTreeDumper();
  
  
private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup );
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
private:
  
  //! name of the output ROOT file
  std::string nameFile_;
  //! name of the tree with the events
  std::string nameTree_;
  //! effectively dump the event tree
  bool dumpTree_;
  //! dump the MC truth (generator) block
  bool dumpMCTruth_;
  //! do the match between reco and truth particles, dump indices into tree
  bool doMCEleMatch_, doMCMuonMatch_;
  //! dump generator informations for the CSA07 soups
  bool dumpGenInfo_;
  //! dump H->WW->2lep 2nu preselection marker
  bool dumpPreselInfo_;
  //! dump H->WW->2lep 2nu gg fusion signal k-factor
  bool dumpSignalKfactor_;
  //! save the LHE event infos, PDF infos
  bool dumpLHE_, dumpPdfWeight_;
  //! dump the basic candidate informations (4-vectors)
  bool saveCand_;
  //! dump specific reco informations in addition to candidate variables
  bool saveTrk_, saveEcal_, saveHcal_, saveDT_, saveCSC_, saveRPC_;
  //! dump more specific reco informations in addition to candidate variables
  bool saveFatTrk_, saveFatEcal_, saveFatHcal_, saveFatDT_, saveFatCSC_, saveFatRPC_;
  //! dump more specific reco informations in addition to candidate variables
  bool savePFTauBasic_, savePFTauDiscriminators_, saveLeadPFCand_;
  //! dump electron ID variables
  bool saveEleID_;
  //! dump the particle flow electron block
  bool savePFEleBasic_; 
  bool savePFEleIsoDep_;
  //! dump the jet alpha parameter
  bool saveJetAlpha_;
  //! dump the jet b-tag output
  bool saveJetBTag_;
  //! dump the electron block
  bool dumpElectrons_;
  //! dump the photon block
  bool dumpPhotons_;
  //! dump the conversion block
  bool dumpConversions_;
  //! dump the particle flow electron block
  bool dumpPFlowElectrons_;
  //! dump the particle flow electron pre-identification block
  bool dumpPFpreId_; 
  //! dump muon block
  bool dumpMuons_;
  //! dump pftau block
  bool dumpPFTaus_, dumphpsPFTaus_, dumphpsTancTaus_;
  //! dump PFCandidates
  bool dumpPFCandidates_;
  //! dump reco / generated / PU corrected jets block
  bool dumpJets_, dumpGenJets_, dumpPUcorrPFJet_;
  //! dump reco / generated MET block
  bool dumpMet_, dumpGenMet_;
  //! dump Super/Basic Clusters block
  bool dumpSCs_, dumpBCs_;
  //! dump tracks
  bool dumpTracks_, dumpGsfTracks_, dumpMuonTracks_;
  //! dump the primary vertices of the event
  bool dumpVertices_;
  //! bool dump the block of the V0 candidates (as K0s)
  bool dumpK0s_;
  //! dump trigger results
  bool dumpTriggerResults_;
  //! dump the Particle Flow objects
  bool dumpParticleFlowObjects_;
  //! dump hcal noise flags
  bool dumpHcalNoiseFlags_;
  bool aodHcalNoiseFlags_;
  //! dump the run info informations
  bool dumpRunInfo_;
  //! PDF collections
  edm::InputTag pdfSet1_, pdfSet2_, pdfSet3_; 
  std::string namePdf1_, namePdf2_, namePdf3_;
  //! save the dE/dx of the tracks (requires the right module to be run)
  bool saveTrackDeDx_;
  //! save the calotowers
  bool dumpCaloTowers_;

  bool usePhotonFix_;

  bool useEnergyRegression_;

  std::string energyRegressionElectronFile_;
  std::string energyRegressionPhotonFile_;

  //! candidate collections in input
  edm::InputTag electronCollection_, muonCollection_,pflowElectronCollection_;
  edm::InputTag photonCollection_;
  edm::InputTag jetCollection1_, genJetCollection_, jetCollection2_, jetCollection3_;
  edm::InputTag PFpuCorrJetCollection1_, PFpuCorrJetCollection2_, PFpuCorrJetCollection3_;
  edm::InputTag metCollection_, TCmetCollection_, genMetCollection_;
  // edm::InputTag corrmetCollection_;
  edm::InputTag vertexCollection_;
  edm::InputTag K0sCollection_;
  edm::InputTag PFjetCollection1_, PFjetCollection2_, PFjetCollection3_, PFmetCollection_, PFChMetCollection_;
  edm::InputTag leptonLinkedPFCandidates_;
  edm::InputTag JPTjetCollection1_, JPTjetCollection2_;
  edm::InputTag chargedMetCollection_;
  //! btag collections (only PF, the other are hardcoded)
  edm::ParameterSet PFJetsBTags_, PFPUcorrJetsBTags_;
  //! jet id mva 
  std::vector<edm::InputTag> PFjetMvaIdCollection_, PFpujetMvaIdCollection_;
  //! track collections
  edm::InputTag trackCollection_, refittedForDeDxTrackCollection_, gsfTrackCollection_;
  edm::InputTag globalMuonTrackCollection_, standAloneMuonTrackCollection_;
  //! taus
  edm::InputTag pfTauCollection_, hpspfTauCollection_, hpsTancTausCollection_;
  //! PF candidates
  edm::InputTag PFCandidateCollection_;
  //! supercluster collections in input
  edm::InputTag ecalSCCollection_; // merged ECAL Superclusters
  edm::InputTag ecalBarrelSCCollection_, ecalEndcapSCCollection_, ecalElePFClusterCollection_, ecalPhoPFClusterCollection_;
  //! basiccluster collections in input
  edm::InputTag ecalBCCollection_;
  //! ECAL rechits to compute the cluster shapes on the fly
  edm::InputTag ecalBarrelRecHits_, ecalEndcapRecHits_, esRecHits_;
  //! coversions collection
  edm::InputTag conversions_;
  //! track and calotowers collection for isolation studies
  edm::InputTag calotowersForIsolationProducer_;
  //! calotowers collections
  edm::InputTag calotowerCollection_, hbheLabel_, hoLabel_, hfLabel_;
  std::vector<edm::InputTag> ecalLabels_;
  //! generator-level particle collection in input
  edm::InputTag mcTruthCollection_;
  //! association map between reco and generated particles
  edm::InputTag electronMatchMap_, muonMatchMap_;
  //! generator-level informations present in the soups
  edm::InputTag hepMcCollection_, genInfoCollection_;
  std::string genWeightCollection_;
  //! PF electrons pre-identification
  edm::InputTag PFpreIdCollection_;
  //! PF Tau Discriminators
  edm::InputTag tauDiscrByLeadingTrackFindingTag_, tauDiscrByLeadingTrackPtCutTag_, tauDiscrByLeadingPionPtCutTag_, 
    tauDiscrByIsolationTag_, tauDiscrByIsolationUsingLeadingPionTag_,
    tauDiscrByTrackIsolationUsingLeadingPionTag_, tauDiscrByTrackIsolationTag_,
    tauDiscrByECALIsolationTag_, tauDiscrByECALIsolationUsingLeadingPionTag_, 
    tauDiscrAgainstMuonTag_, tauDiscrAgainstElectronTag_,
    tauDiscrByTaNCTag_,
    tauDiscrByTaNCfrHalfPercentTag_, tauDiscrByTaNCfrOnePercentTag_,
    tauDiscrByTaNCfrQuarterPercentTag_, tauDiscrByTaNCfrTenthPercentTag_;
  //! HPS PF Tau Discriminators
  edm::InputTag hpsTauDiscrByLooseElectronRejectionTag_, hpsTauDiscrByMediumElectronRejectionTag_, hpsTauDiscrByTightElectronRejectionTag_,
    hpsTauDiscrByLooseMuonRejectionTag_, hpsTauDiscrByTightMuonRejectionTag_,
    hpsTauDiscrByDecayModeFindingTag_,
    hpsTauDiscrByVLooseIsolationTag_, hpsTauDiscrByLooseIsolationTag_,
    hpsTauDiscrByMediumIsolationTag_, hpsTauDiscrByTightIsolationTag_;
  //! HPS Tanc Tau Discriminators
  edm::InputTag hpsTancTausDiscrByLeadingTrackFindingTag_, hpsTancTausDiscrByLeadingTrackPtCutTag_, hpsTancTausDiscrByLeadingPionPtCutTag_,
    hpsTancTausDiscrByTancTag_, hpsTancTausDiscrByTancRawTag_,
    hpsTancTausDiscrByTancVLooseTag_, hpsTancTausDiscrByTancLooseTag_, hpsTancTausDiscrByTancMediumTag_, hpsTancTausDiscrByTancTightTag_,
    hpsTancTausDiscrByLooseElectronRejectionTag_, hpsTancTausDiscrByMediumElectronRejectionTag_, hpsTancTausDiscrByTightElectronRejectionTag_,
    hpsTancTausDiscrByLooseMuonRejectionTag_, hpsTancTausDiscrByTightMuonRejectionTag_,
    hpsTancTausDiscrByDecayModeSelectionTag_,
    hpsTancTausDiscrByVLooseIsolationTag_, hpsTancTausDiscrByLooseIsolationTag_,
    hpsTancTausDiscrByMediumIsolationTag_, hpsTancTausDiscrByTightIsolationTag_;

  //! Hcal noise summary object
  edm::InputTag hcalNoiseSummaryLabel_;
  //! ROOT file with the plain ROOT tree inside
  TFile *fileOut_;
  //! the tree with the events
  CmsTree *tree_;

  //! number of the processed event
  int jevt_;
  int jevtInRun_;

  //! need to keep the HLTObjectDumper to update the trigger configuration on run boundaries
  bool dumpHLTObject_;
  CmsHLTObjectFiller* hltObjectFiller_;
  edm::ParameterSet hltParms_; //parameters for HLTObject filler

  std::vector<std::string>* trgNames_;
  std::vector<std::string>* LHEComments_;

  edm::ParameterSet phFixElePar_;
  edm::ParameterSet phFixPhoPar_;

  edm::ParameterSet posCalcParameters_;

  EGEnergyCorrector* eCorrector_;
  EGEnergyCorrector* pCorrector_;
  PhotonFix* phFixE_;
  PhotonFix* phFixP_;

  EcalClusterFunctionBaseClass* energyCorrectionF;

  PositionCalc* posCalculator_;
};
#endif // HWWTreeDumper_h
