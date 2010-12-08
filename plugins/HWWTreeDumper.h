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
  //! dump the particle flow electron block
  bool dumpPFlowElectrons_;
  //! dump the particle flow electron pre-identification block
  bool dumpPFpreId_; 
  //! dump muon block
  bool dumpMuons_;
  //! dump pftau block
  bool dumpPFTaus_;
  //! dump reco / generated jets block
  bool dumpJets_, dumpGenJets_;
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
  //! dump the run info informations
  bool dumpRunInfo_;
  //! save the dE/dx of the tracks (requires the right module to be run)
  bool saveTrackDeDx_;
  //! save the calotowers
  bool dumpCaloTowers_;

  //! candidate collections in input
  edm::InputTag electronCollection_, muonCollection_,pflowElectronCollection_;
  edm::InputTag photonCollection_;
  edm::InputTag jetCollection1_, genJetCollection_, jetCollection2_;
  edm::InputTag metCollection_, TCmetCollection_, genMetCollection_;
  // edm::InputTag corrmetCollection_;
  edm::InputTag vertexCollection_;
  edm::InputTag K0sCollection_;
  edm::InputTag PFjetCollection1_, PFjetCollection2_, PFmetCollection_;
  edm::InputTag JPTjetCollection1_, JPTjetCollection2_;
  edm::InputTag trackCollection_, refittedForDeDxTrackCollection_, gsfTrackCollection_;
  edm::InputTag globalMuonTrackCollection_, standAloneMuonTrackCollection_;
  edm::InputTag pfTauCollection_;
  //! supercluster collections in input
  edm::InputTag ecalSCCollection_; // merged ECAL Superclusters
  edm::InputTag ecalBarrelSCCollection_, ecalEndcapSCCollection_, ecalPFClusterCollection_;
  //! basiccluster collections in input
  edm::InputTag ecalBCCollection_;
  //! ECAL rechits to compute the cluster shapes on the fly
  edm::InputTag ecalBarrelRecHits_, ecalEndcapRecHits_;
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
  edm::InputTag tauDiscrByLeadTrackFindingTag_, tauDiscrByLeadTrackPtCutTag_,// tauDiscrByNProngsTag_,
    tauDiscrByTrackIsoTag_, tauDiscrByEcalIsoTag_, tauDiscrAgainstMuonsTag_, tauDiscrAgainstElectronsTag_,
    tauDiscrByTaNCTag_,
    tauDiscrByTaNCfrHalfPercentTag_, tauDiscrByTaNCfrOnePercentTag_,
    tauDiscrByTaNCfrQuarterPercentTag_, tauDiscrByTaNCfrTenthPercentTag_;
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

};
#endif // HWWTreeDumper_h
