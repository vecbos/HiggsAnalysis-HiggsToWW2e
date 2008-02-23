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
#include <TFile.h>

class HWWTreeDumper : public edm::EDAnalyzer {
 public:
  explicit HWWTreeDumper(const edm::ParameterSet&);
  ~HWWTreeDumper();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
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
  //! dump the basic candidate informations (4-vectors)
  bool saveCand_;
  //! dump specific reco informations in addition to candidate variables
  bool saveTrk_, saveEcal_, saveHcal_, saveDT_, saveCSC_, saveRPC_;
  //! dump more specific reco informations in addition to candidate variables
  bool saveFatTrk_, saveFatEcal_, saveFatHcal_, saveFatDT_, saveFatCSC_, saveFatRPC_;
  //! dump electron ID variables
  bool saveEleID_;
  //! dump the jet alpha parameter
  bool saveJetAlpha_;
  //! dump the electron block
  bool dumpElectrons_;
  //! dump muon block
  bool dumpMuons_;
  //! dump reco / generated jets block
  bool dumpJets_, dumpGenJets_;
  //! dump reco / generated MET block
  bool dumpMet_, dumpGenMet_;
  //! dump SuperClusters block
  bool dumpSCs_;
  //! dump trigger results
  bool dumpTriggerResults_;
  //! candidate collections in input
  edm::InputTag electronCollection_, muonCollection_;
  edm::InputTag jetCollection_, genJetCollection_;
  edm::InputTag metCollection_, genMetCollection_;
  //! supercluster collections in input
  edm::InputTag hybridSCCollection_, islandSCCollection_;
  //! association maps between the clusters and the cluster shapes
  edm::InputTag ecalBarrelClusterShapes_, ecalEndcapClusterShapes_;
  //! association map between the electron and the electron ID outputs
  edm::InputTag electronIDAssocProducer_;
  //! collection of jet vertices for alpha evaluation
  edm::InputTag jetVertexAlphaCollection_;
  //! generator-level particle collection in input
  edm::InputTag mcTruthCollection_;
  //! association map between reco and generated particles
  edm::InputTag electronMatchMap_, muonMatchMap_;
  //! generator-level informations present in the soups
  edm::InputTag hepMcCollection_, genInfoCollection_, genWeightCollection_;
  //! results of the HLT
  edm::InputTag triggerInputTag_ ;

  //! ROOT file with the plain ROOT tree inside
  TFile *fileOut_;
  //! the tree with the events
  CmsTree *tree_;
  //! the tree with the run conditions (for now only the trigger words)
  CmsTree *treeConditions_;

  //! number of the processed event
  int jevt_;

};
#endif // HWWTreeDumper_h
