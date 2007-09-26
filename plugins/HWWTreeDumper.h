// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    HWWTreeDumper
// Class:      HWWTreeDumper
// 
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
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
  // ----------member data ---------------------------
  std::string nameTree_, nameFile_;
  bool dumpTree_, dumpMCTruth_, dumpMCmatch_, doMCmatch_;
  bool saveTrk_, saveEcal_, saveHcal_, saveDT_, saveCSC_, saveRPC_;
  bool saveFatTrk_, saveFatEcal_, saveFatHcal_, saveFatDT_, saveFatCSC_, saveFatRPC_;
  bool saveEleID_;
  bool saveJetAlpha_;
  bool saveCand_;
  bool dumpElectrons_, dumpMuons_;
  bool dumpJets_, dumpGenJets_, dumpMet_, dumpGenMet_;
  bool dumpTriggerResults_;

  edm::InputTag electronCollection_, muonCollection_;
  edm::InputTag jetCollection_, genJetCollection_, metCollection_, genMetCollection_;
  edm::InputTag mcTruthCollection_, electronMatchMap_, muonMatchMap_;

  edm::InputTag jetVertexAlphaCollection_;

  edm::InputTag triggerInputTag_ ;

  TFile *fileOut_;
  CmsTree *tree_;
  

};
#endif // HWWTreeDumper_h
