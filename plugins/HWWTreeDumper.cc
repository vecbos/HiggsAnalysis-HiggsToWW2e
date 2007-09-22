// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HWWTreeDumper
// Description:
//      Class HWWTreeDumper
//      Analyzer module that takes the Candidate Collections from
//      the analysis producers and dumps an ntuple
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenEvent.h"

#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsElectronFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTriggerTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleTrackerIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleCaloIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTreeDumper.h"


HWWTreeDumper::HWWTreeDumper(const edm::ParameterSet& iConfig)
{
  
  nameFile_     = iConfig.getUntrackedParameter<std::string>("nameFile", "RootOutput.root");
  nameTree_     = iConfig.getUntrackedParameter<std::string>("nameTree", "BaseTree");
  dumpTree_     = iConfig.getUntrackedParameter<bool>("dumpTree", false);
  dumpMCTruth_  = iConfig.getUntrackedParameter<bool>("dumpMCTruth", false);
  doMCmatch_    = iConfig.getUntrackedParameter<bool>("doMCmatch", false);
  
  // control level of Reco Adapters in the tree
  saveTrk_      = iConfig.getUntrackedParameter<bool>("saveTrk", false);
  saveEcal_     = iConfig.getUntrackedParameter<bool>("saveEcal", false);
  saveHcal_     = iConfig.getUntrackedParameter<bool>("saveHcal", false);
  saveDT_       = iConfig.getUntrackedParameter<bool>("saveDT", false);
  saveCSC_      = iConfig.getUntrackedParameter<bool>("saveCSC", false);
  saveRPC_      = iConfig.getUntrackedParameter<bool>("saveRPC", false);
  saveFatTrk_   = iConfig.getUntrackedParameter<bool>("saveFatTrk", false);
  saveFatEcal_  = iConfig.getUntrackedParameter<bool>("saveFatEcal", false);
  saveFatHcal_  = iConfig.getUntrackedParameter<bool>("saveFatHcal", false);
  saveFatDT_    = iConfig.getUntrackedParameter<bool>("saveFatDT", false);
  saveFatCSC_   = iConfig.getUntrackedParameter<bool>("saveFatCSC", false);
  saveFatRPC_   = iConfig.getUntrackedParameter<bool>("saveFatRPC", false);
  saveJetAlpha_ = iConfig.getUntrackedParameter<bool>("saveJetAlpha", false);

  // particle identification
  saveEleID_    = iConfig.getUntrackedParameter<bool>("saveEleID", false);

  // basic kinematic informations
  saveCand_  = iConfig.getUntrackedParameter<bool>("saveCand", true);
  
  // Candidate Collections
  dumpElectrons_  = iConfig.getUntrackedParameter<bool>("dumpElectrons", false);
  dumpMuons_      = iConfig.getUntrackedParameter<bool>("dumpMuons", false);
  dumpJets_       = iConfig.getUntrackedParameter<bool>("dumpJets", false);
  dumpGenJets_    = iConfig.getUntrackedParameter<bool>("dumpGenJets", false);
  dumpMet_        = iConfig.getUntrackedParameter<bool>("dumpMet", false);
  dumpGenMet_     = iConfig.getUntrackedParameter<bool>("dumpGenMet", false);
  
  // jet vertex collections
  jetVertexAlphaCollection_ = iConfig.getParameter<edm::InputTag>("jetVertexAlphaCollection");
  jetVertexBetaCollection_  = iConfig.getParameter<edm::InputTag>("jetVertexBetaCollection");

  electronCollection_ = iConfig.getParameter<edm::InputTag>("electronCollection");
  //   muonCollection_     = iConfig.getParameter<edm::InputTag>("muonCollection");
  jetCollection_      = iConfig.getParameter<edm::InputTag>("jetCollection");
  genJetCollection_   = iConfig.getParameter<edm::InputTag>("genJetCollection");
  metCollection_      = iConfig.getParameter<edm::InputTag>("metCollection");
  genMetCollection_   = iConfig.getParameter<edm::InputTag>("genMetCollection");
  mcTruthCollection_  = iConfig.getParameter<edm::InputTag>("mcTruthCollection");
  electronMatchMap_   = iConfig.getParameter<edm::InputTag>("electronMatchMap");
  //  muonMatchMap_       = iConfig.getParameter<edm::InputTag>("muonMatchMap");

  // trigger Collections
  triggerInputTag_  = iConfig.getParameter<edm::InputTag>("TriggerResultsTag") ;
}


HWWTreeDumper::~HWWTreeDumper()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HWWTreeDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // get MC truth
  CmsMcTruthTreeFiller treeFill(tree_);
  Handle<CandidateCollection> genParticleHandle;
  try { iEvent.getByLabel(mcTruthCollection_, genParticleHandle); }
  catch ( cms::Exception& ex ) { LogWarning("HWWTreeDumper") << "Can't get MC Truth Collection: " << mcTruthCollection_; }
  const CandidateCollection *genParticleCollection = genParticleHandle.product();

  if(dumpMCTruth_) {
    treeFill.writeCollectionToTree(genParticleCollection);
  }

  // fill the trigger paths info
  edm::Handle<edm::TriggerResults> trh;
  try {iEvent.getByLabel(triggerInputTag_,trh);} 
  catch( cms::Exception& ex ) { LogWarning("HWWTreeDumper") << "Trigger results: " << triggerInputTag_ << " not found"; }
  if (!trh.isValid())
    throw cms::Exception("ProductNotValid") << "TriggerResults product not valid";
  CmsTriggerTreeFiller triggerTreeFill (tree_) ;
  std::string prefix ("") ;
  std::string suffix ("Trg") ;
  triggerTreeFill.writeTriggerToTree (trh,prefix,suffix) ;

  // fill Electrons block
  if(dumpElectrons_) {
    CmsElectronFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Ele");
    treeFill.saveCand(saveCand_);
    treeFill.saveTrk(saveTrk_);
    treeFill.saveEcal(saveEcal_);
    treeFill.saveFatTrk(saveFatTrk_);
    treeFill.saveFatEcal(saveFatEcal_);
    treeFill.setMatchMap(electronMatchMap_);
    treeFill.saveEleID(true);

    Handle<CandidateCollection> electronCollectionHandle;
    try { iEvent.getByLabel(electronCollection_, electronCollectionHandle); }
    catch ( cms::Exception& ex ) { LogWarning("HWWTreeDumper") << "Can't get Candidate Collection: " << electronCollection_; }
    const CandidateCollection *electronCollection = electronCollectionHandle.product();
    treeFill.writeCollectionToTree(electronCollection, iEvent, iSetup, prefix, suffix, false);
    if(doMCmatch_) {
      treeFill.writeMcIndicesToTree(electronCollection, iEvent, iSetup, genParticleCollection, prefix, suffix, false);
    }
  }

  
  // fill MET block
  if(dumpMet_) {
    CmsCandidateFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Met");
    treeFill.saveCand(saveCand_);

    Handle<CandidateCollection> metCollectionHandle;
    try { iEvent.getByLabel(metCollection_, metCollectionHandle); }
    catch ( cms::Exception& ex ) { LogWarning("HWWTreeDumper") << "Can't get Candidate Collection: " << metCollection_; }
    const CandidateCollection *metCollection = metCollectionHandle.product();
    LogDebug("HWWTreeDumper") << "Met collection size = " << metCollection->size();
    treeFill.writeCollectionToTree(metCollection, iEvent, iSetup, prefix, suffix, false);

    // dump generated MET
    if(dumpGenMet_) {
      std::string suffix("GenMet");
      Handle<CandidateCollection> genMetCollectionHandle;
      try { iEvent.getByLabel(genMetCollection_, genMetCollectionHandle); }
      catch ( cms::Exception& ex ) { LogWarning("HWWTreeDumper") << "Can't get Candidate Collection: " <<genMetCollection_; }
      const CandidateCollection *genMetCollection = genMetCollectionHandle.product();
      treeFill.writeCollectionToTree(genMetCollection, iEvent, iSetup, prefix, suffix, false);
    }
  }

  // fill JET block (this is temporary... to be replaced with the collection from Vale producer )
  if(dumpJets_) {
    CmsJetFiller treeFill(tree_, jetVertexAlphaCollection_, jetVertexBetaCollection_, true);
    std::string prefix("");
    std::string suffix("Jet");
    treeFill.saveCand(saveCand_);
    treeFill.saveJetExtras(saveJetAlpha_);

    Handle<CandidateCollection> jetCollectionHandle;
    try { iEvent.getByLabel(jetCollection_, jetCollectionHandle); }
    catch ( cms::Exception& ex ) { LogWarning("HWWTreeDumper") << "Can't get Candidate Collection: " << jetCollection_; }
    const CandidateCollection *jetCollection = jetCollectionHandle.product();
    LogDebug("HWWTreeDumper") << "Jet collection size = " << jetCollection->size();
    treeFill.writeCollectionToTree(jetCollection, iEvent, iSetup, prefix, suffix, false);

    // dump generated JETs
    if(dumpGenJets_) {
      std::string suffix("GenJet");
      Handle<CandidateCollection> genJetCollectionHandle;
      try { iEvent.getByLabel(genJetCollection_, genJetCollectionHandle); }
      catch ( cms::Exception& ex ) { LogWarning("HWWTreeDumper") << "Can't get Candidate Collection: " <<genJetCollection_; }
      const CandidateCollection *genJetCollection = genJetCollectionHandle.product();
      treeFill.writeCollectionToTree(genJetCollection, iEvent, iSetup, prefix, suffix, false);
    }
  }

  if(dumpTree_) tree_->dumpData();

}


// ------------ method called once each job just before starting event loop  ------------
void 
HWWTreeDumper::beginJob(const edm::EventSetup&)
{
  fileOut_ = TFile::Open(nameFile_.c_str(), "RECREATE");
  tree_  = new  CmsTree(nameTree_.c_str(),nameTree_.c_str());

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HWWTreeDumper::endJob() {
  fileOut_->cd();
  TTree* treeOut = tree_->getTree();
  treeOut->Write();
  fileOut_->Close();
}


