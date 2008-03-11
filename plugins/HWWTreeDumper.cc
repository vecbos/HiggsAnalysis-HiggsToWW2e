// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class HWWTreeDumper
//      Analyzer module that takes the Candidate Collections from
//      the analysis producers and dumps an ntuple
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
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenEvent.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsElectronFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsSuperClusterFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGenInfoFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConditionsFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTriggerTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleTrackerIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleCaloIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTreeDumper.h"


HWWTreeDumper::HWWTreeDumper(const edm::ParameterSet& iConfig)
{
  
  nameFile_      = iConfig.getUntrackedParameter<std::string>("nameFile", "RootOutput.root");
  nameTree_      = iConfig.getUntrackedParameter<std::string>("nameTree", "BaseTree");
  dumpTree_      = iConfig.getUntrackedParameter<bool>("dumpTree", false);
  dumpMCTruth_   = iConfig.getUntrackedParameter<bool>("dumpMCTruth", false);
  doMCEleMatch_  = iConfig.getUntrackedParameter<bool>("doMCEleMatch", false);
  doMCMuonMatch_ = iConfig.getUntrackedParameter<bool>("doMCMuonMatch", false);
  
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
  dumpPreselInfo_     = iConfig.getUntrackedParameter<bool>("dumpPreselInfo", false);
  dumpSignalKfactor_  = iConfig.getUntrackedParameter<bool>("dumpSignalKfactor", false);
  dumpGenInfoMcAtNlo_ = iConfig.getUntrackedParameter<bool>("dumpGenInfoMcAtNlo", false);
  dumpGenInfo_        = iConfig.getUntrackedParameter<bool>("dumpGenInfo", false);
  dumpElectrons_      = iConfig.getUntrackedParameter<bool>("dumpElectrons", false);
  dumpSCs_            = iConfig.getUntrackedParameter<bool>("dumpSCs", false);
  dumpMuons_          = iConfig.getUntrackedParameter<bool>("dumpMuons", false);
  dumpJets_           = iConfig.getUntrackedParameter<bool>("dumpJets", false);
  dumpGenJets_        = iConfig.getUntrackedParameter<bool>("dumpGenJets", false);
  dumpMet_            = iConfig.getUntrackedParameter<bool>("dumpMet", false);
  dumpGenMet_         = iConfig.getUntrackedParameter<bool>("dumpGenMet", false);
  
  // jet vertex collections
  jetVertexAlphaCollection_ = iConfig.getParameter<edm::InputTag>("jetVertexAlphaCollection");

  electronCollection_      = iConfig.getParameter<edm::InputTag>("electronCollection");
  muonCollection_          = iConfig.getParameter<edm::InputTag>("muonCollection");
  hybridSCCollection_      = iConfig.getParameter<edm::InputTag>("hybridSCCollection");
  islandSCCollection_      = iConfig.getParameter<edm::InputTag>("islandSCCollection");
  ecalBarrelClusterShapes_ = iConfig.getParameter<edm::InputTag>("ecalBarrelClusterShapes");
  ecalEndcapClusterShapes_ = iConfig.getParameter<edm::InputTag>("ecalEndcapClusterShapes");
  electronIDAssocProducer_ = iConfig.getParameter<edm::InputTag>("electronIDAssocProducer");
  tkIsolationProducer_     = iConfig.getParameter<edm::InputTag>("tkIsolationProducer"); 
  towerIsolationProducer_  = iConfig.getParameter<edm::InputTag>("towerIsolationProducer"); 
  jetCollection_           = iConfig.getParameter<edm::InputTag>("jetCollection");
  genJetCollection_        = iConfig.getParameter<edm::InputTag>("genJetCollection");
  metCollection_           = iConfig.getParameter<edm::InputTag>("metCollection");
  genMetCollection_        = iConfig.getParameter<edm::InputTag>("genMetCollection");
  mcTruthCollection_       = iConfig.getParameter<edm::InputTag>("mcTruthCollection");
  electronMatchMap_        = iConfig.getParameter<edm::InputTag>("electronMatchMap");
  muonMatchMap_            = iConfig.getParameter<edm::InputTag>("muonMatchMap");
  hepMcCollection_         = iConfig.getParameter<edm::InputTag>("hepMcCollection");
  genInfoCollection_       = iConfig.getParameter<edm::InputTag>("genInfoCollection");
  genWeightCollection_     = iConfig.getParameter<edm::InputTag>("genWeightCollection");

  // trigger Collections
  triggerInputTag_     = iConfig.getParameter<edm::InputTag>("TriggerResultsTag") ;
  dumpTriggerResults_  = iConfig.getUntrackedParameter<bool>("dumpTriggerResults");
}



HWWTreeDumper::~HWWTreeDumper() { }



// ------------ method called to for each event  ------------
void HWWTreeDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // get MC truth
  CmsMcTruthTreeFiller treeFill(tree_);

  if(dumpMCTruth_) {

    treeFill.writeCollectionToTree( mcTruthCollection_, iEvent );

  }



  // fill the trigger paths info
  if(dumpTriggerResults_) {

    CmsTriggerTreeFiller triggerTreeFill (tree_) ;
    std::string prefix ("") ;
    std::string suffix ("Trg") ;
    triggerTreeFill.writeTriggerToTree (triggerInputTag_, iEvent, prefix, suffix) ;

    /// fill the trigger mask in the Conditions tree
    if ( jevt_ == 1 ) {

      CmsConditionsFiller conditionsFiller( treeConditions_ );
      conditionsFiller.setHLTResults( triggerInputTag_, iEvent );
      conditionsFiller.writeConditionsToTree( iSetup, "", "", true);

      jevt_++;

    }

  }


  // fill preselection output
  if (dumpPreselInfo_) {

    Handle<bool> selected;
    iEvent.getByLabel("preselectionMarker", selected );
    bool isSelected = *selected;
    tree_->column ("evtPresel", isSelected, false, "Pres");

  }


  // fill signal gg fusion k-factor
  if (dumpSignalKfactor_) {
    Handle<double> theFactor;
    iEvent.getByLabel("KFactorProducer", theFactor );
    double theKfactor = *theFactor;
    tree_->column ("evtKfactor", theKfactor, 0., "kFac");
  } 

  // dump infos on MC production in case of MC@NLO production
  if (dumpGenInfoMcAtNlo_) {

    Handle<HepMCProduct> theEvt;
    iEvent.getByType(theEvt);
    HepMC::GenEvent* thisEvt = new  HepMC::GenEvent(*(theEvt->GetEvent()));
    HepMC::WeightContainer allWeights = thisEvt->weights();   
    double theWeight = allWeights.front();
    tree_->column ("evtMcAtNlo", theWeight, 0., "weight");
  }

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
    treeFill.setEcalBarrelClusterShapes(ecalBarrelClusterShapes_);
    treeFill.setEcalEndcapClusterShapes(ecalEndcapClusterShapes_);
    treeFill.setElectronIdProducer(electronIDAssocProducer_);
    treeFill.setTkIsolationProducer(tkIsolationProducer_);
    treeFill.setTowerIsolationProducer(towerIsolationProducer_);
    treeFill.setMatchMap(electronMatchMap_);
    treeFill.saveEleID(true);

    treeFill.writeCollectionToTree(electronCollection_, iEvent, iSetup, prefix, suffix, false);
    if(doMCEleMatch_) {
      treeFill.writeMcIndicesToTree(electronCollection_, iEvent, iSetup, mcTruthCollection_, prefix, suffix, false);
    }

  }



  // fill SC block
  if (dumpSCs_) {

      CmsSuperClusterFiller treeFillHybrid(tree_, 100);
      std::string prefix("");
      std::string hybridSuffix("HybridSCEB");
      treeFillHybrid.writeCollectionToTree(hybridSCCollection_, iEvent, iSetup, prefix, hybridSuffix, false);

      CmsSuperClusterFiller treeFillIsland(tree_, 100);
      std::string islandSuffix("IslandSCEE");
      treeFillIsland.writeCollectionToTree(islandSCCollection_, iEvent, iSetup, prefix, islandSuffix, false);

  }




  // fill muons block
  if(dumpMuons_) {

    CmsCandidateFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Muon");
    treeFill.saveCand(saveCand_);
    treeFill.setMatchMap(muonMatchMap_);

    treeFill.writeCollectionToTree(muonCollection_, iEvent, iSetup, prefix, suffix, false);

    if(doMCMuonMatch_) {
      treeFill.writeMcIndicesToTree(muonCollection_, iEvent, iSetup, mcTruthCollection_, prefix, suffix, false);
    }

  }



  // fill MET block
  if(dumpMet_) {

    CmsCandidateFiller treeRecoFill(tree_, true);
    std::string prefix("");
    std::string suffix("Met");
    treeRecoFill.saveCand(saveCand_);

    treeRecoFill.writeCollectionToTree(metCollection_, iEvent, iSetup, prefix, suffix, false);

    // dump generated MET
    if(dumpGenMet_) {

      CmsCandidateFiller treeGenFill(tree_, true);
      std::string suffix("GenMet");
      treeGenFill.writeCollectionToTree(genMetCollection_, iEvent, iSetup, prefix, suffix, false);

    }

  }



  // fill JET block
  if(dumpJets_) {

    CmsJetFiller treeRecoFill(tree_, jetVertexAlphaCollection_, true);
    std::string prefix("");
    std::string suffix("Jet");
    treeRecoFill.saveCand(saveCand_);
    treeRecoFill.saveJetExtras(saveJetAlpha_);

    treeRecoFill.writeCollectionToTree(jetCollection_, iEvent, iSetup, prefix, suffix, false);

    // dump generated JETs
    if(dumpGenJets_) {

      CmsJetFiller treeGenFill(tree_, jetVertexAlphaCollection_, true);
      std::string suffix("GenJet");
      treeGenFill.saveJetExtras(false);
      treeGenFill.writeCollectionToTree(genJetCollection_, iEvent, iSetup, prefix, suffix, false);

    }

  }
  


  // dump infos on MC production 
  if (dumpGenInfo_) {

      CmsGenInfoFiller treeFill(tree_);
      Handle<int> genProcessID;
      iEvent.getByLabel( "genEventProcID", genProcessID );
      double processID = *genProcessID;
      
      Handle<double> genEventScale;
      iEvent.getByLabel( "genEventScale", genEventScale );
      double pthat = *genEventScale;

      double filter_eff = -99.;
      double cross_section = -99.;
      
      if (processID != 4) {

	  Handle<double> genFilterEff;
	  iEvent.getByLabel( "genEventRunInfo", "FilterEfficiency", genFilterEff);
	  filter_eff = *genFilterEff;
	  
	  Handle<double> genCrossSect;
	  iEvent.getByLabel( "genEventRunInfo", "PreCalculatedCrossSection", genCrossSect); 
	  cross_section = *genCrossSect;

      } else { //ALPGEN case
	  Handle<int> alpgenId;
	  iEvent.getByLabel (genWeightCollection_, alpgenId);
	  processID = *alpgenId;
      }
      
      //       Handle< HepMCProduct > mc;
      //       iEvent.getByLabel( hepMcCollection_, mc );
      //       Handle< GenInfoProduct > gi;
      //       iEvent.getRun().getByLabel( genInfoCollection_, gi);

      Handle< double> weightHandle;
      iEvent.getByLabel (genWeightCollection_, weightHandle);
      double weight = * weightHandle;
      treeFill.writeGenInfoToTree(processID,pthat,filter_eff, cross_section, weight);
  }

 
  if(dumpTree_) tree_->dumpData();

}



// ------------ method called once each job just before starting event loop  ------------
void HWWTreeDumper::beginJob(const edm::EventSetup&) {
  
  fileOut_ = TFile::Open(nameFile_.c_str(), "RECREATE");

  tree_  = new  CmsTree(nameTree_.c_str(),nameTree_.c_str());

  treeConditions_ = new CmsTree("Conditions","Conditions");

  jevt_ = 1;

}



// ------------ method called once each job just after ending the event loop  ------------
void  HWWTreeDumper::endJob() {

  fileOut_->cd();

  TTree* treeEventsOut = tree_->getTree();
  treeEventsOut->Write();

  TTree* treeConditionsOut = treeConditions_->getTree();
  treeConditionsOut->Write();

  fileOut_->Close();

}


