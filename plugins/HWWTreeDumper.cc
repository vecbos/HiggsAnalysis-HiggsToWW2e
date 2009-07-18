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

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "RecoBTag/MCTools/interface/JetFlavourIdentifier.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMuonFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsElectronFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsSuperClusterFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGenInfoFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConditionsFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTrackFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsVertexFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTriggerTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsRunInfoFiller.h"
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
  saveTrk_        = iConfig.getUntrackedParameter<bool>("saveTrk", false);
  saveEcal_       = iConfig.getUntrackedParameter<bool>("saveEcal", false);
  saveHcal_       = iConfig.getUntrackedParameter<bool>("saveHcal", false);
  saveDT_         = iConfig.getUntrackedParameter<bool>("saveDT", false);
  saveCSC_        = iConfig.getUntrackedParameter<bool>("saveCSC", false);
  saveRPC_        = iConfig.getUntrackedParameter<bool>("saveRPC", false);
  saveFatTrk_     = iConfig.getUntrackedParameter<bool>("saveFatTrk", false);
  saveFatEcal_    = iConfig.getUntrackedParameter<bool>("saveFatEcal", false);
  saveFatHcal_    = iConfig.getUntrackedParameter<bool>("saveFatHcal", false);
  saveFatDT_      = iConfig.getUntrackedParameter<bool>("saveFatDT", false);
  saveFatCSC_     = iConfig.getUntrackedParameter<bool>("saveFatCSC", false);
  saveFatRPC_     = iConfig.getUntrackedParameter<bool>("saveFatRPC", false);
  saveJetAlpha_   = iConfig.getUntrackedParameter<bool>("saveJetAlpha", false);
  saveJetFlavour_ = iConfig.getUntrackedParameter<bool>("saveJetFlavour", false);
  saveJet1BTag_    = iConfig.getUntrackedParameter<bool>("saveJet1BTag", false);
  saveJet2BTag_    = iConfig.getUntrackedParameter<bool>("saveJet2BTag", false);

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
  dumpTracks_         = iConfig.getUntrackedParameter<bool>("dumpTracks", false);
  dumpMuons_          = iConfig.getUntrackedParameter<bool>("dumpMuons", false);
  dumpJets_           = iConfig.getUntrackedParameter<bool>("dumpJets", false);
  dumpGenJets_        = iConfig.getUntrackedParameter<bool>("dumpGenJets", false);
  dumpMet_            = iConfig.getUntrackedParameter<bool>("dumpMet", false);
  dumpGenMet_         = iConfig.getUntrackedParameter<bool>("dumpGenMet", false);
  dumpVertices_       = iConfig.getUntrackedParameter<bool>("dumpVertices", false);

  // Particle Flow objects
  dumpParticleFlowObjects_ = iConfig.getUntrackedParameter<bool>("dumpParticleFlowObjects",false);
  
  // data run informations
  dumpRunInfo_ = iConfig.getUntrackedParameter<bool>("dumpRunInfo",false);

  // jet vertex collections
  jetVertexAlphaCollection1_ = iConfig.getParameter<edm::InputTag>("jetVertexAlphaCollection1");
  jetVertexAlphaCollection2_ = iConfig.getParameter<edm::InputTag>("jetVertexAlphaCollection2");
  jetMCFlavourIdentifier_    = iConfig.getParameter<edm::ParameterSet>("jetIdParameters");

  electronCollection_      = iConfig.getParameter<edm::InputTag>("electronCollection");
  muonCollection_          = iConfig.getParameter<edm::InputTag>("muonCollection");
  ecalBarrelSCCollection_  = iConfig.getParameter<edm::InputTag>("ecalBarrelSCCollection");
  ecalEndcapSCCollection_  = iConfig.getParameter<edm::InputTag>("ecalEndcapSCCollection");
  ecalBarrelRecHits_       = iConfig.getParameter<edm::InputTag>("ecalBarrelRecHits");
  ecalEndcapRecHits_       = iConfig.getParameter<edm::InputTag>("ecalEndcapRecHits");
  tracksForIsolationProducer_     = iConfig.getParameter<edm::InputTag>("tracksForIsolationProducer");
  calotowersForIsolationProducer_ = iConfig.getParameter<edm::InputTag>("calotowersForIsolationProducer");
  trackCollection_    = iConfig.getParameter<edm::InputTag>("trackCollection");
  vertexCollection_        = iConfig.getParameter<edm::InputTag>("vertexCollection");
  genJetCollection_       = iConfig.getParameter<edm::InputTag>("genJetCollection");
  jetCollection1_          = iConfig.getParameter<edm::InputTag>("jetCollection1");
  jetCollection2_          = iConfig.getParameter<edm::InputTag>("jetCollection2");
  PFjetCollection1_        = iConfig.getParameter<edm::InputTag>("PFjetCollection1");
  PFjetCollection2_        = iConfig.getParameter<edm::InputTag>("PFjetCollection2");
  metCollection_           = iConfig.getParameter<edm::InputTag>("metCollection");
  PFmetCollection_         = iConfig.getParameter<edm::InputTag>("PFmetCollection");
  genMetCollection_        = iConfig.getParameter<edm::InputTag>("genMetCollection");
  mcTruthCollection_       = iConfig.getParameter<edm::InputTag>("mcTruthCollection");
  electronMatchMap_        = iConfig.getParameter<edm::InputTag>("electronMatchMap");
  muonMatchMap_            = iConfig.getParameter<edm::InputTag>("muonMatchMap");
  hepMcCollection_         = iConfig.getParameter<edm::InputTag>("hepMcCollection");
  genInfoCollection_       = iConfig.getParameter<edm::InputTag>("genInfoCollection");
  genWeightCollection_     = iConfig.getUntrackedParameter<std::string>("genWeightCollection");

  // trigger Collections
  triggerInputTag_     = iConfig.getParameter<edm::InputTag>("TriggerResultsTag") ;
  dumpTriggerResults_  = iConfig.getUntrackedParameter<bool>("dumpTriggerResults");
}



HWWTreeDumper::~HWWTreeDumper() { }



// ------------ method called to for each event  ------------
void HWWTreeDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  /// fill the run info (run number, event, ...)
  if(dumpRunInfo_) {
    CmsRunInfoFiller runFiller( tree_ );
    runFiller.writeRunInfoToTree(iEvent,iSetup,false);
  }

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
    treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
    treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
    // for custom isolation
    treeFill.setTracksProducer(tracksForIsolationProducer_);
    treeFill.setCalotowersProducer(calotowersForIsolationProducer_);
    treeFill.setMatchMap(electronMatchMap_);
    treeFill.saveEleID(true);

    treeFill.writeCollectionToTree(electronCollection_, iEvent, iSetup, prefix, suffix, false);
    if(doMCEleMatch_) {
      treeFill.writeMcIndicesToTree(electronCollection_, iEvent, iSetup, mcTruthCollection_, prefix, suffix, false);
    }

  }



  // fill SC block
  if (dumpSCs_) {

      CmsSuperClusterFiller treeFillBarrel(tree_, 100);
      std::string prefix("");
      std::string barrelSuffix("SCEB");
      treeFillBarrel.writeCollectionToTree(ecalBarrelSCCollection_, iEvent, iSetup, prefix, barrelSuffix, false);

      CmsSuperClusterFiller treeFillEndcap(tree_, 100);
      std::string endcapSuffix("SCEE");
      treeFillEndcap.writeCollectionToTree(ecalEndcapSCCollection_, iEvent, iSetup, prefix, endcapSuffix, false);

  }

  // fill track block
  if(dumpTracks_) {

    CmsTrackFiller treeFiller(tree_, vertexCollection_, true);
    treeFiller.saveFatTrk(saveFatTrk_);

    treeFiller.findPrimaryVertex(iEvent);
    treeFiller.saveVtxTrk(true);

    std::string prefix("");
    std::string suffix("Track");
    treeFiller.saveCand(saveCand_);

    treeFiller.writeCollectionToTree(trackCollection_, iEvent, iSetup, prefix, suffix, false);

  }

  //fill Primary Vertex and associated tracks
  if(dumpVertices_){
    CmsVertexFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("PV");
    treeFill.writeCollectionToTree(vertexCollection_, iEvent, iSetup, prefix, suffix);
  }

  // fill muons block
  if(dumpMuons_) {

    CmsMuonFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Muon");
    treeFill.saveCand(saveCand_);
    treeFill.saveFatTrk(saveFatTrk_);
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

    // particle flow met
    if ( dumpParticleFlowObjects_ ) {
      CmsCandidateFiller pfMetFiller(tree_, true);
      suffix = "PFMet";
      pfMetFiller.saveCand(saveCand_);
      pfMetFiller.writeCollectionToTree(PFmetCollection_, iEvent, iSetup, prefix, suffix, false);
    }

    // dump generated MET
    if(dumpGenMet_) {

      CmsCandidateFiller treeGenFill(tree_, true);
      std::string suffix("GenMet");
      treeGenFill.writeCollectionToTree(genMetCollection_, iEvent, iSetup, prefix, suffix, false);

    }

  }



  // fill JET block
  if(dumpJets_) {

    CmsJetFiller treeRecoFill1(tree_, jetVertexAlphaCollection1_, true);
    std::string prefix("");
    std::string suffix("SisConeCorrJet");
    treeRecoFill1.saveCand(saveCand_);
    treeRecoFill1.saveJetExtras(saveJetAlpha_);
    treeRecoFill1.saveJetFlavour(saveJetFlavour_);
    treeRecoFill1.saveJetBTag(saveJet1BTag_);
    if(saveJetFlavour_) { 
      JetFlavourIdentifier jetMCFlavourIdentifier(jetMCFlavourIdentifier_);
      treeRecoFill1.setJetFlavour(jetMCFlavourIdentifier);
    }
    treeRecoFill1.writeCollectionToTree(jetCollection1_, iEvent, iSetup, prefix, suffix, false);


    CmsJetFiller treeRecoFill2(tree_, jetVertexAlphaCollection2_, true);
    suffix = "SisConeJet";
    treeRecoFill2.saveCand(saveCand_);
    treeRecoFill2.saveJetExtras(saveJetAlpha_);
    treeRecoFill2.saveJetFlavour(saveJetFlavour_);
    treeRecoFill2.saveJetBTag(saveJet2BTag_);
    if(saveJetFlavour_) { 
      JetFlavourIdentifier jetMCFlavourIdentifier(jetMCFlavourIdentifier_);
      treeRecoFill2.setJetFlavour(jetMCFlavourIdentifier);
    }
    treeRecoFill2.writeCollectionToTree(jetCollection2_, iEvent, iSetup, prefix, suffix, false);


    // particle flow jets
    if ( dumpParticleFlowObjects_ ) {  
      CmsPFJetFiller pfJetFiller1(tree_, true);
      suffix = "SisConePFCorrJet";
      pfJetFiller1.saveCand(saveCand_);
      pfJetFiller1.writeCollectionToTree(PFjetCollection1_, iEvent, iSetup, prefix, suffix, false);

      CmsPFJetFiller pfJetFiller2(tree_, true);
      suffix = "SisConePFJet";
      pfJetFiller2.saveCand(saveCand_);
      pfJetFiller2.writeCollectionToTree(PFjetCollection2_, iEvent, iSetup, prefix, suffix, false);
    }

    // dump generated JETs
    if(dumpGenJets_) {

      CmsJetFiller treeGenFill(tree_, jetVertexAlphaCollection1_, true);
      suffix = "SisConeGenJet";
      treeGenFill.saveJetExtras(false);
      treeGenFill.saveJetFlavour(saveJetFlavour_);
      treeGenFill.saveJetBTag(false);
      if(saveJetFlavour_) { 
	JetFlavourIdentifier jetMCFlavourIdentifier(jetMCFlavourIdentifier_);
	treeGenFill.setJetFlavour(jetMCFlavourIdentifier);
      }
      treeGenFill.writeCollectionToTree(genJetCollection_, iEvent, iSetup, prefix, suffix, false);

    }

  }
  


  // dump infos on MC production 
  if (dumpGenInfo_) {

    Handle<double> genEventScale;
    iEvent.getByLabel( "genEventScale", genEventScale );
    double pthat = *genEventScale;

    Handle<double> genEventWeight;
    iEvent.getByLabel( "genEventWeight", genEventWeight );
    double weight = *genEventWeight;

    CmsGenInfoFiller treeFill(tree_);
    treeFill.writeGenInfoToTree(-1.0, pthat, -1.0, -1.0, weight);

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


