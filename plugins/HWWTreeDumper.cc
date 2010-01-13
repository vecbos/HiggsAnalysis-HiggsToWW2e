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

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMuonFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsElectronFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFlowElectronFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsSuperClusterFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsBasicClusterFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGenInfoFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConditionsFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTrackFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGsfTrackFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsVertexFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCaloTowerFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsV0CandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTriggerTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsRunInfoFiller.h"
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
  saveTrackDeDx_  = iConfig.getUntrackedParameter<bool>("saveTrackDeDx", false);
  saveFatEcal_    = iConfig.getUntrackedParameter<bool>("saveFatEcal", false);
  saveFatHcal_    = iConfig.getUntrackedParameter<bool>("saveFatHcal", false);
  saveFatDT_      = iConfig.getUntrackedParameter<bool>("saveFatDT", false);
  saveFatCSC_     = iConfig.getUntrackedParameter<bool>("saveFatCSC", false);
  saveFatRPC_     = iConfig.getUntrackedParameter<bool>("saveFatRPC", false);
  saveJetAlpha_   = iConfig.getUntrackedParameter<bool>("saveJetAlpha", false);
  saveJet1BTag_    = iConfig.getUntrackedParameter<bool>("saveJet1BTag", false);
  saveJet2BTag_    = iConfig.getUntrackedParameter<bool>("saveJet2BTag", false);

  //electron pflow
  savePFEleTrk_   = iConfig.getUntrackedParameter<bool>("savePFEleGsfTrk", false);
  savePFEleBasic_ = iConfig.getUntrackedParameter<bool>("savePFEleBasic", false);

  // particle identification
  saveEleID_    = iConfig.getUntrackedParameter<bool>("saveEleID", false);

  // basic kinematic informations
  saveCand_  = iConfig.getUntrackedParameter<bool>("saveCand", true);
  
  // Candidate Collections
  dumpPreselInfo_     = iConfig.getUntrackedParameter<bool>("dumpPreselInfo", false);
  dumpSignalKfactor_  = iConfig.getUntrackedParameter<bool>("dumpSignalKfactor", false);
  dumpGenInfo_        = iConfig.getUntrackedParameter<bool>("dumpGenInfo", false);
  dumpElectrons_      = iConfig.getUntrackedParameter<bool>("dumpElectrons", false);
  dumpPFlowElectrons_ = iConfig.getUntrackedParameter<bool>("dumpPFlowElectrons", false);
  dumpSCs_            = iConfig.getUntrackedParameter<bool>("dumpSCs", false);
  dumpBCs_            = iConfig.getUntrackedParameter<bool>("dumpBCs", false);
  dumpTracks_         = iConfig.getUntrackedParameter<bool>("dumpTracks", false);
  dumpGsfTracks_      = iConfig.getUntrackedParameter<bool>("dumpGsfTracks", false);
  dumpMuons_          = iConfig.getUntrackedParameter<bool>("dumpMuons", false);
  dumpJets_           = iConfig.getUntrackedParameter<bool>("dumpJets", false);
  dumpGenJets_        = iConfig.getUntrackedParameter<bool>("dumpGenJets", false);
  dumpMet_            = iConfig.getUntrackedParameter<bool>("dumpMet", false);
  dumpGenMet_         = iConfig.getUntrackedParameter<bool>("dumpGenMet", false);
  dumpVertices_       = iConfig.getUntrackedParameter<bool>("dumpVertices", false);
  dumpK0s_            = iConfig.getUntrackedParameter<bool>("dumpK0s", false);
  dumpCaloTowers_     = iConfig.getUntrackedParameter<bool>("dumpCaloTowers", false);

  // Particle Flow objects
  dumpParticleFlowObjects_ = iConfig.getUntrackedParameter<bool>("dumpParticleFlowObjects",false);
  
  // data run informations
  dumpRunInfo_ = iConfig.getUntrackedParameter<bool>("dumpRunInfo",false);

  // jet vertex collections
  jetVertexAlphaCollection1_ = iConfig.getParameter<edm::InputTag>("jetVertexAlphaCollection1");
  jetVertexAlphaCollection2_ = iConfig.getParameter<edm::InputTag>("jetVertexAlphaCollection2");

  electronCollection_      = iConfig.getParameter<edm::InputTag>("electronCollection");
  pflowElectronCollection_ = iConfig.getParameter<edm::InputTag>("pflowElectronCollection");
  muonCollection_          = iConfig.getParameter<edm::InputTag>("muonCollection");
  ecalSCCollection_        = iConfig.getParameter<edm::InputTag>("ecalSCCollection");
  ecalBarrelSCCollection_  = iConfig.getParameter<edm::InputTag>("ecalBarrelSCCollection");
  ecalEndcapSCCollection_  = iConfig.getParameter<edm::InputTag>("ecalEndcapSCCollection");
  ecalPFClusterCollection_ = iConfig.getParameter<edm::InputTag>("ecalPFClusterCollection");
  ecalBCCollection_        = iConfig.getParameter<edm::InputTag>("ecalBCCollection");
  ecalBarrelRecHits_       = iConfig.getParameter<edm::InputTag>("ecalBarrelRecHits");
  ecalEndcapRecHits_       = iConfig.getParameter<edm::InputTag>("ecalEndcapRecHits");
  calotowersForIsolationProducer_ = iConfig.getParameter<edm::InputTag>("calotowersForIsolationProducer");
  trackCollection_         = iConfig.getParameter<edm::InputTag>("trackCollection");
  refittedForDeDxTrackCollection_ = iConfig.getParameter<edm::InputTag>("refittedForDeDxTrackCollection");
  gsfTrackCollection_      = iConfig.getParameter<edm::InputTag>("gsfTrackCollection");
  vertexCollection_        = iConfig.getParameter<edm::InputTag>("vertexCollection");
  K0sCollection_           = iConfig.getParameter<edm::InputTag>("K0sCollection");
  genJetCollection_        = iConfig.getParameter<edm::InputTag>("genJetCollection");
  jetCollection1_          = iConfig.getParameter<edm::InputTag>("jetCollection1");
  jetCollection2_          = iConfig.getParameter<edm::InputTag>("jetCollection2");
  PFjetCollection1_        = iConfig.getParameter<edm::InputTag>("PFjetCollection1");
  PFjetCollection2_        = iConfig.getParameter<edm::InputTag>("PFjetCollection2");
  metCollection_           = iConfig.getParameter<edm::InputTag>("metCollection");
  TCmetCollection_         = iConfig.getParameter<edm::InputTag>("TCmetCollection");
  PFmetCollection_         = iConfig.getParameter<edm::InputTag>("PFmetCollection");
  genMetCollection_        = iConfig.getParameter<edm::InputTag>("genMetCollection");
  mcTruthCollection_       = iConfig.getParameter<edm::InputTag>("mcTruthCollection");
  electronMatchMap_        = iConfig.getParameter<edm::InputTag>("electronMatchMap");
  muonMatchMap_            = iConfig.getParameter<edm::InputTag>("muonMatchMap");
  hepMcCollection_         = iConfig.getParameter<edm::InputTag>("hepMcCollection");
  genInfoCollection_       = iConfig.getParameter<edm::InputTag>("genInfoCollection");
  genWeightCollection_     = iConfig.getUntrackedParameter<std::string>("genWeightCollection");

  // calotowers collections
  calotowerCollection_ = iConfig.getParameter<edm::InputTag>("calotowerCollection");
  hbheLabel_  = iConfig.getParameter<edm::InputTag>("hbheInput");
  hoLabel_    = iConfig.getParameter<edm::InputTag>("hoInput");
  hfLabel_    = iConfig.getParameter<edm::InputTag>("hfInput");
  ecalLabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("ecalInputs");

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
    treeFill.setEcalBarrelSuperClusters(ecalBarrelSCCollection_);
    treeFill.setEcalEndcapSuperClusters(ecalEndcapSCCollection_);
    treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
    treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
    // for custom isolation
    treeFill.setTracksProducer(trackCollection_);
    treeFill.setCalotowersProducer(calotowersForIsolationProducer_);
    treeFill.setMatchMap(electronMatchMap_);
    treeFill.saveEleID(true);
    treeFill.savePFextra(savePFEleTrk_);

    treeFill.writeCollectionToTree(electronCollection_, iEvent, iSetup, prefix, suffix, false);
    if(doMCEleMatch_) {
      treeFill.writeMcIndicesToTree(electronCollection_, iEvent, iSetup, mcTruthCollection_, prefix, suffix, false);
    }

  }

  if(dumpPFlowElectrons_) {
    CmsPFlowElectronFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("PFEle");
    treeFill.saveCand(saveCand_);
    treeFill.savePFEleBasic(savePFEleBasic_);
    treeFill.savePFEleTrk(savePFEleTrk_);
    treeFill.writeCollectionToTree(pflowElectronCollection_, iEvent, iSetup, prefix, suffix, false);
  }


  // fill SC block
  if (dumpSCs_) {

      CmsSuperClusterFiller treeFill(tree_, 1000);
      std::string prefix("");
      std::string suffix("SC");
      treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
      treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
      treeFill.setTracks(trackCollection_);
      treeFill.setGsfTracks(gsfTrackCollection_);
      treeFill.setCalotowers(calotowersForIsolationProducer_);
      treeFill.doTrackBwdPropagation(true);
      treeFill.writeCollectionToTree(ecalSCCollection_, iEvent, iSetup, prefix, suffix, false);

      CmsSuperClusterFiller treeFillPF(tree_, 1000);
      suffix = "PFSC";
      treeFillPF.setEcalBarrelRecHits(ecalBarrelRecHits_);
      treeFillPF.setEcalEndcapRecHits(ecalEndcapRecHits_);
      treeFillPF.setTracks(trackCollection_);
      treeFillPF.setGsfTracks(gsfTrackCollection_);
      treeFillPF.setCalotowers(calotowersForIsolationProducer_);
      treeFillPF.doTrackBwdPropagation(false);
      treeFillPF.writeCollectionToTree(ecalPFClusterCollection_, iEvent, iSetup, prefix, suffix, false);

  }

  // fill BC block
  if (dumpBCs_) {

      CmsBasicClusterFiller treeFill(tree_, 100);
      std::string prefix("");
      std::string barrelSuffix("BC");
      treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
      treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
      treeFill.setEcalSuperClusters(ecalSCCollection_);
      treeFill.writeCollectionToTree(ecalBCCollection_, iEvent, iSetup, prefix, barrelSuffix, false);

  }

  // fill track block
  if(dumpTracks_) {

    CmsTrackFiller treeFiller(tree_, vertexCollection_, true);
    treeFiller.saveFatTrk(saveFatTrk_);
    treeFiller.setRefittedTracksForDeDxProducer(refittedForDeDxTrackCollection_);
    treeFiller.saveDeDx(saveTrackDeDx_);

    treeFiller.findPrimaryVertex(iEvent);
    treeFiller.saveVtxTrk(true);

    std::string prefix("");
    std::string suffix("Track");

    treeFiller.writeCollectionToTree(trackCollection_, iEvent, iSetup, prefix, suffix, false);

  }

  if(dumpGsfTracks_) {

    CmsGsfTrackFiller treeFiller(tree_, vertexCollection_, true);
    treeFiller.saveFatTrk(saveFatTrk_);
    treeFiller.setRefittedTracksForDeDxProducer(gsfTrackCollection_);
    treeFiller.saveDeDx(false);
    treeFiller.saveVtxTrk(false);

    std::string prefix("");
    std::string suffix("GsfTrack");

    treeFiller.writeCollectionToTree(gsfTrackCollection_, iEvent, iSetup, prefix, suffix, false);

  }

  //fill Primary Vertex and associated tracks
  if(dumpVertices_){
    CmsVertexFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("PV");
    treeFill.writeCollectionToTree(vertexCollection_, iEvent, iSetup, prefix, suffix);
  }

  //fill V0 candidates and associated daughter tracks indices
  if(dumpK0s_){
    CmsV0CandidateFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("K0s");
    treeFill.saveCand(saveCand_);
    treeFill.writeCollectionToTree(K0sCollection_, iEvent, iSetup, prefix, suffix);
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

  // fill CaloTower block
  if(dumpCaloTowers_){
    CmsCaloTowerFiller treeFill(tree_, hbheLabel_, hoLabel_, hfLabel_, ecalLabels_, true);
    std::string prefix("");
    std::string suffix("CaloTowers");
    treeFill.saveCand(dumpCaloTowers_);
    treeFill.saveCaloTowerExtras(dumpCaloTowers_);
    treeFill.writeCollectionToTree(calotowerCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  // fill MET block
  if(dumpMet_) {

    // Calo MET
    CmsCandidateFiller treeRecoFill1(tree_, true);
    std::string prefix("");
    std::string suffix("Met");
    treeRecoFill1.saveCand(saveCand_);
    treeRecoFill1.writeCollectionToTree(metCollection_, iEvent, iSetup, prefix, suffix, false);

    // Track-Corrected MET
    CmsCandidateFiller treeRecoFill2(tree_, true);
    suffix = "TCMet";
    treeRecoFill2.saveCand(saveCand_);
    treeRecoFill2.writeCollectionToTree(TCmetCollection_, iEvent, iSetup, prefix, suffix, false);

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

//     CmsJetFiller treeRecoFill1(tree_, jetVertexAlphaCollection1_, true);
    std::string prefix("");
    std::string suffix("SisConeCorrJet");
//     treeRecoFill1.saveCand(saveCand_);
//     treeRecoFill1.saveJetExtras(saveJetAlpha_);
//     treeRecoFill1.saveJetBTag(saveJet1BTag_);
//     treeRecoFill1.writeCollectionToTree(jetCollection1_, iEvent, iSetup, prefix, suffix, false);


    CmsJetFiller treeRecoFill2(tree_, jetVertexAlphaCollection2_, true);
    suffix = "SisConeJet";
    treeRecoFill2.saveCand(saveCand_);
    treeRecoFill2.saveJetExtras(saveJetAlpha_);
    treeRecoFill2.saveJetBTag(saveJet2BTag_);
    treeRecoFill2.writeCollectionToTree(jetCollection2_, iEvent, iSetup, prefix, suffix, false);


    // particle flow jets
    if ( dumpParticleFlowObjects_ ) {  
//       CmsPFJetFiller pfJetFiller1(tree_, true);
//       suffix = "SisConePFCorrJet";
//       pfJetFiller1.saveCand(saveCand_);
//       pfJetFiller1.writeCollectionToTree(PFjetCollection1_, iEvent, iSetup, prefix, suffix, false);

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
      treeGenFill.saveJetBTag(false);
      treeGenFill.writeCollectionToTree(genJetCollection_, iEvent, iSetup, prefix, suffix, false);

    }

  }
  


  // dump infos on MC production 
  if (dumpGenInfo_) {

    Handle<GenEventInfoProduct> gei;
    iEvent.getByLabel( "generator", gei );

    CmsGenInfoFiller treeFill(tree_);
    treeFill.writeGenInfoToTree( gei );

  }

 
  if(dumpTree_) tree_->dumpData();

}



// ------------ method called once each job just before starting event loop  ------------
void HWWTreeDumper::beginJob(const edm::EventSetup&) {
  
  fileOut_ = TFile::Open(nameFile_.c_str(), "RECREATE");

  tree_  = new  CmsTree(nameTree_.c_str(),nameTree_.c_str());

  treeConditions_ = new TTree("Conditions","Conditions");

  jevt_ = 1;

}



// ------------ method called once each job just after ending the event loop  ------------
void  HWWTreeDumper::endJob() {

  fileOut_->cd();

  TTree* treeEventsOut = tree_->getTree();
  treeEventsOut->Write();

  treeConditions_->Write();

  fileOut_->Close();

}


