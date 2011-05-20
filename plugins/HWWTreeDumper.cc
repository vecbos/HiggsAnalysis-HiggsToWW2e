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
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFTauFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFPreIdFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPhotonFiller.h"
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
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsJPTJetFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCaloTowerFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsV0CandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTriggerTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsRunInfoFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsHcalNoiseFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMetFiller.h"
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
  saveJetBTag_    = iConfig.getUntrackedParameter<bool>("saveJetBTag", false);

  //electron pflow
  savePFEleBasic_  = iConfig.getUntrackedParameter<bool>("savePFEleBasic",  true);
  savePFEleIsoDep_ = iConfig.getUntrackedParameter<bool>("savePFEleIsoDep", true);

  //tau pflow                                                                                                                                                                                                                              
  savePFTauBasic_  = iConfig.getUntrackedParameter<bool>("savePFTauBasic", false);
  saveLeadPFCand_  = iConfig.getUntrackedParameter<bool>("saveLeadPFCand", false);
  savePFTauDiscriminators_ = iConfig.getUntrackedParameter<bool>("savePFTauDiscriminators", false);

  // particle identification
  saveEleID_    = iConfig.getUntrackedParameter<bool>("saveEleID", false);

  // basic kinematic informations
  saveCand_  = iConfig.getUntrackedParameter<bool>("saveCand", true);
  
  // Candidate Collections
  dumpPreselInfo_     = iConfig.getUntrackedParameter<bool>("dumpPreselInfo", false);
  dumpSignalKfactor_  = iConfig.getUntrackedParameter<bool>("dumpSignalKfactor", false);
  dumpGenInfo_        = iConfig.getUntrackedParameter<bool>("dumpGenInfo", false);
  dumpElectrons_      = iConfig.getUntrackedParameter<bool>("dumpElectrons", false);
  dumpPhotons_        = iConfig.getUntrackedParameter<bool>("dumpPhotons", false);
  dumpPFlowElectrons_ = iConfig.getUntrackedParameter<bool>("dumpPFlowElectrons", false);
  dumpSCs_            = iConfig.getUntrackedParameter<bool>("dumpSCs", false);
  dumpBCs_            = iConfig.getUntrackedParameter<bool>("dumpBCs", false);
  dumpTracks_         = iConfig.getUntrackedParameter<bool>("dumpTracks", false);
  dumpGsfTracks_      = iConfig.getUntrackedParameter<bool>("dumpGsfTracks", false);
  dumpMuonTracks_     = iConfig.getUntrackedParameter<bool>("dumpMuonTracks", false);
  dumpMuons_          = iConfig.getUntrackedParameter<bool>("dumpMuons", false);
  dumpPFTaus_         = iConfig.getUntrackedParameter<bool>("dumpPFTaus", false);
  dumpJets_           = iConfig.getUntrackedParameter<bool>("dumpJets", false);
  dumpGenJets_        = iConfig.getUntrackedParameter<bool>("dumpGenJets", false);
  dumpPUcorrPFJet_    = iConfig.getUntrackedParameter<bool>("dumpPUcorrPFJet", false);
  dumpMet_            = iConfig.getUntrackedParameter<bool>("dumpMet", false);
  dumpGenMet_         = iConfig.getUntrackedParameter<bool>("dumpGenMet", false);
  dumpVertices_       = iConfig.getUntrackedParameter<bool>("dumpVertices", false);
  dumpK0s_            = iConfig.getUntrackedParameter<bool>("dumpK0s", false);
  dumpCaloTowers_     = iConfig.getUntrackedParameter<bool>("dumpCaloTowers", false);
  dumpHcalNoiseFlags_ = iConfig.getUntrackedParameter<bool>("dumpHcalNoiseFlags", false);
  aodHcalNoiseFlags_  = iConfig.getUntrackedParameter<bool>("AODHcalNoiseFlags", true);

  // Particle Flow objects
  dumpParticleFlowObjects_ = iConfig.getUntrackedParameter<bool>("dumpParticleFlowObjects",false);
  dumpPFpreId_             = iConfig.getUntrackedParameter<bool>("dumpPFpreId", false);  

  // data run informations
  dumpRunInfo_ = iConfig.getUntrackedParameter<bool>("dumpRunInfo",false);

  electronCollection_      = iConfig.getParameter<edm::InputTag>("electronCollection");
  pflowElectronCollection_ = iConfig.getParameter<edm::InputTag>("pflowElectronCollection");
  photonCollection_        = iConfig.getParameter<edm::InputTag>("photonCollection");
  muonCollection_          = iConfig.getParameter<edm::InputTag>("muonCollection");
  pfTauCollection_         = iConfig.getParameter<edm::InputTag>("pfTauCollection");
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
  globalMuonTrackCollection_ = iConfig.getParameter<edm::InputTag>("globalMuonTrackCollection");
  standAloneMuonTrackCollection_ = iConfig.getParameter<edm::InputTag>("standAloneMuonTrackCollection");
  vertexCollection_        = iConfig.getParameter<edm::InputTag>("vertexCollection");
  K0sCollection_           = iConfig.getParameter<edm::InputTag>("K0sCollection");
  genJetCollection_        = iConfig.getParameter<edm::InputTag>("genJetCollection");
  jetCollection1_          = iConfig.getParameter<edm::InputTag>("jetCollection1");
  jetCollection2_          = iConfig.getParameter<edm::InputTag>("jetCollection2");
  PFjetCollection1_        = iConfig.getParameter<edm::InputTag>("PFjetCollection1");
  PFjetCollection2_        = iConfig.getParameter<edm::InputTag>("PFjetCollection2");
  PFpuCorrJetCollection1_  = iConfig.getParameter<edm::InputTag>("PFpuCorrJetCollection1");
  PFpuCorrJetCollection2_  = iConfig.getParameter<edm::InputTag>("PFpuCorrJetCollection2");
  JPTjetCollection1_       = iConfig.getParameter<edm::InputTag>("JPTjetCollection1");
  JPTjetCollection2_       = iConfig.getParameter<edm::InputTag>("JPTjetCollection2");

  // btag collections
  PFJetsBTags_              = iConfig.getUntrackedParameter<edm::ParameterSet>("PFJetsBTags");
  PFPUcorrJetsBTags_        = iConfig.getUntrackedParameter<edm::ParameterSet>("PFPUcorrJetsBTags");

  metCollection_           = iConfig.getParameter<edm::InputTag>("metCollection");
  // corrmetCollection_       = iConfig.getParameter<edm::InputTag>("corrmetCollection");
  TCmetCollection_         = iConfig.getParameter<edm::InputTag>("TCmetCollection");
  PFmetCollection_         = iConfig.getParameter<edm::InputTag>("PFmetCollection");
  genMetCollection_        = iConfig.getParameter<edm::InputTag>("genMetCollection");
  chargedMetCollection_    = iConfig.getParameter<edm::InputTag>("chargedMetCollection");
  mcTruthCollection_       = iConfig.getParameter<edm::InputTag>("mcTruthCollection");
  electronMatchMap_        = iConfig.getParameter<edm::InputTag>("electronMatchMap");
  muonMatchMap_            = iConfig.getParameter<edm::InputTag>("muonMatchMap");
  hepMcCollection_         = iConfig.getParameter<edm::InputTag>("hepMcCollection");
  genInfoCollection_       = iConfig.getParameter<edm::InputTag>("genInfoCollection");
  genWeightCollection_     = iConfig.getUntrackedParameter<std::string>("genWeightCollection");
  PFpreIdCollection_       = iConfig.getParameter<edm::InputTag>("PFpreIdCollection");

  // calotowers collections
  calotowerCollection_ = iConfig.getParameter<edm::InputTag>("calotowerCollection");
  hbheLabel_  = iConfig.getParameter<edm::InputTag>("hbheInput");
  hoLabel_    = iConfig.getParameter<edm::InputTag>("hoInput");
  hfLabel_    = iConfig.getParameter<edm::InputTag>("hfInput");
  ecalLabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("ecalInputs");

  // trigger Collections
  dumpTriggerResults_  = iConfig.getUntrackedParameter<bool>("dumpTriggerResults");
  dumpHLTObject_       = iConfig.getUntrackedParameter<bool>("dumpHLTObjects");
  hltParms_            = iConfig.getUntrackedParameter<edm::ParameterSet>("HLTObjectsInfo");

  // PFTau Discriminators
  tauDiscrByLeadTrackFindingTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByLeadTrackFindingTag");
  tauDiscrByLeadTrackPtCutTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByLeadTrackPtCutTag");
  //  tauDiscrByNProngsTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByNProngsTag");
  tauDiscrByTrackIsoTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTrackIsoTag");
  tauDiscrByEcalIsoTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByEcalIsoTag");
  tauDiscrAgainstMuonsTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrAgainstMuonsTag");
  tauDiscrAgainstElectronsTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrAgainstElectronsTag");
//   tauDiscrByTaNCTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCTag");
//   tauDiscrByTaNCfrHalfPercentTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCfrHalfPercentTag");
//   tauDiscrByTaNCfrOnePercentTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCfrOnePercentTag");
//   tauDiscrByTaNCfrQuarterPercentTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCfrQuarterPercentTag");
//   tauDiscrByTaNCfrTenthPercentTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCfrTenthPercentTag");

  // Hcal collections
  hcalNoiseSummaryLabel_ = iConfig.getParameter<edm::InputTag>("hcalNoiseSummary");
}



HWWTreeDumper::~HWWTreeDumper() { 
  if(hltObjectFiller_) delete(hltObjectFiller_);
}



// ------------ method called to for each event  ------------
void HWWTreeDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  /// fill the run info (run number, event, ...)
  if(dumpRunInfo_) {
    CmsRunInfoFiller runFiller( tree_, dumpMCTruth_ );
    runFiller.writeRunInfoToTree(iEvent,iSetup,false);
  }

  // get MC truth
  CmsMcTruthTreeFiller treeFill(tree_);

  if(dumpMCTruth_) {

    treeFill.writeCollectionToTree( mcTruthCollection_, iEvent, 100 );

  }


  jevtInRun_++;
  // fill the trigger paths info
  if(dumpTriggerResults_) {

    CmsTriggerTreeFiller triggerTreeFill (tree_, hltParms_) ;
    std::string prefix ("") ;
    std::string suffix ("Trg") ;
    triggerTreeFill.writeTriggerToTree (iEvent, prefix, suffix) ;

    /// fill the trigger mask in the Event tree
    std::vector<std::string>* trgNames = new vector<std::string>;
    bool firstEvent = true;
    if(jevt_>1) firstEvent=false; 
    CmsConditionsFiller conditionsFiller( tree_, hltParms_, trgNames );
    conditionsFiller.writeConditionsToTree(iEvent,  firstEvent);
    jevt_++;
  }
  if(dumpHLTObject_) {
    //forward beginRun to HLTObjectFiller
    if(jevtInRun_==1) hltObjectFiller_->beginRun( iEvent, iEvent.getRun(), iSetup);
    hltObjectFiller_->writeHLTObjectToTree(iEvent);
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
    treeFill.setGeneralTracks(trackCollection_);
    treeFill.setEcalBarrelSuperClusters(ecalBarrelSCCollection_);
    treeFill.setEcalEndcapSuperClusters(ecalEndcapSCCollection_);
    treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
    treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
    // for custom isolation
    treeFill.setTracksProducer(trackCollection_);
    treeFill.setCalotowersProducer(calotowersForIsolationProducer_);
    treeFill.setMatchMap(electronMatchMap_);
    treeFill.saveEleID(true);

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
    treeFill.setGeneralTracks(trackCollection_);
    treeFill.savePFEleBasic(savePFEleBasic_);
    treeFill.savePFEleIsoDep(savePFEleIsoDep_);
    treeFill.writeCollectionToTree(pflowElectronCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  if(dumpPFpreId_) {
    CmsPFPreIdFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("PFpreId");
    treeFill.writeCollectionToTree(PFpreIdCollection_, trackCollection_, iEvent, iSetup, prefix, suffix, false);  
  }

  // fill Photons block
  if(dumpPhotons_) {

    CmsPhotonFiller treeFill(tree_, true);
    std::string prefix("");
    std::string suffix("Pho");
    treeFill.saveCand(saveCand_);
    treeFill.setEcalBarrelSuperClusters(ecalBarrelSCCollection_);
    treeFill.setEcalEndcapSuperClusters(ecalEndcapSCCollection_);
    treeFill.writeCollectionToTree(photonCollection_, iEvent, iSetup, prefix, suffix, false);
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
      treeFill.doTrackBwdPropagation(false);
      treeFill.writeCollectionToTree(ecalSCCollection_, iEvent, iSetup, prefix, suffix, false, photonCollection_);

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

    treeFiller.isGsf(false);
    treeFiller.setRefittedTracksForDeDxProducer(refittedForDeDxTrackCollection_);
    treeFiller.saveDeDx(saveTrackDeDx_);

    treeFiller.saveVtxTrk(true);

    std::string prefix("");
    std::string suffix("Track");

    treeFiller.writeCollectionToTree(trackCollection_, iEvent, iSetup, prefix, suffix, false);

  }

  if(dumpGsfTracks_) {

    CmsGsfTrackFiller treeFiller(tree_, vertexCollection_, true);
    treeFiller.saveFatTrk(saveFatTrk_);
    treeFiller.setRefittedTracksForDeDxProducer(gsfTrackCollection_);
    treeFiller.isGsf(true);
    treeFiller.saveDeDx(false);
    treeFiller.saveVtxTrk(true);


    std::string prefix("");
    std::string suffix("GsfTrack");

    treeFiller.writeCollectionToTree(gsfTrackCollection_, iEvent, iSetup, prefix, suffix, false);

  }

  if(dumpMuonTracks_) {

    // dump tracks reconstructed in both tracked and muon detector
    CmsTrackFiller treeFillerGlobalMuonTrack(tree_, vertexCollection_, true);
    treeFillerGlobalMuonTrack.saveFatTrk(saveFatTrk_);
    treeFillerGlobalMuonTrack.setRefittedTracksForDeDxProducer(refittedForDeDxTrackCollection_);
    treeFillerGlobalMuonTrack.saveDeDx(false);
    treeFillerGlobalMuonTrack.isGsf(false);
    treeFillerGlobalMuonTrack.saveVtxTrk(true);

    std::string prefix("");
    std::string suffix1("GlobalMuonTrack");
    treeFillerGlobalMuonTrack.writeCollectionToTree(globalMuonTrackCollection_, iEvent, iSetup, prefix, suffix1, false);

    // dump tracks reconstructed in the muon detector only
    CmsTrackFiller treeFillerSTAMuonTrack(tree_, vertexCollection_, true);
    treeFillerSTAMuonTrack.saveFatTrk(saveFatTrk_);
    treeFillerSTAMuonTrack.setRefittedTracksForDeDxProducer(refittedForDeDxTrackCollection_);
    treeFillerSTAMuonTrack.saveDeDx(false);
    treeFillerSTAMuonTrack.isGsf(false);
    treeFillerSTAMuonTrack.saveVtxTrk(false);

    std::string suffix2("STAMuonTrack");
    treeFillerSTAMuonTrack.writeCollectionToTree(standAloneMuonTrackCollection_, iEvent, iSetup, prefix, suffix2, false);

  }

  //fill Primary Vertex and associated tracks
  if(dumpVertices_){
    CmsVertexFiller treeFillerVertices(tree_, true);
    treeFillerVertices.setChargedMet(chargedMetCollection_);
    std::string prefix("");
    std::string suffix("PV");
    treeFillerVertices.writeCollectionToTree(vertexCollection_, iEvent, iSetup, prefix, suffix);
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

  // fill PFTau block                                                                                                                                                                                                                      
  if(dumpPFTaus_)
    {
      CmsPFTauFiller treeFill(tree_, true);
      std::string prefix("");
      std::string suffix("PFTau");
      treeFill.savePFTauBasic(savePFTauBasic_);
      treeFill.saveLeadPFCand(saveLeadPFCand_);
      treeFill.writeCollectionToTree(pfTauCollection_, iEvent, iSetup, prefix, suffix,
                                     tauDiscrByLeadTrackFindingTag_, tauDiscrByLeadTrackPtCutTag_, //tauDiscrByNProngsTag_,
                                     tauDiscrByTrackIsoTag_, tauDiscrByEcalIsoTag_, tauDiscrAgainstMuonsTag_, tauDiscrAgainstElectronsTag_,
//                                      tauDiscrByTaNCTag_,
//                                      tauDiscrByTaNCfrHalfPercentTag_, tauDiscrByTaNCfrOnePercentTag_,
//                                      tauDiscrByTaNCfrQuarterPercentTag_, tauDiscrByTaNCfrTenthPercentTag_,
                                     false);
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
    CmsMetFiller treeRecoFill1(tree_, true);
    std::string prefix("");
    std::string suffix("Met");
    treeRecoFill1.saveCand(saveCand_);
    treeRecoFill1.writeCollectionToTree(metCollection_, iEvent, iSetup, prefix, suffix, false);

    // Corrected CALO MET
    // CmsCandidateFiller treeRecoFill1bis(tree_, true);
    // suffix = "CorrMet";
    // treeRecoFill1bis.saveCand(saveCand_);
    // treeRecoFill1bis.writeCollectionToTree(corrmetCollection_, iEvent, iSetup, prefix, suffix, false);

    // Track-Corrected MET
    CmsMetFiller treeRecoFill2(tree_, true);
    suffix = "TCMet";
    treeRecoFill2.saveCand(saveCand_);
    treeRecoFill2.writeCollectionToTree(TCmetCollection_, iEvent, iSetup, prefix, suffix, false);

    // particle flow met
    if ( dumpParticleFlowObjects_ ) {
      CmsMetFiller pfMetFiller(tree_, true);
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

    CmsJetFiller caloJetFiller(tree_, true);
    std::string prefix("");
    std::string suffix("AK5Jet");
    caloJetFiller.saveCand(saveCand_);
    caloJetFiller.saveJetExtras(true);
    caloJetFiller.saveJetBTag(saveJetBTag_);
    caloJetFiller.writeCollectionToTree(jetCollection1_, iEvent, iSetup, prefix, suffix, false, jetCollection2_);


    // particle flow jets
    if ( dumpParticleFlowObjects_ ) {  
      CmsPFJetFiller pfJetFiller(tree_, true);
      suffix = "AK5PFJet";
      pfJetFiller.saveCand(saveCand_);
      pfJetFiller.saveJetBTag(false);  // since it is done on the same collection as pfPUcorrJetFiller, do not waste CPU repeating it
      pfJetFiller.setBTags(PFJetsBTags_);
      pfJetFiller.writeCollectionToTree(PFjetCollection1_, iEvent, iSetup, prefix, suffix, false, PFjetCollection2_);

    }

    // particle flow jets with correction for pileup
    if ( dumpParticleFlowObjects_ && dumpPUcorrPFJet_ ) {
      CmsPFJetFiller pfPUcorrJetFiller(tree_, true);
      suffix = "AK5PFPUcorrJet";
      pfPUcorrJetFiller.saveCand(saveCand_);
      pfPUcorrJetFiller.saveJetBTag(saveJetBTag_);
      pfPUcorrJetFiller.setBTags(PFPUcorrJetsBTags_);
      pfPUcorrJetFiller.writeCollectionToTree(PFpuCorrJetCollection1_, iEvent, iSetup, prefix, suffix, false,PFpuCorrJetCollection2_);
    }

    // Jet Plus Tracks jets
    //     CmsJPTJetFiller jptJetFiller(tree_, true);
    //     suffix = "AK5JPTJet";
    //     jptJetFiller.saveCand(saveCand_);
    //     jptJetFiller.saveJetBTag(saveJetBTag_);
    //     jptJetFiller.writeCollectionToTree(JPTjetCollection1_, iEvent, iSetup, prefix, suffix, false, JPTjetCollection2_);

    // dump generated JETs
    if(dumpGenJets_) {

      CmsJetFiller genJetFiller(tree_, true);
      suffix = "AK5GenJet";
      genJetFiller.saveJetExtras(false);
      genJetFiller.saveJetBTag(false);
      genJetFiller.writeCollectionToTree(genJetCollection_, iEvent, iSetup, prefix, suffix, false);

    }

  }
  


  // dump infos on MC production 
  if (dumpGenInfo_) {

    Handle<GenEventInfoProduct> gei;
    iEvent.getByLabel( "generator", gei );

    CmsGenInfoFiller treeFill(tree_);
    treeFill.writeGenInfoToTree( gei );

  }
  
  
  // dump Hcal noise flags
  if(dumpHcalNoiseFlags_) {

    CmsHcalNoiseFiller treeFill(tree_, true);

    treeFill.writeHcalNoiseSummaryToTree(hbheLabel_, hfLabel_, hcalNoiseSummaryLabel_, iEvent, iSetup,
      aodHcalNoiseFlags_);
  }

 
  if(dumpTree_) tree_->dumpData();

}



// ------------ method called once each job just before starting event loop  ------------
void HWWTreeDumper::beginJob() {
  
  fileOut_ = TFile::Open(nameFile_.c_str(), "RECREATE");

  tree_  = new  CmsTree(nameTree_.c_str(),nameTree_.c_str());

  jevt_ = 1;

  //HLTObject Filler needs to exist before beginRun is called
  if(dumpHLTObject_)
    hltObjectFiller_ = new CmsHLTObjectFiller(tree_,hltParms_);

}

void HWWTreeDumper::beginRun( const Run & iRun, const EventSetup & iSetup )
{
  jevtInRun_ = 0;
}



// ------------ method called once each job just after ending the event loop  ------------
void  HWWTreeDumper::endJob() {

  fileOut_->cd();

  TTree* treeEventsOut = tree_->getTree();
  treeEventsOut->Write();

  fileOut_->Close();

}


