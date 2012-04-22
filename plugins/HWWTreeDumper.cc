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
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPdfWeightFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTreeDumper.h"


HWWTreeDumper::HWWTreeDumper(const edm::ParameterSet& iConfig)
{
  
  nameFile_      = iConfig.getUntrackedParameter<std::string>("nameFile", "RootOutput.root");
  nameTree_      = iConfig.getUntrackedParameter<std::string>("nameTree", "BaseTree");
  dumpTree_      = iConfig.getUntrackedParameter<bool>("dumpTree", false);
  dumpMCTruth_   = iConfig.getUntrackedParameter<bool>("dumpMCTruth", false);
  
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
  usePhotonFix_ = iConfig.getUntrackedParameter<bool>("usePhotonFix", true);
  if (usePhotonFix_)
    {
      phFixElePar_  = iConfig.getUntrackedParameter<edm::ParameterSet>("PhotonFix_4_2e");
      phFixPhoPar_  = iConfig.getUntrackedParameter<edm::ParameterSet>("PhotonFix_4_2p");
    }
  useEnergyRegression_ = iConfig.getUntrackedParameter<bool>("useEnergyRegression", true);
  if (useEnergyRegression_)
    {
      energyRegressionElectronFile_ = iConfig.getUntrackedParameter<std::string>("energyRegressionElectronFile");
      energyRegressionPhotonFile_ = iConfig.getUntrackedParameter<std::string>("energyRegressionPhotonFile");
    }

  posCalcParameters_ =
    iConfig.getParameter<edm::ParameterSet>("posCalcParameters");
  
  // Candidate Collections
  dumpPreselInfo_     = iConfig.getUntrackedParameter<bool>("dumpPreselInfo", false);
  dumpSignalKfactor_  = iConfig.getUntrackedParameter<bool>("dumpSignalKfactor", false);
  dumpGenInfo_        = iConfig.getUntrackedParameter<bool>("dumpGenInfo", false);
  dumpLHE_            = iConfig.getUntrackedParameter<bool>("dumpLHE", false);
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
  dumphpsPFTaus_      = iConfig.getUntrackedParameter<bool>("dumphpsPFTaus", false);
  dumphpsTancTaus_    = iConfig.getUntrackedParameter<bool>("dumphpsTancTaus", false);
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

  // eigenvalues for PDF systematics
  dumpPdfWeight_ = iConfig.getUntrackedParameter<bool>("dumpPdfWeight",false);

  electronCollection_      = iConfig.getParameter<edm::InputTag>("electronCollection");
  pflowElectronCollection_ = iConfig.getParameter<edm::InputTag>("pflowElectronCollection");
  photonCollection_        = iConfig.getParameter<edm::InputTag>("photonCollection");
  muonCollection_          = iConfig.getParameter<edm::InputTag>("muonCollection");
  pfTauCollection_         = iConfig.getParameter<edm::InputTag>("pfTauCollection");
  hpspfTauCollection_      = iConfig.getParameter<edm::InputTag>("hpspfTauCollection");
  hpsTancTausCollection_   = iConfig.getParameter<edm::InputTag>("hpsTancTausCollection");
  PFCandidateCollection_   = iConfig.getParameter<edm::InputTag>("PFCandidateCollection");
  ecalSCCollection_        = iConfig.getParameter<edm::InputTag>("ecalSCCollection");
  ecalBarrelSCCollection_  = iConfig.getParameter<edm::InputTag>("ecalBarrelSCCollection");
  ecalEndcapSCCollection_  = iConfig.getParameter<edm::InputTag>("ecalEndcapSCCollection");
  ecalPFClusterCollection_ = iConfig.getParameter<edm::InputTag>("ecalPFClusterCollection");
  ecalBCCollection_        = iConfig.getParameter<edm::InputTag>("ecalBCCollection");
  ecalBarrelRecHits_       = iConfig.getParameter<edm::InputTag>("ecalBarrelRecHits");
  ecalEndcapRecHits_       = iConfig.getParameter<edm::InputTag>("ecalEndcapRecHits");
  esRecHits_               = iConfig.getParameter<edm::InputTag>("esRecHits");
  calotowersForIsolationProducer_ = iConfig.getParameter<edm::InputTag>("calotowersForIsolationProducer");
  conversions_             = iConfig.getParameter<edm::InputTag>("conversionCollection");
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
  jetCollection3_          = iConfig.getParameter<edm::InputTag>("jetCollection3");
  PFjetCollection1_        = iConfig.getParameter<edm::InputTag>("PFjetCollection1");
  PFjetCollection2_        = iConfig.getParameter<edm::InputTag>("PFjetCollection2");
  PFjetCollection3_        = iConfig.getParameter<edm::InputTag>("PFjetCollection3");
  PFpuCorrJetCollection1_  = iConfig.getParameter<edm::InputTag>("PFpuCorrJetCollection1");
  PFpuCorrJetCollection2_  = iConfig.getParameter<edm::InputTag>("PFpuCorrJetCollection2");
  PFpuCorrJetCollection3_  = iConfig.getParameter<edm::InputTag>("PFpuCorrJetCollection3");
  JPTjetCollection1_       = iConfig.getParameter<edm::InputTag>("JPTjetCollection1");
  JPTjetCollection2_       = iConfig.getParameter<edm::InputTag>("JPTjetCollection2");

  // btag collections
  PFJetsBTags_              = iConfig.getUntrackedParameter<edm::ParameterSet>("PFJetsBTags");
  PFPUcorrJetsBTags_        = iConfig.getUntrackedParameter<edm::ParameterSet>("PFPUcorrJetsBTags");

  // MVA based jet id collection
  PFjetMvaIdCollection_     = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("PFjetMvaIdCollection",std::vector<edm::InputTag>());
  PFpujetMvaIdCollection_   = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("PFpujetMvaIdCollection",std::vector<edm::InputTag>());
  metCollection_           = iConfig.getParameter<edm::InputTag>("metCollection");
  // corrmetCollection_       = iConfig.getParameter<edm::InputTag>("corrmetCollection");
  TCmetCollection_         = iConfig.getParameter<edm::InputTag>("TCmetCollection");
  PFmetCollection_         = iConfig.getParameter<edm::InputTag>("PFmetCollection");
  PFChMetCollection_       = iConfig.getParameter<edm::InputTag>("PFChMetCollection");                       
  leptonLinkedPFCandidates_ = iConfig.getParameter<edm::InputTag>("leptonLinkedPFCandidates");
  genMetCollection_        = iConfig.getParameter<edm::InputTag>("genMetCollection");
  chargedMetCollection_    = iConfig.getParameter<edm::InputTag>("chargedMetCollection");
  mcTruthCollection_       = iConfig.getParameter<edm::InputTag>("mcTruthCollection");
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

  // dump PFCandidates
  dumpPFCandidates_  = iConfig.getUntrackedParameter<bool>("dumpPFCandidates");

  // PFTau Discriminators
  tauDiscrByLeadingTrackFindingTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByLeadingTrackFindingTag");
  tauDiscrByLeadingTrackPtCutTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByLeadingTrackPtCutTag");
  tauDiscrByLeadingPionPtCutTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByLeadingPionPtCutTag");
  tauDiscrByIsolationTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByIsolationTag");
  tauDiscrByIsolationUsingLeadingPionTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByIsolationUsingLeadingPionTag");
  tauDiscrByTrackIsolationTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTrackIsolationTag");
  tauDiscrByTrackIsolationUsingLeadingPionTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTrackIsolationUsingLeadingPionTag");
  tauDiscrByECALIsolationTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByECALIsolationTag");
  tauDiscrByECALIsolationUsingLeadingPionTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByECALIsolationUsingLeadingPionTag");
  tauDiscrAgainstMuonTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrAgainstMuonTag");
  tauDiscrAgainstElectronTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrAgainstElectronTag");
  tauDiscrByTaNCTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCTag");
  tauDiscrByTaNCfrHalfPercentTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCfrHalfPercentTag");
  tauDiscrByTaNCfrOnePercentTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCfrOnePercentTag");
  tauDiscrByTaNCfrQuarterPercentTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCfrQuarterPercentTag");
  tauDiscrByTaNCfrTenthPercentTag_ = iConfig.getParameter<edm::InputTag>("tauDiscrByTaNCfrTenthPercentTag");
  // HPS PFTau Discriminators
  hpsTauDiscrByLooseElectronRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByLooseElectronRejectionTag");
  hpsTauDiscrByMediumElectronRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByMediumElectronRejectionTag");
  hpsTauDiscrByTightElectronRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByTightElectronRejectionTag");
  hpsTauDiscrByLooseMuonRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByLooseMuonRejectionTag");
  hpsTauDiscrByTightMuonRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByTightMuonRejectionTag");
  hpsTauDiscrByDecayModeFindingTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByDecayModeFindingTag");
  hpsTauDiscrByVLooseIsolationTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByVLooseIsolationTag");
  hpsTauDiscrByLooseIsolationTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByLooseIsolationTag");
  hpsTauDiscrByMediumIsolationTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByMediumIsolationTag");
  hpsTauDiscrByTightIsolationTag_ = iConfig.getParameter<edm::InputTag>("hpsTauDiscrByTightIsolationTag");
  // HPS Tanc Tau Discriminators
  hpsTancTausDiscrByLeadingTrackFindingTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByLeadingTrackFindingTag");
  hpsTancTausDiscrByLeadingTrackPtCutTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByLeadingTrackPtCutTag");
  hpsTancTausDiscrByLeadingPionPtCutTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByLeadingPionPtCutTag");
  hpsTancTausDiscrByTancTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTancTag");
  hpsTancTausDiscrByTancRawTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTancRawTag");
  hpsTancTausDiscrByTancVLooseTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTancVLooseTag");
  hpsTancTausDiscrByTancLooseTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTancLooseTag");
  hpsTancTausDiscrByTancMediumTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTancMediumTag");
  hpsTancTausDiscrByTancTightTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTancTightTag");
  hpsTancTausDiscrByLooseElectronRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByLooseElectronRejectionTag");
  hpsTancTausDiscrByMediumElectronRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByMediumElectronRejectionTag");
  hpsTancTausDiscrByTightElectronRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTightElectronRejectionTag");
  hpsTancTausDiscrByLooseMuonRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByLooseMuonRejectionTag");
  hpsTancTausDiscrByTightMuonRejectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTightMuonRejectionTag");
  hpsTancTausDiscrByDecayModeSelectionTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByDecayModeSelectionTag");
  hpsTancTausDiscrByVLooseIsolationTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByVLooseIsolationTag");
  hpsTancTausDiscrByLooseIsolationTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByLooseIsolationTag");
  hpsTancTausDiscrByMediumIsolationTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByMediumIsolationTag");
  hpsTancTausDiscrByTightIsolationTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByTightIsolationTag");
  hpsTancTausDiscrByFlightPathTag_ = iConfig.getParameter<edm::InputTag>("hpsTancTausDiscrByFlightPathTag");

  // Hcal collections
  hcalNoiseSummaryLabel_ = iConfig.getParameter<edm::InputTag>("hcalNoiseSummary");

  energyCorrectionF = EcalClusterFunctionFactory::get()->create("EcalClusterEnergyCorrection", iConfig);

  // PDF sets
  pdfSet1_ = iConfig.getParameter<edm::InputTag>("pdfSet1");
  pdfSet2_ = iConfig.getParameter<edm::InputTag>("pdfSet2");
  pdfSet3_ = iConfig.getParameter<edm::InputTag>("pdfSet3");

  namePdf1_ = iConfig.getUntrackedParameter<std::string>("namepdf1", "pdfSet1");
  namePdf2_ = iConfig.getUntrackedParameter<std::string>("namepdf2", "pdfSet2");
  namePdf3_ = iConfig.getUntrackedParameter<std::string>("namepdf3", "pdfSet3");

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
  treeFill.saveLHEComments(dumpLHE_);

  if(dumpMCTruth_) {

    bool firstEvent = (jevt_>1) ? false : true;
    treeFill.writeCollectionToTree( mcTruthCollection_, LHEComments_, iEvent, 100, firstEvent );
  }

  // fill the PDF weights
  if(dumpPdfWeight_) {
    CmsPdfWeightFiller pdfWeightFiller( tree_);
    pdfWeightFiller.writePdfWeightToTree(pdfSet1_, iEvent, iSetup, "", namePdf1_.c_str(), false);
    pdfWeightFiller.writePdfWeightToTree(pdfSet2_, iEvent, iSetup, "", namePdf2_.c_str(), false);
    pdfWeightFiller.writePdfWeightToTree(pdfSet3_, iEvent, iSetup, "", namePdf3_.c_str(), false);
  }

  jevtInRun_++;
  // fill the trigger paths info
  if(dumpTriggerResults_) {

    CmsTriggerTreeFiller triggerTreeFill (tree_, hltParms_) ;
    std::string prefix ("") ;
    std::string suffix ("Trg") ;
    triggerTreeFill.writeTriggerToTree (iEvent, prefix, suffix) ;

    /// fill the trigger mask in the Event tree
    bool firstEvent = true;
    if(jevt_>1) firstEvent=false; 
    CmsConditionsFiller conditionsFiller( tree_, hltParms_, trgNames_ );
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
    treeFill.saveEleID(true);
    // for full vertex fit conversion veto
    treeFill.setConversionsProdcer(conversions_);

    treeFill.writeCollectionToTree(electronCollection_, iEvent, iSetup, prefix, suffix, false);

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
    // for full vertex fit conversion veto
    treeFill.setConversionsProdcer(conversions_);
    treeFill.writeCollectionToTree(photonCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  // fill SC block
  if (dumpSCs_) {

      CmsSuperClusterFiller treeFill(tree_, 1000);
      std::string prefix("");
      std::string suffix("SC");
      treeFill.setEcalBarrelRecHits(ecalBarrelRecHits_);
      treeFill.setEcalEndcapRecHits(ecalEndcapRecHits_);
      treeFill.setESRecHits(esRecHits_);
      treeFill.setCalotowers(calotowersForIsolationProducer_);
      treeFill.setPositionCalc(posCalculator_);
      treeFill.setEnergyCorrectionFunction(energyCorrectionF);
      if (usePhotonFix_)
	{
	  treeFill.setPhotonFixElectron(phFixE_);
	  treeFill.setPhotonFixPhoton(phFixP_);
	}
      if (useEnergyRegression_)
	{
	  treeFill.setEGEnergyCorrectorElectron(eCorrector_);
	  treeFill.setEGEnergyCorrectorPhoton(pCorrector_);
	}
	
      treeFill.writeCollectionToTree(ecalSCCollection_, iEvent, iSetup, prefix, suffix, false);
      

      CmsSuperClusterFiller treeFillPF(tree_, 1000);
      suffix = "PFSC";
      treeFillPF.setEcalBarrelRecHits(ecalBarrelRecHits_);
      treeFillPF.setEcalEndcapRecHits(ecalEndcapRecHits_);
      treeFillPF.setESRecHits(esRecHits_);
      treeFillPF.setCalotowers(calotowersForIsolationProducer_);
      treeFillPF.setPositionCalc(posCalculator_);
      treeFillPF.setEnergyCorrectionFunction(energyCorrectionF);
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
    treeFill.setGeneralTracks(trackCollection_);
    treeFill.saveCand(saveCand_);
    treeFill.saveFatTrk(saveFatTrk_);
    treeFill.writeCollectionToTree(muonCollection_, iEvent, iSetup, prefix, suffix, false);
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
                                     tauDiscrByLeadingTrackFindingTag_,
                                     tauDiscrByLeadingTrackPtCutTag_,
                                     tauDiscrByLeadingPionPtCutTag_,
                                     tauDiscrByIsolationTag_,
                                     tauDiscrByIsolationUsingLeadingPionTag_,
                                     tauDiscrByTrackIsolationTag_,
                                     tauDiscrByTrackIsolationUsingLeadingPionTag_,
                                     tauDiscrByECALIsolationTag_,
                                     tauDiscrByECALIsolationUsingLeadingPionTag_, 
                                     tauDiscrAgainstMuonTag_,
                                     tauDiscrAgainstElectronTag_,
                                     tauDiscrByTaNCTag_,
                                     tauDiscrByTaNCfrHalfPercentTag_,
                                     tauDiscrByTaNCfrOnePercentTag_,
                                     tauDiscrByTaNCfrQuarterPercentTag_,
                                     tauDiscrByTaNCfrTenthPercentTag_,
                                     false);

    }

  if(dumphpsPFTaus_)
    {
      CmsPFTauFiller treeFill(tree_, true);
      std::string prefix("");
      std::string suffix("PFTau");
      treeFill.savePFTauBasic(savePFTauBasic_);
      treeFill.saveLeadPFCand(saveLeadPFCand_);
      treeFill.writeCollectionToTree(hpspfTauCollection_, iEvent, iSetup, prefix, suffix,
                                     hpsTauDiscrByLooseElectronRejectionTag_,
                                     hpsTauDiscrByMediumElectronRejectionTag_,
                                     hpsTauDiscrByTightElectronRejectionTag_,
                                     hpsTauDiscrByLooseMuonRejectionTag_,
                                     hpsTauDiscrByTightMuonRejectionTag_,
                                     hpsTauDiscrByDecayModeFindingTag_,
                                     hpsTauDiscrByVLooseIsolationTag_,
                                     hpsTauDiscrByLooseIsolationTag_,
                                     hpsTauDiscrByMediumIsolationTag_,
                                     hpsTauDiscrByTightIsolationTag_,
                                     false);
    }

  if(dumphpsTancTaus_)
    {
      CmsPFTauFiller treeFill(tree_, true);
      std::string prefix("");
      std::string suffix("PFTau");
      treeFill.savePFTauBasic(savePFTauBasic_);
      treeFill.saveLeadPFCand(saveLeadPFCand_);
      treeFill.writeCollectionToTree(hpsTancTausCollection_, iEvent, iSetup, prefix, suffix,
				     hpsTancTausDiscrByLeadingTrackFindingTag_,
				     hpsTancTausDiscrByLeadingTrackPtCutTag_,
				     hpsTancTausDiscrByLeadingPionPtCutTag_,
				     hpsTancTausDiscrByTancTag_,
				     hpsTancTausDiscrByTancRawTag_,
				     hpsTancTausDiscrByTancVLooseTag_,
				     hpsTancTausDiscrByTancLooseTag_,
				     hpsTancTausDiscrByTancMediumTag_,
				     hpsTancTausDiscrByTancTightTag_,
				     hpsTancTausDiscrByLooseElectronRejectionTag_,
				     hpsTancTausDiscrByMediumElectronRejectionTag_,
				     hpsTancTausDiscrByTightElectronRejectionTag_,
				     hpsTancTausDiscrByLooseMuonRejectionTag_,
				     hpsTancTausDiscrByTightMuonRejectionTag_,
				     hpsTancTausDiscrByDecayModeSelectionTag_,
				     hpsTancTausDiscrByVLooseIsolationTag_,
				     hpsTancTausDiscrByLooseIsolationTag_,
				     hpsTancTausDiscrByMediumIsolationTag_,
				     hpsTancTausDiscrByTightIsolationTag_,
				     hpsTancTausDiscrByFlightPathTag_,
				     false);
    }

  // PF candidates
  if(dumpPFCandidates_) {
    CmsPFCandidateFiller treeFill(tree_);
    std::string prefix("");
    std::string suffix("PFCand");
    treeFill.saveCand(saveCand_);
    treeFill.setGeneralTracks(trackCollection_);
    treeFill.writeCollectionToTree(PFCandidateCollection_, iEvent, iSetup, prefix, suffix, false);
  }

  // PF candidates. Only the ones in a cone 0.1 from leptons to correct ChMET. In Candidate format to save space
  if(dumpParticleFlowObjects_) {
    CmsCandidateFiller treeFill(tree_);
    std::string prefix("");
    std::string suffix("ReducedPFCand");
    treeFill.writeCollectionToTree(leptonLinkedPFCandidates_, iEvent, iSetup, prefix, suffix, false);
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
    treeRecoFill1.isData(!dumpMCTruth_);
    treeRecoFill1.writeCollectionToTree(metCollection_, iEvent, iSetup, prefix, suffix, false);

    // Corrected CALO MET
    // CmsCandidateFiller treeRecoFill1bis(tree_, true);
    // suffix = "CorrMet";
    // treeRecoFill1bis.saveCand(saveCand_);
    // treeRecoFill1bis.writeCollectionToTree(corrmetCollection_, iEvent, iSetup, prefix, suffix, false);

    // Track-Corrected MET
    CmsMetFiller treeRecoFill2(tree_, true);
    suffix = "TCMet";
    treeRecoFill2.isData(false); // the met flags are per event, dumped for caloMET
    treeRecoFill2.saveCand(saveCand_);
    treeRecoFill2.writeCollectionToTree(TCmetCollection_, iEvent, iSetup, prefix, suffix, false);

    // particle flow met
    // [0] = uncorrected PF met
    // [0] = Type-0 corrected PFMET
    // [1] = Type-0 and Type-I corrected PFMET
    // [2] = Type-0, Type-I, and Type-II corrected PFMET
    if ( dumpParticleFlowObjects_ ) {
      CmsMetFiller pfMetFiller(tree_, 4, 4, true);
      suffix = "PFMet";
      pfMetFiller.saveCand(saveCand_);
      pfMetFiller.isData(false); // the met flags are per event, dumped for caloMET
      pfMetFiller.writeCollectionToTree(PFmetCollection_, iEvent, iSetup, prefix, suffix, false);

      // charged PF MET HWW version
      CmsMetFiller treeRecoFill1(tree_, true);
      std::string prefix("");
      std::string suffix("PFChMet");
      treeRecoFill1.isData(false); // the met flags are per event, dumped for caloMET
      treeRecoFill1.saveCand(saveCand_);
      treeRecoFill1.writeCollectionToTree(PFChMetCollection_, iEvent, iSetup, prefix, suffix, false);
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

    // Calo jets: not used for the moment
    CmsJetFiller caloJetFiller(tree_, true);
    std::string prefix("");
    std::string suffix("AK5Jet");
    caloJetFiller.saveCand(saveCand_);
    caloJetFiller.saveJetExtras(true);
    caloJetFiller.saveJetBTag(saveJetBTag_);
    caloJetFiller.writeCollectionToTree(jetCollection1_, iEvent, iSetup, prefix, suffix, false, jetCollection2_, jetCollection3_);

    // particle flow jets
    if ( dumpParticleFlowObjects_ ) {  
      CmsPFJetFiller pfJetFiller(tree_, true);
      suffix = "AK5PFNoPUJet";
      pfJetFiller.saveCand(saveCand_);
      pfJetFiller.saveJetBTag(saveJetBTag_);
      pfJetFiller.setBTags(PFJetsBTags_);
      pfJetFiller.setMvaId(PFjetMvaIdCollection_);
      pfJetFiller.writeCollectionToTree(PFjetCollection1_, iEvent, iSetup, prefix, suffix, false, PFjetCollection2_, PFjetCollection3_);
    }

    // particle flow jets with correction for pileup
    if ( dumpParticleFlowObjects_ && dumpPUcorrPFJet_ ) {
      CmsPFJetFiller pfPUcorrJetFiller(tree_, true);
      suffix = "AK5PFPUcorrJet";
      pfPUcorrJetFiller.saveCand(saveCand_);
      pfPUcorrJetFiller.saveJetBTag(saveJetBTag_);
      pfPUcorrJetFiller.setBTags(PFPUcorrJetsBTags_);
      pfPUcorrJetFiller.setMvaId(PFpujetMvaIdCollection_);
      pfPUcorrJetFiller.writeCollectionToTree(PFpuCorrJetCollection1_, iEvent, iSetup, prefix, suffix, false, PFpuCorrJetCollection2_, PFpuCorrJetCollection3_);
    }

    // Jet Plus Tracks jets: not used for the moment
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
  else hltObjectFiller_ = 0;

  // this pointer MUST survive until tree is closed
  trgNames_ = new vector<std::string>;
  LHEComments_ = new vector<std::string>;

}

void HWWTreeDumper::beginRun( const Run & iRun, const EventSetup & iSetup )
{
  jevtInRun_ = 0;
  if (usePhotonFix_)
    {
      PhotonFix::initialiseGeometry(iSetup);
      phFixE_=new PhotonFix(0.,0.);
      phFixE_->initialiseParameters(phFixElePar_);
      phFixP_=new PhotonFix(0.,0.);
      phFixP_->initialiseParameters(phFixPhoPar_);
    }

  if (useEnergyRegression_)
    {
      eCorrector_=new EGEnergyCorrector();
      pCorrector_=new EGEnergyCorrector();
      if (!eCorrector_->IsInitialized()) eCorrector_->Initialize(iSetup,energyRegressionElectronFile_);
      if (!pCorrector_->IsInitialized()) pCorrector_->Initialize(iSetup,energyRegressionPhotonFile_);
    }

  energyCorrectionF->init(iSetup);
  posCalculator_ = new PositionCalc(posCalcParameters_);
}



// ------------ method called once each job just after ending the event loop  ------------
void  HWWTreeDumper::endJob() {

  fileOut_->cd();

  TTree* treeEventsOut = tree_->getTree();
  treeEventsOut->Write();

  fileOut_->Close();

}


