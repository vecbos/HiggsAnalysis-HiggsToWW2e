//---------------------------------------------------------------------------
//
// Description:
//       Package:   HiggsAnalysis/HiggsToWW2e
//       Class:     CmsPFTauFiller
//
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFTauFiller.h"

#include <TTree.h>
#include <string>

using namespace edm;
using namespace reco;
using namespace std;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsPFTauFiller::CmsPFTauFiller(CmsTree *cmsTree, 
			       int maxTracks, int maxMCTracks,
			       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFTauFillerData)
{
  cmstree=cmsTree;

  savePFTauBasic_  = true;
  saveLeadPFCand_  = true;
  savePFTauDiscriminators_ = true;
  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();  
}

//--------------
// Destructor --
//--------------

CmsPFTauFiller::~CmsPFTauFiller() 
{

  // The following have to be all re-defined according to 
  // DataFormats/TauReco/src/PFTau.cc
  // delete privateData_->tauId;

  delete privateData_->pfJetIndex;
  delete privateData_->isNonNull;
  delete privateData_->charge;
  delete privateData_->energy;
  delete privateData_->et;
  delete privateData_->momentum;
  delete privateData_->theta;
  delete privateData_->eta;
  delete privateData_->phi;
  delete privateData_->x;
  delete privateData_->y;
  delete privateData_->z;
  delete privateData_->vertexX;
  delete privateData_->vertexY;
  delete privateData_->vertexZ;

  // Leading Track: for Electrons
  delete privateData_->isolationPFChargedHadrCandsPtSum;
  delete privateData_->isolationPFGammaCandsEtSum;
  delete privateData_->emFraction;
  delete privateData_->hcalTotOverPLead;
  delete privateData_->hcal3x3OverPLead;
  delete privateData_->ecalStripSumEOverPLead;
  delete privateData_->bremsRecoveryEOverPLead;

  // from the Discriminator for shrinking cone
  delete privateData_->theTauDiscrByLeadingTrackFinding;
  delete privateData_->theTauDiscrByLeadingTrackPtCut;
  delete privateData_->theTauDiscrByLeadingPionPtCut;
  delete privateData_->theTauDiscrByIsolation;
  delete privateData_->theTauDiscrByIsolationUsingLeadingPion;
  delete privateData_->theTauDiscrByTrackIsolation;
  delete privateData_->theTauDiscrByTrackIsolationUsingLeadingPion;
  delete privateData_->theTauDiscrByECALIsolation;
  delete privateData_->theTauDiscrByECALIsolationUsingLeadingPion;
  delete privateData_->theTauDiscrAgainstMuon;
  delete privateData_->theTauDiscrAgainstElectron;
  delete privateData_->theTauDiscrByTaNC;
  delete privateData_->theTauDiscrByTaNCfrHalfPercent;
  delete privateData_->theTauDiscrByTaNCfrOnePercent;
  delete privateData_->theTauDiscrByTaNCfrQuarterPercent;
  delete privateData_->theTauDiscrByTaNCfrTenthPercent;

  // from the Discriminator for HPS
  delete privateData_->thehpsTauDiscrByLooseElectronRejection;
  delete privateData_->thehpsTauDiscrByMediumElectronRejection;
  delete privateData_->thehpsTauDiscrByTightElectronRejection;
  delete privateData_->thehpsTauDiscrByLooseMuonRejection;
  delete privateData_->thehpsTauDiscrByTightMuonRejection;
  delete privateData_->thehpsTauDiscrByDecayModeFinding;
  delete privateData_->thehpsTauDiscrByVLooseIsolation;
  delete privateData_->thehpsTauDiscrByLooseIsolation;
  delete privateData_->thehpsTauDiscrByMediumIsolation;
  delete privateData_->thehpsTauDiscrByTightIsolation;

  // from the Discriminator for HPS TaNC
  delete privateData_->thehpsTancTausDiscrByLeadingTrackFinding;
  delete privateData_->thehpsTancTausDiscrByLeadingTrackPtCut;
  delete privateData_->thehpsTancTausDiscrByLeadingPionPtCut;
  delete privateData_->thehpsTancTausDiscrByTanc;
  delete privateData_->thehpsTancTausDiscrByTancRaw;
  delete privateData_->thehpsTancTausDiscrByTancVLoose;
  delete privateData_->thehpsTancTausDiscrByTancLoose;
  delete privateData_->thehpsTancTausDiscrByTancMedium;
  delete privateData_->thehpsTancTausDiscrByTancTight;
  delete privateData_->thehpsTancTausDiscrByLooseElectronRejection;
  delete privateData_->thehpsTancTausDiscrByMediumElectronRejection;
  delete privateData_->thehpsTancTausDiscrByTightElectronRejection;
  delete privateData_->thehpsTancTausDiscrByLooseMuonRejection;
  delete privateData_->thehpsTancTausDiscrByTightMuonRejection;
  delete privateData_->thehpsTancTausDiscrByDecayModeSelection;
  delete privateData_->thehpsTancTausDiscrByVLooseIsolation;
  delete privateData_->thehpsTancTausDiscrByLooseIsolation;
  delete privateData_->thehpsTancTausDiscrByMediumIsolation;
  delete privateData_->thehpsTancTausDiscrByTightIsolation;

  delete privateData_->ncand;
  delete privateData_;
}


//-------------
// Methods   --
//-------------

void CmsPFTauFiller::savePFTauBasic(bool what) { savePFTauBasic_=what; }

void CmsPFTauFiller::saveLeadPFCand(bool what) { saveLeadPFCand_=what; }

void CmsPFTauFiller::writeCollectionToTree(edm::InputTag collectionTag,
					   const edm::Event& iEvent, const edm::EventSetup& iSetup,
					   const std::string &columnPrefix, const std::string &columnSuffix,
					   edm::InputTag tauDiscrByLeadingTrackFindingTag,
					   edm::InputTag tauDiscrByLeadingTrackPtCutTag,
					   edm::InputTag tauDiscrByLeadingPionPtCutTag,
					   edm::InputTag tauDiscrByIsolationTag,
					   edm::InputTag tauDiscrByIsolationUsingLeadingPionTag,
					   edm::InputTag tauDiscrByTrackIsolationTag,
					   edm::InputTag tauDiscrByTrackIsolationUsingLeadingPionTag,
					   edm::InputTag tauDiscrByECALIsolationTag,
					   edm::InputTag tauDiscrByECALIsolationUsingLeadingPionTag, 
					   edm::InputTag tauDiscrAgainstMuonTag,
					   edm::InputTag tauDiscrAgainstElectronTag,
					   edm::InputTag tauDiscrByTaNCTag,
					   edm::InputTag tauDiscrByTaNCfrHalfPercentTag,
					   edm::InputTag tauDiscrByTaNCfrOnePercentTag,
					   edm::InputTag tauDiscrByTaNCfrQuarterPercentTag,
					   edm::InputTag tauDiscrByTaNCfrTenthPercentTag,
					   bool dumpData) {

  edm::Handle<reco::PFTauCollection> collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau candidate collection: " << collectionTag; }
  const reco::PFTauCollection *collection = collectionHandle.product();

  privateData_->clear();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && (int)collection->size() > maxTracks_) {
      edm::LogInfo("CmsPFTauFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_) {
      edm::LogInfo("CmsPFTauFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ 
				     << ". Collection will be truncated ";
    }

    *(privateData_->ncand) = collection->size();

    // edm::Handle<reco::PFTauCollection> tauJets;

    for ( int iTauJet = 0; iTauJet < (int)collection->size(); ++iTauJet) { //original

      // Maurizio's method: he just took the tauJets, which is the same as collection
      // const reco::PFTau& cand = tauJets.at(iTauJet); // <- maybe this has to be a pointer
      const reco::PFTau& cand = collection->at(iTauJet);
      
      // fill basic kinematics (pT, eta, phi, x, y, z)
      // if(savePFTauBasic_) writePFTauBasicInfo(&(*cand),iEvent,iSetup);
      if(savePFTauBasic_) writePFTauBasicInfo(&(cand),iEvent,iSetup);
    
      // fill Discrminators
      if(savePFTauDiscriminators_) {
	edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadingTrackFindingHandle;
	try { iEvent.getByLabel(tauDiscrByLeadingTrackFindingTag, tauDiscrByLeadingTrackFindingHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByLeadingTrackFindingTag; }
	const reco::PFTauDiscriminator *tauDiscrByLeadingTrackFinding = tauDiscrByLeadingTrackFindingHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadingTrackPtCutHandle;
	try { iEvent.getByLabel(tauDiscrByLeadingTrackPtCutTag, tauDiscrByLeadingTrackPtCutHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByLeadingTrackPtCutTag; }
	const reco::PFTauDiscriminator *tauDiscrByLeadingTrackPtCut = tauDiscrByLeadingTrackPtCutHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadingPionPtCutHandle;
	try { iEvent.getByLabel(tauDiscrByLeadingPionPtCutTag, tauDiscrByLeadingPionPtCutHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByLeadingPionPtCutTag; }
	const reco::PFTauDiscriminator *tauDiscrByLeadingPionPtCut = tauDiscrByLeadingPionPtCutHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByIsolationHandle;
	try { iEvent.getByLabel(tauDiscrByIsolationTag, tauDiscrByIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByIsolationTag; }
	const reco::PFTauDiscriminator *tauDiscrByIsolation = tauDiscrByIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByIsolationUsingLeadingPionHandle;
	try { iEvent.getByLabel(tauDiscrByIsolationUsingLeadingPionTag, tauDiscrByIsolationUsingLeadingPionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByIsolationUsingLeadingPionTag; }
	const reco::PFTauDiscriminator *tauDiscrByIsolationUsingLeadingPion = tauDiscrByIsolationUsingLeadingPionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTrackIsolationHandle;
	try { iEvent.getByLabel(tauDiscrByTrackIsolationTag, tauDiscrByTrackIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTrackIsolationTag; }
	const reco::PFTauDiscriminator *tauDiscrByTrackIsolation = tauDiscrByTrackIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTrackIsolationUsingLeadingPionHandle;
	try { iEvent.getByLabel(tauDiscrByTrackIsolationUsingLeadingPionTag, tauDiscrByTrackIsolationUsingLeadingPionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTrackIsolationUsingLeadingPionTag; }
	const reco::PFTauDiscriminator *tauDiscrByTrackIsolationUsingLeadingPion = tauDiscrByTrackIsolationUsingLeadingPionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByECALIsolationHandle;
	try { iEvent.getByLabel(tauDiscrByECALIsolationTag, tauDiscrByECALIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByECALIsolationTag; }
	const reco::PFTauDiscriminator *tauDiscrByECALIsolation = tauDiscrByECALIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByECALIsolationUsingLeadingPionHandle;
	try { iEvent.getByLabel(tauDiscrByECALIsolationUsingLeadingPionTag, tauDiscrByECALIsolationUsingLeadingPionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByECALIsolationUsingLeadingPionTag; }
	const reco::PFTauDiscriminator *tauDiscrByECALIsolationUsingLeadingPion = tauDiscrByECALIsolationUsingLeadingPionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrAgainstMuonHandle;
	try { iEvent.getByLabel(tauDiscrAgainstMuonTag, tauDiscrAgainstMuonHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrAgainstMuonTag; }
	const reco::PFTauDiscriminator *tauDiscrAgainstMuon = tauDiscrAgainstMuonHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrAgainstElectronHandle;
	try { iEvent.getByLabel(tauDiscrAgainstElectronTag, tauDiscrAgainstElectronHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrAgainstElectronTag; }
	const reco::PFTauDiscriminator *tauDiscrAgainstElectron = tauDiscrAgainstElectronHandle.product();

 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCHandle;
 	try { iEvent.getByLabel(tauDiscrByTaNCTag, tauDiscrByTaNCHandle); }
 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCTag; }
 	const reco::PFTauDiscriminator *tauDiscrByTaNC = tauDiscrByTaNCHandle.product();

 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCfrHalfPercentHandle;
 	try { iEvent.getByLabel(tauDiscrByTaNCfrHalfPercentTag, tauDiscrByTaNCfrHalfPercentHandle); }
 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCfrHalfPercentTag; }
 	const reco::PFTauDiscriminator *tauDiscrByTaNCfrHalfPercent = tauDiscrByTaNCfrHalfPercentHandle.product();

 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCfrOnePercentHandle;
 	try { iEvent.getByLabel(tauDiscrByTaNCfrOnePercentTag, tauDiscrByTaNCfrOnePercentHandle); }
 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCfrOnePercentTag; }
 	const reco::PFTauDiscriminator *tauDiscrByTaNCfrOnePercent = tauDiscrByTaNCfrOnePercentHandle.product();

 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCfrQuarterPercentHandle;
 	try { iEvent.getByLabel(tauDiscrByTaNCfrQuarterPercentTag, tauDiscrByTaNCfrQuarterPercentHandle); }
 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCfrQuarterPercentTag; }
 	const reco::PFTauDiscriminator *tauDiscrByTaNCfrQuarterPercent = tauDiscrByTaNCfrQuarterPercentHandle.product();

 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCfrTenthPercentHandle;
 	try { iEvent.getByLabel(tauDiscrByTaNCfrTenthPercentTag, tauDiscrByTaNCfrTenthPercentHandle); }
 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCfrTenthPercentTag; }
 	const reco::PFTauDiscriminator *tauDiscrByTaNCfrTenthPercent = tauDiscrByTaNCfrTenthPercentHandle.product();

	writePFTauDiscInfo(collectionHandle, iTauJet,
			   tauDiscrByLeadingTrackFinding,
			   tauDiscrByLeadingTrackPtCut,
			   tauDiscrByLeadingPionPtCut,
			   tauDiscrByIsolation,
			   tauDiscrByIsolationUsingLeadingPion,
			   tauDiscrByTrackIsolation,
			   tauDiscrByTrackIsolationUsingLeadingPion,
			   tauDiscrByECALIsolation,
			   tauDiscrByECALIsolationUsingLeadingPion,
			   tauDiscrAgainstMuon,
			   tauDiscrAgainstElectron,
			   tauDiscrByTaNC,
			   tauDiscrByTaNCfrHalfPercent,
			   tauDiscrByTaNCfrOnePercent,
			   tauDiscrByTaNCfrQuarterPercent,
			   tauDiscrByTaNCfrTenthPercent
                           );
      }

      // fill Additional Tau Info
      if(saveLeadPFCand_) writeLeadPFCandInfo(&(cand),iEvent,iSetup); // this method needs to be implemented
      // if(saveLeadPFCand_) writeLeadPFCandInfo(collectionHandle,iTauJet); // this method needs to be implemented

    }
  } else {
    edm::LogWarning("CmsPFTauFiller") << "Warning! The collection seems to be not made by "
				      << "pflow candidate taus, tau-specific infos will be set to default.";       
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(savePFTauBasic_) treePFTauBasicInfo(columnPrefix,columnSuffix);
  if(savePFTauDiscriminators_) treePFTauDiscInfo(columnPrefix,columnSuffix);
  if(saveLeadPFCand_) treeLeadPFCandInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();
	
}

void CmsPFTauFiller::writeCollectionToTree(edm::InputTag collectionTag,
					   const edm::Event& iEvent, const edm::EventSetup& iSetup,
					   const std::string &columnPrefix, const std::string &columnSuffix,
					   edm::InputTag hpsTauDiscrByLooseElectronRejectionTag,
					   edm::InputTag hpsTauDiscrByMediumElectronRejectionTag,
					   edm::InputTag hpsTauDiscrByTightElectronRejectionTag,
					   edm::InputTag hpsTauDiscrByLooseMuonRejectionTag,
					   edm::InputTag hpsTauDiscrByTightMuonRejectionTag,
					   edm::InputTag hpsTauDiscrByDecayModeFindingTag,
					   edm::InputTag hpsTauDiscrByVLooseIsolationTag,
					   edm::InputTag hpsTauDiscrByLooseIsolationTag,
					   edm::InputTag hpsTauDiscrByMediumIsolationTag,
					   edm::InputTag hpsTauDiscrByTightIsolationTag,
					   bool dumpData) {

  edm::Handle<reco::PFTauCollection> collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau candidate collection: " << collectionTag; }
  const reco::PFTauCollection *collection = collectionHandle.product();

  privateData_->clear();

  if(collection) {

    if(hitLimitsMeansNoOutput_ && (int)collection->size() > maxTracks_) {
      edm::LogInfo("CmsPFTauFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_) {
      edm::LogInfo("CmsPFTauFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ 
				     << ". Collection will be truncated ";
    }

    *(privateData_->ncand) = collection->size();

    for ( int iTauJet = 0; iTauJet < (int)collection->size(); ++iTauJet) {

      const reco::PFTau& cand = collection->at(iTauJet);
      if(savePFTauBasic_) writePFTauBasicInfo(&(cand),iEvent,iSetup);

      // fill Discrminators
      if(savePFTauDiscriminators_) {
	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByLooseElectronRejectionHandle;
	try { iEvent.getByLabel(hpsTauDiscrByLooseElectronRejectionTag, hpsTauDiscrByLooseElectronRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByLooseElectronRejectionTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByLooseElectronRejection = hpsTauDiscrByLooseElectronRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByMediumElectronRejectionHandle;
	try { iEvent.getByLabel(hpsTauDiscrByMediumElectronRejectionTag, hpsTauDiscrByMediumElectronRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByMediumElectronRejectionTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByMediumElectronRejection = hpsTauDiscrByMediumElectronRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByTightElectronRejectionHandle;
	try { iEvent.getByLabel(hpsTauDiscrByTightElectronRejectionTag, hpsTauDiscrByTightElectronRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByTightElectronRejectionTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByTightElectronRejection = hpsTauDiscrByTightElectronRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByLooseMuonRejectionHandle;
	try { iEvent.getByLabel(hpsTauDiscrByLooseMuonRejectionTag, hpsTauDiscrByLooseMuonRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByLooseMuonRejectionTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByLooseMuonRejection = hpsTauDiscrByLooseMuonRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByTightMuonRejectionHandle;
	try { iEvent.getByLabel(hpsTauDiscrByTightMuonRejectionTag, hpsTauDiscrByTightMuonRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByTightMuonRejectionTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByTightMuonRejection = hpsTauDiscrByTightMuonRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByDecayModeFindingHandle;
	try { iEvent.getByLabel(hpsTauDiscrByDecayModeFindingTag, hpsTauDiscrByDecayModeFindingHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByDecayModeFindingTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByDecayModeFinding = hpsTauDiscrByDecayModeFindingHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByVLooseIsolationHandle;
	try { iEvent.getByLabel(hpsTauDiscrByVLooseIsolationTag, hpsTauDiscrByVLooseIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByVLooseIsolationTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByVLooseIsolation = hpsTauDiscrByVLooseIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByLooseIsolationHandle;
	try { iEvent.getByLabel(hpsTauDiscrByLooseIsolationTag, hpsTauDiscrByLooseIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByLooseIsolationTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByLooseIsolation = hpsTauDiscrByLooseIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByMediumIsolationHandle;
	try { iEvent.getByLabel(hpsTauDiscrByMediumIsolationTag, hpsTauDiscrByMediumIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByMediumIsolationTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByMediumIsolation = hpsTauDiscrByMediumIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByTightIsolationHandle;
	try { iEvent.getByLabel(hpsTauDiscrByTightIsolationTag, hpsTauDiscrByTightIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTauDiscrByTightIsolationTag; }
	const reco::PFTauDiscriminator *hpsTauDiscrByTightIsolation = hpsTauDiscrByTightIsolationHandle.product();

	writePFTauDiscInfo(collectionHandle, iTauJet,
			   hpsTauDiscrByLooseElectronRejection,
			   hpsTauDiscrByMediumElectronRejection,
			   hpsTauDiscrByTightElectronRejection,
			   hpsTauDiscrByLooseMuonRejection,
			   hpsTauDiscrByTightMuonRejection,
			   hpsTauDiscrByDecayModeFinding,
			   hpsTauDiscrByVLooseIsolation,
			   hpsTauDiscrByLooseIsolation,
			   hpsTauDiscrByMediumIsolation,
			   hpsTauDiscrByTightIsolation
			   );
      }

      // fill Additional Tau Info
      if(saveLeadPFCand_) writeLeadPFCandInfo(&(cand),iEvent,iSetup); // this method needs to be implemented
    }
  } else {
    edm::LogWarning("CmsPFTauFiller") << "Warning! The collection seems to be not made by "
				      << "pflow candidate taus, tau-specific infos will be set to default.";       
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(savePFTauBasic_) treePFTauBasicInfo(columnPrefix,columnSuffix);
  if(savePFTauDiscriminators_) treehpsPFTauDiscInfo(columnPrefix,columnSuffix);
  if(saveLeadPFCand_) treeLeadPFCandInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();
}



void CmsPFTauFiller::writeCollectionToTree(edm::InputTag collectionTag,
					   const edm::Event& iEvent, const edm::EventSetup& iSetup,
					   const std::string &columnPrefix, const std::string &columnSuffix,
					   edm::InputTag hpsTancTausDiscrByLeadingTrackFindingTag,
					   edm::InputTag hpsTancTausDiscrByLeadingTrackPtCutTag,
					   edm::InputTag hpsTancTausDiscrByLeadingPionPtCutTag,
					   edm::InputTag hpsTancTausDiscrByTancTag,
					   edm::InputTag hpsTancTausDiscrByTancRawTag,
					   edm::InputTag hpsTancTausDiscrByTancVLooseTag,
					   edm::InputTag hpsTancTausDiscrByTancLooseTag,
					   edm::InputTag hpsTancTausDiscrByTancMediumTag,
					   edm::InputTag hpsTancTausDiscrByTancTightTag,
					   edm::InputTag hpsTancTausDiscrByLooseElectronRejectionTag,
					   edm::InputTag hpsTancTausDiscrByMediumElectronRejectionTag,
					   edm::InputTag hpsTancTausDiscrByTightElectronRejectionTag,
					   edm::InputTag hpsTancTausDiscrByLooseMuonRejectionTag,
					   edm::InputTag hpsTancTausDiscrByTightMuonRejectionTag,
					   edm::InputTag hpsTancTausDiscrByDecayModeSelectionTag,
					   edm::InputTag hpsTancTausDiscrByVLooseIsolationTag,
					   edm::InputTag hpsTancTausDiscrByLooseIsolationTag,
					   edm::InputTag hpsTancTausDiscrByMediumIsolationTag,
					   edm::InputTag hpsTancTausDiscrByTightIsolationTag,
					   bool dumpData) {

  edm::Handle<reco::PFTauCollection> collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau candidate collection: " << collectionTag; }
  const reco::PFTauCollection *collection = collectionHandle.product();

  privateData_->clear();

  if(collection) {

    if(hitLimitsMeansNoOutput_ && (int)collection->size() > maxTracks_) {
      edm::LogInfo("CmsPFTauFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_) {
      edm::LogInfo("CmsPFTauFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ 
				     << ". Collection will be truncated ";
    }

    *(privateData_->ncand) = collection->size();

    for ( int iTauJet = 0; iTauJet < (int)collection->size(); ++iTauJet) {

      const reco::PFTau& cand = collection->at(iTauJet);
      if(savePFTauBasic_) writePFTauBasicInfo(&(cand),iEvent,iSetup);

      // fill Discrminators
      if(savePFTauDiscriminators_) {
	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByLeadingTrackFindingHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByLeadingTrackFindingTag, hpsTancTausDiscrByLeadingTrackFindingHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByLeadingTrackFindingTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingTrackFinding = hpsTancTausDiscrByLeadingTrackFindingHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByLeadingTrackPtCutHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByLeadingTrackPtCutTag, hpsTancTausDiscrByLeadingTrackPtCutHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByLeadingTrackPtCutTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingTrackPtCut = hpsTancTausDiscrByLeadingTrackPtCutHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByLeadingPionPtCutHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByLeadingPionPtCutTag, hpsTancTausDiscrByLeadingPionPtCutHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByLeadingPionPtCutTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingPionPtCut = hpsTancTausDiscrByLeadingPionPtCutHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTancHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTancTag, hpsTancTausDiscrByTancHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTancTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTanc = hpsTancTausDiscrByTancHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTancRawHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTancRawTag, hpsTancTausDiscrByTancRawHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTancRawTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTancRaw = hpsTancTausDiscrByTancRawHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTancVLooseHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTancVLooseTag, hpsTancTausDiscrByTancVLooseHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTancVLooseTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTancVLoose = hpsTancTausDiscrByTancVLooseHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTancLooseHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTancLooseTag, hpsTancTausDiscrByTancLooseHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTancLooseTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTancLoose = hpsTancTausDiscrByTancLooseHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTancMediumHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTancMediumTag, hpsTancTausDiscrByTancMediumHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTancMediumTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTancMedium = hpsTancTausDiscrByTancMediumHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTancTightHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTancTightTag, hpsTancTausDiscrByTancTightHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTancTightTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTancTight = hpsTancTausDiscrByTancTightHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByLooseElectronRejectionHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByLooseElectronRejectionTag, hpsTancTausDiscrByLooseElectronRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByLooseElectronRejectionTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseElectronRejection = hpsTancTausDiscrByLooseElectronRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByMediumElectronRejectionHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByMediumElectronRejectionTag, hpsTancTausDiscrByMediumElectronRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByMediumElectronRejectionTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByMediumElectronRejection = hpsTancTausDiscrByMediumElectronRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTightElectronRejectionHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTightElectronRejectionTag, hpsTancTausDiscrByTightElectronRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTightElectronRejectionTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTightElectronRejection = hpsTancTausDiscrByTightElectronRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByLooseMuonRejectionHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByLooseMuonRejectionTag, hpsTancTausDiscrByLooseMuonRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByLooseMuonRejectionTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseMuonRejection = hpsTancTausDiscrByLooseMuonRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTightMuonRejectionHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTightMuonRejectionTag, hpsTancTausDiscrByTightMuonRejectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTightMuonRejectionTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTightMuonRejection = hpsTancTausDiscrByTightMuonRejectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByDecayModeSelectionHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByDecayModeSelectionTag, hpsTancTausDiscrByDecayModeSelectionHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByDecayModeSelectionTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByDecayModeSelection = hpsTancTausDiscrByDecayModeSelectionHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByVLooseIsolationHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByVLooseIsolationTag, hpsTancTausDiscrByVLooseIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByVLooseIsolationTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByVLooseIsolation = hpsTancTausDiscrByVLooseIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByLooseIsolationHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByLooseIsolationTag, hpsTancTausDiscrByLooseIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByLooseIsolationTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseIsolation = hpsTancTausDiscrByLooseIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByMediumIsolationHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByMediumIsolationTag, hpsTancTausDiscrByMediumIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByMediumIsolationTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByMediumIsolation = hpsTancTausDiscrByMediumIsolationHandle.product();

	edm::Handle<reco::PFTauDiscriminator> hpsTancTausDiscrByTightIsolationHandle;
	try { iEvent.getByLabel(hpsTancTausDiscrByTightIsolationTag, hpsTancTausDiscrByTightIsolationHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << hpsTancTausDiscrByTightIsolationTag; }
	const reco::PFTauDiscriminator *hpsTancTausDiscrByTightIsolation = hpsTancTausDiscrByTightIsolationHandle.product();

	writePFTauDiscInfo(collectionHandle, iTauJet,
			   hpsTancTausDiscrByLeadingTrackFinding,    
			   hpsTancTausDiscrByLeadingTrackPtCut,	     
			   hpsTancTausDiscrByLeadingPionPtCut,	     
			   hpsTancTausDiscrByTanc,		     
			   hpsTancTausDiscrByTancRaw,		     
			   hpsTancTausDiscrByTancVLoose,	     
			   hpsTancTausDiscrByTancLoose,		     
			   hpsTancTausDiscrByTancMedium,	     
			   hpsTancTausDiscrByTancTight,		     
			   hpsTancTausDiscrByLooseElectronRejection, 
			   hpsTancTausDiscrByMediumElectronRejection,
			   hpsTancTausDiscrByTightElectronRejection,
			   hpsTancTausDiscrByLooseMuonRejection,
			   hpsTancTausDiscrByTightMuonRejection,
			   hpsTancTausDiscrByDecayModeSelection,
			   hpsTancTausDiscrByVLooseIsolation,
			   hpsTancTausDiscrByLooseIsolation, 
			   hpsTancTausDiscrByMediumIsolation,
			   hpsTancTausDiscrByTightIsolation
			   );
      }

      // fill Additional Tau Info
      if(saveLeadPFCand_) writeLeadPFCandInfo(&(cand),iEvent,iSetup); // this method needs to be implemented
    }
  } else {
    edm::LogWarning("CmsPFTauFiller") << "Warning! The collection seems to be not made by "
				      << "pflow candidate taus, tau-specific infos will be set to default.";       
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(savePFTauBasic_) treePFTauBasicInfo(columnPrefix,columnSuffix);
  if(savePFTauDiscriminators_) treehpsTancTausDiscInfo(columnPrefix,columnSuffix);
  if(saveLeadPFCand_) treeLeadPFCandInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();
  
  delete trkIndexName_;

}


// First the write methods **************************
void CmsPFTauFiller::writePFTauBasicInfo(const reco::PFTau *tau, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  privateData_->pfJetIndex->push_back((int)tau->jetRef().key());
  privateData_->isNonNull->push_back((int)tau->leadPFChargedHadrCand().isNonnull());
  privateData_->charge->push_back((int)tau->charge());
  privateData_->energy->push_back(tau->energy());
  privateData_->et->push_back(tau->et());
  privateData_->momentum->push_back(tau->p());
  privateData_->theta->push_back(tau->theta());
  privateData_->eta->push_back(tau->eta());
  privateData_->phi->push_back(tau->phi());
  privateData_->x->push_back(tau->momentum().x());
  privateData_->y->push_back(tau->momentum().y());
  privateData_->z->push_back(tau->momentum().z());
  privateData_->vertexX->push_back(tau->vx());
  privateData_->vertexY->push_back(tau->vy());
  privateData_->vertexZ->push_back(tau->vz());
}

void CmsPFTauFiller::writePFTauDiscInfo(edm::Handle<reco::PFTauCollection> tauCollection, int theTauJetIndex,
					const reco::PFTauDiscriminator *tauDiscrByLeadingTrackFinding,
					const reco::PFTauDiscriminator *tauDiscrByLeadingTrackPtCut,
					const reco::PFTauDiscriminator *tauDiscrByLeadingPionPtCut,
					const reco::PFTauDiscriminator *tauDiscrByIsolation,
					const reco::PFTauDiscriminator *tauDiscrByIsolationUsingLeadingPion,
					const reco::PFTauDiscriminator *tauDiscrByTrackIsolation,
					const reco::PFTauDiscriminator *tauDiscrByTrackIsolationUsingLeadingPion,
					const reco::PFTauDiscriminator *tauDiscrByECALIsolation,
					const reco::PFTauDiscriminator *tauDiscrByECALIsolationUsingLeadingPion,
					const reco::PFTauDiscriminator *tauDiscrAgainstMuon,
					const reco::PFTauDiscriminator *tauDiscrAgainstElectron,
					const reco::PFTauDiscriminator *tauDiscrByTaNC,
					const reco::PFTauDiscriminator *tauDiscrByTaNCfrHalfPercent,
					const reco::PFTauDiscriminator *tauDiscrByTaNCfrOnePercent,
					const reco::PFTauDiscriminator *tauDiscrByTaNCfrQuarterPercent,
					const reco::PFTauDiscriminator *tauDiscrByTaNCfrTenthPercent
                                        )
{
  if ( theTauJetIndex != -1 ) {
    reco::PFTauRef theTauJetRef(tauCollection, theTauJetIndex);
    privateData_->theTauDiscrByLeadingTrackFinding->push_back((*tauDiscrByLeadingTrackFinding)[theTauJetRef]);
    privateData_->theTauDiscrByLeadingTrackPtCut->push_back((*tauDiscrByLeadingTrackPtCut)[theTauJetRef]);
    privateData_->theTauDiscrByLeadingPionPtCut->push_back((*tauDiscrByLeadingPionPtCut)[theTauJetRef]);
    privateData_->theTauDiscrByIsolation->push_back((*tauDiscrByIsolation)[theTauJetRef]);
    privateData_->theTauDiscrByIsolationUsingLeadingPion->push_back((*tauDiscrByIsolationUsingLeadingPion)[theTauJetRef]);
    privateData_->theTauDiscrByTrackIsolation->push_back((*tauDiscrByTrackIsolation)[theTauJetRef]);
    privateData_->theTauDiscrByTrackIsolationUsingLeadingPion->push_back((*tauDiscrByTrackIsolationUsingLeadingPion)[theTauJetRef]);
    privateData_->theTauDiscrByECALIsolation->push_back((*tauDiscrByECALIsolation)[theTauJetRef]);
    privateData_->theTauDiscrByECALIsolationUsingLeadingPion->push_back((*tauDiscrByECALIsolationUsingLeadingPion)[theTauJetRef]);
    privateData_->theTauDiscrAgainstMuon->push_back((*tauDiscrAgainstMuon)[theTauJetRef]);
    privateData_->theTauDiscrAgainstElectron->push_back((*tauDiscrAgainstElectron)[theTauJetRef]);
    privateData_->theTauDiscrByTaNC->push_back((*tauDiscrByTaNC)[theTauJetRef]);
    privateData_->theTauDiscrByTaNCfrHalfPercent->push_back((*tauDiscrByTaNCfrHalfPercent)[theTauJetRef]);
    privateData_->theTauDiscrByTaNCfrOnePercent->push_back((*tauDiscrByTaNCfrOnePercent)[theTauJetRef]);
    privateData_->theTauDiscrByTaNCfrQuarterPercent->push_back((*tauDiscrByTaNCfrQuarterPercent)[theTauJetRef]);
    privateData_->theTauDiscrByTaNCfrTenthPercent->push_back((*tauDiscrByTaNCfrTenthPercent)[theTauJetRef]);
  } else {
    privateData_->theTauDiscrByLeadingTrackFinding->push_back(-1);
    privateData_->theTauDiscrByLeadingTrackPtCut->push_back(-1);
    privateData_->theTauDiscrByLeadingPionPtCut->push_back(-1);
    privateData_->theTauDiscrByIsolation->push_back(-1);
    privateData_->theTauDiscrByIsolationUsingLeadingPion->push_back(-1);
    privateData_->theTauDiscrByTrackIsolation->push_back(-1);
    privateData_->theTauDiscrByTrackIsolationUsingLeadingPion->push_back(-1);
    privateData_->theTauDiscrByECALIsolation->push_back(-1);
    privateData_->theTauDiscrByECALIsolationUsingLeadingPion->push_back(-1);
    privateData_->theTauDiscrAgainstMuon->push_back(-1);
    privateData_->theTauDiscrAgainstElectron->push_back(-1);
    privateData_->theTauDiscrByTaNC->push_back(-1);
    privateData_->theTauDiscrByTaNCfrHalfPercent->push_back(-1);
    privateData_->theTauDiscrByTaNCfrOnePercent->push_back(-1);
    privateData_->theTauDiscrByTaNCfrQuarterPercent->push_back(-1);
    privateData_->theTauDiscrByTaNCfrTenthPercent->push_back(-1);
  }
}

void CmsPFTauFiller::writePFTauDiscInfo(edm::Handle<reco::PFTauCollection> tauCollection, int theTauJetIndex,
					const reco::PFTauDiscriminator *hpsTauDiscrByLooseElectronRejection,
					const reco::PFTauDiscriminator *hpsTauDiscrByMediumElectronRejection,
					const reco::PFTauDiscriminator *hpsTauDiscrByTightElectronRejection,
					const reco::PFTauDiscriminator *hpsTauDiscrByLooseMuonRejection,
					const reco::PFTauDiscriminator *hpsTauDiscrByTightMuonRejection,
					const reco::PFTauDiscriminator *hpsTauDiscrByDecayModeFinding,
					const reco::PFTauDiscriminator *hpsTauDiscrByVLooseIsolation,
					const reco::PFTauDiscriminator *hpsTauDiscrByLooseIsolation,
					const reco::PFTauDiscriminator *hpsTauDiscrByMediumIsolation,
					const reco::PFTauDiscriminator *hpsTauDiscrByTightIsolation
					)
{
  if ( theTauJetIndex != -1 ) {
    reco::PFTauRef theTauJetRef(tauCollection, theTauJetIndex);
    privateData_->thehpsTauDiscrByLooseElectronRejection->push_back((*hpsTauDiscrByLooseElectronRejection)[theTauJetRef]);
    privateData_->thehpsTauDiscrByMediumElectronRejection->push_back((*hpsTauDiscrByMediumElectronRejection)[theTauJetRef]);
    privateData_->thehpsTauDiscrByTightElectronRejection->push_back((*hpsTauDiscrByTightElectronRejection)[theTauJetRef]);
    privateData_->thehpsTauDiscrByLooseMuonRejection->push_back((*hpsTauDiscrByLooseMuonRejection)[theTauJetRef]);
    privateData_->thehpsTauDiscrByTightMuonRejection->push_back((*hpsTauDiscrByTightMuonRejection)[theTauJetRef]);
    privateData_->thehpsTauDiscrByDecayModeFinding->push_back((*hpsTauDiscrByDecayModeFinding)[theTauJetRef]);
    privateData_->thehpsTauDiscrByVLooseIsolation->push_back((*hpsTauDiscrByVLooseIsolation)[theTauJetRef]);
    privateData_->thehpsTauDiscrByLooseIsolation->push_back((*hpsTauDiscrByLooseIsolation)[theTauJetRef]);
    privateData_->thehpsTauDiscrByMediumIsolation->push_back((*hpsTauDiscrByMediumIsolation)[theTauJetRef]);
    privateData_->thehpsTauDiscrByTightIsolation->push_back((*hpsTauDiscrByTightIsolation)[theTauJetRef]);
  } else {
    privateData_->thehpsTauDiscrByLooseElectronRejection->push_back(-1);
    privateData_->thehpsTauDiscrByMediumElectronRejection->push_back(-1);
    privateData_->thehpsTauDiscrByTightElectronRejection->push_back(-1);
    privateData_->thehpsTauDiscrByLooseMuonRejection->push_back(-1);
    privateData_->thehpsTauDiscrByTightMuonRejection->push_back(-1);
    privateData_->thehpsTauDiscrByDecayModeFinding->push_back(-1);
    privateData_->thehpsTauDiscrByVLooseIsolation->push_back(-1);
    privateData_->thehpsTauDiscrByLooseIsolation->push_back(-1);
    privateData_->thehpsTauDiscrByMediumIsolation->push_back(-1);
    privateData_->thehpsTauDiscrByTightIsolation->push_back(-1);
  }
}

void CmsPFTauFiller::writePFTauDiscInfo(edm::Handle<reco::PFTauCollection> tauCollection, int theTauJetIndex,
					const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingTrackFinding,    
					const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingTrackPtCut,	     
					const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingPionPtCut,	     
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTanc,		     
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTancRaw,		     
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTancVLoose,	     
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTancLoose,		     
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTancMedium,	     
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTancTight,		     
					const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseElectronRejection, 
					const reco::PFTauDiscriminator *hpsTancTausDiscrByMediumElectronRejection,
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTightElectronRejection,
					const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseMuonRejection,
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTightMuonRejection,
					const reco::PFTauDiscriminator *hpsTancTausDiscrByDecayModeSelection,
					const reco::PFTauDiscriminator *hpsTancTausDiscrByVLooseIsolation,
					const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseIsolation, 
					const reco::PFTauDiscriminator *hpsTancTausDiscrByMediumIsolation,
					const reco::PFTauDiscriminator *hpsTancTausDiscrByTightIsolation
					)
{
  if ( theTauJetIndex != -1 ) {
    reco::PFTauRef theTauJetRef(tauCollection, theTauJetIndex);
    privateData_->thehpsTancTausDiscrByLeadingTrackFinding->push_back((*hpsTancTausDiscrByLeadingTrackFinding)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByLeadingTrackPtCut->push_back((*hpsTancTausDiscrByLeadingTrackPtCut)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByLeadingPionPtCut->push_back((*hpsTancTausDiscrByLeadingPionPtCut)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTanc->push_back((*hpsTancTausDiscrByTanc)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTancRaw->push_back((*hpsTancTausDiscrByTancRaw)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTancVLoose->push_back((*hpsTancTausDiscrByTancVLoose)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTancLoose->push_back((*hpsTancTausDiscrByTancLoose)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTancMedium->push_back((*hpsTancTausDiscrByTancMedium)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTancTight->push_back((*hpsTancTausDiscrByTancTight)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByLooseElectronRejection->push_back((*hpsTancTausDiscrByLooseElectronRejection)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByMediumElectronRejection->push_back((*hpsTancTausDiscrByMediumElectronRejection)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTightElectronRejection->push_back((*hpsTancTausDiscrByTightElectronRejection)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByLooseMuonRejection->push_back((*hpsTancTausDiscrByLooseMuonRejection)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTightMuonRejection->push_back((*hpsTancTausDiscrByTightMuonRejection)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByDecayModeSelection->push_back((*hpsTancTausDiscrByDecayModeSelection)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByVLooseIsolation->push_back((*hpsTancTausDiscrByVLooseIsolation)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByLooseIsolation->push_back((*hpsTancTausDiscrByLooseIsolation)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByMediumIsolation->push_back((*hpsTancTausDiscrByMediumIsolation)[theTauJetRef]);
    privateData_->thehpsTancTausDiscrByTightIsolation->push_back((*hpsTancTausDiscrByTightIsolation)[theTauJetRef]);
  } else {
    privateData_->thehpsTancTausDiscrByLeadingTrackFinding->push_back(-1);
    privateData_->thehpsTancTausDiscrByLeadingTrackPtCut->push_back(-1);
    privateData_->thehpsTancTausDiscrByLeadingPionPtCut->push_back(-1);
    privateData_->thehpsTancTausDiscrByTanc->push_back(-1);
    privateData_->thehpsTancTausDiscrByTancRaw->push_back(-1);
    privateData_->thehpsTancTausDiscrByTancVLoose->push_back(-1);
    privateData_->thehpsTancTausDiscrByTancLoose->push_back(-1);
    privateData_->thehpsTancTausDiscrByTancMedium->push_back(-1);
    privateData_->thehpsTancTausDiscrByTancTight->push_back(-1);
    privateData_->thehpsTancTausDiscrByLooseElectronRejection->push_back(-1);
    privateData_->thehpsTancTausDiscrByMediumElectronRejection->push_back(-1);
    privateData_->thehpsTancTausDiscrByTightElectronRejection->push_back(-1);
    privateData_->thehpsTancTausDiscrByLooseMuonRejection->push_back(-1);
    privateData_->thehpsTancTausDiscrByTightMuonRejection->push_back(-1);
    privateData_->thehpsTancTausDiscrByDecayModeSelection->push_back(-1);
    privateData_->thehpsTancTausDiscrByVLooseIsolation->push_back(-1);
    privateData_->thehpsTancTausDiscrByLooseIsolation->push_back(-1);
    privateData_->thehpsTancTausDiscrByMediumIsolation->push_back(-1);
    privateData_->thehpsTancTausDiscrByTightIsolation->push_back(-1);
  }
}


void CmsPFTauFiller::writeLeadPFCandInfo(const reco::PFTau *tau, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  privateData_->isolationPFChargedHadrCandsPtSum->push_back( tau->isolationPFChargedHadrCandsPtSum() );
  privateData_->isolationPFGammaCandsEtSum->push_back( tau->isolationPFGammaCandsEtSum() );
  privateData_->emFraction->push_back( tau->emFraction() );
  privateData_->hcalTotOverPLead->push_back( tau->hcalTotOverPLead() );
  privateData_->hcal3x3OverPLead->push_back( tau->hcal3x3OverPLead() );
  privateData_->ecalStripSumEOverPLead->push_back( tau->ecalStripSumEOverPLead() );
  privateData_->bremsRecoveryEOverPLead->push_back( tau->bremsRecoveryEOverPLead() );
}

// Then the tree methods **************************
void CmsPFTauFiller::treePFTauBasicInfo(const std::string &colPrefix, const std::string &colSuffix)
{
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"pfJetIndex"+colSuffix).c_str(), *privateData_->pfJetIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"isNonNull"+colSuffix).c_str(), *privateData_->isNonNull, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"charge"+colSuffix).c_str(), *privateData_->charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"et"+colSuffix).c_str(), *privateData_->et, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"momentum"+colSuffix).c_str(), *privateData_->momentum, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"px"+colSuffix).c_str(), *privateData_->x, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"py"+colSuffix).c_str(), *privateData_->y, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pz"+colSuffix).c_str(), *privateData_->z, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexX"+colSuffix).c_str(), *privateData_->vertexX, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexY"+colSuffix).c_str(), *privateData_->vertexY, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexZ"+colSuffix).c_str(), *privateData_->vertexZ, nCandString.c_str(), 0, "Reco");
}


void CmsPFTauFiller::treePFTauDiscInfo(const std::string &colPrefix, const std::string &colSuffix)
{
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"theTauDiscrByLeadingTrackFinding"+colSuffix).c_str(), *privateData_->theTauDiscrByLeadingTrackFinding, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByLeadingTrackPtCut"+colSuffix).c_str(), *privateData_->theTauDiscrByLeadingTrackPtCut, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByLeadingPionPtCut"+colSuffix).c_str(), *privateData_->theTauDiscrByLeadingPionPtCut, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByIsolation"+colSuffix).c_str(), *privateData_->theTauDiscrByIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByIsolationUsingLeadingPion"+colSuffix).c_str(), *privateData_->theTauDiscrByIsolationUsingLeadingPion, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByTrackIsolation"+colSuffix).c_str(), *privateData_->theTauDiscrByTrackIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByTrackIsolationUsingLeadingPion"+colSuffix).c_str(), *privateData_->theTauDiscrByTrackIsolationUsingLeadingPion, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByECALIsolation"+colSuffix).c_str(), *privateData_->theTauDiscrByECALIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByECALIsolationUsingLeadingPion"+colSuffix).c_str(), *privateData_->theTauDiscrByECALIsolationUsingLeadingPion, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrAgainstMuon"+colSuffix).c_str(), *privateData_->theTauDiscrAgainstMuon, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrAgainstElectron"+colSuffix).c_str(), *privateData_->theTauDiscrAgainstElectron, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByTaNC"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNC, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByTaNCfrHalfPercent"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNCfrHalfPercent, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByTaNCfrOnePercent"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNCfrOnePercent, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByTaNCfrQuarterPercent"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNCfrQuarterPercent, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByTaNCfrTenthPercent"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNCfrTenthPercent, nCandString.c_str(), 0, "Reco");
}

void CmsPFTauFiller::treehpsPFTauDiscInfo(const std::string &colPrefix, const std::string &colSuffix)
{
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"thehpsTauDiscrByLooseElectronRejection"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByLooseElectronRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByMediumElectronRejection"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByMediumElectronRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByTightElectronRejection"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByTightElectronRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByLooseMuonRejection"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByLooseMuonRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByTightMuonRejection"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByTightMuonRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByDecayModeFinding"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByDecayModeFinding, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByVLooseIsolation"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByVLooseIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByLooseIsolation"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByLooseIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByMediumIsolation"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByMediumIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTauDiscrByTightIsolation"+colSuffix).c_str(), *privateData_->thehpsTauDiscrByTightIsolation, nCandString.c_str(), 0, "Reco");
}

void CmsPFTauFiller::treehpsTancTausDiscInfo(const std::string &colPrefix, const std::string &colSuffix)
{
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"thehpsTancTausDiscrByLeadingTrackFinding"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByLeadingTrackFinding, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByLeadingTrackPtCut"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByLeadingTrackPtCut, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByLeadingPionPtCut"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByLeadingPionPtCut, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTanc"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTanc, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTancRaw"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTancRaw, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTancVLoose"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTancVLoose, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTancLoose"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTancLoose, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTancMedium"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTancMedium, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTancTight"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTancTight, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByLooseElectronRejection"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByLooseElectronRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByMediumElectronRejection"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByMediumElectronRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTightElectronRejection"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTightElectronRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByLooseMuonRejection"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByLooseMuonRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTightMuonRejection"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTightMuonRejection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByDecayModeSelection"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByDecayModeSelection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByVLooseIsolation"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByVLooseIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByLooseIsolation"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByLooseIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByMediumIsolation"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByMediumIsolation, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"thehpsTancTausDiscrByTightIsolation"+colSuffix).c_str(), *privateData_->thehpsTancTausDiscrByTightIsolation, nCandString.c_str(), 0, "Reco");
}

void CmsPFTauFiller::treeLeadPFCandInfo(const std::string &colPrefix, const std::string &colSuffix)
{
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"isolationPFChargedHadrCandsPtSum"+colSuffix).c_str(), *privateData_->isolationPFChargedHadrCandsPtSum, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"isolationPFGammaCandsEtSum"+colSuffix).c_str(), *privateData_->isolationPFGammaCandsEtSum, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emFraction"+colSuffix).c_str(), *privateData_->emFraction, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hcalTotOverPLead"+colSuffix).c_str(), *privateData_->hcalTotOverPLead, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hcal3x3OverPLead"+colSuffix).c_str(), *privateData_->hcal3x3OverPLead, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ecalStripSumEOverPLead"+colSuffix).c_str(), *privateData_->ecalStripSumEOverPLead, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"bremsRecoveryEOverPLead"+colSuffix).c_str(), *privateData_->bremsRecoveryEOverPLead, nCandString.c_str(), 0, "Reco");
}

void CmsPFTauFillerData::initialise()
{
  initialiseCandidate();

  // from the basic
  pfJetIndex = new vector<int>;
  isNonNull = new vector<int>;
  // from the Leading Track: for Electrons
  isolationPFChargedHadrCandsPtSum = new vector<float>;
  isolationPFGammaCandsEtSum = new vector<float>;
  emFraction = new vector<float>;
  hcalTotOverPLead = new vector<float>;
  hcal3x3OverPLead = new vector<float>;
  ecalStripSumEOverPLead = new vector<float>;
  bremsRecoveryEOverPLead = new vector<float>;

  // from the Discriminators
  theTauDiscrByLeadingTrackFinding = new vector<float>;
  theTauDiscrByLeadingTrackPtCut = new vector<float>;
  theTauDiscrByLeadingPionPtCut = new vector<float>;
  theTauDiscrByIsolation = new vector<float>;
  theTauDiscrByIsolationUsingLeadingPion = new vector<float>;
  theTauDiscrByTrackIsolation = new vector<float>;
  theTauDiscrByTrackIsolationUsingLeadingPion = new vector<float>;
  theTauDiscrByECALIsolation = new vector<float>;
  theTauDiscrByECALIsolationUsingLeadingPion = new vector<float>;
  theTauDiscrAgainstMuon = new vector<float>;
  theTauDiscrAgainstElectron = new vector<float>;
  theTauDiscrByTaNC = new vector<float>;
  theTauDiscrByTaNCfrHalfPercent = new vector<float>;
  theTauDiscrByTaNCfrOnePercent = new vector<float>;
  theTauDiscrByTaNCfrQuarterPercent = new vector<float>;
  theTauDiscrByTaNCfrTenthPercent = new vector<float>;

  thehpsTauDiscrByLooseElectronRejection = new vector<float>;
  thehpsTauDiscrByMediumElectronRejection = new vector<float>;
  thehpsTauDiscrByTightElectronRejection = new vector<float>;
  thehpsTauDiscrByLooseMuonRejection = new vector<float>;
  thehpsTauDiscrByTightMuonRejection = new vector<float>;
  thehpsTauDiscrByDecayModeFinding = new vector<float>;
  thehpsTauDiscrByVLooseIsolation = new vector<float>;
  thehpsTauDiscrByLooseIsolation = new vector<float>;
  thehpsTauDiscrByMediumIsolation = new vector<float>;
  thehpsTauDiscrByTightIsolation = new vector<float>;

  thehpsTancTausDiscrByLeadingTrackFinding = new vector<float>;
  thehpsTancTausDiscrByLeadingTrackPtCut = new vector<float>;
  thehpsTancTausDiscrByLeadingPionPtCut = new vector<float>;
  thehpsTancTausDiscrByTanc = new vector<float>;
  thehpsTancTausDiscrByTancRaw = new vector<float>;
  thehpsTancTausDiscrByTancVLoose = new vector<float>;
  thehpsTancTausDiscrByTancLoose = new vector<float>;
  thehpsTancTausDiscrByTancMedium = new vector<float>;
  thehpsTancTausDiscrByTancTight = new vector<float>;
  thehpsTancTausDiscrByLooseElectronRejection = new vector<float>;
  thehpsTancTausDiscrByMediumElectronRejection = new vector<float>;
  thehpsTancTausDiscrByTightElectronRejection = new vector<float>;
  thehpsTancTausDiscrByLooseMuonRejection = new vector<float>;
  thehpsTancTausDiscrByTightMuonRejection = new vector<float>;
  thehpsTancTausDiscrByDecayModeSelection = new vector<float>;
  thehpsTancTausDiscrByVLooseIsolation = new vector<float>;
  thehpsTancTausDiscrByLooseIsolation = new vector<float>;
  thehpsTancTausDiscrByMediumIsolation = new vector<float>;
  thehpsTancTausDiscrByTightIsolation = new vector<float>;
}


void CmsPFTauFillerData::clear()
{
  // from the Basic
  pfJetIndex->clear();

  // from the Leading Track: for Electrons
  isolationPFChargedHadrCandsPtSum->clear();
  isolationPFGammaCandsEtSum->clear();
  emFraction->clear();
  hcalTotOverPLead->clear();
  hcal3x3OverPLead->clear();
  ecalStripSumEOverPLead->clear();
  bremsRecoveryEOverPLead->clear();

  // from the Discriminators
  theTauDiscrByLeadingTrackFinding->clear();
  theTauDiscrByLeadingTrackPtCut->clear();
  theTauDiscrByLeadingPionPtCut->clear();
  theTauDiscrByIsolation->clear();
  theTauDiscrByIsolationUsingLeadingPion->clear();
  theTauDiscrByTrackIsolation->clear();
  theTauDiscrByTrackIsolationUsingLeadingPion->clear();
  theTauDiscrByECALIsolation->clear();
  theTauDiscrByECALIsolationUsingLeadingPion->clear();
  theTauDiscrAgainstMuon->clear();
  theTauDiscrAgainstElectron->clear();
  theTauDiscrByTaNC->clear();
  theTauDiscrByTaNCfrHalfPercent->clear();
  theTauDiscrByTaNCfrOnePercent->clear();
  theTauDiscrByTaNCfrQuarterPercent->clear();
  theTauDiscrByTaNCfrTenthPercent->clear();

  thehpsTauDiscrByLooseElectronRejection->clear();
  thehpsTauDiscrByMediumElectronRejection->clear();
  thehpsTauDiscrByTightElectronRejection->clear();
  thehpsTauDiscrByLooseMuonRejection->clear();
  thehpsTauDiscrByTightMuonRejection->clear();
  thehpsTauDiscrByDecayModeFinding->clear();
  thehpsTauDiscrByVLooseIsolation->clear();
  thehpsTauDiscrByLooseIsolation->clear();
  thehpsTauDiscrByMediumIsolation->clear();
  thehpsTauDiscrByTightIsolation->clear();

  thehpsTancTausDiscrByLeadingTrackFinding->clear();
  thehpsTancTausDiscrByLeadingTrackPtCut->clear();
  thehpsTancTausDiscrByLeadingPionPtCut->clear();
  thehpsTancTausDiscrByTanc->clear();
  thehpsTancTausDiscrByTancRaw->clear();
  thehpsTancTausDiscrByTancVLoose->clear();
  thehpsTancTausDiscrByTancLoose->clear();
  thehpsTancTausDiscrByTancMedium->clear();
  thehpsTancTausDiscrByTancTight->clear();
  thehpsTancTausDiscrByLooseElectronRejection->clear();
  thehpsTancTausDiscrByMediumElectronRejection->clear();
  thehpsTancTausDiscrByTightElectronRejection->clear();
  thehpsTancTausDiscrByLooseMuonRejection->clear();
  thehpsTancTausDiscrByTightMuonRejection->clear();
  thehpsTancTausDiscrByDecayModeSelection->clear();
  thehpsTancTausDiscrByVLooseIsolation->clear();
  thehpsTancTausDiscrByLooseIsolation->clear();
  thehpsTancTausDiscrByMediumIsolation->clear();
  thehpsTancTausDiscrByTightIsolation->clear();

}
