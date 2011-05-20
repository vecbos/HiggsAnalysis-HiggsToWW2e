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

  delete privateData_->isNonNull;
  delete privateData_->charge;
  delete privateData_->energy;
  delete privateData_->et;
  delete privateData_->momentum;
  delete privateData_->theta;
  delete privateData_->pt;
  delete privateData_->eta;
  delete privateData_->phi;
  delete privateData_->x;
  delete privateData_->y;
  delete privateData_->z;
  delete privateData_->vertexX;
  delete privateData_->vertexY;
  delete privateData_->vertexZ;
  delete privateData_->mass;
  delete privateData_->mt;
  delete privateData_->pdgId;
  delete privateData_->nDau;


  // Leading Track: for Electrons
  delete privateData_->isolationPFChargedHadrCandsPtSum;
  delete privateData_->isolationPFGammaCandsEtSum;
  delete privateData_->emFraction;
  delete privateData_->hcalTotOverPLead;
  delete privateData_->hcal3x3OverPLead;
  delete privateData_->ecalStripSumEOverPLead;
  delete privateData_->bremsRecoveryEOverPLead;

  delete privateData_->theTauDiscrByLeadTrackFinding;
  delete privateData_->theTauDiscrByLeadTrackPtCut;
  // delete privateData_->theTauDiscrByNProngs;
  delete privateData_->theTauDiscrByTrackIso;
  delete privateData_->theTauDiscrByEcalIso;
  delete privateData_->theTauDiscrAgainstMuons;
  delete privateData_->theTauDiscrAgainstElectrons;
//   delete privateData_->theTauDiscrByTaNC;
//   delete privateData_->theTauDiscrByTaNCfrHalfPercent;
//   delete privateData_->theTauDiscrByTaNCfrOnePercent;
//   delete privateData_->theTauDiscrByTaNCfrQuarterPercent;
//   delete privateData_->theTauDiscrByTaNCfrTenthPercent;

  delete privateData_->ncand;
}


//-------------
// Methods   --
//-------------

void CmsPFTauFiller::savePFTauBasic(bool what) { savePFTauBasic_=what; }

void CmsPFTauFiller::saveLeadPFCand(bool what) { saveLeadPFCand_=what; }

void CmsPFTauFiller::writeCollectionToTree(edm::InputTag collectionTag,
					   const edm::Event& iEvent, const edm::EventSetup& iSetup,
					   const std::string &columnPrefix, const std::string &columnSuffix,
					   edm::InputTag tauDiscrByLeadTrackFindingTag,
					   edm::InputTag tauDiscrByLeadTrackPtCutTag,
					   // edm::InputTag tauDiscrByNProngsTag,
					   edm::InputTag tauDiscrByTrackIsoTag,
					   edm::InputTag tauDiscrByEcalIsoTag, 
					   edm::InputTag tauDiscrAgainstMuonsTag,
					   edm::InputTag tauDiscrAgainstElectronsTag,
// 					   edm::InputTag tauDiscrByTaNCTag,
// 					   edm::InputTag tauDiscrByTaNCfrHalfPercentTag,
// 					   edm::InputTag tauDiscrByTaNCfrOnePercentTag,
// 					   edm::InputTag tauDiscrByTaNCfrQuarterPercentTag,
// 					   edm::InputTag tauDiscrByTaNCfrTenthPercentTag,
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
	edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadTrackFindingHandle;
	try { iEvent.getByLabel(tauDiscrByLeadTrackFindingTag, tauDiscrByLeadTrackFindingHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByLeadTrackFindingTag; }
	const reco::PFTauDiscriminator *tauDiscrByLeadTrackFinding = tauDiscrByLeadTrackFindingHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByLeadTrackPtCutHandle;
	try { iEvent.getByLabel(tauDiscrByLeadTrackPtCutTag, tauDiscrByLeadTrackPtCutHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByLeadTrackPtCutTag; }
	const reco::PFTauDiscriminator *tauDiscrByLeadTrackPtCut = tauDiscrByLeadTrackPtCutHandle.product();

	// edm::Handle<reco::PFTauDiscriminator> tauDiscrByNProngsHandle;
	// try { iEvent.getByLabel(tauDiscrByNProngsTag, tauDiscrByNProngsHandle); }
	// catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByNProngsTag; }
	// const reco::PFTauDiscriminator *tauDiscrByNProngs = tauDiscrByNProngsHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTrackIsoHandle;
	try { iEvent.getByLabel(tauDiscrByTrackIsoTag, tauDiscrByTrackIsoHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTrackIsoTag; }
	const reco::PFTauDiscriminator *tauDiscrByTrackIso = tauDiscrByTrackIsoHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrByEcalIsoHandle;
	try { iEvent.getByLabel(tauDiscrByEcalIsoTag, tauDiscrByEcalIsoHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByEcalIsoTag; }
	const reco::PFTauDiscriminator *tauDiscrByEcalIso = tauDiscrByEcalIsoHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrAgainstMuonsHandle;
	try { iEvent.getByLabel(tauDiscrAgainstMuonsTag, tauDiscrAgainstMuonsHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrAgainstMuonsTag; }
	const reco::PFTauDiscriminator *tauDiscrAgainstMuons = tauDiscrAgainstMuonsHandle.product();

	edm::Handle<reco::PFTauDiscriminator> tauDiscrAgainstElectronsHandle;
	try { iEvent.getByLabel(tauDiscrAgainstElectronsTag, tauDiscrAgainstElectronsHandle); }
	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrAgainstElectronsTag; }
	const reco::PFTauDiscriminator *tauDiscrAgainstElectrons = tauDiscrAgainstElectronsHandle.product();

// 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCHandle;
// 	try { iEvent.getByLabel(tauDiscrByTaNCTag, tauDiscrByTaNCHandle); }
// 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCTag; }
// 	const reco::PFTauDiscriminator *tauDiscrByTaNC = tauDiscrByTaNCHandle.product();

// 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCfrHalfPercentHandle;
// 	try { iEvent.getByLabel(tauDiscrByTaNCfrHalfPercentTag, tauDiscrByTaNCfrHalfPercentHandle); }
// 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCfrHalfPercentTag; }
// 	const reco::PFTauDiscriminator *tauDiscrByTaNCfrHalfPercent = tauDiscrByTaNCfrHalfPercentHandle.product();

// 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCfrOnePercentHandle;
// 	try { iEvent.getByLabel(tauDiscrByTaNCfrOnePercentTag, tauDiscrByTaNCfrOnePercentHandle); }
// 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCfrOnePercentTag; }
// 	const reco::PFTauDiscriminator *tauDiscrByTaNCfrOnePercent = tauDiscrByTaNCfrOnePercentHandle.product();

// 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCfrQuarterPercentHandle;
// 	try { iEvent.getByLabel(tauDiscrByTaNCfrQuarterPercentTag, tauDiscrByTaNCfrQuarterPercentHandle); }
// 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCfrQuarterPercentTag; }
// 	const reco::PFTauDiscriminator *tauDiscrByTaNCfrQuarterPercent = tauDiscrByTaNCfrQuarterPercentHandle.product();

// 	edm::Handle<reco::PFTauDiscriminator> tauDiscrByTaNCfrTenthPercentHandle;
// 	try { iEvent.getByLabel(tauDiscrByTaNCfrTenthPercentTag, tauDiscrByTaNCfrTenthPercentHandle); }
// 	catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFTauFiller") << "Can't get PFTau discriminator: " << tauDiscrByTaNCfrTenthPercentTag; }
// 	const reco::PFTauDiscriminator *tauDiscrByTaNCfrTenthPercent = tauDiscrByTaNCfrTenthPercentHandle.product();

	writePFTauDiscInfo(collectionHandle, iTauJet,
			   tauDiscrByLeadTrackFinding, tauDiscrByLeadTrackPtCut,// tauDiscrByNProngs, 
			   tauDiscrByTrackIso, tauDiscrByEcalIso, tauDiscrAgainstMuons, tauDiscrAgainstElectrons
// 			   tauDiscrByTaNC,
// 			   tauDiscrByTaNCfrHalfPercent, tauDiscrByTaNCfrOnePercent,
// 			   tauDiscrByTaNCfrQuarterPercent, tauDiscrByTaNCfrTenthPercent
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

// First the write methods **************************
void CmsPFTauFiller::writePFTauBasicInfo(const reco::PFTau *tau, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  privateData_->isNonNull->push_back((int)tau->leadPFChargedHadrCand().isNonnull());
  privateData_->charge->push_back((int)tau->charge());
  privateData_->energy->push_back(tau->energy());
  privateData_->et->push_back(tau->et());
  privateData_->momentum->push_back(tau->p());
  privateData_->theta->push_back(tau->theta());
  privateData_->pt->push_back(tau->pt());
  privateData_->eta->push_back(tau->eta());
  privateData_->phi->push_back(tau->phi());
  privateData_->x->push_back(tau->momentum().x());
  privateData_->y->push_back(tau->momentum().y());
  privateData_->z->push_back(tau->momentum().z());
  privateData_->vertexX->push_back(tau->vx());
  privateData_->vertexY->push_back(tau->vy());
  privateData_->vertexZ->push_back(tau->vz());
  privateData_->mass->push_back(tau->mass());
  privateData_->mt->push_back(tau->mt());
  privateData_->pdgId->push_back(tau->pdgId());
  privateData_->nDau->push_back(tau->numberOfDaughters());
}

void CmsPFTauFiller::writePFTauDiscInfo(edm::Handle<reco::PFTauCollection> tauCollection, int theTauJetIndex,
					const reco::PFTauDiscriminator *tauDiscrByLeadTrackFinding,
					const reco::PFTauDiscriminator *tauDiscrByLeadTrackPtCut,
					// const reco::PFTauDiscriminator *tauDiscrByNProngs,
					const reco::PFTauDiscriminator *tauDiscrByTrackIso,
					const reco::PFTauDiscriminator *tauDiscrByEcalIso,
					const reco::PFTauDiscriminator *tauDiscrAgainstMuons,
					const reco::PFTauDiscriminator *tauDiscrAgainstElectrons
// 					const reco::PFTauDiscriminator *tauDiscrByTaNC,
// 					const reco::PFTauDiscriminator *tauDiscrByTaNCfrHalfPercent,
// 					const reco::PFTauDiscriminator *tauDiscrByTaNCfrOnePercent,
// 					const reco::PFTauDiscriminator *tauDiscrByTaNCfrQuarterPercent,
// 					const reco::PFTauDiscriminator *tauDiscrByTaNCfrTenthPercent
                                        )
{
  if ( theTauJetIndex != -1 ) {
    reco::PFTauRef theTauJetRef(tauCollection, theTauJetIndex);
    if(theTauJetRef.isNonnull()) {
    privateData_->theTauDiscrByLeadTrackFinding->push_back((*tauDiscrByLeadTrackFinding)[theTauJetRef]);
    privateData_->theTauDiscrByLeadTrackPtCut->push_back((*tauDiscrByLeadTrackPtCut)[theTauJetRef]);
    // privateData_->theTauDiscrByNProngs->push_back((*tauDiscrByNProngs)[theTauJetRef]);
    privateData_->theTauDiscrByTrackIso->push_back((*tauDiscrByTrackIso)[theTauJetRef]);
    privateData_->theTauDiscrByEcalIso->push_back((*tauDiscrByEcalIso)[theTauJetRef]);
    privateData_->theTauDiscrAgainstMuons->push_back((*tauDiscrAgainstMuons)[theTauJetRef]);
    privateData_->theTauDiscrAgainstElectrons->push_back((*tauDiscrAgainstElectrons)[theTauJetRef]);
//     privateData_->theTauDiscrByTaNC->push_back((*tauDiscrByTaNC)[theTauJetRef]);
//     privateData_->theTauDiscrByTaNCfrHalfPercent->push_back((*tauDiscrByTaNCfrHalfPercent)[theTauJetRef]);
//     privateData_->theTauDiscrByTaNCfrOnePercent->push_back((*tauDiscrByTaNCfrOnePercent)[theTauJetRef]);
//     privateData_->theTauDiscrByTaNCfrQuarterPercent->push_back((*tauDiscrByTaNCfrQuarterPercent)[theTauJetRef]);
//     privateData_->theTauDiscrByTaNCfrTenthPercent->push_back((*tauDiscrByTaNCfrTenthPercent)[theTauJetRef]);
    } else {
      privateData_->theTauDiscrByLeadTrackFinding->push_back(-1);
      privateData_->theTauDiscrByLeadTrackPtCut->push_back(-1);
      // privateData_->theTauDiscrByNProngs->push_back(-1);
      privateData_->theTauDiscrByTrackIso->push_back(-1);
      privateData_->theTauDiscrByEcalIso->push_back(-1);
      privateData_->theTauDiscrAgainstMuons->push_back(-1);
      privateData_->theTauDiscrAgainstElectrons->push_back(-1);
      //     privateData_->theTauDiscrByTaNC->push_back(-1);
      //     privateData_->theTauDiscrByTaNCfrHalfPercent->push_back(-1);
      //     privateData_->theTauDiscrByTaNCfrOnePercent->push_back(-1);
      //     privateData_->theTauDiscrByTaNCfrQuarterPercent->push_back(-1);
      //     privateData_->theTauDiscrByTaNCfrTenthPercent->push_back(-1);      
    }
  } else {
    privateData_->theTauDiscrByLeadTrackFinding->push_back(-1);
    privateData_->theTauDiscrByLeadTrackPtCut->push_back(-1);
    // privateData_->theTauDiscrByNProngs->push_back(-1);
    privateData_->theTauDiscrByTrackIso->push_back(-1);
    privateData_->theTauDiscrByEcalIso->push_back(-1);
    privateData_->theTauDiscrAgainstMuons->push_back(-1);
    privateData_->theTauDiscrAgainstElectrons->push_back(-1);
//     privateData_->theTauDiscrByTaNC->push_back(-1);
//     privateData_->theTauDiscrByTaNCfrHalfPercent->push_back(-1);
//     privateData_->theTauDiscrByTaNCfrOnePercent->push_back(-1);
//     privateData_->theTauDiscrByTaNCfrQuarterPercent->push_back(-1);
//     privateData_->theTauDiscrByTaNCfrTenthPercent->push_back(-1);
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
  cmstree->column((colPrefix+"isNonNull"+colSuffix).c_str(), *privateData_->isNonNull, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"charge"+colSuffix).c_str(), *privateData_->charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"et"+colSuffix).c_str(), *privateData_->et, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"momentum"+colSuffix).c_str(), *privateData_->momentum, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pt"+colSuffix).c_str(), *privateData_->pt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"px"+colSuffix).c_str(), *privateData_->x, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"py"+colSuffix).c_str(), *privateData_->y, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pz"+colSuffix).c_str(), *privateData_->z, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexX"+colSuffix).c_str(), *privateData_->vertexX, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexY"+colSuffix).c_str(), *privateData_->vertexY, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexZ"+colSuffix).c_str(), *privateData_->vertexZ, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mass"+colSuffix).c_str(), *privateData_->mass, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mt"+colSuffix).c_str(), *privateData_->mt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pdgId"+colSuffix).c_str(), *privateData_->pdgId, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nDau"+colSuffix).c_str(), *privateData_->nDau, nCandString.c_str(), 0, "Reco");
}


void CmsPFTauFiller::treePFTauDiscInfo(const std::string &colPrefix, const std::string &colSuffix)
{
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"theTauDiscrByLeadTrackFinding"+colSuffix).c_str(), *privateData_->theTauDiscrByLeadTrackFinding, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByLeadTrackPtCut"+colSuffix).c_str(), *privateData_->theTauDiscrByLeadTrackPtCut, nCandString.c_str(), 0, "Reco");
  // cmstree->column((colPrefix+"theTauDiscrByNProngs"+colSuffix).c_str(), *privateData_->theTauDiscrByNProngs, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByTrackIso"+colSuffix).c_str(), *privateData_->theTauDiscrByTrackIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrByEcalIso"+colSuffix).c_str(), *privateData_->theTauDiscrByEcalIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrAgainstMuons"+colSuffix).c_str(), *privateData_->theTauDiscrAgainstMuons, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theTauDiscrAgainstElectrons"+colSuffix).c_str(), *privateData_->theTauDiscrAgainstElectrons, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"theTauDiscrByTaNC"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNC, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"theTauDiscrByTaNCfrHalfPercent"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNCfrHalfPercent, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"theTauDiscrByTaNCfrOnePercent"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNCfrOnePercent, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"theTauDiscrByTaNCfrQuarterPercent"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNCfrQuarterPercent, nCandString.c_str(), 0, "Reco");
//   cmstree->column((colPrefix+"theTauDiscrByTaNCfrTenthPercent"+colSuffix).c_str(), *privateData_->theTauDiscrByTaNCfrTenthPercent, nCandString.c_str(), 0, "Reco");
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
  isNonNull = new vector<int>;
  charge = new vector<int>;
  energy = new vector<float>;
  et = new vector<float>;
  momentum = new vector<float>;
  theta = new vector<float>;
  pt = new vector<float>;
  eta = new vector<float>;
  phi = new vector<float>;
  x = new vector<float>;
  y = new vector<float>;
  z = new vector<float>;
  vertexX = new vector<float>;
  vertexY = new vector<float>;
  vertexZ = new vector<float>;
  mass = new vector<float>;
  mt = new vector<float>;
  pdgId = new vector<int>;
  ncand = new int;
  nDau = new vector<int>;

  // from the Leading Track: for Electrons
  isolationPFChargedHadrCandsPtSum = new vector<float>;
  isolationPFGammaCandsEtSum = new vector<float>;
  emFraction = new vector<float>;
  hcalTotOverPLead = new vector<float>;
  hcal3x3OverPLead = new vector<float>;
  ecalStripSumEOverPLead = new vector<float>;
  bremsRecoveryEOverPLead = new vector<float>;

  // from the Discriminators
  theTauDiscrByLeadTrackFinding = new vector<float>;
  theTauDiscrByLeadTrackPtCut = new vector<float>;
  // theTauDiscrByNProngs = new vector<float>;
  theTauDiscrByTrackIso = new vector<float>;
  theTauDiscrByEcalIso = new vector<float>;
  theTauDiscrAgainstMuons = new vector<float>;
  theTauDiscrAgainstElectrons = new vector<float>;
//   theTauDiscrByTaNC = new vector<float>;
//   theTauDiscrByTaNCfrHalfPercent = new vector<float>;
//   theTauDiscrByTaNCfrOnePercent = new vector<float>;
//   theTauDiscrByTaNCfrQuarterPercent = new vector<float>;
//   theTauDiscrByTaNCfrTenthPercent = new vector<float>;
}


void CmsPFTauFillerData::clear()
{
  // from the Basic
  isNonNull->clear();
  charge->clear();
  energy->clear();
  et->clear();
  momentum->clear();
  theta->clear();
  pt->clear();
  eta->clear();
  phi->clear();
  x->clear();
  y->clear();
  z->clear();
  vertexX->clear();
  vertexY->clear();
  vertexZ->clear();
  mass->clear();
  mt->clear();
  pdgId->clear();
  nDau->clear();

  // from the Leading Track: for Electrons
  isolationPFChargedHadrCandsPtSum->clear();
  isolationPFGammaCandsEtSum->clear();
  emFraction->clear();
  hcalTotOverPLead->clear();
  hcal3x3OverPLead->clear();
  ecalStripSumEOverPLead->clear();
  bremsRecoveryEOverPLead->clear();

  // from the Discriminators
  theTauDiscrByLeadTrackFinding->clear();
  theTauDiscrByLeadTrackPtCut->clear();
  // theTauDiscrByNProngs->clear();
  theTauDiscrByTrackIso->clear();
  theTauDiscrByEcalIso->clear();
  theTauDiscrAgainstMuons->clear();
  theTauDiscrAgainstElectrons->clear();
//   theTauDiscrByTaNC->clear();
//   theTauDiscrByTaNCfrHalfPercent->clear();
//   theTauDiscrByTaNCfrOnePercent->clear();
//   theTauDiscrByTaNCfrQuarterPercent->clear();
//   theTauDiscrByTaNCfrTenthPercent->clear();
}
