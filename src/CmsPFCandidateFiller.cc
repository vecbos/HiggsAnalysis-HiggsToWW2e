//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsPFCandidateFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  9 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFCandidateFiller.h"
#include <TMath.h>
#include <TTree.h>
#include <Math/VectorUtil.h>

#include <string>

using namespace edm;
using namespace reco;




//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsPFCandidateFiller::CmsPFCandidateFiller(CmsTree *cmsTree, int maxTracks, 
                                           int maxMCTracks,
                                           bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFCandidateFillerData)
{
  cmstree=cmsTree;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();

  minPtPFCand_ = 10.; 
  chargeOnly_ = true;
}

void CmsPFCandidateFiller::setMinPtPFCand(double value) {
  minPtPFCand_ = value;
}

void CmsPFCandidateFiller::dumpChargeOnly(bool value) {
  chargeOnly_ = value;
}

//--------------
// Destructor --
//--------------

CmsPFCandidateFiller::~CmsPFCandidateFiller() {

  // delete here the vector ptr's
  delete privateData_->particleId;
  delete privateData_->pfCandChargedIso03;
  delete privateData_->pfCandNeutralIso03;
  delete privateData_->pfCandPhotonIso03;
  delete privateData_->pfIsolationSumPUPtR03;

  delete privateData_;
  
}


//-------------
// Methods   --
//-------------




void CmsPFCandidateFiller::writeCollectionToTree(edm::InputTag collectionTagNoPU, edm::InputTag collectionTagPU,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {

  // list of PF candidates from PV
  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTagNoPU, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFCandidateFiller") << "Can't get candidate collection: " << collectionTagNoPU; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  // list of PF candidates from other vertices
  edm::Handle< edm::View<reco::Candidate> > collectionHandlePU;
  try { iEvent.getByLabel(collectionTagPU, collectionHandlePU); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFCandidateFiller") << "Can't get candidate collection: " << collectionTagPU; }
  const edm::View<reco::Candidate> *collectionPU = collectionHandlePU.product();

  privateData_->clearTrkVectors();

  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  LogInfo("CmsPFCandidateFiller") << "=== Writing collection " << nCandString << " ===";

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      LogError("CmsPFCandidateFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      LogError("CmsPFCandidateFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ 
				     << ". Collection will be truncated ";
    }
  
    // for track link
    try { iEvent.getByLabel(generalTracks_, h_tracks); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get general track collection: " << generalTracks_; }

    edm::View<reco::Candidate>::const_iterator cand;
    edm::View<reco::Candidate>::const_iterator cand2;

    for(cand=collection->begin(); cand!=collection->end(); cand++) {

      const reco::PFCandidate *pfcand = dynamic_cast< const reco::PFCandidate * > ( &(*cand) );
      // check if particle passes the kinematic requirement
      if(pfcand->pt() < minPtPFCand_) continue; 
      // if requested, ignore neutrals
      if(chargeOnly_ && pfcand->particleId() != reco::PFCandidate::h
	 && pfcand->particleId() != reco::PFCandidate::e
	 && pfcand->particleId() != reco::PFCandidate::mu) continue;

      if(saveCand_) {
	// fill basic kinematics
	writeCandInfo(&(*cand),iEvent,iSetup);
	// fill isolation: loop over other PF candidates and
	// compute variables for PF isolation
	double ptCh=0;
	double ptNeu=0;
	double ptPh=0;
	
	for(cand2=collection->begin(); cand2!=collection->end(); cand2++) {
	  const reco::PFCandidate *pfcand2 = dynamic_cast< const reco::PFCandidate * > ( &(*cand2) );
	  // ignore PF cand if outside the DR=0.3 cone 
	  Double_t dr = ROOT::Math::VectorUtil::DeltaR(pfcand->momentum(), pfcand2->momentum());
	  if(dr >= 0.3) continue;
	  // ignore the PF candidate if it angularly matches the original cand
	  if(dr <=0.0) continue;
	  
	  // photons
	  if (pfcand2->particleId() == reco::PFCandidate::gamma) ptPh += pfcand2->pt();
	  // charged hadrons	  
	  else if(pfcand2->particleId() == reco::PFCandidate::h) ptCh += pfcand2->pt(); 
	  // neutral hadrons
	  else if(pfcand2->particleId() == reco::PFCandidate::h0) ptNeu += pfcand2->pt();
	}
	privateData_->particleId->push_back(pfcand->particleId());
	privateData_->pfCandChargedIso03->push_back(ptCh);
	privateData_->pfCandNeutralIso03->push_back(ptNeu);
	privateData_->pfCandPhotonIso03->push_back(ptPh);

	double ptPU=0;
	for(cand2=collectionPU->begin(); cand2!=collectionPU->end(); cand2++) {
	  const reco::PFCandidate *pfcand2 = dynamic_cast< const reco::PFCandidate * > ( &(*cand2) );
	  // ignore PF cand if outside the DR=0.3 cone 
	  Double_t dr = ROOT::Math::VectorUtil::DeltaR(pfcand->momentum(), pfcand2->momentum());
	  if(dr >= 0.3) continue;
	  // ignore the PF candidate if it angularly matches the original cand
	  if(dr <=0.0) continue;
	  // charged hadrons
	  if(pfcand2->particleId() == reco::PFCandidate::h) ptPU += pfcand2->pt();
	}
	privateData_->pfIsolationSumPUPtR03->push_back(ptPU);
      } 
    }      
    *(privateData_->ncand) = privateData_->particleId->size();
  } else {
    *(privateData_->ncand) = 0;
  }

  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 

  int blockSize = privateData_->particleId->size();
  
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");

  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  treePFInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();

  delete trkIndexName_;

}

void CmsPFCandidateFiller::treePFInfo(const std::string colPrefix, const std::string colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"particleId"+colSuffix).c_str(), *privateData_->particleId, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ChargedIso03"+colSuffix).c_str(), *privateData_->pfCandChargedIso03,nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"NeutralIso03"+colSuffix).c_str(), *privateData_->pfCandNeutralIso03,nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PhotonIso03"+colSuffix).c_str(), *privateData_->pfCandPhotonIso03,nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"SumPUPtR03"+colSuffix).c_str(), *privateData_->pfIsolationSumPUPtR03,nCandString.c_str(), 0, "Reco");

}

void CmsPFCandidateFillerData::initialise() {

  initialiseCandidate();
  particleId = new vector<int>;
  pfCandChargedIso03 = new vector<double>;
  pfCandNeutralIso03 = new vector<double>;
  pfCandPhotonIso03 = new vector<double>;
  pfIsolationSumPUPtR03 = new vector<double>;

}

void CmsPFCandidateFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  particleId->clear();
  pfCandChargedIso03->clear();
  pfCandNeutralIso03->clear();
  pfCandPhotonIso03->clear();
  pfIsolationSumPUPtR03->clear();

}
