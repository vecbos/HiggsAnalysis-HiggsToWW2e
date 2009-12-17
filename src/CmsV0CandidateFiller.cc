//---------------------------------------------------------------------------
//
// Description:
//       Package:   HiggsAnalysis/HiggsToWW2e
//       Class:     CmsMuonFiller
//
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsV0CandidateFiller.h"

#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsV0CandidateFiller::CmsV0CandidateFiller(CmsTree *cmsTree, 
			     int maxTracks, int maxMCTracks,
			     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsV0CandidateFillerData)
{
  cmstree=cmsTree;

  saveTrk_=true;

  trkIndexName_ = new std::string("n");

  maxTracks_=maxTracks;

  privateData_->initialise();
  
}

CmsV0CandidateFiller::CmsV0CandidateFiller(CmsTree *cmsTree, 
			     bool fatTree, 
			     int maxTracks, int maxMCTracks,
			     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsV0CandidateFillerData)
{
  cmstree=cmsTree;

  saveTrk_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsV0CandidateFiller::~CmsV0CandidateFiller() { 
  
  // Index of daughters in the track collection
  delete privateData_->dau1Index;
  delete privateData_->dau2Index;

  delete privateData_->ncand;

}


//-------------
// Methods   --
//-------------

void CmsV0CandidateFiller::writeCollectionToTree(edm::InputTag collectionTag,
					 const edm::Event& iEvent, const edm::EventSetup& iSetup,
					 const std::string &columnPrefix, const std::string &columnSuffix,
					 bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsV0CandidateFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsV0CandidateFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ << " and no output flag is set."
			       << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsV0CandidateFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ 
			       << ". Collection will be truncated ";
    }
  
    *(privateData_->ncand) = collection->size();

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);

      // fill the indices of the daughter tracks
      writeTrkInfo(&(*cand),iEvent,iSetup);

    }
  }
  else {
  *(privateData_->ncand) = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  treeTrkInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();
	
}


void CmsV0CandidateFiller::writeTrkInfo(const Candidate *cand, 
				 const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  TrackRef trk1Ref = (dynamic_cast<const reco::RecoCandidate*>(cand->daughter(0)))->track();
  TrackRef trk2Ref = (dynamic_cast<const reco::RecoCandidate*>(cand->daughter(1)))->track();

  if ( trk1Ref.isNonnull() ) {
    privateData_->dau1Index->push_back( trk1Ref.key() );
    privateData_->dau2Index->push_back( trk2Ref.key() );
  } else {
    privateData_->dau1Index->push_back( -1 );
    privateData_->dau2Index->push_back( -1 );
  }

}

void CmsV0CandidateFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  
  cmstree->column((colPrefix+"dau1Index"+colSuffix).c_str(), *privateData_->dau1Index, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dau2Index"+colSuffix).c_str(), *privateData_->dau2Index, nCandString.c_str(), 0, "Reco");

}

void CmsV0CandidateFillerData::initialise() {

  initialiseCandidate();
  dau1Index = new vector<int>;
  dau2Index = new vector<int>;

}

void CmsV0CandidateFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  dau1Index->clear();
  dau2Index->clear();
  
}
