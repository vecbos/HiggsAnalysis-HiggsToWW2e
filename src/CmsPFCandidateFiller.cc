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
}


//--------------
// Destructor --
//--------------

CmsPFCandidateFiller::~CmsPFCandidateFiller() {

  // delete here the vector ptr's
  delete privateData_->particleId;
  delete privateData_;
}


//-------------
// Methods   --
//-------------




void CmsPFCandidateFiller::writeCollectionToTree(edm::InputTag collectionTag,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFCandidateFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

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
  
    *(privateData_->ncand) = collection->size();

    // for track link
    try { iEvent.getByLabel(generalTracks_, h_tracks); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get general track collection: " << generalTracks_; }

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
      const reco::PFCandidate *pfcand = dynamic_cast< const reco::PFCandidate * > ( &(*cand) );
      privateData_->particleId->push_back(pfcand->particleId());

    }
  }
  else {
    *(privateData_->ncand) = 0;
  }

  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 

  int blockSize = (collection) ? collection->size() : 0;
  
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");

  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  treePFInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();

  delete trkIndexName_;

}







void CmsPFCandidateFiller::treePFInfo(const std::string colPrefix, const std::string colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"particleId"+colSuffix).c_str(), *privateData_->particleId, nCandString.c_str(), 0, "Reco");
}







void CmsPFCandidateFillerData::initialise() {

  initialiseCandidate();
  particleId = new vector<int>;

}

void CmsPFCandidateFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  particleId->clear();

}
