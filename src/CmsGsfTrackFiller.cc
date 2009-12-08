//---------------------------------------------------------------------------
//
// Description:
//       Package:   HiggsAnalysis/HiggsToWW2e
//       Class:     CmsGsfTrackFiller
//
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsGsfTrackFiller.h"

using namespace edm;
using namespace reco;

CmsGsfTrackFiller::CmsGsfTrackFiller(CmsTree *cmsTree, 
                                     edm::InputTag vertexCollection,
                                     int maxTracks, int maxMCTracks,
                                     bool noOutputIfLimitsReached):
  CmsTrackFiller(cmsTree,vertexCollection,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsGsfTrackFillerData)
{
  cmstree=cmsTree;
  vertexCollection_=vertexCollection;
  
  saveTrk_=true;
  saveFatTrk_=true;
  saveVtxTrk_=true; // to change when bug is fixed
  saveDeDx_=false;


  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialiseGsf();
  x0 = 0.;
  y0 = 0.;
  z0 = 0.;
}

CmsGsfTrackFiller::CmsGsfTrackFiller(CmsTree *cmsTree, 
                                     edm::InputTag vertexCollection,
                                     bool fatTree, 
                                     int maxTracks, int maxMCTracks,
                                     bool noOutputIfLimitsReached, bool vtxtrack):
  CmsTrackFiller(cmsTree,vertexCollection,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsGsfTrackFillerData)
{
  cmstree=cmsTree;
  vertexCollection_=vertexCollection;

  saveTrk_ = true;
  saveVtxTrk_ = vtxtrack; 
  saveFatTrk_ = fatTree;
  saveDeDx_=false;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialiseGsf();
  x0 = 0.;
  y0 = 0.;
  z0 = 0.;
}

CmsGsfTrackFiller::~CmsGsfTrackFiller() {

  // Gsf Track information
  delete privateData_->chargeMode;
  delete privateData_->pxMode;
  delete privateData_->pyMode;
  delete privateData_->pzMode;

}

void CmsGsfTrackFiller::writeCollectionToTree(edm::InputTag collectionTag,
                                              const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                              const std::string &columnPrefix, const std::string &columnSuffix,
                                              bool dumpData) {

  edm::Handle< edm::View<reco::Track> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsGsfTrackFiller") << "Can't get GSF track collection: " << collectionTag; }
  const edm::View<reco::Track> *collection = collectionHandle.product();

  privateData_->clearTrkVectorsGsf();

  int blockSize=0;

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsGsfTrackFiller") << "Track length " << collection->size() 
                                        << " is too long for declared max length for tree "
                                        << maxTracks_ << " and no output flag is set."
                                        << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsGsfTrackFiller") << "Track length " << collection->size() 
                                        << " is too long for declared max length for tree "
                                        << maxTracks_ 
                                        << ". Collection will be truncated ";
    }

    try { iEvent.getByLabel(vertexCollection_, primaryVertex_); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsGsfTrackFiller") << "Can't get candidate collection: " << vertexCollection_; }

    if ( saveDeDx_ ) {
      iEvent.getByLabel( "dedxTruncated40", energyLoss_ );
      iEvent.getByLabel(refittedTracksForDeDxTag_,refittedTracksForDeDx_);
      *(privateData_->ncand) = refittedTracksForDeDx_->size();   
      blockSize = (&(*refittedTracksForDeDx_)) ? refittedTracksForDeDx_->size() : 0;
    } else {
      blockSize = (collection) ? blockSize = collection->size() : 0;
      *(privateData_->ncand) = collection->size();
    }

    for(int i=0; i < (int)collection->size(); i++) {
      if( saveDeDx_ ) {
        if ( i != (int)refittedTracksForDeDx_->size() ) {
          RefToBase<reco::Track> refittedTrack(refittedTracksForDeDx_, i);
          if(saveTrk_) writeTrkInfo( refittedTrack );
          writeDeDxInfo( refittedTrack );
        }
      } else {
        RefToBase<reco::Track> trkRef(collectionHandle, i);
        if(saveTrk_) writeTrkInfo( trkRef );
      }

      GsfTrackRef gsfTrackRef = RefToBase<reco::Track>(collectionHandle, i).castTo<GsfTrackRef>();
      writeGsfTrkInfo( gsfTrackRef );
    
    }
    

  } else {
    *(privateData_->ncand) = 0;
  }

  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 

  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");

  if(saveTrk_)  treeTrkInfo(columnPrefix,columnSuffix);
  if(saveDeDx_) treeDeDxInfo(columnPrefix,columnSuffix);
  treeGsfTrkInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();
	
}

void CmsGsfTrackFiller::writeGsfTrkInfo(GsfTrackRef trkRef) {
  
  privateData_->chargeMode->push_back(trkRef->chargeMode());
  privateData_->pxMode->push_back(trkRef->pxMode());
  privateData_->pyMode->push_back(trkRef->pyMode());
  privateData_->pzMode->push_back(trkRef->pzMode());

}

void CmsGsfTrackFiller::treeGsfTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

    cmstree->column((colPrefix+"chargeMode"+colSuffix).c_str(), *privateData_->chargeMode, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pxMode"+colSuffix).c_str(), *privateData_->pxMode, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pyMode"+colSuffix).c_str(), *privateData_->pyMode, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pzMode"+colSuffix).c_str(), *privateData_->pzMode, nCandString.c_str(), 0, "Reco");

}

void CmsGsfTrackFillerData::initialiseGsf() {
  
  initialise();
  chargeMode = new vector<int>;
  pxMode = new vector<float>;
  pyMode = new vector<float>;
  pzMode = new vector<float>;
 
}

void CmsGsfTrackFillerData::clearTrkVectorsGsf() {

  clearTrkVectors();

  chargeMode->clear();
  pxMode->clear();
  pyMode->clear();
  pzMode->clear();

}
