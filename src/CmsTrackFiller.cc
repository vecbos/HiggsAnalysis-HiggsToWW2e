//---------------------------------------------------------------------------
//
// Description:
//       Package:   HiggsAnalysis/HiggsToWW2e
//       Class:     CmsTrackFiller
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

#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTrackFiller.h"
#include "DataFormats/Math/interface/Point3D.h"

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


CmsTrackFiller::CmsTrackFiller(CmsTree *cmsTree, 
			       edm::InputTag vertexCollection,
			       int maxTracks, int maxMCTracks,
			       bool noOutputIfLimitsReached):
  privateData_(new CmsTrackFillerData)
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

  privateData_->initialise();
  x0 = 0.;
  y0 = 0.;
  z0 = 0.;
}

CmsTrackFiller::CmsTrackFiller(CmsTree *cmsTree, 
			       edm::InputTag vertexCollection,
			       bool fatTree, 
			       int maxTracks, int maxMCTracks,
			       bool noOutputIfLimitsReached, bool vtxtrack):
  privateData_(new CmsTrackFillerData)
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

  privateData_->initialise();
  x0 = 0.;
  y0 = 0.;
  z0 = 0.;

}


//--------------
// Destructor --
//--------------

CmsTrackFiller::~CmsTrackFiller() { 

  // Track information
  delete privateData_->pxAtInner;
  delete privateData_->pyAtInner;
  delete privateData_->pzAtInner;
  delete privateData_->xAtInner;
  delete privateData_->yAtInner;
  delete privateData_->zAtInner;
  delete privateData_->pxAtOuter;
  delete privateData_->pyAtOuter;
  delete privateData_->pzAtOuter;
  delete privateData_->xAtOuter;
  delete privateData_->yAtOuter;
  delete privateData_->zAtOuter;
  delete privateData_->px;
  delete privateData_->py;
  delete privateData_->pz;
  delete privateData_->trackNormalizedChi2;
  delete privateData_->qualityMask;
  delete privateData_->trackDxy;
  delete privateData_->trackD0;
  delete privateData_->trackDsz;
  delete privateData_->trackDz;
  delete privateData_->trackDxyError;
  delete privateData_->trackD0Error;
  delete privateData_->trackDszError;
  delete privateData_->trackDzError;
  delete privateData_->trackValidHits;
  delete privateData_->trackLostHits;
  delete privateData_->trackVx;
  delete privateData_->trackVy;
  delete privateData_->trackVz;

  delete privateData_->truncatedDeDx;
  delete privateData_->truncatedDeDxError;
  delete privateData_->truncatedDeDxNoM;
  
  delete privateData_->isPixB1;
  delete privateData_->isPixB2;
  delete privateData_->isPixE1;
  delete privateData_->isPixE2;
  delete privateData_->numberOfValidPixelBarrelHits;
  delete privateData_->numberOfValidPixelEndcapHits;
  delete privateData_->numberOfValidStripTIBHits;
  delete privateData_->numberOfValidStripTIDHits;
  delete privateData_->numberOfValidStripTOBHits;
  delete privateData_->numberOfValidStripTECHits;

  delete privateData_->ncand;

}

//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsTrackFiller::saveTrk(bool what) { saveTrk_=what;}

void CmsTrackFiller::saveFatTrk(bool what) { saveFatTrk_=what;}

void CmsTrackFiller::saveVtxTrk(bool what) { saveVtxTrk_=what;}

void CmsTrackFiller::saveDeDx(bool what ) { saveDeDx_=what; }

void CmsTrackFiller::findPrimaryVertex(const edm::Event& iEvent) {

  edm::Handle< reco::VertexCollection>  primaryVertex  ;
  try { iEvent.getByLabel(vertexCollection_, primaryVertex); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsTrackFiller") << "Can't get candidate collection: " << vertexCollection_; }
  
  if(primaryVertex->size() < 1) { // there is no vertex in the event
    x0 = 0.;
    y0 = 0.;
    z0 = 0.;
  } else {
    float MaxSumPt = -1.;
    VertexCollection::const_iterator vMax = primaryVertex->begin();
    // calculate the vertex pT 
    for(VertexCollection::const_iterator v = primaryVertex->begin();
	v != primaryVertex->end(); ++v){
      float SumPt = 0.0;
      if((*v).tracksSize() > 0){
	std::vector<TrackBaseRef >::const_iterator t;
	for( t = (*v).tracks_begin(); t != (*v).tracks_end(); t++){
	  if((**t).charge() < -1 || (**t).charge() > 1){
	    //illegal charge
	  } else {
	    SumPt += (**t).pt();
	  }
	}
      }
      
      if(SumPt > MaxSumPt) {
	MaxSumPt = SumPt;
	vMax  = v;
      } 
    }
    x0 = vMax->x();
    y0 = vMax->y();
    z0 = vMax->z();

  }
}

void CmsTrackFiller::writeCollectionToTree(edm::InputTag collectionTag,
					   const edm::Event& iEvent, const edm::EventSetup& iSetup,
					   const std::string &columnPrefix, const std::string &columnSuffix,
					   bool dumpData) {

//   edm::Handle< edm::View<reco::Candidate> > collectionHandle;
//   try { iEvent.getByLabel(collectionTag, collectionHandle); }
//   catch ( cms::Exception& ex ) { edm::LogWarning("CmsTrackFiller") << "Can't get candidate collection: " << collectionTag; }
//   const edm::View<reco::Candidate> *collection = collectionHandle.product();

  edm::Handle< edm::View<reco::Track> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsTrackFiller") << "Can't get track collection: " << collectionTag; }
  const edm::View<reco::Track> *collection = collectionHandle.product();

  privateData_->clearTrkVectors();

  int blockSize=0;

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsTrackFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ << " and no output flag is set."
			       << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsTrackFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ 
			       << ". Collection will be truncated ";
    }

    try { iEvent.getByLabel(vertexCollection_, primaryVertex_); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsTrackFiller") << "Can't get candidate collection: " << vertexCollection_; }

    if ( saveDeDx_ ) {
      iEvent.getByLabel( "dedxTruncated40", energyLoss_ );
      iEvent.getByLabel(refittedTracksForDeDxTag_,refittedTracksForDeDx_);
      *(privateData_->ncand) = refittedTracksForDeDx_->size();   
      blockSize = (&(*refittedTracksForDeDx_)) ? refittedTracksForDeDx_->size() : 0;
    } else {
      blockSize = (collection) ? blockSize = collection->size() : 0;
      *(privateData_->ncand) = collection->size();
    }

    //    edm::View<reco::Track>::const_iterator cand;
    //    int index=0;
    //    for(cand=collection->begin(); cand!=collection->end(); cand++,index++) {
    for(int i=0; i < (int)collection->size(); i++) {
      if( saveDeDx_ ) {
        if ( i != (int)refittedTracksForDeDx_->size() ) {
          RefToBase<reco::Track> refittedTrack(refittedTracksForDeDx_, i);
          if(saveTrk_) writeTrkInfo( refittedTrack );
          writeDeDxInfo( refittedTrack );
        }
      } else {
        // if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
        //        TrackRef trkRef = cand->get<TrackRef>();
        RefToBase<reco::Track> trkRef(collectionHandle, i);
        if(saveTrk_) writeTrkInfo( trkRef );
      }
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

  if(dumpData) cmstree->dumpData();
	
}


void CmsTrackFiller::writeTrkInfo(edm::RefToBase<reco::Track> trkRef) {

  if(&trkRef) {
    
    if ( saveVtxTrk_ ) { 
      
      int iVtx = -1;
      int counter = 0;
      double weight = 0.;
      if(saveVtxTrk_) {
	if(primaryVertex_->size() >0 ) { // there is at least one vertex in the event
	  for(VertexCollection::const_iterator v = primaryVertex_->begin();
	      v != primaryVertex_->end(); ++v){
	    double tmpw = v->trackWeight(trkRef);
	    if(tmpw > weight) {
	      if(weight >0) edm::LogWarning("CmsTrackFiller") << "I found this track in two vertices!!!!!!" ;
	      weight = tmpw;
	      iVtx = counter;
	    }
	    counter++;
	  }
	}
      }

      // vertex information
      privateData_->vtxIndex->push_back(iVtx);
      privateData_->vtxWeight->push_back(weight);
      
    }

    privateData_->px->push_back(trkRef->px());
    privateData_->py->push_back(trkRef->py());
    privateData_->pz->push_back(trkRef->pz());

    if ( saveFatTrk_ ) { 

      // Inner Tracker information
      privateData_->pxAtInner->push_back(trkRef->innerMomentum().x());
      privateData_->pyAtInner->push_back(trkRef->innerMomentum().y());
      privateData_->pzAtInner->push_back(trkRef->innerMomentum().z());
      
      privateData_->xAtInner->push_back(trkRef->innerPosition().x());
      privateData_->yAtInner->push_back(trkRef->innerPosition().y());
      privateData_->zAtInner->push_back(trkRef->innerPosition().z());
      
      // Outer Tracker information
      privateData_->pxAtOuter->push_back(trkRef->outerMomentum().x());
      privateData_->pyAtOuter->push_back(trkRef->outerMomentum().y());
      privateData_->pzAtOuter->push_back(trkRef->outerMomentum().z());
      
      privateData_->xAtOuter->push_back(trkRef->outerPosition().x());
      privateData_->yAtOuter->push_back(trkRef->outerPosition().y());
      privateData_->zAtOuter->push_back(trkRef->outerPosition().z());

      privateData_->recHitsSize->push_back(trkRef->recHitsSize());

      const HitPattern trackerPattern = trkRef->hitPattern();
      privateData_->isPixB1->push_back(hasValidHitInNthPixelBarrel(1,trackerPattern));
      privateData_->isPixB2->push_back(hasValidHitInNthPixelBarrel(2,trackerPattern));
      privateData_->isPixE1->push_back(hasValidHitInNthPixelEndcap(1,trackerPattern));
      privateData_->isPixE2->push_back(hasValidHitInNthPixelEndcap(2,trackerPattern));
      privateData_->numberOfValidPixelBarrelHits->push_back(trackerPattern.numberOfValidPixelBarrelHits());
      privateData_->numberOfValidPixelEndcapHits->push_back(trackerPattern.numberOfValidPixelEndcapHits());
      privateData_->numberOfValidStripTIBHits->push_back(trackerPattern.numberOfValidStripTIBHits());
      privateData_->numberOfValidStripTIDHits->push_back(trackerPattern.numberOfValidStripTIDHits());
      privateData_->numberOfValidStripTOBHits->push_back(trackerPattern.numberOfValidStripTOBHits());
      privateData_->numberOfValidStripTECHits->push_back(trackerPattern.numberOfValidStripTECHits());

    }

    else {

      privateData_->pxAtInner->push_back( -1.0 );
      privateData_->pyAtInner->push_back( -1.0 );
      privateData_->pzAtInner->push_back( -1.0 );
      
      privateData_->xAtInner->push_back( -1.0 );
      privateData_->yAtInner->push_back( -1.0 );
      privateData_->zAtInner->push_back( -1.0 );
      
      privateData_->pxAtOuter->push_back( -1.0 );
      privateData_->pyAtOuter->push_back( -1.0 );
      privateData_->pzAtOuter->push_back( -1.0 );
      
      privateData_->xAtOuter->push_back( -1.0 );
      privateData_->yAtOuter->push_back( -1.0 );
      privateData_->zAtOuter->push_back( -1.0 );

      privateData_->recHitsSize->push_back( -1.0 );

      privateData_->isPixB1->push_back(false);
      privateData_->isPixB2->push_back(false);
      privateData_->isPixE1->push_back(false);
      privateData_->isPixE1->push_back(false);
      privateData_->numberOfValidPixelBarrelHits->push_back(-1);
      privateData_->numberOfValidPixelEndcapHits->push_back(-1);
      privateData_->numberOfValidStripTIBHits->push_back(-1);
      privateData_->numberOfValidStripTIDHits->push_back(-1);
      privateData_->numberOfValidStripTOBHits->push_back(-1);
      privateData_->numberOfValidStripTECHits->push_back(-1);

    }

    // track quality
    privateData_->charge->push_back(trkRef->charge());
    privateData_->pterr->push_back(trkRef->ptError());
    privateData_->trackValidHits->push_back(trkRef->numberOfValidHits());
    privateData_->trackLostHits ->push_back(trkRef->numberOfLostHits());
    privateData_->trackNormalizedChi2->push_back(trkRef->normalizedChi2());
    privateData_->qualityMask->push_back(trkRef->qualityMask());

    /// vtx position
    privateData_->trackVx ->push_back(trkRef->vx());
    privateData_->trackVy ->push_back(trkRef->vy());
    privateData_->trackVz ->push_back(trkRef->vz());

    // distance w.r.t. (0,0,0)
    privateData_->trackDxy->push_back(trkRef->dxy());
    privateData_->trackD0 ->push_back(trkRef->d0());
    privateData_->trackDsz->push_back(trkRef->dsz());
    privateData_->trackDz->push_back(trkRef->dz());

    privateData_->trackDxyError->push_back(trkRef->dxyError());
    privateData_->trackD0Error ->push_back(trkRef->d0Error());
    privateData_->trackDszError->push_back(trkRef->dszError());
    privateData_->trackDzError ->push_back(trkRef->dzError());

    if ( saveVtxTrk_ ) { 
      // distance w.r.t. primary vertex
      privateData_->trackDxyPV->push_back(trkRef->dxy(math::XYZPoint(x0,y0,z0)));
      //    privateData_->trackD0PV->push_back(trkRef->d0(math::XYZPoint(x0,y0,z0)));
      privateData_->trackDszPV->push_back(trkRef->dsz(math::XYZPoint(x0,y0,z0)));
      privateData_->trackDzPV->push_back(trkRef->dz(math::XYZPoint(x0,y0,z0)));
      
      //    privateData_->trackDxyErrorPV->push_back(trkRef->dxyError(math::XYZPoint(x0,y0,z0)));
      //    privateData_->trackD0ErrorPV->push_back(trkRef->d0Error(math::XYZPoint(x0,y0,z0)));
      //    privateData_->trackDszErrorPV->push_back(trkRef->dszError(math::XYZPoint(x0,y0,z0)));
      //    privateData_->trackDzErrorPV->push_back(trkRef->dzError(math::XYZPoint(x0,y0,z0)));
    } else {
      privateData_->trackDxyPV->push_back(-1.);
      //    privateData_->trackD0PV->push_back(-1.);
      privateData_->trackDszPV->push_back(-1.);
      privateData_->trackDzPV->push_back(-1.);
    }

  } else {
    
    // vertex information
    privateData_->vtxIndex->push_back(-1);
    privateData_->vtxWeight->push_back(-1.);

    // Inner Tracker information
    privateData_->pxAtInner->push_back(-1.);
    privateData_->pyAtInner->push_back(-1.);
    privateData_->pzAtInner->push_back(-1.);

    privateData_->xAtInner->push_back(-1.);
    privateData_->yAtInner->push_back(-1.);
    privateData_->zAtInner->push_back(-1.);

    // Outer Tracker information
    privateData_->pxAtOuter->push_back(-1.);
    privateData_->pyAtOuter->push_back(-1.);
    privateData_->pzAtOuter->push_back(-1.);

    privateData_->xAtOuter->push_back(-1.);
    privateData_->yAtOuter->push_back(-1.);
    privateData_->zAtOuter->push_back(-1.);

    // track quality
    privateData_->charge->push_back(-1.);
    privateData_->pterr->push_back(-1.);
    privateData_->recHitsSize->push_back(-1.);
    privateData_->trackValidHits->push_back(-1.);
    privateData_->trackLostHits ->push_back(-1.);
    privateData_->trackNormalizedChi2->push_back(-1.);
    privateData_->qualityMask->push_back(-1);

    /// vtx position
    privateData_->trackVx ->push_back(-1.);
    privateData_->trackVy ->push_back(-1.);
    privateData_->trackVz ->push_back(-1.);

    // distance w.r.t. (0,0,0)
    privateData_->trackDxy->push_back(-1.);
    privateData_->trackD0 ->push_back(-1.);
    privateData_->trackDsz->push_back(-1.);
    privateData_->trackDz ->push_back(-1.);

    privateData_->trackDxyError->push_back(-1.);
    privateData_->trackD0Error ->push_back(-1.);
    privateData_->trackDszError->push_back(-1.);
    privateData_->trackDzError ->push_back(-1.);

    // distance w.r.t. primary vertex
    privateData_->trackDxyPV->push_back(-1.);
    //    privateData_->trackD0PV->push_back(-1.);
    privateData_->trackDszPV->push_back(-1.);
    privateData_->trackDzPV->push_back(-1.);

//     privateData_->trackDxyErrorPV->push_back(-1.);
//     privateData_->trackD0ErrorPV->push_back(-1.);
//     privateData_->trackDszErrorPV->push_back(-1.);
//     privateData_->trackDzErrorPV->push_back(-1.);

  }

}

void CmsTrackFiller::writeDeDxInfo( edm::RefToBase<reco::Track> refittedTrack ) {

  const DeDxDataValueMap & dedxTruncated40Val = *energyLoss_;
  
  privateData_->truncatedDeDx->push_back( dedxTruncated40Val[refittedTrack].dEdx() );
  privateData_->truncatedDeDxError->push_back( dedxTruncated40Val[refittedTrack].dEdxError() );
  privateData_->truncatedDeDxNoM->push_back( dedxTruncated40Val[refittedTrack].numberOfMeasurements() );
  
}

void CmsTrackFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"px"+colSuffix).c_str(), *privateData_->px, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"py"+colSuffix).c_str(), *privateData_->py, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pz"+colSuffix).c_str(), *privateData_->pz, nCandString.c_str(), 0, "Reco");

  if ( saveVtxTrk_ ) {
    cmstree->column((colPrefix+"vtxIndex"+colSuffix).c_str(), *privateData_->vtxIndex, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"vtxWeight"+colSuffix).c_str(), *privateData_->vtxWeight, nCandString.c_str(), 0, "Reco");
  }

  cmstree->column((colPrefix+"charge"+colSuffix).c_str(), *privateData_->charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ptError"+colSuffix).c_str(), *privateData_->pterr, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackValidHits"+colSuffix).c_str(),  *privateData_->trackValidHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackLostHits"+colSuffix).c_str(),   *privateData_->trackLostHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackNormalizedChi2"+colSuffix).c_str(), *privateData_->trackNormalizedChi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"qualityMask"+colSuffix).c_str(), *privateData_->qualityMask, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"trackDxy"+colSuffix).c_str(), *privateData_->trackDxy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackD0"+colSuffix).c_str(),  *privateData_->trackD0, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDsz"+colSuffix).c_str(), *privateData_->trackDsz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDz"+colSuffix).c_str(),  *privateData_->trackDz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDxyError"+colSuffix).c_str(), *privateData_->trackDxyError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackD0Error"+colSuffix).c_str(),  *privateData_->trackD0Error, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDszError"+colSuffix).c_str(), *privateData_->trackDszError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDzError"+colSuffix).c_str(),  *privateData_->trackDzError, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"trackDxyPV"+colSuffix).c_str(), *privateData_->trackDxyPV, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDszPV"+colSuffix).c_str(), *privateData_->trackDszPV, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackDzPV"+colSuffix).c_str(),  *privateData_->trackDzPV, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"trackVx"+colSuffix).c_str(),  *privateData_->trackVx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackVy"+colSuffix).c_str(),  *privateData_->trackVy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackVz"+colSuffix).c_str(),  *privateData_->trackVz, nCandString.c_str(), 0, "Reco");

  if ( saveFatTrk_ ) {

    cmstree->column((colPrefix+"pxAtOuter"+colSuffix).c_str(), *privateData_->pxAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pyAtOuter"+colSuffix).c_str(), *privateData_->pyAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pzAtOuter"+colSuffix).c_str(), *privateData_->pzAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"xAtOuter"+colSuffix).c_str(), *privateData_->xAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"yAtOuter"+colSuffix).c_str(), *privateData_->yAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"zAtOuter"+colSuffix).c_str(), *privateData_->zAtOuter, nCandString.c_str(), 0, "Reco");
    
    cmstree->column((colPrefix+"pxAtInner"+colSuffix).c_str(), *privateData_->pxAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pyAtInner"+colSuffix).c_str(), *privateData_->pyAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pzAtInner"+colSuffix).c_str(), *privateData_->pzAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"xAtInner"+colSuffix).c_str(), *privateData_->xAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"yAtInner"+colSuffix).c_str(), *privateData_->yAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"zAtInner"+colSuffix).c_str(), *privateData_->zAtInner, nCandString.c_str(), 0, "Reco");
    
    cmstree->column((colPrefix+"recHitsSize"+colSuffix).c_str(), *privateData_->recHitsSize, nCandString.c_str(), 0, "Reco");

    cmstree->column((colPrefix+"isPixB1"+colSuffix).c_str(), *privateData_->isPixB1, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"isPixB2"+colSuffix).c_str(), *privateData_->isPixB2, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"isPixE1"+colSuffix).c_str(), *privateData_->isPixE1, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"isPixE2"+colSuffix).c_str(), *privateData_->isPixE2, nCandString.c_str(), 0, "Reco");

    cmstree->column((colPrefix+"numberOfValidPixelBarrelHits"+colSuffix).c_str(), *privateData_->numberOfValidPixelBarrelHits, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"numberOfValidPixelEndcapHits"+colSuffix).c_str(), *privateData_->numberOfValidPixelEndcapHits, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"numberOfValidStripTIBHits"+colSuffix).c_str(), *privateData_->numberOfValidStripTIBHits, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"numberOfValidStripTIDHits"+colSuffix).c_str(), *privateData_->numberOfValidStripTIDHits, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"numberOfValidStripTOBHits"+colSuffix).c_str(), *privateData_->numberOfValidStripTOBHits, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"numberOfValidStripTECHits"+colSuffix).c_str(), *privateData_->numberOfValidStripTECHits, nCandString.c_str(), 0, "Reco");

  }

}

void CmsTrackFiller::treeDeDxInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"truncatedDeDx"+colSuffix).c_str(),  *privateData_->truncatedDeDx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"truncatedDeDxError"+colSuffix).c_str(),  *privateData_->truncatedDeDxError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"truncatedDeDxNoM"+colSuffix).c_str(),  *privateData_->truncatedDeDxNoM, nCandString.c_str(), 0, "Reco");

}

void CmsTrackFillerData::initialise() {
  
  ncand = new int;
  vtxIndex  = new vector<int>;
  vtxWeight = new vector<float>;
  pxAtInner = new vector<float>;
  pyAtInner = new vector<float>;
  pzAtInner = new vector<float>;
  xAtInner = new vector<float>;
  yAtInner = new vector<float>;
  zAtInner = new vector<float>;
  pxAtOuter = new vector<float>;
  pyAtOuter = new vector<float>;
  pzAtOuter = new vector<float>;
  xAtOuter = new vector<float>;
  yAtOuter = new vector<float>;
  zAtOuter = new vector<float>;
  px = new vector<float>;
  py = new vector<float>;
  pz = new vector<float>;
  charge = new vector<float>;
  pterr = new vector<float>;
  recHitsSize  = new vector<float>;
  trackValidHits  = new vector<float>;
  trackLostHits  = new vector<float>;
  trackNormalizedChi2 = new vector<float>;
  qualityMask = new vector<int>;
  trackDxy = new vector<float>;
  trackD0 = new vector<float>;
  trackDsz = new vector<float>;
  trackDz = new vector<float>;
  trackDxyError = new vector<float>;
  trackD0Error = new vector<float>;
  trackDszError = new vector<float>;
  trackDzError = new vector<float>;
  trackDxyPV = new vector<float>;
  trackDszPV = new vector<float>;
  trackDzPV = new vector<float>;
  trackVx = new vector<float>;
  trackVy = new vector<float>;
  trackVz = new vector<float>;
  truncatedDeDx = new vector<float>;
  truncatedDeDxError = new vector<float>;
  truncatedDeDxNoM = new vector<float>;
  isPixB1 = new vector<bool>;
  isPixB2 = new vector<bool>;
  isPixE1 = new vector<bool>;
  isPixE2 = new vector<bool>;
  numberOfValidPixelBarrelHits = new vector<int>;
  numberOfValidPixelEndcapHits = new vector<int>;
  numberOfValidStripTIBHits = new vector<int>;
  numberOfValidStripTIDHits = new vector<int>;
  numberOfValidStripTOBHits = new vector<int>;
  numberOfValidStripTECHits = new vector<int>;

}

void CmsTrackFillerData::clearTrkVectors() {

  vtxIndex->clear();
  vtxWeight->clear();
  pxAtOuter->clear();
  pyAtOuter->clear();
  pzAtOuter->clear();
  xAtOuter->clear();
  yAtOuter->clear();
  zAtOuter->clear();
  pxAtInner->clear();
  pyAtInner->clear();
  pzAtInner->clear();
  xAtInner->clear();
  yAtInner->clear();
  zAtInner->clear();
  px->clear();
  py->clear();
  pz->clear();
  charge->clear();
  pterr->clear();
  recHitsSize->clear();
  trackValidHits->clear();
  trackLostHits->clear();
  trackNormalizedChi2->clear();
  qualityMask->clear();
  trackDxy->clear();
  trackD0 ->clear();
  trackDsz->clear();
  trackDz ->clear();
  trackDxyError->clear();
  trackD0Error ->clear();
  trackDszError->clear();
  trackDzError ->clear();
  trackValidHits->clear();
  trackLostHits ->clear();
  trackVx->clear();
  trackVy->clear();
  trackVz->clear();
  truncatedDeDx->clear();
  truncatedDeDxError->clear();
  truncatedDeDxNoM->clear();
  isPixB1->clear();
  isPixB2->clear();
  isPixE1->clear();
  isPixE2->clear();
  numberOfValidPixelBarrelHits->clear();
  numberOfValidPixelEndcapHits->clear();
  numberOfValidStripTIBHits->clear();
  numberOfValidStripTIDHits->clear();
  numberOfValidStripTOBHits->clear();
  numberOfValidStripTECHits->clear();

}

bool CmsTrackFiller::hasValidHitInNthPixelBarrel(uint32_t nlayer, HitPattern hitpattern) {
  for (int i=0; i<(PatternSize * 32) / HitSize; i++) {
    uint32_t pattern = hitpattern.getHitPattern(i);
    if (hitpattern.pixelBarrelHitFilter(pattern)) {
      if (hitpattern.getLayer(pattern) == nlayer) {
        if (hitpattern.validHitFilter(pattern)) {
          return true;
        }
      }
    }
  }
  return false;
}

bool CmsTrackFiller::hasValidHitInNthPixelEndcap(uint32_t nlayer, HitPattern hitpattern) {
  for (int i=0; i<(PatternSize * 32) / HitSize; i++) {
    uint32_t pattern = hitpattern.getHitPattern(i);
    if (hitpattern.pixelEndcapHitFilter(pattern)) {
      if (hitpattern.getLayer(pattern) == nlayer) {
        if (hitpattern.validHitFilter(pattern)) {
          return true;
        }
      }
    }
  }
  return false;
}
