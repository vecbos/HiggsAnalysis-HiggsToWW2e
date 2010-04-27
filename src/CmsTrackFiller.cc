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

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTrackFiller.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"

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
  isGsf_=false;
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
  delete privateData_->trackLostHits;
  delete privateData_->trackVx;
  delete privateData_->trackVy;
  delete privateData_->trackVz;

  delete privateData_->truncatedDeDx;
  delete privateData_->truncatedDeDxError;
  delete privateData_->truncatedDeDxNoM;
  delete privateData_->medianDeDx;
  delete privateData_->medianDeDxError;
  delete privateData_->medianDeDxNoM;
  delete privateData_->harmonic2DeDx;
  delete privateData_->harmonic2DeDxError;
  delete privateData_->harmonic2DeDxNoM;

  delete privateData_->pixelHits;
  delete privateData_->expInnerLayers;
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

void CmsTrackFiller::isGsf(bool what) { isGsf_=what;}

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
    bestPrimaryVertex_ = *vMax;
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

  try { iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder_); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsTrackFiller") << "Can't get TransientTrackBuilder from Event Setup."; }

  privateData_->clearTrkVectors();

  edm::ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get( magfield );     
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);       
  
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
      iEvent.getByLabel( "dedxTruncated40", truncatedEnergyLoss_ );
      iEvent.getByLabel( "dedxMedian", medianEnergyLoss_ );
      iEvent.getByLabel( "dedxHarmonic2", harmonic2EnergyLoss_ );
      iEvent.getByLabel(refittedTracksForDeDxTag_,refittedTracksForDeDx_);
      *(privateData_->ncand) = refittedTracksForDeDx_->size();   
      blockSize = (&(*refittedTracksForDeDx_)) ? refittedTracksForDeDx_->size() : 0;
    } else {
      blockSize = (collection) ? collection->size() : 0;
      *(privateData_->ncand) = collection->size();
    }

    //    edm::View<reco::Track>::const_iterator cand;
    //    int index=0;
    //    for(cand=collection->begin(); cand!=collection->end(); cand++,index++) {
    for(int i=0; i < (int)collection->size(); i++) {
      if( saveDeDx_ ) {
        if ( i != (int)refittedTracksForDeDx_->size() ) {
          RefToBase<reco::Track> refittedTrack(refittedTracksForDeDx_, i);
          if(saveTrk_) writeTrkInfo( refittedTrack , magfield, theTrackingGeometry);
          writeDeDxInfo( refittedTrack );
        }
      } else {
        // if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
        //        TrackRef trkRef = cand->get<TrackRef>();
        RefToBase<reco::Track> trkRef(collectionHandle, i);
        if(saveTrk_) writeTrkInfo( trkRef , magfield, theTrackingGeometry );
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


void CmsTrackFiller::writeTrkInfo(edm::RefToBase<reco::Track> trkRef, const edm::ESHandle<MagneticField>& magfield, const edm::ESHandle<GlobalTrackingGeometry>& theTrackingGeometry) 
{

  if(trkRef.isNonnull()) {
    
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
    
    // conversions rejection
    const HitPattern trackerPattern = trkRef->hitPattern();

    int isPixB1 = (hasValidHitInNthPixelBarrel(1,trackerPattern)) ? 1 : 0;
    int isPixB2 = (hasValidHitInNthPixelBarrel(2,trackerPattern)) ? 1 : 0;
    int isPixE1 = (hasValidHitInNthPixelEndcap(1,trackerPattern)) ? 1 : 0;
    int isPixE2 = (hasValidHitInNthPixelEndcap(2,trackerPattern)) ? 1 : 0;

    int packed_sel = (isPixB1 << 3) | (isPixB2 << 2) | (isPixE1 << 1) | isPixE2;

    privateData_->pixelHits->push_back(packed_sel);

    privateData_->numberOfValidPixelBarrelHits->push_back(trackerPattern.numberOfValidPixelBarrelHits());
    privateData_->numberOfValidPixelEndcapHits->push_back(trackerPattern.numberOfValidPixelEndcapHits());
    privateData_->numberOfValidStripTIBHits->push_back(trackerPattern.numberOfValidStripTIBHits());
    privateData_->numberOfValidStripTIDHits->push_back(trackerPattern.numberOfValidStripTIDHits());
    privateData_->numberOfValidStripTOBHits->push_back(trackerPattern.numberOfValidStripTOBHits());
    privateData_->numberOfValidStripTECHits->push_back(trackerPattern.numberOfValidStripTECHits());
    
    const HitPattern innerPattern = trkRef->trackerExpectedHitsInner();
    int expInnerLayers = innerPattern.numberOfHits();
    privateData_->expInnerLayers->push_back(expInnerLayers);

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

    if ( saveVtxTrk_ && trkRef.isNonnull() ) 
      { 

	GlobalVector direction(trkRef->px(), trkRef->py(), trkRef->pz());		    

	TransientTrack tt;
	if (isGsf_)
	  {
	    GsfTrackRef gsfTrackRef = trkRef.castTo<GsfTrackRef>();
	    if (gsfTrackRef.isNonnull())
		tt = TransientTrack(new GsfTransientTrack(gsfTrackRef,&(*magfield),&(*theTrackingGeometry)));
	  }
	else
	  {
	    tt = TransientTrack(trkRef.castTo<TrackRef>(),&(*magfield),&(*theTrackingGeometry));
	  }
	
	std::pair<bool,Measurement1D> sgnImpPar3D = std::pair<bool,Measurement1D>(false,Measurement1D());
	
	try {
	  if (bestPrimaryVertex_.isValid() && tt.isValid() )
	    sgnImpPar3D = IPTools::signedImpactParameter3D(tt,direction,bestPrimaryVertex_);
	}
	catch ( cms::Exception& ex ) {
	}
	
	if( sgnImpPar3D.first ) {
	  privateData_->impactPar3D->push_back(sgnImpPar3D.second.value());
	  privateData_->impactPar3DError->push_back(sgnImpPar3D.second.error());
	}
	else {
	  privateData_->impactPar3D->push_back(-1.);
	  privateData_->impactPar3DError->push_back(-1.);
	}
	
	std::pair<bool,Measurement1D> sgnTransvImpPar = std::pair<bool,Measurement1D>(false,Measurement1D());
	
	try {
	  if (bestPrimaryVertex_.isValid() && tt.isValid())
	    sgnTransvImpPar  = IPTools::signedTransverseImpactParameter(tt,direction,bestPrimaryVertex_);
	}
	catch ( cms::Exception& ex ) {
	}
	
	if( sgnTransvImpPar.first ) {
	  privateData_->transvImpactPar->push_back(sgnTransvImpPar.second.value());
	  privateData_->transvImpactParError->push_back(sgnTransvImpPar.second.error());
	} else {
	  privateData_->transvImpactPar->push_back(-1.);
	  privateData_->transvImpactParError->push_back(-1.);
	}
	
      } else {
      privateData_->impactPar3D->push_back(-1.);
      privateData_->impactPar3DError->push_back(-1.);
      privateData_->transvImpactPar->push_back(-1.);
      privateData_->transvImpactParError->push_back(-1.);
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

    // distance w.r.t. primary vertex
    privateData_->impactPar3D->push_back(-1.);
    privateData_->impactPar3DError->push_back(-1.);
    privateData_->transvImpactPar->push_back(-1.);
    privateData_->transvImpactParError->push_back(-1.);
  }

}

void CmsTrackFiller::writeDeDxInfo( edm::RefToBase<reco::Track> refittedTrack ) {

  const DeDxDataValueMap & dedxTruncated40Val = *truncatedEnergyLoss_;
  
  privateData_->truncatedDeDx->push_back( dedxTruncated40Val[refittedTrack].dEdx() );
  privateData_->truncatedDeDxError->push_back( dedxTruncated40Val[refittedTrack].dEdxError() );
  privateData_->truncatedDeDxNoM->push_back( dedxTruncated40Val[refittedTrack].numberOfMeasurements() );

  const DeDxDataValueMap & dedxMedianVal = *medianEnergyLoss_;
  
  privateData_->medianDeDx->push_back( dedxMedianVal[refittedTrack].dEdx() );
  privateData_->medianDeDxError->push_back( dedxMedianVal[refittedTrack].dEdxError() );
  privateData_->medianDeDxNoM->push_back( dedxMedianVal[refittedTrack].numberOfMeasurements() );
  
  const DeDxDataValueMap & dedxHarmonic2Val = *harmonic2EnergyLoss_;
  
  privateData_->harmonic2DeDx->push_back( dedxHarmonic2Val[refittedTrack].dEdx() );
  privateData_->harmonic2DeDxError->push_back( dedxHarmonic2Val[refittedTrack].dEdxError() );
  privateData_->harmonic2DeDxNoM->push_back( dedxHarmonic2Val[refittedTrack].numberOfMeasurements() );

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

  cmstree->column((colPrefix+"impactPar3D"+colSuffix).c_str(), *privateData_->impactPar3D, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"impactPar3DError"+colSuffix).c_str(), *privateData_->impactPar3DError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"transvImpactPar"+colSuffix).c_str(), *privateData_->transvImpactPar, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"transvImpactParError"+colSuffix).c_str(), *privateData_->transvImpactParError, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"trackVx"+colSuffix).c_str(),  *privateData_->trackVx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackVy"+colSuffix).c_str(),  *privateData_->trackVy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackVz"+colSuffix).c_str(),  *privateData_->trackVz, nCandString.c_str(), 0, "Reco");

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

  cmstree->column((colPrefix+"pixelHits"+colSuffix).c_str(), *privateData_->pixelHits, nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"expInnerLayers"+colSuffix).c_str(), *privateData_->expInnerLayers, nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"numberOfValidPixelBarrelHits"+colSuffix).c_str(), *privateData_->numberOfValidPixelBarrelHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"numberOfValidPixelEndcapHits"+colSuffix).c_str(), *privateData_->numberOfValidPixelEndcapHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"numberOfValidStripTIBHits"+colSuffix).c_str(), *privateData_->numberOfValidStripTIBHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"numberOfValidStripTIDHits"+colSuffix).c_str(), *privateData_->numberOfValidStripTIDHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"numberOfValidStripTOBHits"+colSuffix).c_str(), *privateData_->numberOfValidStripTOBHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"numberOfValidStripTECHits"+colSuffix).c_str(), *privateData_->numberOfValidStripTECHits, nCandString.c_str(), 0, "Reco");

}

void CmsTrackFiller::treeDeDxInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"truncatedDeDx"+colSuffix).c_str(),  *privateData_->truncatedDeDx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"truncatedDeDxError"+colSuffix).c_str(),  *privateData_->truncatedDeDxError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"truncatedDeDxNoM"+colSuffix).c_str(),  *privateData_->truncatedDeDxNoM, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"medianDeDx"+colSuffix).c_str(),  *privateData_->medianDeDx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"medianDeDxError"+colSuffix).c_str(),  *privateData_->medianDeDxError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"medianDeDxNoM"+colSuffix).c_str(),  *privateData_->medianDeDxNoM, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"harmonic2DeDx"+colSuffix).c_str(),  *privateData_->harmonic2DeDx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"harmonic2DeDxError"+colSuffix).c_str(),  *privateData_->harmonic2DeDxError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"harmonic2DeDxNoM"+colSuffix).c_str(),  *privateData_->harmonic2DeDxNoM, nCandString.c_str(), 0, "Reco");

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
  impactPar3D = new vector<float>;
  impactPar3DError = new vector<float>;
  transvImpactPar = new vector<float>;
  transvImpactParError = new vector<float>;
  trackVx = new vector<float>;
  trackVy = new vector<float>;
  trackVz = new vector<float>;
  truncatedDeDx = new vector<float>;
  truncatedDeDxError = new vector<float>;
  truncatedDeDxNoM = new vector<float>;
  medianDeDx = new vector<float>;
  medianDeDxError = new vector<float>;
  medianDeDxNoM = new vector<float>;
  harmonic2DeDx = new vector<float>;
  harmonic2DeDxError = new vector<float>;
  harmonic2DeDxNoM = new vector<float>;
  pixelHits = new vector<int>;
  expInnerLayers = new vector<int>;
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
  trackValidHits->clear();
  trackLostHits ->clear();
  impactPar3D->clear();
  impactPar3DError->clear();
  transvImpactPar->clear();
  transvImpactParError->clear();
  trackVx->clear();
  trackVy->clear();
  trackVz->clear();
  truncatedDeDx->clear();
  truncatedDeDxError->clear();
  truncatedDeDxNoM->clear();
  medianDeDx->clear();
  medianDeDxError->clear();
  medianDeDxNoM->clear();
  harmonic2DeDx->clear();
  harmonic2DeDxError->clear();
  harmonic2DeDxNoM->clear();
  pixelHits->clear();
  expInnerLayers->clear();
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
