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
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMuonFiller.h"

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


CmsMuonFiller::CmsMuonFiller(CmsTree *cmsTree, 
			     int maxTracks, int maxMCTracks,
			     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsMuonFillerData)
{
  cmstree=cmsTree;

  saveMuonExtras_=true;
  saveTrk_=true;
  saveFatTrk_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();

  
}

CmsMuonFiller::CmsMuonFiller(CmsTree *cmsTree, 
			     bool fatTree, 
			     int maxTracks, int maxMCTracks,
			     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsMuonFillerData)
{
  cmstree=cmsTree;

  saveMuonExtras_=true;
  saveTrk_=true;
  saveFatTrk_=fatTree;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsMuonFiller::~CmsMuonFiller() { 
  
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
  delete privateData_->muTrackNormalizedChi2;
  delete privateData_->muTrackDxy;
  delete privateData_->muTrackD0;
  delete privateData_->muTrackDsz;
  delete privateData_->muTrackDz;
  delete privateData_->muTrackDxyError;
  delete privateData_->muTrackD0Error;
  delete privateData_->muTrackDszError;
  delete privateData_->muTrackDzError;
  delete privateData_->muTrackValidHits;
  delete privateData_->muTrackLostHits;
  delete privateData_->muTrackVx;
  delete privateData_->muTrackVy;
  delete privateData_->muTrackVz;

  delete privateData_->isGlobal;
  delete privateData_->isTracker;
  delete privateData_->isStandAlone;
  delete privateData_->isCalo;
  
  delete privateData_->  sumPt03;
  delete privateData_->  emEt03;
  delete privateData_->  hadEt03;
  delete privateData_->  hoEt03;
  delete privateData_->  nTrk03;
  delete privateData_->  nJets03;

  delete privateData_->  sumPt05;
  delete privateData_->  emEt05;
  delete privateData_->  hadEt05;
  delete privateData_->  hoEt05;
  delete privateData_->  nTrk05;
  delete privateData_->  nJets05;

  delete privateData_->EcalExpDepo;
  delete privateData_->HcalExpDepo;
  delete privateData_->HoExpDepo;
  delete privateData_->emS9;
  delete privateData_->hadS9;
  delete privateData_->hoS9;
  delete privateData_->CaloComp;

  
  delete privateData_->ncand;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsMuonFiller::saveMuonExtras(bool what) { saveMuonExtras_=what; }

void CmsMuonFiller::saveTrk(bool what) { saveTrk_=what;}

void CmsMuonFiller::saveFatTrk(bool what) { saveFatTrk_=what;}

void CmsMuonFiller::writeCollectionToTree(edm::InputTag collectionTag,
					 const edm::Event& iEvent, const edm::EventSetup& iSetup,
					 const std::string &columnPrefix, const std::string &columnSuffix,
					 bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMuonFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsMuonFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ << " and no output flag is set."
			       << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsMuonFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ 
			       << ". Collection will be truncated ";
    }
  
    *(privateData_->ncand) = collection->size();

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);

      // fill muon extra information
      const reco::Muon *muon = dynamic_cast< const reco::Muon *> ( &(*cand));
      if(saveMuonExtras_) writeMuonInfo(&(*cand),iEvent,iSetup,&(*muon));

      // fill tracks extra informations (only if the muon has a tracker track)
      if(saveTrk_) writeTrkInfo(&(*cand),iEvent,iSetup,&(*muon));

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
  if(saveTrk_) treeTrkInfo(columnPrefix,columnSuffix);
  if(saveMuonExtras_) treeMuonInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();
	
}


void CmsMuonFiller::writeTrkInfo(const Candidate *cand, 
				 const edm::Event& iEvent, const edm::EventSetup& iSetup,
				 const Muon *muon) {

  
  TrackRef trkRef;
  bool hasTrackerTrack = false;

  if( & muon ) {  
    if ( muon->track().isNonnull() ) {
      hasTrackerTrack = true;
      trkRef = cand->get<TrackRef>();
    }
  }

  if( hasTrackerTrack && &trkRef!=0 ) {

    if ( saveFatTrk_ ) { 
      
      privateData_->pxAtInner->push_back(trkRef->innerMomentum().x());
      privateData_->pyAtInner->push_back(trkRef->innerMomentum().y());
      privateData_->pzAtInner->push_back(trkRef->innerMomentum().z());
      
      privateData_->xAtInner->push_back(trkRef->innerPosition().x());
      privateData_->yAtInner->push_back(trkRef->innerPosition().y());
      privateData_->zAtInner->push_back(trkRef->innerPosition().z());
      
      privateData_->pxAtOuter->push_back(trkRef->outerMomentum().x());
      privateData_->pyAtOuter->push_back(trkRef->outerMomentum().y());
      privateData_->pzAtOuter->push_back(trkRef->outerMomentum().z());
      
      privateData_->xAtOuter->push_back(trkRef->outerPosition().x());
      privateData_->yAtOuter->push_back(trkRef->outerPosition().y());
      privateData_->zAtOuter->push_back(trkRef->outerPosition().z());
      
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

    }

    privateData_->muTrackNormalizedChi2->push_back(trkRef->normalizedChi2());
    
    privateData_->muTrackDxy->push_back(trkRef->dxy());
    privateData_->muTrackD0 ->push_back(trkRef->d0());
    privateData_->muTrackDsz->push_back(trkRef->dsz());
    privateData_->muTrackDz ->push_back(trkRef->dz());

    privateData_->muTrackDxyError->push_back(trkRef->dxyError());
    privateData_->muTrackD0Error ->push_back(trkRef->d0Error());
    privateData_->muTrackDszError->push_back(trkRef->dszError());
    privateData_->muTrackDzError ->push_back(trkRef->dzError());

    privateData_->muTrackValidHits->push_back(trkRef->numberOfValidHits());
    privateData_->muTrackLostHits ->push_back(trkRef->numberOfLostHits());

    privateData_->muTrackVx ->push_back(trkRef->vx());
    privateData_->muTrackVy ->push_back(trkRef->vy());
    privateData_->muTrackVz ->push_back(trkRef->vz());

  }

  else {
    privateData_->pxAtInner->push_back(-1.);
    privateData_->pyAtInner->push_back(-1.);
    privateData_->pzAtInner->push_back(-1.);

    privateData_->xAtInner->push_back(-1.);
    privateData_->yAtInner->push_back(-1.);
    privateData_->zAtInner->push_back(-1.);

    privateData_->pxAtOuter->push_back(-1.);
    privateData_->pyAtOuter->push_back(-1.);
    privateData_->pzAtOuter->push_back(-1.);

    privateData_->xAtOuter->push_back(-1.);
    privateData_->yAtOuter->push_back(-1.);
    privateData_->zAtOuter->push_back(-1.);

    privateData_->muTrackNormalizedChi2->push_back(-1.);

    privateData_->muTrackDxy->push_back(-1.);
    privateData_->muTrackD0 ->push_back(-1.);
    privateData_->muTrackDsz->push_back(-1.);
    privateData_->muTrackDz ->push_back(-1.);

    privateData_->muTrackDxyError->push_back(-1.);
    privateData_->muTrackD0Error ->push_back(-1.);
    privateData_->muTrackDszError->push_back(-1.);
    privateData_->muTrackDzError ->push_back(-1.);

    privateData_->muTrackValidHits->push_back(-1.);					       
    privateData_->muTrackLostHits ->push_back(-1.);

    privateData_->muTrackVx ->push_back(-1.);
    privateData_->muTrackVy ->push_back(-1.);
    privateData_->muTrackVz ->push_back(-1.);

  }
}

void CmsMuonFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  
  cmstree->column((colPrefix+"muTrackNormalizedChi2"+colSuffix).c_str(), *privateData_->muTrackNormalizedChi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackDxy"+colSuffix).c_str(), *privateData_->muTrackDxy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackD0"+colSuffix).c_str(),  *privateData_->muTrackD0, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackDsz"+colSuffix).c_str(), *privateData_->muTrackDsz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackDz"+colSuffix).c_str(),  *privateData_->muTrackDz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackDxyError"+colSuffix).c_str(), *privateData_->muTrackDxyError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackD0Error"+colSuffix).c_str(),  *privateData_->muTrackD0Error, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackDszError"+colSuffix).c_str(), *privateData_->muTrackDszError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackDzError"+colSuffix).c_str(),  *privateData_->muTrackDzError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackValidHits"+colSuffix).c_str(),  *privateData_->muTrackValidHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackLostHits"+colSuffix).c_str(),   *privateData_->muTrackLostHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackVx"+colSuffix).c_str(),  *privateData_->muTrackVx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackVy"+colSuffix).c_str(),  *privateData_->muTrackVy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"muTrackVz"+colSuffix).c_str(),  *privateData_->muTrackVz, nCandString.c_str(), 0, "Reco");

  if(saveFatTrk_) {

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

  }
}

void CmsMuonFiller::writeMuonInfo(const Candidate *cand, const edm::Event& iEvent, 
				  const edm::EventSetup& iSetup, const Muon *muon) {
  if(&muon) {

    privateData_->isGlobal->push_back(muon->isGlobalMuon());
    privateData_->isTracker->push_back(muon->isTrackerMuon());
    privateData_->isStandAlone->push_back(muon->isStandAloneMuon());
    privateData_->isCalo->push_back(muon->isCaloMuon());

    // default isolation variables 0.3
    MuonIsolation Iso03  = muon->isolationR03();
    privateData_->sumPt03->push_back(Iso03.sumPt);
    privateData_->emEt03->push_back(Iso03.emEt);
    privateData_->hadEt03->push_back(Iso03.hadEt);
    privateData_->hoEt03->push_back(Iso03.hoEt);
    privateData_->nTrk03->push_back(Iso03.nTracks);
    privateData_->nJets03->push_back(Iso03.nJets);

    // default isolation variables 0.5
    MuonIsolation Iso05  = muon->isolationR05();
    privateData_->sumPt05->push_back(Iso05.sumPt);
    privateData_->emEt05->push_back(Iso05.emEt);
    privateData_->hadEt05->push_back(Iso05.hadEt);
    privateData_->hoEt05->push_back(Iso05.hoEt);
    privateData_->nTrk05->push_back(Iso05.nTracks);
    privateData_->nJets05->push_back(Iso05.nJets);

    // Expected deposits in CALO
    privateData_->EcalExpDepo->push_back(muon->calEnergy().em);
    privateData_->HcalExpDepo->push_back(muon->calEnergy().had);
    privateData_->HoExpDepo->push_back(muon->calEnergy().ho);
    privateData_->emS9->push_back(muon->calEnergy().emS9);
    privateData_->hadS9->push_back(muon->calEnergy().hadS9);
    privateData_->hoS9->push_back(muon->calEnergy().hoS9);
    privateData_->CaloComp->push_back(muon->caloCompatibility());

  } else {

    // default isolation variables 0.3
    privateData_->sumPt03->push_back(-1.);
    privateData_->emEt03->push_back(-1.);
    privateData_->hadEt03->push_back(-1.);
    privateData_->hoEt03->push_back(-1.);
    privateData_->nTrk03->push_back(-1.);
    privateData_->nJets03->push_back(-1.);
    // default isolation variables 0.5
     privateData_->sumPt05->push_back(-1.);
    privateData_->emEt05->push_back(-1.);
    privateData_->hadEt05->push_back(-1.);
    privateData_->hoEt05->push_back(-1.);
    privateData_->nTrk05->push_back(-1.);
    privateData_->nJets05->push_back(-1.);

    // Expected deposits in CALO
    privateData_->EcalExpDepo->push_back(-1.);
    privateData_->HcalExpDepo->push_back(-1.);
    privateData_->HoExpDepo->push_back(-1.);
    privateData_->emS9->push_back(-1.);
    privateData_->hadS9->push_back(-1.);
    privateData_->hoS9->push_back(-1.);
    privateData_->CaloComp->push_back(-1.);

  }
}

void CmsMuonFiller::treeMuonInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
   
  // type of muon
  cmstree->column((colPrefix+"isGlobal"+colSuffix).c_str(), *privateData_->isGlobal, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"isTracker"+colSuffix).c_str(), *privateData_->isTracker, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"isStandAlone"+colSuffix).c_str(), *privateData_->isStandAlone, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"isCalo"+colSuffix).c_str(), *privateData_->isCalo, nCandString.c_str(), 0, "Reco");

  // isolation R=0.3
  cmstree->column((colPrefix+"sumPt03"+colSuffix).c_str(), *privateData_->sumPt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emEt03"+colSuffix).c_str(), *privateData_->emEt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadEt03"+colSuffix).c_str(), *privateData_->hadEt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hoEt03"+colSuffix).c_str(), *privateData_->hoEt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nTrk03"+colSuffix).c_str(), *privateData_->nTrk03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nJets03"+colSuffix).c_str(), *privateData_->nJets03, nCandString.c_str(), 0, "Reco");

  // isolation R=0.5
  cmstree->column((colPrefix+"sumPt05"+colSuffix).c_str(), *privateData_->sumPt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emEt05"+colSuffix).c_str(), *privateData_->emEt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadEt05"+colSuffix).c_str(), *privateData_->hadEt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hoEt05"+colSuffix).c_str(), *privateData_->hoEt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nTrk05"+colSuffix).c_str(), *privateData_->nTrk05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nJets05"+colSuffix).c_str(), *privateData_->nJets05, nCandString.c_str(), 0, "Reco");

  //  Expected deposits in CALO
  cmstree->column((colPrefix+"EcalExpDepo"+colSuffix).c_str(), *privateData_->EcalExpDepo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HcalExpDepo"+colSuffix).c_str(), *privateData_->HcalExpDepo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HoExpDepo"+colSuffix).c_str(), *privateData_->HoExpDepo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emS9"+colSuffix).c_str(), *privateData_->emS9, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadS9"+colSuffix).c_str(), *privateData_->hadS9, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hoS9"+colSuffix).c_str(), *privateData_->hoS9, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"CaloComp"+colSuffix).c_str(), *privateData_->CaloComp, nCandString.c_str(), 0, "Reco");

}

void CmsMuonFillerData::initialise() {

  initialiseCandidate();
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
  muTrackNormalizedChi2 = new vector<float>;
  muTrackDxy = new vector<float>;
  muTrackD0 = new vector<float>;
  muTrackDsz = new vector<float>;
  muTrackDz = new vector<float>;
  muTrackDxyError = new vector<float>;
  muTrackD0Error = new vector<float>;
  muTrackDszError = new vector<float>;
  muTrackDzError = new vector<float>;
  muTrackValidHits = new vector<float>;
  muTrackLostHits = new vector<float>;
  muTrackVx = new vector<float>;
  muTrackVy = new vector<float>;
  muTrackVz = new vector<float>;

  isGlobal = new vector<int>;
  isTracker = new vector<int>;
  isStandAlone = new vector<int>;
  isCalo = new vector<int>;
  
  sumPt03 = new vector<float>;
  emEt03 = new vector<float>;
  hadEt03 = new vector<float>;
  hoEt03 = new vector<float>;
  nTrk03 = new vector<float>;
  nJets03 = new vector<float>;
  sumPt05 = new vector<float>;
  emEt05 = new vector<float>;
  hadEt05 = new vector<float>;
  hoEt05 = new vector<float>;
  nTrk05 = new vector<float>;
  nJets05 = new vector<float>;

  EcalExpDepo = new vector<float>;
  HcalExpDepo = new vector<float>;
  HoExpDepo = new vector<float>;
  emS9 = new vector<float>;
  hadS9 = new vector<float>;
  hoS9 = new vector<float>;
  CaloComp = new vector<float>;
  
}

void CmsMuonFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

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

  muTrackNormalizedChi2->clear();
  muTrackDxy->clear();
  muTrackD0 ->clear();
  muTrackDsz->clear();
  muTrackDz ->clear();
  muTrackDxyError->clear();
  muTrackD0Error ->clear();
  muTrackDszError->clear();
  muTrackDzError ->clear();
  muTrackValidHits->clear();
  muTrackLostHits ->clear();
  muTrackVx->clear();
  muTrackVy->clear();
  muTrackVz->clear();
 
  isGlobal->clear();
  isTracker->clear();
  isStandAlone->clear();
  isCalo->clear();

  sumPt03->clear();
  emEt03->clear();
  hadEt03->clear();
  hoEt03->clear();
  nTrk03->clear();
  nJets03->clear();

  sumPt05->clear();
  emEt05->clear();
  hadEt05->clear();
  hoEt05->clear();
  nTrk05->clear();
  nJets05->clear();

  EcalExpDepo->clear();
  HcalExpDepo->clear();
  HoExpDepo->clear();
  emS9->clear();
  hadS9->clear();
  hoS9->clear();
  CaloComp->clear();
  
}
