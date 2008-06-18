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
#include "PhysicsTools/HepMCCandAlgos/interface/MCCandMatcher.h"

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
  
  delete privateData_->  pxECAL;
  delete privateData_->  pyECAL;
  delete privateData_->  pzECAL;
  delete privateData_->  xECAL;
  delete privateData_->  yECAL;
  delete privateData_->  zECAL;
  delete privateData_->  EexpECAL;
  
  delete privateData_->  pxHCAL;
  delete privateData_->  pyHCAL;
  delete privateData_->  pzHCAL;
  delete privateData_->  xHCAL;
  delete privateData_->  yHCAL;
  delete privateData_->  zHCAL;
  delete privateData_->  EexpHCAL;

  delete privateData_->  pxHO;
  delete privateData_->  pyHO;
  delete privateData_->  pzHO;
  delete privateData_->  xHO;
  delete privateData_->  yHO;
  delete privateData_->  zHO;
  delete privateData_->  EexpHO;
  
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

      // fill tracks extra informations
      TrackRef trkRef = cand->get<TrackRef>();
      if(saveTrk_) writeTrkInfo(&(*cand),iEvent,iSetup,trkRef);

      // fill muon extra information 
      //      const MuonRef muonRef = collection->refAt(index).castTo<MuonRef>();
      //      MuonRef muonRef = cand->get<MuonRef>();

      const reco::Muon *muon = dynamic_cast< const reco::Muon *> ( &(*cand));
      if(saveMuonExtras_) writeMuonInfo(&(*cand),iEvent,iSetup,&(*muon));

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
				     TrackRef trkRef) {
  if(&trkRef) {
    
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

    // commented out for 184
    // Create the FreeTrajectoryState from the track
//     GlobalVector vector( trkRef->momentum().x(), trkRef->momentum().y(), trkRef->momentum().z() );
//     GlobalPoint point( trkRef->vertex().x(), trkRef->vertex().y(),  trkRef->vertex().z() );
//     GlobalTrajectoryParameters tPars(point, vector, trkRef->charge(), &*bField);    
//     HepSymMatrix covT(6,1); covT *= 1e-6; // initialize to sigma=1e-3
//     CartesianTrajectoryError tCov(covT);
    
//     // Create the Helix and propagate it
//     FreeTrajectoryState* innerState = new FreeTrajectoryState(tPars, tCov);
//     SteppingHelixStateInfo trackOrigin(*innerState);
//     cachedTrajectory_.setStateAtIP(trackOrigin);

//     cachedTrajectory_.propagateAll(trackOrigin);

//     // get trajectory in calorimeters
//     cachedTrajectory_.findEcalTrajectory(ecalDetIdAssociator_.volume() );
//     cachedTrajectory_.findHcalTrajectory(hcalDetIdAssociator_.volume() );
//     cachedTrajectory_.findHOTrajectory( hoDetIdAssociator_.volume() );
    
    //    XYZPoint ECALPoint = cachedTrajectory_.getStateAtEcal().position(); 

//     // 3Momentum at ECAL inner surface
//     privateData_->xECAL->push_back(cachedTrajectory_.getStateAtEcal().position().x());
//     privateData_->yECAL->push_back(cachedTrajectory_.getStateAtEcal().position().y());
//     privateData_->zECAL->push_back(cachedTrajectory_.getStateAtEcal().position().z());
    
//     // 3D Position at ECAL inner surface
//     privateData_->pxECAL->push_back(cachedTrajectory_.getStateAtEcal().momentum().x());
//     privateData_->pyECAL->push_back(cachedTrajectory_.getStateAtEcal().momentum().y());
//     privateData_->pzECAL->push_back(cachedTrajectory_.getStateAtEcal().momentum().z());
    
//     // 3Momentum at HCAL inner surface
//     privateData_->xHCAL->push_back(cachedTrajectory_.getStateAtHcal().position().x());
//     privateData_->yHCAL->push_back(cachedTrajectory_.getStateAtHcal().position().y());
//     privateData_->zHCAL->push_back(cachedTrajectory_.getStateAtHcal().position().z());
    
//     // 3D Position at HCAL inner surface
//     privateData_->pxHCAL->push_back(cachedTrajectory_.getStateAtHcal().momentum().x());
//     privateData_->pyHCAL->push_back(cachedTrajectory_.getStateAtHcal().momentum().y());
//     privateData_->pzHCAL->push_back(cachedTrajectory_.getStateAtHcal().momentum().z());
    
//     // 3Momentum at HO inner surface
//     privateData_->xHO->push_back(cachedTrajectory_.getStateAtHO().position().x());
//     privateData_->yHO->push_back(cachedTrajectory_.getStateAtHO().position().y());
//     privateData_->zHO->push_back(cachedTrajectory_.getStateAtHO().position().z());
    
//     // 3D Position at HO inner surface
//     privateData_->pxHO->push_back(cachedTrajectory_.getStateAtHO().momentum().x());
//     privateData_->pyHO->push_back(cachedTrajectory_.getStateAtHO().momentum().y());
//     privateData_->pzHO->push_back(cachedTrajectory_.getStateAtHO().momentum().z());
    
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

    // 3Momentum at ECAL inner surface
//     privateData_->pxECAL->push_back(-1.);
//     privateData_->pyECAL->push_back(-1.);
//     privateData_->pzECAL->push_back(-1.);
    
//     // 3D Position at ECAL inner surface
//     privateData_->xECAL->push_back(-1.);
//     privateData_->yECAL->push_back(-1.);
//     privateData_->zECAL->push_back(-1.);
    
//     // 3Momentum at HCAL inner surface
//     privateData_->pxHCAL->push_back(-1.);
//     privateData_->pyHCAL->push_back(-1.);
//     privateData_->pzHCAL->push_back(-1.);
    
//     // 3D Position at HCAL inner surface
//     privateData_->xHCAL->push_back(-1.);
//     privateData_->yHCAL->push_back(-1.);
//     privateData_->zHCAL->push_back(-1.);
    
//     // 3Momentum at HO inner surface
//     privateData_->pxHO->push_back(-1.);
//     privateData_->pyHO->push_back(-1.);
//     privateData_->pzHO->push_back(-1.);
    
//     // 3D Position at HO inner surface
//     privateData_->xHO->push_back(-1.);
//     privateData_->yHO->push_back(-1.);
//     privateData_->zHO->push_back(-1.);
    
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

  // 3Momentum at ECAL inner surface
  //     cmstree->column((colPrefix+"pxECAL"+colSuffix).c_str(), *privateData_->pxECAL, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"pyECAL"+colSuffix).c_str(), *privateData_->pyECAL, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"pzECAL"+colSuffix).c_str(), *privateData_->pzECAL, nCandString.c_str(), 0, "Reco");  
    
  //     // 3D Position at ECAL inner surface
  //     cmstree->column((colPrefix+"xECAL"+colSuffix).c_str(), *privateData_->xECAL, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"yECAL"+colSuffix).c_str(), *privateData_->yECAL, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"zECAL"+colSuffix).c_str(), *privateData_->zECAL, nCandString.c_str(), 0, "Reco");
    
  //     // 3Momentum at HCAL inner surface
  //     cmstree->column((colPrefix+"pxHCAL"+colSuffix).c_str(), *privateData_->pxHCAL, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"pyHCAL"+colSuffix).c_str(), *privateData_->pyHCAL, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"pzHCAL"+colSuffix).c_str(), *privateData_->pzHCAL, nCandString.c_str(), 0, "Reco");  
    
  //     // 3D Position at HCAL inner surface
  //     cmstree->column((colPrefix+"xHCAL"+colSuffix).c_str(), *privateData_->xHCAL, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"yHCAL"+colSuffix).c_str(), *privateData_->yHCAL, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"zHCAL"+colSuffix).c_str(), *privateData_->zHCAL, nCandString.c_str(), 0, "Reco");
    
  //     // 3Momentum at HO inner surface
  //     cmstree->column((colPrefix+"pxHO"+colSuffix).c_str(), *privateData_->pxHO, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"pyHO"+colSuffix).c_str(), *privateData_->pyHO, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"pzHO"+colSuffix).c_str(), *privateData_->pzHO, nCandString.c_str(), 0, "Reco");  
    
  //     // 3D Position at HO inner surface
  //     cmstree->column((colPrefix+"xHO"+colSuffix).c_str(), *privateData_->xHO, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"yHO"+colSuffix).c_str(), *privateData_->yHO, nCandString.c_str(), 0, "Reco");
  //     cmstree->column((colPrefix+"zHO"+colSuffix).c_str(), *privateData_->zHO, nCandString.c_str(), 0, "Reco");

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

  pxECAL = new vector<float>;
  pyECAL = new vector<float>;
  pzECAL = new vector<float>;
  xECAL = new vector<float>;
  yECAL = new vector<float>;
  zECAL = new vector<float>;
  EexpECAL = new vector<float>;

  pxHCAL = new vector<float>;
  pyHCAL = new vector<float>;
  pzHCAL = new vector<float>;
  xHCAL = new vector<float>;
  yHCAL = new vector<float>;
  zHCAL = new vector<float>;
  EexpHCAL = new vector<float>;

  pxHO = new vector<float>;
  pyHO = new vector<float>;
  pzHO = new vector<float>;
  xHO = new vector<float>;
  yHO = new vector<float>;
  zHO = new vector<float>;
  EexpHO = new vector<float>;

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
 
  pxECAL->clear();
  pyECAL->clear();
  pzECAL->clear();
  xECAL->clear();
  yECAL->clear();
  zECAL->clear();
  EexpECAL->clear();

  pxHCAL->clear();
  pyHCAL->clear();
  pzHCAL->clear();
  xHCAL->clear();
  yHCAL->clear();
  zHCAL->clear();
  EexpHCAL->clear();

  pxHO->clear();
  pyHO->clear();
  pzHO->clear();
  xHO->clear();
  yHO->clear();
  zHO->clear();
  EexpHO->clear();

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


void CmsMuonFiller::SetGeometry(const edm::EventSetup& iSetup) {

  // SEE http://cmslxr.fnal.gov/lxr/source/TrackingTools/TrackAssociator/src/TrackDetectorAssociator.cc?v=CMSSW_1_6_9

  // setup propagator
//   iSetup.get<IdealMagneticFieldRecord>().get(bField);

//   SteppingHelixPropagator* prop  = new SteppingHelixPropagator(&*bField,anyDirection);
//   prop->setMaterialMode(false);
//   prop->applyRadX0Correction(true);
//   // prop->setDebug(true); // tmp
//   cachedTrajectory_.setPropagator(prop);

//   // access the calorimeter geometry
//   iSetup.get<IdealGeometryRecord>().get(theCaloGeometry_);
//   if (!theCaloGeometry_.isValid()) 
//     throw cms::Exception("FatalError") << "Unable to find IdealGeometryRecord in event!\n";
//   // get the tracking Geometry
//   iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry_);
//   if (!theTrackingGeometry_.isValid()) 
//     throw cms::Exception("FatalError") << "Unable to find GlobalTrackingGeometryRecord in event!\n";
  

//   ecalDetIdAssociator_.setGeometry(&*theCaloGeometry_);
//   caloDetIdAssociator_.setGeometry(&*theCaloGeometry_);
//   hcalDetIdAssociator_.setGeometry(&*theCaloGeometry_);
//   hoDetIdAssociator_.setGeometry(&*theCaloGeometry_);
//   muonDetIdAssociator_.setGeometry(&*theTrackingGeometry_);

  //   iSetup.get<DetIdAssociatorRecord>().get("EcalDetIdAssociator", ecalDetIdAssociator_);
  //   iSetup.get<DetIdAssociatorRecord>().get("HcalDetIdAssociator", hcalDetIdAssociator_);
  //   iSetup.get<DetIdAssociatorRecord>().get("HODetIdAssociator", hoDetIdAssociator_);
  //   iSetup.get<DetIdAssociatorRecord>().get("CaloDetIdAssociator", caloDetIdAssociator_);
  //   iSetup.get<DetIdAssociatorRecord>().get("MuonDetIdAssociator", muonDetIdAssociator_);  

  //! CHECK! changed setDetectorRadius and setDetectorLength to the following in 1_8_4 migration
  // crashes: map is not valid
//   cachedTrajectory_.setMaxDetectorRadius(muonDetIdAssociator_.volume().maxR());
//   cachedTrajectory_.setMaxDetectorLength(2.*muonDetIdAssociator_.volume().maxZ());

//     double HOmaxR = hoDetIdAssociator_.volume().maxR();
//     double HOmaxZ = hoDetIdAssociator_.volume().maxZ();
//     double minR = ecalDetIdAssociator_.volume().minR();
//     double minZ = ecalDetIdAssociator_.volume().minZ();
//     cachedTrajectory_.setMaxHORadius(HOmaxR);
//     cachedTrajectory_.setMaxHOLength(HOmaxZ*2.);
//     cachedTrajectory_.setMinDetectorRadius(minR);
//     cachedTrajectory_.setMinDetectorLength(minZ*2.);
//     double maxR = muonDetIdAssociator_.volume().maxR();
//     double maxZ = muonDetIdAssociator_.volume().maxZ();
//     cachedTrajectory_.setMaxDetectorRadius(maxR);
//     cachedTrajectory_.setMaxDetectorLength(maxZ*2.);
  
  

}
