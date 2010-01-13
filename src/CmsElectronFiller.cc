//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsTreeFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  9 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsElectronFiller.h"

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


CmsElectronFiller::CmsElectronFiller(CmsTree *cmsTree, 
				     int maxTracks, int maxMCTracks,
				     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsElectronFillerData)
{
  cmstree=cmsTree;

  saveTrk_=true;
  saveEcal_=true;
  saveFatTrk_=true;
  saveFatEcal_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsElectronFiller::CmsElectronFiller(CmsTree *cmsTree, bool fatTree, 
				     int maxTracks, int maxMCTracks,
				     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsElectronFillerData)
{
  cmstree=cmsTree;

  saveTrk_=true;
  saveEcal_=true;
  saveFatTrk_=fatTree;
  saveFatEcal_=fatTree;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();

}


//--------------
// Destructor --
//--------------

CmsElectronFiller::~CmsElectronFiller() {

  // delete here the vector ptr's
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
  delete privateData_->modePx;
  delete privateData_->modePy;
  delete privateData_->modePz;
  delete privateData_->eleTrackNormalizedChi2;
  delete privateData_->eleTrackDxy;
  delete privateData_->eleTrackD0;
  delete privateData_->eleTrackDsz;
  delete privateData_->eleTrackDz;
  delete privateData_->eleTrackDxyError;
  delete privateData_->eleTrackD0Error;
  delete privateData_->eleTrackDszError;
  delete privateData_->eleTrackDzError;
  delete privateData_->eleTrackValidHits;
  delete privateData_->eleTrackLostHits;
  delete privateData_->eleTrackVx;
  delete privateData_->eleTrackVy;
  delete privateData_->eleTrackVz;

  delete privateData_->fiducialFlags;
  delete privateData_->recoFlags;
  delete privateData_->esEnergy;
  delete privateData_->energyCorrections;

  delete privateData_->superClusterIndex;
  delete privateData_->PFsuperClusterIndex;
  delete privateData_->trackIndex;
  delete privateData_->gsfTrackIndex;

  delete privateData_->ncand;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsElectronFiller::saveTrk(bool what) { saveTrk_=what;}

void CmsElectronFiller::saveEcal(bool what) { saveEcal_=what;}

void CmsElectronFiller::saveFatTrk(bool what) { saveFatTrk_=what;}

void CmsElectronFiller::saveFatEcal(bool what) { saveFatEcal_=what;}

void CmsElectronFiller::saveEleID(bool what) { saveEleID_=what;}



void CmsElectronFiller::writeCollectionToTree(edm::InputTag collectionTag,
					      const edm::Event& iEvent, const edm::EventSetup& iSetup,
					      const std::string &columnPrefix, const std::string &columnSuffix,
					      bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get electron candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsElectronFiller") << "Track length " << collection->size() 
				   << " is too long for declared max length for tree "
				   << maxTracks_ << " and no output flag is set."
				   << " No tracks written to tuple for this event ";
      return;
    }
  
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsElectronFiller") << "Track length " << collection->size() 
				   << " is too long for declared max length for tree "
				   << maxTracks_ 
				   << ". Collection will be truncated ";
    }

    *(privateData_->ncand) = collection->size();

    // superclusters
    Handle<SuperClusterCollection> EcalBarrelSuperClusters;
    try { iEvent.getByLabel(EcalBarrelSuperClusters_, EcalBarrelSuperClusters); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL barrel supercluster Collection" << EcalBarrelSuperClusters_; }
    
    Handle<SuperClusterCollection> EcalEndcapSuperClusters;
    try { iEvent.getByLabel(EcalEndcapSuperClusters_, EcalEndcapSuperClusters); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL endcap supercluster Collection" << EcalEndcapSuperClusters_; }
    
    barrelSuperClustersSize = EcalBarrelSuperClusters->size();

    // for cluster shape variables
    Handle< EcalRecHitCollection > EcalBarrelRecHits;
    try { iEvent.getByLabel(EcalBarrelRecHits_, EcalBarrelRecHits); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL barrel rec hits Collection" << EcalBarrelRecHits_; }
    const EcalRecHitCollection *EBRecHits = EcalBarrelRecHits.product();

    Handle< EcalRecHitCollection > EcalEndcapRecHits;
    try { iEvent.getByLabel(EcalEndcapRecHits_, EcalEndcapRecHits); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL endcap rec hits Collection" << EcalEndcapRecHits_; }
    const EcalRecHitCollection *EERecHits = EcalEndcapRecHits.product();

    for(int index = 0; index < (int)collection->size(); index++) {

      // fill basic kinematics
      const Candidate *cand = &(collection->at(index));
      if(saveCand_) writeCandInfo(cand,iEvent,iSetup);

      const GsfElectronRef electronRef = collection->refAt(index).castTo<GsfElectronRef>();

      if ( !(electronRef.isNull()) ) {

        // fill Cluster Adapter
        SuperClusterRef sclusRef = electronRef->superCluster();
        SuperClusterRef pfclusRef = electronRef->pflowSuperCluster();
        if(saveEcal_) writeEcalInfo(electronRef,iEvent,iSetup,sclusRef,pfclusRef,EBRecHits,EERecHits );
        // fill (GSF) Track Adapter
        GsfTrackRef trkRef = cand->get<GsfTrackRef>();
        if(saveTrk_) writeTrkInfo(electronRef,iEvent,iSetup,trkRef);

      } else {
        edm::LogWarning("CmsElectronFiller") << "Warning! The collection seems to be not made by "
                                             << "electrons, electron-specific infos will be set to default.";
      }

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
  if(saveEcal_) treeEcalInfo(columnPrefix,columnSuffix);
  if(saveTrk_) treeTrkInfo(columnPrefix,columnSuffix);
  if(saveEleID_) {
    CmsEleIDTreeFiller eIDFiller(cmstree);
    eIDFiller.setStandalone(false);
    eIDFiller.setEcalBarrelRecHits(EcalBarrelRecHits_);
    eIDFiller.setEcalEndcapRecHits(EcalEndcapRecHits_);
    // those for egamma official isolations
    eIDFiller.setTkIsolationProducer(tkIsolationProducer_);
    eIDFiller.setTowerIsolationProducer(towerIsolationProducer_);
    // those for private H->WW isolations
    eIDFiller.setTracksProducer(tracksProducer_);
    eIDFiller.setCalotowersProducer(calotowersProducer_);
    eIDFiller.writeCollectionToTree(collectionTag,iEvent,iSetup,columnPrefix,columnSuffix,false);
  }
  
  
  if(dumpData) cmstree->dumpData();

}




void CmsElectronFiller::writeTrkInfo(const GsfElectronRef electronRef, 
				     const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				     GsfTrackRef trkRef) {
  if(&trkRef) {

    privateData_->gsfTrackIndex->push_back(trkRef.key());

    reco::TrackRef closeCtfTrack = electronRef->closestCtfTrackRef();
    privateData_->trackIndex->push_back(closeCtfTrack.key());

    privateData_->modePx->push_back(electronRef->trackMomentumAtVtx().x());
    privateData_->modePy->push_back(electronRef->trackMomentumAtVtx().y());
    privateData_->modePz->push_back(electronRef->trackMomentumAtVtx().z());

    if ( saveFatTrk_ ) {

      privateData_->eleTrackNormalizedChi2->push_back(trkRef->normalizedChi2());
      
      privateData_->eleTrackDxy->push_back(trkRef->dxy());
      privateData_->eleTrackD0 ->push_back(trkRef->d0());
      privateData_->eleTrackDsz->push_back(trkRef->dsz());
      privateData_->eleTrackDz ->push_back(trkRef->dz());
      
      privateData_->eleTrackDxyError->push_back(trkRef->dxyError());
      privateData_->eleTrackD0Error ->push_back(trkRef->d0Error());
      privateData_->eleTrackDszError->push_back(trkRef->dszError());
      privateData_->eleTrackDzError ->push_back(trkRef->dzError());
      
      privateData_->eleTrackValidHits->push_back(trkRef->numberOfValidHits());
      privateData_->eleTrackLostHits ->push_back(trkRef->numberOfLostHits());
      
      privateData_->eleTrackVx ->push_back(trkRef->vx());
      privateData_->eleTrackVy ->push_back(trkRef->vy());
      privateData_->eleTrackVz ->push_back(trkRef->vz());
      
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

      privateData_->eleTrackNormalizedChi2->push_back(-1.);

      privateData_->eleTrackDxy->push_back(-1.);
      privateData_->eleTrackD0 ->push_back(-1.);
      privateData_->eleTrackDsz->push_back(-1.);
      privateData_->eleTrackDz ->push_back(-1.);

      privateData_->eleTrackDxyError->push_back(-1.);
      privateData_->eleTrackD0Error ->push_back(-1.);
      privateData_->eleTrackDszError->push_back(-1.);
      privateData_->eleTrackDzError ->push_back(-1.);

      privateData_->eleTrackValidHits->push_back(-1.);					       
      privateData_->eleTrackLostHits ->push_back(-1.);

      privateData_->eleTrackVx ->push_back(-1.);
      privateData_->eleTrackVy ->push_back(-1.);
      privateData_->eleTrackVz ->push_back(-1.);

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
    
  }
  else {


  }
}




void CmsElectronFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"modePx"+colSuffix).c_str(),  *privateData_->modePx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"modePy"+colSuffix).c_str(),  *privateData_->modePy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"modePz"+colSuffix).c_str(),  *privateData_->modePz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(),  *privateData_->trackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsfTrackIndex"+colSuffix).c_str(),  *privateData_->gsfTrackIndex, nCandString.c_str(), 0, "Reco");

  if(saveFatTrk_) {
    cmstree->column((colPrefix+"trackNormalizedChi2"+colSuffix).c_str(), *privateData_->eleTrackNormalizedChi2, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackDxy"+colSuffix).c_str(), *privateData_->eleTrackDxy, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackD0"+colSuffix).c_str(),  *privateData_->eleTrackD0, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackDsz"+colSuffix).c_str(), *privateData_->eleTrackDsz, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackDz"+colSuffix).c_str(),  *privateData_->eleTrackDz, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackDxyError"+colSuffix).c_str(), *privateData_->eleTrackDxyError, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackD0Error"+colSuffix).c_str(),  *privateData_->eleTrackD0Error, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackDszError"+colSuffix).c_str(), *privateData_->eleTrackDszError, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackDzError"+colSuffix).c_str(),  *privateData_->eleTrackDzError, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackValidHits"+colSuffix).c_str(),  *privateData_->eleTrackValidHits, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackLostHits"+colSuffix).c_str(),   *privateData_->eleTrackLostHits, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackVx"+colSuffix).c_str(),  *privateData_->eleTrackVx, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackVy"+colSuffix).c_str(),  *privateData_->eleTrackVy, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"trackVz"+colSuffix).c_str(),  *privateData_->eleTrackVz, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pxAtInner"+colSuffix).c_str(), *privateData_->pxAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pyAtInner"+colSuffix).c_str(), *privateData_->pyAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pzAtInner"+colSuffix).c_str(), *privateData_->pzAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"xAtInner"+colSuffix).c_str(), *privateData_->xAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"yAtInner"+colSuffix).c_str(), *privateData_->yAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"zAtInner"+colSuffix).c_str(), *privateData_->zAtInner, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pxAtOuter"+colSuffix).c_str(), *privateData_->pxAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pyAtOuter"+colSuffix).c_str(), *privateData_->pyAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"pzAtOuter"+colSuffix).c_str(), *privateData_->pzAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"xAtOuter"+colSuffix).c_str(), *privateData_->xAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"yAtOuter"+colSuffix).c_str(), *privateData_->yAtOuter, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"zAtOuter"+colSuffix).c_str(), *privateData_->zAtOuter, nCandString.c_str(), 0, "Reco");
  }
}





void CmsElectronFiller::writeEcalInfo(const GsfElectronRef electronRef, 
				      const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				      SuperClusterRef sclusRef, SuperClusterRef pfclusRef,
				      const EcalRecHitCollection *EBRecHits,
				      const EcalRecHitCollection *EERecHits) {


  bool validTopologyAndGeometry = false;
  if(&electronRef) {

    // fiducial flags in ECAL
    int packed_sel = -1;
    int isEB = ( electronRef->isEB() ) ? 1 : 0;
    int isEE = ( electronRef->isEE() ) ? 1 : 0;
    int isGap = ( electronRef->isGap() ) ? 1 : 0;
    int isEBEEGap = ( electronRef->isEBEEGap() ) ? 1 : 0;
    int isEBGap = ( electronRef->isEBGap() ) ? 1 : 0;
    int isEBEtaGap = ( electronRef->isEBEtaGap() ) ? 1 : 0;
    int isEBPhiGap = ( electronRef->isEBPhiGap() ) ? 1 : 0;
    int isEEGap = ( electronRef->isEEGap() ) ? 1 : 0;
    int isEEDeeGap = ( electronRef->isEEDeeGap() ) ? 1 : 0;
    int isEERingGap = ( electronRef->isEERingGap() ) ? 1 : 0;
    
    packed_sel = ( isEB << 9 ) | ( isEE << 8 ) | ( isGap << 7 ) |
      ( isEBEEGap << 6 ) | ( isEBGap << 5 ) | ( isEBEtaGap << 4 ) | ( isEBPhiGap << 3 ) |
      ( isEEGap << 2 ) | ( isEEDeeGap << 1 ) | isEERingGap;

    privateData_->fiducialFlags->push_back(packed_sel);

    int packed_reco;
    int isEcalDriven = ( electronRef->isEcalDriven() ) ? 1 : 0;
    int isTrackerDriven = ( electronRef->isTrackerDriven() ) ? 1 : 0;
    packed_reco = ( isEcalDriven << 1 ) | isTrackerDriven;
    privateData_->recoFlags->push_back( packed_reco );

    // link to the supercluster (collections are merged: barrel + endcap in this order)
    if ( isEcalDriven && sclusRef.isNonnull() ) {
      int offset = ( fabs(sclusRef->eta() ) < 1.479 ) ? 0 : barrelSuperClustersSize;
      privateData_->superClusterIndex->push_back( sclusRef.key() + offset );
    } else {
      privateData_->superClusterIndex->push_back( -1 );
    }

    if ( isTrackerDriven && pfclusRef.isNonnull() ) {
      privateData_->PFsuperClusterIndex->push_back( pfclusRef.key() );
    } else {
      privateData_->PFsuperClusterIndex->push_back( -1 );
    }
    
    int packed_corr;
    int isEcalEnergyCorrected = ( electronRef->isEcalEnergyCorrected() ) ? 1 : 0;
    int isMomentumCorrected = ( electronRef->isMomentumCorrected() ) ? 1 : 0;
    packed_corr = ( isEcalEnergyCorrected << 1 ) | isMomentumCorrected;
    privateData_->energyCorrections->push_back( packed_corr );

    // preshower energy
    privateData_->esEnergy->push_back(sclusRef->preshowerEnergy());

  } else {
    privateData_->fiducialFlags->push_back(-1);
    privateData_->recoFlags->push_back(-1);
    privateData_->superClusterIndex->push_back( -1 );
    privateData_->PFsuperClusterIndex->push_back( -1 );
    privateData_->energyCorrections->push_back( -1 );
    privateData_->esEnergy->push_back(-1.);
  }

}

void CmsElectronFiller::treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"fiducialFlags"+colSuffix).c_str(), *privateData_->fiducialFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"recoFlags"+colSuffix).c_str(), *privateData_->recoFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energyCorrections"+colSuffix).c_str(), *privateData_->energyCorrections, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esEnergy"+colSuffix).c_str(), *privateData_->esEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"superClusterIndex"+colSuffix).c_str(), *privateData_->superClusterIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFsuperClusterIndex"+colSuffix).c_str(), *privateData_->PFsuperClusterIndex, nCandString.c_str(), 0, "Reco");

}




void CmsElectronFillerData::initialise() {
  
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
  modePx = new vector<float>;
  modePy = new vector<float>;
  modePz = new vector<float>;
  eleTrackNormalizedChi2 = new vector<float>;
  eleTrackDxy = new vector<float>;
  eleTrackD0 = new vector<float>;
  eleTrackDsz = new vector<float>;
  eleTrackDz = new vector<float>;
  eleTrackDxyError = new vector<float>;
  eleTrackD0Error = new vector<float>;
  eleTrackDszError = new vector<float>;
  eleTrackDzError = new vector<float>;
  eleTrackValidHits = new vector<float>;
  eleTrackLostHits = new vector<float>;
  eleTrackVx = new vector<float>;
  eleTrackVy = new vector<float>;
  eleTrackVz = new vector<float>;

  fiducialFlags = new vector<int>;
  recoFlags = new vector<int>;
  esEnergy = new vector<float>;
  energyCorrections = new vector<int>;

  superClusterIndex = new vector<int>;
  PFsuperClusterIndex = new vector<int>;
  trackIndex = new vector<int>;
  gsfTrackIndex = new vector<int>;

}

void CmsElectronFillerData::clearTrkVectors() {

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
  modePx->clear();
  modePy->clear();
  modePz->clear();  
  eleTrackNormalizedChi2->clear();
  eleTrackDxy->clear();
  eleTrackD0 ->clear();
  eleTrackDsz->clear();
  eleTrackDz ->clear();
  eleTrackDxyError->clear();
  eleTrackD0Error ->clear();
  eleTrackDszError->clear();
  eleTrackDzError ->clear();
  eleTrackValidHits->clear();
  eleTrackLostHits ->clear();
  eleTrackVx->clear();
  eleTrackVy->clear();
  eleTrackVz->clear();

  fiducialFlags->clear();
  recoFlags->clear();
  esEnergy->clear();
  energyCorrections->clear();

  superClusterIndex->clear();
  PFsuperClusterIndex->clear();
  trackIndex->clear();
  gsfTrackIndex->clear();

}
