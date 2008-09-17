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
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

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

  delete privateData_->ecal;
  delete privateData_->eraw;
  delete privateData_->caloEta;
  delete privateData_->caloPhi;
  delete privateData_->nClu;
  delete privateData_->nCry;
  delete privateData_->e2x2;
  delete privateData_->e3x3;
  delete privateData_->e5x5;
  delete privateData_->eMax;
  delete privateData_->e2nd;
  delete privateData_->s1s9;
  delete privateData_->s9s25;
  delete privateData_->covEtaEta;
  delete privateData_->covEtaPhi;
  delete privateData_->covPhiPhi;
  delete privateData_->lat;
  delete privateData_->phiLat;
  delete privateData_->etaLat;
  delete privateData_->a20;
  delete privateData_->a42;

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

    // for cluster shape variables
    Handle< EcalRecHitCollection > EcalBarrelRecHits;
    try { iEvent.getByLabel(EcalBarrelRecHits_, EcalBarrelRecHits); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL barrel rec hits Collection" << EcalBarrelRecHits_; }
    const EcalRecHitCollection *EBRecHits = EcalBarrelRecHits.product();

    Handle< EcalRecHitCollection > EcalEndcapRecHits;
    try { iEvent.getByLabel(EcalEndcapRecHits_, EcalEndcapRecHits); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL endcap rec hits Collection" << EcalEndcapRecHits_; }
    const EcalRecHitCollection *EERecHits = EcalEndcapRecHits.product();

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
      // fill Cluster Adapter
      SuperClusterRef sclusRef = cand->get<SuperClusterRef>();
      if(saveEcal_) writeEcalInfo(&(*cand),iEvent,iSetup,sclusRef,EBRecHits,EERecHits );
      // fill (GSF) Track Adapter
      GsfTrackRef trkRef = cand->get<GsfTrackRef>();
      if(saveTrk_) writeTrkInfo(&(*cand),iEvent,iSetup,trkRef);
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
    eIDFiller.setElectronIdCutsLabel(electronIdCutsLabel_);
    eIDFiller.setElectronIdLikelihoodLabel(electronIdLikelihoodLabel_);
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




void CmsElectronFiller::writeTrkInfo(const Candidate *cand, 
				     const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				     GsfTrackRef trkRef) {
  if(&trkRef) {

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

  }
}




void CmsElectronFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"eleTrackNormalizedChi2"+colSuffix).c_str(), *privateData_->eleTrackNormalizedChi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackDxy"+colSuffix).c_str(), *privateData_->eleTrackDxy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackD0"+colSuffix).c_str(),  *privateData_->eleTrackD0, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackDsz"+colSuffix).c_str(), *privateData_->eleTrackDsz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackDz"+colSuffix).c_str(),  *privateData_->eleTrackDz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackDxyError"+colSuffix).c_str(), *privateData_->eleTrackDxyError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackD0Error"+colSuffix).c_str(),  *privateData_->eleTrackD0Error, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackDszError"+colSuffix).c_str(), *privateData_->eleTrackDszError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackDzError"+colSuffix).c_str(),  *privateData_->eleTrackDzError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackValidHits"+colSuffix).c_str(),  *privateData_->eleTrackValidHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackLostHits"+colSuffix).c_str(),   *privateData_->eleTrackLostHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackVx"+colSuffix).c_str(),  *privateData_->eleTrackVx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackVy"+colSuffix).c_str(),  *privateData_->eleTrackVy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackVz"+colSuffix).c_str(),  *privateData_->eleTrackVz, nCandString.c_str(), 0, "Reco");

  if(saveFatTrk_) {
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





void CmsElectronFiller::writeEcalInfo(const Candidate *cand, 
				      const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				      SuperClusterRef sclusRef,
				      const EcalRecHitCollection *EBRecHits,
				      const EcalRecHitCollection *EERecHits) {


  bool validTopologyAndGeometry = false;
  if(&sclusRef) {
    // Cluster related variables
    privateData_->ecal->push_back(sclusRef->energy());
    privateData_->nClu->push_back(sclusRef->clustersSize());

    int ncry=0;
    reco::basicCluster_iterator bcItr;
    for(bcItr=sclusRef->clustersBegin(); bcItr!=sclusRef->clustersEnd(); bcItr++){
      reco::BasicClusterRef bclusRef = *bcItr;
      ncry+=bclusRef->getHitsByDetId().size();
    }
    privateData_->nCry->push_back(ncry);
    
    if(saveFatEcal_) {
      privateData_->eraw->push_back(sclusRef->rawEnergy());
      privateData_->caloEta->push_back(sclusRef->eta());
      privateData_->caloPhi->push_back(sclusRef->phi());
    }

    edm::ESHandle<CaloTopology> pTopology;
    iSetup.get<CaloTopologyRecord>().get(pTopology);

    edm::ESHandle<CaloGeometry> pGeometry;
    iSetup.get<CaloGeometryRecord>().get(pGeometry);

    if ( pTopology.isValid() && pGeometry.isValid() ) {

      validTopologyAndGeometry = true;

      const CaloTopology *topology = pTopology.product();
      const CaloGeometry *geometry = pGeometry.product();

      BasicClusterRef theSeed = sclusRef->seed();

      const EcalRecHitCollection *rechits = 0;

      float seedEta = theSeed->position().eta();

      if( fabs(seedEta) < 1.479 ) rechits = EBRecHits;
      else rechits = EERecHits; 

      float eMax = EcalClusterTools::eMax( *theSeed, &(*rechits) );
      float e3x3 = EcalClusterTools::e3x3( *theSeed, &(*rechits), topology );
      float e5x5 = EcalClusterTools::e5x5( *theSeed, &(*rechits), topology );

      privateData_->e3x3->push_back(e3x3);
      privateData_->e5x5->push_back(e5x5);
      privateData_->eMax->push_back(eMax);

      std::vector<float> vLat = EcalClusterTools::lat( *theSeed, &(*rechits), geometry );
      float etaLat = vLat[0];
      float phiLat = vLat[1];
      float lat = vLat[2];

      privateData_->lat->push_back(lat);
      privateData_->phiLat->push_back(phiLat);
      privateData_->etaLat->push_back(etaLat);

      if(saveFatEcal_) {

	float e2x2 = EcalClusterTools::e2x2( *theSeed, &(*rechits), topology );
	float e2nd = EcalClusterTools::e2nd( *theSeed, &(*rechits) );
	float s1s9 = eMax/e3x3;
	float s9s25 = e3x3/e5x5;
	std::vector<float> vCov = EcalClusterTools::covariances( *theSeed, &(*rechits), topology, geometry );

	float covEtaEta = vCov[0];
	float covEtaPhi = vCov[1];
	float covPhiPhi = vCov[2];

	privateData_->e2x2->push_back(e2x2);
	privateData_->e2nd->push_back(e2nd);
	privateData_->s1s9->push_back(s1s9);
	privateData_->s9s25->push_back(s9s25);
	privateData_->covEtaEta->push_back(covEtaEta);
	privateData_->covEtaPhi->push_back(covEtaPhi);
	privateData_->covPhiPhi->push_back(covPhiPhi);

	float zernike20 = EcalClusterTools::zernike20( *theSeed, &(*rechits), geometry );
	float zernike42 = EcalClusterTools::zernike42( *theSeed, &(*rechits), geometry );

 	privateData_->a20->push_back(zernike20);
 	privateData_->a42->push_back(zernike42);

      }
    }
    else { edm::LogWarning("CmsElectronFiller") << "ECAL topology or geometry not valid, not filling the ECAL cluster shapes"; }
  }
  if(!(&sclusRef) || ((&sclusRef) && (!validTopologyAndGeometry)) ) {
    privateData_->ecal->push_back(-1.);
    privateData_->nClu->push_back(-1);
    privateData_->nCry->push_back(-1);
    privateData_->e3x3->push_back(-1.);
    privateData_->e5x5->push_back(-1.);
    privateData_->eMax->push_back(-1.);
    privateData_->lat->push_back(-1.);
    privateData_->phiLat->push_back(-1.);
    privateData_->etaLat->push_back(-1.);
    if(saveFatEcal_) {
      privateData_->eraw->push_back(-1.);
      privateData_->caloEta->push_back(-1.);
      privateData_->caloPhi->push_back(-1.);
      privateData_->e2x2->push_back(-1.);
      privateData_->e2nd->push_back(-1.);
      privateData_->s1s9->push_back(-1.);
      privateData_->s9s25->push_back(-1.);
      privateData_->covEtaEta->push_back(-1.);
      privateData_->covEtaPhi->push_back(-1.);
      privateData_->covPhiPhi->push_back(-1.);
      privateData_->a20->push_back(-1.);
      privateData_->a42->push_back(-1.);
    }
  }
}





void CmsElectronFiller::treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"ecal"+colSuffix).c_str(), *privateData_->ecal, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nClu"+colSuffix).c_str(), *privateData_->nClu, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nCry"+colSuffix).c_str(), *privateData_->nCry, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e3x3"+colSuffix).c_str(), *privateData_->e3x3, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"e5x5"+colSuffix).c_str(), *privateData_->e5x5, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eMax"+colSuffix).c_str(), *privateData_->eMax, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"lat"+colSuffix).c_str(), *privateData_->lat, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phiLat"+colSuffix).c_str(), *privateData_->phiLat, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"etaLat"+colSuffix).c_str(), *privateData_->etaLat, nCandString.c_str(), 0, "Reco");
  
  if(saveFatEcal_) {
    cmstree->column((colPrefix+"eraw"+colSuffix).c_str(), *privateData_->eraw, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"caloEta"+colSuffix).c_str(), *privateData_->caloEta, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"caloPhi"+colSuffix).c_str(), *privateData_->caloPhi, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"e2x2"+colSuffix).c_str(), *privateData_->e2x2, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"e2nd"+colSuffix).c_str(), *privateData_->e2nd, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"s1s9"+colSuffix).c_str(), *privateData_->s1s9, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"s9s25"+colSuffix).c_str(), *privateData_->s9s25, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"covEtaEta"+colSuffix).c_str(), *privateData_->covEtaEta, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"covEtaPhi"+colSuffix).c_str(), *privateData_->covEtaPhi, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"covPhiPhi"+colSuffix).c_str(), *privateData_->covPhiPhi, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"a20"+colSuffix).c_str(), *privateData_->a20, nCandString.c_str(), 0, "Reco");
    cmstree->column((colPrefix+"a42"+colSuffix).c_str(), *privateData_->a42, nCandString.c_str(), 0, "Reco");
  }
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

  ecal = new vector<float>;
  nClu = new vector<int>;
  nCry = new vector<int>;
  eraw = new vector<float>;
  caloEta = new vector<float>;
  caloPhi = new vector<float>;
  eMax = new vector<float>;
  e2nd = new vector<float>;
  s1s9 = new vector<float>;
  s9s25 = new vector<float>;
  e2x2 = new vector<float>;
  e3x3 = new vector<float>;
  e5x5 = new vector<float>;
  covEtaEta = new vector<float>;
  covEtaPhi = new vector<float>;
  covPhiPhi = new vector<float>;
  lat = new vector<float>;
  phiLat = new vector<float>;
  etaLat = new vector<float>;
  a20 = new vector<float>;
  a42 = new vector<float>;

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

  ecal->clear();
  nClu->clear();
  nCry->clear();
  e3x3->clear();
  e5x5->clear();
  eMax->clear();
  lat->clear();
  phiLat->clear();
  etaLat->clear();
  eraw->clear();
  caloEta->clear();
  caloPhi->clear();
  e2x2->clear();
  e2nd->clear();
  s1s9->clear();
  s9s25->clear();
  covEtaEta->clear();
  covEtaPhi->clear();
  covPhiPhi->clear();
  a20->clear();
  a42->clear();

}
