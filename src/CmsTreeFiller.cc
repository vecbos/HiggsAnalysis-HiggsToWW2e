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
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectronFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTreeFiller.h"

#include "PhysicsTools/HepMCCandAlgos/interface/MCCandMatcher.h"
#include "CLHEP/HepMC/GenParticle.h"


#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;

struct CmsTreeFillerData {

  CmsTree *cmstree;
  bool saveTrk;
  bool saveEcal;
  bool saveHcal;
  bool saveDT;
  bool saveCSC;
  bool saveRPC;

  bool saveFatTrk;
  bool saveFatEcal;
  bool saveFatHcal;
  bool saveFatDT;
  bool saveFatCSC;
  bool saveFatRPC;

  bool saveEleID;

  bool saveCand;
  bool doMcMatch;

  bool hitLimitsMeansNoOutput;
  int maxTracks;
  int maxNeutrals;
  int maxMCTracks;

  std::string *trkIndexName;

  // All the vectors that will store the stuff
  // going into the tree.

  vector<int> *charge;
  vector<float> *energy, *et, *momentum;
  vector<float> *vertexX, *vertexY, *vertexZ;
  vector<float> *theta, *eta, *phi;
  vector<float> *x, *y, *z;
  vector<float> *mass, *mt;
  vector<int> *pdgId;

  vector<int> *mcIndex;

  vector<float> *pxAtInner, *pyAtInner, *pzAtInner, *xAtInner, *yAtInner, *zAtInner;
  vector<float> *pxAtOuter, *pyAtOuter, *pzAtOuter, *xAtOuter, *yAtOuter, *zAtOuter;

  vector<float> *ecal, *eraw, *caloEta, *caloPhi;
  vector<int> *nClu, *nCry;

  vector<float> *e2x2, *e3x3, *e5x5;
  vector<float> *eMax, *e2nd;
  vector<float> *s1s9, *s9s25;
  vector<float> *covEtaEta, *covEtaPhi, *covPhiPhi;
  vector<float> *lat, *phiLat, *etaLat, *a20, *a42;
  
  int *ncand;

};


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsTreeFiller::CmsTreeFiller(CmsTree *cmstree, int maxTracks, 
					   int maxNeutrals,int maxMCTracks,
					   bool noOutputIfLimitsReached):
  privateData_(new CmsTreeFillerData)
{
  privateData_->cmstree=cmstree;

  privateData_->saveFatTrk=true;
  privateData_->saveFatEcal=true;
  privateData_->saveFatHcal=true;
  privateData_->saveFatDT=true;
  privateData_->saveFatCSC=true;
  privateData_->saveFatRPC=true;

  privateData_->saveCand=true;
  privateData_->doMcMatch=false;

  privateData_->trkIndexName = new std::string("n");

  privateData_->hitLimitsMeansNoOutput = noOutputIfLimitsReached;
  privateData_->maxTracks=maxTracks;
  privateData_->maxNeutrals=maxNeutrals;
  privateData_->maxMCTracks=maxMCTracks;

  initialise();
}

CmsTreeFiller::CmsTreeFiller(CmsTree *cmstree, bool fatTree, int maxTracks,
					   int maxNeutrals,int maxMCTracks,
					   bool noOutputIfLimitsReached):
  privateData_(new CmsTreeFillerData)
{
  privateData_->cmstree=cmstree;

  privateData_->saveFatTrk=fatTree;
  privateData_->saveFatEcal=fatTree;
  privateData_->saveFatHcal=fatTree;
  privateData_->saveFatDT=fatTree;
  privateData_->saveFatCSC=fatTree;
  privateData_->saveFatRPC=fatTree;

  privateData_->trkIndexName = new std::string("n");

  privateData_->hitLimitsMeansNoOutput = noOutputIfLimitsReached;
  privateData_->maxTracks=maxTracks;
  privateData_->maxNeutrals=maxNeutrals;
  privateData_->maxMCTracks=maxMCTracks;

  initialise();
}


//--------------
// Destructor --
//--------------

CmsTreeFiller::~CmsTreeFiller() {

  // delete here the vector ptr's
  delete privateData_->charge;
  delete privateData_->energy;
  delete privateData_->et;
  delete privateData_->momentum;
  delete privateData_->theta;
  delete privateData_->eta;
  delete privateData_->phi;
  delete privateData_->x;
  delete privateData_->y;
  delete privateData_->z;
  delete privateData_->vertexX;
  delete privateData_->vertexY;
  delete privateData_->vertexZ;
  delete privateData_->mass;
  delete privateData_->mt;
  delete privateData_->pdgId;

  delete privateData_->mcIndex;

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

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsTreeFiller::saveTrk(bool what) { privateData_->saveTrk=what;}

void CmsTreeFiller::saveEcal(bool what) { privateData_->saveEcal=what;}

void CmsTreeFiller::saveHcal(bool what) { privateData_->saveHcal=what;}

void CmsTreeFiller::saveDT(bool what) { privateData_->saveDT=what;}

void CmsTreeFiller::saveCSC(bool what) { privateData_->saveCSC=what;}

void CmsTreeFiller::saveRPC(bool what) { privateData_->saveRPC=what;}

void CmsTreeFiller::saveFatTrk(bool what) { privateData_->saveFatTrk=what;}

void CmsTreeFiller::saveFatEcal(bool what) { privateData_->saveFatEcal=what;}

void CmsTreeFiller::saveFatHcal(bool what) { privateData_->saveFatHcal=what;}

void CmsTreeFiller::saveFatDT(bool what) { privateData_->saveFatDT=what;}

void CmsTreeFiller::saveFatCSC(bool what) { privateData_->saveFatCSC=what;}

void CmsTreeFiller::saveFatRPC(bool what) { privateData_->saveFatRPC=what;}

void CmsTreeFiller::saveCand(bool what) { privateData_->saveCand=what;}

void CmsTreeFiller::saveEleID(bool what) { privateData_->saveEleID=what;}

void CmsTreeFiller::doMcMatch(bool what) { privateData_->doMcMatch=what;}

void CmsTreeFiller::writeCollectionToTree(const CandidateCollection *collection,
			   const edm::Event& iEvent, const edm::EventSetup& iSetup,
			   const std::string &columnPrefix, const std::string &columnSuffix,
			   bool dumpData) {
  
  if(privateData_->hitLimitsMeansNoOutput && 
     (int)collection->size() > privateData_->maxTracks){
    LogError("CmsTreeFiller") << "Track length " << collection->size() 
				  << " is too long for declared max length for tree "
				  << privateData_-> maxTracks << " and no output flag is set."
				  << " No tracks written to tuple for this event ";
    return;
  }

  if((int)collection->size() > privateData_->maxTracks){
    LogError("CmsTreeFiller") << "Track length " << collection->size() 
				  << " is too long for declared max length for tree "
				  << privateData_-> maxTracks 
				  << ". Collection will be truncated ";
  }
  
  *(privateData_->ncand) = collection->size();

  LogDebug("CmsTreeFiller") << "Filling candidate vectors";
  this->clearTrkVectors();
  CandidateCollection::const_iterator cand;
  for(cand=collection->begin(); cand!=collection->end(); cand++) {
    // fill basic kinematics
    if(privateData_->saveCand) writeCandInfo(&(*cand),iEvent,iSetup);
    // fill Cluster Adapter
    SuperClusterRef sclusRef = cand->get<SuperClusterRef>();
    if(privateData_->saveEcal) writeEcalInfo(&(*cand),iEvent,iSetup,sclusRef);
    // fill (GSF) Track Adapter
    GsfTrackRef trkRef = cand->get<GsfTrackRef>();
    if(privateData_->saveTrk) writeTrkInfo(&(*cand),iEvent,iSetup,trkRef);
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  
  std::string nCandString = columnPrefix+(*privateData_->trkIndexName)+columnSuffix; 
  privateData_->cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  
  if(privateData_->saveCand) treeCandInfo(columnPrefix,columnSuffix);
  if(privateData_->saveEcal) treeEcalInfo(columnPrefix,columnSuffix);
  if(privateData_->saveTrk) treeTrkInfo(columnPrefix,columnSuffix);
  if(privateData_->saveEleID) {
    CmsEleIDTreeFiller eIDFiller(privateData_->cmstree);
    eIDFiller.setStandalone(false);
    eIDFiller.writeCollectionToTree(collection,iEvent,iSetup,columnPrefix,columnSuffix,false);
  }

  if(dumpData) privateData_->cmstree->dumpData();

}

void CmsTreeFiller::writeMcIndicesToTree(const CandidateCollection *recoCollection,
					     const edm::Event& iEvent, const edm::EventSetup& iSetup,
					     const CandidateCollection *genCollection,
					     const std::string &columnPrefix, const std::string &columnSuffix,
					     bool dumpData) {
  LogDebug("CmsTreeFiller") << "Filling the tree with references to the MC truth particles";
  writeMcMatchInfo(recoCollection, iEvent, iSetup, genCollection);
  treeMcMatchInfo(columnPrefix, columnSuffix);

}

void CmsTreeFiller::writeTrkInfo(const Candidate *cand, const edm::Event& iEvent, const edm::EventSetup& iSetup, GsfTrackRef trkRef) {
  if(&trkRef) {

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

void CmsTreeFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*privateData_->trkIndexName)+colSuffix;

  privateData_->cmstree->column((colPrefix+"pxAtOuter"+colSuffix).c_str(), *privateData_->pxAtOuter, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"pyAtOuter"+colSuffix).c_str(), *privateData_->pyAtOuter, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"pzAtOuter"+colSuffix).c_str(), *privateData_->pzAtOuter, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"xAtOuter"+colSuffix).c_str(), *privateData_->xAtOuter, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"yAtOuter"+colSuffix).c_str(), *privateData_->yAtOuter, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"zAtOuter"+colSuffix).c_str(), *privateData_->zAtOuter, nCandString.c_str(), 0, "Reco");

  if(privateData_->saveFatTrk) {
    privateData_->cmstree->column((colPrefix+"pxAtInner"+colSuffix).c_str(), *privateData_->pxAtInner, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"pyAtInner"+colSuffix).c_str(), *privateData_->pyAtInner, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"pzAtInner"+colSuffix).c_str(), *privateData_->pzAtInner, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"xAtInner"+colSuffix).c_str(), *privateData_->xAtInner, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"yAtInner"+colSuffix).c_str(), *privateData_->yAtInner, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"zAtInner"+colSuffix).c_str(), *privateData_->zAtInner, nCandString.c_str(), 0, "Reco");
  }
}

void CmsTreeFiller::writeEcalInfo(const Candidate *cand, const edm::Event& iEvent, const edm::EventSetup& iSetup, SuperClusterRef sclusRef) {
  bool hasBarrel=true;
  bool hasEndcap=true;
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
    
    if(privateData_->saveFatEcal) {
      privateData_->eraw->push_back(sclusRef->rawEnergy());
      privateData_->caloEta->push_back(sclusRef->eta());
      privateData_->caloPhi->push_back(sclusRef->phi());
    }

    // Cluster shape variables
    Handle<BasicClusterShapeAssociationCollection> barrelClShpHandle;
    try { iEvent.getByLabel("hybridSuperClusters","hybridShapeAssoc", barrelClShpHandle); }
    catch ( cms::Exception& ex ) { LogWarning("CmsTreeFiller") << "Can't get ECAL barrel Cluster Shape Collection"; }
    const reco::BasicClusterShapeAssociationCollection& barrelClShpMap = *barrelClShpHandle;

    Handle<BasicClusterShapeAssociationCollection> endcapClShpHandle;
    try { iEvent.getByLabel("islandBasicClusters","islandEndcapShapeAssoc", endcapClShpHandle); }
    catch ( cms::Exception& ex ) { LogWarning("CmsTreeFiller") << "Can't get ECAL endcap Cluster Shape Collection"; }
    const reco::BasicClusterShapeAssociationCollection& endcapClShpMap = *endcapClShpHandle;

    reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
    seedShpItr = barrelClShpMap.find(sclusRef->seed());
    if(seedShpItr==barrelClShpMap.end()) {
      LogDebug("CmsTreeFiller") << "This ECAL cluster has not barrel hit, looking for encap ones";
      hasBarrel=false;
      seedShpItr=endcapClShpMap.find(sclusRef->seed());
      if(seedShpItr==endcapClShpMap.end()) hasEndcap=false;
    }
    if(hasBarrel || hasEndcap) {
      const ClusterShapeRef& sClShape = seedShpItr->val;  

      privateData_->e3x3->push_back(sClShape->e3x3());
      privateData_->e5x5->push_back(sClShape->e5x5());
      privateData_->eMax->push_back(sClShape->eMax());
      privateData_->lat->push_back(sClShape->lat());
      privateData_->phiLat->push_back(sClShape->phiLat());
      privateData_->etaLat->push_back(sClShape->etaLat());
      if(privateData_->saveFatEcal) {
	privateData_->e2x2->push_back(sClShape->e2x2());
	privateData_->e2nd->push_back(sClShape->e2nd());
	privateData_->s1s9->push_back(sClShape->eMax()/sClShape->e3x3());
	privateData_->s9s25->push_back(sClShape->e3x3()/sClShape->e5x5());
	privateData_->covEtaEta->push_back(sClShape->covEtaEta());
	privateData_->covEtaPhi->push_back(sClShape->covEtaPhi());
	privateData_->covPhiPhi->push_back(sClShape->covPhiPhi());
 	privateData_->a20->push_back(sClShape->zernike20());
 	privateData_->a42->push_back(sClShape->zernike42());
      }
    }
    else { LogWarning("CmsTreeFiller") << "Cannot find hits in ECAL barrel or ECAL encap. Why are you requesting filling ECAL infos?";}
  }
  if(!(&sclusRef) || ((&sclusRef) && (!hasBarrel & !hasEndcap)) ) {
    privateData_->ecal->push_back(-1.);
    privateData_->nClu->push_back(-1);
    privateData_->nCry->push_back(-1);
    privateData_->e3x3->push_back(-1.);
    privateData_->e5x5->push_back(-1.);
    privateData_->eMax->push_back(-1.);
    privateData_->lat->push_back(-1.);
    privateData_->phiLat->push_back(-1.);
    privateData_->etaLat->push_back(-1.);
    if(privateData_->saveFatEcal) {
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

void CmsTreeFiller::treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*privateData_->trkIndexName)+colSuffix;
  privateData_->cmstree->column((colPrefix+"ecal"+colSuffix).c_str(), *privateData_->ecal, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"nClu"+colSuffix).c_str(), *privateData_->nClu, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"nCry"+colSuffix).c_str(), *privateData_->nCry, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"e3x3"+colSuffix).c_str(), *privateData_->e3x3, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"e5x5"+colSuffix).c_str(), *privateData_->e5x5, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eMax"+colSuffix).c_str(), *privateData_->eMax, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"lat"+colSuffix).c_str(), *privateData_->lat, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"phiLat"+colSuffix).c_str(), *privateData_->phiLat, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"etaLat"+colSuffix).c_str(), *privateData_->etaLat, nCandString.c_str(), 0, "Reco");
  
  if(privateData_->saveFatEcal) {
    privateData_->cmstree->column((colPrefix+"eraw"+colSuffix).c_str(), *privateData_->eraw, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"caloEta"+colSuffix).c_str(), *privateData_->caloEta, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"caloPhi"+colSuffix).c_str(), *privateData_->caloPhi, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"e2x2"+colSuffix).c_str(), *privateData_->e2x2, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"e2nd"+colSuffix).c_str(), *privateData_->e2nd, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"s1s9"+colSuffix).c_str(), *privateData_->s1s9, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"s9s25"+colSuffix).c_str(), *privateData_->s9s25, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"covEtaEta"+colSuffix).c_str(), *privateData_->covEtaEta, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"covEtaPhi"+colSuffix).c_str(), *privateData_->covEtaPhi, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"covPhiPhi"+colSuffix).c_str(), *privateData_->covPhiPhi, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"a20"+colSuffix).c_str(), *privateData_->a20, nCandString.c_str(), 0, "Reco");
    privateData_->cmstree->column((colPrefix+"a42"+colSuffix).c_str(), *privateData_->a42, nCandString.c_str(), 0, "Reco");
  }
}

void CmsTreeFiller::writeHcalInfo(const Candidate *cand, const edm::Event& iEvent, const edm::EventSetup& iSetup) {
}

void CmsTreeFiller::treeHcalInfo(const std::string &colPrefix, const std::string &colSuffix) {
}

void CmsTreeFiller::writeMuonInfo(const Candidate *cand, const edm::Event& iEvent, const edm::EventSetup& iSetup) {
}

void CmsTreeFiller::treeMuonInfo(const std::string &colPrefix, const std::string &colSuffix) {
}

void CmsTreeFiller::writeCandInfo(const Candidate *cand, const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  privateData_->charge->push_back((int)cand->charge());
  privateData_->energy->push_back(cand->energy());
  privateData_->et->push_back(cand->et());
  privateData_->momentum->push_back(cand->p());
  privateData_->theta->push_back(cand->theta());
  privateData_->eta->push_back(cand->eta());
  privateData_->phi->push_back(cand->phi());
  privateData_->x->push_back(cand->momentum().x());
  privateData_->y->push_back(cand->momentum().y());
  privateData_->z->push_back(cand->momentum().z());
  privateData_->vertexX->push_back(cand->vx());
  privateData_->vertexY->push_back(cand->vy());
  privateData_->vertexZ->push_back(cand->vz());
  privateData_->mass->push_back(cand->mass());
  privateData_->mt->push_back(cand->mt());
  privateData_->pdgId->push_back(cand->pdgId());

}

void CmsTreeFiller::treeCandInfo(const std::string colPrefix, const std::string colSuffix) {

  std::string nCandString = colPrefix+(*privateData_->trkIndexName)+colSuffix;
  privateData_->cmstree->column((colPrefix+"charge"+colSuffix).c_str(), *privateData_->charge, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"et"+colSuffix).c_str(), *privateData_->et, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"momentum"+colSuffix).c_str(), *privateData_->momentum, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"px"+colSuffix).c_str(), *privateData_->x, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"py"+colSuffix).c_str(), *privateData_->y, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"pz"+colSuffix).c_str(), *privateData_->z, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"vertexX"+colSuffix).c_str(), *privateData_->vertexX, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"vertexY"+colSuffix).c_str(), *privateData_->vertexY, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"vertexZ"+colSuffix).c_str(), *privateData_->vertexZ, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"mass"+colSuffix).c_str(), *privateData_->mass, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"mt"+colSuffix).c_str(), *privateData_->mt, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"pdgId"+colSuffix).c_str(), *privateData_->pdgId, nCandString.c_str(), 0, "Reco");

}

void CmsTreeFiller::writeMcMatchInfo(const CandidateCollection *recoCollection, const edm::Event& iEvent, const edm::EventSetup& iSetup,
				     const CandidateCollection *genCollection) {
  
  edm::Handle<reco::CandMatchMap> mcMatchMap;
  try { iEvent.getByLabel( matchMap_, mcMatchMap ); }
  catch( cms::Exception& ex ) { edm::LogWarning("CmsMcTruthTreeFiller") << "Can't get MC match map " << matchMap_; }
  MCCandMatcher match(*mcMatchMap);
  //  MCCandMatcher<reco::CandidateCollection> match(*mcMatchMap); // this for releases > 15X, in <=14X defined with no template argument 
  
  CandidateCollection::const_iterator recoCand;
  for(recoCand=recoCollection->begin(); recoCand!=recoCollection->end(); recoCand++) {
    CandidateRef mcRef = match(*recoCand);
    
    // find the index in the MC collection
    int indMatched=-1;
    bool matched=false;
    int idx=0;
    reco::CandidateCollection::const_iterator genCandIter;
    for(genCandIter=genCollection->begin(); genCandIter!=genCollection->end(); genCandIter++) {
      const Candidate *iCand=&(*genCandIter);
      if(&(*mcRef)==&(*iCand)) {
	indMatched=idx;
	matched=true;
	break;
      }
      idx++;
    }

    if(matched) privateData_->mcIndex->push_back(indMatched);
    else privateData_->mcIndex->push_back(-1);
  }
}

void CmsTreeFiller::treeMcMatchInfo(const std::string colPrefix, const std::string colSuffix) {
  std::string nCandString = colPrefix+(*privateData_->trkIndexName)+colSuffix;
  privateData_->cmstree->column((colPrefix+"index"+colSuffix).c_str(), *privateData_->mcIndex, nCandString.c_str(), 0, "Reco");
}

void CmsTreeFiller::initialise() {

  privateData_->charge = new vector<int>;
  privateData_->energy = new vector<float>;
  privateData_->et = new vector<float>;
  privateData_->momentum = new vector<float>;
  privateData_->theta = new vector<float>;
  privateData_->eta = new vector<float>;
  privateData_->phi = new vector<float>;
  privateData_->x = new vector<float>;
  privateData_->y = new vector<float>;
  privateData_->z = new vector<float>;
  privateData_->vertexX = new vector<float>;
  privateData_->vertexY = new vector<float>;
  privateData_->vertexZ = new vector<float>;
  privateData_->mass = new vector<float>;
  privateData_->mt = new vector<float>;
  privateData_->pdgId = new vector<int>;
  privateData_->ncand = new int;

  privateData_->mcIndex = new vector<int>;

  privateData_->pxAtInner = new vector<float>;
  privateData_->pyAtInner = new vector<float>;
  privateData_->pzAtInner = new vector<float>;
  privateData_->xAtInner = new vector<float>;
  privateData_->yAtInner = new vector<float>;
  privateData_->zAtInner = new vector<float>;
  privateData_->pxAtOuter = new vector<float>;
  privateData_->pyAtOuter = new vector<float>;
  privateData_->pzAtOuter = new vector<float>;
  privateData_->xAtOuter = new vector<float>;
  privateData_->yAtOuter = new vector<float>;
  privateData_->zAtOuter = new vector<float>;

  privateData_->ecal = new vector<float>;
  privateData_->nClu = new vector<int>;
  privateData_->nCry = new vector<int>;
  privateData_->eraw = new vector<float>;
  privateData_->caloEta = new vector<float>;
  privateData_->caloPhi = new vector<float>;
  privateData_->eMax = new vector<float>;
  privateData_->e2nd = new vector<float>;
  privateData_->s1s9 = new vector<float>;
  privateData_->s9s25 = new vector<float>;
  privateData_->e2x2 = new vector<float>;
  privateData_->e3x3 = new vector<float>;
  privateData_->e5x5 = new vector<float>;
  privateData_->covEtaEta = new vector<float>;
  privateData_->covEtaPhi = new vector<float>;
  privateData_->covPhiPhi = new vector<float>;
  privateData_->lat = new vector<float>;
  privateData_->phiLat = new vector<float>;
  privateData_->etaLat = new vector<float>;
  privateData_->a20 = new vector<float>;
  privateData_->a42 = new vector<float>;

  privateData_->saveTrk=true;
  privateData_->saveEcal=true;
  privateData_->saveHcal=true;
  privateData_->saveDT=true;
  privateData_->saveCSC=true;
  privateData_->saveRPC=true;
  
  privateData_->saveCand=true;

  privateData_->saveEleID=false;
}

void CmsTreeFiller::clearTrkVectors() {

  if(privateData_->saveTrk) {
    privateData_->pxAtOuter->clear();
    privateData_->pyAtOuter->clear();
    privateData_->pzAtOuter->clear();
    privateData_->xAtOuter->clear();
    privateData_->yAtOuter->clear();
    privateData_->zAtOuter->clear();
    if(privateData_->saveFatTrk) {
      privateData_->pxAtInner->clear();
      privateData_->pyAtInner->clear();
      privateData_->pzAtInner->clear();
      privateData_->xAtInner->clear();
      privateData_->yAtInner->clear();
      privateData_->zAtInner->clear();
    }
  }

  if(privateData_->saveEcal) {
    privateData_->ecal->clear();
    privateData_->nClu->clear();
    privateData_->nCry->clear();
    privateData_->e3x3->clear();
    privateData_->e5x5->clear();
    privateData_->eMax->clear();
    privateData_->lat->clear();
    privateData_->phiLat->clear();
    privateData_->etaLat->clear();
    if(privateData_->saveFatEcal) {
      privateData_->eraw->clear();
      privateData_->caloEta->clear();
      privateData_->caloPhi->clear();
      privateData_->e2x2->clear();
      privateData_->e2nd->clear();
      privateData_->s1s9->clear();
      privateData_->s9s25->clear();
      privateData_->covEtaEta->clear();
      privateData_->covEtaPhi->clear();
      privateData_->covPhiPhi->clear();
      privateData_->a20->clear();
      privateData_->a42->clear();
    }
  }
  
  privateData_->mcIndex->clear();
  if(privateData_->saveCand) {
    privateData_->charge->clear();
    privateData_->energy->clear();
    privateData_->et->clear();
    privateData_->momentum->clear();
    privateData_->theta->clear();
    privateData_->eta->clear();
    privateData_->phi->clear();
    privateData_->x->clear();
    privateData_->y->clear();
    privateData_->z->clear();
    privateData_->vertexX->clear();
    privateData_->vertexY->clear();
    privateData_->vertexZ->clear();
    privateData_->mass->clear();
    privateData_->mt->clear();
    privateData_->pdgId->clear();
  }
}
