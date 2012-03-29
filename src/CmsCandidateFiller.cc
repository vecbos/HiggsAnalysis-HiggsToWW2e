//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsCandidateFiller
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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


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


CmsCandidateFiller::CmsCandidateFiller(CmsTree *cmsTree, int maxTracks, 
				       int maxMCTracks,
				       bool noOutputIfLimitsReached):
  privateData_(new CmsCandidateFillerData)
{
  cmstree=cmsTree;

  saveCand_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialiseCandidate();
}

CmsCandidateFiller::CmsCandidateFiller(CmsTree *cmsTree, bool fatTree, int maxTracks,
				       int maxMCTracks,
				       bool noOutputIfLimitsReached):
  privateData_(new CmsCandidateFillerData)
{
  cmstree=cmsTree;

  saveCand_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_ = maxTracks;
  maxMCTracks_ = maxMCTracks;

  privateData_->initialiseCandidate();
}


//--------------
// Destructor --
//--------------

CmsCandidateFiller::~CmsCandidateFiller() {

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
  delete privateData_->ncand;
  delete privateData_;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out



void CmsCandidateFiller::saveCand(bool what) { saveCand_=what;}


void CmsCandidateFiller::writeCollectionToTree(edm::InputTag collectionTag,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsCandidateFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  privateData_->clearTrkVectorsCandidate();

  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  LogInfo("CmsCandidateFiller") << "=== Writing collection " << nCandString << " ===";

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      LogError("CmsCandidateFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      LogError("CmsCandidateFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ 
				     << ". Collection will be truncated ";
    }
  
    *(privateData_->ncand) = collection->size();

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
    }
  }
  else {
    *(privateData_->ncand) = 0;
  }

  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 

  int blockSize = (collection) ? collection->size() : 0;
  
  //  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");

  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();

  delete trkIndexName_;

}


void CmsCandidateFiller::writeCandInfo(const Candidate *cand, 
				       const edm::Event& iEvent, 
				       const edm::EventSetup& iSetup) {

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
}







void CmsCandidateFiller::treeCandInfo(const std::string colPrefix, const std::string colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"charge"+colSuffix).c_str(), *privateData_->charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"px"+colSuffix).c_str(), *privateData_->x, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"py"+colSuffix).c_str(), *privateData_->y, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pz"+colSuffix).c_str(), *privateData_->z, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexX"+colSuffix).c_str(), *privateData_->vertexX, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexY"+colSuffix).c_str(), *privateData_->vertexY, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexZ"+colSuffix).c_str(), *privateData_->vertexZ, nCandString.c_str(), 0, "Reco");
}


void CmsCandidateFillerData::initialiseCandidate() {

  charge = new vector<int>;
  energy = new vector<float>;
  et = new vector<float>;
  momentum = new vector<float>;
  theta = new vector<float>;
  eta = new vector<float>;
  phi = new vector<float>;
  x = new vector<float>;
  y = new vector<float>;
  z = new vector<float>;
  vertexX = new vector<float>;
  vertexY = new vector<float>;
  vertexZ = new vector<float>;
  ncand = new int;
}

void CmsCandidateFillerData::clearTrkVectorsCandidate() {
  charge->clear();
  energy->clear();
  et->clear();
  momentum->clear();
  theta->clear();
  eta->clear();
  phi->clear();
  x->clear();
  y->clear();
  z->clear();
  vertexX->clear();
  vertexY->clear();
  vertexZ->clear();
}
