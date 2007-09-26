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
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectronFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include "PhysicsTools/HepMCCandAlgos/interface/MCCandMatcher.h"
#include "CLHEP/HepMC/GenParticle.h"


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
  doMcMatch_=false;

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
  delete privateData_->mass;
  delete privateData_->mt;
  delete privateData_->pdgId;

  delete privateData_->mcIndex;
  delete privateData_->ncand;
}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out



void CmsCandidateFiller::saveCand(bool what) { saveCand_=what;}

void CmsCandidateFiller::doMcMatch(bool what) { doMcMatch_=what;}




void CmsCandidateFiller::writeCollectionToTree(const CandidateCollection *collection,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {

  privateData_->clearTrkVectorsCandidate();
  
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

    CandidateCollection::const_iterator cand;
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
  
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(collection) {
    if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  }

  if(dumpData) cmstree->dumpData();

}






void CmsCandidateFiller::writeMcIndicesToTree(const CandidateCollection *recoCollection,
					      const edm::Event& iEvent, const edm::EventSetup& iSetup,
					      const CandidateCollection *genCollection,
					      const std::string &columnPrefix, const std::string &columnSuffix,
					      bool dumpData) {

  writeMcMatchInfo(recoCollection, iEvent, iSetup, genCollection);
  treeMcMatchInfo(columnPrefix, columnSuffix);

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
  privateData_->mass->push_back(cand->mass());
  privateData_->mt->push_back(cand->mt());
  privateData_->pdgId->push_back(cand->pdgId());
}







void CmsCandidateFiller::treeCandInfo(const std::string colPrefix, const std::string colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"charge"+colSuffix).c_str(), *privateData_->charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energy"+colSuffix).c_str(), *privateData_->energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"et"+colSuffix).c_str(), *privateData_->et, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"momentum"+colSuffix).c_str(), *privateData_->momentum, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"theta"+colSuffix).c_str(), *privateData_->theta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(), *privateData_->eta, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phi"+colSuffix).c_str(), *privateData_->phi, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"px"+colSuffix).c_str(), *privateData_->x, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"py"+colSuffix).c_str(), *privateData_->y, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pz"+colSuffix).c_str(), *privateData_->z, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexX"+colSuffix).c_str(), *privateData_->vertexX, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexY"+colSuffix).c_str(), *privateData_->vertexY, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"vertexZ"+colSuffix).c_str(), *privateData_->vertexZ, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mass"+colSuffix).c_str(), *privateData_->mass, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mt"+colSuffix).c_str(), *privateData_->mt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pdgId"+colSuffix).c_str(), *privateData_->pdgId, nCandString.c_str(), 0, "Reco");

}



void CmsCandidateFiller::writeMcMatchInfo(const CandidateCollection *recoCollection, 
					  const edm::Event& iEvent, const edm::EventSetup& iSetup,
					  const CandidateCollection *genCollection) {
  
  edm::Handle<reco::CandMatchMap> mcMatchMap;
  try { iEvent.getByLabel( matchMap_, mcMatchMap ); }
  catch( cms::Exception& ex ) { edm::LogWarning("CmsMcTruthTreeFiller") << "Can't get MC match map " << matchMap_; }
  //  MCCandMatcher match(*mcMatchMap);
  MCCandMatcher<reco::CandidateCollection> match(*mcMatchMap); // this for releases > 15X, in <=14X defined with no template argument 

  if(recoCollection) {
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

      //    if(mcRef) privateData_->mcIndex->push_back(mcRef.key());
      if(matched) privateData_->mcIndex->push_back(indMatched);
      else privateData_->mcIndex->push_back(-1);
    }
  }
}




void CmsCandidateFiller::treeMcMatchInfo(const std::string colPrefix, const std::string colSuffix) {
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"index"+colSuffix).c_str(), *privateData_->mcIndex, nCandString.c_str(), 0, "Reco");
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
  mass = new vector<float>;
  mt = new vector<float>;
  pdgId = new vector<int>;
  ncand = new int;

  mcIndex = new vector<int>;

}

void CmsCandidateFillerData::clearTrkVectorsCandidate() {

  mcIndex->clear();
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
  mass->clear();
  mt->clear();
  pdgId->clear();

}
