//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsPFlowPreIdFiller
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

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFlowPreIdFiller.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TTree.h>
#include <TVector3.h>
#include <string>

using namespace edm;
using namespace reco;


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsPFlowPreIdFiller::CmsPFlowPreIdFiller(CmsTree *cmsTree, 
					 int maxTracks, int maxMCTracks,
					 bool noOutputIfLimitsReached):
  privateData_(new CmsPFlowPreIdFillerData)
{
  cmstree=cmsTree;

  trkIndexName_ = new std::string("n");  

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsPFlowPreIdFiller::CmsPFlowPreIdFiller(CmsTree *cmsTree, 
					 bool fatTree, int maxTracks, int maxMCTracks,					 
					 bool noOutputIfLimitsReached):
  privateData_(new CmsPFlowPreIdFillerData)
{
  cmstree=cmsTree;
  
  trkIndexName_ = new std::string("n");
  
  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;
  
  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsPFlowPreIdFiller::~CmsPFlowPreIdFiller() {
  
  // delete here the vector ptr's
  delete privateData_->deltaEtaMatch;
  delete privateData_->deltaPhiMatch;
  delete privateData_->chiEtaMatch;
  delete privateData_->chiPhiMatch;
  delete privateData_->chi2Match;
  delete privateData_->eopMatch;
  delete privateData_->pt;
  delete privateData_->eta;
  delete privateData_->kfChi2;
  delete privateData_->kfNHits;
  delete privateData_->ecalPosX;
  delete privateData_->ecalPosY;
  delete privateData_->ecalPosZ;
  delete privateData_->meanShowerX;
  delete privateData_->meanShowerY;
  delete privateData_->meanShowerZ;
  delete privateData_->gsfChi2Ratio;
  delete privateData_->kfNewChi2;
  delete privateData_->mva;
  delete privateData_->ecalMatching;
  delete privateData_->psMatching;
  delete privateData_->trackFiltered;
  delete privateData_->preided;
  delete privateData_->trackIndex;
  delete privateData_->ncand;
}


//-------------
// Methods   --
//-------------

void CmsPFlowPreIdFiller::writeCollectionToTree(edm::InputTag collectionTag,
						const edm::Event& iEvent, const edm::EventSetup& iSetup,
						const std::string &columnPrefix, const std::string &columnSuffix,
						bool dumpData) {
  
  edm::Handle<PreIdCollection> collection;
  try { iEvent.getByLabel(collectionTag, collection); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFlowPreIdFiller") << "Can't get preId collection: " << collectionTag; }
  
  privateData_->clearTrkVectors();
  
  int blockSize=0;

  if(&(*collection)) {

    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFlowPreIdFiller") << "Track length " << collection->size() 
					  << " is too long for declared max length for tree "
					  << maxTracks_ << " and no output flag is set."
					  << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFlowPreIdFiller") << "Track length " << collection->size() 
					  << " is too long for declared max length for tree "
					  << maxTracks_ 
					  << ". Collection will be truncated ";
    }
    
    blockSize = collection->size();
    *(privateData_->ncand) = collection->size();
  
    for(PreIdCollection::const_iterator cand=collection->begin(); cand!=collection->end(); cand++) {
      
      const reco::TrackRef theTrackRef = cand->trackRef();    


      // chiara check
      // TVector3 p3ChiaraTrack(theTrackRef->px(),theTrackRef->py(),theTrackRef->pz()); 
      // float chiaraTrackPt  = p3ChiaraTrack.Pt();
      // float chiaraTrackEta = p3ChiaraTrack.Eta();
      // std::cout << "TEST: " 
      // 	<< chiaraTrackPt  << " " << cand->pt()
      //	<< chiaraTrackEta << " " << cand->eta()
      //	<< std::endl;
			     

      
      privateData_->deltaEtaMatch->push_back( cand->deltaEtaMatch() );
      privateData_->deltaPhiMatch->push_back( cand->deltaPhiMatch() );
      privateData_->chiEtaMatch->push_back( cand->chiEtaMatch() );
      privateData_->chiPhiMatch->push_back( cand->chiPhiMatch() );
      privateData_->chi2Match->push_back( cand->chi2Match() );
      privateData_->eopMatch->push_back( cand->eopMatch() );
      privateData_->pt->push_back( cand->pt() );
      privateData_->eta->push_back( cand->eta() );
      privateData_->kfChi2->push_back( cand->kfChi2() );
      privateData_->kfNHits->push_back( cand->kfNHits() );
      privateData_->ecalPosX->push_back( cand->ecalPos().X() );
      privateData_->ecalPosY->push_back( cand->ecalPos().Y() );
      privateData_->ecalPosZ->push_back( cand->ecalPos().Z() );
      privateData_->meanShowerX->push_back( cand->meanShower().X() );
      privateData_->meanShowerY->push_back( cand->meanShower().Y() );
      privateData_->meanShowerZ->push_back( cand->meanShower().Z() );
      privateData_->gsfChi2Ratio->push_back( cand->gsfChi2Ratio() );
      privateData_->kfNewChi2->push_back( cand->kfNewChi2() );
      privateData_->mva->push_back( cand->mva() );
      privateData_->ecalMatching->push_back( cand->ecalMatching() );
      privateData_->psMatching->push_back( cand->psMatching() );
      privateData_->trackFiltered->push_back( cand->trackFiltered() );
      privateData_->preided->push_back( cand->preided() );
      privateData_->trackIndex->push_back(theTrackRef.key());
    }
  }
  else {
    *(privateData_->ncand) = 0;
    blockSize = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the tree
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  treePFPreIdInfo(columnPrefix,columnSuffix);
  
  if(dumpData) cmstree->dumpData();	
}

void CmsPFlowPreIdFiller::treePFPreIdInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"deltaEtaMatch"+colSuffix).c_str(), *privateData_->deltaEtaMatch, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiMatch"+colSuffix).c_str(), *privateData_->deltaPhiMatch, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chiEtaMatch"+colSuffix).c_str(),   *privateData_->chiEtaMatch,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chiPhiMatch"+colSuffix).c_str(),   *privateData_->chiPhiMatch,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chi2Match"+colSuffix).c_str(),     *privateData_->chi2Match,     nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eopMatch"+colSuffix).c_str(),      *privateData_->eopMatch,      nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pt"+colSuffix).c_str(),            *privateData_->pt,            nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"eta"+colSuffix).c_str(),           *privateData_->eta,           nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"kfChi2"+colSuffix).c_str(),        *privateData_->kfChi2,        nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"kfNHits"+colSuffix).c_str(),       *privateData_->kfNHits,       nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"ecalPosX"+colSuffix).c_str(),      *privateData_->ecalPosX,      nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"ecalPosY"+colSuffix).c_str(),      *privateData_->ecalPosY,      nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"ecalPosZ"+colSuffix).c_str(),      *privateData_->ecalPosZ,      nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"meanShowerX"+colSuffix).c_str(),   *privateData_->meanShowerX,   nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"meanShowerY"+colSuffix).c_str(),   *privateData_->meanShowerY,   nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"meanShowerZ"+colSuffix).c_str(),   *privateData_->meanShowerZ,   nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"gsfChi2Ratio"+colSuffix).c_str(),  *privateData_->gsfChi2Ratio,  nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"kfNewChi2"+colSuffix).c_str(),     *privateData_->kfNewChi2,     nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"mva"+colSuffix).c_str(),           *privateData_->mva,           nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"ecalMatching"+colSuffix).c_str(),  *privateData_->ecalMatching,  nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"psMatching"+colSuffix).c_str(),    *privateData_->psMatching,    nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"trackFiltered"+colSuffix).c_str(), *privateData_->trackFiltered, nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"preided"+colSuffix).c_str(),       *privateData_->preided,       nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(),    *privateData_->trackIndex,    nCandString.c_str(), 0, "Reco");  
}

void CmsPFlowPreIdFillerData::initialise() {
  
  ncand = new int;
  trackIndex = new vector<int>;
  deltaEtaMatch = new vector<float>;
  deltaPhiMatch = new vector<float>;
  chiEtaMatch = new vector<float>;
  chiPhiMatch = new vector<float>;
  chi2Match = new vector<float>;
  eopMatch = new vector<float>;
  pt = new vector<float>;
  eta = new vector<float>;
  kfChi2 = new vector<float>;
  kfNHits = new vector<float>;
  ecalPosX = new vector<float>; 
  ecalPosY = new vector<float>; 
  ecalPosZ = new vector<float>; 
  meanShowerX = new vector<float>; 
  meanShowerY = new vector<float>; 
  meanShowerZ = new vector<float>; 
  gsfChi2Ratio = new vector<float>; 
  kfNewChi2 = new vector<float>; 
  mva = new vector<float>;
  ecalMatching = new vector<bool>; 
  psMatching = new vector<bool>; 
  trackFiltered = new vector<bool>;
  preided = new vector<bool>;  
}

void CmsPFlowPreIdFillerData::clearTrkVectors() {

  deltaEtaMatch -> clear();
  deltaPhiMatch -> clear();
  chiEtaMatch   -> clear();
  chiPhiMatch   -> clear();
  chi2Match     -> clear();
  eopMatch      -> clear();
  pt            -> clear();
  eta           -> clear();
  kfChi2        -> clear();
  kfNHits       -> clear();
  ecalPosX      -> clear();
  ecalPosY      -> clear();
  ecalPosZ      -> clear();
  meanShowerX   -> clear();
  meanShowerY   -> clear();
  meanShowerZ   -> clear(); 
  ecalMatching  -> clear();
  psMatching    -> clear();
  trackFiltered -> clear();
  preided       -> clear();
  gsfChi2Ratio  -> clear();
  kfNewChi2     -> clear();
  mva           -> clear();
  trackIndex    -> clear();
}
