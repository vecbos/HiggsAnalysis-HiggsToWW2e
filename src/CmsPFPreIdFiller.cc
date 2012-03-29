//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsPFPreIdFiller
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
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFPreIdFiller.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TTree.h>
#include <TVector3.h>
#include <string>

using namespace edm;
using namespace reco;


//----------------------------------------
// -- Public Function Member Definitions --
//----------------------------------------

CmsPFPreIdFiller::CmsPFPreIdFiller(CmsTree *cmsTree, 
				   int maxTracks, int maxMCTracks,
				   bool noOutputIfLimitsReached):
  privateData_(new CmsPFPreIdFillerData)
{
  cmstree=cmsTree;
  
  trkIndexName_ = new std::string("n");  
  
  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_   = maxTracks;
  maxMCTracks_ = maxMCTracks;
  
  privateData_->initialise();
}

CmsPFPreIdFiller::CmsPFPreIdFiller(CmsTree *cmsTree, 
				   bool fatTree, int maxTracks, int maxMCTracks, 
				   bool noOutputIfLimitsReached):
  privateData_(new CmsPFPreIdFillerData)
{
  cmstree=cmsTree;
  
  trkIndexName_ = new std::string("n");
  
  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_   = maxTracks;
  maxMCTracks_ = maxMCTracks;
  
  privateData_->initialise();
}


CmsPFPreIdFiller::~CmsPFPreIdFiller() {
  
  delete privateData_->deltaEtaMatch;
  delete privateData_->deltaPhiMatch;
  delete privateData_->chiEtaMatch;
  delete privateData_->chiPhiMatch;
  delete privateData_->chi2Match;
  delete privateData_->eopMatch;
  delete privateData_->kfChi2;
  delete privateData_->kfNHits;
  delete privateData_->gsfChi2;
  delete privateData_->chi2Ratio;
  delete privateData_->ecalMatching;
  delete privateData_->psMatching;
  delete privateData_->trackFiltered;
  delete privateData_->preided;
  delete privateData_->trackIndex;
  delete privateData_->ncand;
  delete privateData_;
}

void CmsPFPreIdFiller::writeCollectionToTree(edm::InputTag PreIdMapLabel_, edm::InputTag TrackLabel_,
					     const edm::Event& iEvent, const edm::EventSetup& iSetup,
					     const std::string &columnPrefix, const std::string &columnSuffix,
					     bool dumpData) {
  
  edm::Handle<reco::TrackCollection> trackh;   
  try { iEvent.getByLabel(TrackLabel_,trackh); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFPreIdFiller") << "Can't get track collection: " << TrackLabel_; }
  const reco::TrackCollection & tracks = *(trackh.product());

  edm::Handle<edm::ValueMap<reco::PreIdRef> > vmaph;
  try { iEvent.getByLabel(PreIdMapLabel_,vmaph); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFPreIdFiller") << "Can't get preID collection: " << PreIdMapLabel_; }
  const edm::ValueMap<reco::PreIdRef> & preidMap = *(vmaph.product());

  privateData_->clearTrkVectors();
    
  int savedElements = 0;    

  if(&(*vmaph) && &(*trackh)) {

    for(unsigned itrack=0; itrack<tracks.size(); ++itrack) {
      
      reco::TrackRef theTrackRef(trackh,itrack);

      if(preidMap[theTrackRef].isNull()) continue;
      savedElements++;

      const reco::PreId & myPreId(*(preidMap[theTrackRef]));

      privateData_->deltaEtaMatch->push_back( myPreId.geomMatching()[0] );  
      privateData_->deltaPhiMatch->push_back( myPreId.geomMatching()[1] );  
      privateData_->chiEtaMatch->push_back( myPreId.geomMatching()[2] );  
      privateData_->chiPhiMatch->push_back( myPreId.geomMatching()[3] );  
      privateData_->chi2Match->push_back( myPreId.geomMatching()[4] );  
      privateData_->eopMatch->push_back( myPreId.eopMatch() );          
      privateData_->kfChi2->push_back( myPreId.kfChi2() );      
      privateData_->kfNHits->push_back( myPreId.kfNHits() );    
      privateData_->gsfChi2->push_back( myPreId.gsfChi2() );  
      privateData_->chi2Ratio->push_back( myPreId.chi2Ratio() );  
      privateData_->ecalMatching->push_back( myPreId.ecalMatching() ); 
      privateData_->psMatching->push_back( myPreId.esMatching() );
      privateData_->trackFiltered->push_back( myPreId.trackFiltered() ); 
      privateData_->preided->push_back( myPreId.preIded() );  
      privateData_->trackIndex->push_back(theTrackRef.key()); 
    }
  } // existing collections

  if(hitLimitsMeansNoOutput_ && savedElements>maxTracks_) {
    edm::LogInfo("CmsPFPreIdFiller") << "Track length " << vmaph->size() 
				     << " is too long for declared max length for tree "
				     <<  maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
    return;
  }
    
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the tree
  int blockSize = &(*vmaph) ? savedElements : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  treePFPreIdInfo(columnPrefix,columnSuffix);
  
  if(dumpData) cmstree->dumpData();
}

void CmsPFPreIdFiller::treePFPreIdInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"deltaEtaMatch"+colSuffix).c_str(), *privateData_->deltaEtaMatch, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiMatch"+colSuffix).c_str(), *privateData_->deltaPhiMatch, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chiEtaMatch"+colSuffix).c_str(),   *privateData_->chiEtaMatch,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chiPhiMatch"+colSuffix).c_str(),   *privateData_->chiPhiMatch,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chi2Match"+colSuffix).c_str(),     *privateData_->chi2Match,     nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eopMatch"+colSuffix).c_str(),      *privateData_->eopMatch,      nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"kfChi2"+colSuffix).c_str(),        *privateData_->kfChi2,        nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"kfNHits"+colSuffix).c_str(),       *privateData_->kfNHits,       nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"gsfChi2"+colSuffix).c_str(),       *privateData_->gsfChi2,       nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"chi2Ratio"+colSuffix).c_str(),     *privateData_->chi2Ratio,     nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"ecalMatching"+colSuffix).c_str(),  *privateData_->ecalMatching,  nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"psMatching"+colSuffix).c_str(),    *privateData_->psMatching,    nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"trackFiltered"+colSuffix).c_str(), *privateData_->trackFiltered, nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"preided"+colSuffix).c_str(),       *privateData_->preided,       nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(),    *privateData_->trackIndex,    nCandString.c_str(), 0, "Reco");  
}

void CmsPFPreIdFillerData::initialise() {
  
  ncand = new int;
  trackIndex = new vector<int>;
  deltaEtaMatch = new vector<float>;
  deltaPhiMatch = new vector<float>;
  chiEtaMatch = new vector<float>;
  chiPhiMatch = new vector<float>;
  chi2Match = new vector<float>;
  eopMatch = new vector<float>;
  kfChi2 = new vector<float>;
  kfNHits = new vector<float>;
  gsfChi2 = new vector<float>; 
  chi2Ratio = new vector<float>; 
  ecalMatching = new vector<bool>; 
  psMatching = new vector<bool>; 
  trackFiltered = new vector<bool>;
  preided = new vector<bool>;  
}

void CmsPFPreIdFillerData::clearTrkVectors() {

  deltaEtaMatch -> clear();
  deltaPhiMatch -> clear();
  chiEtaMatch   -> clear();
  chiPhiMatch   -> clear();
  chi2Match     -> clear();
  eopMatch      -> clear();
  kfChi2        -> clear();
  kfNHits       -> clear();
  ecalMatching  -> clear();
  psMatching    -> clear();
  trackFiltered -> clear();
  preided       -> clear();
  gsfChi2       -> clear();
  chi2Ratio     -> clear();
  trackIndex    -> clear();
}
