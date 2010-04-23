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

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFlowElectronFiller.h"

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


CmsPFlowElectronFiller::CmsPFlowElectronFiller(CmsTree *cmsTree, 
					       int maxTracks, int maxMCTracks,
					       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFlowElectronFillerData)
{
  cmstree=cmsTree;

  savePFEleTrk_=true;     
  savePFEleBasic_=true;
  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsPFlowElectronFiller::CmsPFlowElectronFiller(CmsTree *cmsTree, bool fatTree, 
					       int maxTracks, int maxMCTracks,
					       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFlowElectronFillerData)
{
  cmstree=cmsTree;

  savePFEleTrk_=true;     
  savePFEleBasic_=true;
  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();

}


//--------------
// Destructor --
//--------------

CmsPFlowElectronFiller::~CmsPFlowElectronFiller() {

  delete privateData_->ncand;
  
  delete privateData_->trackIndex;
  delete privateData_->gsfTrackIndex;

  delete privateData_->MvaOutput;
  delete privateData_->PS1Energy;
  delete privateData_->PS2Energy;
  delete privateData_->EcalEnergy;  
  delete privateData_->RawEcalEnergy;
  delete privateData_->HcalEnergy;  
  delete privateData_->RawHcalEnergy;  
  delete privateData_->PositionAtEcalX;
  delete privateData_->PositionAtEcalY;
  delete privateData_->PositionAtEcalZ;
}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsPFlowElectronFiller::writeCollectionToTree(edm::InputTag collectionTag,
						   const edm::Event& iEvent, const edm::EventSetup& iSetup,
						   const std::string &columnPrefix, const std::string &columnSuffix,
						   bool dumpData) {
  
  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFlowElectronFiller") << "Can't get electron candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  privateData_->clearTrkVectors();
  
  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFlowElectronFiller") << "Track length " << collection->size() 
					     << " is too long for declared max length for tree "
					     << maxTracks_ << " and no output flag is set."
					     << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFlowElectronFiller") << "Track length " << collection->size() 
					     << " is too long for declared max length for tree "
					     << maxTracks_ 
					     << ". Collection will be truncated ";
    }
    
    *(privateData_->ncand) = collection->size();
    
    
    for(int index = 0; index < (int)collection->size(); index++) {
      
      // fill basic kinematics
      const Candidate *cand = &(collection->at(index));
      const PFCandidateRef pflowCandRef = collection->refAt(index).castTo<PFCandidateRef>();
      if(saveCand_) writeCandInfo(cand,iEvent,iSetup);
      
      if ( !(pflowCandRef.isNull()) ) {
	
	// basic PF infos
	if(savePFEleBasic_) writePFEleBasicInfo(pflowCandRef);
	
	// tracker based infos
	GsfTrackRef gsfRef  = pflowCandRef->gsfTrackRef();
	TrackRef kfTrackRef = pflowCandRef->trackRef();
	if(savePFEleTrk_) writePFEleTrkInfo(gsfRef,kfTrackRef);
      }
      else {
	edm::LogWarning("CmsPFlowElectronFiller") << "Warning! The collection seems to be not made by "
						  << "pflow candidates electrons, electron-specific infos will be set to default.";
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
  
  if(saveCand_)       treeCandInfo(columnPrefix,columnSuffix);
  if(savePFEleBasic_) treePFEleBasicInfo(columnPrefix,columnSuffix);
  if(savePFEleTrk_)   treePFEleTrkInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();

}

void CmsPFlowElectronFiller::writePFEleTrkInfo(reco::GsfTrackRef gsfRef, reco::TrackRef kfTrackRef ) {
  
  if(gsfRef.isNonnull()) {
    privateData_->gsfTrackIndex->push_back(gsfRef.key());
  }
  else {
    privateData_->gsfTrackIndex->push_back( -1 );
  }
  
  if(kfTrackRef.isNonnull()) {
    privateData_->trackIndex->push_back(kfTrackRef.key());
  } else {
    privateData_->trackIndex->push_back( -1 );
  }
}

void CmsPFlowElectronFiller::writePFEleBasicInfo(const reco::PFCandidateRef pflowCandRef) {

  if (pflowCandRef.isNonnull()) {
    
    privateData_->MvaOutput->push_back(pflowCandRef->mva_e_pi());
    privateData_->PS1Energy->push_back(pflowCandRef->pS1Energy());
    privateData_->PS2Energy->push_back(pflowCandRef->pS2Energy());
    privateData_->EcalEnergy->push_back(pflowCandRef->ecalEnergy());
    privateData_->HcalEnergy->push_back(pflowCandRef->hcalEnergy());    
    privateData_->RawEcalEnergy->push_back(pflowCandRef->rawEcalEnergy());
    privateData_->RawHcalEnergy->push_back(pflowCandRef->rawHcalEnergy());    
    privateData_->PositionAtEcalX->push_back(pflowCandRef->positionAtECALEntrance().x());
    privateData_->PositionAtEcalY->push_back(pflowCandRef->positionAtECALEntrance().y());
    privateData_->PositionAtEcalZ->push_back(pflowCandRef->positionAtECALEntrance().z());

  } else {  
    
    privateData_->MvaOutput->push_back(-1.);
    privateData_->PS1Energy->push_back(-1.);
    privateData_->PS2Energy->push_back(-1.);
    privateData_->EcalEnergy->push_back(-1.);
    privateData_->HcalEnergy->push_back(-1.);
    privateData_->RawEcalEnergy->push_back(-1.);
    privateData_->RawHcalEnergy->push_back(-1.);
    privateData_->PositionAtEcalX->push_back(-1.);
    privateData_->PositionAtEcalY->push_back(-1.);
    privateData_->PositionAtEcalZ->push_back(-1.);
  }
}



void CmsPFlowElectronFiller::treePFEleTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"gsfTrackIndex"+colSuffix).c_str(), *privateData_->gsfTrackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(), *privateData_->trackIndex, nCandString.c_str(), 0, "Reco");
}

void CmsPFlowElectronFiller::treePFEleBasicInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"MvaOutput"+colSuffix).c_str(),  *privateData_->MvaOutput, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PS1Energy"+colSuffix).c_str(),  *privateData_->PS1Energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PS2Energy"+colSuffix).c_str(),  *privateData_->PS2Energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"EcalEnergy"+colSuffix).c_str(),  *privateData_->EcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HcalEnergy"+colSuffix).c_str(),  *privateData_->HcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"RawEcalEnergy"+colSuffix).c_str(),  *privateData_->RawEcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"RawHcalEnergy"+colSuffix).c_str(),  *privateData_->RawHcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PositionAtEcalX"+colSuffix).c_str(),  *privateData_->PositionAtEcalX, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PositionAtEcalY"+colSuffix).c_str(),  *privateData_->PositionAtEcalY, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PositionAtEcalZ"+colSuffix).c_str(),  *privateData_->PositionAtEcalZ, nCandString.c_str(), 0, "Reco");
}



void CmsPFlowElectronFillerData::initialise() {  
  initialiseCandidate();

  // track infos
  trackIndex = new vector<int>;
  gsfTrackIndex = new vector<int>;  

  // basic pf candidates info
  MvaOutput = new vector<float>;
  PS1Energy = new vector<float>;
  PS2Energy = new vector<float>;
  EcalEnergy = new vector<float>;
  HcalEnergy = new vector<float>;
  RawEcalEnergy = new vector<float>;
  RawHcalEnergy = new vector<float>;
  PositionAtEcalX = new vector<float>;
  PositionAtEcalY = new vector<float>;
  PositionAtEcalZ = new vector<float>;
}
void CmsPFlowElectronFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  trackIndex->clear();
  gsfTrackIndex->clear();

  // basic pf candidates info
  MvaOutput->clear();
  PS1Energy->clear();
  PS2Energy->clear();
  EcalEnergy->clear();
  HcalEnergy->clear();
  RawEcalEnergy->clear();
  RawHcalEnergy->clear();
  PositionAtEcalX->clear();
  PositionAtEcalY->clear();
  PositionAtEcalZ->clear();
}
