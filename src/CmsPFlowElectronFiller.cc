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


CmsPFLowElectronFiller::CmsPFLowElectronFiller(CmsTree *cmsTree, 
					       int maxTracks, int maxMCTracks,
					       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFLowElectronFillerData)
{
  cmstree=cmsTree;

  // savePFEleKfTrk_=true;
  savePFEleGsfTrk_=true;
  savePFEleBasic_=true;
  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsPFLowElectronFiller::CmsPFLowElectronFiller(CmsTree *cmsTree, bool fatTree, 
					       int maxTracks, int maxMCTracks,
					       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFLowElectronFillerData)
{
  cmstree=cmsTree;

  savePFEleGsfTrk_=true;
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

CmsPFLowElectronFiller::~CmsPFLowElectronFiller() {

  

  delete privateData_->ncand;

  // gsf tracks
  delete privateData_->gsf_pxAtInner;
  delete privateData_->gsf_pyAtInner;
  delete privateData_->gsf_pzAtInner;
  delete privateData_->gsf_pxAtInnerMode;
  delete privateData_->gsf_pyAtInnerMode;
  delete privateData_->gsf_pzAtInnerMode;
  delete privateData_->gsf_xAtInner;
  delete privateData_->gsf_yAtInner;
  delete privateData_->gsf_zAtInner;
  delete privateData_->gsf_pxAtOuter;
  delete privateData_->gsf_pyAtOuter;
  delete privateData_->gsf_pzAtOuter;
  delete privateData_->gsf_xAtOuter;
  delete privateData_->gsf_yAtOuter;
  delete privateData_->gsf_zAtOuter;
  delete privateData_->gsf_TrackNormalizedChi2;
  delete privateData_->gsf_TrackDxy;
  delete privateData_->gsf_TrackD0;
  delete privateData_->gsf_TrackDsz;
  delete privateData_->gsf_TrackDz;
  delete privateData_->gsf_TrackDxyError;
  delete privateData_->gsf_TrackD0Error;
  delete privateData_->gsf_TrackDszError;
  delete privateData_->gsf_TrackDzError;
  delete privateData_->gsf_TrackValidHits;
  delete privateData_->gsf_TrackLostHits;
  delete privateData_->gsf_TrackVx;
  delete privateData_->gsf_TrackVy;
  delete privateData_->gsf_TrackVz;
  delete privateData_->gsf_charge;
  delete privateData_->gsf_chargeMode;
  delete privateData_->gsf_EcalDriven;
  delete privateData_->gsf_TrackerDriven; 


  // basic pf candidates info
  delete privateData_->MvaOutput;
  delete privateData_->PS1Energy;
  delete privateData_->PS2Energy;
  delete privateData_->EcalEnergy;
  delete privateData_->EcalElectronEnergy;
}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

//void CmsPFlowElectronFiller::savePFEleGsfTrk(bool what) { savePFEleGsfTrk_=what;}


void CmsPFLowElectronFiller::writeCollectionToTree(edm::InputTag collectionTag,
						   const edm::Event& iEvent, const edm::EventSetup& iSetup,
						   const std::string &columnPrefix, const std::string &columnSuffix,
						   bool dumpData) {
  
  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFLowElectronFiller") << "Can't get electron candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  privateData_->clearTrkVectors();
  
  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFLowElectronFiller") << "Track length " << collection->size() 
					     << " is too long for declared max length for tree "
					     << maxTracks_ << " and no output flag is set."
					     << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFLowElectronFiller") << "Track length " << collection->size() 
					     << " is too long for declared max length for tree "
					     << maxTracks_ 
					     << ". Collection will be truncated ";
    }
    
    *(privateData_->ncand) = collection->size();
    
    
    for(int index = 0; index < (int)collection->size(); index++) {
      
      // fill basic kinematics
      const Candidate *cand = &(collection->at(index));
      
      const PFCandidateRef pflowCandRef = collection->refAt(index).castTo<PFCandidateRef>();
      
      if(saveCand_) 
	writeCandInfo(cand,iEvent,iSetup);

      if ( !(pflowCandRef.isNull()) ) {
	GsfTrackRef gsfRef = pflowCandRef->gsfTrackRef();
	if(savePFEleBasic_)
	  writePFEleBasicInfo(pflowCandRef);

	if(savePFEleGsfTrk_) 
	  writePFEleGsfTrkInfo(gsfRef);
      }
      else {
	edm::LogWarning("CmsPFLowElectronFiller") << "Warning! The collection seems to be not made by "
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
  
  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  if(savePFEleBasic_) treePFEleBasicInfo(columnPrefix,columnSuffix);
  if(savePFEleGsfTrk_) treePFEleGsfTrkInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();

}

void CmsPFLowElectronFiller::writePFEleGsfTrkInfo(reco::GsfTrackRef gsfRef) {
  if(&gsfRef) {
    
    privateData_->gsf_TrackNormalizedChi2->push_back(gsfRef->normalizedChi2());
    
    privateData_->gsf_TrackDxy->push_back(gsfRef->dxy());
    privateData_->gsf_TrackD0->push_back(gsfRef->d0());
    privateData_->gsf_TrackDsz->push_back(gsfRef->dsz());
    privateData_->gsf_TrackDz->push_back(gsfRef->dz());
    
    privateData_->gsf_TrackDxyError->push_back(gsfRef->dxyError());
    privateData_->gsf_TrackD0Error->push_back(gsfRef->d0Error());
    privateData_->gsf_TrackDszError->push_back(gsfRef->dszError());
    privateData_->gsf_TrackDzError->push_back(gsfRef->dzError());
    
    privateData_->gsf_TrackValidHits->push_back(gsfRef->numberOfValidHits());
    privateData_->gsf_TrackLostHits->push_back(gsfRef->numberOfLostHits());
    
    privateData_->gsf_TrackVx->push_back(gsfRef->vx());
    privateData_->gsf_TrackVy->push_back(gsfRef->vy());
    privateData_->gsf_TrackVz->push_back(gsfRef->vz());
    
    privateData_->gsf_pxAtInner->push_back(gsfRef->innerMomentum().x());
    privateData_->gsf_pyAtInner->push_back(gsfRef->innerMomentum().y());
    privateData_->gsf_pzAtInner->push_back(gsfRef->innerMomentum().z());
   

    privateData_->gsf_pxAtInnerMode->push_back(gsfRef->pxMode());
    privateData_->gsf_pyAtInnerMode->push_back(gsfRef->pyMode());
    privateData_->gsf_pzAtInnerMode->push_back(gsfRef->pzMode());


    privateData_->gsf_charge->push_back(gsfRef->charge());
    privateData_->gsf_chargeMode->push_back(gsfRef->chargeMode());


    privateData_->gsf_xAtInner->push_back(gsfRef->innerPosition().x());
    privateData_->gsf_yAtInner->push_back(gsfRef->innerPosition().y());
    privateData_->gsf_zAtInner->push_back(gsfRef->innerPosition().z());
    
    privateData_->gsf_pxAtOuter->push_back(gsfRef->outerMomentum().x());
    privateData_->gsf_pyAtOuter->push_back(gsfRef->outerMomentum().y());
    privateData_->gsf_pzAtOuter->push_back(gsfRef->outerMomentum().z());
    
    privateData_->gsf_xAtOuter->push_back(gsfRef->outerPosition().x());
    privateData_->gsf_yAtOuter->push_back(gsfRef->outerPosition().y());
    privateData_->gsf_zAtOuter->push_back(gsfRef->outerPosition().z());
    
    //if (&(gsfRef.seedRef())) {
    ElectronSeedRef seedRef= gsfRef->extra()->seedRef().castTo<ElectronSeedRef>();
    if(seedRef.isNonnull()) {
      if(seedRef->caloCluster().isNonnull()) {
	privateData_->gsf_EcalDriven->push_back(1);
      }
      else {
	privateData_->gsf_EcalDriven->push_back(0);
      }
      if(seedRef->ctfTrack().isNonnull()) {
	privateData_->gsf_TrackerDriven->push_back(1);
      }
      else {
	privateData_->gsf_TrackerDriven->push_back(0);
      }
    }
    else {
      privateData_->gsf_EcalDriven->push_back(-1);
      privateData_->gsf_TrackerDriven->push_back(-1);
    }

    
  }
  else {
    
    privateData_->gsf_pxAtInner->push_back( -1.0 );
    privateData_->gsf_pyAtInner->push_back( -1.0 );
    privateData_->gsf_pzAtInner->push_back( -1.0 );
    
    privateData_->gsf_pxAtInnerMode->push_back(-1.0);
    privateData_->gsf_pyAtInnerMode->push_back(-1.0);
    privateData_->gsf_pzAtInnerMode->push_back(-1.0);


    privateData_->gsf_charge->push_back(-1.0);
    privateData_->gsf_chargeMode->push_back(-1.0);


    privateData_->gsf_xAtInner->push_back( -1.0 );
    privateData_->gsf_yAtInner->push_back( -1.0 );
    privateData_->gsf_zAtInner->push_back( -1.0 );
    
    privateData_->gsf_pxAtOuter->push_back( -1.0 );
    privateData_->gsf_pyAtOuter->push_back( -1.0 );
    privateData_->gsf_pzAtOuter->push_back( -1.0 );
    
    privateData_->gsf_xAtOuter->push_back( -1.0 );
    privateData_->gsf_yAtOuter->push_back( -1.0 );
    privateData_->gsf_zAtOuter->push_back( -1.0 );

    privateData_->gsf_TrackNormalizedChi2->push_back(-1.);

    privateData_->gsf_TrackDxy->push_back(-1.);
    privateData_->gsf_TrackD0->push_back(-1.);
    privateData_->gsf_TrackDsz->push_back(-1.);
    privateData_->gsf_TrackDz->push_back(-1.);

    privateData_->gsf_TrackDxyError->push_back(-1.);
    privateData_->gsf_TrackD0Error->push_back(-1.);
    privateData_->gsf_TrackDszError->push_back(-1.);
    privateData_->gsf_TrackDzError->push_back(-1.);

    privateData_->gsf_TrackValidHits->push_back(-1.);					       
    privateData_->gsf_TrackLostHits->push_back(-1.);

    privateData_->gsf_TrackVx->push_back(-1.);
    privateData_->gsf_TrackVy->push_back(-1.);
    privateData_->gsf_TrackVz->push_back(-1.);

    privateData_->gsf_pxAtInner->push_back(-1.);
    privateData_->gsf_pyAtInner->push_back(-1.);
    privateData_->gsf_pzAtInner->push_back(-1.);

    privateData_->gsf_xAtInner->push_back(-1.);
    privateData_->gsf_yAtInner->push_back(-1.);
    privateData_->gsf_zAtInner->push_back(-1.);

    privateData_->gsf_pxAtOuter->push_back(-1.);
    privateData_->gsf_pyAtOuter->push_back(-1.);
    privateData_->gsf_pzAtOuter->push_back(-1.);

    privateData_->gsf_xAtOuter->push_back(-1.);
    privateData_->gsf_yAtOuter->push_back(-1.);
    privateData_->gsf_zAtOuter->push_back(-1.);

    
    privateData_->gsf_EcalDriven->push_back(-1);
    privateData_->gsf_TrackerDriven->push_back(-1);

  }
}

void CmsPFLowElectronFiller::writePFEleBasicInfo(const reco::PFCandidateRef pflowCandRef) {
  
  privateData_->MvaOutput->push_back(pflowCandRef->mva_e_pi());
  privateData_->PS1Energy->push_back(pflowCandRef->pS1Energy());
  privateData_->PS2Energy->push_back(pflowCandRef->pS2Energy());
  privateData_->EcalEnergy->push_back(pflowCandRef->ecalEnergy());
  privateData_->EcalElectronEnergy->push_back(pflowCandRef->rawEcalEnergy());

}



void CmsPFLowElectronFiller::treePFEleGsfTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"gsf_TrackNormalizedChi2"+colSuffix).c_str(), *privateData_->gsf_TrackNormalizedChi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackDxy"+colSuffix).c_str(), *privateData_->gsf_TrackDxy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackD0"+colSuffix).c_str(),  *privateData_->gsf_TrackD0, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackDsz"+colSuffix).c_str(), *privateData_->gsf_TrackDsz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackDz"+colSuffix).c_str(),  *privateData_->gsf_TrackDz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackDxyError"+colSuffix).c_str(), *privateData_->gsf_TrackDxyError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackD0Error"+colSuffix).c_str(),  *privateData_->gsf_TrackD0Error, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackDszError"+colSuffix).c_str(), *privateData_->gsf_TrackDszError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackDzError"+colSuffix).c_str(),  *privateData_->gsf_TrackDzError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackValidHits"+colSuffix).c_str(),  *privateData_->gsf_TrackValidHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackLostHits"+colSuffix).c_str(),   *privateData_->gsf_TrackLostHits, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackVx"+colSuffix).c_str(),  *privateData_->gsf_TrackVx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackVy"+colSuffix).c_str(),  *privateData_->gsf_TrackVy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackVz"+colSuffix).c_str(),  *privateData_->gsf_TrackVz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pxAtInner"+colSuffix).c_str(), *privateData_->gsf_pxAtInner, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pyAtInner"+colSuffix).c_str(), *privateData_->gsf_pyAtInner, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pzAtInner"+colSuffix).c_str(), *privateData_->gsf_pzAtInner, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pxAtInnerMode"+colSuffix).c_str(), *privateData_->gsf_pxAtInnerMode, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pyAtInnerMode"+colSuffix).c_str(), *privateData_->gsf_pyAtInnerMode, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pzAtInnerMode"+colSuffix).c_str(), *privateData_->gsf_pzAtInnerMode, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_xAtInner"+colSuffix).c_str(), *privateData_->gsf_xAtInner, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_yAtInner"+colSuffix).c_str(), *privateData_->gsf_yAtInner, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_zAtInner"+colSuffix).c_str(), *privateData_->gsf_zAtInner, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pxAtOuter"+colSuffix).c_str(), *privateData_->gsf_pxAtOuter, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pyAtOuter"+colSuffix).c_str(), *privateData_->gsf_pyAtOuter, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_pzAtOuter"+colSuffix).c_str(), *privateData_->gsf_pzAtOuter, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_xAtOuter"+colSuffix).c_str(), *privateData_->gsf_xAtOuter, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_yAtOuter"+colSuffix).c_str(), *privateData_->gsf_yAtOuter, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_zAtOuter"+colSuffix).c_str(), *privateData_->gsf_zAtOuter, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_charge"+colSuffix).c_str(), *privateData_->gsf_charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_chargeMode"+colSuffix).c_str(), *privateData_->gsf_chargeMode, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_EcalDriven"+colSuffix).c_str(), *privateData_->gsf_EcalDriven, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsf_TrackerDriven"+colSuffix).c_str(), *privateData_->gsf_TrackerDriven, nCandString.c_str(), 0, "Reco");
  //  cmstree->column((colPrefix+""+colSuffix).c_str(), *privateData_->, nCandString.c_str(), 0, "Reco");

}

void CmsPFLowElectronFiller::treePFEleBasicInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"MvaOutput"+colSuffix).c_str(),  *privateData_->MvaOutput, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PS1Energy"+colSuffix).c_str(),  *privateData_->PS1Energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PS2Energy"+colSuffix).c_str(),  *privateData_->PS2Energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"EcalEnergy"+colSuffix).c_str(),  *privateData_->EcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"EcalElectronEnergy"+colSuffix).c_str(),  *privateData_->EcalElectronEnergy, nCandString.c_str(), 0, "Reco");

}



void CmsPFLowElectronFillerData::initialise() {  
  initialiseCandidate();

  // gsf tracks
  gsf_pxAtInner = new vector<float>;
  gsf_pyAtInner = new vector<float>;
  gsf_pzAtInner = new vector<float>;
  gsf_pxAtInnerMode = new vector<float>;
  gsf_pyAtInnerMode = new vector<float>;
  gsf_pzAtInnerMode = new vector<float>;
  gsf_xAtInner = new vector<float>;
  gsf_yAtInner = new vector<float>;
  gsf_zAtInner = new vector<float>;
  gsf_pxAtOuter = new vector<float>;
  gsf_pyAtOuter = new vector<float>;
  gsf_pzAtOuter = new vector<float>;
  gsf_xAtOuter = new vector<float>;
  gsf_yAtOuter = new vector<float>;
  gsf_zAtOuter = new vector<float>;
  gsf_TrackNormalizedChi2 = new vector<float>;
  gsf_TrackDxy = new vector<float>;
  gsf_TrackD0 = new vector<float>;
  gsf_TrackDsz = new vector<float>;
  gsf_TrackDz = new vector<float>;
  gsf_TrackDxyError = new vector<float>;
  gsf_TrackD0Error = new vector<float>;
  gsf_TrackDszError = new vector<float>;
  gsf_TrackDzError = new vector<float>;
  gsf_TrackValidHits = new vector<float>;
  gsf_TrackLostHits = new vector<float>;
  gsf_TrackVx = new vector<float>;
  gsf_TrackVy = new vector<float>;
  gsf_TrackVz = new vector<float>;
  gsf_charge = new vector<float>;
  gsf_chargeMode = new vector<float>;
  gsf_EcalDriven = new vector<int>;
  gsf_TrackerDriven = new vector<int>;


  // basic pf candidates info
  MvaOutput = new vector<float>;
  PS1Energy = new vector<float>;
  PS2Energy = new vector<float>;
  EcalEnergy = new vector<float>;
  EcalElectronEnergy = new vector<float>;


}
void CmsPFLowElectronFillerData::clearTrkVectors() {

  // gsf track
  clearTrkVectorsCandidate();
  gsf_pxAtOuter->clear();
  gsf_pyAtOuter->clear();
  gsf_pzAtOuter->clear();
  gsf_xAtOuter->clear();
  gsf_yAtOuter->clear();
  gsf_zAtOuter->clear();
  gsf_pxAtInner->clear();
  gsf_pyAtInner->clear();
  gsf_pzAtInner->clear();
  gsf_pxAtInnerMode->clear();
  gsf_pyAtInnerMode->clear();
  gsf_pzAtInnerMode->clear();
  gsf_xAtInner->clear();
  gsf_yAtInner->clear();
  gsf_zAtInner->clear();
  gsf_TrackNormalizedChi2->clear();
  gsf_TrackDxy->clear();
  gsf_TrackD0->clear();
  gsf_TrackDsz->clear();
  gsf_TrackDz->clear();
  gsf_TrackDxyError->clear();
  gsf_TrackD0Error->clear();
  gsf_TrackDszError->clear();
  gsf_TrackDzError->clear();
  gsf_TrackValidHits->clear();
  gsf_TrackLostHits->clear();
  gsf_TrackVx->clear();
  gsf_TrackVy->clear();
  gsf_TrackVz->clear();
  gsf_charge->clear();
  gsf_chargeMode->clear();
  gsf_EcalDriven->clear();
  gsf_TrackerDriven->clear();

  // basic pf candidates info
  MvaOutput->clear();
  PS1Energy->clear();
  PS2Energy->clear();
  EcalEnergy->clear();
  EcalElectronEnergy->clear();
}
