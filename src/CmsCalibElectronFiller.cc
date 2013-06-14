//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsCalibElectronFiller
//
// Original Author:  Chiara Rovelli
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCalibElectronFiller.h"

#include <TTree.h>

#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <string>

using namespace edm;
using namespace reco;


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------

CmsCalibElectronFiller::CmsCalibElectronFiller(CmsTree *cmsTree, 
					       int maxTracks, int maxMCTracks,
					       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsCalibElectronFillerData)
{
  cmstree=cmsTree;
  
  trkIndexName_ = new std::string("n"); 
  
  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsCalibElectronFiller::~CmsCalibElectronFiller() {

  delete privateData_->energyError;  
  delete privateData_->ncand;
  delete privateData_;
}

void CmsCalibElectronFiller::writeCollectionToTree(edm::InputTag collectionTag,
						   const edm::Event& iEvent, const edm::EventSetup& iSetup,
						   const std::string &columnPrefix, const std::string &columnSuffix,
						   bool dumpData) {
  
  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsCalibElectronFiller") << "Can't get calibrated electron candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  privateData_->clearTrkVectors();
  
  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsCalibElectronFiller") << "Track length " << collection->size() 
					     << " is too long for declared max length for tree "
					     << maxTracks_ << " and no output flag is set."
					     << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsCalibElectronFiller") << "Track length " << collection->size() 
					     << " is too long for declared max length for tree "
					     << maxTracks_ 
					     << ". Collection will be truncated ";
    }
    
    *(privateData_->ncand) = collection->size();

    for(int index = 0; index < (int)collection->size(); index++) {
      
      // fill basic kinematics
      const Candidate *cand = &(collection->at(index));
      if(saveCand_) writeCandInfo(cand,iEvent,iSetup);

      // getting final energy error
      const GsfElectronRef electronRef = collection->refAt(index).castTo<GsfElectronRef>();
      
      if ( !(electronRef.isNull()) ) {
	privateData_->energyError->push_back( electronRef->p4Error( (electronRef->candidateP4Kind()) ) );
      } else {
	privateData_->energyError->push_back( -999. );
        edm::LogWarning("CmsCalibElectronFiller") << "Warning! The collection seems to be not made by "
						  << "electrons, electron-specific infos will be set to default.";
      }
    }
  }
  else {
    *(privateData_->ncand) = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  treeCandInfo(columnPrefix,columnSuffix);
  treeExtraEleInfo(columnPrefix,columnSuffix);

  cmstree->dumpData();

  delete trkIndexName_;
}

void CmsCalibElectronFiller::treeExtraEleInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"energyError"+colSuffix).c_str(),  *privateData_->energyError, nCandString.c_str(), 0, "Reco");
}

void CmsCalibElectronFillerData::initialise() {
  
  initialiseCandidate();
  energyError = new vector<float>;
}

void CmsCalibElectronFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  energyError->clear();
}
