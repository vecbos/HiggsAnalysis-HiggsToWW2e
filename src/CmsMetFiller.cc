//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsMetFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  9 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMetFiller.h"
#include "DataFormats/METReco/interface/AnomalousECALVariables.h"

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


CmsMetFiller::CmsMetFiller(CmsTree *cmsTree, int maxTracks, 
				       int maxMCTracks,
				       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsMetFillerData)
{
  cmstree=cmsTree;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  isData_ = false;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsMetFiller::~CmsMetFiller() {

  // delete here the vector ptr's
  delete privateData_->sumEt;
  delete privateData_->mEtSig;
  delete privateData_->significance;
  delete privateData_->filterBits;
  delete privateData_;
}


//-------------
// Methods   --
//-------------




void CmsMetFiller::writeCollectionToTree(edm::InputTag collectionTag,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  privateData_->clearTrkVectors();

  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  LogInfo("CmsMetFiller") << "=== Writing collection " << nCandString << " ===";

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      LogError("CmsMetFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ << " and no output flag is set."
				     << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      LogError("CmsMetFiller") << "Track length " << collection->size() 
				     << " is too long for declared max length for tree "
				     << maxTracks_ 
				     << ". Collection will be truncated ";
    }

    if(isData_) {
      // take the noise filter bits
      // ECAL MET optional filters
      edm::Handle< bool > ECALTPFilter;
      try { iEvent.getByLabel("EcalDeadCellEventFlagProducer", ECALTPFilter); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << ECALTPFilter; }
      bool ECALTPFilterFlag = *ECALTPFilter;

      edm::Handle< bool > HBHENoiseFilterResult;
      try { iEvent.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"), HBHENoiseFilterResult); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << HBHENoiseFilterResult; }
      bool HBHENoiseFilterResultFlag = *HBHENoiseFilterResult;

      edm::Handle< bool > hcalLaserEventFilter;
      try { iEvent.getByLabel("hcalLaserEventFilter", hcalLaserEventFilter); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << hcalLaserEventFilter; }
      bool hcalLaserEventFilterFlag = *hcalLaserEventFilter;

      edm::Handle< bool > eeBadScFilter;
      try { iEvent.getByLabel("eeBadScFilter", eeBadScFilter); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << eeBadScFilter; }
      bool eeBadScFilterFlag = *eeBadScFilter;

      edm::Handle< int > ECALDeadDRFilter;
      try { iEvent.getByLabel("simpleDRFlagProducer","deadCellStatus", ECALDeadDRFilter); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get int: " << ECALDeadDRFilter; }
      int ECALDeadDRFilterFlag = *ECALDeadDRFilter;
    
      edm::Handle< int > ECALBoundaryDRFilter;
      try { iEvent.getByLabel("simpleDRFlagProducer","boundaryStatus", ECALBoundaryDRFilter); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get int: " << ECALBoundaryDRFilter; }
      int ECALBoundaryDRFilterFlag = *ECALBoundaryDRFilter;

      int drDead = (ECALDeadDRFilterFlag>0) ? 1 : 0;
      int drBoundary = (ECALBoundaryDRFilterFlag>0) ? 1 : 0;

      edm::Handle< BeamHaloSummary > beamHaloH;
      iEvent.getByLabel("BeamHaloSummary", beamHaloH);
      bool CSCHaloFilterFlag = ! beamHaloH->CSCTightHaloId();

      edm::Handle< bool > trackerFailureFilter;
      try { iEvent.getByLabel("trackingFailureFilter", trackerFailureFilter); }
      catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << trackerFailureFilter; }
      bool trackerFailureFilterFlag = *trackerFailureFilter;    

      edm::InputTag ecalAnomalousFilterTag("EcalAnomalousEventFilter","anomalousECALVariables");
      Handle<AnomalousECALVariables> anomalousECALvarsHandle;
      iEvent.getByLabel(ecalAnomalousFilterTag, anomalousECALvarsHandle);
      AnomalousECALVariables anomalousECALvars;
      if (anomalousECALvarsHandle.isValid()) {
        anomalousECALvars = *anomalousECALvarsHandle;
      } else {
        edm::LogWarning("anomalous ECAL Vars not valid/found");
      } 

      bool isNotDeadEcalCluster = !(anomalousECALvars.isDeadEcalCluster());
      
      *(privateData_->filterBits) = 
        (eeBadScFilterFlag << 8) | (hcalLaserEventFilterFlag << 7) | (HBHENoiseFilterResultFlag << 6) | 
        (isNotDeadEcalCluster << 5) | (trackerFailureFilterFlag << 4) | (CSCHaloFilterFlag << 3) | 
        ( drDead << 2 ) | ( drBoundary << 1 ) | ECALTPFilterFlag;
    }

    *(privateData_->ncand) = collection->size();

    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);
      const MET *thisMET = dynamic_cast< const MET * > ( &(*cand) );
      privateData_->sumEt->push_back(thisMET->sumEt());
      privateData_->mEtSig->push_back(thisMET->mEtSig());
      privateData_->significance->push_back(thisMET->significance());
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
  treeMetInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();

  delete trkIndexName_;

}







void CmsMetFiller::treeMetInfo(const std::string colPrefix, const std::string colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"sumEt"+colSuffix).c_str(), *privateData_->sumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mEtSig"+colSuffix).c_str(), *privateData_->mEtSig, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"significance"+colSuffix).c_str(), *privateData_->significance, nCandString.c_str(), 0, "Reco");
  
  if(isData_) cmstree->column("METFlags", *privateData_->filterBits, 0, "Reco");
  
}







void CmsMetFillerData::initialise() {

  initialiseCandidate();
  sumEt = new vector<float>;
  mEtSig = new vector<float>;
  significance = new vector<float>;
  filterBits = new int;

}

void CmsMetFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();
  sumEt->clear();
  mEtSig->clear();
  significance->clear();

}
