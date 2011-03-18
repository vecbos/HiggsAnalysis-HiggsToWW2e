//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsVertexFiller
//
// Original Author:  Chris Rogan
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"


#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsVertexFiller.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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


CmsVertexFiller::CmsVertexFiller(CmsTree *cmsTree,
				 bool fatTree, 
				 int maxTracks, int maxMCTracks,
				 bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,fatTree,maxTracks,maxMCTracks,noOutputIfLimitsReached), 
  privateData_(new CmsVertexFillerData)
{
  cmstree=cmsTree;
  trkIndexName_ = new std::string("n");
  privateData_->initialise();
}


// ---------------------------------------------------------------


CmsVertexFiller::~CmsVertexFiller() {
  delete privateData_->PVx;
  delete privateData_->PVy;
  delete privateData_->PVz;
  delete privateData_->PVErrx;
  delete privateData_->PVErry;
  delete privateData_->PVErrz;
  delete privateData_->SumPt;
  delete privateData_->ndof;
  delete privateData_->chi2;
  delete privateData_->isFake;  
}


// ---------------------------------------------------------------
void
CmsVertexFiller::writeCollectionToTree (edm::InputTag vtxcollectionTag,
					const edm::Event& iEvent, 
					const edm::EventSetup& iSetup,
					const std::string &columnPrefix,
					const std::string &columnSuffix)
{
  privateData_->clearTrkVectors();
  
  int nVtx = 0;

  Handle<reco::VertexCollection> primaryVertex;
  try { iEvent.getByLabel(vtxcollectionTag, primaryVertex); }
  catch(cms::Exception& ex ) {edm::LogWarning("CmsVertexFiller") << "Can't get candidate collection: " << vtxcollectionTag; }

  if(primaryVertex->size() > 0) {
    for(VertexCollection::const_iterator v = primaryVertex->begin();
	v != primaryVertex->end(); ++v){
      float SumPt = 0.0;
      if((*v).tracksSize() > 0){
	std::vector<TrackBaseRef >::const_iterator t;
	for( t = (*v).tracks_begin(); t != (*v).tracks_end(); t++){
	  if((**t).charge() < -1 || (**t).charge() > 1){
	    //illegal charge
	  } else {
	    SumPt += (**t).pt();
	  }
	}
	privateData_->PVx->push_back((*v).x());
	privateData_->PVy->push_back((*v).y());
	privateData_->PVz->push_back((*v).z());
	privateData_->PVErrx->push_back((*v).xError());
	privateData_->PVErry->push_back((*v).yError());
	privateData_->PVErrz->push_back((*v).zError());
	privateData_->SumPt->push_back(SumPt);
	privateData_->ndof->push_back((*v).ndof());
	privateData_->chi2->push_back((*v).chi2());
	if ((*v).isFake())   privateData_->isFake->push_back(1);
	if (!(*v).isFake())  privateData_->isFake->push_back(0);
	nVtx++;
      }
    }
  } else {
    privateData_->PVx->push_back(-1.);
    privateData_->PVy->push_back(-1.);
    privateData_->PVz->push_back(-1.);
    privateData_->PVErrx->push_back(-1.);
    privateData_->PVErry->push_back(-1.);
    privateData_->PVErrz->push_back(-1.);
    privateData_->SumPt->push_back(-1.);
    privateData_->ndof->push_back(-1.);
    privateData_->chi2->push_back(-1.);
    privateData_->isFake->push_back(-1.);
  }

  *(privateData_->ncand) = nVtx;

  int blockSize = nVtx;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix;
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");

  treeVertexInfo(columnPrefix,columnSuffix);
  
}

// ---------------------------------------------------------------

void CmsVertexFiller::treeVertexInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  
  cmstree->column((colPrefix+"PVx"+colSuffix).c_str(), *privateData_->PVx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PVy"+colSuffix).c_str(), *privateData_->PVy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PVz"+colSuffix).c_str(), *privateData_->PVz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PVErrx"+colSuffix).c_str(), *privateData_->PVErrx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PVErry"+colSuffix).c_str(), *privateData_->PVErry, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PVErrz"+colSuffix).c_str(), *privateData_->PVErrz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"SumPt"+colSuffix).c_str(), *privateData_->SumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ndof"+colSuffix).c_str(), *privateData_->ndof, nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"chi2"+colSuffix).c_str(), *privateData_->chi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"isFake"+colSuffix).c_str(), *privateData_->isFake, nCandString.c_str(), 0, "Reco");
  
}

void CmsVertexFillerData::initialise() {
  initialiseCandidate();
  PVx = new vector<float>;
  PVy = new vector<float>;
  PVz = new vector<float>;
  PVErrx = new vector<float>;
  PVErry = new vector<float>;
  PVErrz = new vector<float>;
  SumPt = new vector<float>;
  ndof = new vector<float>;
  chi2 = new vector<float>;
  isFake = new vector<int>;
}

void CmsVertexFillerData::clearTrkVectors() {
  clearTrkVectorsCandidate();
  PVx->clear();
  PVy->clear();
  PVz->clear();
  PVErrx->clear();
  PVErry->clear();
  PVErrz->clear();
  SumPt->clear();
  ndof->clear();
  chi2->clear();
  isFake->clear();

}

