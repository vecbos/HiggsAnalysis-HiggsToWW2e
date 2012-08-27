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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsVertexFiller.h"

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
  delete privateData_->rho;
  delete privateData_->ndof;
  delete privateData_->chi2;
  delete privateData_->normChi2;
  delete privateData_->pxChMet;
  delete privateData_->pyChMet;
  delete privateData_->pzChMet;
  delete privateData_->isFake;  
  delete privateData_->isValid;  
  delete privateData_->trackSize;  
  delete privateData_;
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

  edm::Handle< edm::ValueMap<reco::PFMET> > chargedMets;
  if(ChargedMets_.label().size()!=0) {
    try { iEvent.getByLabel(ChargedMets_, chargedMets); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsVertexFiller") << "Can't get charged PFMET candidate collection: " << ChargedMets_; }
  }
  if(primaryVertex->size() > 0) {
    int ivtx=0;
    for(VertexCollection::const_iterator v = primaryVertex->begin();
	v != primaryVertex->end(); ++v){
      float SumPt = 0.0;
      if(tracksCollection_.label().size()!=0) {
        if((*v).tracksSize() > 0){
          std::vector<TrackBaseRef >::const_iterator t;
          for( t = (*v).tracks_begin(); t != (*v).tracks_end(); t++){
            if((**t).charge() < -1 || (**t).charge() > 1){
              //illegal charge
            } else {
              SumPt += (**t).pt();
            }
          }
        }

        reco::VertexRef vertexRef(primaryVertex, ivtx);
        
	privateData_->PVx->push_back((*v).x());
	privateData_->PVy->push_back((*v).y());
	privateData_->PVz->push_back((*v).z());
	privateData_->PVErrx->push_back((*v).xError());
	privateData_->PVErry->push_back((*v).yError());
	privateData_->PVErrz->push_back((*v).zError());
	privateData_->SumPt->push_back(SumPt);
	privateData_->rho->push_back((*v).position().Rho());
	privateData_->ndof->push_back((*v).ndof());
	privateData_->chi2->push_back((*v).chi2());
	privateData_->normChi2->push_back((*v).normalizedChi2());
        if(vertexRef.isNonnull() && ChargedMets_.label().size()!=0) {
          privateData_->pxChMet->push_back(((*chargedMets)[vertexRef]).px());
          privateData_->pyChMet->push_back(((*chargedMets)[vertexRef]).py());
          privateData_->pzChMet->push_back(((*chargedMets)[vertexRef]).pz());
        } else {
          privateData_->pxChMet->push_back(-1000.);
          privateData_->pyChMet->push_back(-1000.);
          privateData_->pzChMet->push_back(-1000.);          
        }
	if ((*v).isFake())   privateData_->isFake->push_back(1);
	if (!(*v).isFake())  privateData_->isFake->push_back(0);
	privateData_->isValid->push_back((*v).isValid());
	privateData_->trackSize->push_back((*v).tracksSize());
	
	nVtx++;
      }
      ivtx++;
    }
  } else {
    privateData_->PVx->push_back(-1000.);
    privateData_->PVy->push_back(-1000.);
    privateData_->PVz->push_back(-1000.);
    privateData_->PVErrx->push_back(-1000.);
    privateData_->PVErry->push_back(-1000.);
    privateData_->PVErrz->push_back(-1000.);
    privateData_->SumPt->push_back(-1000.);
    privateData_->rho->push_back(-1000.);
    privateData_->ndof->push_back(-1000.);
    privateData_->chi2->push_back(-1000.);
    privateData_->normChi2->push_back(-1000.);
    privateData_->pxChMet->push_back(-1000.);
    privateData_->pyChMet->push_back(-1000.);
    privateData_->pzChMet->push_back(-1000.);
    privateData_->isFake->push_back(-1000.);
    privateData_->isValid->push_back(-1000.);
    privateData_->trackSize->push_back(-1000.);
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
  cmstree->column((colPrefix+"rho"+colSuffix).c_str(), *privateData_->rho, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ndof"+colSuffix).c_str(), *privateData_->ndof, nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"chi2"+colSuffix).c_str(), *privateData_->chi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"normChi2"+colSuffix).c_str(), *privateData_->normChi2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pxChMet"+colSuffix).c_str(), *privateData_->pxChMet, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pyChMet"+colSuffix).c_str(), *privateData_->pyChMet, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pzChMet"+colSuffix).c_str(), *privateData_->pzChMet, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"isFake"+colSuffix).c_str(), *privateData_->isFake, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"isValid"+colSuffix).c_str(), *privateData_->isValid, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackSize"+colSuffix).c_str(), *privateData_->trackSize, nCandString.c_str(), 0, "Reco");
  
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
  rho = new vector<float>;
  ndof = new vector<float>;
  chi2 = new vector<float>;
  normChi2 = new vector<float>;
  pxChMet = new vector<float>;
  pyChMet = new vector<float>;
  pzChMet = new vector<float>;
  isFake = new vector<int>;
  isValid = new vector<int>;
  trackSize = new vector<int>;
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
  rho->clear();
  ndof->clear();
  chi2->clear();
  normChi2->clear();
  pxChMet->clear();
  pyChMet->clear();
  pzChMet->clear();
  isFake->clear();
  isValid->clear();
  trackSize->clear();
}

