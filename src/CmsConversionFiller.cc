//---------------------------------------------------------------------------
//
// Description:
//       Package:   HiggsAnalysis/HiggsToWW2e
//       Class:     CmsConversionFiller
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

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsConversionFiller.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

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


CmsConversionFiller::CmsConversionFiller(CmsTree *cmsTree, 
					 int maxConv, int maxMCConv,
					 bool noOutputIfLimitsReached):
  privateData_(new CmsConversionFillerData)
{
  cmstree=cmsTree;

  convIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxConv_=maxConv;
  maxMCConv_=maxMCConv;

  privateData_->initialise();
}



//--------------
// Destructor --
//--------------

CmsConversionFiller::~CmsConversionFiller() { 
  delete privateData_->pxPair;
  delete privateData_->pyPair;
  delete privateData_->pzPair;

  delete privateData_->pxRefittedPair;
  delete privateData_->pyRefittedPair;
  delete privateData_->pzRefittedPair;
  
  delete privateData_->etaRefittedPair;
  delete privateData_->phiRefittedPair;
  delete privateData_->ptRefittedPair;
  delete privateData_->energyRefittedPair;
  delete privateData_->eOverPRefittedPair;
  delete privateData_->zOfPVFromTracks;
  
  delete privateData_->xVtx;
  delete privateData_->yVtx;
  delete privateData_->zVtx;
  delete privateData_->chi2Vtx;
  delete privateData_->chi2ProbVtx;
  delete privateData_->mvaOutVtx;
  delete privateData_->isValidVtx;
  delete privateData_->nTracksVtx;

  delete privateData_->trk1Dz;
  delete privateData_->trk1DzError;
  delete privateData_->trk1Charge;
  delete privateData_->trk1Algo;
  delete privateData_->trk1Pt;
  delete privateData_->trk1D0;
  delete privateData_->trk1Pout;
  delete privateData_->trk1Pin;

  delete privateData_->trk2Dz;
  delete privateData_->trk2DzError;
  delete privateData_->trk2Charge;
  delete privateData_->trk2Algo;
  delete privateData_->trk2Pt;
  delete privateData_->trk2D0;
  delete privateData_->trk2Pout;
  delete privateData_->trk2Pin;

  delete privateData_->ncand;
  delete privateData_;
}

//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsConversionFiller::writeCollectionToTree(edm::InputTag collectionTag,
					   const edm::Event& iEvent, const edm::EventSetup& iSetup,
					   const std::string &columnPrefix, const std::string &columnSuffix,
					   bool dumpData) {

//   edm::Handle< edm::View<reco::Candidate> > collectionHandle;
//   try { iEvent.getByLabel(collectionTag, collectionHandle); }
//   catch ( cms::Exception& ex ) { edm::LogWarning("CmsConversionFiller") << "Can't get candidate collection: " << collectionTag; }
//   const edm::View<reco::Candidate> *collection = collectionHandle.product();

  edm::Handle< edm::View<reco::Conversion> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsConversionFiller") << "Can't get conversion collection: " << collectionTag; }
  const edm::View<reco::Conversion> *collection = collectionHandle.product();

  privateData_->clearConvVectors();

  edm::ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get( magfield );     
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);       
  
  int blockSize=0;

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxConv_){
      edm::LogInfo("CmsConversionFiller") << "Conversion length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxConv_ << " and no output flag is set."
			       << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxConv_){
      edm::LogInfo("CmsConversionFiller") << "Conversion length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxConv_ 
			       << ". Collection will be truncated ";
    }

    for(int i=0; i < (int)collection->size(); i++) {
      RefToBase<reco::Conversion> convRef(collectionHandle, i);
      writeConvInfo(convRef , magfield, theTrackingGeometry );
    }
  } else {
    edm::LogWarning("CmsConversionFiller") << "No Conversions!";
    *(privateData_->ncand) = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 

  std::string nCandString = columnPrefix+(*convIndexName_)+columnSuffix; 
  if( (int)collection->size() ) blockSize = collection->size();
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");

  treeConvInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();

  delete convIndexName_;

}


void CmsConversionFiller::writeConvInfo(edm::RefToBase<reco::Conversion> convRef, const edm::ESHandle<MagneticField>& magfield, const edm::ESHandle<GlobalTrackingGeometry>& theTrackingGeometry) 
{

  if(convRef.isNonnull()) {
    privateData_->pxPair->push_back(convRef->pairMomentum().x());
    privateData_->pyPair->push_back(convRef->pairMomentum().y());
    privateData_->pzPair->push_back(convRef->pairMomentum().z());

    privateData_->pxRefittedPair->push_back(convRef->refittedPairMomentum().x());
    privateData_->pyRefittedPair->push_back(convRef->refittedPairMomentum().y());
    privateData_->pzRefittedPair->push_back(convRef->refittedPairMomentum().z());

    privateData_->etaRefittedPair->push_back(convRef->refittedPair4Momentum().eta());
    privateData_->phiRefittedPair->push_back(convRef->refittedPair4Momentum().phi());
    privateData_->ptRefittedPair->push_back(convRef->refittedPair4Momentum().pt());
    privateData_->energyRefittedPair->push_back(convRef->refittedPair4Momentum().energy());
    privateData_->eOverPRefittedPair->push_back(convRef->EoverPrefittedTracks());
    privateData_->zOfPVFromTracks->push_back(convRef->zOfPrimaryVertexFromTracks());

    privateData_->isValidVtx->push_back(convRef->conversionVertex().isValid());
    if(convRef->conversionVertex().isValid()){
      privateData_->xVtx->push_back(convRef->conversionVertex().x());
      privateData_->yVtx->push_back(convRef->conversionVertex().y());
      privateData_->zVtx->push_back(convRef->conversionVertex().z());

      privateData_->chi2Vtx->push_back(convRef->conversionVertex().chi2());
      privateData_->chi2ProbVtx->push_back(ChiSquaredProbability(convRef->conversionVertex().chi2(),convRef->conversionVertex().ndof()));

      privateData_->nTracksVtx->push_back(convRef->nTracks());
      privateData_->mvaOutVtx->push_back(convRef->MVAout());      
    }else{
      privateData_->xVtx->push_back(-1.);
      privateData_->yVtx->push_back(-1.);
      privateData_->zVtx->push_back(-1);

      privateData_->chi2Vtx->push_back(-1);
      privateData_->chi2ProbVtx->push_back(-1);

      privateData_->nTracksVtx->push_back(-1);
      privateData_->mvaOutVtx->push_back(-1);
    }

    //track info
    const std::vector<edm::RefToBase<reco::Track> > tracks = convRef->tracks();
    if(tracks.size()>0){ //at least 1 track
      privateData_->trk1Dz->push_back(tracks[0]->dz());
      privateData_->trk1DzError->push_back(tracks[0]->dzError());
      privateData_->trk1Charge->push_back(tracks[0]->charge());
      privateData_->trk1Algo->push_back(tracks[0]->algo());
      privateData_->trk1Pt->push_back(tracks[0]->pt());
      privateData_->trk1D0->push_back(convRef->tracksSigned_d0()[0]);
      privateData_->trk1Pout->push_back(sqrt(convRef->tracksPout()[0].Mag2()));
      privateData_->trk1Pin->push_back(sqrt(convRef->tracksPin()[0].Mag2()));
    }else{
      privateData_->trk1Dz->push_back(-99);
      privateData_->trk1DzError->push_back(-99);
      privateData_->trk1Charge->push_back(-99);
      privateData_->trk1Algo->push_back(-99);
      privateData_->trk1Pt->push_back(-99);
      privateData_->trk1D0->push_back(-99);
      privateData_->trk1Pout->push_back(-99);
      privateData_->trk1Pin->push_back(-99);
    }
    if(tracks.size()>1){ //at least 2 track
      privateData_->trk2Dz->push_back(tracks[1]->dz());
      privateData_->trk2DzError->push_back(tracks[1]->dzError());
      privateData_->trk2Charge->push_back(tracks[1]->charge());
      privateData_->trk2Algo->push_back(tracks[1]->algo());
      privateData_->trk2Pt->push_back(tracks[1]->pt());
      privateData_->trk2D0->push_back(convRef->tracksSigned_d0()[1]);
      privateData_->trk2Pout->push_back(sqrt(convRef->tracksPout()[1].Mag2()));
      privateData_->trk2Pin->push_back(sqrt(convRef->tracksPin()[1].Mag2()));
    }else{
      privateData_->trk2Dz->push_back(-99);
      privateData_->trk2DzError->push_back(-99);
      privateData_->trk2Charge->push_back(-99);
      privateData_->trk2Algo->push_back(-99);
      privateData_->trk2Pt->push_back(-99);
      privateData_->trk2D0->push_back(-99);
      privateData_->trk2Pout->push_back(-99);
      privateData_->trk2Pin->push_back(-99);
    }
  }
}


void CmsConversionFiller::treeConvInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*convIndexName_)+colSuffix;

  cmstree->column((colPrefix+"pxPair"+colSuffix).c_str(), *privateData_->pxPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pyPair"+colSuffix).c_str(), *privateData_->pyPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pzPair"+colSuffix).c_str(), *privateData_->pzPair, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"pxRefittedPair"+colSuffix).c_str(), *privateData_->pxRefittedPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pyRefittedPair"+colSuffix).c_str(), *privateData_->pyRefittedPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pzRefittedPair"+colSuffix).c_str(), *privateData_->pzRefittedPair, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"etaRefittedPair"+colSuffix).c_str(), *privateData_->etaRefittedPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phiRefittedPair"+colSuffix).c_str(), *privateData_->phiRefittedPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ptRefittedPair"+colSuffix).c_str(), *privateData_->ptRefittedPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energyRefittedPair"+colSuffix).c_str(), *privateData_->energyRefittedPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eOverPRefittedPair"+colSuffix).c_str(), *privateData_->eOverPRefittedPair, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"zOfPVFromTracks"+colSuffix).c_str(), *privateData_->zOfPVFromTracks, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"xVtx"+colSuffix).c_str(), *privateData_->xVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"yVtx"+colSuffix).c_str(), *privateData_->yVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"zVtx"+colSuffix).c_str(), *privateData_->zVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chi2Vtx"+colSuffix).c_str(), *privateData_->chi2Vtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chi2ProbVtx"+colSuffix).c_str(), *privateData_->chi2ProbVtx, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"isValidVtx"+colSuffix).c_str(), *privateData_->isValidVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nTracksVtx"+colSuffix).c_str(), *privateData_->nTracksVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mvaOutVtx"+colSuffix).c_str(), *privateData_->mvaOutVtx, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"trk1Dz"+colSuffix).c_str(), *privateData_->trk1Dz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk1DzError"+colSuffix).c_str(), *privateData_->trk1DzError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk1Charge"+colSuffix).c_str(), *privateData_->trk1Charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk1Algo"+colSuffix).c_str(), *privateData_->trk1Algo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk1Pt"+colSuffix).c_str(), *privateData_->trk1Pt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk1D0"+colSuffix).c_str(), *privateData_->trk1D0, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk1Pout"+colSuffix).c_str(), *privateData_->trk1Pout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk1Pin"+colSuffix).c_str(), *privateData_->trk1Pin, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"trk2Dz"+colSuffix).c_str(), *privateData_->trk2Dz, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk2DzError"+colSuffix).c_str(), *privateData_->trk2DzError, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk2Charge"+colSuffix).c_str(), *privateData_->trk2Charge, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk2Algo"+colSuffix).c_str(), *privateData_->trk2Algo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk2Pt"+colSuffix).c_str(), *privateData_->trk2Pt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk2D0"+colSuffix).c_str(), *privateData_->trk2D0, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk2Pout"+colSuffix).c_str(), *privateData_->trk2Pout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trk2Pin"+colSuffix).c_str(), *privateData_->trk2Pin, nCandString.c_str(), 0, "Reco");

}


void CmsConversionFillerData::initialise() {
  
  ncand = new int;
  pxPair = new vector<float>;
  pyPair = new vector<float>;
  pzPair = new vector<float>;

  pxRefittedPair = new vector<float>;
  pyRefittedPair = new vector<float>;
  pzRefittedPair = new vector<float>;

  etaRefittedPair = new vector<float>;
  phiRefittedPair = new vector<float>;
  ptRefittedPair = new vector<float>;
  energyRefittedPair = new vector<float>;
  eOverPRefittedPair = new vector<float>;
  zOfPVFromTracks = new vector<float>;

  xVtx = new vector<float>;
  yVtx = new vector<float>;
  zVtx = new vector<float>;
  chi2Vtx = new vector<float>;
  chi2ProbVtx = new vector<float>;
  mvaOutVtx = new vector<float>;
  isValidVtx = new vector<int>;
  nTracksVtx = new vector<int>;

  trk1Dz = new vector<float>;
  trk1DzError = new vector<float>;
  trk1Charge = new vector<float>;
  trk1Algo = new vector<float>;
  trk1Pt = new vector<float>;
  trk1D0 = new vector<float>;
  trk1Pout = new vector<float>;
  trk1Pin = new vector<float>;

  trk2Dz = new vector<float>;
  trk2DzError = new vector<float>;
  trk2Charge = new vector<float>;
  trk2Algo = new vector<float>;
  trk2Pt = new vector<float>;
  trk2D0 = new vector<float>;
  trk2Pout = new vector<float>;
  trk2Pin = new vector<float>;

}
void CmsConversionFillerData::clearConvVectors() {
  pxPair->clear();
  pyPair->clear();
  pzPair->clear();

  pxRefittedPair->clear();
  pyRefittedPair->clear();
  pzRefittedPair->clear();

  etaRefittedPair->clear();
  phiRefittedPair->clear();
  ptRefittedPair->clear();
  energyRefittedPair->clear();
  eOverPRefittedPair->clear();
  zOfPVFromTracks->clear();

  xVtx->clear();
  yVtx->clear();
  zVtx->clear();
  chi2Vtx->clear();
  chi2ProbVtx->clear();
  mvaOutVtx->clear();
  isValidVtx->clear();
  nTracksVtx->clear();
  trk1Dz->clear();
  trk1DzError->clear();
  trk1Charge->clear();
  trk1Algo->clear();
  trk1Pt->clear();
  trk1D0->clear();
  trk1Pout->clear();
  trk1Pin->clear();

  trk2Dz->clear();
  trk2DzError->clear();
  trk2Charge->clear();
  trk2Algo->clear();
  trk2Pt->clear();
  trk2D0->clear();
  trk2Pout->clear();
  trk2Pin->clear();
}

