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

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPhotonFiller.h"

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


CmsPhotonFiller::CmsPhotonFiller(CmsTree *cmsTree, 
                                 int maxTracks, int maxMCTracks,
                                 bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPhotonFillerData)
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

CmsPhotonFiller::~CmsPhotonFiller() {

  // delete here the vector ptr's
  delete privateData_->fiducialFlags;
  delete privateData_->recoFlags;
  delete privateData_->superClusterIndex;
  delete privateData_->PFsuperClusterIndex;
  
  delete privateData_->hOverE;
  delete privateData_->dr03HollowTkSumPt;
  delete privateData_->dr03TkSumPt;
  delete privateData_->dr03EcalRecHitSumEt;
  delete privateData_->dr03HcalTowerSumEt;
  delete privateData_->dr04HollowTkSumPt;
  delete privateData_->dr04TkSumPt;
  delete privateData_->dr04EcalRecHitSumEt;
  delete privateData_->dr04HcalTowerSumEt;
  delete privateData_->chargedHadronIso;
  delete privateData_->neutralHadronIso;
  delete privateData_->photonIso;
  delete privateData_->hasPixelSeed;
  delete privateData_->hasMatchedConversion;

  delete privateData_->ncand;

}


//-------------
// Methods   --
//-------------

void CmsPhotonFiller::writeCollectionToTree(edm::InputTag collectionTag,
					      const edm::Event& iEvent, const edm::EventSetup& iSetup,
					      const std::string &columnPrefix, const std::string &columnSuffix,
					      bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPhotonFiller") << "Can't get electron candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPhotonFiller") << "Track length " << collection->size() 
				   << " is too long for declared max length for tree "
				   << maxTracks_ << " and no output flag is set."
				   << " No tracks written to tuple for this event ";
      return;
    }
  
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPhotonFiller") << "Track length " << collection->size() 
				   << " is too long for declared max length for tree "
				   << maxTracks_ 
				   << ". Collection will be truncated ";
    }

    *(privateData_->ncand) = collection->size();

    // superclusters
    Handle<SuperClusterCollection> EcalBarrelSuperClusters;
    try { iEvent.getByLabel(EcalBarrelSuperClusters_, EcalBarrelSuperClusters); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsPhotonFiller") << "Can't get ECAL barrel supercluster Collection" << EcalBarrelSuperClusters_; }
    
    Handle<SuperClusterCollection> EcalEndcapSuperClusters;
    try { iEvent.getByLabel(EcalEndcapSuperClusters_, EcalEndcapSuperClusters); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsPhotonFiller") << "Can't get ECAL endcap supercluster Collection" << EcalEndcapSuperClusters_; }
    
    barrelSuperClustersSize = EcalBarrelSuperClusters->size();

    // for conversions with full vertex fit
    iEvent.getByLabel("offlineBeamSpot", bsHandle);

    try { iEvent.getByLabel(conversionsProducer_, hConversions); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsPhotonFiller") << "Can't get conversions collection " << conversionsProducer_; }

    for(int index = 0; index < (int)collection->size(); index++) {

      // fill basic kinematics
      const Candidate *cand = &(collection->at(index));
      if(saveCand_) writeCandInfo(cand,iEvent,iSetup);

      const PhotonRef photonRef = collection->refAt(index).castTo<PhotonRef>();

      if ( !(photonRef.isNull()) ) {

        // fill Cluster Adapter
        SuperClusterRef sclusRef = photonRef->superCluster();
        //        SuperClusterRef pfclusRef = photonRef->pflowSuperCluster();
        SuperClusterRef pfclusRef = photonRef->superCluster(); // placeholder
        writeEcalInfo(photonRef,iEvent,iSetup,sclusRef,pfclusRef);

      } else {
        edm::LogWarning("CmsPhotonFiller") << "Warning! The collection seems to be not made by "
                                             << "electrons, electron-specific infos will be set to default.";
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
  treeEcalInfo(columnPrefix,columnSuffix);
  
  if(dumpData) cmstree->dumpData();

  delete trkIndexName_;

}




void CmsPhotonFiller::writeEcalInfo(const PhotonRef photonRef, 
				      const edm::Event& iEvent, const edm::EventSetup& iSetup, 
                                    SuperClusterRef sclusRef, SuperClusterRef pfclusRef) {

  if(photonRef.isNonnull()) {

    // fiducial flags in ECAL
    int packed_sel = -1;
    int isEB = ( photonRef->isEB() ) ? 1 : 0;
    int isEE = ( photonRef->isEE() ) ? 1 : 0;
    int isGap = 0; // photon has not this method, but mantain the synch with ele
    int isEBEEGap = ( photonRef->isEBEEGap() ) ? 1 : 0;
    int isEBGap = ( photonRef->isEBGap() ) ? 1 : 0;
    int isEBEtaGap = ( photonRef->isEBEtaGap() ) ? 1 : 0;
    int isEBPhiGap = ( photonRef->isEBPhiGap() ) ? 1 : 0;
    int isEEGap = ( photonRef->isEEGap() ) ? 1 : 0;
    int isEEDeeGap = ( photonRef->isEEDeeGap() ) ? 1 : 0;
    int isEERingGap = ( photonRef->isEERingGap() ) ? 1 : 0;
    
    packed_sel = ( isEB << 9 ) | ( isEE << 8 ) | ( isGap << 7 ) |
      ( isEBEEGap << 6 ) | ( isEBGap << 5 ) | ( isEBEtaGap << 4 ) | ( isEBPhiGap << 3 ) |
      ( isEEGap << 2 ) | ( isEEDeeGap << 1 ) | isEERingGap;

    privateData_->fiducialFlags->push_back(packed_sel);

    //    int packed_reco; // the following available only in > 420
    //    int isStdPho = ( photonRef->isStandardPhoton() ) ? 1 : 0;
    //    int isPFPho = ( photonRef->isPFlowPhoton() ) ? 1 : 0;
    //    packed_reco = ( isStdPho << 1 ) | isPFPho;
    //    privateData_->recoFlags->push_back( packed_reco );
    privateData_->recoFlags->push_back( -1 );

    // link to the supercluster (collections are merged: barrel + endcap in this order)
    //    if ( isStdPho && sclusRef.isNonnull() ) {
    if ( sclusRef.isNonnull() ) {
      int offset = ( fabs(sclusRef->eta() ) < 1.479 ) ? 0 : barrelSuperClustersSize;
      privateData_->superClusterIndex->push_back( sclusRef.key() + offset );
    } else {
      privateData_->superClusterIndex->push_back( -1 );
    }

    //    if ( isPFPho && pfclusRef.isNonnull() ) {
//     if ( pfclusRef.isNonnull() ) { 
//       privateData_->PFsuperClusterIndex->push_back( pfclusRef.key() );
//     } else {
//       privateData_->PFsuperClusterIndex->push_back( -1 );
//     }
    privateData_->PFsuperClusterIndex->push_back( -1 );
    
    // isolations
    privateData_->hOverE->push_back(photonRef->hadronicOverEm());
    privateData_->dr03TkSumPt->push_back(photonRef->trkSumPtSolidConeDR03());
    privateData_->dr03HollowTkSumPt->push_back(photonRef->trkSumPtHollowConeDR03());
    privateData_->dr03EcalRecHitSumEt->push_back(photonRef->ecalRecHitSumEtConeDR03());
    privateData_->dr03HcalTowerSumEt->push_back(photonRef->hcalTowerSumEtConeDR03());
    privateData_->dr04TkSumPt->push_back(photonRef->trkSumPtSolidConeDR04());
    privateData_->dr04HollowTkSumPt->push_back(photonRef->trkSumPtHollowConeDR04());
    privateData_->dr04EcalRecHitSumEt->push_back(photonRef->ecalRecHitSumEtConeDR04());
    privateData_->dr04HcalTowerSumEt->push_back(photonRef->hcalTowerSumEtConeDR04());
    privateData_->chargedHadronIso->push_back(photonRef->chargedHadronIso());
    privateData_->neutralHadronIso->push_back(photonRef->neutralHadronIso());
    privateData_->photonIso->push_back(photonRef->photonIso());
    privateData_->hasPixelSeed->push_back(int(photonRef->hasPixelSeed()));

    // conversion rejection
    bool matchesConv = ConversionTools::hasMatchedConversion(*sclusRef,hConversions,bsHandle->position());
    privateData_->hasMatchedConversion->push_back(matchesConv);

  } else {
    privateData_->fiducialFlags->push_back(-1);
    privateData_->recoFlags->push_back(-1);
    privateData_->superClusterIndex->push_back( -1 );
    privateData_->PFsuperClusterIndex->push_back( -1 );

    privateData_->hOverE->push_back(-999.);
    privateData_->dr03TkSumPt->push_back(-999.);
    privateData_->dr03HollowTkSumPt->push_back(-999.);
    privateData_->dr03EcalRecHitSumEt->push_back(-999.);
    privateData_->dr03HcalTowerSumEt->push_back(-999.);
    privateData_->dr04TkSumPt->push_back(-999.);
    privateData_->dr04HollowTkSumPt->push_back(-999.);
    privateData_->dr04EcalRecHitSumEt->push_back(-999.);
    privateData_->dr04HcalTowerSumEt->push_back(-999.);
    privateData_->chargedHadronIso->push_back(-999.);
    privateData_->neutralHadronIso->push_back(-999.);
    privateData_->photonIso->push_back(-999.);
    privateData_->hasPixelSeed->push_back(0);
    privateData_->hasMatchedConversion->push_back(false);

  }

}

void CmsPhotonFiller::treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"fiducialFlags"+colSuffix).c_str(), *privateData_->fiducialFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"recoFlags"+colSuffix).c_str(), *privateData_->recoFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"superClusterIndex"+colSuffix).c_str(), *privateData_->superClusterIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFsuperClusterIndex"+colSuffix).c_str(), *privateData_->PFsuperClusterIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hOverE"+colSuffix).c_str(), *privateData_->hOverE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03TkSumPt"+colSuffix).c_str(), *privateData_->dr03TkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03HollowTkSumPt"+colSuffix).c_str(), *privateData_->dr03HollowTkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03EcalRecHitSumEt"+colSuffix).c_str(), *privateData_->dr03EcalRecHitSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03HcalTowerSumEt"+colSuffix).c_str(), *privateData_->dr03HcalTowerSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04TkSumPt"+colSuffix).c_str(), *privateData_->dr04TkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04HollowTkSumPt"+colSuffix).c_str(), *privateData_->dr04HollowTkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04EcalRecHitSumEt"+colSuffix).c_str(), *privateData_->dr04EcalRecHitSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04HcalTowerSumEt"+colSuffix).c_str(), *privateData_->dr04HcalTowerSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chargedHadronIso"+colSuffix).c_str(), *privateData_->chargedHadronIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"neutralHadronIso"+colSuffix).c_str(), *privateData_->neutralHadronIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"photonIso"+colSuffix).c_str(), *privateData_->photonIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hasPixelSeed"+colSuffix).c_str(), *privateData_->hasPixelSeed, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hasMatchedConversion"+colSuffix).c_str(),  *privateData_->hasMatchedConversion, nCandString.c_str(), 0, "Reco");

}




void CmsPhotonFillerData::initialise() {
  
  initialiseCandidate();

  fiducialFlags = new vector<int>;
  recoFlags = new vector<int>;

  superClusterIndex = new vector<int>;
  PFsuperClusterIndex = new vector<int>;

  hOverE                   = new vector<float>;
  dr03TkSumPt              = new vector<float>;
  dr03HollowTkSumPt        = new vector<float>;
  dr03EcalRecHitSumEt      = new vector<float>;
  dr03HcalTowerSumEt       = new vector<float>;
  dr04TkSumPt              = new vector<float>;
  dr04HollowTkSumPt        = new vector<float>;
  dr04EcalRecHitSumEt      = new vector<float>;
  dr04HcalTowerSumEt       = new vector<float>;
  chargedHadronIso         = new vector<float>;
  neutralHadronIso         = new vector<float>;
  photonIso                = new vector<float>;
  hasPixelSeed             = new vector<int>;
  hasMatchedConversion     = new vector<bool>;

}

void CmsPhotonFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  fiducialFlags->clear();
  recoFlags->clear();

  superClusterIndex->clear();
  PFsuperClusterIndex->clear();

  dr03TkSumPt              ->clear();
  dr03HollowTkSumPt        ->clear();
  dr03EcalRecHitSumEt      ->clear();
  dr03HcalTowerSumEt       ->clear();
  dr04TkSumPt              ->clear();
  dr04HollowTkSumPt        ->clear();
  dr04EcalRecHitSumEt      ->clear();
  dr04HcalTowerSumEt       ->clear();
  
  chargedHadronIso         ->clear();
  neutralHadronIso         ->clear();
  photonIso                ->clear();
  hasPixelSeed             ->clear();
  hasMatchedConversion     ->clear();

}
