#include <memory>
#include <iostream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleTrackerIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleCaloIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"

#include "EgammaAnalysis/ElectronIDAlgos/interface/ElectronLikelihood.h"
#include "EgammaAnalysis/ElectronIDESSources/interface/ElectronLikelihoodESSource.h"


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------

CmsEleIDTreeFiller::CmsEleIDTreeFiller(CmsTree *cmsTree, int maxTracks, bool noOutputIfLimitsReached ):
  CmsCandidateFiller(cmsTree,maxTracks,500,noOutputIfLimitsReached),
  privateData_(new CmsEleIDTreeFillerData)
{
  cmstree=cmsTree;
  maxTracks_=maxTracks;
  trkIndexName_ = new std::string("n");
  
  standalone_ = true;

  privateData_->initialise();
}

//--------------
// Destructor --
//--------------
CmsEleIDTreeFiller::~CmsEleIDTreeFiller() {
  // delete here the vector ptr's
  delete privateData_->eleClass;
  delete privateData_->eleHoE;
  delete privateData_->eleNotCorrEoP;
  delete privateData_->eleCorrEoP;
  delete privateData_->eleNotCorrEoPout;
  delete privateData_->eleCorrEoPout;
  delete privateData_->eleDeltaEtaAtVtx;
  delete privateData_->eleDeltaEtaAtCalo;
  delete privateData_->eleDeltaPhiAtVtx;
  delete privateData_->eleDeltaPhiAtCalo;
  delete privateData_->eleFullCorrE;
  delete privateData_->eleCaloCorrE;
  delete privateData_->eleNxtalCorrE;
  delete privateData_->eleRawE;
  delete privateData_->eleTrackerP;
  delete privateData_->eleTrackerIso_minDR;
  delete privateData_->eleTrackerIso_minDR_veto;
  delete privateData_->eleTrackerIso_sumPt;
  delete privateData_->eleCaloIso_minDR;
  delete privateData_->eleCaloIso_sumPt;
  delete privateData_->eleLik;
  delete privateData_->eleTip;
}


//-------------
// Methods   --
//-------------
void CmsEleIDTreeFiller::setStandalone(bool what) { standalone_=what; }


void CmsEleIDTreeFiller::writeCollectionToTree(const CandidateCollection *collection,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {
  privateData_->clearTrkVectors();
  
  if(collection) {
    CandidateCollection::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      if ( cand->hasMasterClone() ) {   
	CandidateBaseRef master = cand->masterClone();
	PixelMatchGsfElectronRef electronRef = cand->masterClone().castTo<PixelMatchGsfElectronRef>();
	const PixelMatchGsfElectron &electron = *electronRef;
	writeEleInfo(&electron,iEvent,iSetup);
      }
    }
  }

  // if used standalone, it is necessary to initialize the size of the block event by event
  if(standalone_) {
    std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
    cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  }

  treeEleInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();
}





void CmsEleIDTreeFiller::writeCollectionToTree(const PixelMatchGsfElectronCollection *collection,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {
  privateData_->clearTrkVectors();

  if(collection) {
    PixelMatchGsfElectronCollection::const_iterator electron;
    for(electron=collection->begin(); electron!=collection->end(); electron++) {
      writeEleInfo(&(*electron),iEvent,iSetup);
    }
  }

  // if used standalone, it is necessary to initialize the size of the block event by event
  if(standalone_) {
    std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
    cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  }

  treeEleInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();
}




void CmsEleIDTreeFiller::writeEleInfo(const PixelMatchGsfElectron *electron, 
				      const edm::Event& iEvent, const edm::EventSetup& iSetup) { 

  // --------------------------------------
  // collections needed for isolation
  //
  // get reconstructed tracks
  edm::Handle<TrackCollection> tracks;
  try { iEvent.getByLabel("ctfWithMaterialTracks","",tracks); }
  catch ( cms::Exception& ex ) { cout << "Can't get collection: " << "ctfWithMaterialTracks" << endl; }

  /// get hcal cells
  edm::Handle<HBHERecHitCollection> hcalrhits;
  try { iEvent.getByLabel("hbhereco", hcalrhits); }
  catch ( cms::Exception& ex ) { cout << "Can't get collection: " << "hbhereco" << endl; }

  // taking the calo geometry
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<IdealGeometryRecord>().get(pG);
  caloGeo = pG.product();
  //--------------------------------------

  SuperClusterRef sclusRef = electron->superCluster();
  // ele corr - notcorr energy
  float myEleTrackerP   = electron->trackMomentumAtVtx().r();
  float myEleFullCorrE  = electron->energy();
  float myEleCaloCorrE  = electron->caloEnergy();
  float myEleNxtalCorrE = 0.;
  float myEleRawE       = 0.;
  float mySeedE         = 0.;
  float mySeedCorrE     = 0.;            
  if(&sclusRef) {
    myEleNxtalCorrE = sclusRef->energy();
    myEleRawE       = sclusRef->rawEnergy();
    mySeedE         = sclusRef->seed()->energy();
    if ((int(electron->classification()/10) == 3) || (int(electron->classification()/10) == 13) ){
      double Ebrem = 0.;  
      basicCluster_iterator bc;
      for(bc = sclusRef->clustersBegin(); bc!=sclusRef->clustersEnd(); bc++) {
	Ebrem = Ebrem + (*bc)->energy();
      }
      Ebrem = Ebrem - mySeedE;
      mySeedCorrE = myEleNxtalCorrE - Ebrem;
    }
    else {mySeedCorrE = myEleFullCorrE;}
  }

  // transverse impact parameter
  GsfTrackRef trRef = electron->gsfTrack();
  float myTip = sqrt((trRef->vertex().x())*(trRef->vertex().x()) + (trRef->vertex().y())*(trRef->vertex().y()));

  // eleID
  privateData_->eleFullCorrE     ->push_back(myEleFullCorrE);  
  privateData_->eleCaloCorrE     ->push_back(myEleCaloCorrE);  
  privateData_->eleNxtalCorrE    ->push_back(myEleNxtalCorrE); 
  privateData_->eleRawE          ->push_back(myEleRawE); 
  privateData_->eleTrackerP      ->push_back(myEleTrackerP);
  privateData_->eleClass         ->push_back(electron->classification());
  privateData_->eleHoE           ->push_back(electron->hadronicOverEm());
  privateData_->eleCorrEoP       ->push_back(electron->eSuperClusterOverP());
  privateData_->eleNotCorrEoP    ->push_back(myEleRawE/myEleTrackerP);
  privateData_->eleCorrEoPout    ->push_back(electron->eSeedClusterOverPout());
  privateData_->eleNotCorrEoPout ->push_back(electron->eSeedClusterOverPout()*(mySeedCorrE/mySeedE));
  privateData_->eleDeltaEtaAtVtx ->push_back(electron->deltaEtaSuperClusterTrackAtVtx());
  privateData_->eleDeltaPhiAtVtx ->push_back(electron->deltaPhiSuperClusterTrackAtVtx());
  privateData_->eleDeltaEtaAtCalo->push_back(electron->deltaEtaSeedClusterTrackAtCalo());
  privateData_->eleDeltaPhiAtCalo->push_back(electron->deltaPhiSeedClusterTrackAtCalo());
  privateData_->eleTip           ->push_back(myTip);

  // electron likelihood
  // FIX: reorganize this, cluster shape collections have to be passed
  // in the cfg file as parameters
  bool hasBarrel=true;
  bool hasEndcap=true;

  Handle<BasicClusterShapeAssociationCollection> barrelClShpHandle;
  try { iEvent.getByLabel("hybridSuperClusters","hybridShapeAssoc", barrelClShpHandle); }
  catch ( cms::Exception& ex ) { LogWarning("CmsTreeFiller") << "Can't get ECAL barrel Cluster Shape Collection"; }
  const reco::BasicClusterShapeAssociationCollection& barrelClShpMap = *barrelClShpHandle;

  Handle<BasicClusterShapeAssociationCollection> endcapClShpHandle;
  try { iEvent.getByLabel("islandBasicClusters","islandEndcapShapeAssoc", endcapClShpHandle); }
  catch ( cms::Exception& ex ) { LogWarning("CmsTreeFiller") << "Can't get ECAL endcap Cluster Shape Collection"; }
  const reco::BasicClusterShapeAssociationCollection& endcapClShpMap = *endcapClShpHandle;

  reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
  seedShpItr = barrelClShpMap.find(sclusRef->seed());
  if(seedShpItr==barrelClShpMap.end()) {
    hasBarrel=false;
    seedShpItr=endcapClShpMap.find(sclusRef->seed());
    if(seedShpItr==endcapClShpMap.end()) hasEndcap=false;
  }
  if(hasBarrel || hasEndcap) {
    const ClusterShapeRef& sClShape = seedShpItr->val;  
    edm::ESHandle<ElectronLikelihood> likelihood;
    iSetup.getData( likelihood );
    privateData_->eleLik->push_back( likelihood->result(*electron,*sClShape) );

  }
  else {
    edm::LogWarning("CmsEleIDTreeFiller") << "cannot find cluster shapes in ECAL barrel or endcap "
					  << " setting likelihood value to -1";
    privateData_->eleLik->push_back( -1.);
  }

  // tracker isolation
  const TrackCollection tracksC = *(tracks.product());
  hwwEleTrackerIsolation trackIsolation(electron, tracksC);
  trackIsolation.setExtRadius(0.2);    
  trackIsolation.setIntRadius(0.015);    
  float minDR_tracker     = trackIsolation.minDeltaR(0.15);  
  float minDRveto_tracker = trackIsolation.minDeltaR_withVeto(0.15);  
  float sumPt_tracker     = trackIsolation.getPtTracks();  
  privateData_->eleTrackerIso_minDR->push_back(minDR_tracker);
  privateData_->eleTrackerIso_minDR_veto->push_back(minDRveto_tracker);
  privateData_->eleTrackerIso_sumPt->push_back(sumPt_tracker);

  // calo isolation
  const HBHERecHitCollection hcalRecHits = *(hcalrhits.product());
  hwwEleCaloIsolation caloIsolation(electron, hcalRecHits, caloGeo);
  //float minDR_calo = caloIsolation.minDeltaR(0.15);  
  float minDR_calo = -1000.;
  caloIsolation.setExtRadius(0.2);    
  float sumEt_calo = caloIsolation.getEtHcal();  
  privateData_->eleCaloIso_minDR->push_back(minDR_calo);
  privateData_->eleCaloIso_sumPt->push_back(sumEt_calo);

}




void CmsEleIDTreeFiller::treeEleInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"eleFullCorrE"+colSuffix).c_str(), *privateData_->eleFullCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCaloCorrE"+colSuffix).c_str(), *privateData_->eleCaloCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNxtalCorrE"+colSuffix).c_str(), *privateData_->eleNxtalCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleRawE"+colSuffix).c_str(), *privateData_->eleRawE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackerP"+colSuffix).c_str(), *privateData_->eleTrackerP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleClass"+colSuffix).c_str(), *privateData_->eleClass, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleHoE"+colSuffix).c_str(), *privateData_->eleHoE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCorrEoP"+colSuffix).c_str(), *privateData_->eleCorrEoP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNotCorrEoP"+colSuffix).c_str(), *privateData_->eleNotCorrEoP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCorrEoPout"+colSuffix).c_str(), *privateData_->eleCorrEoPout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNotCorrEoPout"+colSuffix).c_str(), *privateData_->eleNotCorrEoPout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaEtaAtVtx"+colSuffix).c_str(), *privateData_->eleDeltaEtaAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaPhiAtVtx"+colSuffix).c_str(), *privateData_->eleDeltaPhiAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaEtaAtCalo"+colSuffix).c_str(), *privateData_->eleDeltaEtaAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaPhiAtCalo"+colSuffix).c_str(), *privateData_->eleDeltaPhiAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackerIso_minDR"+colSuffix).c_str(), *privateData_->eleTrackerIso_minDR, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackerIso_minDR_veto"+colSuffix).c_str(), *privateData_->eleTrackerIso_minDR_veto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackerIso_sumPt"+colSuffix).c_str(), *privateData_->eleTrackerIso_sumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCaloIso_minDR"+colSuffix).c_str(), *privateData_->eleCaloIso_minDR, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCaloIso_sumPt"+colSuffix).c_str(), *privateData_->eleCaloIso_sumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleLikelihood"+colSuffix).c_str(), *privateData_->eleLik, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTip"+colSuffix).c_str(), *privateData_->eleTip, nCandString.c_str(), 0, "Reco");
}




void CmsEleIDTreeFillerData::initialise() {
  initialiseCandidate();
  eleClass                 = new vector<int>;
  eleHoE                   = new vector<float>;
  eleNotCorrEoP            = new vector<float>;
  eleCorrEoP               = new vector<float>;
  eleNotCorrEoPout         = new vector<float>;
  eleCorrEoPout            = new vector<float>;
  eleDeltaEtaAtVtx         = new vector<float>;
  eleDeltaEtaAtCalo        = new vector<float>;
  eleDeltaPhiAtVtx         = new vector<float>;
  eleDeltaPhiAtCalo        = new vector<float>;
  eleFullCorrE             = new vector<float>;
  eleCaloCorrE             = new vector<float>;
  eleNxtalCorrE            = new vector<float>;
  eleRawE                  = new vector<float>;
  eleTrackerP              = new vector<float>;
  eleTrackerIso_minDR      = new vector<float>;
  eleTrackerIso_minDR_veto = new vector<float>;
  eleTrackerIso_sumPt      = new vector<float>;
  eleCaloIso_minDR         = new vector<float>;
  eleCaloIso_sumPt         = new vector<float>;
  eleLik                   = new vector<float>;
  eleTip                   = new vector<float>;
}




void CmsEleIDTreeFillerData::clearTrkVectors() {
  clearTrkVectorsCandidate();
  eleClass            ->clear();
  eleHoE              ->clear();
  eleNotCorrEoP       ->clear();
  eleCorrEoP          ->clear();
  eleNotCorrEoPout    ->clear();
  eleCorrEoPout       ->clear();
  eleDeltaEtaAtVtx    ->clear();
  eleDeltaEtaAtCalo   ->clear();
  eleDeltaPhiAtVtx    ->clear();
  eleDeltaPhiAtCalo   ->clear();
  eleFullCorrE        ->clear();
  eleCaloCorrE        ->clear();
  eleNxtalCorrE       ->clear();
  eleRawE             ->clear();
  eleTrackerP         ->clear();
  eleTrackerIso_minDR ->clear();
  eleTrackerIso_minDR_veto ->clear();
  eleTrackerIso_sumPt ->clear();
  eleCaloIso_minDR    ->clear();
  eleCaloIso_sumPt    ->clear();
  eleLik              ->clear();
  eleTip              ->clear();
}

