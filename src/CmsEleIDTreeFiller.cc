#include <memory>
#include <iostream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "MyAnalysis/IsolationTools/interface/SuperClusterHitsEcalIsolation.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleTrackerIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleCalotowerIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"



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
  delete privateData_->eleStandardClass;
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
  delete privateData_->minDR03;
  delete privateData_->minDRveto03;
  delete privateData_->sumPt03;
  delete privateData_->sumPtSquared03;
  delete privateData_->sumN03;
  delete privateData_->sumPt04;
  delete privateData_->sumPt05;
  delete privateData_->sumPtPreselection;
  delete privateData_->sumHadEt04;
  delete privateData_->sumEmEt04;
  delete privateData_->sumHadEt05;
  delete privateData_->sumEmEt05;
  delete privateData_->isoFromDepsTk;
  delete privateData_->isoFromDepsEcal;
  delete privateData_->isoFromDepsHcal;
  delete privateData_->scBasedEcalSum04;
  delete privateData_->scBasedEcalSum05;
  delete privateData_->scHaloBasedEcalSum04;
  delete privateData_->scHaloBasedEcalSum05;
  delete privateData_->eleLik;
  delete privateData_->eleIdCutsLoose;
  delete privateData_->eleIdStandardCutsRobust;
  delete privateData_->eleIdStandardCutsLoose;
  delete privateData_->eleIdStandardCutsTight;
  delete privateData_->eleTip;
}


//-------------
// Methods   --
//-------------
void CmsEleIDTreeFiller::setStandalone(bool what) { standalone_=what; }


void CmsEleIDTreeFiller::writeCollectionToTree(edm::InputTag collectionTag,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {

  // used for the general tree dump
  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  privateData_->clearTrkVectors();
  
  if(collection) {

    // for cluster shape variables
    Handle< EcalRecHitCollection > EcalBarrelRecHits;
    try { iEvent.getByLabel(EcalBarrelRecHits_, EcalBarrelRecHits); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get ECAL barrel rec hits Collection" << EcalBarrelRecHits_; }
    const EcalRecHitCollection *EBRecHits = EcalBarrelRecHits.product();

    Handle< EcalRecHitCollection > EcalEndcapRecHits;
    try { iEvent.getByLabel(EcalEndcapRecHits_, EcalEndcapRecHits); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get ECAL endcap rec hits Collection" << EcalEndcapRecHits_; }
    const EcalRecHitCollection *EERecHits = EcalEndcapRecHits.product();

    
    eleIdResults_ = new eleIdContainer(5);

    iEvent.getByLabel( "egammaIDCutsLoose", (*eleIdResults_)[0] );
    iEvent.getByLabel( "egammaIDLikelihood", (*eleIdResults_)[1] );
    iEvent.getByLabel( "egammaIDStandardCutsRobust", (*eleIdResults_)[2] ); 
    iEvent.getByLabel( "egammaIDStandardCutsLoose", (*eleIdResults_)[3] );
    iEvent.getByLabel( "egammaIDStandardCutsTight", (*eleIdResults_)[4] );

    // Read the tracks and calotowers collections for isolation
    try { iEvent.getByLabel(tracksProducer_, m_tracks); }
    catch (  cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get tracks product" << tracksProducer_; }
    
    try { iEvent.getByLabel(calotowersProducer_, m_calotowers); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get calotowers product" << calotowersProducer_; }

    //read the isolation with iso deposits
    eIsoFromDepsValueMap_ = new isoContainer(3);
    iEvent.getByLabel( "eleIsoFromDepsTk", (*eIsoFromDepsValueMap_)[0] ); 
    iEvent.getByLabel( "eleIsoFromDepsEcalFromHits", (*eIsoFromDepsValueMap_)[1] ); 
    iEvent.getByLabel( "eleIsoFromDepsHcalFromHits", (*eIsoFromDepsValueMap_)[2] ); 

    for(int index = 0; index < (int)collection->size(); index++) {
	  
      const GsfElectronRef electronRef = collection->refAt(index).castTo<GsfElectronRef>();
      if ( !(electronRef.isNull()) )
	writeEleInfo(electronRef,iEvent,iSetup,EBRecHits,EERecHits);
      else edm::LogInfo("CmsEleIDTreeFiller") << "Warning! The collection seems to be not made by "
					      << "electrons, electron-specific infos will be set to default.";

    }

    delete eleIdResults_;
    delete eIsoFromDepsValueMap_;

  }

  // if used standalone, it is necessary to initialize the size of the block event by event
  if(standalone_) {
    std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
    cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  }

  treeEleInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();
}



void CmsEleIDTreeFiller::writeEleInfo(const GsfElectronRef electronRef,
				      const edm::Event& iEvent, const edm::EventSetup& iSetup,
				      const EcalRecHitCollection *EBRecHits,
				      const EcalRecHitCollection *EERecHits) {

  SuperClusterRef sclusRef = electronRef->superCluster();
  // ele corr - notcorr energy
  float myEleTrackerP   = electronRef->trackMomentumAtVtx().r();
  float myEleFullCorrE  = electronRef->energy();
  float myEleCaloCorrE  = electronRef->caloEnergy();
  float myEleNxtalCorrE = 0.;
  float myEleRawE       = 0.;
  float mySeedE         = 0.;
  float mySeedCorrE     = 0.;            
  if(&sclusRef) {
    myEleNxtalCorrE = sclusRef->energy();
    myEleRawE       = sclusRef->rawEnergy();
    mySeedE         = sclusRef->seed()->energy();
    if ((int(electronRef->classification()/10) == 3) || (int(electronRef->classification()/10) == 13) ){
      double Ebrem = 0.;  
      basicCluster_iterator bc;
      for(bc = sclusRef->clustersBegin(); bc!=sclusRef->clustersEnd(); bc++) {
	Ebrem = Ebrem +(*bc)->energy();
      }
      Ebrem = Ebrem - mySeedE;
      mySeedCorrE = myEleNxtalCorrE - Ebrem;
    }
    else {mySeedCorrE = myEleFullCorrE;}
  }
  
  // transverse impact parameter
  GsfTrackRef trRef = electronRef->gsfTrack();
  float myTip = sqrt((trRef->vertex().x())*(trRef->vertex().x()) + (trRef->vertex().y())*(trRef->vertex().y()));

  // eleID
  privateData_->eleFullCorrE     ->push_back(myEleFullCorrE);  
  privateData_->eleCaloCorrE     ->push_back(myEleCaloCorrE);  
  privateData_->eleNxtalCorrE    ->push_back(myEleNxtalCorrE); 
  privateData_->eleRawE          ->push_back(myEleRawE); 
  privateData_->eleTrackerP      ->push_back(myEleTrackerP);
  privateData_->eleClass         ->push_back(electronRef->classification());
  privateData_->eleStandardClass ->push_back(stdEleIdClassify(&(*electronRef)));
  privateData_->eleHoE           ->push_back(electronRef->hadronicOverEm());
  privateData_->eleCorrEoP       ->push_back(electronRef->eSuperClusterOverP());
  privateData_->eleNotCorrEoP    ->push_back(myEleRawE/myEleTrackerP);
  privateData_->eleCorrEoPout    ->push_back(electronRef->eSeedClusterOverPout());
  privateData_->eleNotCorrEoPout ->push_back(electronRef->eSeedClusterOverPout()*(mySeedCorrE/mySeedE));
  privateData_->eleDeltaEtaAtVtx ->push_back(electronRef->deltaEtaSuperClusterTrackAtVtx());
  privateData_->eleDeltaPhiAtVtx ->push_back(electronRef->deltaPhiSuperClusterTrackAtVtx());
  privateData_->eleDeltaEtaAtCalo->push_back(electronRef->deltaEtaSeedClusterTrackAtCalo());
  privateData_->eleDeltaPhiAtCalo->push_back(electronRef->deltaPhiSeedClusterTrackAtCalo());
  privateData_->eleTip           ->push_back(myTip);

  // results of standard electron ID sequences
  const eleIdMap & eleIdCutsLooseVal = *( (*eleIdResults_)[0] );
  const eleIdMap & eleIdLikelihoodVal = *( (*eleIdResults_)[1] );
  const eleIdMap & eleIdStandardCutsRobustVal = *( (*eleIdResults_)[2] );
  const eleIdMap & eleIdStandardCutsLooseVal = *( (*eleIdResults_)[3] );
  const eleIdMap & eleIdStandardCutsTightVal = *( (*eleIdResults_)[4] );

  privateData_->eleIdCutsLoose->push_back( eleIdCutsLooseVal[electronRef] );
  privateData_->eleLik->push_back( eleIdLikelihoodVal[electronRef] );  
  privateData_->eleIdStandardCutsRobust->push_back( eleIdStandardCutsRobustVal[electronRef] );  
  privateData_->eleIdStandardCutsLoose->push_back( eleIdStandardCutsLooseVal[electronRef] );  
  privateData_->eleIdStandardCutsTight->push_back( eleIdStandardCutsTightVal[electronRef] );  

  const isoFromDepositsMap & eIsoFromDepsTkVal = *( (*eIsoFromDepsValueMap_)[0] );
  const isoFromDepositsMap & eIsoFromDepsEcalVal = *( (*eIsoFromDepsValueMap_)[1] );
  const isoFromDepositsMap & eIsoFromDepsHcalVal = *( (*eIsoFromDepsValueMap_)[2] );

  privateData_->isoFromDepsTk->push_back( eIsoFromDepsTkVal[electronRef] );
  privateData_->isoFromDepsEcal->push_back( eIsoFromDepsEcalVal[electronRef] );
  privateData_->isoFromDepsHcal->push_back( eIsoFromDepsHcalVal[electronRef] );

  // --- isolations ---
  
  // for tracker isolation studies
  const TrackCollection tracksC = *(m_tracks.product());
  hwwEleTrackerIsolation trackIsolation(&(*electronRef), tracksC);
  trackIsolation.setExtRadius(0.3);    
  trackIsolation.setIntRadius(0.02);    
  float minDR03     = trackIsolation.minDeltaR(0.15);  
  float minDRveto03 = trackIsolation.minDeltaR_withVeto(0.15);  
  float sumPt03  = trackIsolation.getPtTracks(true,false);
  float sumPtSquared03  = trackIsolation.getPtTracks(true,true);
  float sumN03 = trackIsolation.getNTracks(1.0);

  trackIsolation.setExtRadius(0.4);    
  float sumPt04     = trackIsolation.getPtTracks();  

  trackIsolation.setExtRadius(0.5);
  float sumPt05     = trackIsolation.getPtTracks();  

  trackIsolation.setIntRadius(0.015);
  trackIsolation.setExtRadius(0.2);
  float sumPtPreselection =  trackIsolation.getPtTracks();

  privateData_->minDR03->push_back(minDR03);
  privateData_->minDRveto03->push_back(minDRveto03);
  privateData_->sumPt03->push_back(sumPt03);
  privateData_->sumPtSquared03->push_back(sumPtSquared03);
  privateData_->sumN03->push_back(sumN03);
  privateData_->sumPt04->push_back(sumPt04);
  privateData_->sumPt05->push_back(sumPt05);
  privateData_->sumPtPreselection->push_back(sumPtPreselection);

  // calo isolation - calotower based: can be used on both RECO and AOD
  const CaloTowerCollection *calotowersC = m_calotowers.product();
  hwwEleCalotowerIsolation calotowerIsolation(&(*electronRef), calotowersC);
  calotowerIsolation.setIntRadius(0.1);
  calotowerIsolation.setExtRadius(0.40);
  
  bool relative = true;
  float sumHadEt04 = calotowerIsolation.getEtHcal(relative);
  float sumEmEt04 = calotowerIsolation.getEtEcal(relative);

  calotowerIsolation.setExtRadius(0.50);
  float sumHadEt05 = calotowerIsolation.getEtHcal(relative);
  float sumEmEt05 = calotowerIsolation.getEtEcal(relative);

  privateData_->sumHadEt04->push_back(sumHadEt04);
  privateData_->sumEmEt04->push_back(sumEmEt04);
  privateData_->sumHadEt05->push_back(sumHadEt05);
  privateData_->sumEmEt05->push_back(sumEmEt05);

  // ecal isolation with SC rechits removal
  SuperClusterHitsEcalIsolation scBasedIsolation(EBRecHits,EERecHits);
  reco::SuperClusterRef sc = electronRef->get<reco::SuperClusterRef>();

  scBasedIsolation.setExtRadius(0.4);
  scBasedIsolation.excludeHalo(false);
  float scBasedEcalSum04 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scBasedEcalSum04->push_back(scBasedEcalSum04);

  scBasedIsolation.excludeHalo(true);
  float scHaloBasedEcalSum04 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scHaloBasedEcalSum04->push_back(scHaloBasedEcalSum04);

  scBasedIsolation.setExtRadius(0.5);
  scBasedIsolation.excludeHalo(false);
  float scBasedEcalSum05 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scBasedEcalSum05->push_back(scBasedEcalSum05);
  
  scBasedIsolation.excludeHalo(true);
  float scHaloBasedEcalSum05 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scHaloBasedEcalSum05->push_back(scHaloBasedEcalSum05);

}




void CmsEleIDTreeFiller::treeEleInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"eleFullCorrE"+colSuffix).c_str(), *privateData_->eleFullCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCaloCorrE"+colSuffix).c_str(), *privateData_->eleCaloCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNxtalCorrE"+colSuffix).c_str(), *privateData_->eleNxtalCorrE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleRawE"+colSuffix).c_str(), *privateData_->eleRawE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTrackerP"+colSuffix).c_str(), *privateData_->eleTrackerP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleClass"+colSuffix).c_str(), *privateData_->eleClass, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleStandardClass"+colSuffix).c_str(), *privateData_->eleStandardClass, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleHoE"+colSuffix).c_str(), *privateData_->eleHoE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCorrEoP"+colSuffix).c_str(), *privateData_->eleCorrEoP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNotCorrEoP"+colSuffix).c_str(), *privateData_->eleNotCorrEoP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleCorrEoPout"+colSuffix).c_str(), *privateData_->eleCorrEoPout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleNotCorrEoPout"+colSuffix).c_str(), *privateData_->eleNotCorrEoPout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaEtaAtVtx"+colSuffix).c_str(), *privateData_->eleDeltaEtaAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaPhiAtVtx"+colSuffix).c_str(), *privateData_->eleDeltaPhiAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaEtaAtCalo"+colSuffix).c_str(), *privateData_->eleDeltaEtaAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleDeltaPhiAtCalo"+colSuffix).c_str(), *privateData_->eleDeltaPhiAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleMinDR03"+colSuffix).c_str(), *privateData_->minDR03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleMinDRveto03"+colSuffix).c_str(), *privateData_->minDRveto03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumPt03"+colSuffix).c_str(), *privateData_->sumPt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumPtSquared03"+colSuffix).c_str(), *privateData_->sumPtSquared03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumN03"+colSuffix).c_str(), *privateData_->sumN03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumPt04"+colSuffix).c_str(), *privateData_->sumPt04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumPt05"+colSuffix).c_str(), *privateData_->sumPt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumPtPreselection"+colSuffix).c_str(), *privateData_->sumPtPreselection, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumHadEt04"+colSuffix).c_str(), *privateData_->sumHadEt04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumEmEt04"+colSuffix).c_str(), *privateData_->sumEmEt04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumHadEt05"+colSuffix).c_str(), *privateData_->sumHadEt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleSumEmEt05"+colSuffix).c_str(), *privateData_->sumEmEt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIsoFromDepsTk"+colSuffix).c_str(), *privateData_->isoFromDepsTk, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIsoFromDepsEcal"+colSuffix).c_str(), *privateData_->isoFromDepsEcal, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIsoFromDepsHcal"+colSuffix).c_str(), *privateData_->isoFromDepsHcal, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleScBasedEcalSum04"+colSuffix).c_str(), *privateData_->scBasedEcalSum04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleScBasedEcalSum05"+colSuffix).c_str(), *privateData_->scBasedEcalSum05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleScHaloBasedEcalSum04"+colSuffix).c_str(), *privateData_->scHaloBasedEcalSum04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleScHaloBasedEcalSum05"+colSuffix).c_str(), *privateData_->scHaloBasedEcalSum05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIdCutBased"+colSuffix).c_str(), *privateData_->eleIdCutsLoose, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleLikelihood"+colSuffix).c_str(), *privateData_->eleLik, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIdStandardCutsRobust"+colSuffix).c_str(), *privateData_->eleIdStandardCutsRobust, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIdStandardCutsLoose"+colSuffix).c_str(), *privateData_->eleIdStandardCutsLoose, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIdStandardCutsTight"+colSuffix).c_str(), *privateData_->eleIdStandardCutsTight, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleTip"+colSuffix).c_str(), *privateData_->eleTip, nCandString.c_str(), 0, "Reco");
}




void CmsEleIDTreeFillerData::initialise() {
  initialiseCandidate();
  eleClass                 = new vector<int>;
  eleStandardClass         = new vector<int>;
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
  minDR03                  = new vector<float>;
  minDRveto03              = new vector<float>;
  sumPt03                  = new vector<float>;
  sumPtSquared03           = new vector<float>;
  sumN03                   = new vector<float>;
  sumPt04                  = new vector<float>;
  sumPt05                  = new vector<float>;
  sumPtPreselection        = new vector<float>;
  sumHadEt04               = new vector<float>;
  sumEmEt04                = new vector<float>;
  sumHadEt05               = new vector<float>;
  sumEmEt05                = new vector<float>;
  isoFromDepsTk            = new vector<float>;
  isoFromDepsEcal          = new vector<float>;
  isoFromDepsHcal          = new vector<float>;
  scBasedEcalSum04         = new vector<float>;
  scBasedEcalSum05         = new vector<float>;
  scHaloBasedEcalSum04     = new vector<float>;
  scHaloBasedEcalSum05     = new vector<float>;
  eleIdCutsLoose           = new vector<bool>;
  eleLik                   = new vector<float>;
  eleIdStandardCutsRobust  = new vector<bool>;
  eleIdStandardCutsLoose   = new vector<bool>;
  eleIdStandardCutsTight   = new vector<bool>;
  eleTip                   = new vector<float>;
}

void CmsEleIDTreeFillerData::clearTrkVectors() {
  clearTrkVectorsCandidate();
  eleClass                 ->clear();
  eleStandardClass         ->clear();
  eleHoE                   ->clear();
  eleNotCorrEoP            ->clear();
  eleCorrEoP               ->clear();
  eleNotCorrEoPout         ->clear();
  eleCorrEoPout            ->clear();
  eleDeltaEtaAtVtx         ->clear();
  eleDeltaEtaAtCalo        ->clear();
  eleDeltaPhiAtVtx         ->clear();
  eleDeltaPhiAtCalo        ->clear();
  eleFullCorrE             ->clear();
  eleCaloCorrE             ->clear();
  eleNxtalCorrE            ->clear();
  eleRawE                  ->clear();
  eleTrackerP              ->clear();
  minDR03                  ->clear();
  minDRveto03              ->clear();
  sumPt03                  ->clear();
  sumPtSquared03           ->clear();
  sumN03                   ->clear();
  sumPt04                  ->clear();
  sumPt05                  ->clear();
  sumPtPreselection        ->clear();
  sumHadEt04               ->clear();
  sumEmEt04                ->clear();
  sumHadEt05               ->clear();
  sumEmEt05                ->clear();
  isoFromDepsTk            ->clear();
  isoFromDepsEcal          ->clear();
  isoFromDepsHcal          ->clear();
  scBasedEcalSum04         ->clear();
  scBasedEcalSum05         ->clear();
  scHaloBasedEcalSum04     ->clear();
  scHaloBasedEcalSum05     ->clear();
  eleIdCutsLoose           ->clear();
  eleLik                   ->clear();
  eleIdStandardCutsRobust  ->clear();
  eleIdStandardCutsLoose   ->clear();
  eleIdStandardCutsTight   ->clear();
  eleTip                   ->clear();
}

int CmsEleIDTreeFiller::stdEleIdClassify(const GsfElectron* electron) {
  
  double eta = electron->p4().Eta();
  double eOverP = electron->eSuperClusterOverP();
  double pin  = electron->trackMomentumAtVtx().R(); 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  
  int cat;
  if((fabs(eta)<1.479 && fBrem<0.06) || (fabs(eta)>1.479 && fBrem<0.1)) 
    cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) 
    cat=0;
  else 
    cat=2;
  
  return cat;
}
