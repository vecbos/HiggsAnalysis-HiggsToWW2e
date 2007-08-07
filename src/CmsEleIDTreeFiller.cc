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

// #include "HtoWWElectrons/HtoWWLHRecord/interface/PidLHElectronRecord.h"
// #include "HtoWWElectrons/HtoWWLHAlgo/interface/PidLHElectronAlgo.h"
// #include "HtoWWElectrons/HtoWWLHESSource/interface/PidLHElectronESSource.h"

struct CmsEleIDTreeFillerData {
  CmsTree *cmstree;
  std::string *trkIndexName;
  bool standalone;
  int maxTracks;

  vector<int>   *eleClass;
  vector<float> *eleHoE;
  vector<float> *eleNotCorrEoP,    *eleCorrEoP;
  vector<float> *eleNotCorrEoPout, *eleCorrEoPout;
  vector<float> *eleDeltaEtaAtVtx, *eleDeltaEtaAtCalo;
  vector<float> *eleDeltaPhiAtVtx, *eleDeltaPhiAtCalo;
  vector<float> *eleTrackerP;
  vector<float> *eleTrackerIso_minDR, *eleTrackerIso_sumPt;
  vector<float> *eleTrackerIso_minDR_veto;
  vector<float> *eleCaloIso_minDR,    *eleCaloIso_sumPt;
  vector<float> *eleFullCorrE,     *eleCaloCorrE;
  vector<float> *eleNxtalCorrE,    *eleRawE;
  vector<float> *eleLik;
  vector<float> *eleTip;

};

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------

CmsEleIDTreeFiller::CmsEleIDTreeFiller(CmsTree *cmstree, int maxTracks, bool noOutputIfLimitsReached ):
  privateData_(new CmsEleIDTreeFillerData)
{
  privateData_->cmstree=cmstree;
  privateData_->maxTracks=maxTracks;
  privateData_->trkIndexName = new std::string("n");
  
  privateData_->standalone = true;

  initialise();
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
void CmsEleIDTreeFiller::setStandalone(bool what) { privateData_->standalone=what; }

void CmsEleIDTreeFiller::writeCollectionToTree(const CandidateCollection *collection,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {
  this->clearTrkVectors();
  CandidateCollection::const_iterator cand;
  for(cand=collection->begin(); cand!=collection->end(); cand++) {
    const PixelMatchGsfElectron *electron = dynamic_cast<const PixelMatchGsfElectron*>(&(*cand));
    writeEleInfo(electron,iEvent,iSetup);
  }

  // if used standalone, it is necessary to initialize the size of the block event by event
  if(privateData_->standalone) {
    std::string nCandString = columnPrefix+(*privateData_->trkIndexName)+columnSuffix; 
    privateData_->cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  }

  treeEleInfo(columnPrefix,columnSuffix);
  if(dumpData) privateData_->cmstree->dumpData();
}

void CmsEleIDTreeFiller::writeCollectionToTree(const PixelMatchGsfElectronCollection *collection,
					       const edm::Event& iEvent, const edm::EventSetup& iSetup,
					       const std::string &columnPrefix, const std::string &columnSuffix,
					       bool dumpData) {
  this->clearTrkVectors();
  PixelMatchGsfElectronCollection::const_iterator electron;
  for(electron=collection->begin(); electron!=collection->end(); electron++) {
    writeEleInfo(&(*electron),iEvent,iSetup);
  }

  // if used standalone, it is necessary to initialize the size of the block event by event
  if(privateData_->standalone) {
    std::string nCandString = columnPrefix+(*privateData_->trkIndexName)+columnSuffix; 
    privateData_->cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  }

  treeEleInfo(columnPrefix,columnSuffix);
  if(dumpData) privateData_->cmstree->dumpData();
}

void CmsEleIDTreeFiller::writeEleInfo(const PixelMatchGsfElectron *electron, const edm::Event& iEvent, const edm::EventSetup& iSetup) { 

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


  SuperClusterRef sclusRef = electron->get<SuperClusterRef>();
  
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
//   GsfLHSelector eleIDlikelihood;
//   eleIDlikelihood.setup("HtoWWElectrons/HtoWWEleProducer/data/EBpdfs.root","HtoWWElectrons/HtoWWEleProducer/data/EEpdfs.root");
//   eleIDlikelihood.getInputVar(electron,iEvent,iSetup);
//   privateData_->eleLik->push_back(eleIDlikelihood.getLHRatio());

//   edm::ESHandle<PidLHElectronAlgo> likelihood;
//   iSetup.getData( likelihood );
//   privateData_->eleLik->push_back(likelihood->getLHRatio (electron,iEvent,iSetup) );
  privateData_->eleLik->push_back( -1. );

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
  float minDR_calo = -1000.
  caloIsolation.setExtRadius(0.2);    
  float sumEt_calo = caloIsolation.getEtHcal();  
  privateData_->eleCaloIso_minDR->push_back(minDR_calo);
  privateData_->eleCaloIso_sumPt->push_back(sumEt_calo);
}

void CmsEleIDTreeFiller::treeEleInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString = colPrefix+(*privateData_->trkIndexName)+colSuffix;
  privateData_->cmstree->column((colPrefix+"eleFullCorrE"+colSuffix).c_str(), *privateData_->eleFullCorrE, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleCaloCorrE"+colSuffix).c_str(), *privateData_->eleCaloCorrE, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleNxtalCorrE"+colSuffix).c_str(), *privateData_->eleNxtalCorrE, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleRawE"+colSuffix).c_str(), *privateData_->eleRawE, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleTrackerP"+colSuffix).c_str(), *privateData_->eleTrackerP, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleClass"+colSuffix).c_str(), *privateData_->eleClass, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleHoE"+colSuffix).c_str(), *privateData_->eleHoE, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleCorrEoP"+colSuffix).c_str(), *privateData_->eleCorrEoP, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleNotCorrEoP"+colSuffix).c_str(), *privateData_->eleNotCorrEoP, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleCorrEoPout"+colSuffix).c_str(), *privateData_->eleCorrEoPout, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleNotCorrEoPout"+colSuffix).c_str(), *privateData_->eleNotCorrEoPout, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleDeltaEtaAtVtx"+colSuffix).c_str(), *privateData_->eleDeltaEtaAtVtx, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleDeltaPhiAtVtx"+colSuffix).c_str(), *privateData_->eleDeltaPhiAtVtx, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleDeltaEtaAtCalo"+colSuffix).c_str(), *privateData_->eleDeltaEtaAtCalo, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleDeltaPhiAtCalo"+colSuffix).c_str(), *privateData_->eleDeltaPhiAtCalo, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleTrackerIso_minDR"+colSuffix).c_str(), *privateData_->eleTrackerIso_minDR, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleTrackerIso_minDR_veto"+colSuffix).c_str(), *privateData_->eleTrackerIso_minDR_veto, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleTrackerIso_sumPt"+colSuffix).c_str(), *privateData_->eleTrackerIso_sumPt, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleCaloIso_minDR"+colSuffix).c_str(), *privateData_->eleCaloIso_minDR, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleCaloIso_sumPt"+colSuffix).c_str(), *privateData_->eleCaloIso_sumPt, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleLikelihood"+colSuffix).c_str(), *privateData_->eleLik, nCandString.c_str(), 0, "Reco");
  privateData_->cmstree->column((colPrefix+"eleTip"+colSuffix).c_str(), *privateData_->eleTip, nCandString.c_str(), 0, "Reco");
}

void CmsEleIDTreeFiller::initialise() {
  privateData_->eleClass                 = new vector<int>;
  privateData_->eleHoE                   = new vector<float>;
  privateData_->eleNotCorrEoP            = new vector<float>;
  privateData_->eleCorrEoP               = new vector<float>;
  privateData_->eleNotCorrEoPout         = new vector<float>;
  privateData_->eleCorrEoPout            = new vector<float>;
  privateData_->eleDeltaEtaAtVtx         = new vector<float>;
  privateData_->eleDeltaEtaAtCalo        = new vector<float>;
  privateData_->eleDeltaPhiAtVtx         = new vector<float>;
  privateData_->eleDeltaPhiAtCalo        = new vector<float>;
  privateData_->eleFullCorrE             = new vector<float>;
  privateData_->eleCaloCorrE             = new vector<float>;
  privateData_->eleNxtalCorrE            = new vector<float>;
  privateData_->eleRawE                  = new vector<float>;
  privateData_->eleTrackerP              = new vector<float>;
  privateData_->eleTrackerIso_minDR      = new vector<float>;
  privateData_->eleTrackerIso_minDR_veto = new vector<float>;
  privateData_->eleTrackerIso_sumPt      = new vector<float>;
  privateData_->eleCaloIso_minDR         = new vector<float>;
  privateData_->eleCaloIso_sumPt         = new vector<float>;
  privateData_->eleLik                   = new vector<float>;
  privateData_->eleTip                   = new vector<float>;
}
void CmsEleIDTreeFiller::clearTrkVectors() {
  privateData_->eleClass            ->clear();
  privateData_->eleHoE              ->clear();
  privateData_->eleNotCorrEoP       ->clear();
  privateData_->eleCorrEoP          ->clear();
  privateData_->eleNotCorrEoPout    ->clear();
  privateData_->eleCorrEoPout       ->clear();
  privateData_->eleDeltaEtaAtVtx    ->clear();
  privateData_->eleDeltaEtaAtCalo   ->clear();
  privateData_->eleDeltaPhiAtVtx    ->clear();
  privateData_->eleDeltaPhiAtCalo   ->clear();
  privateData_->eleFullCorrE        ->clear();
  privateData_->eleCaloCorrE        ->clear();
  privateData_->eleNxtalCorrE       ->clear();
  privateData_->eleRawE             ->clear();
  privateData_->eleTrackerP         ->clear();
  privateData_->eleTrackerIso_minDR ->clear();
  privateData_->eleTrackerIso_minDR_veto ->clear();
  privateData_->eleTrackerIso_sumPt ->clear();
  privateData_->eleCaloIso_minDR    ->clear();
  privateData_->eleCaloIso_sumPt    ->clear();
  privateData_->eleLik              ->clear();
  privateData_->eleTip              ->clear();
}

