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

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"

using namespace std;
using namespace edm;
using namespace reco;

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
  delete privateData_->classification;
  delete privateData_->standardClassification;
  delete privateData_->fbrem;
  delete privateData_->nbrems;
  delete privateData_->hOverE;
  delete privateData_->eSuperClusterOverP;
  delete privateData_->eSeedOverPout;
  delete privateData_->deltaEtaAtVtx;
  delete privateData_->deltaEtaAtCalo;
  delete privateData_->deltaPhiAtVtx;
  delete privateData_->deltaPhiAtCalo;
  delete privateData_->tip;
  delete privateData_->dr03TkSumPt;
  delete privateData_->dr03EcalRecHitSumEt;
  delete privateData_->dr03HcalTowerSumEt;
  delete privateData_->dr04TkSumPt;
  delete privateData_->dr04EcalRecHitSumEt;
  delete privateData_->dr04HcalTowerSumEt;
  delete privateData_->scBasedEcalSum03;
  delete privateData_->scBasedEcalSum04;
  delete privateData_->scHaloBasedEcalSum03;
  delete privateData_->scHaloBasedEcalSum04;
  delete privateData_->eleIdCuts;
  //  delete privateData_->eleLik;
  delete privateData_->pflowMVA;
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

    
    eleIdResults_ = new eleIdContainer(6);

    iEvent.getByLabel( "eidLoose", (*eleIdResults_)[0] );
    iEvent.getByLabel( "eidTight", (*eleIdResults_)[1] );
    iEvent.getByLabel( "eidRobustLoose", (*eleIdResults_)[2] );
    iEvent.getByLabel( "eidRobustTight", (*eleIdResults_)[3] );
    iEvent.getByLabel( "eidRobustHighEnergy", (*eleIdResults_)[4] );
    //    iEvent.getByLabel( "egammaIDLikelihood", (*eleIdResults_)[5] );

    // Read the tracks and calotowers collections for isolation
    try { iEvent.getByLabel(tracksProducer_, m_tracks); }
    catch (  cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get tracks product" << tracksProducer_; }
    
    try { iEvent.getByLabel(calotowersProducer_, m_calotowers); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get calotowers product" << calotowersProducer_; }

    for(int index = 0; index < (int)collection->size(); index++) {
	  
      const GsfElectronRef electronRef = collection->refAt(index).castTo<GsfElectronRef>();
      if ( !(electronRef.isNull()) )
	writeEleInfo(electronRef,iEvent,iSetup,EBRecHits,EERecHits);
      else edm::LogWarning("CmsEleIDTreeFiller") << "Warning! The collection seems to be not made by "
                                                 << "electrons, electron-specific infos will be set to default.";

    }

    delete eleIdResults_;

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
  float mySeedE         = 0.;
  if(sclusRef.isNonnull()) mySeedE = sclusRef->seed()->energy();
  
  // transverse impact parameter
  GsfTrackRef trRef = electronRef->gsfTrack();
  float myTip = sqrt((trRef->vertex().x())*(trRef->vertex().x()) + (trRef->vertex().y())*(trRef->vertex().y()));

  // eleID
  privateData_->classification->push_back(electronRef->classification());
  privateData_->standardClassification->push_back(stdEleIdClassify(&(*electronRef)));
  privateData_->fbrem->push_back(electronRef->fbrem());
  privateData_->nbrems->push_back(electronRef->numberOfBrems());
  privateData_->hOverE->push_back(electronRef->hadronicOverEm());
  privateData_->eSuperClusterOverP->push_back(electronRef->eSuperClusterOverP());
  privateData_->eSeedOverPout->push_back(electronRef->eSeedClusterOverPout());
  privateData_->deltaEtaAtVtx->push_back(electronRef->deltaEtaSuperClusterTrackAtVtx());
  privateData_->deltaPhiAtVtx->push_back(electronRef->deltaPhiSuperClusterTrackAtVtx());
  privateData_->deltaEtaAtCalo->push_back(electronRef->deltaEtaSeedClusterTrackAtCalo());
  privateData_->deltaPhiAtCalo->push_back(electronRef->deltaPhiSeedClusterTrackAtCalo());
  privateData_->tip->push_back(myTip);

  // results of standard electron ID sequences
  const eleIdMap & eleIdLooseVal = *( (*eleIdResults_)[0] );
  const eleIdMap & eleIdTightVal = *( (*eleIdResults_)[1] );
  const eleIdMap & eleIdRobustLooseVal = *( (*eleIdResults_)[2] );
  const eleIdMap & eleIdRobustTightVal = *( (*eleIdResults_)[3] );
  const eleIdMap & eleIdRobustHighEnergyVal = *( (*eleIdResults_)[4] );
  //  const eleIdMap & eleIdLikelihoodVal = *( (*eleIdResults_)[5] );

  int packed_sel = -1;
  int eleIdLoose = ( eleIdLooseVal[electronRef] ) ? 1 : 0;
  int eleIdTight = ( eleIdTightVal[electronRef] ) ? 1 : 0;
  int eleIdRobustLoose = ( eleIdRobustLooseVal[electronRef] ) ? 1 : 0;
  int eleIdRobustTight = ( eleIdRobustTightVal[electronRef] ) ? 1 : 0;
  int eleIdRobustHighEnergy = ( eleIdRobustHighEnergyVal[electronRef] ) ? 1 : 0;

  packed_sel = ( eleIdLoose << 4 ) | ( eleIdTight << 3 ) |
    ( eleIdRobustLoose << 2 ) | ( eleIdRobustTight << 1 ) | eleIdRobustHighEnergy;

  privateData_->eleIdCuts->push_back( packed_sel );
  //  privateData_->eleLik->push_back( eleIdLikelihoodVal[electronRef] );  
  privateData_->pflowMVA->push_back( electronRef->mva() );

  // --- isolations ---

  privateData_->dr03TkSumPt->push_back( electronRef->dr03TkSumPt() );
  privateData_->dr03EcalRecHitSumEt->push_back( electronRef->dr03EcalRecHitSumEt() );
  privateData_->dr03HcalTowerSumEt->push_back( electronRef->dr03HcalTowerSumEt() );

  privateData_->dr04TkSumPt->push_back( electronRef->dr04TkSumPt() );
  privateData_->dr04EcalRecHitSumEt->push_back( electronRef->dr04EcalRecHitSumEt() );
  privateData_->dr04HcalTowerSumEt->push_back( electronRef->dr04HcalTowerSumEt() );

  // ecal isolation with SC rechits removal
  SuperClusterHitsEcalIsolation scBasedIsolation(EBRecHits,EERecHits);
  reco::SuperClusterRef sc = electronRef->get<reco::SuperClusterRef>();

  scBasedIsolation.setExtRadius(0.3);
  scBasedIsolation.excludeHalo(false);
  float scBasedEcalSum03 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scBasedEcalSum03->push_back(scBasedEcalSum03);
  
  scBasedIsolation.excludeHalo(true);
  float scHaloBasedEcalSum03 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scHaloBasedEcalSum03->push_back(scHaloBasedEcalSum03);

  scBasedIsolation.setExtRadius(0.4);
  scBasedIsolation.excludeHalo(false);
  float scBasedEcalSum04 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scBasedEcalSum04->push_back(scBasedEcalSum04);

  scBasedIsolation.excludeHalo(true);
  float scHaloBasedEcalSum04 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scHaloBasedEcalSum04->push_back(scHaloBasedEcalSum04);

}




void CmsEleIDTreeFiller::treeEleInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"classification"+colSuffix).c_str(), *privateData_->classification, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"standardClassification"+colSuffix).c_str(), *privateData_->standardClassification, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hOverE"+colSuffix).c_str(), *privateData_->hOverE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eSuperClusterOverP"+colSuffix).c_str(), *privateData_->eSuperClusterOverP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eSeedOverPout"+colSuffix).c_str(), *privateData_->eSeedOverPout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaEtaAtVtx"+colSuffix).c_str(), *privateData_->deltaEtaAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiAtVtx"+colSuffix).c_str(), *privateData_->deltaPhiAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaEtaAtCalo"+colSuffix).c_str(), *privateData_->deltaEtaAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiAtCalo"+colSuffix).c_str(), *privateData_->deltaPhiAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"tip"+colSuffix).c_str(), *privateData_->tip, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03TkSumPt"+colSuffix).c_str(), *privateData_->dr03TkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03EcalRecHitSumEt"+colSuffix).c_str(), *privateData_->dr03EcalRecHitSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03HcalTowerSumEt"+colSuffix).c_str(), *privateData_->dr03HcalTowerSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04TkSumPt"+colSuffix).c_str(), *privateData_->dr04TkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04EcalRecHitSumEt"+colSuffix).c_str(), *privateData_->dr04EcalRecHitSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04HcalTowerSumEt"+colSuffix).c_str(), *privateData_->dr04HcalTowerSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"scBasedEcalSum03"+colSuffix).c_str(), *privateData_->scBasedEcalSum03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"scBasedEcalSum04"+colSuffix).c_str(), *privateData_->scBasedEcalSum04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"scHaloBasedEcalSum03"+colSuffix).c_str(), *privateData_->scHaloBasedEcalSum03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"scHaloBasedEcalSum04"+colSuffix).c_str(), *privateData_->scHaloBasedEcalSum04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIdCuts"+colSuffix).c_str(), *privateData_->eleIdCuts, nCandString.c_str(), 0, "Reco");
  //  cmstree->column((colPrefix+"eleIdLikelihood"+colSuffix).c_str(), *privateData_->eleLik, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pflowMVA"+colSuffix).c_str(), *privateData_->pflowMVA, nCandString.c_str(), 0, "Reco");
}




void CmsEleIDTreeFillerData::initialise() {
  initialiseCandidate();
  classification           = new vector<int>;
  standardClassification   = new vector<int>;
  fbrem                    = new vector<float>;
  nbrems                   = new vector<int>;
  hOverE                   = new vector<float>;
  eSuperClusterOverP       = new vector<float>;
  eSeedOverPout            = new vector<float>;
  deltaEtaAtVtx            = new vector<float>;
  deltaEtaAtCalo           = new vector<float>;
  deltaPhiAtVtx            = new vector<float>;
  deltaPhiAtCalo           = new vector<float>;
  tip                      = new vector<float>;
  dr03TkSumPt              = new vector<float>;
  dr03EcalRecHitSumEt      = new vector<float>;
  dr03HcalTowerSumEt       = new vector<float>;
  dr04TkSumPt              = new vector<float>;
  dr04EcalRecHitSumEt      = new vector<float>;
  dr04HcalTowerSumEt       = new vector<float>;
  scBasedEcalSum03         = new vector<float>;
  scBasedEcalSum04         = new vector<float>;
  scHaloBasedEcalSum03     = new vector<float>;
  scHaloBasedEcalSum04     = new vector<float>;
  eleIdCuts                = new vector<int>;
  //  eleLik                   = new vector<float>;
  pflowMVA                 = new vector<bool>;
}

void CmsEleIDTreeFillerData::clearTrkVectors() {
  clearTrkVectorsCandidate();
  classification           ->clear();
  standardClassification   ->clear();
  fbrem                    ->clear();
  nbrems                   ->clear();
  hOverE                   ->clear();
  eSuperClusterOverP       ->clear();
  eSeedOverPout            ->clear();
  deltaEtaAtVtx            ->clear();
  deltaEtaAtCalo           ->clear();
  deltaPhiAtVtx            ->clear();
  deltaPhiAtCalo           ->clear();
  tip                      ->clear();
  dr03TkSumPt              ->clear();
  dr03EcalRecHitSumEt      ->clear();
  dr03HcalTowerSumEt       ->clear();
  dr04TkSumPt              ->clear();
  dr04EcalRecHitSumEt      ->clear();
  dr04HcalTowerSumEt       ->clear();
  scBasedEcalSum03         ->clear();
  scBasedEcalSum04         ->clear();
  scHaloBasedEcalSum03     ->clear();
  scHaloBasedEcalSum04     ->clear();
  eleIdCuts                ->clear();
  //  eleLik                   ->clear();
  pflowMVA                 ->clear();
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
