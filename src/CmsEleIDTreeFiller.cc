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
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h" 

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
  delete privateData_->ambiguousGsfTracksSize;
  delete privateData_->hOverE;
  delete privateData_->eSuperClusterOverP;
  delete privateData_->eSeedOverPout;
  delete privateData_->eEleClusterOverPout;
  delete privateData_->deltaEtaEleClusterTrackAtCalo;
  delete privateData_->deltaPhiEleClusterTrackAtCalo;
  delete privateData_->deltaEtaAtVtx;
  delete privateData_->deltaEtaAtCalo;
  delete privateData_->deltaPhiAtVtx;
  delete privateData_->deltaPhiAtCalo;
  delete privateData_->dr03TkSumPt;
  delete privateData_->dr03EcalRecHitSumEt;
  delete privateData_->dr03HcalTowerSumEt;
  delete privateData_->dr04TkSumPt;
  delete privateData_->dr04EcalRecHitSumEt;
  delete privateData_->dr04HcalTowerSumEt;
  delete privateData_->scBasedEcalSum03;
  delete privateData_->scBasedEcalSum04;
  delete privateData_->dr03HcalTowerSumEtFullCone;
  delete privateData_->dr04HcalTowerSumEtFullCone;
  delete privateData_->eleLik;
  delete privateData_->pflowMVA;
  delete privateData_->mvaidtrig;
  delete privateData_->mvaidisotrig;
  delete privateData_->mvaidnontrig;
  delete privateData_->pfCombinedIso;
  delete privateData_->pfCandChargedIso01;
  delete privateData_->pfCandNeutralIso01;
  delete privateData_->pfCandPhotonIso01;
  delete privateData_->pfCandChargedIso02;
  delete privateData_->pfCandNeutralIso02;
  delete privateData_->pfCandPhotonIso02;
  delete privateData_->pfCandChargedIso03;
  delete privateData_->pfCandNeutralIso03;
  delete privateData_->pfCandPhotonIso03;
  delete privateData_->pfCandChargedIso04;
  delete privateData_->pfCandNeutralIso04;
  delete privateData_->pfCandPhotonIso04;
  delete privateData_->pfCandChargedIso05;
  delete privateData_->pfCandNeutralIso05;
  delete privateData_->pfCandPhotonIso05;
  delete privateData_->pfCandChargedIso06;
  delete privateData_->pfCandNeutralIso06;
  delete privateData_->pfCandPhotonIso06;
  delete privateData_->pfCandChargedIso07;
  delete privateData_->pfCandNeutralIso07;
  delete privateData_->pfCandPhotonIso07;
  delete privateData_->pfCandChargedDirIso04;
  delete privateData_->pfCandNeutralDirIso04;
  delete privateData_->pfCandPhotonDirIso04;

  delete privateData_;
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

    
    eleIdResults_ = new eleIdContainer(1);

    iEvent.getByLabel( "egammaIDLikelihood", (*eleIdResults_)[0] );

    eIsoFromPFCandsValueMap_ = new isoContainer(25);
    iEvent.getByLabel( "electronCombinedPFIsoMapProducer", (*eIsoFromPFCandsValueMap_)[0] ); 
    iEvent.getByLabel( "electronPFIsoChHad01", (*eIsoFromPFCandsValueMap_)[1] );
    iEvent.getByLabel( "electronPFIsoNHad01", (*eIsoFromPFCandsValueMap_)[2] );
    iEvent.getByLabel( "electronPFIsoPhoton01", (*eIsoFromPFCandsValueMap_)[3] );
    iEvent.getByLabel( "electronPFIsoChHad02", (*eIsoFromPFCandsValueMap_)[4] );
    iEvent.getByLabel( "electronPFIsoNHad02", (*eIsoFromPFCandsValueMap_)[5] );
    iEvent.getByLabel( "electronPFIsoPhoton02", (*eIsoFromPFCandsValueMap_)[6] );
    iEvent.getByLabel( "electronPFIsoChHad03", (*eIsoFromPFCandsValueMap_)[7] );
    iEvent.getByLabel( "electronPFIsoNHad03", (*eIsoFromPFCandsValueMap_)[8] );
    iEvent.getByLabel( "electronPFIsoPhoton03", (*eIsoFromPFCandsValueMap_)[9] );
    iEvent.getByLabel( "electronPFIsoChHad04", (*eIsoFromPFCandsValueMap_)[10] );
    iEvent.getByLabel( "electronPFIsoNHad04", (*eIsoFromPFCandsValueMap_)[11] );
    iEvent.getByLabel( "electronPFIsoPhoton04", (*eIsoFromPFCandsValueMap_)[12] );
    iEvent.getByLabel( "electronPFIsoChHad05", (*eIsoFromPFCandsValueMap_)[13] );
    iEvent.getByLabel( "electronPFIsoNHad05", (*eIsoFromPFCandsValueMap_)[14] );
    iEvent.getByLabel( "electronPFIsoPhoton05", (*eIsoFromPFCandsValueMap_)[15] );
    iEvent.getByLabel( "electronPFIsoChHad06", (*eIsoFromPFCandsValueMap_)[16] );
    iEvent.getByLabel( "electronPFIsoNHad06", (*eIsoFromPFCandsValueMap_)[17] );
    iEvent.getByLabel( "electronPFIsoPhoton06", (*eIsoFromPFCandsValueMap_)[18] );
    iEvent.getByLabel( "electronPFIsoChHad07", (*eIsoFromPFCandsValueMap_)[19] );
    iEvent.getByLabel( "electronPFIsoNHad07", (*eIsoFromPFCandsValueMap_)[20] );
    iEvent.getByLabel( "electronPFIsoPhoton07", (*eIsoFromPFCandsValueMap_)[21] );
    iEvent.getByLabel( "electronDirPFIsoChHad04", (*eIsoFromPFCandsValueMap_)[22] );
    iEvent.getByLabel( "electronDirPFIsoNHad04", (*eIsoFromPFCandsValueMap_)[23] );
    iEvent.getByLabel( "electronDirPFIsoPhoton04", (*eIsoFromPFCandsValueMap_)[24] );

    // Read the tracks and calotowers collections for isolation
    try { iEvent.getByLabel(tracksProducer_, m_tracks); }
    catch (  cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get tracks product" << tracksProducer_; }
    
    try { iEvent.getByLabel(calotowersProducer_, m_calotowers); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsEleIDTreeFiller") << "Can't get calotowers product" << calotowersProducer_; }

    // vertices and transient track builder needed for ele ID MVA
    try { iEvent.getByLabel(m_vxtCollectionTag, primaryVertex); }
    catch(cms::Exception& ex ) {edm::LogWarning("CmsVertexFiller") << "Can't get candidate collection: " << m_vxtCollectionTag; }

    try { iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder_); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsTrackFiller") << "Can't get TransientTrackBuilder from Event Setup."; }

    // the full pf candidates to make the combined ID+ISO MVA
    try { iEvent.getByLabel(m_pfcandCollectionTag, pfcandidates_); }
    catch(cms::Exception& ex ) {edm::LogWarning("CmsVertexFiller") << "Can't get candidate collection: " << m_pfcandCollectionTag; }
    
    iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"),rho_);

    // HCAL isolation in the full cone (to be used w/o H/E)
    float egHcalIsoConeSizeOutSmall=0.3, egHcalIsoConeSizeOutLarge=0.4;
    float egHcalIsoConeSizeIn=0.0, egHcalIsoPtMin=0.0;
    int egHcalDepth1=1, egHcalDepth2=2;
    hadDepth1Isolation03_ = new EgammaTowerIsolation(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,m_calotowers.product()) ;
    hadDepth2Isolation03_ = new EgammaTowerIsolation(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,m_calotowers.product()) ;
    hadDepth1Isolation04_ = new EgammaTowerIsolation(egHcalIsoConeSizeOutLarge,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,m_calotowers.product()) ;
    hadDepth2Isolation04_ = new EgammaTowerIsolation(egHcalIsoConeSizeOutLarge,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,m_calotowers.product()) ;

    for(int index = 0; index < (int)collection->size(); index++) {
	  
      const GsfElectronRef electronRef = collection->refAt(index).castTo<GsfElectronRef>();
      if ( !(electronRef.isNull()) )
	writeEleInfo(electronRef,iEvent,iSetup,EBRecHits,EERecHits);
      else edm::LogWarning("CmsEleIDTreeFiller") << "Warning! The collection seems to be not made by "
                                                 << "electrons, electron-specific infos will be set to default.";

    }

    delete eleIdResults_;
    delete eIsoFromPFCandsValueMap_;

  }

  // if used standalone, it is necessary to initialize the size of the block event by event
  if(standalone_) {
    std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
    cmstree->column(nCandString.c_str(),collection->size(),0,"Reco");
  }

  treeEleInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();
  
  delete hadDepth1Isolation03_;
  delete hadDepth2Isolation03_;
  delete hadDepth1Isolation04_;
  delete hadDepth2Isolation04_;
  delete trkIndexName_;
}



void CmsEleIDTreeFiller::writeEleInfo(const GsfElectronRef electronRef,
				      const edm::Event& iEvent, const edm::EventSetup& iSetup,
				      const EcalRecHitCollection *EBRecHits,
				      const EcalRecHitCollection *EERecHits) {

  SuperClusterRef sclusRef = electronRef->superCluster();
  
  // eleID
  privateData_->classification->push_back(electronRef->classification());
  privateData_->standardClassification->push_back(stdEleIdClassify(&(*electronRef)));
  privateData_->fbrem->push_back(electronRef->fbrem());
  privateData_->nbrems->push_back(electronRef->numberOfBrems());
  privateData_->ambiguousGsfTracksSize->push_back(electronRef->ambiguousGsfTracksSize());
  privateData_->hOverE->push_back(electronRef->hadronicOverEm());
  privateData_->eSuperClusterOverP->push_back(electronRef->eSuperClusterOverP());
  privateData_->eSeedOverPout->push_back(electronRef->eSeedClusterOverPout());
  privateData_->eEleClusterOverPout->push_back(electronRef->eEleClusterOverPout());
  privateData_->deltaEtaAtVtx->push_back(electronRef->deltaEtaSuperClusterTrackAtVtx());
  privateData_->deltaPhiAtVtx->push_back(electronRef->deltaPhiSuperClusterTrackAtVtx());
  privateData_->deltaEtaAtCalo->push_back(electronRef->deltaEtaSeedClusterTrackAtCalo());
  privateData_->deltaPhiAtCalo->push_back(electronRef->deltaPhiSeedClusterTrackAtCalo());
  privateData_->deltaEtaEleClusterTrackAtCalo->push_back(electronRef->deltaEtaEleClusterTrackAtCalo());
  privateData_->deltaPhiEleClusterTrackAtCalo->push_back(electronRef->deltaPhiEleClusterTrackAtCalo());

  // results of standard electron ID sequences
  const eleIdMap & eleIdLikelihoodVal = *( (*eleIdResults_)[0] );

  privateData_->eleLik->push_back( eleIdLikelihoodVal[electronRef] );  
  privateData_->pflowMVA->push_back( electronRef->mva() );

  EcalClusterLazyTools lazyTools( iEvent, iSetup, EcalBarrelRecHits_, EcalEndcapRecHits_ );
  const TransientTrackBuilder thebuilder = *(trackBuilder_.product());

  const reco::GsfElectron *ele = &(*electronRef);
  double mvaidnontrig = myMVANonTrig->mvaValue(*ele,
                                               primaryVertex->front(),
                                               thebuilder,
                                               lazyTools,
                                               false);
  privateData_->mvaidnontrig->push_back(mvaidnontrig);

  double mvaidtrig = myMVATrig->mvaValue(*ele,
                                         primaryVertex->front(),
                                         thebuilder,
                                         lazyTools,
                                         false);
  privateData_->mvaidtrig->push_back(mvaidtrig);

  double mvaidisotrig = myMVATrigIdIsoCombined->IDIsoCombinedMvaValue(*ele,
                                                                      primaryVertex->front(),
                                                                      thebuilder,
                                                                      lazyTools,
                                                                      *pfcandidates_,
                                                                      *rho_,
                                                                      ElectronEffectiveArea::kEleEAData2012,
                                                                      false);
  privateData_->mvaidisotrig->push_back(mvaidisotrig);

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
  
  scBasedIsolation.setExtRadius(0.4);
  scBasedIsolation.excludeHalo(false);
  float scBasedEcalSum04 = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
  privateData_->scBasedEcalSum04->push_back(scBasedEcalSum04);

  float hcalDepth1TowerSumEt03 = hadDepth1Isolation03_->getTowerEtSum(&(*electronRef));
  float hcalDepth2TowerSumEt03 = hadDepth2Isolation03_->getTowerEtSum(&(*electronRef));
  float hcalDepth1TowerSumEt04 = hadDepth1Isolation04_->getTowerEtSum(&(*electronRef));
  float hcalDepth2TowerSumEt04 = hadDepth2Isolation04_->getTowerEtSum(&(*electronRef));

  privateData_->dr03HcalTowerSumEtFullCone->push_back(hcalDepth1TowerSumEt03+hcalDepth2TowerSumEt03);
  privateData_->dr04HcalTowerSumEtFullCone->push_back(hcalDepth1TowerSumEt04+hcalDepth2TowerSumEt04);

  const isoFromPFCandsMap & electronsPfCombinedIsoVal = *( (*eIsoFromPFCandsValueMap_)[0] );
  const isoFromPFCandsMap & electronsPfCandChHad01IsoVal = *( (*eIsoFromPFCandsValueMap_)[1] );
  const isoFromPFCandsMap & electronsPfCandNHad01IsoVal = *( (*eIsoFromPFCandsValueMap_)[2] );
  const isoFromPFCandsMap & electronsPfCandPhoton01IsoVal = *( (*eIsoFromPFCandsValueMap_)[3] );
  const isoFromPFCandsMap & electronsPfCandChHad02IsoVal = *( (*eIsoFromPFCandsValueMap_)[4] );
  const isoFromPFCandsMap & electronsPfCandNHad02IsoVal = *( (*eIsoFromPFCandsValueMap_)[5] );
  const isoFromPFCandsMap & electronsPfCandPhoton02IsoVal = *( (*eIsoFromPFCandsValueMap_)[6] );
  const isoFromPFCandsMap & electronsPfCandChHad03IsoVal = *( (*eIsoFromPFCandsValueMap_)[7] );
  const isoFromPFCandsMap & electronsPfCandNHad03IsoVal = *( (*eIsoFromPFCandsValueMap_)[8] );
  const isoFromPFCandsMap & electronsPfCandPhoton03IsoVal = *( (*eIsoFromPFCandsValueMap_)[9] );
  const isoFromPFCandsMap & electronsPfCandChHad04IsoVal = *( (*eIsoFromPFCandsValueMap_)[10] );
  const isoFromPFCandsMap & electronsPfCandNHad04IsoVal = *( (*eIsoFromPFCandsValueMap_)[11] );
  const isoFromPFCandsMap & electronsPfCandPhoton04IsoVal = *( (*eIsoFromPFCandsValueMap_)[12] );
  const isoFromPFCandsMap & electronsPfCandChHad05IsoVal = *( (*eIsoFromPFCandsValueMap_)[13] );
  const isoFromPFCandsMap & electronsPfCandNHad05IsoVal = *( (*eIsoFromPFCandsValueMap_)[14] );
  const isoFromPFCandsMap & electronsPfCandPhoton05IsoVal = *( (*eIsoFromPFCandsValueMap_)[15] );
  const isoFromPFCandsMap & electronsPfCandChHad06IsoVal = *( (*eIsoFromPFCandsValueMap_)[16] );
  const isoFromPFCandsMap & electronsPfCandNHad06IsoVal = *( (*eIsoFromPFCandsValueMap_)[17] );
  const isoFromPFCandsMap & electronsPfCandPhoton06IsoVal = *( (*eIsoFromPFCandsValueMap_)[18] );
  const isoFromPFCandsMap & electronsPfCandChHad07IsoVal = *( (*eIsoFromPFCandsValueMap_)[19] );
  const isoFromPFCandsMap & electronsPfCandNHad07IsoVal = *( (*eIsoFromPFCandsValueMap_)[20] );
  const isoFromPFCandsMap & electronsPfCandPhoton07IsoVal = *( (*eIsoFromPFCandsValueMap_)[21] );
  const isoFromPFCandsMap & electronsPfCandChHad04DirIsoVal = *( (*eIsoFromPFCandsValueMap_)[22] );
  const isoFromPFCandsMap & electronsPfCandNHad04DirIsoVal = *( (*eIsoFromPFCandsValueMap_)[23] );
  const isoFromPFCandsMap & electronsPfCandPhoton04DirIsoVal = *( (*eIsoFromPFCandsValueMap_)[24] );

  privateData_->pfCombinedIso->push_back( electronsPfCombinedIsoVal[electronRef] );
  privateData_->pfCandChargedIso01->push_back( electronsPfCandChHad01IsoVal[electronRef] );
  privateData_->pfCandNeutralIso01->push_back( electronsPfCandNHad01IsoVal[electronRef] );
  privateData_->pfCandPhotonIso01->push_back( electronsPfCandPhoton01IsoVal[electronRef] );
  privateData_->pfCandChargedIso02->push_back( electronsPfCandChHad02IsoVal[electronRef] );
  privateData_->pfCandNeutralIso02->push_back( electronsPfCandNHad02IsoVal[electronRef] );
  privateData_->pfCandPhotonIso02->push_back( electronsPfCandPhoton02IsoVal[electronRef] );
  privateData_->pfCandChargedIso03->push_back( electronsPfCandChHad03IsoVal[electronRef] );
  privateData_->pfCandNeutralIso03->push_back( electronsPfCandNHad03IsoVal[electronRef] );
  privateData_->pfCandPhotonIso03->push_back( electronsPfCandPhoton03IsoVal[electronRef] );
  privateData_->pfCandChargedIso04->push_back( electronsPfCandChHad04IsoVal[electronRef] );
  privateData_->pfCandNeutralIso04->push_back( electronsPfCandNHad04IsoVal[electronRef] );
  privateData_->pfCandPhotonIso04->push_back( electronsPfCandPhoton04IsoVal[electronRef] );
  privateData_->pfCandChargedIso05->push_back( electronsPfCandChHad05IsoVal[electronRef] );
  privateData_->pfCandNeutralIso05->push_back( electronsPfCandNHad05IsoVal[electronRef] );
  privateData_->pfCandPhotonIso05->push_back( electronsPfCandPhoton05IsoVal[electronRef] );
  privateData_->pfCandChargedIso06->push_back( electronsPfCandChHad06IsoVal[electronRef] );
  privateData_->pfCandNeutralIso06->push_back( electronsPfCandNHad06IsoVal[electronRef] );
  privateData_->pfCandPhotonIso06->push_back( electronsPfCandPhoton06IsoVal[electronRef] );
  privateData_->pfCandChargedIso07->push_back( electronsPfCandChHad07IsoVal[electronRef] );
  privateData_->pfCandNeutralIso07->push_back( electronsPfCandNHad07IsoVal[electronRef] );
  privateData_->pfCandPhotonIso07->push_back( electronsPfCandPhoton07IsoVal[electronRef] );
  privateData_->pfCandChargedDirIso04->push_back( electronsPfCandChHad04DirIsoVal[electronRef] );
  privateData_->pfCandNeutralDirIso04->push_back( electronsPfCandNHad04DirIsoVal[electronRef] );
  privateData_->pfCandPhotonDirIso04->push_back( electronsPfCandPhoton04DirIsoVal[electronRef] );

}




void CmsEleIDTreeFiller::treeEleInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"classification"+colSuffix).c_str(), *privateData_->classification, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"standardClassification"+colSuffix).c_str(), *privateData_->standardClassification, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"fbrem"+colSuffix).c_str(), *privateData_->fbrem, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nbrems"+colSuffix).c_str(), *privateData_->nbrems, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"ambiguousGsfTracksSize"+colSuffix).c_str(), *privateData_->ambiguousGsfTracksSize, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hOverE"+colSuffix).c_str(), *privateData_->hOverE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eSuperClusterOverP"+colSuffix).c_str(), *privateData_->eSuperClusterOverP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eSeedOverPout"+colSuffix).c_str(), *privateData_->eSeedOverPout, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eEleClusterOverPout"+colSuffix).c_str(), *privateData_->eEleClusterOverPout, nCandString.c_str(), 0, "Reco");  
  cmstree->column((colPrefix+"deltaEtaAtVtx"+colSuffix).c_str(), *privateData_->deltaEtaAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiAtVtx"+colSuffix).c_str(), *privateData_->deltaPhiAtVtx, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaEtaAtCalo"+colSuffix).c_str(), *privateData_->deltaEtaAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiAtCalo"+colSuffix).c_str(), *privateData_->deltaPhiAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaEtaEleClusterTrackAtCalo"+colSuffix).c_str(), *privateData_->deltaEtaEleClusterTrackAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiEleClusterTrackAtCalo"+colSuffix).c_str(), *privateData_->deltaPhiEleClusterTrackAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03TkSumPt"+colSuffix).c_str(), *privateData_->dr03TkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03EcalRecHitSumEt"+colSuffix).c_str(), *privateData_->dr03EcalRecHitSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03HcalTowerSumEt"+colSuffix).c_str(), *privateData_->dr03HcalTowerSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04TkSumPt"+colSuffix).c_str(), *privateData_->dr04TkSumPt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04EcalRecHitSumEt"+colSuffix).c_str(), *privateData_->dr04EcalRecHitSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04HcalTowerSumEt"+colSuffix).c_str(), *privateData_->dr04HcalTowerSumEt, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"scBasedEcalSum03"+colSuffix).c_str(), *privateData_->scBasedEcalSum03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"scBasedEcalSum04"+colSuffix).c_str(), *privateData_->scBasedEcalSum04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr03HcalTowerSumEtFullCone"+colSuffix).c_str(), *privateData_->dr03HcalTowerSumEtFullCone, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"dr04HcalTowerSumEtFullCone"+colSuffix).c_str(), *privateData_->dr04HcalTowerSumEtFullCone, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"eleIdLikelihood"+colSuffix).c_str(), *privateData_->eleLik, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pflowMVA"+colSuffix).c_str(), *privateData_->pflowMVA, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mvaidnontrig"+colSuffix).c_str(), *privateData_->mvaidnontrig, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mvaidtrig"+colSuffix).c_str(), *privateData_->mvaidtrig, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"mvaidisotrig"+colSuffix).c_str(), *privateData_->mvaidisotrig, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCombinedIso"+colSuffix).c_str(),  *privateData_->pfCombinedIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso01"+colSuffix).c_str(),  *privateData_->pfCandChargedIso01, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso01"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso01, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso01"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso01, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso02"+colSuffix).c_str(),  *privateData_->pfCandChargedIso02, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso02"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso02, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso02"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso02, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso03"+colSuffix).c_str(),  *privateData_->pfCandChargedIso03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso03"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso03"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso04"+colSuffix).c_str(),  *privateData_->pfCandChargedIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso04"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso04"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso05"+colSuffix).c_str(),  *privateData_->pfCandChargedIso05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso05"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso05"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso06"+colSuffix).c_str(),  *privateData_->pfCandChargedIso06, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso06"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso06, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso06"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso06, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso07"+colSuffix).c_str(),  *privateData_->pfCandChargedIso07, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso07"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso07, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso07"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso07, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedDirIso04"+colSuffix).c_str(),  *privateData_->pfCandChargedDirIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralDirIso04"+colSuffix).c_str(),  *privateData_->pfCandNeutralDirIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonDirIso04"+colSuffix).c_str(),  *privateData_->pfCandPhotonDirIso04, nCandString.c_str(), 0, "Reco");
  
}




void CmsEleIDTreeFillerData::initialise() {
  initialiseCandidate();
  classification           = new vector<int>;
  standardClassification   = new vector<int>;
  fbrem                    = new vector<float>;
  nbrems                   = new vector<int>;
  ambiguousGsfTracksSize   = new vector<int>;
  hOverE                   = new vector<float>;
  eSuperClusterOverP       = new vector<float>;
  eSeedOverPout            = new vector<float>;
  eEleClusterOverPout      = new vector<float>;
  deltaEtaAtVtx            = new vector<float>;
  deltaEtaAtCalo           = new vector<float>;
  deltaPhiAtVtx            = new vector<float>;
  deltaPhiAtCalo           = new vector<float>;
  deltaEtaEleClusterTrackAtCalo = new vector<float>;
  deltaPhiEleClusterTrackAtCalo = new vector<float>;
  dr03TkSumPt              = new vector<float>;
  dr03EcalRecHitSumEt      = new vector<float>;
  dr03HcalTowerSumEt       = new vector<float>;
  dr04TkSumPt              = new vector<float>;
  dr04EcalRecHitSumEt      = new vector<float>;
  dr04HcalTowerSumEt       = new vector<float>;
  scBasedEcalSum03         = new vector<float>;
  scBasedEcalSum04         = new vector<float>;
  dr03HcalTowerSumEtFullCone = new vector<float>;
  dr04HcalTowerSumEtFullCone = new vector<float>;
  eleLik                   = new vector<float>;
  pflowMVA                 = new vector<float>;
  mvaidnontrig             = new vector<float>;
  mvaidtrig                = new vector<float>;
  mvaidisotrig             = new vector<float>;
  pfCombinedIso            = new vector<float>;
  pfCandChargedIso01       = new vector<float>;
  pfCandNeutralIso01       = new vector<float>;
  pfCandPhotonIso01        = new vector<float>;
  pfCandChargedIso02       = new vector<float>;
  pfCandNeutralIso02       = new vector<float>;
  pfCandPhotonIso02        = new vector<float>;
  pfCandChargedIso03       = new vector<float>;
  pfCandNeutralIso03       = new vector<float>;
  pfCandPhotonIso03        = new vector<float>;
  pfCandChargedIso04       = new vector<float>;
  pfCandNeutralIso04       = new vector<float>;
  pfCandPhotonIso04        = new vector<float>;
  pfCandChargedIso05       = new vector<float>;
  pfCandNeutralIso05       = new vector<float>;
  pfCandPhotonIso05        = new vector<float>;
  pfCandChargedIso06       = new vector<float>;
  pfCandNeutralIso06       = new vector<float>;
  pfCandPhotonIso06        = new vector<float>;
  pfCandChargedIso07       = new vector<float>;
  pfCandNeutralIso07       = new vector<float>;
  pfCandPhotonIso07        = new vector<float>;
  pfCandChargedDirIso04    = new vector<float>;
  pfCandNeutralDirIso04    = new vector<float>;
  pfCandPhotonDirIso04     = new vector<float>;
}

void CmsEleIDTreeFillerData::clearTrkVectors() {
  clearTrkVectorsCandidate();
  classification           ->clear();
  standardClassification   ->clear();
  fbrem                    ->clear();
  nbrems                   ->clear();
  ambiguousGsfTracksSize   ->clear();
  hOverE                   ->clear();
  eSuperClusterOverP       ->clear();
  eSeedOverPout            ->clear();
  eEleClusterOverPout      ->clear();
  deltaEtaAtVtx            ->clear();
  deltaEtaAtCalo           ->clear();
  deltaPhiAtVtx            ->clear();
  deltaPhiAtCalo           ->clear();
  deltaEtaEleClusterTrackAtCalo->clear();
  deltaPhiEleClusterTrackAtCalo->clear();
  dr03TkSumPt              ->clear();
  dr03EcalRecHitSumEt      ->clear();
  dr03HcalTowerSumEt       ->clear();
  dr04TkSumPt              ->clear();
  dr04EcalRecHitSumEt      ->clear();
  dr04HcalTowerSumEt       ->clear();
  scBasedEcalSum03         ->clear();
  scBasedEcalSum04         ->clear();
  dr03HcalTowerSumEtFullCone -> clear();
  dr04HcalTowerSumEtFullCone -> clear();
  eleLik                   ->clear();
  pflowMVA                 ->clear();
  mvaidnontrig -> clear();
  mvaidtrig    -> clear();
  mvaidisotrig -> clear();
  pfCombinedIso    ->clear();
  pfCandChargedIso01 ->clear();
  pfCandNeutralIso01 ->clear();
  pfCandPhotonIso01  ->clear();
  pfCandChargedIso02 ->clear();
  pfCandNeutralIso02 ->clear();
  pfCandPhotonIso02  ->clear();
  pfCandChargedIso03 ->clear();
  pfCandNeutralIso03 ->clear();
  pfCandPhotonIso03  ->clear();
  pfCandChargedIso04 ->clear();
  pfCandNeutralIso04 ->clear();
  pfCandPhotonIso04  ->clear();
  pfCandChargedIso05 ->clear();
  pfCandNeutralIso05 ->clear();
  pfCandPhotonIso05  ->clear();
  pfCandChargedIso06 ->clear();
  pfCandNeutralIso06 ->clear();
  pfCandPhotonIso06  ->clear();
  pfCandChargedIso07 ->clear();
  pfCandNeutralIso07 ->clear();
  pfCandPhotonIso07  ->clear();
  pfCandChargedDirIso04 ->clear();
  pfCandNeutralDirIso04 ->clear();
  pfCandPhotonDirIso04  ->clear();
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
