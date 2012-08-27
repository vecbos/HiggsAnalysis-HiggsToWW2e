// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsEleIDTreeFiller
//      Simple class for dumping electron ID contents to a ROOT tree
//      
// Original Authors:  Chiara Ilaria Rovelli, Emanuele Di Marco
//          Created:  Fri May  18 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsEleIDTreeFiller_h
#define CmsEleIDTreeFiller_h

#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/Common/interface/ValueMap.h"

struct CmsEleIDTreeFillerData : public CmsCandidateFillerData {

  vector<int>   *classification, *standardClassification;
  vector<float> *fbrem;
  vector<int>   *nbrems, *ambiguousGsfTracksSize;
  vector<float> *hOverE, *eSuperClusterOverP, *eSeedOverPout, *eEleClusterOverPout;
  vector<float> *deltaEtaAtVtx, *deltaEtaAtCalo, *deltaPhiAtVtx, *deltaPhiAtCalo, *deltaEtaEleClusterTrackAtCalo, *deltaPhiEleClusterTrackAtCalo;
  vector<float> *dr03TkSumPt, *dr03EcalRecHitSumEt, *dr03HcalTowerSumEt;
  vector<float> *dr04TkSumPt, *dr04EcalRecHitSumEt, *dr04HcalTowerSumEt;
  vector<float> *scBasedEcalSum03, *scBasedEcalSum04;
  vector<float> *dr03HcalTowerSumEtFullCone, *dr04HcalTowerSumEtFullCone;
  vector<float> *pfCombinedIso, 
    *pfCandChargedIso01, *pfCandNeutralIso01, *pfCandPhotonIso01,
    *pfCandChargedIso02, *pfCandNeutralIso02, *pfCandPhotonIso02,
    *pfCandChargedIso03, *pfCandNeutralIso03, *pfCandPhotonIso03,
    *pfCandChargedIso04, *pfCandNeutralIso04, *pfCandPhotonIso04,
    *pfCandChargedIso05, *pfCandNeutralIso05, *pfCandPhotonIso05,
    *pfCandChargedIso06, *pfCandNeutralIso06, *pfCandPhotonIso06,
    *pfCandChargedIso07, *pfCandNeutralIso07, *pfCandPhotonIso07,
    *pfCandChargedDirIso04, *pfCandNeutralDirIso04, *pfCandPhotonDirIso04;

  vector<float> *eleLik, *pflowMVA, *mvaidtrig, *mvaidisotrig, *mvaidnontrig;

public:
  void initialise();
  void clearTrkVectors();

};


class CmsEleIDTreeFiller : public CmsCandidateFiller {

 public:

  //! Constructor
  CmsEleIDTreeFiller(CmsTree *, int maxTracks=500, 
		     bool noOutputIfLimitsReached=false );

  //! Destructor
  virtual ~CmsEleIDTreeFiller();

  //! write the electron ID starting from a candidate collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);
  
  //! set the rechits for ECAL barrel (needed for cluster shapes)
  void setEcalBarrelRecHits( edm::InputTag EcalBarrelRecHits ) { EcalBarrelRecHits_ = EcalBarrelRecHits; }
  //! set the rechits for ECAL endcap (needed for cluster shapes)
  void setEcalEndcapRecHits( edm::InputTag EcalEndcapRecHits ) { EcalEndcapRecHits_ = EcalEndcapRecHits; }
  //! set the tracker isolation producer (egamma)
  void setTkIsolationProducer( edm::InputTag tkIsolationProducer ) { tkIsolationProducer_ = tkIsolationProducer; }
  //! set the HCAL isolation producer with calo towers (egamma)
  void setTowerIsolationProducer( edm::InputTag towerIsolationProducer ) { towerIsolationProducer_ = towerIsolationProducer; }
  //! set the track producer for tracker isolation
  void setTracksProducer( edm::InputTag tracksProducer ) { tracksProducer_ = tracksProducer; }
  //! set the calotower producer for calorimetric isolation
  void setCalotowersProducer( edm::InputTag calotowersProducer ) { calotowersProducer_ = calotowersProducer; }
  //! set the vertex collection
  void setVertexCollection(edm::InputTag collectionTag) { m_vxtCollectionTag = collectionTag; }
  //! set the PF candidates collection
  void setPFCandidateCollection(edm::InputTag collectionTag) { m_pfcandCollectionTag = collectionTag; }
  //! set the eleID MVA algos
  void setEleIdMVAs(EGammaMvaEleEstimator* algotrig, EGammaMvaEleEstimator* algotrigidiso, EGammaMvaEleEstimator* algonontrig) { myMVATrig = algotrig; myMVATrigIdIsoCombined = algotrigidiso; myMVANonTrig = algonontrig; }
  //! dump pflow isolation related variables
  void savePFlowIsolations(bool what) { savePFlowIsolation_ = what; }

  //! set to false if the column with the block size is set by another object
  void setStandalone(bool );

 private:

  void writeEleInfo(const reco::GsfElectronRef, const edm::Event&, const edm::EventSetup&,
		    const EcalRecHitCollection *EBRecHits,
		    const EcalRecHitCollection *EERecHits);
  void treeEleInfo(const std::string &colPrefix, const std::string &colSuffix);

  int stdEleIdClassify(const reco::GsfElectron* electron);

  int maxTracks_;
  std::string *trkIndexName_;
  bool standalone_;

  edm::InputTag EcalBarrelRecHits_;
  edm::InputTag EcalEndcapRecHits_;
  
  edm::InputTag electronIdCutsLabel_;
  edm::InputTag electronIdLikelihoodLabel_;
  edm::InputTag tkIsolationProducer_;
  edm::InputTag towerIsolationProducer_;

  edm::InputTag tracksProducer_;
  edm::InputTag calotowersProducer_;
  edm::InputTag m_vxtCollectionTag;
  edm::InputTag m_pfcandCollectionTag;

  bool savePFlowIsolation_;

  CmsTree *cmstree;
  CmsEleIDTreeFillerData *privateData_;

  typedef edm::ValueMap<float> eleIdMap;
  typedef std::vector< edm::Handle<eleIdMap> > eleIdContainer;
  eleIdContainer *eleIdResults_;

  typedef edm::ValueMap<float> isoFromPFCandsMap;
  typedef std::vector< edm::Handle<isoFromPFCandsMap> > isoContainer;
  isoContainer *eIsoFromPFCandsValueMap_;

  edm::Handle<CaloTowerCollection> m_calotowers;
  edm::Handle<reco::TrackCollection> m_tracks;

  EgammaTowerIsolation *hadDepth1Isolation03_, *hadDepth2Isolation03_, 
    *hadDepth1Isolation04_, *hadDepth2Isolation04_;

  EGammaMvaEleEstimator* myMVANonTrig, *myMVATrig, *myMVATrigIdIsoCombined;
  edm::Handle<reco::VertexCollection> primaryVertex;
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;
  edm::Handle<reco::PFCandidateCollection> pfcandidates_;

  edm::Handle<double> rho_;

  const CaloGeometry* caloGeo;
};

#endif // CmsEleIDTreeFiller_h
