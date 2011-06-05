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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/Common/interface/ValueMap.h"

struct CmsEleIDTreeFillerData : public CmsCandidateFillerData {

  vector<int>   *classification, *standardClassification;
  vector<float> *fbrem;
  vector<int>   *nbrems;
  vector<float> *hOverE, *eSuperClusterOverP, *eSeedOverPout;
  vector<float> *deltaEtaAtVtx, *deltaEtaAtCalo, *deltaPhiAtVtx, *deltaPhiAtCalo;
  vector<float> *dr03TkSumPt, *dr03EcalRecHitSumEt, *dr03HcalTowerSumEt;
  vector<float> *dr04TkSumPt, *dr04EcalRecHitSumEt, *dr04HcalTowerSumEt;
  vector<float> *scBasedEcalSum03, *scBasedEcalSum04;
  vector<float> *dr03HcalTowerSumEtFullCone, *dr04HcalTowerSumEtFullCone;
  vector<float> *pfChargedIso, *pfNeutralIso, *pfPhotonIso, *pfGenericChargedIso, *pfGenericNeutralIso, *pfGenericPhotonIso,
    *pfGenericNoOverChargedIso, *pfGenericNoOverNeutralIso, *pfGenericNoOverPhotonIso;

  vector<float> *eleLik, *pflowMVA;

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

  CmsTree *cmstree;
  CmsEleIDTreeFillerData *privateData_;

  typedef edm::ValueMap<float> eleIdMap;
  typedef std::vector< edm::Handle<eleIdMap> > eleIdContainer;
  eleIdContainer *eleIdResults_;

  typedef edm::ValueMap<double> isoFromDepositsMap;
  typedef std::vector< edm::Handle<isoFromDepositsMap> > isoContainer;
  isoContainer *eIsoFromDepsValueMap_;

  edm::Handle<CaloTowerCollection> m_calotowers;
  edm::Handle<reco::TrackCollection> m_tracks;

  EgammaTowerIsolation *hadDepth1Isolation03_, *hadDepth2Isolation03_, 
    *hadDepth1Isolation04_, *hadDepth2Isolation04_;

  const CaloGeometry* caloGeo;
};

#endif // CmsEleIDTreeFiller_h
