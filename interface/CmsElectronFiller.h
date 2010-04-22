// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsElectronFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsElectronFiller_h
#define CmsElectronFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include <TTree.h>

struct CmsElectronFillerData : public CmsCandidateFillerData {

  // All the vectors that will store the stuff
  // going into the tree.
  vector<int> *trackIndex, *gsfTrackIndex;

  vector<int> *fiducialFlags, *recoFlags, *scPixCharge;

  vector<int> *superClusterIndex, *PFsuperClusterIndex;
  vector<float> *ecal, *eraw, *esEnergy, *caloEta, *caloPhi;
  vector<int> *energyCorrections;

  vector<float> *convDist, *convDcot, *convRadius, *convX, *convY, *convZ, *convChi2Prob;
  vector<int> *convTrackIndex;

public:
  void initialise();
  void clearTrkVectors();


};

class CmsElectronFiller : public CmsCandidateFiller {

 public:

  //! Dump everything
  CmsElectronFiller(CmsTree *, int maxTracks=500,
		int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  //! Dump  everything if fatTree is true and less informations otherwise
  CmsElectronFiller(CmsTree *, bool fatTree, int maxTracks=500,
		int maxMCTracks=2000, bool noOutputIfLimitsReached=false );


  //! Destructor
  virtual ~CmsElectronFiller();

  //! dump tracking related variables
  void saveTrk(bool );
  //! dump ECAL related variables
  void saveEcal(bool );
  //! dump more tracking related variables
  void saveFatTrk(bool );
  //! dump more ECAL related variables
  void saveFatEcal(bool );
  //! dump electron ID variables
  void saveEleID(bool );

  //! write the electron related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  //! set the general tracks
  void setGeneralTracks( edm::InputTag generalTracks) { generalTracks_ = generalTracks; }
  //! set the superclusters for ECAL barrel
  void setEcalBarrelSuperClusters( edm::InputTag EcalBarrelSuperClusters) { EcalBarrelSuperClusters_ = EcalBarrelSuperClusters; }
  //! set the superclusters for ECAL endcap
  void setEcalEndcapSuperClusters( edm::InputTag EcalEndcapSuperClusters) { EcalEndcapSuperClusters_ = EcalEndcapSuperClusters; }
  //! set the rechits for ECAL barrel (needed for cluster shapes)
  void setEcalBarrelRecHits( edm::InputTag EcalBarrelRecHits ) { EcalBarrelRecHits_ = EcalBarrelRecHits; }
  //! set the rechits for ECAL endcap (needed for cluster shapes)
  void setEcalEndcapRecHits( edm::InputTag EcalEndcapRecHits ) { EcalEndcapRecHits_ = EcalEndcapRecHits; }
  //! set the tracker isolation producer
  void setTkIsolationProducer( edm::InputTag tkIsolationProducer ) { tkIsolationProducer_ = tkIsolationProducer; }
  //! set the HCAL isolation producer with calo towers
  void setTowerIsolationProducer( edm::InputTag towerIsolationProducer ) { towerIsolationProducer_ = towerIsolationProducer; }
  //! set the track producer for tracker isolation
  void setTracksProducer( edm::InputTag tracksProducer ) { tracksProducer_ = tracksProducer; }
  //! set the calotower producer for calorimetric isolation
  void setCalotowersProducer( edm::InputTag calotowersProducer ) { calotowersProducer_ = calotowersProducer; }


 private:
  
  void writeTrkInfo(const reco::GsfElectronRef, const edm::Event&, const edm::EventSetup&, reco::GsfTrackRef trkRef);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeEcalInfo(const reco::GsfElectronRef, const edm::Event&, const edm::EventSetup&, 
                     reco::SuperClusterRef sclusRef, reco::SuperClusterRef pfclusRef,
		     const EcalRecHitCollection *EBRecHits, const EcalRecHitCollection *EERecHits);
  void treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix);


  bool saveTrk_;
  bool saveEcal_;
  bool saveFatTrk_;
  bool saveFatEcal_;
  bool saveEleID_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsElectronFillerData *privateData_;
  edm::InputTag matchMap_;

  edm::InputTag generalTracks_;
  edm::InputTag EcalBarrelSuperClusters_, EcalEndcapSuperClusters_;

  edm::InputTag EcalBarrelRecHits_;
  edm::InputTag EcalEndcapRecHits_;

  edm::InputTag electronIdCutsLabel_;
  edm::InputTag electronIdLikelihoodLabel_;
  edm::InputTag tkIsolationProducer_;
  edm::InputTag towerIsolationProducer_;

  edm::InputTag tracksProducer_;
  edm::InputTag calotowersProducer_;

  edm::Handle< reco::TrackCollection > h_tracks;

  int barrelSuperClustersSize;

  CmsTree *cmstree;


};

#endif // CmsElectronFiller_h
