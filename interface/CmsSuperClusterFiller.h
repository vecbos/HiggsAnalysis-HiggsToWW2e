// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsSuperClusterFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsSuperClusterFiller_h
#define CmsSuperClusterFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

struct CmsSuperClusterFillerData {
  // All the vectors that will store the stuff
  // going into the tree.
  vector<int>  *nBC, *nCrystals, *iAlgo, *recoFlag, *sevClosProbl, *idClosProbl;
  vector<float> *rawEnergy, *energy, *seedEnergy, *eta, *theta, *phi, *time, *chi2Prob, *fracClosProbl;
  vector<float> *e3x3, * e5x5, *eMax, *e2x2, *e2nd, *covIEtaIEta, *covIEtaIPhi, *covIPhiIPhi;
  vector<int> *trackIndex;
  vector<float> *deltaR, *deltaPhi, *deltaEta;
  int *nSC;

public:
  void initialiseCandidate();
  void clear();
};


class CmsSuperClusterFiller {

public:

  // Constructors

  // Dump everything
  CmsSuperClusterFiller(CmsTree *, int maxSC=500);


  // Destructor
  virtual ~CmsSuperClusterFiller();

  // Operators

  // run number and all of that --- to implement

  virtual void writeCollectionToTree(edm::InputTag collectionTag,
				     const edm::Event&, const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);

  //! set the track collection (to match the superclusters as hand-made electrons)
  void setTracks( edm::InputTag Tracks ) { Tracks_ = Tracks; }
  //! set the rechits for ECAL barrel (needed for cluster shapes)
  void setEcalBarrelRecHits( edm::InputTag EcalBarrelRecHits ) { EcalBarrelRecHits_ = EcalBarrelRecHits; }
  //! set the rechits for ECAL endcap (needed for cluster shapes)
  void setEcalEndcapRecHits( edm::InputTag EcalEndcapRecHits ) { EcalEndcapRecHits_ = EcalEndcapRecHits; }


protected:
  
  // fraction of SC energy around closest problematic
  float fractionAroundClosestProblematic( const reco::CaloCluster & , const EcalRecHitCollection &,  const EcalChannelStatus &, const CaloTopology* topology );
  // retrieve closest problematic channel and its severity wrt seed crystal using as distance sqrt(ieta^2+ieta^2+iphi^2+iphi^2). Return a null detId in case not found within a search region of 11 (ieta) x 51 (iphi)  
  std::pair<DetId,int> closestProblematic( const reco::CaloCluster & , const EcalRecHitCollection &,  const EcalChannelStatus &, const CaloTopology* topology );


  //return the distance in eta units between two EBDetId
  static int distanceEta(const EBDetId& a,const EBDetId& b); 
  //return the distance in phi units between two EBDetId
  static int distancePhi(const EBDetId& a,const EBDetId& b);

  virtual void writeSCInfo(const reco::SuperCluster *cand, 
			   const edm::Event&, const edm::EventSetup&,
                           const EcalRecHitCollection *EBRecHits, const EcalRecHitCollection *EERecHits, 
                           const reco::TrackCollection *Tracks);
  virtual void treeSCInfo(const std::string colPrefix, const std::string colSuffix);
  
  // Friends

  CmsSuperClusterFillerData *privateData_;

  CmsTree *cmstree;

  edm::InputTag Tracks_;
  edm::InputTag EcalBarrelRecHits_;
  edm::InputTag EcalEndcapRecHits_;

  DetId closestProb_;
  int severityClosestProb_;

  int maxSC_;
  std::string *trkIndexName_;
};

#endif // CmsSuperClusterFiller_h
