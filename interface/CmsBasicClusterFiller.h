// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsBasicClusterFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsBasicClusterFiller_h
#define CmsBasicClusterFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

struct CmsBasicClusterFillerData {
  // All the vectors that will store the stuff
  // going into the tree.
  vector<int>  *nCrystals, *recoFlag, *sevClosProbl, *idClosProbl, *indexSC;
  vector<float> *energy, *seedEnergy, *eta, *theta, *phi, *time, *chi2, *fracClosProbl;
  vector<float> *e3x3, * e5x5, *eMax, *e2x2, *e2nd, *covIEtaIEta, *covIEtaIPhi, *covIPhiIPhi;
  vector<float> *eTop, *eBottom, *eLeft, *eRight;
  vector<float> *etaCrystal,*phiCrystal,*thetaTilt,*phiTilt;
  vector<int> *iEta,*iPhi;
  int *nBC;

public:
  void initialiseCandidate();
  void clear();
};


class CmsBasicClusterFiller {

public:

  // Constructors

  // Dump everything
  CmsBasicClusterFiller(CmsTree *, int maxBC=500);


  // Destructor
  virtual ~CmsBasicClusterFiller();

  // Operators

  // run number and all of that --- to implement

  virtual void writeCollectionToTree(edm::InputTag collectionTag,
				     const edm::Event&, const edm::EventSetup&,
				     const std::string &columnPrefix, const std::string &columnSuffix,
				     bool dumpData=false);

  //! set the calotowers collection (to calculate H/E)
  void setCalotowers( edm::InputTag Calotowers ) { Calotowers_ = Calotowers; }
  //! set the rechits for ECAL barrel (needed for cluster shapes)
  void setEcalBarrelRecHits( edm::InputTag EcalBarrelRecHits ) { EcalBarrelRecHits_ = EcalBarrelRecHits; }
  //! set the rechits for ECAL endcap (needed for cluster shapes)
  void setEcalEndcapRecHits( edm::InputTag EcalEndcapRecHits ) { EcalEndcapRecHits_ = EcalEndcapRecHits; }
  //! set SC  
  void setEcalSuperClusters( edm::InputTag a_tag ) { EcalSC_ = a_tag; }



protected:
  
  // fraction of SC energy around closest problematic
  float fractionAroundClosestProblematic( const edm::EventSetup & , const reco::CaloCluster & , const EcalRecHitCollection &,  const EcalChannelStatus &, const CaloTopology* topology );
  // retrieve closest problematic channel and its severity wrt seed crystal using as distance sqrt(ieta^2+ieta^2+iphi^2+iphi^2). Return a null detId in case not found within a search region of 11 (ieta) x 51 (iphi)  
  std::pair<DetId,int> closestProblematic( const edm::EventSetup & , const reco::CaloCluster & , const EcalRecHitCollection &,  const EcalChannelStatus &, const CaloTopology* topology );


  //return the distance in eta units between two EBDetId
  static int distanceEta(const EBDetId& a,const EBDetId& b); 
  //return the distance in phi units between two EBDetId
  static int distancePhi(const EBDetId& a,const EBDetId& b);

  virtual void writeBCInfo(const reco::BasicCluster *cand, 
			   const edm::Event&, const edm::EventSetup&,
                           const EcalRecHitCollection *EBRecHits, const EcalRecHitCollection *EERecHits, 
			   const reco::SuperClusterCollection *ESCCollectionEB );
  virtual void treeBCInfo(const std::string colPrefix, const std::string colSuffix);
  

  // Friends
  CmsTree *cmstree;

  CmsBasicClusterFillerData *privateData_;

  // FC added SC link 
  edm::InputTag EcalSC_; // not mandatory 


  edm::InputTag EcalBarrelRecHits_;
  edm::InputTag EcalEndcapRecHits_;
  edm::InputTag Calotowers_;

  edm::Handle<CaloTowerCollection> calotowers_;

  int nBCCandOfSC;

  DetId closestProb_;
  int severityClosestProb_;

  int maxBC_;
  std::string *trkIndexName_;
};

#endif // CmsBasicClusterFiller_h
