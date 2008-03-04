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

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "CLHEP/HepMC/GenEvent.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;

struct CmsElectronFillerData : public CmsCandidateFillerData {

  // All the vectors that will store the stuff
  // going into the tree.

  vector<float> *pxAtInner, *pyAtInner, *pzAtInner, *xAtInner, *yAtInner, *zAtInner;
  vector<float> *pxAtOuter, *pyAtOuter, *pzAtOuter, *xAtOuter, *yAtOuter, *zAtOuter;
  vector<float> *eleTrackNormalizedChi2;
  vector<float> *eleTrackDxy, *eleTrackD0, *eleTrackDsz, *eleTrackDz;
  vector<float> *eleTrackDxyError, *eleTrackD0Error, *eleTrackDszError, *eleTrackDzError;
  vector<float> *eleTrackValidHits, *eleTrackLostHits;
  vector<float> *eleTrackVx, *eleTrackVy, *eleTrackVz; 

  vector<float> *ecal, *eraw, *caloEta, *caloPhi;
  vector<int> *nClu, *nCry;

  vector<float> *e2x2, *e3x3, *e5x5;
  vector<float> *eMax, *e2nd;
  vector<float> *s1s9, *s9s25;
  vector<float> *covEtaEta, *covEtaPhi, *covPhiPhi;
  //  vector<float> *lat, *phiLat, *etaLat, *a20, *a42;
  
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

  //! set the cluster shape association map for ECAL barrel
  void setEcalBarrelClusterShapes( edm::InputTag EcalBarrelClusterShapes ) { EcalBarrelClusterShapes_ = EcalBarrelClusterShapes; }
  //! set the cluster shape association map for ECAL endcap
  void setEcalEndcapClusterShapes( edm::InputTag EcalEndcapClusterShapes ) { EcalEndcapClusterShapes_ = EcalEndcapClusterShapes; }
  //! set the electron ID association map
  void setElectronIdProducer( edm::InputTag electronIDAssocProducer ) { electronIDAssocProducer_ = electronIDAssocProducer; }
  //! set the tracker isolation producer
  void setTkIsolationProducer( edm::InputTag tkIsolationProducer ) { tkIsolationProducer_ = tkIsolationProducer; }
  //! set the HCAL isolation producer with calo towers
  void setTowerIsolationProducer( edm::InputTag towerIsolationProducer ) { towerIsolationProducer_ = towerIsolationProducer; }

 private:
  
  void writeTrkInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&, GsfTrackRef trkRef);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeEcalInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&, SuperClusterRef sclusRef,
		     const reco::BasicClusterShapeAssociationCollection& barrelClShpMap, const reco::BasicClusterShapeAssociationCollection& endcapClShpMap);
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

  edm::InputTag EcalBarrelClusterShapes_;
  edm::InputTag EcalEndcapClusterShapes_;

  edm::InputTag electronIDAssocProducer_;
  edm::InputTag tkIsolationProducer_;
  edm::InputTag towerIsolationProducer_;

  CmsTree *cmstree;


};

#endif // CmsElectronFiller_h
