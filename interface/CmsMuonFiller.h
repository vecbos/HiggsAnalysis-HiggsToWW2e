// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      ZjetsAnalysis/ZllProducer
// Description:
//      Class CmsMuonFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
//-----------------------------------------------------------------------

#ifndef CmsMuonFiller_h
#define CmsMuonFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/TrackAssociator/interface/CachedTrajectory.h"
#include "TrackingTools/TrackAssociator/interface/CaloDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/EcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/MuonDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/HcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/HODetIdAssociator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "CLHEP/HepMC/GenEvent.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;

struct CmsMuonFillerData : public CmsCandidateFillerData {

  vector<float> *pxAtInner, *pyAtInner, *pzAtInner, *xAtInner, *yAtInner, *zAtInner;
  vector<float> *pxAtOuter, *pyAtOuter, *pzAtOuter, *xAtOuter, *yAtOuter, *zAtOuter;
  vector<float> *muTrackNormalizedChi2;
  vector<float> *muTrackDxy, *muTrackD0, *muTrackDsz, *muTrackDz;
  vector<float> *muTrackDxyError, *muTrackD0Error, *muTrackDszError, *muTrackDzError;
  vector<float> *muTrackValidHits, *muTrackLostHits;
  vector<float> *muTrackVx, *muTrackVy, *muTrackVz; 

  vector<float> *pxECAL, *pyECAL, *pzECAL, *xECAL, *yECAL, *zECAL, *EexpECAL; 
  vector<float> *pxHCAL, *pyHCAL, *pzHCAL, *xHCAL, *yHCAL, *zHCAL, *EexpHCAL; 
  vector<float> *pxHO, *pyHO, *pzHO, *xHO, *yHO, *zHO, *EexpHO; 

  vector<float> *sumPt03, *emEt03, *hadEt03, *hoEt03, *nTrk03, *nJets03;
  vector<float> *sumPt05, *emEt05, *hadEt05, *hoEt05, *nTrk05, *nJets05;

  vector<float> *EcalExpDepo, *HcalExpDepo, *HoExpDepo, *emS9, *hadS9, *hoS9, *CaloComp;

public:
  void initialise();
  void clearTrkVectors();

};

class CmsMuonFiller : public CmsCandidateFiller {

 public:

  // Constructors

  // Dump everything
  CmsMuonFiller(CmsTree *, 
		int maxTracks=500, int maxMCTracks=2000, 
		bool noOutputIfLimitsReached=false );

  // Dump  everything if fatTree is true and less informations otherwise
  CmsMuonFiller(CmsTree *, 
		bool fatTree, 
		int maxTracks=500, int maxMCTracks=2000, 
		bool noOutputIfLimitsReached=false );
  
  // Destructor
  virtual ~CmsMuonFiller();

  //! dump more  variables
  void saveMuonExtras(bool );
  //! dump tracking related variables
  void saveTrk(bool );
  //! dump more ECAL related variables
  void saveFatTrk(bool );

  // Operators

  //! write the muon related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  
  void SetGeometry(const edm::EventSetup& iSetup);

 private:
  
  void writeTrkInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&, TrackRef trkRef);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeMuonInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&, const Muon *muon);
  void treeMuonInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool saveMuonExtras_;
  bool saveTrk_;
  bool saveFatTrk_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsMuonFillerData *privateData_;
  edm::InputTag matchMap_;

  CmsTree *cmstree;

  // Geometry
  EcalDetIdAssociator ecalDetIdAssociator_;
  HcalDetIdAssociator hcalDetIdAssociator_;
  HODetIdAssociator   hoDetIdAssociator_;
  CaloDetIdAssociator caloDetIdAssociator_;
  MuonDetIdAssociator muonDetIdAssociator_;


  CachedTrajectory cachedTrajectory_;
  edm::ESHandle<MagneticField> bField;
  edm::ESHandle<CaloGeometry> theCaloGeometry_;
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry_;

};

#endif // CmsMuonFiller_h
