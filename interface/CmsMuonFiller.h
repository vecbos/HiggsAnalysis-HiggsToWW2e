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

#ifndef CmsMuonFiller_h
#define CmsMuonFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/Common/interface/ValueMap.h"

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
#include <TTree.h>

struct CmsMuonFillerData : public CmsCandidateFillerData {

  vector<int> *trackIndex, *standAloneTrackIndex, *combinedTrackIndex;

  vector<float> *sumPt03, *emEt03, *hadEt03, *hoEt03, *nTrk03, *nJets03;
  vector<float> *sumPt05, *emEt05, *hadEt05, *hoEt05, *nTrk05, *nJets05;
  vector<int> *muonId, *type, *numberOfMatches;
  vector<float> *pfChargedIso, *pfNeutralIso, *pfPhotonIso, *pfGenericChargedIso, *pfGenericNeutralIso, *pfGenericPhotonIso,
    *pfGenericNoOverChargedIso, *pfGenericNoOverNeutralIso, *pfGenericNoOverPhotonIso, *pfCombinedIso;
  vector<float> *kink;

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

  //! set the general tracks (reduced collection)
  void setGeneralTracks( edm::InputTag generalTracks) { generalTracks_ = generalTracks; }
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
  
  void writeTrkInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&, const reco::Muon *muon);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeMuonInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&, const reco::Muon *muon, const reco::MuonRef muonRef);
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

  edm::InputTag generalTracks_;

  CmsTree *cmstree;

  // Geometry
  EcalDetIdAssociator ecalDetIdAssociator_;
  HcalDetIdAssociator hcalDetIdAssociator_;
  HODetIdAssociator   hoDetIdAssociator_;
  CaloDetIdAssociator caloDetIdAssociator_;
  MuonDetIdAssociator muonDetIdAssociator_;


  CachedTrajectory cachedTrajectory_;
  edm::Handle< reco::TrackRefVector > h_tracks;
  edm::ESHandle<MagneticField> bField;
  edm::ESHandle<CaloGeometry> theCaloGeometry_;
  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry_;

  typedef edm::ValueMap<double> isoFromDepositsMap;
  typedef std::vector< edm::Handle<isoFromDepositsMap> > isoContainer;
  isoContainer *eIsoFromDepsValueMap_;

  typedef edm::ValueMap<float> isoFromPFCandsMap;
  typedef std::vector< edm::Handle<isoFromPFCandsMap> > isoContainer2;
  isoContainer2 *eIsoFromPFCandsValueMap_;

};

#endif // CmsMuonFiller_h
