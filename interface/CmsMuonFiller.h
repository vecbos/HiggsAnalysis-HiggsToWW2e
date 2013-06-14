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
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"

#include <TTree.h>

struct CmsMuonFillerData : public CmsCandidateFillerData {

  vector<int> *trackIndex, *standAloneTrackIndex, *combinedTrackIndex;

  vector<float> *sumPt03, *emEt03, *hadEt03, *hoEt03, *nTrk03, *nJets03;
  vector<float> *sumPt05, *emEt05, *hadEt05, *hoEt05, *nTrk05, *nJets05;
  vector<int> *muonId, *type, *numberOfMatches;
  vector<bool> *pfmuonId;
  vector<float> *pfCombinedIso;
  vector<float> *pfCandChargedIso01, *pfCandNeutralIso01, *pfCandPhotonIso01,
    *pfCandChargedIso02, *pfCandNeutralIso02, *pfCandPhotonIso02,
    *pfCandChargedIso03, *pfCandNeutralIso03, *pfCandPhotonIso03,
    *pfCandChargedIso04, *pfCandNeutralIso04, *pfCandPhotonIso04,
    *pfCandChargedIso05, *pfCandNeutralIso05, *pfCandPhotonIso05,
    *pfCandChargedIso06, *pfCandNeutralIso06, *pfCandPhotonIso06,
    *pfCandChargedIso07, *pfCandNeutralIso07, *pfCandPhotonIso07,
    *pfCandChargedDirIso04, *pfCandNeutralDirIso04, *pfCandPhotonDirIso04,
    *pfIsolationSumPUPtR03, *pfIsolationSumPUPtR04,
    *mvaiso,
    *pfCandChargedPUIso03, *pfCandChargedPUIso04;
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
  //! set the vertices collection (needed for MVA isolation)
  void setVertexCollection(edm::InputTag collectionTag);
  //! set the muon MVA algo
  void setMuonIsoMVA(MuonMVAEstimator* algo) { fMuonIsoMVA = algo; }

  // Operators

  //! write the muon related informations for the given collection
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  
 private:
  
  void writeTrkInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&, const reco::Muon *muon);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeMuonInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&, const reco::Muon *muon, const reco::MuonRef muonRef, double rho);
  void treeMuonInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool saveMuonExtras_;
  bool saveTrk_;
  bool saveFatTrk_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsMuonFillerData *privateData_;

  edm::InputTag generalTracks_;
  edm::InputTag m_vxtCollectionTag;

  CmsTree *cmstree;

  edm::Handle< reco::TrackRefVector > h_tracks;

  typedef edm::ValueMap<float> isoFromPFCandsMap;
  typedef std::vector< edm::Handle<isoFromPFCandsMap> > isoContainer;
  isoContainer *eIsoFromPFCandsValueMap_;

  MuonMVAEstimator *fMuonIsoMVA;
  reco::Vertex *firstVtx;
  edm::Handle<reco::PFCandidateCollection> pfCands;
  edm::Handle<reco::VertexCollection> primaryVertex;

};

#endif // CmsMuonFiller_h
