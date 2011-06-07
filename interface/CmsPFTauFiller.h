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

#ifndef CmsPFTauFiller_h
#define CmsPFTauFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/TrackAssociator/interface/CachedTrajectory.h"
#include "TrackingTools/TrackAssociator/interface/CaloDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/EcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/MuonDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/HcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/HODetIdAssociator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

// #include "DataFormats/TauReco/interface/BaseTau.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminatorByIsolation.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"
#include "RecoTauTag/TauTagTools/interface/PFTauElementsOperators.h"
#include "RecoTauTag/TauTagTools/interface/TauTagTools.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include <TTree.h>

struct CmsPFTauFillerData : public CmsCandidateFillerData {

  // from the basic
  vector<int> *isNonNull;
  vector<int> *charge;
  vector<float> *energy, *et, *momentum, *pt;
  vector<float> *vertexX, *vertexY, *vertexZ;
  vector<float> *theta, *eta, *phi;
  vector<float> *x, *y, *z;
  vector<float> *mass, *mt;
  vector<int> *pdgId;
  vector<int> *nDau;

  // from the Leading Track: for Electrons
  vector<float> *isolationPFChargedHadrCandsPtSum;
  vector<float> *isolationPFGammaCandsEtSum;
  vector<float> *emFraction;
  vector<float> *hcalTotOverPLead;
  vector<float> *hcal3x3OverPLead;
  vector<float> *ecalStripSumEOverPLead;
  vector<float> *bremsRecoveryEOverPLead;


  // from the Discriminator
  vector<float> *theTauDiscrByLeadTrackFinding;
  vector<float> *theTauDiscrByLeadTrackPtCut;
  // vector<float> *theTauDiscrByNProngs;
  vector<float> *theTauDiscrByTrackIso;
  vector<float> *theTauDiscrByEcalIso;
  vector<float> *theTauDiscrAgainstMuons;
  vector<float> *theTauDiscrAgainstElectrons;
//   vector<float> *theTauDiscrByTaNC;
//   vector<float> *theTauDiscrByTaNCfrHalfPercent;
//   vector<float> *theTauDiscrByTaNCfrOnePercent;
//   vector<float> *theTauDiscrByTaNCfrQuarterPercent;
//   vector<float> *theTauDiscrByTaNCfrTenthPercent;


public:
  void initialise();
  void clear();
};

class CmsPFTauFiller : public CmsCandidateFiller {

public:

  // Constructor
  CmsPFTauFiller(CmsTree *, 
		 int maxTracks=500, int maxMCTracks=2000, 
		 bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsPFTauFiller();

  //! dump more  variables
  void savePFTauBasic(bool );
  void saveLeadPFCand(bool );

  // Operators
  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     edm::InputTag tauDiscrByLeadTrackFindingTag,
			     edm::InputTag tauDiscrByLeadTrackPtCutTag,
			     // edm::InputTag tauDiscrByNProngsTag,
			     edm::InputTag tauDiscrByTrackIsoTag,
			     edm::InputTag tauDiscrByEcalIsoTag, 
			     edm::InputTag tauDiscrAgainstMuonsTag,
			     edm::InputTag tauDiscrAgainstElectronsTag,
// 			     edm::InputTag tauDiscrByTaNCTag,
// 			     edm::InputTag tauDiscrByTaNCfrHalfPercentTag,
// 			     edm::InputTag tauDiscrByTaNCfrOnePercentTag,
// 			     edm::InputTag tauDiscrByTaNCfrQuarterPercentTag,
// 			     edm::InputTag tauDiscrByTaNCfrTenthPercentTag,
			     bool dumpData=false);
  
private:
  
  // void writePFTauBasicInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&, const reco::PFTau *tau);
  void writePFTauBasicInfo(const reco::PFTau *cand, const edm::Event&, const edm::EventSetup&);

  void writeLeadPFCandInfo(const reco::PFTau *cand, const edm::Event&, const edm::EventSetup&);

  void writePFTauDiscInfo(edm::Handle<reco::PFTauCollection> tauCollection, int theTauJetIndex,
			  const reco::PFTauDiscriminator *tauDiscrByLeadTrackFinding,
			  const reco::PFTauDiscriminator *tauDiscrByLeadTrackPtCut,
			  // const reco::PFTauDiscriminator *tauDiscrByNProngs,
			  const reco::PFTauDiscriminator *tauDiscrByTrackIso,
			  const reco::PFTauDiscriminator *tauDiscrByEcalIso,
			  const reco::PFTauDiscriminator *tauDiscrAgainstMuons,
			  const reco::PFTauDiscriminator *tauDiscrAgainstElectrons
// 			  const reco::PFTauDiscriminator *tauDiscrByTaNC,
// 			  const reco::PFTauDiscriminator *tauDiscrByTaNCfrHalfPercent,
// 			  const reco::PFTauDiscriminator *tauDiscrByTaNCfrOnePercent,
// 			  const reco::PFTauDiscriminator *tauDiscrByTaNCfrQuarterPercent,
// 			  const reco::PFTauDiscriminator *tauDiscrByTaNCfrTenthPercent
                          );

  void treePFTauBasicInfo(const std::string &colPrefix, const std::string &colSuffix);
  void treePFTauDiscInfo(const std::string &colPrefix, const std::string &colSuffix);
  void treeLeadPFCandInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool savePFTauBasic_;
  bool saveLeadPFCand_;
  bool savePFTauDiscriminators_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsPFTauFillerData *privateData_;

  CmsTree *cmstree;

};
#endif // CmsPFTauFiller_h
