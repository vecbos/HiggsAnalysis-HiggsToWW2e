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
#include "MagneticField/Engine/interface/MagneticField.h"

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
  vector<int> *pfJetIndex;
  vector<int> *isNonNull;

  // from the Leading Track: for Electrons
  vector<float> *isolationPFChargedHadrCandsPtSum;
  vector<float> *isolationPFGammaCandsEtSum;
  vector<float> *emFraction;
  vector<float> *hcalTotOverPLead;
  vector<float> *hcal3x3OverPLead;
  vector<float> *ecalStripSumEOverPLead;
  vector<float> *bremsRecoveryEOverPLead;

  // from the Discriminator for shrinking cone
  vector<float> *theTauDiscrByLeadingTrackFinding;
  vector<float> *theTauDiscrByLeadingTrackPtCut;
  vector<float> *theTauDiscrByLeadingPionPtCut;
  vector<float> *theTauDiscrByIsolation;
  vector<float> *theTauDiscrByIsolationUsingLeadingPion;
  vector<float> *theTauDiscrByTrackIsolationUsingLeadingPion;
  vector<float> *theTauDiscrByTrackIsolation;
  vector<float> *theTauDiscrByECALIsolation;
  vector<float> *theTauDiscrByECALIsolationUsingLeadingPion;
  vector<float> *theTauDiscrAgainstMuon;
  vector<float> *theTauDiscrAgainstElectron;
  vector<float> *theTauDiscrByTaNC;
  vector<float> *theTauDiscrByTaNCfrHalfPercent;
  vector<float> *theTauDiscrByTaNCfrOnePercent;
  vector<float> *theTauDiscrByTaNCfrQuarterPercent;
  vector<float> *theTauDiscrByTaNCfrTenthPercent;

  // from the Discriminator for HPS
  vector<float> *thehpsTauDiscrByLooseElectronRejection;
  vector<float> *thehpsTauDiscrByMediumElectronRejection;
  vector<float> *thehpsTauDiscrByTightElectronRejection;
  vector<float> *thehpsTauDiscrByLooseMuonRejection;
  vector<float> *thehpsTauDiscrByTightMuonRejection;
  vector<float> *thehpsTauDiscrByDecayModeFinding;
  vector<float> *thehpsTauDiscrByVLooseIsolation;
  vector<float> *thehpsTauDiscrByLooseIsolation;
  vector<float> *thehpsTauDiscrByMediumIsolation;
  vector<float> *thehpsTauDiscrByTightIsolation;

  // from the Discriminator for HPS TaNC
  vector<float> *thehpsTancTausDiscrByLeadingTrackFinding;
  vector<float> *thehpsTancTausDiscrByLeadingTrackPtCut;
  vector<float> *thehpsTancTausDiscrByLeadingPionPtCut;
  vector<float> *thehpsTancTausDiscrByTanc;
  vector<float> *thehpsTancTausDiscrByTancRaw;
  vector<float> *thehpsTancTausDiscrByTancVLoose;
  vector<float> *thehpsTancTausDiscrByTancLoose;
  vector<float> *thehpsTancTausDiscrByTancMedium;
  vector<float> *thehpsTancTausDiscrByTancTight;
  vector<float> *thehpsTancTausDiscrByLooseElectronRejection;
  vector<float> *thehpsTancTausDiscrByMediumElectronRejection;
  vector<float> *thehpsTancTausDiscrByTightElectronRejection;
  vector<float> *thehpsTancTausDiscrByLooseMuonRejection;
  vector<float> *thehpsTancTausDiscrByTightMuonRejection;
  vector<float> *thehpsTancTausDiscrByDecayModeSelection;
  vector<float> *thehpsTancTausDiscrByVLooseIsolation;
  vector<float> *thehpsTancTausDiscrByLooseIsolation;
  vector<float> *thehpsTancTausDiscrByMediumIsolation;
  vector<float> *thehpsTancTausDiscrByTightIsolation;

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
			     edm::InputTag tauDiscrByLeadingTrackFindingTag,
			     edm::InputTag tauDiscrByLeadingTrackPtCutTag,
			     edm::InputTag tauDiscrByLeadingPionPtCutTag,
			     edm::InputTag tauDiscrByIsolationTag,
			     edm::InputTag tauDiscrByIsolationUsingLeadingPionTag,
			     edm::InputTag tauDiscrByTrackIsolationUsingLeadingPionTag,
			     edm::InputTag tauDiscrByTrackIsolationTag,
			     edm::InputTag tauDiscrByECALIsolationTag,
			     edm::InputTag tauDiscrByECALIsolationUsingLeadingPionTag, 
			     edm::InputTag tauDiscrAgainstMuonTag,
			     edm::InputTag tauDiscrAgainstElectronTag,
 			     edm::InputTag tauDiscrByTaNCTag,
 			     edm::InputTag tauDiscrByTaNCfrHalfPercentTag,
 			     edm::InputTag tauDiscrByTaNCfrOnePercentTag,
 			     edm::InputTag tauDiscrByTaNCfrQuarterPercentTag,
 			     edm::InputTag tauDiscrByTaNCfrTenthPercentTag,
			     bool dumpData=false);

  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     edm::InputTag hpsTauDiscrByLooseElectronRejectionTag,
			     edm::InputTag hpsTauDiscrByMediumElectronRejectionTag,
			     edm::InputTag hpsTauDiscrByTightElectronRejectionTag,
			     edm::InputTag hpsTauDiscrByLooseMuonRejectionTag,
			     edm::InputTag hpsTauDiscrByTightMuonRejectionTag,
			     edm::InputTag hpsTauDiscrByDecayModeFindingTag,
			     edm::InputTag hpsTauDiscrByVLooseIsolationTag,
			     edm::InputTag hpsTauDiscrByLooseIsolationTag,
			     edm::InputTag hpsTauDiscrByMediumIsolationTag,
			     edm::InputTag hpsTauDiscrByTightIsolationTag,
			     bool dumpData=false);

  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event& iEvent, const edm::EventSetup& iSetup,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     edm::InputTag hpsTancTausDiscrByLeadingTrackFindingTag,
			     edm::InputTag hpsTancTausDiscrByLeadingTrackPtCutTag,
			     edm::InputTag hpsTancTausDiscrByLeadingPionPtCutTag,
			     edm::InputTag hpsTancTausDiscrByTancTag,
			     edm::InputTag hpsTancTausDiscrByTancRawTag,
			     edm::InputTag hpsTancTausDiscrByTancVLooseTag,
			     edm::InputTag hpsTancTausDiscrByTancLooseTag,
			     edm::InputTag hpsTancTausDiscrByTancMediumTag,
			     edm::InputTag hpsTancTausDiscrByTancTightTag,
			     edm::InputTag hpsTancTausDiscrByLooseElectronRejectionTag,
			     edm::InputTag hpsTancTausDiscrByMediumElectronRejectionTag,
			     edm::InputTag hpsTancTausDiscrByTightElectronRejectionTag,
			     edm::InputTag hpsTancTausDiscrByLooseMuonRejectionTag,
			     edm::InputTag hpsTancTausDiscrByTightMuonRejectionTag,
			     edm::InputTag hpsTancTausDiscrByDecayModeSelectionTag,
			     edm::InputTag hpsTancTausDiscrByVLooseIsolationTag,
			     edm::InputTag hpsTancTausDiscrByLooseIsolationTag,
			     edm::InputTag hpsTancTausDiscrByMediumIsolationTag,
			     edm::InputTag hpsTancTausDiscrByTightIsolationTag,
			     bool dumpData=false);
private:
  
  // void writePFTauBasicInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&, const reco::PFTau *tau);
  void writePFTauBasicInfo(const reco::PFTau *cand, const edm::Event&, const edm::EventSetup&);

  void writeLeadPFCandInfo(const reco::PFTau *cand, const edm::Event&, const edm::EventSetup&);

  void writePFTauDiscInfo(edm::Handle<reco::PFTauCollection> tauCollection, int theTauJetIndex,
			  const reco::PFTauDiscriminator *tauDiscrByLeadingTrackFinding,
			  const reco::PFTauDiscriminator *tauDiscrByLeadingTrackPtCut,
			  const reco::PFTauDiscriminator *tauDiscrByLeadingPionPtCut,
			  const reco::PFTauDiscriminator *tauDiscrByIsolation,
			  const reco::PFTauDiscriminator *tauDiscrByIsolationUsingLeadingPion,
			  const reco::PFTauDiscriminator *tauDiscrByTrackIsolationUsingLeadingPion,
			  const reco::PFTauDiscriminator *tauDiscrByTrackIsolation,
			  const reco::PFTauDiscriminator *tauDiscrByECALIsolation,
			  const reco::PFTauDiscriminator *tauDiscrByECALIsolationUsingLeadingPion,
			  const reco::PFTauDiscriminator *tauDiscrAgainstMuon,
			  const reco::PFTauDiscriminator *tauDiscrAgainstElectron,
 			  const reco::PFTauDiscriminator *tauDiscrByTaNC,
 			  const reco::PFTauDiscriminator *tauDiscrByTaNCfrHalfPercent,
 			  const reco::PFTauDiscriminator *tauDiscrByTaNCfrOnePercent,
 			  const reco::PFTauDiscriminator *tauDiscrByTaNCfrQuarterPercent,
 			  const reco::PFTauDiscriminator *tauDiscrByTaNCfrTenthPercent
                          );

  void writePFTauDiscInfo(edm::Handle<reco::PFTauCollection> tauCollection, int theTauJetIndex,
			 const reco::PFTauDiscriminator *hpsTauDiscrByLooseElectronRejection,
			 const reco::PFTauDiscriminator *hpsTauDiscrByMediumElectronRejection,
			 const reco::PFTauDiscriminator *hpsTauDiscrByTightElectronRejection,
			 const reco::PFTauDiscriminator *hpsTauDiscrByLooseMuonRejection,
			 const reco::PFTauDiscriminator *hpsTauDiscrByTightMuonRejection,
			 const reco::PFTauDiscriminator *hpsTauDiscrByDecayModeFinding,
			 const reco::PFTauDiscriminator *hpsTauDiscrByVLooseIsolation,
			 const reco::PFTauDiscriminator *hpsTauDiscrByLooseIsolation,
			 const reco::PFTauDiscriminator *hpsTauDiscrByMediumIsolation,
			 const reco::PFTauDiscriminator *hpsTauDiscrByTightIsolation
			 );

  void writePFTauDiscInfo(edm::Handle<reco::PFTauCollection> tauCollection, int theTauJetIndex,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingTrackFinding,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingTrackPtCut,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByLeadingPionPtCut,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTanc,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTancRaw,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTancVLoose,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTancLoose,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTancMedium,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTancTight,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseElectronRejection,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByMediumElectronRejection,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTightElectronRejection,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseMuonRejection,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTightMuonRejection,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByDecayModeSelection,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByVLooseIsolation,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByLooseIsolation,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByMediumIsolation,
			  const reco::PFTauDiscriminator *hpsTancTausDiscrByTightIsolation
			  );

  void treePFTauBasicInfo(const std::string &colPrefix, const std::string &colSuffix);
  void treePFTauDiscInfo(const std::string &colPrefix, const std::string &colSuffix);
  void treehpsPFTauDiscInfo(const std::string &colPrefix, const std::string &colSuffix);
  void treehpsTancTausDiscInfo(const std::string &colPrefix, const std::string &colSuffix);
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
