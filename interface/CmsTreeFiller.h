// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsTreeFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsTreeFiller_h
#define CmsTreeFiller_h

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

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "CLHEP/HepMC/GenEvent.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;

class CmsTreeFillerData;

class CmsTreeFiller {

 public:

  // Constructors

  // Dump everything
  CmsTreeFiller(CmsTree *, int maxTracks=500, int maxNeutrals=500,
		    int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  // Dump  everything if fatTree is true and less informations otherwise
  CmsTreeFiller(CmsTree *, bool fatTree, int maxTracks=500, int maxNeutrals=500,
		    int maxMCTracks=2000, bool noOutputIfLimitsReached=false );


  // Destructor
  virtual ~CmsTreeFiller();

  // Modifiers

  void saveTrk(bool );
  void saveEcal(bool );
  void saveHcal(bool );
  void saveDT(bool );
  void saveCSC(bool );
  void saveRPC(bool );
  void saveJetAlpha(bool );

  void saveCand(bool );

  void saveFatTrk(bool );
  void saveFatEcal(bool );
  void saveFatHcal(bool );
  void saveFatDT(bool );
  void saveFatCSC(bool );
  void saveFatRPC(bool );

  void saveEleID(bool );

  void doMcMatch(bool );
  void setMatchMap(edm::InputTag matchMap) { matchMap_ = matchMap;};

  // Operators

  // run number and all of that --- to implement

  void writeCollectionToTree(const CandidateCollection *,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);
  void writeMcIndicesToTree(const CandidateCollection *,
			    const edm::Event&, const edm::EventSetup&,
			    const CandidateCollection *,
			    const std::string &columnPrefix, const std::string &columnSuffix,
			    bool dumpData=false);

 private:
  
  void writeTrkInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&, GsfTrackRef trkRef);
  void treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeEcalInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&, SuperClusterRef);
  void treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeHcalInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeHcalInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeMuonInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeMuonInfo(const std::string &colPrefix, const std::string &colSuffix);

  void writeCandInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeCandInfo(const std::string colPrefix, const std::string colSuffix);
  
  void writeMcMatchInfo(const CandidateCollection *, const edm::Event&, const edm::EventSetup&, const CandidateCollection *);
  void treeMcMatchInfo(const std::string colPrefix, const std::string colSuffix);

  void writeJetAlphaInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeJetAlphaInfo(const std::string &colPrefix, const std::string &colSuffix);

  void initialise();
  void clearTrkVectors();

  // Friends

  CmsTreeFillerData *privateData_;
  edm::InputTag matchMap_;

};

#endif // CmsTreeFiller_h
