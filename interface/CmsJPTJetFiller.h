// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsJPTJetFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Mon Sep  29 18:05:34 CEST 2008
//
//-----------------------------------------------------------------------

#ifndef CmsJPTJetFiller_h
#define CmsJPTJetFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include <vector>
#include <TTree.h>

struct CmsJPTJetFillerData : public CmsCandidateFillerData {

  std::vector<float> *chargedHadronEnergy, *neutralHadronEnergy, *chargedEmEnergy, *neutralEmEnergy;
  std::vector<float> *chargedMultiplicity, *muonMultiplicity, *elecMultiplicity;
  std::vector<float> *ZSPCor, *Eta2momtr, *Phi2momtr;
  std::vector<float> *uncorrEnergy;
  std::vector<float> *combinedSecondaryVertexBJetTags, 
    *combinedSecondaryVertexMVABJetTags,
    *jetBProbabilityBJetTags,
    *jetProbabilityBJetTags,
    *simpleSecondaryVertexBJetTags,
    *softMuonBJetTags,
    *softMuonByIP3dBJetTags,
    *softMuonByPtBJetTags,
    *softElectronBJetTags,
    *softElectronByIP3dBJetTags,
    *softElectronByPtBJetTags,
    *trackCountingHighPurBJetTags,
    *trackCountingHighEffBJetTags;

public:
  void initialise();
  void clearTrkVectors();

};

class CmsJPTJetFiller : public CmsCandidateFiller {

 public:

  // Constructors

  // Dump everything
  CmsJPTJetFiller(CmsTree *, 
	       int maxTracks=500, int maxMCTracks=2000, 
	       bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsJPTJetFiller();

  // Modifiers
  void saveJetBTag(bool );

  // Operators

  // run number and all of that --- to implement

  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false,
                             edm::InputTag uncorrectedCollectionTag=edm::InputTag("",""));

 private:
  
  void writeJetInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeJetInfo(const std::string &colPrefix, const std::string &colSuffix);

  // Friends
  bool saveJetBTag_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  bool dumpUncorrEnergy_;

  std::string *trkIndexName_;

  CmsJPTJetFillerData *privateData_;

  CmsTree *cmstree;

};

#endif // CmsJPTJetFiller_h
