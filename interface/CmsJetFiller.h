// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsJetFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsJetFiller_h
#define CmsJetFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

struct CmsJetFillerData : public CmsCandidateFillerData {

  std::vector<float> *emFrac, *hadFrac;
  std::vector<float> *alpha;
  std::vector<float> *combinedSecondaryVertexBJetTags, 
    *combinedSecondaryVertexMVABJetTags,
    *jetBProbabilityBJetTags,
    *jetProbabilityBJetTags,
    *simpleSecondaryVertexBJetTags,
    *softMuonBJetTags,
    *trackCountingHighPurBJetTags,
    *trackCountingHighEffBJetTags;
  std::vector<float> *uncorrEnergy;

public:
  void initialise();
  void clearTrkVectors();

};

class CmsJetFiller : public CmsCandidateFiller {

 public:

  // Constructors

  // Dump everything
  CmsJetFiller(CmsTree *, 
	       edm::InputTag jetVertexAlphaCollection,
	       int maxTracks=500, int maxMCTracks=2000, 
	       bool noOutputIfLimitsReached=false );

  // Dump  everything if fatTree is true and less informations otherwise
  CmsJetFiller(CmsTree *, 
	       edm::InputTag jetVertexAlphaCollection,
	       bool fatTree, 
	       int maxTracks=500, int maxMCTracks=2000, 
	       bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsJetFiller();

  // Modifiers
  void saveJetExtras(bool );

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
  bool saveJetExtras_;
  bool saveJetBTag_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  bool dumpUncorrEnergy_;

  std::string *trkIndexName_;

  CmsJetFillerData *privateData_;
  edm::InputTag jetVertexAlphaCollection_;

  CmsTree *cmstree;


};

#endif // CmsJetFiller_h
