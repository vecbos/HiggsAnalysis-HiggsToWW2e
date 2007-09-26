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

struct CmsJetFillerData : public CmsCandidateFillerData {

  std::vector<float> *emFrac, *hadFrac;
  std::vector<float> *alpha;

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

  // Operators

  // run number and all of that --- to implement

  void writeCollectionToTree(const CandidateCollection *,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

 private:
  
  void writeJetInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeJetInfo(const std::string &colPrefix, const std::string &colSuffix);

  // Friends
  bool saveJetExtras_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsJetFillerData *privateData_;
  edm::InputTag jetVertexAlphaCollection_;

  CmsTree *cmstree;

};

#endif // CmsJetFiller_h
