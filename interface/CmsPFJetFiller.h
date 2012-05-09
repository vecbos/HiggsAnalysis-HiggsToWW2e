// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsPFJetFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Mon Sep  29 18:05:34 CEST 2008
//
//-----------------------------------------------------------------------

#ifndef CmsPFJetFiller_h
#define CmsPFJetFiller_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include <vector>
#include <TTree.h>


struct QGLikelihoodVars {

  float ptD;
  float rmsCand;

};


struct CmsPFJetFillerData : public CmsCandidateFillerData {

  std::vector<float> *chargedHadronEnergy, *neutralHadronEnergy, *photonEnergy, *electronEnergy, *muonEnergy, 
    *HFHadronEnergy, *HFEMEnergy;
  std::vector<int> *chargedHadronMultiplicity, *neutralHadronMultiplicity, 
    *photonMultiplicity, *electronMultiplicity, *muonMultiplicity,
    *HFHadronMultiplicity, *HFEMMultiplicity;
  std::vector<float> *uncorrEnergy, *area;
  std::vector<float> *ptD, *rmsCand;
  std::vector<float> *combinedSecondaryVertexBJetTags, 
    *simpleSecondaryVertexHighEffBJetTags,
    *simpleSecondaryVertexHighPurBJetTags,
    *trackCountingHighPurBJetTags,
    *trackCountingHighEffBJetTags,
    *trackCountingVeryHighEffBJetTags;
  std::vector<float> *weightedDz1, *weightedDz2;

  std::vector<float> *betastar, *rmsCandsHand;
  std::vector<float> *jetIdMva;
  std::vector<float> *nChargedIdMva, *nNeutralsIdMva;
  std::vector<float> *dZIdMva, *nParticlesIdMva;
  std::vector<float> *dR2MeanIdMva, *dRMeanIdMva;
  std::vector<float> *frac01IdMva, *frac02IdMva, *frac03IdMva, *frac04IdMva, *frac05IdMva; 
  std::vector<float> *betaIdMva, *betastarIdMva;
  std::vector<float> *betastarclassicIdMva;

  // for backward compatibility with existing trees
  std::vector<float> *chargedEmEnergy, *neutralEmEnergy;

public:
  void initialise();
  void clearTrkVectors();

};

class CmsPFJetFiller : public CmsCandidateFiller {

 public:

  // Constructors

  // Dump everything
  CmsPFJetFiller(CmsTree *, 
	       int maxTracks=500, int maxMCTracks=2000, 
	       bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsPFJetFiller();

  // Modifiers
  void saveJetBTag(bool );

  void setBTags(edm::ParameterSet btagcollections);
  void setJetCorrectionService(std::string jcs) { m_jcs = jcs; }

  // Operators

  // run number and all of that --- to implement

  void writeCollectionToTree(edm::InputTag collectionTag,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  void setMvaId( std::vector<edm::InputTag> PFjetMvaIdCollection) { PFjetMvaIdCollection_ = PFjetMvaIdCollection; }

 private:
  
  void writeJetInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeJetInfo(const std::string &colPrefix, const std::string &colSuffix);
  float DzCorrected(reco::TrackRef trk, reco::Vertex vtx);

  // Friends
  bool saveJetBTag_;

  edm::ParameterSet BTagCollections_;
  std::vector<edm::InputTag> PFjetMvaIdCollection_;

  std::vector<edm::ValueMap<float> > mvas_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsPFJetFillerData *privateData_;
  std::string m_jcs;

  CmsTree *cmstree;

  reco::Vertex bestPrimaryVertex_;

};

#endif // CmsPFJetFiller_h
