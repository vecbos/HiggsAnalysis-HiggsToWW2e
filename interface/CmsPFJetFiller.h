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
#include "CMGTools/External/interface/PileupJetIdAlgo.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"

#include <vector>
#include <TTree.h>


struct QGLikelihoodVars {

  float ptD;
  float rmsCand;
  float axis1;
  float axis2;
  float pull;
  float r_ch;
  float tana;

};


struct CmsPFJetFillerData {

  // candidate. Repeated because we want to store the corrected jets energies
  vector<int> *charge;
  vector<float> *energy, *et, *momentum;
  vector<float> *vertexX, *vertexY, *vertexZ;
  vector<float> *theta, *eta, *phi;
  vector<float> *x, *y, *z;
  vector<float> *uncorrx, *uncorry, *uncorrz;
  int *ncand;

  std::vector<float> *chargedHadronEnergy, *neutralHadronEnergy, *photonEnergy, *electronEnergy, *muonEnergy, 
    *HFHadronEnergy, *HFEMEnergy;
  std::vector<int> *chargedHadronMultiplicity, *neutralHadronMultiplicity, 
    *photonMultiplicity, *electronMultiplicity, *muonMultiplicity,
    *HFHadronMultiplicity, *HFEMMultiplicity;
  std::vector<float> *uncorrEnergy, *area;
  std::vector<float> *ptD, *rmsCand, *axis1, *axis2, *pull, *r_ch, *tana;
  std::vector<float> *combinedSecondaryVertexBJetTags, 
    *combinedSecondaryVertexMVABJetTags, 
    *jetBProbabilityBJetTags, 
    *jetProbabilityBJetTags,
    *simpleSecondaryVertexHighEffBJetTags,
    *simpleSecondaryVertexHighPurBJetTags,
    *trackCountingHighPurBJetTags,
    *trackCountingHighEffBJetTags,
    *trackCountingVeryHighEffBJetTags;
  std::vector<float> *weightedDz1, *weightedDz2;

  std::vector<float> *betastar, *rmsCandsHand;
  std::vector<float> *jetIdMvaSimple, *jetIdMvaFull, *jetIdMvaPhilV1;
  std::vector<float> *nChargedIdMva, *nNeutralsIdMva;
  std::vector<float> *dZIdMva;
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
  void setVertexCollection ( edm::InputTag verticesTag ) { _verticesTag = verticesTag; }
  void setjetMVAAlgos( std::vector<PileupJetIdAlgo* > jetId_algos ) { _jetId_algos = jetId_algos; }
  
 private:
  
  void writeJetInfo(const reco::Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeCandInfo(const std::string colPrefix, const std::string colSuffix);
  void treeJetInfo(const std::string &colPrefix, const std::string &colSuffix);
  float DzCorrected(reco::TrackRef trk, reco::Vertex vtx);

  // Friends
  bool saveJetBTag_;

  edm::ParameterSet BTagCollections_;
  std::vector<edm::InputTag> PFjetMvaIdCollection_;
  edm::InputTag _verticesTag;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;
  std::vector<PileupJetIdAlgo* > _jetId_algos;

  CmsPFJetFillerData *privateData_;
  std::string m_jcs;

  CmsTree *cmstree;

  reco::Vertex bestPrimaryVertex_;

};

#endif // CmsPFJetFiller_h
