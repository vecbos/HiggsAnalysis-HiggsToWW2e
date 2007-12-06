// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HiggsAnalysis/HiggsToWW2e
// Description:
//      Class CmsCandidateFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
// Original Author:  Emanuele Di Marco
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsCandidateFiller_h
#define CmsCandidateFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

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

struct CmsCandidateFillerData {
  // All the vectors that will store the stuff
  // going into the tree.
  vector<int> *charge;
  vector<float> *energy, *et, *momentum;
  vector<float> *vertexX, *vertexY, *vertexZ;
  vector<float> *theta, *eta, *phi;
  vector<float> *x, *y, *z;
  vector<float> *mass, *mt;
  vector<int> *pdgId;
  vector<int> *nDau;
  vector<int> *d1Index;
  vector<int> *d2Index;

  vector<int> *mcIndex;
  int *ncand;

public:
  void initialiseCandidate();
  void clearTrkVectorsCandidate();
};


class CmsCandidateFiller {

 public:

  // Constructors

  // Dump everything
  CmsCandidateFiller(CmsTree *, int maxTracks=500, 
		    int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  // Dump  everything if fatTree is true and less informations otherwise
  CmsCandidateFiller(CmsTree *, bool fatTree, int maxTracks=500, 
		     int maxMCTracks=2000, bool noOutputIfLimitsReached=false );

  // Destructor
  virtual ~CmsCandidateFiller();

  // Modifiers


  void saveCand(bool );
  void doMcMatch(bool );

  void setMatchMap(edm::InputTag matchMap) { matchMap_ = matchMap;};

  // Operators

  // run number and all of that --- to implement

  virtual void writeCollectionToTree(const CandidateCollection *,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  virtual void writeCollectionToTree(const CompositeCandidateCollection *,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

  virtual void writeMcIndicesToTree(const CandidateCollection *,
			    const edm::Event&, const edm::EventSetup&,
			    const CandidateCollection *,
			    const std::string &columnPrefix, const std::string &columnSuffix,
			    bool dumpData=false);


  // Add daughter list
  virtual void addDaughterList(const CandidateCollection *);

 protected:
  

  virtual void writeCandInfo(const Candidate *cand, 
		     const edm::Event&, const edm::EventSetup&);
  virtual void treeCandInfo(const std::string colPrefix, const std::string colSuffix);
  
  virtual void writeMcMatchInfo(const CandidateCollection *, 
			const edm::Event&, const edm::EventSetup&, 
			const CandidateCollection *);
  virtual void treeMcMatchInfo(const std::string colPrefix, const std::string colSuffix);

  virtual bool candOverlap(const Candidate *cand1, const Candidate *cand2);

  // Friends

  CmsCandidateFillerData *privateData_;
  edm::InputTag matchMap_;
  std::vector<const CandidateCollection*> daugCollectionList_;

  CmsTree *cmstree;

  bool saveCand_;
  bool doMcMatch_;

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;
};

#endif // CmsCandidateFiller_h
