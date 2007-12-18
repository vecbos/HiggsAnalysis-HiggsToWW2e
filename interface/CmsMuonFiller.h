// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      ZjetsAnalysis/ZllProducer
// Description:
//      Class CmsMuonFiller
//      Simple class for dumping RECO (or AOD) contents to a ROOT tree
//      
//-----------------------------------------------------------------------

#ifndef CmsMuonFiller_h
#define CmsMuonFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "CLHEP/HepMC/GenEvent.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;

struct CmsMuonFillerData : public CmsCandidateFillerData {

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

  // Operators

  // run number and all of that --- to implement

  void writeCollectionToTree(const edm::View<reco::Candidate> *,
			     const edm::Event&, const edm::EventSetup&,
			     const std::string &columnPrefix, const std::string &columnSuffix,
			     bool dumpData=false);

 private:
  
  void writeMuonInfo(const Candidate *cand, const edm::Event&, const edm::EventSetup&);
  void treeMuonInfo(const std::string &colPrefix, const std::string &colSuffix);

  bool hitLimitsMeansNoOutput_;
  int maxTracks_;
  int maxMCTracks_;

  std::string *trkIndexName_;

  CmsMuonFillerData *privateData_;
  edm::InputTag matchMap_;

  CmsTree *cmstree;

};

#endif // CmsMuonFiller_h
