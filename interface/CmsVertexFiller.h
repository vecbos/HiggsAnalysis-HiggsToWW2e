// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Package:    
//      HtoWWTreeDumper
// Description:
//      Class CmsGenInfoFiller
//      Simple class for dumping RECO (or AOD) contents to an ntuple
//      
// Original Author:  Alessio Ghezzi, Pietro Govoni
//         Created:  Fri Apr  6 18:05:34 CEST 2007
//
//-----------------------------------------------------------------------

#ifndef CmsVertexFiller_h
#define CmsVertexFiller_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include <TTree.h>

using namespace cms;
using namespace edm;
using namespace reco;

struct CmsVertexFillerData : public CmsCandidateFillerData {
  std::vector<float> *PVx, *PVy, *PVz;
  std::vector<float> *PVErrx, *PVErry, *PVErrz;
  std::vector<float> *SumPt, *ndof, *chi2, *normChi2, *rho;
  std::vector<float> *pxChMet, *pyChMet, *pzChMet;
  std::vector<int> *isFake, *isValid, *trackSize;

public:
  void initialise();
  void clearTrkVectors();
};

struct CmsVertexTracksFillerData : public CmsCandidateFillerData {
  //  std::vector<float> *pt, *eta, *phi, *d0, *Chi2, *pterr;
  //  std::vector<int> *charge, *VtxIndex, *recHitsSize, *nvalhits;
  
public:
  void initialise();
  void clearTrkVectors();
};

class CmsVertexFiller : public CmsCandidateFiller {

public:

  //! Constructors
  CmsVertexFiller(CmsTree *, bool,
		  int maxTracks=500, int maxMCTracks=2000, 
		  bool noOutputIfLimitsReached=false );

  //! Destructor
  virtual ~CmsVertexFiller();

  void writeCollectionToTree (edm::InputTag vtxcollectionTag,
  			      const edm::Event& iEvent, 
			      const edm::EventSetup& iSetup,
			      const std::string &columnPrefix,
			      const std::string &columnSuffix);

  //! set the charged met collection associated to the PVs
  void setChargedMet( edm::InputTag ChargedMets ) { ChargedMets_ = ChargedMets; }

private:
  void treeVertexInfo(const std::string &colPrefix,
		      const std::string &colSuffix);

  CmsVertexFillerData *privateData_;
  //    CmsVertexTracksFillerData *privateDatat_;
  edm::InputTag ChargedMets_;
  std::string *trkIndexName_;
  CmsTree *cmstree;

};

#endif // CmsVertexFiller_h
