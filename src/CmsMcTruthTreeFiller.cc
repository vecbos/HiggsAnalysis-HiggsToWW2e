//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsMcTruthTreeFiller
//      Simple class for dumping MC truth info to an ntuple. 
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  21 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "CLHEP/HepMC/GenParticle.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMcTruthTreeFiller.h"


#include <string>

struct CmsMcTruthTreeFillerData {
  CmsTree *cmstree;
  
};

//----------------
// Constructors --
//----------------
CmsMcTruthTreeFiller::CmsMcTruthTreeFiller(CmsTree *cmstree):
  privateData_(new CmsMcTruthTreeFillerData)
{
  privateData_->cmstree=cmstree; 
}

//--------------
// Destructor --
//--------------
CmsMcTruthTreeFiller::~CmsMcTruthTreeFiller() {
  delete privateData_;
}

//-------------
// Methods   --
//-------------
void CmsMcTruthTreeFiller::writeCollectionToTree(edm::InputTag mcTruthCollection, 
						 const edm::Event& iEvent, int range) {

  // prepared for 200, when CandidateCollection will be substituted by GenParticleCollection
  //  edm::Handle< reco::GenParticleCollection > genParticleHandle;
  edm::Handle< reco::CandidateCollection > genParticleHandle;
  try { iEvent.getByLabel(mcTruthCollection, genParticleHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("HWWTreeDumper") << "Can't get MC Truth Collection: " << mcTruthCollection; }
  //  const reco::GenParticleCollection *genParticleCollection = genParticleHandle.product();
  const reco::CandidateCollection *genParticleCollection = genParticleHandle.product(); 

  vector<float> pMC,massMC,thetaMC,etaMC,phiMC,energyMC;
  vector<float> xMC,yMC,zMC;
  vector<int> idMC,mothMC,nDauMC,statusMC;

  //  reco::GenParticleCollection::const_iterator genPart;
  reco::CandidateCollection::const_iterator genPart;
  for(genPart=genParticleCollection->begin(); genPart!=genParticleCollection->end(); genPart++) {
    const reco::Candidate & cand = *genPart;
    if((int)pMC.size()>range) break;
    // write only stables or documentation particles
    if(cand.status() == 1 || cand.status() == 3) { 
      pMC.push_back(cand.p());
      massMC.push_back(cand.mass());
      thetaMC.push_back(cand.theta());
      etaMC.push_back(cand.eta());
      phiMC.push_back(cand.phi());
      energyMC.push_back(cand.energy());
      idMC.push_back(cand.pdgId());
      nDauMC.push_back(cand.numberOfDaughters());
      statusMC.push_back(cand.status());
    
      int indMom=-1;
      int idx=0;
      //    reco::GenParticleCollection::const_iterator candIter;
      reco::CandidateCollection::const_iterator candIter;
      for(candIter=genParticleCollection->begin(); candIter!=genParticleCollection->end(); candIter++) {
	if(cand.status() == 1 || cand.status() == 3) {
	  const reco::Candidate *mom = cand.mother();
	  if(&(*candIter)==&(*mom)) {
	    indMom=idx;
	    break;
	  }
	  idx++;
	}
      }
      mothMC.push_back(indMom);
    
      // decay vertex
      xMC.push_back(cand.vx());
      yMC.push_back(cand.vy());
      zMC.push_back(cand.vz());
    }
  }
  
  std::string theName = "Mc";
  std::string indName = "n"+theName;
  
  privateData_->cmstree->column(indName.c_str(),pMC.size(),0,theName.c_str());
  privateData_->cmstree->column(("p"+theName).c_str(), pMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("mass"+theName).c_str(), massMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("theta"+theName).c_str(), thetaMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("eta"+theName).c_str(), etaMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("phi"+theName).c_str(), phiMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("energy"+theName).c_str(), energyMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("id"+theName).c_str(), idMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("moth"+theName).c_str(), mothMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("nDau"+theName).c_str(), nDauMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("status"+theName).c_str(), statusMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("x"+theName).c_str(), xMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("y"+theName).c_str(), yMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("z"+theName).c_str(), zMC, indName.c_str(), 0, theName.c_str());
    
}
  
