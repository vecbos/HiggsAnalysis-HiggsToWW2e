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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

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
  saveLHE_ = false;
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
void CmsMcTruthTreeFiller::writeCollectionToTree(edm::InputTag mcTruthCollection, std::vector<std::string>* lheComments,
						 const edm::Event& iEvent, int range, bool firstEvent) {

  edm::Handle< reco::GenParticleCollection > genParticleHandle;
  try { iEvent.getByLabel(mcTruthCollection, genParticleHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("HWWTreeDumper") << "Can't get MC Truth Collection: " << mcTruthCollection; }
  const reco::GenParticleCollection *genParticleCollection = genParticleHandle.product();

  vector<float> pMC,massMC,thetaMC,etaMC,phiMC,energyMC,vxMC,vyMC,vzMC;
  vector<float> xMC,yMC,zMC;
  vector<int> idMC,mothMC,nDauMC,statusMC;

  reco::GenParticleCollection::const_iterator genPart;
  for(genPart=genParticleCollection->begin(); genPart!=genParticleCollection->end(); genPart++) {
    const reco::Candidate & cand = *genPart;
    if((int)pMC.size()>range) break;
    // write only stables or documentation particles
    // if(cand.status() == 1 || cand.status() == 3) { 
    pMC.push_back(cand.p());
    massMC.push_back(cand.mass());
    thetaMC.push_back(cand.theta());
    etaMC.push_back(cand.eta());
    phiMC.push_back(cand.phi());
    energyMC.push_back(cand.energy());
    vxMC.push_back(cand.vx());
    vyMC.push_back(cand.vy());
    vzMC.push_back(cand.vz());
    idMC.push_back(cand.pdgId());
    nDauMC.push_back(cand.numberOfDaughters());
    statusMC.push_back(cand.status());
    
    int indMom=-1;
    int idx=0;
    reco::GenParticleCollection::const_iterator candIter;
    for(candIter=genParticleCollection->begin(); candIter!=genParticleCollection->end(); candIter++) {
      // if(cand.status() == 1 || cand.status() == 3) {
      const reco::Candidate *mom = cand.mother();
      if(&(*candIter)==&(*mom)) {
	indMom=idx;
	break;
      }
      idx++;
      //}
    }
    mothMC.push_back(indMom);
    
    // decay vertex
    xMC.push_back(cand.vx());
    yMC.push_back(cand.vy());
    zMC.push_back(cand.vz());
    //    }
  }

  // LHE Event
  // to be written in the ntuple
  if(saveLHE_) {
    lheComments->clear();
    edm::Handle<LHEEventProduct> product;
    iEvent.getByLabel("source", product);
    std::vector<std::string>::const_iterator c_begin = product->comments_begin();
    std::vector<std::string>::const_iterator c_end = product->comments_end();
    for( std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
      lheComments->push_back(*cit);
    }
  }
  
  std::string theName = "Mc";
  std::string indName = "n"+theName;
  
  privateData_->cmstree->column(indName.c_str(),(int)pMC.size(),0,theName.c_str());
  privateData_->cmstree->column(("p"+theName).c_str(), pMC, indName.c_str(), 0, theName.c_str());
  //  privateData_->cmstree->column(("mass"+theName).c_str(), massMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("theta"+theName).c_str(), thetaMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("eta"+theName).c_str(), etaMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("phi"+theName).c_str(), phiMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("energy"+theName).c_str(), energyMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("vx"+theName).c_str(), vxMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("vy"+theName).c_str(), vyMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("vz"+theName).c_str(), vzMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("id"+theName).c_str(), idMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("moth"+theName).c_str(), mothMC, indName.c_str(), 0, theName.c_str());
  //  privateData_->cmstree->column(("nDau"+theName).c_str(), nDauMC, indName.c_str(), 0, theName.c_str());
  privateData_->cmstree->column(("status"+theName).c_str(), statusMC, indName.c_str(), 0, theName.c_str());
  //  privateData_->cmstree->column(("x"+theName).c_str(), xMC, indName.c_str(), 0, theName.c_str());
  //  privateData_->cmstree->column(("y"+theName).c_str(), yMC, indName.c_str(), 0, theName.c_str());
  //  privateData_->cmstree->column(("z"+theName).c_str(), zMC, indName.c_str(), 0, theName.c_str());

  if(firstEvent && saveLHE_) privateData_->cmstree->getTree()->Branch("commentLHE", &(*lheComments));

}
  
