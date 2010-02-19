//---------------------------------------------------------------------------
//
// Description:
//       Package:   HtoWWTreeDumper
//       Class:     CmsTreeFiller
//
// Original Author:  Emanuele Di Marco
//         Created:  Mon Apr  9 11:01:00 CEST 2007
//-----------------------------------------------------------------------

// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementGsfTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFlowElectronFiller.h"

#include <TTree.h>
#include <string>

using namespace edm;
using namespace reco;


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsPFlowElectronFiller::CmsPFlowElectronFiller(CmsTree *cmsTree, 
					       int maxTracks, int maxMCTracks,
					       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFlowElectronFillerData)
{
  cmstree=cmsTree;

  savePFEleTrk_=true;     
  savePFEleBasic_=true;
  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}

CmsPFlowElectronFiller::CmsPFlowElectronFiller(CmsTree *cmsTree, bool fatTree, 
					       int maxTracks, int maxMCTracks,
					       bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsPFlowElectronFillerData)
{
  cmstree=cmsTree;

  savePFEleTrk_=true;     
  savePFEleBasic_=true;
  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();

}


//--------------
// Destructor --
//--------------

CmsPFlowElectronFiller::~CmsPFlowElectronFiller() {

  delete privateData_->ncand;
  
  delete privateData_->trackIndex;
  delete privateData_->gsfTrackIndex;

  delete privateData_->MvaOutput;
  delete privateData_->PS1Energy;
  delete privateData_->PS2Energy;
  delete privateData_->EcalEnergy;  
  delete privateData_->EcalElectronEnergy;

  delete privateData_->EoPMode;
  delete privateData_->EoPoutMode;
  delete privateData_->EBremoDeltaP;
  delete privateData_->deltaEtaAtCalo;
  delete privateData_->deltaPhiAtCalo;
  delete privateData_->hOverHE;
  delete privateData_->hOverP;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsPFlowElectronFiller::writeCollectionToTree(edm::InputTag collectionTag,
						   const edm::Event& iEvent, const edm::EventSetup& iSetup,
						   const std::string &columnPrefix, const std::string &columnSuffix,
						   bool dumpData) {
  
  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFlowElectronFiller") << "Can't get electron candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();
  
  privateData_->clearTrkVectors();
  
  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFlowElectronFiller") << "Track length " << collection->size() 
					     << " is too long for declared max length for tree "
					     << maxTracks_ << " and no output flag is set."
					     << " No tracks written to tuple for this event ";
      return;
    }
    
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsPFlowElectronFiller") << "Track length " << collection->size() 
					     << " is too long for declared max length for tree "
					     << maxTracks_ 
					     << ". Collection will be truncated ";
    }
    
    *(privateData_->ncand) = collection->size();
    
    
    for(int index = 0; index < (int)collection->size(); index++) {
      
      // fill basic kinematics
      const Candidate *cand = &(collection->at(index));
      const PFCandidateRef pflowCandRef = collection->refAt(index).castTo<PFCandidateRef>();
      if(saveCand_) writeCandInfo(cand,iEvent,iSetup);
      
      if ( !(pflowCandRef.isNull()) ) {
	
	// basic PF infos
        if(savePFEleBasic_) writePFEleBasicInfo(pflowCandRef);
	
	// tracker based infos
	GsfTrackRef gsfRef  = pflowCandRef->gsfTrackRef();
	TrackRef kfTrackRef = pflowCandRef->trackRef();
	if(savePFEleTrk_) writePFEleTrkInfo(gsfRef,kfTrackRef);
      }
      else {
	edm::LogWarning("CmsPFlowElectronFiller") << "Warning! The collection seems to be not made by "
						  << "pflow candidates electrons, electron-specific infos will be set to default.";
      }
    }
    
  }
  else {
    *(privateData_->ncand) = 0;
  }
  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  int blockSize = (collection) ? collection->size() : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(saveCand_)       treeCandInfo(columnPrefix,columnSuffix);
  if(savePFEleBasic_) treePFEleBasicInfo(columnPrefix,columnSuffix);
  if(savePFEleTrk_)   treePFEleTrkInfo(columnPrefix,columnSuffix);
  if(dumpData) cmstree->dumpData();

}

void CmsPFlowElectronFiller::writePFEleTrkInfo(reco::GsfTrackRef gsfRef, reco::TrackRef kfTrackRef ) {
  
  if(gsfRef.isNonnull()) {
    privateData_->gsfTrackIndex->push_back(gsfRef.key());
  }
  else {
    privateData_->gsfTrackIndex->push_back( -1 );
  }
  
  if(kfTrackRef.isNonnull()) {
    privateData_->trackIndex->push_back(kfTrackRef.key());
  } else {
    privateData_->trackIndex->push_back( -1 );
  }
}

void CmsPFlowElectronFiller::writePFEleBasicInfo(const reco::PFCandidateRef pflowCandRef) {

  float DEtaGsfEcalClust, DPhiGsfEcalClust, sigmaEtaEta, HOverHE, HOverPin;
  DEtaGsfEcalClust = DPhiGsfEcalClust = sigmaEtaEta = HOverHE = HOverPin = -999;
  
  float deta_gsfecal, dphi_gsfecal;
  deta_gsfecal = dphi_gsfecal = -999;

  float Ene_ecalgsf, Ene_ecalbrem, Ene_tot, Eout_gsf;
  Ene_ecalgsf = Ene_ecalbrem = Ene_tot = Eout_gsf = 0.;

  float EtotPinMode, EGsfPoutMode, EtotBremPinPoutMode;
  EtotPinMode = EGsfPoutMode = EtotBremPinPoutMode = -1;

  if (pflowCandRef.isNonnull()) {
    
    privateData_->MvaOutput->push_back(pflowCandRef->mva_e_pi());
    privateData_->PS1Energy->push_back(pflowCandRef->pS1Energy());
    privateData_->PS2Energy->push_back(pflowCandRef->pS2Energy());
    privateData_->EcalEnergy->push_back(pflowCandRef->ecalEnergy());
    privateData_->EcalElectronEnergy->push_back(pflowCandRef->rawEcalEnergy());

    // the electron ID variables like dEta/dPhi are done with the PF-linked cluster,
    // not the seed. So we need to loop over the elements and find the linked cluster
    const PFCandidate::ElementsInBlocks& theElements = pflowCandRef->elementsInBlocks();
    typedef PFCandidate::ElementsInBlocks::const_iterator IEB; 

    int firstEcal = 0;

    // the closest PF cluster is attached as the first cluster to the PF electron
    for (IEB ieb=theElements.begin(); ieb<theElements.end(); ++ieb) {

      const PFBlock& block = *(ieb->first);
      const PFBlockElement& pfbe = block.elements()[ieb->second];

      PFBlockElement::Type typeassoc = pfbe.type();

      if(typeassoc == reco::PFBlockElement::ECAL) {
        PFClusterRef clusterRef = pfbe.clusterRef();

        if(firstEcal == 0) {
          GsfTrackRef gsfRef  = pflowCandRef->gsfTrackRef();

          deta_gsfecal = clusterRef->position().eta() - pflowCandRef->positionAtECALEntrance().eta();
          dphi_gsfecal = clusterRef->position().phi() - pflowCandRef->positionAtECALEntrance().phi();
          
          vector< const reco::PFCluster * > pfClust_vec(0);
          pfClust_vec.clear();
          pfClust_vec.push_back(&(*clusterRef));

//           PFClusterWidthAlgo pfwidth(pfClust_vec);
//           sigmaEtaEta = pfwidth.pflowSigmaEtaEta();

        } 

        firstEcal++;
       }

      if(typeassoc == reco::PFBlockElement::GSF) {

        // take the mode pout from the PFBlockElement
        const reco::PFBlockElementGsfTrack *  GsfEl =  dynamic_cast<const reco::PFBlockElementGsfTrack*>(&pfbe); 
        Eout_gsf = GsfEl->Pout().t();

      }

    }
    
    // Eraw = energy of cluster associated to GSF track
    // E = all the supercluster energy
    Ene_ecalgsf = pflowCandRef->rawEcalEnergy();
    Ene_tot = pflowCandRef->ecalEnergy();
    Ene_ecalbrem = Ene_tot - Ene_ecalgsf;

    float Ein_gsf   = 0.;
    GsfTrackRef gsfRef  = pflowCandRef->gsfTrackRef();
    if (gsfRef.isNonnull()) {
      float m_el=0.00051;
      Ein_gsf =sqrt(gsfRef->pMode()*
                    gsfRef->pMode()+m_el*m_el);
    }

    float Ene_hcalgsf = pflowCandRef->hcalEnergy();

    EtotPinMode        =  Ene_tot / Ein_gsf; 
    EGsfPoutMode       =  Ene_ecalgsf/Eout_gsf;
    EtotBremPinPoutMode = Ene_ecalbrem /(Ein_gsf - Eout_gsf);
    DEtaGsfEcalClust  = deta_gsfecal;
    DPhiGsfEcalClust = dphi_gsfecal;
    HOverHE = Ene_hcalgsf/(Ene_hcalgsf + Ene_ecalgsf);
    HOverPin = Ene_hcalgsf / Ein_gsf;
  
    privateData_->EoPMode->push_back(EtotPinMode);
    privateData_->EoPoutMode->push_back(EGsfPoutMode);
    privateData_->EBremoDeltaP->push_back(EtotBremPinPoutMode);
    privateData_->deltaEtaAtCalo->push_back(DEtaGsfEcalClust);
    privateData_->deltaPhiAtCalo->push_back(DPhiGsfEcalClust);
    privateData_->hOverHE->push_back(HOverHE);
    privateData_->hOverP->push_back(HOverPin);
    
  } else {  
    
    privateData_->MvaOutput->push_back(-1.);
    privateData_->PS1Energy->push_back(-1.);
    privateData_->PS2Energy->push_back(-1.);
    privateData_->EcalEnergy->push_back(-1.);
    privateData_->EcalElectronEnergy->push_back(-1.);
  }
}



void CmsPFlowElectronFiller::treePFEleTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"gsfTrackIndex"+colSuffix).c_str(), *privateData_->gsfTrackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(), *privateData_->trackIndex, nCandString.c_str(), 0, "Reco");
}

void CmsPFlowElectronFiller::treePFEleBasicInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"MvaOutput"+colSuffix).c_str(),  *privateData_->MvaOutput, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PS1Energy"+colSuffix).c_str(),  *privateData_->PS1Energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PS2Energy"+colSuffix).c_str(),  *privateData_->PS2Energy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"EcalEnergy"+colSuffix).c_str(),  *privateData_->EcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"EcalElectronEnergy"+colSuffix).c_str(),  *privateData_->EcalElectronEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"EoPMode"+colSuffix).c_str(),  *privateData_->EoPMode, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"EoPoutMode"+colSuffix).c_str(),  *privateData_->EoPoutMode, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"EBremoDeltaP"+colSuffix).c_str(),  *privateData_->EBremoDeltaP, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaEtaAtCalo"+colSuffix).c_str(),  *privateData_->deltaEtaAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiAtCalo"+colSuffix).c_str(),  *privateData_->deltaPhiAtCalo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hOverHE"+colSuffix).c_str(),  *privateData_->hOverHE, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hOverP"+colSuffix).c_str(),  *privateData_->hOverP, nCandString.c_str(), 0, "Reco");

}



void CmsPFlowElectronFillerData::initialise() {  
  initialiseCandidate();

  // track infos
  trackIndex = new vector<int>;
  gsfTrackIndex = new vector<int>;  

  // basic pf candidates info
  MvaOutput = new vector<float>;
  PS1Energy = new vector<float>;
  PS2Energy = new vector<float>;
  EcalEnergy = new vector<float>;
  EcalElectronEnergy = new vector<float>;
  EoPMode = new vector<float>;  
  EoPoutMode = new vector<float>;
  EBremoDeltaP = new vector<float>;
  deltaEtaAtCalo = new vector<float>;
  deltaPhiAtCalo = new vector<float>;
  hOverHE = new vector<float>;
  hOverP = new vector<float>;

}
void CmsPFlowElectronFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  trackIndex->clear();
  gsfTrackIndex->clear();

  // basic pf candidates info
  MvaOutput->clear();
  PS1Energy->clear();
  PS2Energy->clear();
  EcalEnergy->clear();
  EcalElectronEnergy->clear();
  EoPMode->clear();
  EoPoutMode->clear();
  EBremoDeltaP->clear();
  deltaEtaAtCalo->clear();
  deltaPhiAtCalo->clear();
  hOverHE->clear();
  hOverP->clear();

}
