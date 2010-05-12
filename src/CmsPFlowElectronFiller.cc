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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFlowElectronFiller.h"

#include <TTree.h>
#include <string>

using namespace edm;
using namespace reco;
using namespace std;
using namespace math;
using namespace reco::isodeposit;


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
  
  savePFEleTrk_    = true;     
  savePFEleBasic_  = true;
  savePFEleIsoDep_ = true;
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

  savePFEleTrk_    = true;     
  savePFEleBasic_  = true;
  savePFEleIsoDep_ = true;
  trkIndexName_    = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsPFlowElectronFiller::~CmsPFlowElectronFiller() {

  delete privateData_->trackIndex;
  delete privateData_->gsfTrackIndex;

  delete privateData_->MvaOutput;
  delete privateData_->PS1Energy;
  delete privateData_->PS2Energy;
  delete privateData_->EcalEnergy;  
  delete privateData_->RawEcalEnergy;
  delete privateData_->HcalEnergy;  
  delete privateData_->RawHcalEnergy;  
  delete privateData_->PositionAtEcalX;
  delete privateData_->PositionAtEcalY;
  delete privateData_->PositionAtEcalZ;
  
  delete privateData_->chIso03veto;
  delete privateData_->chIso04veto;
  delete privateData_->chIso05veto;
  delete privateData_->chIso03noVeto;
  delete privateData_->chIso04noVeto;
  delete privateData_->chIso05noVeto;
  delete privateData_->nhIso03veto;
  delete privateData_->nhIso04veto;
  delete privateData_->nhIso05veto;
  delete privateData_->nhIso03noVeto;
  delete privateData_->nhIso04noVeto;
  delete privateData_->nhIso05noVeto;
  delete privateData_->phIso03veto;
  delete privateData_->phIso04veto;
  delete privateData_->phIso05veto;
  delete privateData_->phIso03noVeto;
  delete privateData_->phIso04noVeto;
  delete privateData_->phIso05noVeto;
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
  
  // counting the number of PF electrons candidates
  int howManyEle = 0;
  
  if(collection) {
    
    for(int index = 0; index < (int)collection->size(); index++) {
      
      const Candidate *cand = &(collection->at(index));
      const PFCandidateRef pflowCandRef = collection->refAt(index).castTo<PFCandidateRef>();

      if ( (!(pflowCandRef.isNull()) ) ) {

	const PFCandidate::ParticleType ptype = pflowCandRef->particleId();	
	if ( ptype==PFCandidate::e ) {  // only for PF electrons
	  howManyEle++; // counting
	  
	  // fill basic kinematics
	  if(saveCand_) writeCandInfo(cand,iEvent,iSetup);	  

	  // basic PF infos
	  if(savePFEleBasic_) writePFEleBasicInfo(pflowCandRef);
	  
	  // tracker based infos
	  GsfTrackRef gsfRef  = pflowCandRef->gsfTrackRef();
	  TrackRef kfTrackRef = pflowCandRef->trackRef();
	  if(savePFEleTrk_) writePFEleTrkInfo(gsfRef,kfTrackRef);

	  // computing isolation: building iso-deposits
	  if (savePFEleIsoDep_) {

	    const reco::PFCandidate & pfCandidate = *pflowCandRef;   
	    Direction pfDir(pfCandidate.eta(), pfCandidate.phi());
	    IsoDeposit photonIsoDep(pfDir);
	    IsoDeposit chIsoDep(pfDir);
	    IsoDeposit nhIsoDep(pfDir);
	    
	    // here I loop over the full PF candidates collection	  
	    for(int jindex = 0; jindex < (int)collection->size(); jindex++) {  
	      
	      const PFCandidateRef c2Ref = collection->refAt(jindex).castTo<PFCandidateRef>();
	      const PFCandidate::ParticleType ptype2 = c2Ref->particleId();
	      
	      // check I'm not checking the same  candidate as my electron
	      // if( index == jindex ) continue; 
	      if ( ptype2==PFCandidate::e ) {  // can be the same only if reconstructed as ele
		if (  fabs(pflowCandRef->energy() - c2Ref->energy()) < 1e-6 ) { // matched energy
		  if (  fabs(pflowCandRef->eta() - c2Ref->eta()) < 1e-6 ) {     // matched eta
		    if (  fabs(pflowCandRef->phi() - c2Ref->phi()) < 1e-6 ) {   // matched phi
		      continue;  // they are actually the same particle! 
		    }}}}
	      
	      // selecting deposits in dR < 1 wrt the candidate
	      Direction dirPfc(c2Ref->eta(), c2Ref->phi());
	      double dR = pfDir.deltaR(dirPfc);
	      if(dR > 1) continue;        // hardcoded!!
	      
	      double pt = c2Ref->pt();
	      // only for charged hadrons, consider only those coming from the same vertex as the electron
	      if( c2Ref->particleId() == 1) {   // charged 
		if ( fabs(c2Ref->vertex().z() - pflowCandRef->vertex().z()) < 0.2 ) {    // hardcoded!!
		  chIsoDep.addDeposit(dirPfc, pt); 
		}
	      }
	      if( c2Ref->particleId() == 5) nhIsoDep.addDeposit(dirPfc, pt);       // neutral
	      if( c2Ref->particleId() == 4) photonIsoDep.addDeposit(dirPfc, pt);   // gamma
	    }
	    
	    // building vetos
	    AbsVetos emptyVetos;
	    AbsVetos photonsVetos, chargedVetos, neutralVetos; 
	    
	    AbsVeto *myConeThreshVetoCH = new ConeThresholdVeto(pfDir, 0.04, 0.7);
	    chargedVetos.push_back( myConeThreshVetoCH );
	    
	    AbsVeto *myConeVetoNH = new ConeVeto(pfDir, 0.07);
	    neutralVetos.push_back( myConeVetoNH );
	    
	    AbsVeto *myRectVetoPH = new RectangularEtaPhiVeto( pfDir, -0.02, 0.02, -0.3, 0.3 );   // could be studied
	    photonsVetos.push_back( myRectVetoPH );
	    
	    double chIso03veto_   = chIsoDep.depositAndCountWithin( 0.3, chargedVetos, false ).first;
	    double chIso04veto_   = chIsoDep.depositAndCountWithin( 0.4, chargedVetos, false ).first;
	    double chIso05veto_   = chIsoDep.depositAndCountWithin( 0.5, chargedVetos, false ).first;
	    double chIso03noVeto_ = chIsoDep.depositAndCountWithin( 0.3, emptyVetos, true ).first;
	    double chIso04noVeto_ = chIsoDep.depositAndCountWithin( 0.4, emptyVetos, true ).first;
	    double chIso05noVeto_ = chIsoDep.depositAndCountWithin( 0.5, emptyVetos, true ).first;
	    
	    double nhIso03veto_   = nhIsoDep.depositAndCountWithin( 0.3, neutralVetos, false ).first;
	    double nhIso04veto_   = nhIsoDep.depositAndCountWithin( 0.4, neutralVetos, false ).first;
	    double nhIso05veto_   = nhIsoDep.depositAndCountWithin( 0.5, neutralVetos, false ).first;
	    double nhIso03noVeto_ = nhIsoDep.depositAndCountWithin( 0.3, emptyVetos, true ).first;
	    double nhIso04noVeto_ = nhIsoDep.depositAndCountWithin( 0.4, emptyVetos, true ).first;
	    double nhIso05noVeto_ = nhIsoDep.depositAndCountWithin( 0.5, emptyVetos, true ).first;
	    
	    double phIso03veto_   = photonIsoDep.depositAndCountWithin( 0.3, photonsVetos, false ).first;
	    double phIso04veto_   = photonIsoDep.depositAndCountWithin( 0.4, photonsVetos, false ).first;
	    double phIso05veto_   = photonIsoDep.depositAndCountWithin( 0.5, photonsVetos, false ).first;
	    double phIso03noVeto_ = photonIsoDep.depositAndCountWithin( 0.3, emptyVetos, true ).first;
	    double phIso04noVeto_ = photonIsoDep.depositAndCountWithin( 0.4, emptyVetos, true ).first;
	    double phIso05noVeto_ = photonIsoDep.depositAndCountWithin( 0.5, emptyVetos, true ).first;
	    
	    privateData_->chIso03veto  ->push_back(chIso03veto_);
	    privateData_->chIso04veto  ->push_back(chIso04veto_);
	    privateData_->chIso05veto  ->push_back(chIso05veto_);
	    privateData_->chIso03noVeto->push_back(chIso03noVeto_);
	    privateData_->chIso04noVeto->push_back(chIso04noVeto_);
	    privateData_->chIso05noVeto->push_back(chIso05noVeto_);
	    privateData_->nhIso03veto  ->push_back(nhIso03veto_);
	    privateData_->nhIso04veto  ->push_back(nhIso04veto_);
	    privateData_->nhIso05veto  ->push_back(nhIso05veto_);
	    privateData_->nhIso03noVeto->push_back(nhIso03noVeto_);
	    privateData_->nhIso04noVeto->push_back(nhIso04noVeto_);
	    privateData_->nhIso05noVeto->push_back(nhIso05noVeto_);
	    privateData_->phIso03veto  ->push_back(phIso03veto_);
	    privateData_->phIso04veto  ->push_back(phIso04veto_);
	    privateData_->phIso05veto  ->push_back(phIso05veto_);
	    privateData_->phIso03noVeto->push_back(phIso03noVeto_);
	    privateData_->phIso04noVeto->push_back(phIso04noVeto_);
	    privateData_->phIso05noVeto->push_back(phIso05noVeto_);

	  } // save iso infos

	} // cand is an electron
      } // null candidate
      else {
	edm::LogWarning("CmsPFlowElectronFiller") << "Warning! The collection seems to be not made by "
						  << "pflow candidates electrons, electron-specific infos will be set to default.";       
      }
    }
    
  }

  if(hitLimitsMeansNoOutput_ && howManyEle > maxTracks_){
    edm::LogInfo("CmsPFlowElectronFiller") << "Track length " << howManyEle
					   << " is too long for declared max length for tree "
					   << maxTracks_ << " and no output flag is set."
					   << " No tracks written to tuple for this event ";
    return;
  }
  
  if(howManyEle > maxTracks_){
    edm::LogInfo("CmsPFlowElectronFiller") << "Track length " << howManyEle  
					   << " is too long for declared max length for tree "
					   << maxTracks_ 
					   << ". Collection will be truncated ";
  }
  
    

  
  // The class member vectors containing the relevant quantities 
  // have all been filled. Now transfer those we want into the 
  // tree 
  int blockSize = (collection) ? howManyEle : 0;
  std::string nCandString = columnPrefix+(*trkIndexName_)+columnSuffix; 
  cmstree->column(nCandString.c_str(),blockSize,0,"Reco");
  
  if(saveCand_)        treeCandInfo(columnPrefix,columnSuffix);
  if(savePFEleBasic_)  treePFEleBasicInfo(columnPrefix,columnSuffix);
  if(savePFEleTrk_)    treePFEleTrkInfo(columnPrefix,columnSuffix);
  if(savePFEleIsoDep_) treePFEleIsoInfo(columnPrefix,columnSuffix);   
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

  if (pflowCandRef.isNonnull()) {    
    privateData_->MvaOutput->push_back(pflowCandRef->mva_e_pi());
    privateData_->PS1Energy->push_back(pflowCandRef->pS1Energy());
    privateData_->PS2Energy->push_back(pflowCandRef->pS2Energy());
    privateData_->EcalEnergy->push_back(pflowCandRef->ecalEnergy());
    privateData_->HcalEnergy->push_back(pflowCandRef->hcalEnergy());    
    privateData_->RawEcalEnergy->push_back(pflowCandRef->rawEcalEnergy());
    privateData_->RawHcalEnergy->push_back(pflowCandRef->rawHcalEnergy());    
    privateData_->PositionAtEcalX->push_back(pflowCandRef->positionAtECALEntrance().x());
    privateData_->PositionAtEcalY->push_back(pflowCandRef->positionAtECALEntrance().y());
    privateData_->PositionAtEcalZ->push_back(pflowCandRef->positionAtECALEntrance().z());
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
  cmstree->column((colPrefix+"HcalEnergy"+colSuffix).c_str(),  *privateData_->HcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"RawEcalEnergy"+colSuffix).c_str(),  *privateData_->RawEcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"RawHcalEnergy"+colSuffix).c_str(),  *privateData_->RawHcalEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PositionAtEcalX"+colSuffix).c_str(),  *privateData_->PositionAtEcalX, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PositionAtEcalY"+colSuffix).c_str(),  *privateData_->PositionAtEcalY, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PositionAtEcalZ"+colSuffix).c_str(),  *privateData_->PositionAtEcalZ, nCandString.c_str(), 0, "Reco");
}

void CmsPFlowElectronFiller::treePFEleIsoInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"chIso03veto"+colSuffix).c_str(),    *privateData_->chIso03veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso04veto"+colSuffix).c_str(),    *privateData_->chIso04veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso05veto"+colSuffix).c_str(),    *privateData_->chIso05veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso03noVeto"+colSuffix).c_str(),  *privateData_->chIso03noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso04noVeto"+colSuffix).c_str(),  *privateData_->chIso04noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso05noVeto"+colSuffix).c_str(),  *privateData_->chIso05noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso03veto"+colSuffix).c_str(),    *privateData_->nhIso03veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso04veto"+colSuffix).c_str(),    *privateData_->nhIso04veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso05veto"+colSuffix).c_str(),    *privateData_->nhIso05veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso03noVeto"+colSuffix).c_str(),  *privateData_->nhIso03noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso04noVeto"+colSuffix).c_str(),  *privateData_->nhIso04noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso05noVeto"+colSuffix).c_str(),  *privateData_->nhIso05noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso03veto"+colSuffix).c_str(),    *privateData_->phIso03veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso04veto"+colSuffix).c_str(),    *privateData_->phIso04veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso05veto"+colSuffix).c_str(),    *privateData_->phIso05veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso03noVeto"+colSuffix).c_str(),  *privateData_->phIso03noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso04noVeto"+colSuffix).c_str(),  *privateData_->phIso04noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso05noVeto"+colSuffix).c_str(),  *privateData_->phIso05noVeto, nCandString.c_str(), 0, "Reco");
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
  HcalEnergy = new vector<float>;
  RawEcalEnergy = new vector<float>;
  RawHcalEnergy = new vector<float>;
  PositionAtEcalX = new vector<float>;
  PositionAtEcalY = new vector<float>;
  PositionAtEcalZ = new vector<float>;
  chIso03veto   = new vector<float>;
  chIso04veto   = new vector<float>;
  chIso05veto   = new vector<float>;
  chIso03noVeto = new vector<float>;
  chIso04noVeto = new vector<float>;
  chIso05noVeto = new vector<float>;
  nhIso03veto   = new vector<float>;
  nhIso04veto   = new vector<float>;
  nhIso05veto   = new vector<float>;
  nhIso03noVeto = new vector<float>;
  nhIso04noVeto = new vector<float>;
  nhIso05noVeto = new vector<float>;
  phIso03veto   = new vector<float>;
  phIso04veto   = new vector<float>;
  phIso05veto   = new vector<float>;
  phIso03noVeto = new vector<float>;
  phIso04noVeto = new vector<float>;
  phIso05noVeto = new vector<float>;
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
  HcalEnergy->clear();
  RawEcalEnergy->clear();
  RawHcalEnergy->clear();
  PositionAtEcalX->clear();
  PositionAtEcalY->clear();
  PositionAtEcalZ->clear();
  chIso03veto   ->clear();
  chIso04veto   ->clear();
  chIso05veto   ->clear();
  chIso03noVeto ->clear();
  chIso04noVeto ->clear();
  chIso05noVeto ->clear();
  nhIso03veto   ->clear();
  nhIso04veto   ->clear();
  nhIso05veto   ->clear();
  nhIso03noVeto ->clear();
  nhIso04noVeto ->clear();
  nhIso05noVeto ->clear();
  phIso03veto   ->clear();
  phIso04veto   ->clear();
  phIso05veto   ->clear();
  phIso03noVeto ->clear();
  phIso04noVeto ->clear();
  phIso05noVeto ->clear();
}
