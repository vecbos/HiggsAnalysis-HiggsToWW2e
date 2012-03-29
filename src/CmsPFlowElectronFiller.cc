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

#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementGsfTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementBrem.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h" 
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsPFlowElectronFiller.h"

#include <TTree.h>
#include <string>
#include "TMath.h"

using namespace edm;
using namespace reco;
using namespace std;
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
  
  savePFEleBasic_  = true;
  savePFEleIsoDep_ = true;
  trkIndexName_    = new std::string("n");
  
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
  delete privateData_->convDist;
  delete privateData_->convDcot;
  delete privateData_->convRadius;
  delete privateData_->convTrackIndex;
  delete privateData_->deltaEtaAtVtxEg;
  delete privateData_->deltaPhiAtVtxEg;
  // delete privateData_->sigmaEtaEtaEg;
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
  delete privateData_->chIso03vetoNVC;
  delete privateData_->chIso04vetoNVC;
  delete privateData_->chIso05vetoNVC;
  delete privateData_->chIso03noVetoNVC;
  delete privateData_->chIso04noVetoNVC;
  delete privateData_->chIso05noVetoNVC;
  delete privateData_->nhIso03veto;
  delete privateData_->nhIso04veto;
  delete privateData_->nhIso05veto;
  delete privateData_->nhIso03vetoTJ;
  delete privateData_->nhIso04vetoTJ;
  delete privateData_->nhIso05vetoTJ;
  delete privateData_->nhIso03noVeto;
  delete privateData_->nhIso04noVeto;
  delete privateData_->nhIso05noVeto;
  delete privateData_->phIso03veto;
  delete privateData_->phIso04veto;
  delete privateData_->phIso05veto;
  delete privateData_->phIso03vetoTJ;
  delete privateData_->phIso04vetoTJ;
  delete privateData_->phIso05vetoTJ;
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
  

  // for conversion rejection  
  try { iEvent.getByLabel(generalTracks_, h_tracks); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsPFElectronFiller") << "Can't get general track collection: " << generalTracks_; }

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

	  // eleID a la egamma infos
	  if(savePFEleBasic_) writePFEleEgammaIDInfo(pflowCandRef);
	  
	  // tracker based infos
	  GsfTrackRef gsfRef  = pflowCandRef->gsfTrackRef();
	  TrackRef kfTrackRef = pflowCandRef->trackRef();
	  if(savePFEleBasic_) writePFEleTrkInfo(pflowCandRef,gsfRef,kfTrackRef,iSetup);  

	  // computing isolation: building iso-deposits
	  if (savePFEleIsoDep_) {

	    const reco::PFCandidate & pfCandidate = *pflowCandRef;   
	    Direction pfDir(pfCandidate.eta(), pfCandidate.phi());
	    IsoDeposit photonIsoDep(pfDir);
	    IsoDeposit chIsoDep(pfDir);
	    IsoDeposit chIsoDepNoVtxComp(pfDir);
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
	      if( c2Ref->particleId() == 1) chIsoDepNoVtxComp.addDeposit(dirPfc, pt);   // charged no vtx compatibility
	      if( c2Ref->particleId() == 5) nhIsoDep.addDeposit(dirPfc, pt);            // neutral
	      if( c2Ref->particleId() == 4) photonIsoDep.addDeposit(dirPfc, pt);        // gamma
	    }
	    
	    // building vetos
	    AbsVetos emptyVetos;
	    AbsVetos photonsVetos, photonsVetosTJ, chargedVetos, neutralVetos, neutralVetosTJ; 
	    
	    // charged hadrons, as from Mihal's studies
	    AbsVeto *myConeThreshVetoCH = new ConeThresholdVeto(pfDir, 0.04, 0.7);
	    chargedVetos.push_back( myConeThreshVetoCH );

	    // neutral hadrons, as from Mihal's studies
	    AbsVeto *myConeVetoNH = new ConeVeto(pfDir, 0.07);
	    neutralVetos.push_back( myConeVetoNH );

	    // neutral hadrons, as from TJ's studies
	    AbsVeto *myConeThreshVetoNH_TJ = new ConeThresholdVeto(pfDir, 0.1, 0.5);  
	    neutralVetosTJ.push_back( myConeThreshVetoNH_TJ );
	    
	    // photons, as from Mihal's studies
	    AbsVeto *myRectVetoPH = new RectangularEtaPhiVeto( pfDir, -0.02, 0.02, -0.3, 0.3 );  
	    photonsVetos.push_back( myRectVetoPH );
	    
	    // photons, as from TJ's studies
	    AbsVeto *myRectVetoPH_TJ   = new RectangularEtaPhiVeto( pfDir, -0.025, 0.025, -0.3, 0.3 );  
	    AbsVeto *myThreshVetoPH_TJ = new ThresholdVeto(0.5);
	    photonsVetosTJ.push_back( myRectVetoPH_TJ );    
	    photonsVetosTJ.push_back( myThreshVetoPH_TJ ); 

	    double chIso03veto_   = chIsoDep.depositAndCountWithin( 0.3, chargedVetos, false ).first;
	    double chIso04veto_   = chIsoDep.depositAndCountWithin( 0.4, chargedVetos, false ).first;
	    double chIso05veto_   = chIsoDep.depositAndCountWithin( 0.5, chargedVetos, false ).first;
	    double chIso03noVeto_ = chIsoDep.depositAndCountWithin( 0.3, emptyVetos, true ).first;
	    double chIso04noVeto_ = chIsoDep.depositAndCountWithin( 0.4, emptyVetos, true ).first;
	    double chIso05noVeto_ = chIsoDep.depositAndCountWithin( 0.5, emptyVetos, true ).first;

	    double chIso03vetoNVC_   = chIsoDepNoVtxComp.depositAndCountWithin( 0.3, chargedVetos, false ).first;
	    double chIso04vetoNVC_   = chIsoDepNoVtxComp.depositAndCountWithin( 0.4, chargedVetos, false ).first;
	    double chIso05vetoNVC_   = chIsoDepNoVtxComp.depositAndCountWithin( 0.5, chargedVetos, false ).first;
	    double chIso03noVetoNVC_ = chIsoDepNoVtxComp.depositAndCountWithin( 0.3, emptyVetos, true ).first;
	    double chIso04noVetoNVC_ = chIsoDepNoVtxComp.depositAndCountWithin( 0.4, emptyVetos, true ).first;
	    double chIso05noVetoNVC_ = chIsoDepNoVtxComp.depositAndCountWithin( 0.5, emptyVetos, true ).first;
	    
	    double nhIso03veto_   = nhIsoDep.depositAndCountWithin( 0.3, neutralVetos, false ).first;
	    double nhIso04veto_   = nhIsoDep.depositAndCountWithin( 0.4, neutralVetos, false ).first;
	    double nhIso05veto_   = nhIsoDep.depositAndCountWithin( 0.5, neutralVetos, false ).first;
	    double nhIso03vetoTJ_ = nhIsoDep.depositAndCountWithin( 0.3, neutralVetosTJ, false ).first;
	    double nhIso04vetoTJ_ = nhIsoDep.depositAndCountWithin( 0.4, neutralVetosTJ, false ).first;
	    double nhIso05vetoTJ_ = nhIsoDep.depositAndCountWithin( 0.5, neutralVetosTJ, false ).first;
	    double nhIso03noVeto_ = nhIsoDep.depositAndCountWithin( 0.3, emptyVetos, true ).first;
	    double nhIso04noVeto_ = nhIsoDep.depositAndCountWithin( 0.4, emptyVetos, true ).first;
	    double nhIso05noVeto_ = nhIsoDep.depositAndCountWithin( 0.5, emptyVetos, true ).first;
	    
	    double phIso03veto_   = photonIsoDep.depositAndCountWithin( 0.3, photonsVetos, false ).first;
	    double phIso04veto_   = photonIsoDep.depositAndCountWithin( 0.4, photonsVetos, false ).first;
	    double phIso05veto_   = photonIsoDep.depositAndCountWithin( 0.5, photonsVetos, false ).first;
	    double phIso03vetoTJ_ = photonIsoDep.depositAndCountWithin( 0.3, photonsVetosTJ, false ).first;
	    double phIso04vetoTJ_ = photonIsoDep.depositAndCountWithin( 0.4, photonsVetosTJ, false ).first;
	    double phIso05vetoTJ_ = photonIsoDep.depositAndCountWithin( 0.5, photonsVetosTJ, false ).first;
	    double phIso03noVeto_ = photonIsoDep.depositAndCountWithin( 0.3, emptyVetos, true ).first;
	    double phIso04noVeto_ = photonIsoDep.depositAndCountWithin( 0.4, emptyVetos, true ).first;
	    double phIso05noVeto_ = photonIsoDep.depositAndCountWithin( 0.5, emptyVetos, true ).first;
	    
	    privateData_->chIso03veto  ->push_back(chIso03veto_);
	    privateData_->chIso04veto  ->push_back(chIso04veto_);
	    privateData_->chIso05veto  ->push_back(chIso05veto_);
	    privateData_->chIso03noVeto->push_back(chIso03noVeto_);
	    privateData_->chIso04noVeto->push_back(chIso04noVeto_);
	    privateData_->chIso05noVeto->push_back(chIso05noVeto_);
	    privateData_->chIso03vetoNVC  ->push_back(chIso03vetoNVC_);
	    privateData_->chIso04vetoNVC  ->push_back(chIso04vetoNVC_);
	    privateData_->chIso05vetoNVC  ->push_back(chIso05vetoNVC_);
	    privateData_->chIso03noVetoNVC->push_back(chIso03noVetoNVC_);
	    privateData_->chIso04noVetoNVC->push_back(chIso04noVetoNVC_);
	    privateData_->chIso05noVetoNVC->push_back(chIso05noVetoNVC_);
	    privateData_->nhIso03veto  ->push_back(nhIso03veto_);
	    privateData_->nhIso04veto  ->push_back(nhIso04veto_);
	    privateData_->nhIso05veto  ->push_back(nhIso05veto_);
	    privateData_->nhIso03vetoTJ  ->push_back(nhIso03vetoTJ_);
	    privateData_->nhIso04vetoTJ  ->push_back(nhIso04vetoTJ_);
	    privateData_->nhIso05vetoTJ  ->push_back(nhIso05vetoTJ_);
	    privateData_->nhIso03noVeto->push_back(nhIso03noVeto_);
	    privateData_->nhIso04noVeto->push_back(nhIso04noVeto_);
	    privateData_->nhIso05noVeto->push_back(nhIso05noVeto_);
	    privateData_->phIso03veto  ->push_back(phIso03veto_);
	    privateData_->phIso04veto  ->push_back(phIso04veto_);
	    privateData_->phIso05veto  ->push_back(phIso05veto_);
	    privateData_->phIso03vetoTJ->push_back(phIso03vetoTJ_);
	    privateData_->phIso04vetoTJ->push_back(phIso04vetoTJ_);
	    privateData_->phIso05vetoTJ->push_back(phIso05vetoTJ_);
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
  if(savePFEleBasic_)  { 
    treePFEleBasicInfo(columnPrefix,columnSuffix);
    treePFEleTrkInfo(columnPrefix,columnSuffix);
    treePFEleEgammaIDInfo(columnPrefix,columnSuffix);
  }
  if(savePFEleIsoDep_) treePFEleIsoInfo(columnPrefix,columnSuffix);   
  if(dumpData) cmstree->dumpData();

  delete trkIndexName_;

}

void CmsPFlowElectronFiller::writePFEleTrkInfo(const PFCandidateRef pflowCandRef, reco::GsfTrackRef gsfRef, reco::TrackRef kfTrackRef, const edm::EventSetup& iSetup ) {
  
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
  
  // for conversion rejection
  if( h_tracks.isValid() ) {
    
    edm::ESHandle<MagneticField> theMagField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
    const MagneticField *magField=&(*theMagField);
    
    // if track not found these values will not be changed
    theCFDist_       = -999.;
    theCFDcot_       = -999.;
    theCFRconv_      = -999.;
    theCFTrackIndex_ = -999;
    
    GlobalPoint origin(0.,0.,0.);
    getConversionInfo(pflowCandRef,h_tracks,(magField->inTesla(origin)).mag());
    
    privateData_->convDist->push_back(theCFDist_);
    privateData_->convDcot->push_back(theCFDcot_);
    privateData_->convRadius->push_back(theCFRconv_);
    privateData_->convTrackIndex->push_back(theCFTrackIndex_);
  } else {
    privateData_->convDist->push_back(999);
    privateData_->convDcot->push_back(999);
    privateData_->convRadius->push_back(999);
    privateData_->convTrackIndex->push_back(-1);
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

void CmsPFlowElectronFiller::writePFEleEgammaIDInfo(const reco::PFCandidateRef pflowCandRef) {

  const reco::PFCandidate & pfCandidate = *pflowCandRef;   

  if (pflowCandRef.isNonnull()) {      
    
    const PFCandidate::ElementsInBlocks& theElements = pfCandidate.elementsInBlocks();

    typedef PFCandidate::ElementsInBlocks::const_iterator IEB; 
    
    int firstEcal = 0;
    unsigned int gsf_index = 10000;
    // float sigmaEtaEta  = -999.;
    float deta_gsfecal = -999.;
    float dphi_gsfecal = -999.;
    float clusterEta = -999.;
    float clusterPhi = -999.;
    float trackEta   = -999.;
    float trackPhi   = -999.;
    
    for (IEB ieb=theElements.begin(); ieb<theElements.end(); ++ieb) {
      const PFBlock& block = *(ieb->first);
      const PFBlockElement& pfbe = block.elements()[ieb->second];
      
      if(pfbe.type()==reco::PFBlockElement::ECAL) {
	if(firstEcal == 0) {
	  reco::PFClusterRef clusterRef = pfbe.clusterRef();
	  // vector< const reco::PFCluster * > pfClust_vec(0);
	  // pfClust_vec.clear();
	  // pfClust_vec.push_back(&(*clusterRef));
	  
	  // PFClusterWidthAlgo pfwidth(pfClust_vec);
	  // sigmaEtaEta = pfwidth.pflowSigmaEtaEta();	
	  clusterEta = clusterRef->position().eta();  
	  clusterPhi = clusterRef->position().phi();  
	}
	firstEcal++;
      }
      
      if(pfbe.type()==reco::PFBlockElement::GSF) { 
	gsf_index = ieb->second;
	
	const reco::PFBlockElementGsfTrack * GsfEl  =  
	  dynamic_cast<const reco::PFBlockElementGsfTrack*>(&pfbe);
	trackEta = GsfEl->positionAtECALEntrance().eta();   
	trackPhi = GsfEl->positionAtECALEntrance().phi();
      }
    }

    deta_gsfecal = clusterEta - trackEta;
    dphi_gsfecal = clusterPhi - trackPhi;
    if ( dphi_gsfecal < -3.14159265358979323846 ) 
      dphi_gsfecal = dphi_gsfecal + 2.*3.14159265358979323846;
    else if ( dphi_gsfecal > 3.14159265358979323846 ) 
      dphi_gsfecal = dphi_gsfecal - 2.*3.14159265358979323846;
    
    privateData_->deltaEtaAtVtxEg->push_back( deta_gsfecal );
    privateData_->deltaPhiAtVtxEg->push_back( dphi_gsfecal );
    // privateData_->sigmaEtaEtaEg->push_back( sigmaEtaEta );
  }
  else {
    privateData_->deltaEtaAtVtxEg->push_back( -1. );
    privateData_->deltaPhiAtVtxEg->push_back( -1. );
    // privateData_->sigmaEtaEtaEg->push_back( -1. );
  }
}


void CmsPFlowElectronFiller::treePFEleTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"gsfTrackIndex"+colSuffix).c_str(), *privateData_->gsfTrackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(), *privateData_->trackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"convDist"+colSuffix).c_str(), *privateData_->convDist, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"convDcot"+colSuffix).c_str(), *privateData_->convDcot, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"convRadius"+colSuffix).c_str(), *privateData_->convRadius, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"convTrackIndex"+colSuffix).c_str(), *privateData_->convTrackIndex, nCandString.c_str(), 0, "Reco");
}

void CmsPFlowElectronFiller::treePFEleEgammaIDInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"deltaEtaAtVtxEg"+colSuffix).c_str(),  *privateData_->deltaEtaAtVtxEg, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"deltaPhiAtVtxEg"+colSuffix).c_str(),  *privateData_->deltaPhiAtVtxEg, nCandString.c_str(), 0, "Reco");
  // cmstree->column((colPrefix+"sigmaEtaEtaEg"+colSuffix).c_str(),  *privateData_->sigmaEtaEtaEg, nCandString.c_str(), 0, "Reco");
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
  cmstree->column((colPrefix+"chIso03vetoNVC"+colSuffix).c_str(),    *privateData_->chIso03vetoNVC,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso04vetoNVC"+colSuffix).c_str(),    *privateData_->chIso04vetoNVC,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso05vetoNVC"+colSuffix).c_str(),    *privateData_->chIso05vetoNVC,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso03noVetoNVC"+colSuffix).c_str(),  *privateData_->chIso03noVetoNVC, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso04noVetoNVC"+colSuffix).c_str(),  *privateData_->chIso04noVetoNVC, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"chIso05noVetoNVC"+colSuffix).c_str(),  *privateData_->chIso05noVetoNVC, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso03veto"+colSuffix).c_str(),    *privateData_->nhIso03veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso04veto"+colSuffix).c_str(),    *privateData_->nhIso04veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso05veto"+colSuffix).c_str(),    *privateData_->nhIso05veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso03vetoTJ"+colSuffix).c_str(),  *privateData_->nhIso03vetoTJ, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso04vetoTJ"+colSuffix).c_str(),  *privateData_->nhIso04vetoTJ, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso05vetoTJ"+colSuffix).c_str(),  *privateData_->nhIso05vetoTJ, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso03noVeto"+colSuffix).c_str(),  *privateData_->nhIso03noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso04noVeto"+colSuffix).c_str(),  *privateData_->nhIso04noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nhIso05noVeto"+colSuffix).c_str(),  *privateData_->nhIso05noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso03veto"+colSuffix).c_str(),    *privateData_->phIso03veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso04veto"+colSuffix).c_str(),    *privateData_->phIso04veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso05veto"+colSuffix).c_str(),    *privateData_->phIso05veto,   nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso03vetoTJ"+colSuffix).c_str(),  *privateData_->phIso03vetoTJ, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso04vetoTJ"+colSuffix).c_str(),  *privateData_->phIso04vetoTJ, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso05vetoTJ"+colSuffix).c_str(),  *privateData_->phIso05vetoTJ, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso03noVeto"+colSuffix).c_str(),  *privateData_->phIso03noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso04noVeto"+colSuffix).c_str(),  *privateData_->phIso04noVeto, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"phIso05noVeto"+colSuffix).c_str(),  *privateData_->phIso05noVeto, nCandString.c_str(), 0, "Reco");
}


void CmsPFlowElectronFillerData::initialise() {  
  initialiseCandidate();

  // track infos
  trackIndex = new vector<int>;
  gsfTrackIndex = new vector<int>;  
  convDist = new vector<float>;
  convDcot = new vector<float>;
  convRadius = new vector<float>;
  convTrackIndex = new vector<int>;

  // id a la egamma
  deltaEtaAtVtxEg = new vector<float>;
  deltaPhiAtVtxEg = new vector<float>;
  // sigmaEtaEtaEg   = new vector<float>;

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

  // isolation
  chIso03veto   = new vector<float>;
  chIso04veto   = new vector<float>;
  chIso05veto   = new vector<float>;
  chIso03noVeto = new vector<float>;
  chIso04noVeto = new vector<float>;
  chIso05noVeto = new vector<float>;
  chIso03vetoNVC   = new vector<float>;
  chIso04vetoNVC   = new vector<float>;
  chIso05vetoNVC   = new vector<float>;
  chIso03noVetoNVC = new vector<float>;
  chIso04noVetoNVC = new vector<float>;
  chIso05noVetoNVC = new vector<float>;
  nhIso03veto   = new vector<float>;
  nhIso04veto   = new vector<float>;
  nhIso05veto   = new vector<float>;
  nhIso03vetoTJ = new vector<float>;
  nhIso04vetoTJ = new vector<float>;
  nhIso05vetoTJ = new vector<float>;
  nhIso03noVeto = new vector<float>;
  nhIso04noVeto = new vector<float>;
  nhIso05noVeto = new vector<float>;
  phIso03veto   = new vector<float>;
  phIso04veto   = new vector<float>;
  phIso05veto   = new vector<float>;
  phIso03vetoTJ = new vector<float>;
  phIso04vetoTJ = new vector<float>;
  phIso05vetoTJ = new vector<float>;
  phIso03noVeto = new vector<float>;
  phIso04noVeto = new vector<float>;
  phIso05noVeto = new vector<float>;
}

void CmsPFlowElectronFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  trackIndex     -> clear();
  gsfTrackIndex  -> clear();
  convDist       -> clear(); 
  convDcot       -> clear(); 
  convRadius     -> clear();
  convTrackIndex -> clear(); 

  deltaEtaAtVtxEg -> clear();
  deltaPhiAtVtxEg -> clear();
  // sigmaEtaEtaEg   -> clear();

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
  chIso03vetoNVC   ->clear();
  chIso04vetoNVC   ->clear();
  chIso05vetoNVC   ->clear();
  chIso03noVetoNVC ->clear();
  chIso04noVetoNVC ->clear();
  chIso05noVetoNVC ->clear();
  nhIso03veto   ->clear();
  nhIso04veto   ->clear();
  nhIso05veto   ->clear();
  nhIso03vetoTJ ->clear();
  nhIso04vetoTJ ->clear();
  nhIso05vetoTJ ->clear();
  nhIso03noVeto ->clear();
  nhIso04noVeto ->clear();
  nhIso05noVeto ->clear();
  phIso03veto   ->clear();
  phIso04veto   ->clear();
  phIso05veto   ->clear();
  phIso03vetoTJ ->clear();
  phIso04vetoTJ ->clear();
  phIso05vetoTJ ->clear();
  phIso03noVeto ->clear();
  phIso04noVeto ->clear();
  phIso05noVeto ->clear();
}

void CmsPFlowElectronFiller::getConversionInfo(reco::PFCandidateRef pflowCandRef,
						   const edm::Handle<reco::TrackCollection>& track_h, 
						   const double bFieldAtOrigin) {
  
  const TrackCollection *ctftracks = track_h.product();
  
  if (pflowCandRef.isNonnull()) {     // if null, keep the default value as outside  

    if ((pflowCandRef->trackRef()).isNonnull()) {  // if null, keep the default value as outside  

      int ctfidx = -999.;  

      TrackRef el_ctftrack = pflowCandRef->trackRef();
      ctfidx = static_cast<int>(el_ctftrack.key());
    
      math::XYZTLorentzVector el_tk_p4(pflowCandRef->trackRef()->px(), pflowCandRef->trackRef()->py(), pflowCandRef->trackRef()->pz(), pflowCandRef->trackRef()->p());
    
      int tk_i = 0;
      double mindcot = 9999.;
      //make a null Track Ref
      TrackRef candCtfTrackRef = TrackRef() ;
    
      for(TrackCollection::const_iterator tk = ctftracks->begin();
	  tk != ctftracks->end(); tk++, tk_i++) {
	//if the general Track is the same one as made by the electron, skip it
	if((tk_i == ctfidx))
	  continue;
      
	math::XYZTLorentzVector tk_p4 = math::XYZTLorentzVector(tk->px(), tk->py(),tk->pz(), tk->p());
	
	// look only in a cone of 0.5
	double dR = deltaR(el_tk_p4, tk_p4);
	if(dR>0.5) continue;
      
	//require opp. sign -> Should we use the majority logic??
	if(tk->charge() + pflowCandRef->trackRef()->charge() != 0) continue;
	
	double dcot = fabs(1./tan(tk_p4.theta()) - 1./tan(el_tk_p4.theta()));
	if(dcot < mindcot) {
	  mindcot = dcot;
	  candCtfTrackRef = reco::TrackRef(track_h, tk_i);
	}
      } //track loop
      
      if(candCtfTrackRef.isNonnull()) { 
	
	// now calculate the conversion related information
	double elCurvature = -0.3*bFieldAtOrigin*(pflowCandRef->trackRef()->charge()/el_tk_p4.pt())/100.;
	double rEl = fabs(1./elCurvature);
	double xEl = -1*(1./elCurvature - pflowCandRef->trackRef()->d0())*sin(el_tk_p4.phi());
	double yEl = (1./elCurvature - pflowCandRef->trackRef()->d0())*cos(el_tk_p4.phi());
	
	math::XYZTLorentzVector cand_p4 = math::XYZTLorentzVector(candCtfTrackRef->px(), candCtfTrackRef->py(),candCtfTrackRef->pz(), candCtfTrackRef->p());
	double candCurvature = -0.3*bFieldAtOrigin*(candCtfTrackRef->charge()/cand_p4.pt())/100.;
	double rCand = fabs(1./candCurvature);
	double xCand = -1*(1./candCurvature - candCtfTrackRef->d0())*sin(cand_p4.phi());
	double yCand = (1./candCurvature - candCtfTrackRef->d0())*cos(cand_p4.phi());
	
	double d = sqrt(pow(xEl-xCand, 2) + pow(yEl-yCand , 2));
	double dist = d - (rEl + rCand);
	double dcot = 1./tan(el_tk_p4.theta()) - 1./tan(cand_p4.theta());
      
	//get the point of conversion
	double xa1 = xEl   + (xCand-xEl) * rEl/d;
	double xa2 = xCand + (xEl-xCand) * rCand/d;
	double ya1 = yEl   + (yCand-yEl) * rEl/d;
	double ya2 = yCand + (yEl-yCand) * rCand/d;
	
	double x=.5*(xa1+xa2);
	double y=.5*(ya1+ya2);
	double rconv = sqrt(pow(x,2) + pow(y,2));
	double z = pflowCandRef->trackRef()->dz() + rEl*pflowCandRef->trackRef()->pz()*TMath::ACos(1-pow(rconv,2)/(2.*pow(rEl,2)))/pflowCandRef->trackRef()->pt();
	
	math::XYZPoint convPoint(x, y, z);
	
	// now assign a sign to the radius of conversion
	float tempsign = pflowCandRef->trackRef()->px()*x + pflowCandRef->trackRef()->py()*y;
	tempsign = tempsign/fabs(tempsign);
	rconv = tempsign*rconv;
	
	// saving to bring out of here
	theCFDist_       = dist;
	theCFDcot_       = dcot;
	theCFRconv_      = rconv;
	theCFTrackIndex_ = candCtfTrackRef.key();
      }
    }
  }
}
