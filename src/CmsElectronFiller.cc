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

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"

#include "RecoParticleFlow/PFProducer/interface/PFGeometry.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsElectronFiller.h"

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


CmsElectronFiller::CmsElectronFiller(CmsTree *cmsTree, 
				     int maxTracks, int maxMCTracks,
				     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsElectronFillerData)
{
  cmstree=cmsTree;

  saveTrk_=true;
  saveEcal_=true;
  saveFatTrk_=true;
  saveFatEcal_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  string etamap("RecoParticleFlow/PFBlockProducer/data/resmap_ECAL_eta.dat"); 
  string phimap("RecoParticleFlow/PFBlockProducer/data/resmap_ECAL_phi.dat"); 
  FileInPath ecalEtaMap(etamap);
  FileInPath ecalPhiMap(phimap);
  resMapEtaECAL_ = new PFResolutionMap("ECAL_eta",ecalEtaMap.fullPath().c_str());
  resMapPhiECAL_ = new PFResolutionMap("ECAL_phi",ecalPhiMap.fullPath().c_str());

  privateData_->initialise();
}

CmsElectronFiller::CmsElectronFiller(CmsTree *cmsTree, bool fatTree, 
				     int maxTracks, int maxMCTracks,
				     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsElectronFillerData)
{
  cmstree=cmsTree;

  saveTrk_=true;
  saveEcal_=true;
  saveFatTrk_=fatTree;
  saveFatEcal_=fatTree;
  savePFextra_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  string etamap("RecoParticleFlow/PFBlockProducer/data/resmap_ECAL_eta.dat"); 
  string phimap("RecoParticleFlow/PFBlockProducer/data/resmap_ECAL_phi.dat"); 
  FileInPath ecalEtaMap(etamap);
  FileInPath ecalPhiMap(phimap);
  resMapEtaECAL_ = new PFResolutionMap("ECAL_eta",ecalEtaMap.fullPath().c_str());
  resMapPhiECAL_ = new PFResolutionMap("ECAL_phi",ecalPhiMap.fullPath().c_str());

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsElectronFiller::~CmsElectronFiller() {

  // delete here the vector ptr's
  delete privateData_->fiducialFlags;
  delete privateData_->recoFlags;
  delete privateData_->esEnergy;
  delete privateData_->energyCorrections;

  delete privateData_->superClusterIndex;
  delete privateData_->PFsuperClusterIndex;
  delete privateData_->trackIndex;
  delete privateData_->gsfTrackIndex;

  delete privateData_->PFChi2EcalTrack; 
  delete privateData_->PFesEneL1;
  delete privateData_->PFesEneL2;
  delete privateData_->PFesChi2L1;
  delete privateData_->PFesChi2L2;

  delete privateData_->ncand;

  delete resMapEtaECAL_;
  delete resMapPhiECAL_;
}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsElectronFiller::saveTrk(bool what) { saveTrk_=what;}

void CmsElectronFiller::saveEcal(bool what) { saveEcal_=what;}

void CmsElectronFiller::saveFatTrk(bool what) { saveFatTrk_=what;}

void CmsElectronFiller::saveFatEcal(bool what) { saveFatEcal_=what;}

void CmsElectronFiller::saveEleID(bool what) { saveEleID_=what;}

void CmsElectronFiller::savePFextra(bool what) { savePFextra_=what;}

void CmsElectronFiller::writeCollectionToTree(edm::InputTag collectionTag,
					      const edm::Event& iEvent, const edm::EventSetup& iSetup,
					      const std::string &columnPrefix, const std::string &columnSuffix,
					      bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get electron candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsElectronFiller") << "Track length " << collection->size() 
				   << " is too long for declared max length for tree "
				   << maxTracks_ << " and no output flag is set."
				   << " No tracks written to tuple for this event ";
      return;
    }
  
    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsElectronFiller") << "Track length " << collection->size() 
				   << " is too long for declared max length for tree "
				   << maxTracks_ 
				   << ". Collection will be truncated ";
    }

    *(privateData_->ncand) = collection->size();

    // Magnetic Field
    ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    B_=magneticField->inTesla(GlobalPoint(0,0,0));

    // superclusters
    Handle<SuperClusterCollection> EcalBarrelSuperClusters;
    try { iEvent.getByLabel(EcalBarrelSuperClusters_, EcalBarrelSuperClusters); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL barrel supercluster Collection" << EcalBarrelSuperClusters_; }
    
    Handle<SuperClusterCollection> EcalEndcapSuperClusters;
    try { iEvent.getByLabel(EcalEndcapSuperClusters_, EcalEndcapSuperClusters); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL endcap supercluster Collection" << EcalEndcapSuperClusters_; }
    
    barrelSuperClustersSize = EcalBarrelSuperClusters->size();

    // for cluster shape variables
    Handle< EcalRecHitCollection > EcalBarrelRecHits;
    try { iEvent.getByLabel(EcalBarrelRecHits_, EcalBarrelRecHits); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL barrel rec hits Collection" << EcalBarrelRecHits_; }
    const EcalRecHitCollection *EBRecHits = EcalBarrelRecHits.product();

    Handle< EcalRecHitCollection > EcalEndcapRecHits;
    try { iEvent.getByLabel(EcalEndcapRecHits_, EcalEndcapRecHits); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get ECAL endcap rec hits Collection" << EcalEndcapRecHits_; }
    const EcalRecHitCollection *EERecHits = EcalEndcapRecHits.product();

    for(int index = 0; index < (int)collection->size(); index++) {
      
      // fill basic kinematics
      const Candidate *cand = &(collection->at(index));
      if(saveCand_) writeCandInfo(cand,iEvent,iSetup);

      const GsfElectronRef electronRef = collection->refAt(index).castTo<GsfElectronRef>();

      if ( !(electronRef.isNull()) ) {

        // fill Cluster Adapter
        SuperClusterRef sclusRef = electronRef->superCluster();
        SuperClusterRef pfclusRef = electronRef->pflowSuperCluster();
        if(saveEcal_) writeEcalInfo(electronRef,iEvent,iSetup,sclusRef,pfclusRef,EBRecHits,EERecHits );
        // fill (GSF) Track Adapter
        GsfTrackRef trkRef = cand->get<GsfTrackRef>();
        if(saveTrk_) writeTrkInfo(electronRef,iEvent,iSetup,trkRef);
	// fill extra info for PF seeding
	if(savePFextra_) writePFextraInfo(electronRef,pfclusRef);

      } else {
        edm::LogWarning("CmsElectronFiller") << "Warning! The collection seems to be not made by "
                                             << "electrons, electron-specific infos will be set to default.";
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
  
  if(saveCand_) treeCandInfo(columnPrefix,columnSuffix);
  if(saveEcal_) treeEcalInfo(columnPrefix,columnSuffix);
  if(saveTrk_)  treeTrkInfo(columnPrefix,columnSuffix);
  if(savePFextra_) treePFextraInfo(columnPrefix,columnSuffix); 
  if(saveEleID_) {
    CmsEleIDTreeFiller eIDFiller(cmstree);
    eIDFiller.setStandalone(false);
    eIDFiller.setEcalBarrelRecHits(EcalBarrelRecHits_);
    eIDFiller.setEcalEndcapRecHits(EcalEndcapRecHits_);
    // those for egamma official isolations
    eIDFiller.setTkIsolationProducer(tkIsolationProducer_);
    eIDFiller.setTowerIsolationProducer(towerIsolationProducer_);
    // those for private H->WW isolations
    eIDFiller.setTracksProducer(tracksProducer_);
    eIDFiller.setCalotowersProducer(calotowersProducer_);
    eIDFiller.writeCollectionToTree(collectionTag,iEvent,iSetup,columnPrefix,columnSuffix,false);
  }
  
  if(dumpData) cmstree->dumpData();

}




void CmsElectronFiller::writeTrkInfo(const GsfElectronRef electronRef, 
				     const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				     GsfTrackRef trkRef) {
  if( trkRef.isNonnull() ) {

    privateData_->gsfTrackIndex->push_back(trkRef.key());

    reco::TrackRef closeCtfTrack = electronRef->closestCtfTrackRef();
    if ( closeCtfTrack.isNonnull() ) {
      privateData_->trackIndex->push_back(closeCtfTrack.key());
    } else {
      privateData_->trackIndex->push_back( -1 );
    }

  } else {
    privateData_->gsfTrackIndex->push_back( -1 );
    privateData_->trackIndex->push_back( -1 );
  }
    
}

void CmsElectronFiller::writePFextraInfo(const GsfElectronRef electronRef, SuperClusterRef pfclusRef) {

  // variables to be put in the tree
  ps2En =-1.;
  ps1En =-1.;
  ps2chi=-1.; 
  ps1chi=-1.;
  float chichi = -1.;

  if( electronRef.isNonnull() ) {

    bool isTrackerDriven = electronRef->isTrackerDriven();
    if (isTrackerDriven) {
      
      reco::TrackRef closeCtfTrack = electronRef->closestCtfTrackRef();
      
      if (closeCtfTrack.isNonnull() && pfclusRef.isNonnull()) {
	
	reco::TrackBase::TrackQuality trackQuality_;
	trackQuality_=TrackBase::qualityByName("highPurity");
	if (closeCtfTrack->quality(trackQuality_)){ 

	  // clusters/track matching

	  // fixme! we do not have trajectories in reco. We use tracks instead 
	  // float PTOB=Tj[i].lastMeasurement().updatedState().globalMomentum().mag();
	  // int nhitpi=Tj[i].foundHits();
	  float PTOB=closeCtfTrack->outerMomentum().R();

	  float pfmass = 0.0005;
	  float pfoutenergy=sqrt((pfmass*pfmass)+closeCtfTrack->outerMomentum().Mag2());
	  math::XYZTLorentzVector mom = math::XYZTLorentzVector(closeCtfTrack->outerMomentum().x(),
								closeCtfTrack->outerMomentum().y(),
								closeCtfTrack->outerMomentum().z(),
								pfoutenergy);
	  math::XYZTLorentzVector pos = math::XYZTLorentzVector(closeCtfTrack->outerPosition().x(),
								closeCtfTrack->outerPosition().y(),
								closeCtfTrack->outerPosition().z(),
								0.);
	  BaseParticlePropagator theOutParticle = BaseParticlePropagator( RawParticle(mom,pos),0,0,B_.z());
	  theOutParticle.setCharge(closeCtfTrack->charge());
	  theOutParticle.propagateToEcalEntrance(false);
	  
	  float toteta=1000;
	  float totphi=1000;
	  float EE=0;
	  float feta=0;
	  math::XYZPointF ElecTrkEcalPos(0,0,0);
	  if(theOutParticle.getSuccess()!=0){
	    ElecTrkEcalPos=math::XYZPointF(theOutParticle.vertex().x(),
					   theOutParticle.vertex().y(),
					   theOutParticle.vertex().z());
	    bool isBelowPS=(fabs(theOutParticle.vertex().eta())>1.65) ? true :false;

	    double ecalShowerDepth = PFCluster::getDepthCorrection(pfclusRef->energy(), isBelowPS, false);

	    math::XYZPoint meanShower=math::XYZPoint(theOutParticle.vertex())+
	      math::XYZTLorentzVector(theOutParticle.momentum()).Vect().Unit()*ecalShowerDepth;
	      
	    float etarec=meanShower.eta();
	    float phirec=meanShower.phi();
	    float tmp_ep=pfclusRef->energy()/PTOB;
	    totphi=fabs(pfclusRef->position().phi()-phirec);
	    if (totphi>TMath::TwoPi()) totphi-= TMath::TwoPi();
	    
	    if ((tmp_ep>0.3)&&(tmp_ep<3.)){
	      toteta=pfclusRef->position().eta()-etarec;
	      EE=pfclusRef->energy();
	      feta= pfclusRef->position().eta();
	    }
	  }
	  
	  // resolution maps  
	  double ecaletares = resMapEtaECAL_->GetBinContent(resMapEtaECAL_->FindBin(feta,EE));
	  double ecalphires = resMapPhiECAL_->GetBinContent(resMapPhiECAL_->FindBin(feta,EE)); 
	    
	  // geometrical compatibility
	  float chieta=(toteta!=1000)? toteta/ecaletares : toteta;
	  float chiphi=(totphi!=1000)? totphi/ecalphires : totphi;
	  chichi= sqrt(chieta*chieta + chiphi*chiphi);

	  // preshower in endcap
	  if ( (fabs(closeCtfTrack->eta())>1.68) ) PSforTMVA(mom,pos,pfclusRef);

	}  // track quality 
      } // ok track reference 
    } // tracker driven
  } // ok electron 

  privateData_->PFChi2EcalTrack->push_back(chichi);
  privateData_->PFesEneL1->push_back(ps1En);
  privateData_->PFesEneL2->push_back(ps2En);
  privateData_->PFesChi2L1->push_back(ps1chi);
  privateData_->PFesChi2L2->push_back(ps2chi);
}

void CmsElectronFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(),  *privateData_->trackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsfTrackIndex"+colSuffix).c_str(),  *privateData_->gsfTrackIndex, nCandString.c_str(), 0, "Reco");
}

void CmsElectronFiller::treePFextraInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"PFChi2EcalTrack"+colSuffix).c_str(),  *privateData_->PFChi2EcalTrack, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFesEneL1"+colSuffix).c_str(),   *privateData_->PFesEneL1, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFesEneL2"+colSuffix).c_str(),   *privateData_->PFesEneL2, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFesChi2L1"+colSuffix).c_str(),  *privateData_->PFesChi2L1, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFesChi2L2"+colSuffix).c_str(),  *privateData_->PFesChi2L2, nCandString.c_str(), 0, "Reco");
}

void CmsElectronFiller::writeEcalInfo(const GsfElectronRef electronRef, 
				      const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				      SuperClusterRef sclusRef, SuperClusterRef pfclusRef,
				      const EcalRecHitCollection *EBRecHits,
				      const EcalRecHitCollection *EERecHits) {

  if(&electronRef) {

    // fiducial flags in ECAL
    int packed_sel = -1;
    int isEB = ( electronRef->isEB() ) ? 1 : 0;
    int isEE = ( electronRef->isEE() ) ? 1 : 0;
    int isGap = ( electronRef->isGap() ) ? 1 : 0;
    int isEBEEGap = ( electronRef->isEBEEGap() ) ? 1 : 0;
    int isEBGap = ( electronRef->isEBGap() ) ? 1 : 0;
    int isEBEtaGap = ( electronRef->isEBEtaGap() ) ? 1 : 0;
    int isEBPhiGap = ( electronRef->isEBPhiGap() ) ? 1 : 0;
    int isEEGap = ( electronRef->isEEGap() ) ? 1 : 0;
    int isEEDeeGap = ( electronRef->isEEDeeGap() ) ? 1 : 0;
    int isEERingGap = ( electronRef->isEERingGap() ) ? 1 : 0;
    
    packed_sel = ( isEB << 9 ) | ( isEE << 8 ) | ( isGap << 7 ) |
      ( isEBEEGap << 6 ) | ( isEBGap << 5 ) | ( isEBEtaGap << 4 ) | ( isEBPhiGap << 3 ) |
      ( isEEGap << 2 ) | ( isEEDeeGap << 1 ) | isEERingGap;

    privateData_->fiducialFlags->push_back(packed_sel);

    int packed_reco;
    int isEcalDriven = ( electronRef->isEcalDriven() ) ? 1 : 0;
    int isTrackerDriven = ( electronRef->isTrackerDriven() ) ? 1 : 0;
    packed_reco = ( isEcalDriven << 1 ) | isTrackerDriven;
    privateData_->recoFlags->push_back( packed_reco );

    // link to the supercluster (collections are merged: barrel + endcap in this order)
    if ( isEcalDriven && sclusRef.isNonnull() ) {
      int offset = ( fabs(sclusRef->eta() ) < 1.479 ) ? 0 : barrelSuperClustersSize;
      privateData_->superClusterIndex->push_back( sclusRef.key() + offset );
    } else {
      privateData_->superClusterIndex->push_back( -1 );
    }

    if ( isTrackerDriven && pfclusRef.isNonnull() ) {
      privateData_->PFsuperClusterIndex->push_back( pfclusRef.key() );
    } else {
      privateData_->PFsuperClusterIndex->push_back( -1 );
    }
    
    int packed_corr;
    int isEcalEnergyCorrected = ( electronRef->isEcalEnergyCorrected() ) ? 1 : 0;
    int isMomentumCorrected = ( electronRef->isMomentumCorrected() ) ? 1 : 0;
    packed_corr = ( isEcalEnergyCorrected << 1 ) | isMomentumCorrected;
    privateData_->energyCorrections->push_back( packed_corr );

    // preshower energy
    privateData_->esEnergy->push_back(sclusRef->preshowerEnergy());

  } else {
    privateData_->fiducialFlags->push_back(-1);
    privateData_->recoFlags->push_back(-1);
    privateData_->superClusterIndex->push_back( -1 );
    privateData_->PFsuperClusterIndex->push_back( -1 );
    privateData_->energyCorrections->push_back( -1 );
    privateData_->esEnergy->push_back(-1.);
  }

}

void CmsElectronFiller::treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"fiducialFlags"+colSuffix).c_str(), *privateData_->fiducialFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"recoFlags"+colSuffix).c_str(), *privateData_->recoFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energyCorrections"+colSuffix).c_str(), *privateData_->energyCorrections, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"esEnergy"+colSuffix).c_str(), *privateData_->esEnergy, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"superClusterIndex"+colSuffix).c_str(), *privateData_->superClusterIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFsuperClusterIndex"+colSuffix).c_str(), *privateData_->PFsuperClusterIndex, nCandString.c_str(), 0, "Reco");

}




void CmsElectronFillerData::initialise() {
  
  initialiseCandidate();

  fiducialFlags = new vector<int>;
  recoFlags = new vector<int>;
  esEnergy = new vector<float>;
  energyCorrections = new vector<int>;

  superClusterIndex = new vector<int>;
  PFsuperClusterIndex = new vector<int>;
  trackIndex = new vector<int>;
  gsfTrackIndex = new vector<int>;

  PFChi2EcalTrack = new vector<float>;
  PFesEneL1 = new vector<float>;
  PFesEneL2 = new vector<float>;
  PFesChi2L1 = new vector<float>;
  PFesChi2L2 = new vector<float>;
}

void CmsElectronFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  fiducialFlags->clear();
  recoFlags->clear();
  esEnergy->clear();
  energyCorrections->clear();

  superClusterIndex->clear();
  PFsuperClusterIndex->clear();
  trackIndex->clear();
  gsfTrackIndex->clear();

  PFChi2EcalTrack->clear();
  PFesEneL1->clear();
  PFesEneL2->clear();
  PFesChi2L1->clear();
  PFesChi2L2->clear();
}

void CmsElectronFiller::PSforTMVA(math::XYZTLorentzVector mom, math::XYZTLorentzVector pos, SuperClusterRef pfclusRef ){
  
  BaseParticlePropagator OutParticle(RawParticle(mom,pos),0.,0.,B_.z());
				     
  OutParticle.propagateToPreshowerLayer1(false);
  if (OutParticle.getSuccess()!=0){

    math::XYZPoint v1=math::XYZPoint(OutParticle.vertex());
    //if ((v1.Rho() >=
    // PFGeometry::innerRadius(PFGeometry::PS1)) &&
    //(v1.Rho() <=
    // PFGeometry::outerRadius(PFGeometry::PS1))) {
    if ((v1.Rho() >= 45.0) && (v1.Rho() <= 125.0) ) {

      float enPScl1=0;
      float chi1=100;
      for (CaloCluster_iterator ips=pfclusRef->preshowerClustersBegin(); ips!=pfclusRef->preshowerClustersEnd(); ips++){
	
	ESDetId firstId = (ESDetId)((*ips)->hitsAndFractions()[0].first);

	if (firstId.plane()==2) continue;   // chiara check

	float ax=((*ips)->position().x()-v1.x())/0.114;
	float ay=((*ips)->position().y()-v1.y())/2.43;
	float pschi= sqrt(ax*ax+ay*ay);
	if (pschi<chi1){
	  chi1=pschi;
	  enPScl1=(*ips)->energy();
	}
      }
      ps1En=enPScl1;
      ps1chi=chi1;

      
      OutParticle.propagateToPreshowerLayer2(false);
      if (OutParticle.getSuccess()!=0){
	math::XYZPoint v2=math::XYZPoint(OutParticle.vertex());
	//if ((v2.Rho() >=
	//     PFGeometry::innerRadius(PFGeometry::PS2)) &&
	//    (v2.Rho() <=
	//     PFGeometry::outerRadius(PFGeometry::PS2))){
	if ((v2.Rho() >= 45.0) && (v2.Rho() <= 125.0) ) {
	  float enPScl2=0;
	  float chi2=100;

	  for (CaloCluster_iterator ips=pfclusRef->preshowerClustersBegin(); ips!=pfclusRef->preshowerClustersEnd(); ips++){

	    ESDetId firstId = (ESDetId)((*ips)->hitsAndFractions()[0].first);
	    if (firstId.plane()==1) continue;   // check
	    
	    float ax=((*ips)->position().x()-v2.x())/1.88;
	    float ay=((*ips)->position().y()-v2.y())/0.1449;
	    float pschi= sqrt(ax*ax+ay*ay);
	    if (pschi<chi2){
	      chi2=pschi;
	      enPScl2=(*ips)->energy();
	    }
	  }
	  
	  ps2En=enPScl2;
	  ps2chi=chi2;
	}
      }
    }
  }
}
