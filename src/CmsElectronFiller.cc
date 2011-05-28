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
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

//Kinematic constraint vertex fitter
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/ColinearityKinematicConstraint.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsEleIDTreeFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsElectronFiller.h"

#include <TTree.h>

#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/Exception.h"

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

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();

}


//--------------
// Destructor --
//--------------

CmsElectronFiller::~CmsElectronFiller() {

  // delete here the vector ptr's
  delete privateData_->fiducialFlags;
  delete privateData_->recoFlags;
  delete privateData_->energyCorrections;
  delete privateData_->scPixCharge;

  delete privateData_->superClusterIndex;
  delete privateData_->PFsuperClusterIndex;
  delete privateData_->trackIndex;
  delete privateData_->gsfTrackIndex;

  delete privateData_->convDist;
  delete privateData_->convDcot;
  delete privateData_->convRadius;
  delete privateData_->convTrackIndex;

  delete privateData_->ncand;

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

    // for track link
    try { iEvent.getByLabel(generalTracks_, h_tracks); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get general track collection: " << generalTracks_; }

    // for conversions
    try { iEvent.getByLabel("generalTracks", h_tracksTot); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsElectronFiller") << "Can't get general track collection: generalTracks"; }

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
  if(saveTrk_) treeTrkInfo(columnPrefix,columnSuffix);
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

    int ctfTrackLink=-1;
    reco::TrackRef closeCtfTrack = electronRef->closestCtfTrackRef();
    if ( closeCtfTrack.isNonnull() ) {
      if(h_tracks.isValid()) {
        for(int itrk=0; itrk<(int)h_tracks->size(); itrk++) {
          reco::TrackRef iTrackRef = h_tracks->at(itrk);
          if(closeCtfTrack.key() == iTrackRef.key()) ctfTrackLink = itrk;
        }
      }
      privateData_->trackIndex->push_back(ctfTrackLink);
      privateData_->scPixCharge->push_back(electronRef->scPixCharge());
    } else {
      privateData_->trackIndex->push_back( -1 );
      privateData_->scPixCharge->push_back( -999. );
    }

    // for conversion rejection
    if( h_tracks.isValid() ) {
      edm::ESHandle<MagneticField> theMagField;
      iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
      const MagneticField *magField=&(*theMagField);

      GlobalPoint origin(0.,0.,0.);
      ConversionFinder finder;
      ConversionInfo convInfo = finder.getConversionInfo(*electronRef,
                                                         h_tracksTot,
                                                         (magField->inTesla(origin)).mag());
      
      privateData_->convDist->push_back(convInfo.dist());
      privateData_->convDcot->push_back(convInfo.dcot());
      privateData_->convRadius->push_back(convInfo.radiusOfConversion());

      if( convInfo.conversionPartnerCtfTk().isNonnull() )
        privateData_->convTrackIndex->push_back(convInfo.conversionPartnerCtfTk().key());
      else privateData_->convTrackIndex->push_back(-1);
      
    } else {
      privateData_->convDist->push_back(999);
      privateData_->convDcot->push_back(999);
      privateData_->convRadius->push_back(999);
      privateData_->convTrackIndex->push_back(-1);
    }

  } else {
    privateData_->gsfTrackIndex->push_back( -1 );
    privateData_->trackIndex->push_back( -1 );
    privateData_->convDist->push_back(999);
    privateData_->convDcot->push_back(999);
    privateData_->convRadius->push_back(999);
    privateData_->convTrackIndex->push_back(-1);
  }
    
}




void CmsElectronFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;

  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(),  *privateData_->trackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"gsfTrackIndex"+colSuffix).c_str(),  *privateData_->gsfTrackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"convDist"+colSuffix).c_str(),  *privateData_->convDist, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"convDcot"+colSuffix).c_str(),  *privateData_->convDcot, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"convRadius"+colSuffix).c_str(),  *privateData_->convRadius, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"convTrackIndex"+colSuffix).c_str(),  *privateData_->convTrackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"scPixCharge"+colSuffix).c_str(),  *privateData_->scPixCharge, nCandString.c_str(), 0, "Reco");

}





void CmsElectronFiller::writeEcalInfo(const GsfElectronRef electronRef, 
				      const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				      SuperClusterRef sclusRef, SuperClusterRef pfclusRef,
				      const EcalRecHitCollection *EBRecHits,
				      const EcalRecHitCollection *EERecHits) {

  if(electronRef.isNonnull()) {

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
    int isEcalDriven = ( electronRef->ecalDrivenSeed() ) ? 1 : 0;
    int isTrackerDriven = ( electronRef->trackerDrivenSeed() ) ? 1 : 0;
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
    packed_corr = isEcalEnergyCorrected;
    privateData_->energyCorrections->push_back( packed_corr );

  } else {
    privateData_->fiducialFlags->push_back(-1);
    privateData_->recoFlags->push_back(-1);
    privateData_->superClusterIndex->push_back( -1 );
    privateData_->PFsuperClusterIndex->push_back( -1 );
    privateData_->energyCorrections->push_back( -1 );
  }

}

void CmsElectronFiller::treeEcalInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString = colPrefix+(*trkIndexName_)+colSuffix;
  cmstree->column((colPrefix+"fiducialFlags"+colSuffix).c_str(), *privateData_->fiducialFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"recoFlags"+colSuffix).c_str(), *privateData_->recoFlags, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"energyCorrections"+colSuffix).c_str(), *privateData_->energyCorrections, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"superClusterIndex"+colSuffix).c_str(), *privateData_->superClusterIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"PFsuperClusterIndex"+colSuffix).c_str(), *privateData_->PFsuperClusterIndex, nCandString.c_str(), 0, "Reco");

}




void CmsElectronFillerData::initialise() {
  
  initialiseCandidate();

  fiducialFlags = new vector<int>;
  recoFlags = new vector<int>;
  energyCorrections = new vector<int>;

  superClusterIndex = new vector<int>;
  PFsuperClusterIndex = new vector<int>;
  trackIndex = new vector<int>;
  gsfTrackIndex = new vector<int>;
  convDist = new vector<float>;
  convDcot = new vector<float>;
  convRadius = new vector<float>;
  convTrackIndex = new vector<int>;

  scPixCharge = new vector<int>;

}

void CmsElectronFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  fiducialFlags->clear();
  recoFlags->clear();
  energyCorrections->clear();

  superClusterIndex->clear();
  PFsuperClusterIndex->clear();
  trackIndex->clear();
  gsfTrackIndex->clear();
  convDist->clear();
  convDcot->clear();
  convRadius->clear();
  convTrackIndex->clear();
  
  scPixCharge->clear();

}
