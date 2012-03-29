//---------------------------------------------------------------------------
//
// Description:
//       Package:   HiggsAnalysis/HiggsToWW2e
//       Class:     CmsMuonFiller
//
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
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsCandidateFiller.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/CmsMuonFiller.h"

#include <TTree.h>

#include <string>

using namespace edm;
using namespace reco;
using namespace muon;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


CmsMuonFiller::CmsMuonFiller(CmsTree *cmsTree, 
			     int maxTracks, int maxMCTracks,
			     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsMuonFillerData)
{
  cmstree=cmsTree;

  saveMuonExtras_=true;
  saveTrk_=true;
  saveFatTrk_=true;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();

  
}

CmsMuonFiller::CmsMuonFiller(CmsTree *cmsTree, 
			     bool fatTree, 
			     int maxTracks, int maxMCTracks,
			     bool noOutputIfLimitsReached):
  CmsCandidateFiller(cmsTree,maxTracks,fatTree,maxMCTracks,noOutputIfLimitsReached),
  privateData_(new CmsMuonFillerData)
{
  cmstree=cmsTree;

  saveMuonExtras_=true;
  saveTrk_=true;
  saveFatTrk_=fatTree;

  trkIndexName_ = new std::string("n");

  hitLimitsMeansNoOutput_ = noOutputIfLimitsReached;
  maxTracks_=maxTracks;
  maxMCTracks_=maxMCTracks;

  privateData_->initialise();
}


//--------------
// Destructor --
//--------------

CmsMuonFiller::~CmsMuonFiller() { 
  
  // Track information
  delete privateData_->trackIndex;
  delete privateData_->standAloneTrackIndex;
  delete privateData_->combinedTrackIndex;

  delete privateData_->muonId;
  delete privateData_->type;
  delete privateData_->numberOfMatches;

  delete privateData_->  sumPt03;
  delete privateData_->  emEt03;
  delete privateData_->  hadEt03;
  delete privateData_->  hoEt03;
  delete privateData_->  nTrk03;
  delete privateData_->  nJets03;

  delete privateData_->  sumPt05;
  delete privateData_->  emEt05;
  delete privateData_->  hadEt05;
  delete privateData_->  hoEt05;
  delete privateData_->  nTrk05;
  delete privateData_->  nJets05;
  delete privateData_->  kink;

  delete privateData_->EcalExpDepo;
  delete privateData_->HcalExpDepo;
  delete privateData_->HoExpDepo;
  delete privateData_->emS9;
  delete privateData_->hadS9;
  delete privateData_->hoS9;
  delete privateData_->CaloComp;

  delete privateData_->pfCombinedIso;
  delete privateData_->pfCandChargedIso01;
  delete privateData_->pfCandNeutralIso01;
  delete privateData_->pfCandPhotonIso01;
  delete privateData_->pfCandChargedIso02;
  delete privateData_->pfCandNeutralIso02;
  delete privateData_->pfCandPhotonIso02;
  delete privateData_->pfCandChargedIso03;
  delete privateData_->pfCandNeutralIso03;
  delete privateData_->pfCandPhotonIso03;
  delete privateData_->pfCandChargedIso04;
  delete privateData_->pfCandNeutralIso04;
  delete privateData_->pfCandPhotonIso04;
  delete privateData_->pfCandChargedIso05;
  delete privateData_->pfCandNeutralIso05;
  delete privateData_->pfCandPhotonIso05;
  delete privateData_->pfCandChargedIso06;
  delete privateData_->pfCandNeutralIso06;
  delete privateData_->pfCandPhotonIso06;
  delete privateData_->pfCandChargedIso07;
  delete privateData_->pfCandNeutralIso07;
  delete privateData_->pfCandPhotonIso07;
  delete privateData_->pfCandChargedDirIso04;
  delete privateData_->pfCandNeutralDirIso04;
  delete privateData_->pfCandPhotonDirIso04;

  delete privateData_->ncand;

}


//-------------
// Methods   --
//-------------

// Set boolean control options for quantities that are written out

void CmsMuonFiller::saveMuonExtras(bool what) { saveMuonExtras_=what; }

void CmsMuonFiller::saveTrk(bool what) { saveTrk_=what;}

void CmsMuonFiller::saveFatTrk(bool what) { saveFatTrk_=what;}

void CmsMuonFiller::writeCollectionToTree(edm::InputTag collectionTag,
					 const edm::Event& iEvent, const edm::EventSetup& iSetup,
					 const std::string &columnPrefix, const std::string &columnSuffix,
					 bool dumpData) {

  edm::Handle< edm::View<reco::Candidate> > collectionHandle;
  try { iEvent.getByLabel(collectionTag, collectionHandle); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMuonFiller") << "Can't get candidate collection: " << collectionTag; }
  const edm::View<reco::Candidate> *collection = collectionHandle.product();

  privateData_->clearTrkVectors();

  if(collection) {
    if(hitLimitsMeansNoOutput_ && 
       (int)collection->size() > maxTracks_){
      edm::LogInfo("CmsMuonFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ << " and no output flag is set."
			       << " No tracks written to tuple for this event ";
      return;
    }

    if((int)collection->size() > maxTracks_){
      edm::LogInfo("CmsMuonFiller") << "Track length " << collection->size() 
			       << " is too long for declared max length for tree "
			       << maxTracks_ 
			       << ". Collection will be truncated ";
    }
  
    *(privateData_->ncand) = collection->size();

    // for track link
    try { iEvent.getByLabel(generalTracks_, h_tracks); }
    catch ( cms::Exception& ex ) { edm::LogWarning("CmsMuonFiller") << "Can't get general track collection: " << generalTracks_; }

    //read the PF isolation with PF candidates
    eIsoFromPFCandsValueMap_ = new isoContainer(25);
    iEvent.getByLabel( "muonCombinedPFIsoMapProducer", (*eIsoFromPFCandsValueMap_)[0] ); 
    iEvent.getByLabel( "muonPFIsoChHad01", (*eIsoFromPFCandsValueMap_)[1] );
    iEvent.getByLabel( "muonPFIsoNHad01", (*eIsoFromPFCandsValueMap_)[2] );
    iEvent.getByLabel( "muonPFIsoPhoton01", (*eIsoFromPFCandsValueMap_)[3] );
    iEvent.getByLabel( "muonPFIsoChHad02", (*eIsoFromPFCandsValueMap_)[4] );
    iEvent.getByLabel( "muonPFIsoNHad02", (*eIsoFromPFCandsValueMap_)[5] );
    iEvent.getByLabel( "muonPFIsoPhoton02", (*eIsoFromPFCandsValueMap_)[6] );
    iEvent.getByLabel( "muonPFIsoChHad03", (*eIsoFromPFCandsValueMap_)[7] );
    iEvent.getByLabel( "muonPFIsoNHad03", (*eIsoFromPFCandsValueMap_)[8] );
    iEvent.getByLabel( "muonPFIsoPhoton03", (*eIsoFromPFCandsValueMap_)[9] );
    iEvent.getByLabel( "muonPFIsoChHad04", (*eIsoFromPFCandsValueMap_)[10] );
    iEvent.getByLabel( "muonPFIsoNHad04", (*eIsoFromPFCandsValueMap_)[11] );
    iEvent.getByLabel( "muonPFIsoPhoton04", (*eIsoFromPFCandsValueMap_)[12] );
    iEvent.getByLabel( "muonPFIsoChHad05", (*eIsoFromPFCandsValueMap_)[13] );
    iEvent.getByLabel( "muonPFIsoNHad05", (*eIsoFromPFCandsValueMap_)[14] );
    iEvent.getByLabel( "muonPFIsoPhoton05", (*eIsoFromPFCandsValueMap_)[15] );
    iEvent.getByLabel( "muonPFIsoChHad06", (*eIsoFromPFCandsValueMap_)[16] );
    iEvent.getByLabel( "muonPFIsoNHad06", (*eIsoFromPFCandsValueMap_)[17] );
    iEvent.getByLabel( "muonPFIsoPhoton06", (*eIsoFromPFCandsValueMap_)[18] );
    iEvent.getByLabel( "muonPFIsoChHad07", (*eIsoFromPFCandsValueMap_)[19] );
    iEvent.getByLabel( "muonPFIsoNHad07", (*eIsoFromPFCandsValueMap_)[20] );
    iEvent.getByLabel( "muonPFIsoPhoton07", (*eIsoFromPFCandsValueMap_)[21] );
    iEvent.getByLabel( "muonDirPFIsoChHad04", (*eIsoFromPFCandsValueMap_)[22] );
    iEvent.getByLabel( "muonDirPFIsoNHad04", (*eIsoFromPFCandsValueMap_)[23] );
    iEvent.getByLabel( "muonDirPFIsoPhoton04", (*eIsoFromPFCandsValueMap_)[24] );

    int imu=0;
    edm::View<reco::Candidate>::const_iterator cand;
    for(cand=collection->begin(); cand!=collection->end(); cand++) {
      // fill basic kinematics
      if(saveCand_) writeCandInfo(&(*cand),iEvent,iSetup);

      // fill muon extra information
      const reco::Muon *muon = dynamic_cast< const reco::Muon *> ( &(*cand));
      const reco::MuonRef muonRef = collection->refAt(imu).castTo<reco::MuonRef>();

      if(saveMuonExtras_) writeMuonInfo(&(*cand),iEvent,iSetup,&(*muon),muonRef);

      // fill tracks extra informations (only if the muon has a tracker track)
      if(saveTrk_) writeTrkInfo(&(*cand),iEvent,iSetup,&(*muon));
      imu++;
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
  if(saveTrk_) treeTrkInfo(columnPrefix,columnSuffix);
  if(saveMuonExtras_) treeMuonInfo(columnPrefix,columnSuffix);

  if(dumpData) cmstree->dumpData();
	
  delete trkIndexName_;

}


void CmsMuonFiller::writeTrkInfo(const Candidate *cand, 
				 const edm::Event& iEvent, const edm::EventSetup& iSetup,
				 const Muon *muon) {

  if( muon ) {
    TrackRef track = muon->track();
    TrackRef standAloneTrack = muon->standAloneMuon();
    TrackRef combinedTrack = muon->combinedMuon();
    
    int ctfTrackLink=-1;
    if ( track.isNonnull() && h_tracks.isValid() ) {
      for(int itrk=0; itrk<(int)h_tracks->size(); itrk++) {
        reco::TrackRef iTrackRef = h_tracks->at(itrk);
        if(track.key() == iTrackRef.key()) {
          ctfTrackLink = itrk;
          break;
        }
      }
    }
    
    privateData_->trackIndex->push_back(ctfTrackLink);
    privateData_->standAloneTrackIndex->push_back(standAloneTrack.key());
    privateData_->combinedTrackIndex->push_back(combinedTrack.key());
  }

}

void CmsMuonFiller::treeTrkInfo(const std::string &colPrefix, const std::string &colSuffix) {

  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
  
  cmstree->column((colPrefix+"trackIndex"+colSuffix).c_str(), *privateData_->trackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"standAloneTrackIndex"+colSuffix).c_str(), *privateData_->standAloneTrackIndex, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"combinedTrackIndex"+colSuffix).c_str(), *privateData_->combinedTrackIndex, nCandString.c_str(), 0, "Reco");

}

void CmsMuonFiller::writeMuonInfo(const Candidate *cand, const edm::Event& iEvent, 
				  const edm::EventSetup& iSetup, const Muon *muon, const reco::MuonRef muonRef) {
  if(muon) {

    // now put muon ID codified:
    int AllGlobalMuons = ( muon::isGoodMuon( *muon, muon::AllGlobalMuons ) ) ? 1 : 0; 
    int AllStandAloneMuons = ( muon::isGoodMuon( *muon, muon::AllStandAloneMuons ) ) ? 1 : 0;
    int AllTrackerMuons = ( muon::isGoodMuon( *muon, muon::AllTrackerMuons ) ) ? 1 : 0;
    int TrackerMuonArbitrated = ( muon::isGoodMuon( *muon, muon::TrackerMuonArbitrated ) ) ? 1 : 0;
    int AllArbitrated = ( muon::isGoodMuon( *muon, muon::AllArbitrated ) ) ? 1 : 0;
    int GlobalMuonPromptTight = ( muon::isGoodMuon( *muon, muon::GlobalMuonPromptTight ) ) ? 1 : 0;
    int TMLastStationLoose = ( muon::isGoodMuon( *muon, muon::TMLastStationLoose ) ) ? 1 : 0;
    int TMLastStationTight = ( muon::isGoodMuon( *muon, muon::TMLastStationTight ) ) ? 1 : 0;
    int TM2DCompatibilityLoose = ( muon::isGoodMuon( *muon, muon::TM2DCompatibilityLoose ) ) ? 1 : 0;
    int TM2DCompatibilityTight = ( muon::isGoodMuon( *muon, muon::TM2DCompatibilityTight ) ) ? 1 : 0;
    int TMOneStationLoose = ( muon::isGoodMuon( *muon, muon::TMOneStationLoose ) ) ? 1 : 0;
    int TMOneStationTight = ( muon::isGoodMuon( *muon, muon::TMOneStationTight ) ) ? 1 : 0;
    int TMLastStationOptimizedLowPtLoose = ( muon::isGoodMuon( *muon, muon::TMLastStationOptimizedLowPtLoose ) ) ? 1 : 0;
    int TMLastStationOptimizedLowPtTight = ( muon::isGoodMuon( *muon, muon::TMLastStationOptimizedLowPtTight ) ) ? 1 : 0;

    int packed_sel = ( AllGlobalMuons << 13 ) | ( AllStandAloneMuons << 12 ) | ( AllTrackerMuons << 11 ) |
      ( TrackerMuonArbitrated << 10 ) | ( AllArbitrated << 9 ) | 
      ( GlobalMuonPromptTight << 8 ) | 
      ( TMLastStationLoose << 7 ) | ( TMLastStationTight << 6 ) |
      ( TM2DCompatibilityLoose << 5 ) | ( TM2DCompatibilityTight << 4 ) |
      ( TMOneStationLoose << 3 ) | ( TMOneStationTight << 2 ) |
      ( TMLastStationOptimizedLowPtLoose << 1 ) | TMLastStationOptimizedLowPtTight;

    privateData_->muonId->push_back(packed_sel);
    privateData_->type->push_back(muon->type());
    privateData_->numberOfMatches->push_back(muon->numberOfMatches());

    // default isolation variables 0.3
    MuonIsolation Iso03  = muon->isolationR03();
    privateData_->sumPt03->push_back(Iso03.sumPt);
    privateData_->emEt03->push_back(Iso03.emEt);
    privateData_->hadEt03->push_back(Iso03.hadEt);
    privateData_->hoEt03->push_back(Iso03.hoEt);
    privateData_->nTrk03->push_back(Iso03.nTracks);
    privateData_->nJets03->push_back(Iso03.nJets);

    // default isolation variables 0.5
    MuonIsolation Iso05  = muon->isolationR05();
    privateData_->sumPt05->push_back(Iso05.sumPt);
    privateData_->emEt05->push_back(Iso05.emEt);
    privateData_->hadEt05->push_back(Iso05.hadEt);
    privateData_->hoEt05->push_back(Iso05.hoEt);
    privateData_->nTrk05->push_back(Iso05.nTracks);
    privateData_->nJets05->push_back(Iso05.nJets);

    // PF isolation
    const isoFromPFCandsMap & muonsPfCombinedIsoVal = *( (*eIsoFromPFCandsValueMap_)[0] );
    const isoFromPFCandsMap & muonsPfCandChHad01IsoVal = *( (*eIsoFromPFCandsValueMap_)[1] );
    const isoFromPFCandsMap & muonsPfCandNHad01IsoVal = *( (*eIsoFromPFCandsValueMap_)[2] );
    const isoFromPFCandsMap & muonsPfCandPhoton01IsoVal = *( (*eIsoFromPFCandsValueMap_)[3] );
    const isoFromPFCandsMap & muonsPfCandChHad02IsoVal = *( (*eIsoFromPFCandsValueMap_)[4] );
    const isoFromPFCandsMap & muonsPfCandNHad02IsoVal = *( (*eIsoFromPFCandsValueMap_)[5] );
    const isoFromPFCandsMap & muonsPfCandPhoton02IsoVal = *( (*eIsoFromPFCandsValueMap_)[6] );
    const isoFromPFCandsMap & muonsPfCandChHad03IsoVal = *( (*eIsoFromPFCandsValueMap_)[7] );
    const isoFromPFCandsMap & muonsPfCandNHad03IsoVal = *( (*eIsoFromPFCandsValueMap_)[8] );
    const isoFromPFCandsMap & muonsPfCandPhoton03IsoVal = *( (*eIsoFromPFCandsValueMap_)[9] );
    const isoFromPFCandsMap & muonsPfCandChHad04IsoVal = *( (*eIsoFromPFCandsValueMap_)[10] );
    const isoFromPFCandsMap & muonsPfCandNHad04IsoVal = *( (*eIsoFromPFCandsValueMap_)[11] );
    const isoFromPFCandsMap & muonsPfCandPhoton04IsoVal = *( (*eIsoFromPFCandsValueMap_)[12] );
    const isoFromPFCandsMap & muonsPfCandChHad05IsoVal = *( (*eIsoFromPFCandsValueMap_)[13] );
    const isoFromPFCandsMap & muonsPfCandNHad05IsoVal = *( (*eIsoFromPFCandsValueMap_)[14] );
    const isoFromPFCandsMap & muonsPfCandPhoton05IsoVal = *( (*eIsoFromPFCandsValueMap_)[15] );
    const isoFromPFCandsMap & muonsPfCandChHad06IsoVal = *( (*eIsoFromPFCandsValueMap_)[16] );
    const isoFromPFCandsMap & muonsPfCandNHad06IsoVal = *( (*eIsoFromPFCandsValueMap_)[17] );
    const isoFromPFCandsMap & muonsPfCandPhoton06IsoVal = *( (*eIsoFromPFCandsValueMap_)[18] );
    const isoFromPFCandsMap & muonsPfCandChHad07IsoVal = *( (*eIsoFromPFCandsValueMap_)[19] );
    const isoFromPFCandsMap & muonsPfCandNHad07IsoVal = *( (*eIsoFromPFCandsValueMap_)[20] );
    const isoFromPFCandsMap & muonsPfCandPhoton07IsoVal = *( (*eIsoFromPFCandsValueMap_)[21] );
    const isoFromPFCandsMap & muonsPfCandChHad04DirIsoVal = *( (*eIsoFromPFCandsValueMap_)[22] );
    const isoFromPFCandsMap & muonsPfCandNHad04DirIsoVal = *( (*eIsoFromPFCandsValueMap_)[23] );
    const isoFromPFCandsMap & muonsPfCandPhoton04DirIsoVal = *( (*eIsoFromPFCandsValueMap_)[24] );
    
    privateData_->pfCombinedIso->push_back( muonsPfCombinedIsoVal[muonRef] );
    privateData_->pfCandChargedIso01->push_back( muonsPfCandChHad01IsoVal[muonRef] );
    privateData_->pfCandNeutralIso01->push_back( muonsPfCandNHad01IsoVal[muonRef] );
    privateData_->pfCandPhotonIso01->push_back( muonsPfCandPhoton01IsoVal[muonRef] );
    privateData_->pfCandChargedIso02->push_back( muonsPfCandChHad02IsoVal[muonRef] );
    privateData_->pfCandNeutralIso02->push_back( muonsPfCandNHad02IsoVal[muonRef] );
    privateData_->pfCandPhotonIso02->push_back( muonsPfCandPhoton02IsoVal[muonRef] );
    privateData_->pfCandChargedIso03->push_back( muonsPfCandChHad03IsoVal[muonRef] );
    privateData_->pfCandNeutralIso03->push_back( muonsPfCandNHad03IsoVal[muonRef] );
    privateData_->pfCandPhotonIso03->push_back( muonsPfCandPhoton03IsoVal[muonRef] );
    privateData_->pfCandChargedIso04->push_back( muonsPfCandChHad04IsoVal[muonRef] );
    privateData_->pfCandNeutralIso04->push_back( muonsPfCandNHad04IsoVal[muonRef] );
    privateData_->pfCandPhotonIso04->push_back( muonsPfCandPhoton04IsoVal[muonRef] );
    privateData_->pfCandChargedIso05->push_back( muonsPfCandChHad05IsoVal[muonRef] );
    privateData_->pfCandNeutralIso05->push_back( muonsPfCandNHad05IsoVal[muonRef] );
    privateData_->pfCandPhotonIso05->push_back( muonsPfCandPhoton05IsoVal[muonRef] );
    privateData_->pfCandChargedIso06->push_back( muonsPfCandChHad06IsoVal[muonRef] );
    privateData_->pfCandNeutralIso06->push_back( muonsPfCandNHad06IsoVal[muonRef] );
    privateData_->pfCandPhotonIso06->push_back( muonsPfCandPhoton06IsoVal[muonRef] );
    privateData_->pfCandChargedIso07->push_back( muonsPfCandChHad07IsoVal[muonRef] );
    privateData_->pfCandNeutralIso07->push_back( muonsPfCandNHad07IsoVal[muonRef] );
    privateData_->pfCandPhotonIso07->push_back( muonsPfCandPhoton07IsoVal[muonRef] );
    privateData_->pfCandChargedDirIso04->push_back( muonsPfCandChHad04DirIsoVal[muonRef] );
    privateData_->pfCandNeutralDirIso04->push_back( muonsPfCandNHad04DirIsoVal[muonRef] );
    privateData_->pfCandPhotonDirIso04->push_back( muonsPfCandPhoton04DirIsoVal[muonRef] );

    // track kinks
    privateData_->kink->push_back(muonRef->combinedQuality().trkKink) ;

    // Expected deposits in CALO
    privateData_->EcalExpDepo->push_back(muon->calEnergy().em);
    privateData_->HcalExpDepo->push_back(muon->calEnergy().had);
    privateData_->HoExpDepo->push_back(muon->calEnergy().ho);
    privateData_->emS9->push_back(muon->calEnergy().emS9);
    privateData_->hadS9->push_back(muon->calEnergy().hadS9);
    privateData_->hoS9->push_back(muon->calEnergy().hoS9);
    privateData_->CaloComp->push_back(muon->caloCompatibility());

  } else {

    privateData_->muonId->push_back(0);
    privateData_->type->push_back(-1);
    privateData_->numberOfMatches->push_back(-1);
    // default isolation variables 0.3
    privateData_->sumPt03->push_back(-1.);
    privateData_->emEt03->push_back(-1.);
    privateData_->hadEt03->push_back(-1.);
    privateData_->hoEt03->push_back(-1.);
    privateData_->nTrk03->push_back(-1.);
    privateData_->nJets03->push_back(-1.);
    // default isolation variables 0.5
     privateData_->sumPt05->push_back(-1.);
    privateData_->emEt05->push_back(-1.);
    privateData_->hadEt05->push_back(-1.);
    privateData_->hoEt05->push_back(-1.);
    privateData_->nTrk05->push_back(-1.);
    privateData_->nJets05->push_back(-1.);

    privateData_->pfCombinedIso->push_back( -1 );
    privateData_->pfCandChargedIso01->push_back( -1 );
    privateData_->pfCandNeutralIso01->push_back( -1 );
    privateData_->pfCandPhotonIso01->push_back( -1 );
    privateData_->pfCandChargedIso02->push_back( -1 );
    privateData_->pfCandNeutralIso02->push_back( -1 );
    privateData_->pfCandPhotonIso02->push_back( -1 );
    privateData_->pfCandChargedIso03->push_back( -1 );
    privateData_->pfCandNeutralIso03->push_back( -1 );
    privateData_->pfCandPhotonIso03->push_back( -1 );
    privateData_->pfCandChargedIso04->push_back( -1 );
    privateData_->pfCandNeutralIso04->push_back( -1 );
    privateData_->pfCandPhotonIso04->push_back( -1 );
    privateData_->pfCandChargedIso05->push_back( -1 );
    privateData_->pfCandNeutralIso05->push_back( -1 );
    privateData_->pfCandPhotonIso05->push_back( -1 );
    privateData_->pfCandChargedIso06->push_back( -1 );
    privateData_->pfCandNeutralIso06->push_back( -1 );
    privateData_->pfCandPhotonIso06->push_back( -1 );
    privateData_->pfCandChargedIso07->push_back( -1 );
    privateData_->pfCandNeutralIso07->push_back( -1 );
    privateData_->pfCandPhotonIso07->push_back( -1 );
    privateData_->pfCandChargedDirIso04->push_back( -1 );
    privateData_->pfCandNeutralDirIso04->push_back( -1 );
    privateData_->pfCandPhotonDirIso04->push_back( -1 );

    privateData_->kink->push_back( -1 );

    // Expected deposits in CALO
    privateData_->EcalExpDepo->push_back(-1.);
    privateData_->HcalExpDepo->push_back(-1.);
    privateData_->HoExpDepo->push_back(-1.);
    privateData_->emS9->push_back(-1.);
    privateData_->hadS9->push_back(-1.);
    privateData_->hoS9->push_back(-1.);
    privateData_->CaloComp->push_back(-1.);

  }
}

void CmsMuonFiller::treeMuonInfo(const std::string &colPrefix, const std::string &colSuffix) {
  std::string nCandString=colPrefix+(*trkIndexName_)+colSuffix;
   
  // muon id 
  cmstree->column((colPrefix+"muonId"+colSuffix).c_str(), *privateData_->muonId, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"type"+colSuffix).c_str(), *privateData_->type, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"numberOfMatches"+colSuffix).c_str(), *privateData_->numberOfMatches, nCandString.c_str(), 0, "Reco");

  // isolation R=0.3
  cmstree->column((colPrefix+"sumPt03"+colSuffix).c_str(), *privateData_->sumPt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emEt03"+colSuffix).c_str(), *privateData_->emEt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadEt03"+colSuffix).c_str(), *privateData_->hadEt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hoEt03"+colSuffix).c_str(), *privateData_->hoEt03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nTrk03"+colSuffix).c_str(), *privateData_->nTrk03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nJets03"+colSuffix).c_str(), *privateData_->nJets03, nCandString.c_str(), 0, "Reco");

  // isolation R=0.5
  cmstree->column((colPrefix+"sumPt05"+colSuffix).c_str(), *privateData_->sumPt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emEt05"+colSuffix).c_str(), *privateData_->emEt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadEt05"+colSuffix).c_str(), *privateData_->hadEt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hoEt05"+colSuffix).c_str(), *privateData_->hoEt05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nTrk05"+colSuffix).c_str(), *privateData_->nTrk05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"nJets05"+colSuffix).c_str(), *privateData_->nJets05, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"pfCombinedIso"+colSuffix).c_str(),  *privateData_->pfCombinedIso, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso01"+colSuffix).c_str(),  *privateData_->pfCandChargedIso01, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso01"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso01, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso01"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso01, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso02"+colSuffix).c_str(),  *privateData_->pfCandChargedIso02, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso02"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso02, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso02"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso02, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso03"+colSuffix).c_str(),  *privateData_->pfCandChargedIso03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso03"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso03"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso03, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso04"+colSuffix).c_str(),  *privateData_->pfCandChargedIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso04"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso04"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso05"+colSuffix).c_str(),  *privateData_->pfCandChargedIso05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso05"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso05"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso05, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso06"+colSuffix).c_str(),  *privateData_->pfCandChargedIso06, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso06"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso06, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso06"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso06, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedIso07"+colSuffix).c_str(),  *privateData_->pfCandChargedIso07, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralIso07"+colSuffix).c_str(),  *privateData_->pfCandNeutralIso07, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonIso07"+colSuffix).c_str(),  *privateData_->pfCandPhotonIso07, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandChargedDirIso04"+colSuffix).c_str(),  *privateData_->pfCandChargedDirIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandNeutralDirIso04"+colSuffix).c_str(),  *privateData_->pfCandNeutralDirIso04, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"pfCandPhotonDirIso04"+colSuffix).c_str(),  *privateData_->pfCandPhotonDirIso04, nCandString.c_str(), 0, "Reco");

  cmstree->column((colPrefix+"kink"+colSuffix).c_str(),  *privateData_->kink, nCandString.c_str(), 0, "Reco");

  //  Expected deposits in CALO
  cmstree->column((colPrefix+"EcalExpDepo"+colSuffix).c_str(), *privateData_->EcalExpDepo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HcalExpDepo"+colSuffix).c_str(), *privateData_->HcalExpDepo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"HoExpDepo"+colSuffix).c_str(), *privateData_->HoExpDepo, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"emS9"+colSuffix).c_str(), *privateData_->emS9, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hadS9"+colSuffix).c_str(), *privateData_->hadS9, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"hoS9"+colSuffix).c_str(), *privateData_->hoS9, nCandString.c_str(), 0, "Reco");
  cmstree->column((colPrefix+"CaloComp"+colSuffix).c_str(), *privateData_->CaloComp, nCandString.c_str(), 0, "Reco");

}

void CmsMuonFillerData::initialise() {

  initialiseCandidate();

  trackIndex = new vector<int>;
  standAloneTrackIndex = new vector<int>;
  combinedTrackIndex = new vector<int>;

  muonId = new vector<int>;
  type = new vector<int>;
  numberOfMatches = new vector<int>;

  sumPt03 = new vector<float>;
  emEt03 = new vector<float>;
  hadEt03 = new vector<float>;
  hoEt03 = new vector<float>;
  nTrk03 = new vector<float>;
  nJets03 = new vector<float>;
  sumPt05 = new vector<float>;
  emEt05 = new vector<float>;
  hadEt05 = new vector<float>;
  hoEt05 = new vector<float>;
  nTrk05 = new vector<float>;
  nJets05 = new vector<float>;

  pfCombinedIso      = new vector<float>;
  pfCandChargedIso01 = new vector<float>;
  pfCandNeutralIso01 = new vector<float>;
  pfCandPhotonIso01  = new vector<float>;
  pfCandChargedIso02 = new vector<float>;
  pfCandNeutralIso02 = new vector<float>;
  pfCandPhotonIso02  = new vector<float>;
  pfCandChargedIso03 = new vector<float>;
  pfCandNeutralIso03 = new vector<float>;
  pfCandPhotonIso03  = new vector<float>;
  pfCandChargedIso04 = new vector<float>;
  pfCandNeutralIso04 = new vector<float>;
  pfCandPhotonIso04  = new vector<float>;
  pfCandChargedIso05 = new vector<float>;
  pfCandNeutralIso05 = new vector<float>;
  pfCandPhotonIso05  = new vector<float>;
  pfCandChargedIso06 = new vector<float>;
  pfCandNeutralIso06 = new vector<float>;
  pfCandPhotonIso06  = new vector<float>;
  pfCandChargedIso07 = new vector<float>;
  pfCandNeutralIso07 = new vector<float>;
  pfCandPhotonIso07  = new vector<float>;
  pfCandChargedDirIso04 = new vector<float>;
  pfCandNeutralDirIso04 = new vector<float>;
  pfCandPhotonDirIso04  = new vector<float>;

  kink = new vector<float>;

  EcalExpDepo = new vector<float>;
  HcalExpDepo = new vector<float>;
  HoExpDepo = new vector<float>;
  emS9 = new vector<float>;
  hadS9 = new vector<float>;
  hoS9 = new vector<float>;
  CaloComp = new vector<float>;
  
}

void CmsMuonFillerData::clearTrkVectors() {

  clearTrkVectorsCandidate();

  trackIndex->clear();
  standAloneTrackIndex->clear();
  combinedTrackIndex->clear();
 
  muonId->clear();
  type->clear();
  numberOfMatches->clear();

  sumPt03->clear();
  emEt03->clear();
  hadEt03->clear();
  hoEt03->clear();
  nTrk03->clear();
  nJets03->clear();

  sumPt05->clear();
  emEt05->clear();
  hadEt05->clear();
  hoEt05->clear();
  nTrk05->clear();
  nJets05->clear();

  pfCombinedIso    ->clear();
  pfCandChargedIso01 ->clear();
  pfCandNeutralIso01 ->clear();
  pfCandPhotonIso01  ->clear();
  pfCandChargedIso02 ->clear();
  pfCandNeutralIso02 ->clear();
  pfCandPhotonIso02  ->clear();
  pfCandChargedIso03 ->clear();
  pfCandNeutralIso03 ->clear();
  pfCandPhotonIso03  ->clear();
  pfCandChargedIso04 ->clear();
  pfCandNeutralIso04 ->clear();
  pfCandPhotonIso04  ->clear();
  pfCandChargedIso05 ->clear();
  pfCandNeutralIso05 ->clear();
  pfCandPhotonIso05  ->clear();
  pfCandChargedIso06 ->clear();
  pfCandNeutralIso06 ->clear();
  pfCandPhotonIso06  ->clear();
  pfCandChargedIso07 ->clear();
  pfCandNeutralIso07 ->clear();
  pfCandPhotonIso07  ->clear();
  pfCandChargedDirIso04 ->clear();
  pfCandNeutralDirIso04 ->clear();
  pfCandPhotonDirIso04  ->clear();

  kink->clear();

  EcalExpDepo->clear();
  HcalExpDepo->clear();
  HoExpDepo->clear();
  emS9->clear();
  hadS9->clear();
  hoS9->clear();
  CaloComp->clear();
  
}
