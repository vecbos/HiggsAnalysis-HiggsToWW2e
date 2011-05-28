// my includes
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/LeptonTrackFilter.h"

#include <iostream>
#include <algorithm>

LeptonTrackFilter::LeptonTrackFilter (const edm::ParameterSet& conf)
{
  m_ElectronLabel = conf.getParameter<edm::InputTag>("ElectronLabel");
  m_MuonLabel = conf.getParameter<edm::InputTag>("MuonLabel");
  m_JetLabel = conf.getParameter<edm::InputTag>("JetLabel");
}


// ------------------------------------------------------------


void 
LeptonTrackFilter::select (edm::Handle<collection> input, 
                               const edm::Event& iEvent,  const edm::EventSetup& evtSetup ) 
{

  m_selected.clear() ;

  edm::LogInfo("LeptonTrackFilter") << "LeptonTrackFilter starting, original collection size = " 
					 << input->size();

  // get the GsfElectronCollection
  edm::Handle< edm::View<reco::Candidate> > electrons;
  try { iEvent.getByLabel(m_ElectronLabel, electrons); }
  catch ( cms::Exception& ex ) { edm::LogWarning("LeptonTrackFilter") << "Can't get electron candidate collection: " << m_ElectronLabel; }

  // get the MuonCollection
  edm::Handle< edm::View<reco::Candidate> > muons;
  try { iEvent.getByLabel(m_MuonLabel, muons); }
  catch ( cms::Exception& ex ) { edm::LogWarning("LeptonTrackFilter") << "Can't get electron candidate collection: " << m_MuonLabel; }

  // get the JetCollection
  edm::Handle< edm::View<reco::Candidate> > jets;
  try { iEvent.getByLabel(m_JetLabel, jets); }
  catch ( cms::Exception& ex ) { edm::LogWarning("LeptonTrackFilter") << "Can't get electron candidate collection: " << m_JetLabel; }

  for(unsigned int itrk=0; itrk<input->size(); itrk++) {

    reco::TrackRef trackRef(input, itrk);

    // look if it is linked to an electron
    for(unsigned int iele = 0; iele < electrons->size(); iele++) {
      const reco::GsfElectronRef electronRef = electrons->refAt(iele).castTo<reco::GsfElectronRef>();
      if ( !(electronRef.isNull()) ) {
        reco::TrackRef closeCtfTrack = electronRef->closestCtfTrackRef();
        if ( closeCtfTrack.isNonnull() ) { 
          if( closeCtfTrack.key() == trackRef.key() ) m_selected.push_back(trackRef);
        }
      }
    }

    // look if it is linked to a muon
    for(unsigned int imu = 0; imu < muons->size(); imu++) {
      const reco::MuonRef muonRef = muons->refAt(imu).castTo<reco::MuonRef>();
      if ( !(muonRef.isNull()) ) {
        reco::TrackRef linkCtfTrack = muonRef->track();
        if ( linkCtfTrack.isNonnull() ) { 
          if( linkCtfTrack.key() == trackRef.key() ) m_selected.push_back(trackRef);
        }
      }
    }

    // look if it is in a cone DR < 0.5 wrt jet (temporary for calojet ID)
    edm::View<reco::Candidate>::const_iterator jet;
    for(jet=jets->begin(); jet!=jets->end(); jet++) {
      if(jet->pt()>15. && reco::deltaR(jet->eta(),jet->phi(),trackRef->eta(),trackRef->phi())<0.5) m_selected.push_back(trackRef);
    }

  }

  edm::LogInfo("LeptonTrackFilter") << "LeptonTrackFilter ending, final collection size = " 
                                    << m_selected.size();
  
  return ;
}       



