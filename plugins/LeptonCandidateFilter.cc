// my includes
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/VectorUtil.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/LeptonCandidateFilter.h"

#include <iostream>
#include <algorithm>

LeptonCandidateFilter::LeptonCandidateFilter (const edm::ParameterSet& conf)
{
  m_ElectronLabel = conf.getParameter<edm::InputTag>("ElectronLabel");
  m_MuonLabel = conf.getParameter<edm::InputTag>("MuonLabel");
  m_PhotonLabel = conf.getParameter<edm::InputTag>("PhotonLabel");
}


// ------------------------------------------------------------


void 
LeptonCandidateFilter::select (edm::Handle<collection> input, 
                               const edm::Event& iEvent,  const edm::EventSetup& evtSetup ) 
{

  m_selected.clear() ;

  edm::LogInfo("LeptonCandidateFilter") << "LeptonCandidateFilter starting, original collection size = " 
					 << input->size();

  // get the GsfElectronCollection
  edm::Handle< edm::View<reco::Candidate> > electrons;
  try { iEvent.getByLabel(m_ElectronLabel, electrons); }
  catch ( cms::Exception& ex ) { edm::LogWarning("LeptonCandidateFilter") << "Can't get electron candidate collection: " << m_ElectronLabel; }

  // get the MuonCollection
  edm::Handle< edm::View<reco::Candidate> > muons;
  try { iEvent.getByLabel(m_MuonLabel, muons); }
  catch ( cms::Exception& ex ) { edm::LogWarning("LeptonCandidateFilter") << "Can't get muon candidate collection: " << m_MuonLabel; }

  // get the PhotonCollection
  edm::Handle< edm::View<reco::Candidate> > photons;
  try { iEvent.getByLabel(m_PhotonLabel, photons); }
  catch ( cms::Exception& ex ) { edm::LogWarning("LeptonCandidateFilter") << "Can't get photon candidate collection: " << m_PhotonLabel; }

  for(unsigned int itrk=0; itrk<input->size(); itrk++) {

    bool used=false;
    reco::PFCandidateRef candRef(input, itrk);

    // look if it is in a cone dr=0.5 from an electron
    for(unsigned int iele = 0; iele < electrons->size(); iele++) {
      const reco::GsfElectronRef electronRef = electrons->refAt(iele).castTo<reco::GsfElectronRef>();
      if ( !(electronRef.isNull()) && fabs(ROOT::Math::VectorUtil::DeltaR(electronRef->p4(),candRef->p4())) <= 0.5 ) {
        m_selected.push_back(candRef);
        used=true;
      }
    }

    // look if it is in a cone dr=0.5 from a muon
    for(unsigned int imu = 0; imu < muons->size(); imu++) {
      const reco::MuonRef muonRef = muons->refAt(imu).castTo<reco::MuonRef>();
      if ( !(muonRef.isNull()) && fabs(ROOT::Math::VectorUtil::DeltaR(muonRef->p4(),candRef->p4())) <= 0.5 ) {
        if(!used) {
          m_selected.push_back(candRef);
          used=true;
        }
      }
    }

    // look if it is in a cone dr=0.5 from a photon
    for(unsigned int ipho = 0; ipho < photons->size(); ipho++) {
      const reco::PhotonRef photonRef = photons->refAt(ipho).castTo<reco::PhotonRef>();
      if ( !(photonRef.isNull()) && fabs(ROOT::Math::VectorUtil::DeltaR(photonRef->p4(),candRef->p4())) <= 0.5 ) {
        if(!used) {
          m_selected.push_back(candRef);
          used=true;
        }
      }
    }

  }

  edm::LogInfo("LeptonCandidateFilter") << "LeptonCandidateFilter ending, final collection size = " 
                                          << m_selected.size();
  
  return ;
}       



