#ifndef LeptonCandidateFilter_h
#define LeptonCandidateFilter_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <functional>
#include <vector>
#include <map>


class LeptonCandidateFilter{

 public:

  typedef const reco::PFCandidate * track ;
  typedef reco::PFCandidateCollection collection ;
  typedef reco::PFCandidateRefVector container;
  typedef container::const_iterator const_iterator ;
  
  //! ctor
  LeptonCandidateFilter (const edm::ParameterSet& conf) ;
  //!dtor
  ~LeptonCandidateFilter () {};

  //! iterator to the begin of the selected collection
  const_iterator begin () const { return m_selected.begin () ; }
  
  //! iterator to the end of the selected collection
  const_iterator end () const { return m_selected.end () ; }

  //! the actual selector method 
  void select (edm::Handle<collection>, const edm::Event&, const edm::EventSetup& ) ;
     

 private:

  //! the selected collection
  container m_selected ;

  //! the configurable inputs
  edm::InputTag m_ElectronLabel, m_MuonLabel, m_PhotonLabel;

};  

#endif

