#ifndef LeptonTrackFilter_h
#define LeptonTrackFilter_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <functional>
#include <vector>
#include <map>


class LeptonTrackFilter{

 public:

  typedef const reco::Track * track ;
  typedef reco::TrackCollection collection ;
  typedef reco::TrackRefVector container;
  typedef container::const_iterator const_iterator ;
  
  //! ctor
  LeptonTrackFilter (const edm::ParameterSet& conf) ;
  //!dtor
  ~LeptonTrackFilter () {};

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
  edm::InputTag m_ElectronLabel, m_MuonLabel, m_JetLabel;

};  

#endif

