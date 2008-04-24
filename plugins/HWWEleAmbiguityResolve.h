#ifndef HWWEleAmbiguityResolve_h
#define HWWEleAmbiguityResolve_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <functional>
#include <vector>
#include <map>


class HWWEleAmbiguityResolve{

 public:

  typedef const reco::PixelMatchGsfElectron * electron ;
  typedef reco::PixelMatchGsfElectronCollection collection ;
  typedef reco::PixelMatchGsfElectronRefVector container;
  typedef container::const_iterator const_iterator ;
  
  //! ctor
  HWWEleAmbiguityResolve (const edm::ParameterSet& conf) ;
  //!dtor
  ~HWWEleAmbiguityResolve () ;

  //! iterator to the begin of the selected collection
  const_iterator begin () const { return m_selected.begin () ; }
  
  //! iterator to the end of the selected collection
  const_iterator end () const { return m_selected.end () ; }

  //! the actual selector method 
  void select (edm::Handle<collection>, const edm::Event& );
  // for CMSSW >= 1_6_10
  //  void select (edm::Handle<collection>, const edm::Event&, const edm::EventSetup& ) ;
     

 private:

  //! find the non ambiguous electrons, initialise the ambiguity map
  void Init(const edm::Handle<collection> & input);
  
  //! ambiguity resolution
  void ResolveByEoverP(const edm::Handle<collection> & input);

  //! the selected collection
  container m_selected ;
 
  //! map between the ambiguous electrons
  std::vector<std::pair<unsigned int, unsigned int> > ambEle;

  //! if doRefCheck, select only among the elements present in the reducedElectron collection
  bool m_doRefCheck;
  edm::InputTag m_reducedElectronsRefCollectionLabel;
  edm::Handle< edm::RefVector<collection> > m_reduced;

};  

#endif

