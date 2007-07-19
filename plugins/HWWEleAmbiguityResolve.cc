// my includes
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWEleAmbiguityResolve.h"

//CLHEP
#include <CLHEP/Vector/LorentzVector.h>

#include <iostream>

HWWEleAmbiguityResolve::HWWEleAmbiguityResolve (const edm::ParameterSet& conf)
{}


// ------------------------------------------------------------


void 
HWWEleAmbiguityResolve::select (edm::Handle<collection> inputHandle, 
                        const edm::Event& evt, const edm::EventSetup& evtStp ) 
{
  m_selected.clear () ;
  //PG get the actual product

  dump (&m_selected, inputHandle) ;

  std::vector<const reco::PixelMatchGsfElectron *>::iterator ShorterEnd;

  sort (m_selected.begin (), m_selected.end (), TKeEoPSorting ()) ; 
  ShorterEnd = unique (m_selected.begin (), m_selected.end (), equalTk ()) ; 
  m_selected.erase(ShorterEnd, m_selected.end());
  sort (m_selected.begin (), m_selected.end (), SCEeEoPSorting ()) ; 
  ShorterEnd = unique (m_selected.begin (), m_selected.end (), equalSC ()) ; 
  m_selected.erase (ShorterEnd, m_selected.end());

  return ;
}       


// ------------------------------------------------------------


HWWEleAmbiguityResolve::~HWWEleAmbiguityResolve()
{}


// ------------------------------------------------------------


void 
HWWEleAmbiguityResolve::dump (container * output, 
                              const edm::Handle<collection> & input) 
{
  for (collection::const_iterator it = input->begin () ;
       it != input->end () ;
       ++it)
    { output->push_back (&(*it)) ; }   
  return ;
}

