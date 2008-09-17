#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTkIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/TkIsolation.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"

#include <iostream>
#include <algorithm>

//ROOT includes
#include <Math/VectorUtil.h>

using namespace reco;
using namespace std;
//using namespace ROOT::Math::VectorUtil ;

HWWTkIsolation::HWWTkIsolation(const edm::ParameterSet& conf)
{

  tracksLabel_ = conf.getParameter<string>("tracks");
  radiusConeExt_  = conf.getParameter<double>("radiusConeExt");
  radiusConeInt_ = conf.getParameter<double>("radiusConeInt");
  pTMin_ = conf.getParameter<double>("pTMin");
  lip_ = conf.getParameter<double>("lip");
  cut_ = conf.getParameter<double>("cut");
  eleIn = 0 ;
  eleOut = 0 ;
}  
  
HWWTkIsolation::~HWWTkIsolation()
{
    std::cout << "eleIn=" << eleIn << std::endl;
    std::cout << "eleOut=" << eleOut << std::endl;
    std::cout << "efficiency=" << double(eleOut)/double(eleIn) << std::endl;
}

void HWWTkIsolation::select(edm::Handle<reco::GsfElectronCollection> c,
                            const edm::Event& e, const edm::EventSetup& evtStp)
{

  std::cout << "[HWWTkIsolation::select] nbr of initial electrons : " << c->size() << std::endl;
  for( reco::GsfElectronCollection::const_iterator el = c->begin();  el != c->end(); el++ ) {
    std::cout << "[HWWTkIsolation::select] new initial electron with pt = " << (el)->pt() << " and eta = " <<
    (el)->eta() << std::endl;
  }  

  selected_.clear();

  // get track collection
  edm::Handle<TrackCollection> tracks;
  e.getByLabel(tracksLabel_, tracks); 
  const reco::TrackCollection* trackCollection = tracks.product () ;
  std::cout << "[HWWTkIsolation::select] trackCollection->size()=" << trackCollection->size() << std::endl;
  
  for( reco::GsfElectronCollection::const_iterator el = c->begin(); 
   el != c->end(); el++ ) {

    TkIsolation myTkIsolation(&(*el),trackCollection) ;
    myTkIsolation.setExtRadius (radiusConeExt_) ; 
    myTkIsolation.setIntRadius (radiusConeInt_) ;
    myTkIsolation.setPtLow (pTMin_) ;
    myTkIsolation.setLip (lip_) ;
    double isoValue = myTkIsolation.getPtTracks()/(el)->pt() ;
    eleIn++;
    std::cout << "isoValue="<< isoValue << std::endl;
    std::cout << "myTkIsolation.getPtTracks()=" << myTkIsolation.getPtTracks() << std::endl ;
    if ( isoValue < cut_ ) 
      {
	eleOut++;
    	selected_.push_back( & * el );
	std::cout << "selected" << std::endl;
      }
  }  
    
  std::cout << "[HWWTkIsolation::select] nbr of selected electrons : " << selected_.size() << std::endl;
  for( const_iterator el = begin();  el != end(); el++ ) {
    std::cout << "[HWWTkIsolation::select] new selected electron with pt = " << (*el)->pt() << " and eta = " <<
    (*el)->eta() << std::endl;
  }  

}

