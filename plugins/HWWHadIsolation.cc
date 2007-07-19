#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWHadIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HadIsolation.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"

#include <iostream>
#include <algorithm>

using namespace reco;
using namespace std;

HWWHadIsolation::HWWHadIsolation(const edm::ParameterSet& conf)
{
  hcalrhitsLabel_ = conf.getParameter<string>("hcalrhits");
  radiusConeExt_  = conf.getParameter<double>("radiusConeExt");
  radiusConeInt_ = conf.getParameter<double>("radiusConeInt");
  eTMin_ = conf.getParameter<double>("eTMin");
  cut_ = conf.getParameter<double>("cut");
}  
  
HWWHadIsolation::~HWWHadIsolation()
{
}

void HWWHadIsolation::select(edm::Handle<reco::PixelMatchGsfElectronCollection> c, 
                             const edm::Event& e, const edm::EventSetup& es )
{

  std::cout << "[HWWHadIsolation::select] nbr of initial electrons : " << c->size() << std::endl;
  for( reco::PixelMatchGsfElectronCollection::const_iterator el = c->begin();  el != c->end(); el++ ) {
    std::cout << "[HWWHadIsolation::select] new initial electron with pt = " << (el)->pt() << " and eta = " <<
    (el)->eta() << std::endl;
  }  

  selected_.clear();

  //get the HCAL rechit collection
  edm::Handle<HBHERecHitCollection> hbhe;
  HBHERecHitMetaCollection *mhbhe = 0;
  e.getByLabel(hcalrhitsLabel_, hbhe);//getByType(hbhe);  
  mhbhe =  &HBHERecHitMetaCollection(*hbhe);  //FIXME, generates warning

  //services
  //get calo geometry
  es.get<IdealGeometryRecord>().get(theCaloGeom_);
  //product the geometry
  theCaloGeom_.product() ;

  for( reco::PixelMatchGsfElectronCollection::const_iterator el = c->begin(); 
   el != c->end(); el++ ) {

    HadIsolation myHadIsolation(theCaloGeom_,mhbhe,&(*el)) ;
    myHadIsolation.setExtRadius (radiusConeExt_) ; 
    myHadIsolation.setEtLow (eTMin_) ;
    double isoValue = myHadIsolation.getEtHadClusters ()/(el)->pt() ;
    if ( isoValue < cut_ ) 
    	selected_.push_back( & * el );
  }  
    
  std::cout << "[HWWHadIsolation::select] nbr of selected electrons : " << selected_.size() << std::endl;
  for( const_iterator el = begin();  el != end(); el++ ) {
    std::cout << "[HWWHadIsolation::select] new selected electron with pt = " << (*el)->pt() << " and eta = " <<
    (*el)->eta() << std::endl;
  }  

}

