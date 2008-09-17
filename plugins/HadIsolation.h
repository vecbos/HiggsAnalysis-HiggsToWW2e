#ifndef HadIsolation_h
#define HadIsolation_h

//C++ includes
#include <vector>
#include <functional>

//CMSSW includes
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"

class HadIsolation {
 public:
  
  //constructors
  HadIsolation ( ) ;
  HadIsolation (edm::ESHandle<CaloGeometry> ,
                HBHERecHitMetaCollection*  ,
		const reco::GsfElectron* ) ;
  HadIsolation (edm::ESHandle<CaloGeometry> , 
                HBHERecHitMetaCollection*  ,
		const reco::GsfElectron* ,
		const reco::GsfElectronCollection* ) ;
  
  //methods
  void setExtRadius (double extRadius) ;
  void setIntRadius (double intRadius) ;
  void setEtLow (double etLow) ;

  double getEtHadClusters () const ;
  double getHoE () const ;

  //destructor 
  ~HadIsolation() ;
  
 private:
  
  edm::ESHandle<CaloGeometry>  theCaloGeom_ ;
  HBHERecHitMetaCollection* mhbhe_ ;
  const reco::GsfElectron  *electron_ ;
  const reco::GsfElectronCollection *electronCollection_ ;
  
  double extRadius_ ;
  double intRadius_ ;
  double etLow_ ;

};

#endif
