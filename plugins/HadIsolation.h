#ifndef HadIsolation_h
#define HadIsolation_h

//C++ includes
#include <vector>
#include <functional>

//CMSSW includes
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"

class HadIsolation {
 public:
  
  //constructors
  HadIsolation ( ) ;
  HadIsolation (edm::ESHandle<CaloGeometry> ,
                HBHERecHitMetaCollection*  ,
		const reco::PixelMatchGsfElectron* ) ;
  HadIsolation (edm::ESHandle<CaloGeometry> , 
                HBHERecHitMetaCollection*  ,
		const reco::PixelMatchGsfElectron* ,
		const reco::PixelMatchGsfElectronCollection* ) ;
  
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
  const reco::PixelMatchGsfElectron  *electron_ ;
  const reco::PixelMatchGsfElectronCollection *electronCollection_ ;
  
  double extRadius_ ;
  double intRadius_ ;
  double etLow_ ;

};

#endif
