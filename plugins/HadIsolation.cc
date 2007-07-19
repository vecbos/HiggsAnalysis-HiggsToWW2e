//C++ includes
#include <vector>
#include <functional>

//ROOT includes
#include <Math/VectorUtil.h>

//CMSSW includes
#include "HiggsAnalysis/HiggsToWW2e/plugins/HadIsolation.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"


using namespace std;

HadIsolation::HadIsolation ()
{
}

HadIsolation::HadIsolation ( edm::ESHandle<CaloGeometry> theCaloGeom ,
                             HBHERecHitMetaCollection*  mhbhe,
                             const reco::PixelMatchGsfElectron* electron ) :
   theCaloGeom_(theCaloGeom) ,  
   mhbhe_(mhbhe) ,
   electron_(electron) 
{
  electronCollection_ = 0 ;
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  etLow_ = 1.5 ; 
}

HadIsolation::HadIsolation (edm::ESHandle<CaloGeometry> theCaloGeom ,
                            HBHERecHitMetaCollection*  mhbhe,
			    const reco::PixelMatchGsfElectron* electron , 
			    const reco::PixelMatchGsfElectronCollection* electronCollection ) : 
  theCaloGeom_(theCaloGeom) , 
  mhbhe_(mhbhe) ,
  electron_(electron) ,
  electronCollection_(electronCollection)  
{
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  etLow_ = 1.5 ; 
}  

HadIsolation::~HadIsolation ()
{
}

void HadIsolation::setExtRadius (double extRadius)
{
  extRadius_ = extRadius ;
}

void HadIsolation::setIntRadius (double intRadius)
{  
  intRadius_ = intRadius ;
}

void HadIsolation::setEtLow (double etLow)
{  
  etLow_ = etLow ;
}

double HadIsolation::getEtHadClusters () const
{

  double hcalEt = 0.;
  if (mhbhe_) 
   {
      //Take the SC position
      const CaloGeometry* caloGeom = theCaloGeom_.product();
      CaloConeSelector sel(extRadius_ , caloGeom, DetId::Hcal);
      math::XYZPoint theCaloPosition = electron_->caloPosition () ;
      GlobalPoint pclu (theCaloPosition.x () ,
                	theCaloPosition.y () ,
			theCaloPosition.z () );
      //Compute the HCAL energy behind ECAL
      std::auto_ptr<CaloRecHitMetaCollectionV> chosen = sel.select(pclu,*mhbhe_);
      for (CaloRecHitMetaCollectionV::const_iterator i = chosen->begin () ; 
                                                     i!= chosen->end () ; 
						     ++i) 
       {
	 double hcalHit_eta = caloGeom->getPosition(i->detid()).eta();
	 double hcalHit_Et = i->energy()*sin(2*atan(exp(-hcalHit_eta)));
	 if ( hcalHit_Et > etLow_)
	      hcalEt += i->energy();
       }
    } 
  return hcalEt ;
}

double HadIsolation::getHoE () const
{

  double HoE ;
  if (mhbhe_) 
   {
     //Take the SC position
     const CaloGeometry* caloGeom = theCaloGeom_.product();
     CaloConeSelector sel(extRadius_ , caloGeom, DetId::Hcal);
     math::XYZPoint theCaloPosition = electron_->caloPosition () ;
     GlobalPoint pclu (theCaloPosition.x () ,
                       theCaloPosition.y () ,
		       theCaloPosition.z () );
     //Compute the HCAL energy behind ECAL
     double hcalEnergy = 0. ;
     std::auto_ptr<CaloRecHitMetaCollectionV> chosen = sel.select(pclu,*mhbhe_);
     for (CaloRecHitMetaCollectionV::const_iterator i = chosen->begin () ; 
                                                    i!= chosen->end () ; 
						    ++i) 
     {
       hcalEnergy += i->energy();
     }
     //Take the SC energy
     double ecalEnergy = electron_->caloEnergy () ;
     //Compute HoE
     HoE = hcalEnergy/ecalEnergy ;
   } 
  else HoE = 0. ;

  return HoE ;
}

