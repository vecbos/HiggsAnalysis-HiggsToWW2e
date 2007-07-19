#ifndef TkIsolation_h
#define TkIsolation_h

//C++ includes
#include <vector>
#include <functional>

//Root includes
#include "TObjArray.h"

//CMSSW includes 
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"


class TkIsolation {
 public:
  
  //constructors
  TkIsolation () ;
  TkIsolation (const reco::PixelMatchGsfElectron* , 
	       const reco::TrackCollection* ) ;
  TkIsolation (const reco::PixelMatchGsfElectron* , 
	       const reco::TrackCollection* ,
	       const reco::PixelMatchGsfElectronCollection* ) ;

  //methods
  void setExtRadius (double extRadius) ;
  void setIntRadius (double intRadius) ;
  void setPtLow (double ptLow) ;
  void setLip (double lip) ;

  int getNumberTracks() const ;
  double getPtTracks () const ;
/*
  int getNumberECinsideCone () const ;
  double getPtMax () const ;
  bool isIsolated (double ptCut = 0.2) const ;
 */
  //destructor 
  ~TkIsolation() ;
  
 private:
  
  const reco::PixelMatchGsfElectron  *electron_ ;
  const reco::TrackCollection *trackCollection_ ;
  const reco::PixelMatchGsfElectronCollection *electronCollection_ ;
  
  double extRadius_ ;
  double intRadius_ ;
  double ptLow_ ;
  double lip_ ;
};

#endif
