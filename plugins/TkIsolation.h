#ifndef TkIsolation_h
#define TkIsolation_h

//C++ includes
#include <vector>
#include <functional>

//Root includes
#include "TObjArray.h"

//CMSSW includes 
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

class TkIsolation {
 public:
  
  //constructors
  TkIsolation () ;
  TkIsolation (const reco::GsfElectron* , 
	       const reco::TrackCollection* ) ;
  TkIsolation (const reco::GsfElectron* , 
	       const reco::TrackCollection* ,
	       const reco::GsfElectronCollection* ) ;

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
  
  const reco::GsfElectron  *electron_ ;
  const reco::TrackCollection *trackCollection_ ;
  const reco::GsfElectronCollection *electronCollection_ ;
  
  double extRadius_ ;
  double intRadius_ ;
  double ptLow_ ;
  double lip_ ;
};

#endif
