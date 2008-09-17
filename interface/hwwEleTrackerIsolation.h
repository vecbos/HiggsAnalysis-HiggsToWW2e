// -*- C++ -*-
//
// Package:    HtoWWElectrons
// Class:      hwwEleTrackerIsolation
// 
/*
   Description: <one line class summary>
   Electron isolation using tracker info

   Implementation:
 
*/
//
// Original Author:  Chiara Rovelli
//
//


#ifndef hwwEleTrackerIsolation_h
#define hwwEleTrackerIsolation_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

using namespace edm;
using namespace std;
using namespace reco;

class hwwEleTrackerIsolation{
 public:
  
  //constructors
  hwwEleTrackerIsolation();
  hwwEleTrackerIsolation(const GsfElectron *gsfEle, const TrackCollection trackColl);

  //methods
  void setExtRadius (float extRadius);
  void setIntRadius (float intRadius);
  
  float getPtTracks (bool relative=true, bool squared=false) const;
  float getNTracks (float minPtTrack=1.0) const;
  float minDeltaR (float minPt) const;
  float minDeltaR_withVeto (float minPt) const;
  bool isIsolated (float ptCut = 0.05) const;

  //destructor 
  ~hwwEleTrackerIsolation();
  
 private:

  const GsfElectron*_myGsfEle;  	  
  const TrackCollection _tracks;
  
  float _extRadius;
  float _intRadius;
};

#endif
