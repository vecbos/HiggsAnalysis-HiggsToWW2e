// my includes
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "HiggsAnalysis/HiggsToWW2e/interface/hwwEleTrackerIsolation.h"

// temporary version of tracker isolation
// we'd like to compare parameters at vertex, but no way.... take them at innermost state
// also missing: cut on longitudinal impact parameter

//CLHEP
#include <CLHEP/Vector/LorentzVector.h>

using namespace std;
hwwEleTrackerIsolation::hwwEleTrackerIsolation (){}

hwwEleTrackerIsolation::hwwEleTrackerIsolation (const PixelMatchGsfElectron *gsfEle, const TrackCollection trackColl) : 
  _myGsfEle(gsfEle),
  _tracks(trackColl)   
{
  _extRadius = 0.20;
  _intRadius = 0.015;
}

hwwEleTrackerIsolation::~hwwEleTrackerIsolation (){}


void hwwEleTrackerIsolation::setExtRadius (float extRadius){_extRadius = extRadius; }

void hwwEleTrackerIsolation::setIntRadius (float intRadius){_intRadius = intRadius; }

// sum of pt of tracks (within a given cone) / electron pt
float hwwEleTrackerIsolation::getPtTracks () const
{
  float dummyPt = 0 ;

  Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  float ele_pt  = _myGsfEle->pt();
  // float ele_lip = _myGsfEle->lip();    // to be fixed

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    
    Hep3Vector trackPAtVtx(this_track->innerMomentum().x(),this_track->innerMomentum().y(),this_track->innerMomentum().z());
    float this_pt  = this_track->innerMomentum().rho();

    // only tracks from the same vertex as the electron
    // float this_lip = TrackLip[iso];    
    // if ( fabs(this_lip - my_ele_lip) > 0.2 ){ continue; }

    double dr = elePAtVtx.deltaR(trackPAtVtx);
    if ( fabs(dr) < _extRadius && fabs(dr) > _intRadius ){ dummyPt += this_pt; } 
    
  } //end loop over tracks		       
  
  // sum tracks pt / ele pt
  dummyPt = dummyPt/ele_pt; 

  return dummyPt;
}

// minimum distance from tracks upper a given pt cut - without veto
float hwwEleTrackerIsolation::minDeltaR(float minPt) const
{
  float minDR = 100000. ;

  Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  // float ele_lip = _myGsfEle->lip();    // to be fixed

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    
    Hep3Vector trackPAtVtx(this_track->innerMomentum().x(),this_track->innerMomentum().y(),this_track->innerMomentum().z());
    float this_pt  = this_track->innerMomentum().rho();
    if (this_pt < minPt){ continue;} 

    // only tracks from the same vertex as the electron
    // float this_lip = TrackLip[iso];      // to be fixed  
    // if ( fabs(this_lip - my_ele_lip) > 0.2 ){ continue; }
    
    double dr = elePAtVtx.deltaR(trackPAtVtx);
    if ( fabs(dr) < minDR ){ minDR = dr; }
 
  } //end loop over tracks		       
  
  return minDR;
}

// minimum distance from tracks upper a given pt cut - with veto
float hwwEleTrackerIsolation::minDeltaR_withVeto(float minPt) const
{
  float minDR = 100000. ;

  Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  // float ele_lip = _myGsfEle->lip();    // to be fixed

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    
    Hep3Vector trackPAtVtx(this_track->innerMomentum().x(),this_track->innerMomentum().y(),this_track->innerMomentum().z());
    float this_pt  = this_track->innerMomentum().rho();
    if (this_pt < minPt){ continue;} 

    // only tracks from the same vertex as the electron
    // float this_lip = TrackLip[iso];      // to be fixed  
    // if ( fabs(this_lip - my_ele_lip) > 0.2 ){ continue; }
    
    double dr = elePAtVtx.deltaR(trackPAtVtx);
    if ( (fabs(dr) < minDR) && (fabs(dr)>_intRadius) ){ minDR = dr; }
 
  } //end loop over tracks		       
  
  return minDR;
}

bool hwwEleTrackerIsolation::isIsolated (float ptCut) const
{
  bool dummyIsolation = true ;
  
  if (hwwEleTrackerIsolation::getPtTracks() > ptCut )   // 0.05 
    dummyIsolation = false ;
  
  return dummyIsolation ;
}
