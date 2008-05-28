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
float hwwEleTrackerIsolation::getPtTracks (bool relative, bool squared) const
{
  float dummyPt = 0 ;

  Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  float ele_pt  = _myGsfEle->pt();
  float ele_lip = _myGsfEle->vz();   

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    
    Hep3Vector trackPAtVtx(this_track->px(),this_track->py(),this_track->pz());
    float this_pt  = trackPAtVtx.perp();

    // only tracks from the same vertex as the electron
    float this_lip = this_track->vz();
    if ( fabs(this_lip - ele_lip) > 0.2 ){ continue; }

    // usually the < 1GeV tracks are fakes 
    if ( this_pt < 1.0 ) continue;

    double dr = elePAtVtx.deltaR(trackPAtVtx);
    if ( fabs(dr) < _extRadius && fabs(dr) > _intRadius ){ 
      if ( relative ) {
	if ( squared ) dummyPt += pow(this_pt/ele_pt,2);
	else dummyPt += this_pt/ele_pt;
      }
      else {
	if ( squared ) dummyPt += pow(this_pt,2);
	else dummyPt += this_pt;
      }
    } 
    
  } //end loop over tracks		       
  
  return dummyPt;
}

float hwwEleTrackerIsolation::getNTracks (float minPtTrack) const
{
  int dummyN = 0 ;

  Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  float ele_lip = _myGsfEle->vz();   

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    
    Hep3Vector trackPAtVtx(this_track->px(),this_track->py(),this_track->pz());
    float this_pt  = trackPAtVtx.perp();

    // only tracks from the same vertex as the electron
    float this_lip = this_track->vz();
    if ( fabs(this_lip - ele_lip) > 0.2 ){ continue; }

    // only tracks with not negligible pT enter
    if ( this_pt < minPtTrack ) continue;

    double dr = elePAtVtx.deltaR(trackPAtVtx);
    if ( fabs(dr) < _extRadius && fabs(dr) > _intRadius ){ dummyN++; } 
    
  } //end loop over tracks		       
  
  return dummyN;
}

// minimum distance from tracks upper a given pt cut - without veto
float hwwEleTrackerIsolation::minDeltaR(float minPt) const
{
  float minDR = 100000. ;

  Hep3Vector elePAtVtx(_myGsfEle->px(), _myGsfEle->py(), _myGsfEle->pz()); 
  float ele_lip = _myGsfEle->vz();   

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    Hep3Vector trackPAtVtx(this_track->px(),this_track->py(),this_track->pz());
    float this_pt  = trackPAtVtx.perp();

    if (this_pt < minPt){ continue;} 

    // only tracks from the same vertex as the electron
    float this_lip = this_track->vz();
    if ( fabs(this_lip - ele_lip) > 0.2 ){ continue; }
    
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
  float ele_lip = _myGsfEle->vz();   

  for(TrackCollection::const_iterator this_track = _tracks.begin(); this_track != _tracks.end(); this_track++ ){ 
    
    Hep3Vector trackPAtVtx(this_track->px(),this_track->py(),this_track->pz());
    float this_pt  = trackPAtVtx.perp();
    if (this_pt < minPt){ continue;} 

    // only tracks from the same vertex as the electron
    float this_lip = this_track->vz();
    if ( fabs(this_lip - ele_lip) > 0.2 ){ continue; }
    
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
