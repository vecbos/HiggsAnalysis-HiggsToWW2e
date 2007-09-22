#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWEleAmbiguityResolve.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWEleId.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTkIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWHadIsolation.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWMetProducer.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTreeDumper.h"
//PG temp obj sel for 131
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWObjectSelector.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
//#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "PhysicsTools/Utilities/interface/PtComparator.h"
#include "PhysicsTools/UtilAlgos/interface/PtMinSelector.h"
//#include "PhysicsTools/UtilAlgos/interface/SingleElementCollectionSelector.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectCountFilter.h"

typedef HWWObjectSelector<HWWEleAmbiguityResolve> AmbResolver ;
typedef HWWObjectSelector<HWWEleId> EleIdentifier ;
typedef HWWObjectSelector<HWWTkIsolation> TkIso;  
typedef HWWObjectSelector<HWWHadIsolation> HadIso;

DEFINE_SEAL_MODULE () ;
DEFINE_ANOTHER_FWK_MODULE (AmbResolver) ;
DEFINE_ANOTHER_FWK_MODULE (EleIdentifier) ;
DEFINE_ANOTHER_FWK_MODULE (TkIso) ;
DEFINE_ANOTHER_FWK_MODULE (HadIso);  
DEFINE_ANOTHER_FWK_MODULE (HWWMetProducer);  
DEFINE_ANOTHER_FWK_MODULE (HWWTreeDumper);  
