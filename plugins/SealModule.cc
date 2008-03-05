#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWEleAmbiguityResolve.h"
#include "HiggsAnalysis/HiggsToWW2e/plugins/HWWTreeDumper.h"


typedef ObjectSelector< HWWEleAmbiguityResolve, 
			reco::PixelMatchGsfElectronRefVector
			> AmbResolver ;

DEFINE_SEAL_MODULE () ;
DEFINE_ANOTHER_FWK_MODULE (AmbResolver) ;
DEFINE_ANOTHER_FWK_MODULE (HWWTreeDumper);  
