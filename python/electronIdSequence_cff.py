import FWCore.ParameterSet.Config as cms

# electron ID cuts egamma robust
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
egammaIDStandardCutsRobust = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
egammaIDStandardCutsRobust.electronQuality = "robust"

# electron ID cuts egamma loose
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
egammaIDStandardCutsLoose = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
egammaIDStandardCutsLoose.electronQuality = "loose"

# electron ID cuts egamma tight
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
egammaIDStandardCutsTight = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
egammaIDStandardCutsTight.electronQuality = "tight"

eIdSequence = cms.Sequence( egammaIDStandardCutsRobust +
                            egammaIDStandardCutsLoose +
                            egammaIDStandardCutsTight )
