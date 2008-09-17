import FWCore.ParameterSet.Config as cms

from RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi import *

import RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi
egammaIDCutsLoose = RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi.eidCutBasedClassesExt.clone()

from RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi import *

import RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi
egammaIDLikelihood = RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi.eidLikelihoodExt.clone()

eIdSequence = cms.Sequence( egammaIDCutsLoose + egammaIDLikelihood )
