import FWCore.ParameterSet.Config as cms

from RecoEgamma.ElectronIdentification.electronIdSequence_cff import *

from RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi
egammaIDLikelihood = RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi.eidLikelihoodExt.clone()

eIdSequence = cms.Sequence( eIdSequence * egammaIDLikelihood )
