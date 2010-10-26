import FWCore.ParameterSet.Config as cms

# to run likelihood and standard CIC
from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentification_cfi import *
eIdSequence = cms.Sequence(eidVeryLoose+eidLoose+eidMedium+eidTight+eidSuperTight+eidHyperTight1+eidHyperTight2+eidHyperTight3)

# to run the HWW V04 CIC
from HiggsAnalysis.HiggsToWW2e.cutsInCategoriesHWWElectronIdentificationV04_cfi import *
eIdHWWSequence = cms.Sequence(eidHWWVeryLoose+eidHWWLoose+eidHWWMedium+eidHWWTight+eidHWWSuperTight+eidHWWHyperTight1+eidHWWHyperTight2+eidHWWHyperTight3)

from RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi
egammaIDLikelihood = RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi.eidLikelihoodExt.clone()

eIdSequence = cms.Sequence( eIdSequence * eIdHWWSequence * egammaIDLikelihood )
