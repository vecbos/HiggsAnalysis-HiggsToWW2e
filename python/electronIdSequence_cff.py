import FWCore.ParameterSet.Config as cms

# to run likelihood and standard CIC
from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi import *
eIdSequence = cms.Sequence(eidVeryLoose+eidLoose+eidMedium+eidTight+eidSuperTight+eidHyperTight1)

# to run the HWW V04 CIC
from HiggsAnalysis.HiggsToWW2e.cutsInCategoriesHWWElectronIdentificationV04_cfi import *
eIdHWWSequence = cms.Sequence(eidHWWVeryLoose+eidHWWLoose+eidHWWMedium+eidHWWTight+eidHWWSuperTight+eidHWWHyperTight1+eidHWWHyperTight2+eidHWWHyperTight3)

from RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi
egammaIDLikelihood = RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi.eidLikelihoodExt.clone()

eIdSequence = cms.Sequence( eIdSequence * eIdHWWSequence * egammaIDLikelihood )


# to compute FastJet rho to correct isolation (note: EtaMax restricted to 2.5)
from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

FastjetForIsolation = cms.Sequence( kt6PFJetsForIsolation )
