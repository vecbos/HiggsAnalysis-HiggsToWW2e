#jetCorrectionsData = 1

import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *
#from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *

# ESSources for JPT corrections from file (until JPT corrections are not in DB)
from HiggsAnalysis.HiggsToWW2e.jptL2L3Corrections_cff import *

from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.Configuration.RecoJets_cff   import *

from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *

##-------------------- Turn-on the FastJet density calculation -----------------------
from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJets= kt4PFJets.clone( rParam = 0.6, doRhoFastjet = cms.bool(True), doAreaFastjet = cms.bool(True) )

#kt6PFJets.doRhoFastjet = True
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
ak5PFJets.doAreaFastjet = True
ak5CaloJets.doAreaFastjet = True

offsetCorrection = cms.Sequence(kt6PFJets)
offsetCaloCorrection = cms.Sequence(kt6PFJets)

############## PF jets PF no PU
# produce PFnoPU jets
from RecoJets.JetProducers.ak5PFJets_cfi import *
ak5PFNoPUJets = ak5PFJets.clone( src = 'pfNoPileUp' )

# calculate rho from this
kt6PFJetsNoPU = kt4PFJets.clone( src = 'pfNoPileUp', rParam = 0.6, doRhoFastjet = cms.bool(True), doAreaFastjet = cms.bool(True) )

# uncorrected jet sequence
FastjetForPFNoPU = cms.Sequence( kt6PFJetsNoPU * ak5PFNoPUJets )


#ak5PFNoPUL1Fastjet = ak5PFL1Fastjet.clone(
#    algorithm = "AK5PFchs",
#    srcRho = cms.InputTag('kt6PFJetsNoPU', 'rho')
#    )


# correct the PFnoPU jets
ak5PFNoPUL1Fastjet = ak5PFL1Fastjet.clone( srcRho = cms.InputTag('kt6PFJets','rho') )
ak5PFNoPUL1FastL2L3 = ak5PFL2L3.clone()
ak5PFNoPUL1FastL2L3.correctors.insert(0,'ak5PFNoPUL1Fastjet')

ak5PFNoPUL1FastL2L3Residual = ak5PFL2L3Residual.clone()
ak5PFNoPUL1FastL2L3Residual.correctors.insert(0,'ak5PFNoPUL1Fastjet')

ak5PFNoPUJetsL2L3 = ak5PFJetsL2L3.clone( src = 'ak5PFNoPUJets', correctors = ['ak5PFL2L3'] )
ak5PFNoPUJetsL2L3Residual = ak5PFJetsL2L3.clone( src = 'ak5PFNoPUJets', correctors = ['ak5PFL2L3Residual'] )

ak5PFNoPUJetsL1FastL2L3 = ak5PFJetsL1FastL2L3.clone( src = 'ak5PFNoPUJets', correctors = ['ak5PFNoPUL1FastL2L3'] )
ak5PFNoPUJetsL1FastL2L3Residual = ak5PFJetsL1FastL2L3.clone( src = 'ak5PFNoPUJets', correctors = ['ak5PFNoPUL1FastL2L3Residual'] )
#ak5PFNoPUJetsL1FastL2L3 = ak5PFJetsL1FastL2L3.clone( srcRho = cms.InputTag('kt6PFJetsNoPU','rho'), src = 'ak5PFNoPUJets', correctors = ['ak5PFL1FastL2L3'] )
#ak5PFNoPUJetsL1FastL2L3Residual = ak5PFJetsL1FastL2L3.clone( srcRho = cms.InputTag('kt6PFJetsNoPU','rho'), src = 'ak5PFNoPUJets', correctors = ['ak5PFL1FastL2L3Residual'] )
####################################

# data sequences use residual corrections
# ak5CaloJetsL1FastL2L3Residual not yet in GlobalTag
#CaloJetSequenceData = cms.Sequence( ak5CaloJetsL2L3Residual* offsetCaloCorrection* ak5CaloJetsL1FastL2L3Residual)                   
#PFJetAK5SequenceData = cms.Sequence( ak5PFJetsL2L3Residual * offsetCorrection * ak5PFJetsL1FastL2L3Residual)
#PFNoPUJetAK5SequenceData = cms.Sequence( FastjetForPFNoPU * ak5PFNoPUJetsL1FastL2L3Residual)
#JPTjetsAK5SequenceData = cms.Sequence( ak5JPTJetsL2L3Residual ) # not run for the moment

CaloJetSequenceData = cms.Sequence( ak5CaloJets * kt6PFJets * ak5CaloJetsL1FastL2L3Residual)   # not run for the moment                 
PFJetAK5SequenceData = cms.Sequence( ak5PFJets * kt6PFJets * ak5PFJetsL1FastL2L3Residual)
PFNoPUJetAK5SequenceData = cms.Sequence( FastjetForPFNoPU * ak5PFNoPUJetsL1FastL2L3Residual)
JPTjetsAK5SequenceData = cms.Sequence( ak5JPTJetsL2L3Residual ) # not run for the moment

ourJetSequenceData = cms.Sequence( PFJetAK5SequenceData * PFNoPUJetAK5SequenceData)
ourJetSequenceDataReduced = cms.Sequence( PFJetAK5SequenceData * PFNoPUJetAK5SequenceData * CaloJetSequenceData )

# MC sequeces use only L2L3 corrections
CaloJetSequenceMC = cms.Sequence( ak5CaloJets * kt6PFJets * ak5CaloJetsL1FastL2L3)  # not run for the moment
PFJetAK5SequenceMC = cms.Sequence( ak5PFJets * kt6PFJets * ak5PFJetsL1FastL2L3 )
PFNoPUJetAK5SequenceMC = cms.Sequence( FastjetForPFNoPU * ak5PFNoPUJetsL1FastL2L3)
JPTjetsAK5SequenceMC = cms.Sequence( ak5JPTJetsL2L3 ) # not run for the moment
ourJetSequenceMC = cms.Sequence( PFJetAK5SequenceMC * PFNoPUJetAK5SequenceMC)
ourJetSequenceMCReduced = cms.Sequence( PFJetAK5SequenceMC * PFNoPUJetAK5SequenceMC * CaloJetSequenceMC )

