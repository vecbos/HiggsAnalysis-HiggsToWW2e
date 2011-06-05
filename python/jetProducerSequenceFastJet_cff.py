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
kt6PFJets.doRhoFastjet = True
kt6PFJets.Rho_EtaMax=cms.double(4.5)
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
ak5PFJets.doAreaFastjet = True
ak5PFJets.Rho_EtaMax=cms.double(4.5)
ak5CaloJets.doAreaFastjet = True
ak5CaloJets.Rho_EtaMax=cms.double(4.5)

offsetCorrection = cms.Sequence(kt6PFJets)
offsetCaloCorrection = cms.Sequence(kt6PFJets)

############## PF jets PF no PU
# produce PFnoPU jets
from CommonTools.ParticleFlow.pfNoPileUp_cff import *
from RecoJets.JetProducers.ak5PFJets_cfi import *
ak5PFJetsNoPU = ak5PFJets.clone( src = 'pfNoPileUp' )

# calculate rho from this
from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJetsNoPu = kt4PFJets.clone( src = 'pfNoPileUp', rParam = 0.6, doRhoFastjet = True )
kt6PFJetsNoPu.Rho_EtaMax = cms.double(4.5)

# uncorrected jet sequence
FastjetForPFNoPU = cms.Sequence( pfPileUp * pfNoPileUp * kt6PFJetsNoPu * ak5PFJetsNoPU )

# correct the PFnoPU jets
ak5PFJetsNoPUL1FastL2L3 = ak5PFJetsL2L3.clone( src = 'ak5PFJetsNoPU', correctors = ['ak5PFL1FastL2L3'] )
ak5PFJetsNoPUL1FastL2L3Residual = ak5PFJetsL2L3.clone( src = 'ak5PFJetsNoPU', correctors = ['ak5PFL1FastL2L3Residual'] )
####################################

# data sequences use residual corrections
# ak5CaloJetsL1FastL2L3Residual not yet in GlobalTag
#CaloJetSequenceData = cms.Sequence( ak5CaloJetsL2L3Residual* offsetCaloCorrection* ak5CaloJetsL1FastL2L3Residual)                   
#PFJetAK5SequenceData = cms.Sequence( ak5PFJetsL2L3Residual * offsetCorrection * ak5PFJetsL1FastL2L3Residual)
#PFNoPUJetAK5SequenceData = cms.Sequence( FastjetForPFNoPU * ak5PFJetsNoPUL1FastL2L3Residual)
#JPTjetsAK5SequenceData = cms.Sequence( ak5JPTJetsL2L3Residual ) # not run for the moment

CaloJetSequenceData = cms.Sequence( ak5CaloJetsL2L3* offsetCaloCorrection* ak5CaloJetsL1FastL2L3)                   
PFJetAK5SequenceData = cms.Sequence( ak5PFJetsL2L3 * offsetCorrection * ak5PFJetsL1FastL2L3)
PFNoPUJetAK5SequenceData = cms.Sequence( FastjetForPFNoPU * ak5PFJetsNoPUL1FastL2L3)
JPTjetsAK5SequenceData = cms.Sequence( ak5JPTJetsL2L3 ) # not run for the moment

ourJetSequenceData = cms.Sequence( CaloJetSequenceData * PFJetAK5SequenceData * PFNoPUJetAK5SequenceData)
ourJetSequenceDataReduced = cms.Sequence( CaloJetSequenceData * PFJetAK5SequenceData)

# MC sequeces use only L2L3 corrections
CaloJetSequenceMC = cms.Sequence( ak5CaloJetsL2L3* offsetCaloCorrection* ak5CaloJetsL1FastL2L3)
PFJetAK5SequenceMC = cms.Sequence( ak5PFJetsL2L3 * offsetCorrection * ak5PFJetsL1FastL2L3 )
PFNoPUJetAK5SequenceMC = cms.Sequence( FastjetForPFNoPU * ak5PFJetsNoPUL1FastL2L3)
JPTjetsAK5SequenceMC = cms.Sequence( ak5JPTJetsL2L3 ) # not run for the moment
ourJetSequenceMC = cms.Sequence( CaloJetSequenceMC * PFJetAK5SequenceMC * PFNoPUJetAK5SequenceMC)
ourJetSequenceMCReduced = cms.Sequence( CaloJetSequenceMC * PFJetAK5SequenceMC)

