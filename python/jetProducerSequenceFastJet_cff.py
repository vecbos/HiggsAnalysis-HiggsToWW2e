#jetCorrectionsData = 1

import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *
#from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *

from JetMETCorrections.Configuration.JetCorrectionServices_cff  import ak5PFL1FastL2L3
from JetMETCorrections.Configuration.JetCorrectionServices_cff  import ak5PFL1FastL2L3Residual

from JetMETCorrections.Configuration.JetCorrectionServices_cff  import ak5CaloL1FastL2L3
from JetMETCorrections.Configuration.JetCorrectionServices_cff  import ak5CaloL1FastL2L3Residual

# ESSources for JPT corrections from file (until JPT corrections are not in DB)
from HiggsAnalysis.HiggsToWW2e.jptL2L3Corrections_cff import *

from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.Configuration.RecoJets_cff   import *

from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *

##-------------------- Turn-on the FastJet density calculation -----------------------
from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJets = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = cms.bool(True), doAreaFastjet = cms.bool(True) )

##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
ak5PFJets.doAreaFastjet = True
ak5CaloJets.doAreaFastjet = True

############## PF jets PF no PU
# produce PFnoPU jets
from RecoJets.JetProducers.ak5PFJets_cfi import *
ak5PFNoPUJets = ak5PFJets.clone( src = 'pfNoPileUp' )

# uncorrected jet sequence
FastjetForPFNoPU = cms.Sequence( ak5PFNoPUJets )


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

CaloJetSequenceData = cms.Sequence( ak5CaloJets)   # not run for the moment
PFJetAK5SequenceData = cms.Sequence( ak5PFJets * kt6PFJets )
PFNoPUJetAK5SequenceData = cms.Sequence( FastjetForPFNoPU )

ourJetSequenceData = cms.Sequence( PFJetAK5SequenceData * PFNoPUJetAK5SequenceData)
ourJetSequenceDataReduced = cms.Sequence( PFJetAK5SequenceData * PFNoPUJetAK5SequenceData * CaloJetSequenceData )

# MC sequeces use only L2L3 corrections
CaloJetSequenceMC = cms.Sequence( ak5CaloJets )  # not run for the moment
PFJetAK5SequenceMC = cms.Sequence( ak5PFJets * kt6PFJets )
PFNoPUJetAK5SequenceMC = cms.Sequence( FastjetForPFNoPU )
ourJetSequenceMC = cms.Sequence( PFJetAK5SequenceMC * PFNoPUJetAK5SequenceMC)
ourJetSequenceMCReduced = cms.Sequence( PFJetAK5SequenceMC * PFNoPUJetAK5SequenceMC * CaloJetSequenceMC )

