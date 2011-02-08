#jetCorrectionsData = 1

import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *
#from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *

# ESSources for JPT corrections from file (until JPT corrections are not in DB)
from HiggsAnalysis.HiggsToWW2e.jptL2L3Corrections_cff import *

from RecoJets.Configuration.RecoPFJets_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *

##-------------------- Disable the CondDB for the L1FastJet (until they are included in a new global tag) -------
ak5PFL1Fastjet.useCondDB = False
##-------------------- Turn-on the FastJet density calculation -----------------------
kt6PFJets.doRhoFastjet = True
kt6PFJets.Rho_EtaMax=cms.double(4.5)
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
ak5PFJets.doAreaFastjet = True
ak5PFJets.Rho_EtaMax=cms.double(4.5)

offsetCorrection = cms.Sequence(kt6PFJets * ak5PFJets)

# data sequences use residual corrections
CaloJetSequenceData = cms.Sequence( ak5CaloJetsL2L3Residual )                   
PFJetAK5SequenceData = cms.Sequence( ak5PFJetsL2L3Residual * offsetCorrection * ak5PFJetsL1FastL2L3Residual)
JPTjetsAK5SequenceData = cms.Sequence( ak5JPTJetsL2L3Residual ) # not run for the moment
ourJetSequenceData = cms.Sequence( CaloJetSequenceData * PFJetAK5SequenceData )

# MC sequeces use only L2L3 corrections
CaloJetSequenceMC = cms.Sequence( ak5CaloJetsL2L3 )
PFJetAK5SequenceMC = cms.Sequence( ak5PFJetsL2L3 * offsetCorrection * ak5PFJetsL1FastL2L3 )
JPTjetsAK5SequenceMC = cms.Sequence( ak5JPTJetsL2L3 ) # not run for the moment
ourJetSequenceMC = cms.Sequence( CaloJetSequenceMC * PFJetAK5SequenceMC )

