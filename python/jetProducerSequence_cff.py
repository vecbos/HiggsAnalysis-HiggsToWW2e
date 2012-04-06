#jetCorrectionsData = 1

import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *
#from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *

# ESSources for JPT corrections from file (until JPT corrections are not in DB)
from HiggsAnalysis.HiggsToWW2e.jptL2L3Corrections_cff import *

from RecoJets.Configuration.RecoPFJets_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *

# data sequences use residual corrections
CaloJetSequenceData = cms.Sequence( ak5CaloJetsL2L3Residual )   # not run for the moment                 
PFJetAK5SequenceData = cms.Sequence( ak5PFJetsL1L2L3Residual )    
JPTjetsAK5SequenceData = cms.Sequence( ak5JPTJetsL2L3Residual ) # not run for the moment
ourJetSequenceData = cms.Sequence( PFJetAK5SequenceData )

# MC sequeces use only L2L3 corrections
CaloJetSequenceMC = cms.Sequence( ak5CaloJetsL2L3 )   # not run for the moment
PFJetAK5SequenceMC = cms.Sequence( ak5PFJetsL1L2L3 )
JPTjetsAK5SequenceMC = cms.Sequence( ak5JPTJetsL2L3 ) # not run for the moment
ourJetSequenceMC = cms.Sequence( PFJetAK5SequenceMC )

