import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *
#from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *

# ESSources for JPT corrections from file (until JPT corrections are not in DB)
from HiggsAnalysis.HiggsToWW2e.jptL2L3Corrections_cff import *

CaloJetSequence = cms.Sequence( ak5CaloJetsL2L3Residual )
PFJetAK5Sequence = cms.Sequence( ak5PFJetsL2L3Residual )
JPTjetsAK5Sequence = cms.Sequence( ak5JPTJetsL2L3Residual )

ourJetSequence = cms.Sequence( CaloJetSequence * PFJetAK5Sequence * JPTjetsAK5Sequence )

