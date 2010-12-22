#jetCorrectionsData = 1

import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *
#from JetMETCorrections.Configuration.JetCorrectionCondDB_cff import *

# ESSources for JPT corrections from file (until JPT corrections are not in DB)
from HiggsAnalysis.HiggsToWW2e.jptL2L3Corrections_cff import *

from RecoJets.Configuration.RecoPFJets_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *

# to run the offset corrections on PFJets and L2L3 on top of them
PUcorrAk5PFJets=ak5PFJets.clone()
PUcorrAk5PFJets.doAreaFastjet = True
PUcorrKt6PFJets=kt6PFJets.clone()
PUcorrKt6PFJets.doRhoFastjet = True
PUcorrL1Fastjet=L1Fastjet.clone();
PUcorrL1Fastjet.algorithm = cms.string('AK5Calo') #DUMMY THESE are DUMMY,
PUcorrL1Fastjet.era = 'Spring10'                  #DUMMY
PUcorrL1Fastjet.level = cms.string('L2Relative')  #DUMMY
PUcorrL1Fastjet.useCondDB = cms.untracked.bool(False)
PUcorrL1Fastjet.srcMedianPt = 'PUcorrKt6PFJets'
PUcorrAk5PFJetsL1 = ak5PFJetsL1.clone()
PUcorrAk5PFJetsL1.src = 'PUcorrAk5PFJets'
PUcorrAk5PFJetsL1.correctors = ['PUcorrL1Fastjet']

PUcorrAk5PFJetsL1L2L3Residual   = ak5PFJetsL2L3.clone(src = 'PUcorrAk5PFJetsL1', correctors = ['ak5PFL2L3Residual'])
PUcorrAk5PFJetsL1L2L3  = ak5PFJetsL2L3.clone(src = 'PUcorrAk5PFJetsL1', correctors = ['ak5PFL2L3'])

offsetCorrection = cms.Sequence(PUcorrAk5PFJets * PUcorrKt6PFJets * PUcorrAk5PFJetsL1)

# data sequences use residual corrections
CaloJetSequenceData = cms.Sequence( ak5CaloJetsL2L3Residual )
PFJetAK5SequenceData = cms.Sequence( ak5PFJetsL2L3Residual * offsetCorrection * PUcorrAk5PFJetsL1L2L3Residual )
JPTjetsAK5SequenceData = cms.Sequence( ak5JPTJetsL2L3Residual ) # not run for the moment
ourJetSequenceData = cms.Sequence( CaloJetSequenceData * PFJetAK5SequenceData )

# MC sequeces use only L2L3 corrections
CaloJetSequenceMC = cms.Sequence( ak5CaloJetsL2L3 )
PFJetAK5SequenceMC = cms.Sequence( ak5PFJetsL2L3 * offsetCorrection * PUcorrAk5PFJetsL1L2L3 )
JPTjetsAK5SequenceMC = cms.Sequence( ak5JPTJetsL2L3 ) # not run for the moment
ourJetSequenceMC = cms.Sequence( CaloJetSequenceMC * PFJetAK5SequenceMC )

