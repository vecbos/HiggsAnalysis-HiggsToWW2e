from MyAnalysis.METFlags.EcalDeadCellEventFlagProducer_cfi import *
from MyAnalysis.METFlags.simpleDRFlagProducer_cfi import *

metOptionalFilterSequence = cms.Sequence( EcalDeadCellEventFlagProducer * simpleDRFlagProducer )
