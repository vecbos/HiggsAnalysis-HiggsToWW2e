from MyAnalysis.METFlags.EcalDeadCellEventFlagProducer_cfi import *
from MyAnalysis.METFlags.simpleDRFlagProducer_cfi import *
# the CSC Halo filter only works in the RECO
#from MyAnalysis.METFlags.CSCHaloFlagProducer_cfi import *
from MyAnalysis.METFlags.trackingFailureFlagProducer_cfi import *

metOptionalFilterSequence = cms.Sequence( EcalDeadCellEventFlagProducer * simpleDRFlagProducer * trackingFailureFlagProducer )
