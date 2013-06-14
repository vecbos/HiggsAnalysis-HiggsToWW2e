import FWCore.ParameterSet.Config as cms

# this is the simplest PF isolation (combined)
import WWAnalysis.Tools.electronPFIsoMapProd_cfi
electronCombinedPFIsoMapProducer = WWAnalysis.Tools.electronPFIsoMapProd_cfi.electronPFIsoMapProd.clone()
electronCombinedPFIsoMapProducer.vtxLabel = 'offlinePrimaryVertices' # if the event has the first vertex bad, will be discarded offline.

import WWAnalysis.Tools.muonPFIsoMapProd_cfi
muonCombinedPFIsoMapProducer = WWAnalysis.Tools.muonPFIsoMapProd_cfi.muonPFIsoMapProd.clone()
muonCombinedPFIsoMapProducer.vtxLabel = 'offlinePrimaryVertices' # if the event has the first vertex bad, will be discarded offline.

# this is the PF candidate isolation with pfnopu input and custom vetoes
from MyAnalysis.IsolationTools.electronPFIsolations_cff import *
from MyAnalysis.IsolationTools.muonPFIsolations_cff import *
from MyAnalysis.IsolationTools.electronDirectionalPFIsolations_cff import *
from MyAnalysis.IsolationTools.muonDirectionalPFIsolations_cff import *
from MyAnalysis.IsolationTools.electronPFPUIsolations_cff import *
from MyAnalysis.IsolationTools.muonPFPUIsolations_cff import *

from CommonTools.ParticleFlow.pfNoPileUp_cff import *
pfPileUp.PFCandidates = "particleFlow"
pfNoPileUp.bottomCollection = "particleFlow"

pfPUSequence = cms.Sequence( pfPileUp * pfNoPileUp )

pfIsolationCombined = cms.Sequence( electronCombinedPFIsoMapProducer * muonCombinedPFIsoMapProducer )

pfIsolationSingleType = cms.Sequence( electronPFIsoChHad01 * electronPFIsoNHad01 * electronPFIsoPhoton01 * muonPFIsoChHad01 * muonPFIsoNHad01 * muonPFIsoPhoton01 *
                                      electronPFIsoChHad02 * electronPFIsoNHad02 * electronPFIsoPhoton02 * muonPFIsoChHad02 * muonPFIsoNHad02 * muonPFIsoPhoton02 *
                                      electronPFIsoChHad03 * electronPFIsoNHad03 * electronPFIsoPhoton03 * muonPFIsoChHad03 * muonPFIsoNHad03 * muonPFIsoPhoton03 *
                                      electronPFIsoChHad04 * electronPFIsoNHad04 * electronPFIsoPhoton04 * muonPFIsoChHad04 * muonPFIsoNHad04 * muonPFIsoPhoton04 *
                                      electronPFIsoChHad05 * electronPFIsoNHad05 * electronPFIsoPhoton05 * muonPFIsoChHad05 * muonPFIsoNHad05 * muonPFIsoPhoton05 *
                                      electronPFIsoChHad06 * electronPFIsoNHad06 * electronPFIsoPhoton06 * muonPFIsoChHad06 * muonPFIsoNHad06 * muonPFIsoPhoton06 *
                                      electronPFIsoChHad07 * electronPFIsoNHad07 * electronPFIsoPhoton07 * muonPFIsoChHad07 * muonPFIsoNHad07 * muonPFIsoPhoton07 )

# just one cone size
pfDirectionalIsolationSingleType = cms.Sequence( electronDirPFIsoChHad04 * electronDirPFIsoNHad04 * electronDirPFIsoPhoton04 * muonDirPFIsoChHad04 * muonDirPFIsoNHad04 * muonDirPFIsoPhoton04 )

# just for charged, only 2 cones
pfPUIsolationSingleType = cms.Sequence( electronPUPFIsoChHad03 * electronPUPFIsoChHad04 * muonPUPFIsoChHad03 * muonPUPFIsoChHad04 )

pfIsolationAllSequence = cms.Sequence( pfPUSequence * pfIsolationCombined * pfIsolationSingleType * pfDirectionalIsolationSingleType * pfPUIsolationSingleType )

