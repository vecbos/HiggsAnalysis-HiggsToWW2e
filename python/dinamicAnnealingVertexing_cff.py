# re-reconstructing the primary vertices with the Deterministic Annealing (DA) vertex finder
# from B. Mangano studies
import FWCore.ParameterSet.Config as cms

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi import *
import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi

offlinePrimaryVertices = RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi.offlinePrimaryVerticesDA.clone()
offlinePrimaryVertices.useBeamConstraint = cms.bool(False)
offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01)

vertexingSequence = cms.Sequence( offlinePrimaryVertices )
