import FWCore.ParameterSet.Config as cms

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
import RecoJets.JetAssociationProducers.ak5JTA_cff
newJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
newJetTracksAssociatorAtVertex.jets = "ak5CaloJetsL2L3Residual"
newJetTracksAssociatorAtVertex.tracks = "generalTracks"

# impact parameter b-tag
import RecoBTag.Configuration.RecoBTag_cff
newImpactParameterTagInfos = RecoBTag.Configuration.RecoBTag_cff.impactParameterTagInfos.clone()
newImpactParameterTagInfos.jetTracks = "newJetTracksAssociatorAtVertex"
newTrackCountingHighEffBJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighEffBJetTags.clone()
newTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
newTrackCountingHighPurBJetTags = RecoBTag.Configuration.RecoBTag_cff.trackCountingHighPurBJetTags.clone()
newTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
newJetProbabilityBJetTags = RecoBTag.Configuration.RecoBTag_cff.jetProbabilityBJetTags.clone()
newJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
newJetBProbabilityBJetTags = RecoBTag.Configuration.RecoBTag_cff.jetBProbabilityBJetTags.clone()
newJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
newImpactParameterMVABJetTags = RecoBTag.Configuration.RecoBTag_cff.impactParameterMVABJetTags.clone()
newImpactParameterMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )

# secondary vertex b-tag
newSecondaryVertexTagInfos = RecoBTag.Configuration.RecoBTag_cff.secondaryVertexTagInfos.clone()
newSecondaryVertexTagInfos.trackIPTagInfos = "newImpactParameterTagInfos"
newSimpleSecondaryVertexHighEffBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighEffBJetTags.clone()
newSimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSecondaryVertexTagInfos") )
newSimpleSecondaryVertexHighPurBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexHighPurBJetTags.clone()
newSimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSecondaryVertexTagInfos") )
newCombinedSecondaryVertexBJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
newCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )
newCombinedSecondaryVertexMVABJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
newCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )

# soft electron b-tag
newSoftElectronTagInfos = RecoBTag.Configuration.RecoBTag_cff.softElectronTagInfos.clone()
newSoftElectronTagInfos.jets = "ak5CaloJetsL2L3Residual"
newSoftElectronBJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronBJetTags.clone()
newSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftElectronTagInfos") )
newSoftElectronByIP3dBJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByIP3dBJetTags.clone()
newSoftElectronByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftElectronTagInfos") )
newSoftElectronByPtBJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronByPtBJetTags.clone()
newSoftElectronByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftElectronTagInfos") )

# soft muon b-tag
newSoftMuonTagInfos = RecoBTag.Configuration.RecoBTag_cff.softMuonTagInfos.clone()
newSoftMuonTagInfos.jets = "ak5CaloJetsL2L3Residual"
newSoftMuonBJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonBJetTags.clone()
newSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )
newSoftMuonByIP3dBJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByIP3dBJetTags.clone()
newSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )
newSoftMuonByPtBJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonByPtBJetTags.clone()
newSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )

# prepare a path running the new modules
newJetTracksAssociator = cms.Sequence(
    newJetTracksAssociatorAtVertex
    )

newJetBtaggingIP = cms.Sequence(
    newImpactParameterTagInfos * (
       newTrackCountingHighEffBJetTags +
       newTrackCountingHighPurBJetTags +
       newJetProbabilityBJetTags +
       newJetBProbabilityBJetTags )
    )

newJetBtaggingSV = cms.Sequence(
    newImpactParameterTagInfos *
    newSecondaryVertexTagInfos * (
       newSimpleSecondaryVertexHighEffBJetTags +
       newSimpleSecondaryVertexHighPurBJetTags +
       newCombinedSecondaryVertexBJetTags +
       newCombinedSecondaryVertexMVABJetTags )
    )

newJetBtaggingEle = cms.Sequence(
    softElectronCands * 
    newSoftElectronTagInfos *
       newSoftElectronBJetTags +
       newSoftElectronByIP3dBJetTags +
       newSoftElectronByPtBJetTags
    )

newJetBtaggingMu = cms.Sequence(
    newSoftMuonTagInfos * (
       newSoftMuonBJetTags +
       newSoftMuonByIP3dBJetTags +
       newSoftMuonByPtBJetTags )
    )

newJetBtagging = cms.Sequence(
    newJetBtaggingIP +
    newJetBtaggingSV +
    newJetBtaggingEle +
    newJetBtaggingMu )

newBtaggingSequence = cms.Sequence(
    newJetTracksAssociator *
       newJetBtagging )
