import FWCore.ParameterSet.Config as cms

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
import RecoJets.JetAssociationProducers.ak5JTA_cff
newJetTracksAssociatorAtVertex = RecoJets.JetAssociationProducers.ak5JTA_cff.ak5JetTracksAssociatorAtVertex.clone()
newJetTracksAssociatorAtVertex.jets = "ak5CaloJetsL2L3"
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
newSimpleSecondaryVertexBJetTags = RecoBTag.Configuration.RecoBTag_cff.simpleSecondaryVertexBJetTags.clone()
newSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSecondaryVertexTagInfos") )
newCombinedSecondaryVertexBJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexBJetTags.clone()
newCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )
newCombinedSecondaryVertexMVABJetTags = RecoBTag.Configuration.RecoBTag_cff.combinedSecondaryVertexMVABJetTags.clone()
newCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )

# soft electron b-tag
newSoftElectronTagInfos = RecoBTag.Configuration.RecoBTag_cff.softElectronTagInfos.clone()
newSoftElectronTagInfos.jets = "ak5CaloJetsL2L3"
newSoftElectronBJetTags = RecoBTag.Configuration.RecoBTag_cff.softElectronBJetTags.clone()
newSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftElectronTagInfos") )

# soft muon b-tag
newSoftMuonTagInfos = RecoBTag.Configuration.RecoBTag_cff.softMuonTagInfos.clone()
newSoftMuonTagInfos.jets = "ak5CaloJetsL2L3"
newSoftMuonBJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonBJetTags.clone()
newSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )
newSoftMuonNoIPBJetTags = RecoBTag.Configuration.RecoBTag_cff.softMuonNoIPBJetTags.clone()
newSoftMuonNoIPBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )

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
       newSimpleSecondaryVertexBJetTags +
       newCombinedSecondaryVertexBJetTags +
       newCombinedSecondaryVertexMVABJetTags )
    )

newJetBtaggingEle = cms.Sequence(
    newSoftElectronTagInfos *
       newSoftElectronBJetTags
    )

newJetBtaggingMu = cms.Sequence(
    newSoftMuonTagInfos * (
       newSoftMuonBJetTags +
       newSoftMuonNoIPBJetTags )
    )

newJetBtagging = cms.Sequence(
    newJetBtaggingIP +
    newJetBtaggingSV +
    newJetBtaggingEle +
    newJetBtaggingMu )

newBtaggingSequence = cms.Sequence(
    newJetTracksAssociator *
       newJetBtagging )
