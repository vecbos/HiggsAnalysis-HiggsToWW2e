import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff import *
import RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff

alternativeIsolationSequence = cms.Sequence( egammaIsolationSequence )
