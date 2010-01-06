import FWCore.ParameterSet.Config as cms

import RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi
lowThrCaloTowers = RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi.calotowermaker.clone()

lowThrCaloTowers.EBSumThreshold = 0.0
lowThrCaloTowers.EESumThreshold = 0.0
lowThrCaloTowers.HBThreshold = 0.3
lowThrCaloTowers.HF1Threshold = 0.6
lowThrCaloTowers.EBThreshold = 0.02
lowThrCaloTowers.HF2Threshold = 0.9
lowThrCaloTowers.EEThreshold = 0.2
lowThrCaloTowers.HESThreshold = 0.3
lowThrCaloTowers.HEDThreshold = 0.3
