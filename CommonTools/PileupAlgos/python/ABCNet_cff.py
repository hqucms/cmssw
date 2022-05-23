import FWCore.ParameterSet.Config as cms

abc = cms.EDProducer("ABCNetProducer",
                     candName = cms.InputTag("packedPFCandidates"),
                     graph_path = cms.FileInPath("toy/path/for/now/model.pb")
)
