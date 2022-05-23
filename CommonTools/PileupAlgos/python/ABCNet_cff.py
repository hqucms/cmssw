import FWCore.ParameterSet.Config as cms

abc = cms.EDProducer("ABCNetProducer",
                     candName = cms.InputTag("packedPFCandidates"),
                     graph_path = cms.FileInPath("CommonTools/PileupAlgos/plugins/AttentionBasedPileupRejectionModel_Run2.pb")
)
