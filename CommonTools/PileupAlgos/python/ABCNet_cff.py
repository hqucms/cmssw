import FWCore.ParameterSet.Config as cms

abc = cms.EDProducer("ABCNetProducer",
                     candName = cms.InputTag("packedPFCandidates"),
                     graph_path = cms.FileInPath("CommonTools/PileupAlgos/plugins/AttentionBasedPileupRejectionModel_Run2.pb"),
                     preprocess_json = cms.FileInPath("CommonTools/PileupAlgos/plugins/preprocessing_info.json"),
                     input_tensor_name = cms.string("input_1"),
                     output_tensor_name = cms.string("Identity"),
                     debug = cms.bool(False)
)
