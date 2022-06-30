import FWCore.ParameterSet.Config as cms

abc = cms.EDProducer("ABCNetProducer",
                     candName = cms.InputTag("packedPFCandidates"),
                     graph_path = cms.FileInPath("CommonTools/PileupAlgos/data/AttentionBasedPileupRejectionModel_Run2.pb"),
                     preprocess_json = cms.FileInPath("CommonTools/PileupAlgos/data/preprocessing_info.json"),
                     input_tensor_name = cms.string("input_1"),
                     output_tensor_name = cms.string("Identity"),
                     n_pf_cands = cms.int32(4000),
                     n_feats = cms.int32(19),
                     debug = cms.bool(False)
)
