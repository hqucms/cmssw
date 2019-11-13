import FWCore.ParameterSet.Config as cms

from RecoBTag.FeatureTools.pfParticleNetTagInfos_cfi import pfParticleNetTagInfos
from RecoBTag.MXNet.pfParticleNetJetTags_cfi import pfParticleNetJetTags as _pfParticleNetJetTags
from RecoBTag.MXNet.Parameters.ParticleNet.V01.pfParticleNetPreprocessParams_cfi import pfParticleNetPreprocessParams
from RecoBTag.MXNet.pfParticleNetDiscriminatorsJetTags_cfi import pfParticleNetDiscriminatorsJetTags

pfParticleNetJetTags = _pfParticleNetJetTags.clone(
    preprocessParams = pfParticleNetPreprocessParams,
    model_path = 'RecoBTag/Combined/data/ParticleNet/V01/ParticleNet-symbol.json',
    param_path = 'RecoBTag/Combined/data/ParticleNet/V01/ParticleNet-0000.params',
    debugMode  = False, # debug
)

from CommonTools.PileupAlgos.Puppi_cff import puppi
from PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi import primaryVertexAssociation

# This task is not used, useful only if we run it from RECO jets (RECO/AOD)
pfParticleNetTask = cms.Task(puppi, primaryVertexAssociation,
                             pfParticleNetTagInfos, pfParticleNetJetTags, pfParticleNetDiscriminatorsJetTags)

# declare all the discriminators
_pfParticleNetJetTagsProbs = ['pfParticleNetJetTags:' + flav_name
                              for flav_name in pfParticleNetJetTags.flav_names]
# meta-taggers
_pfParticleNetJetTagsMetaDiscrs = ['pfParticleNetDiscriminatorsJetTags:' + disc.name.value()
                                   for disc in pfParticleNetDiscriminatorsJetTags.discriminators]

_pfParticleNetJetTagsAll = _pfParticleNetJetTagsProbs + _pfParticleNetJetTagsMetaDiscrs
