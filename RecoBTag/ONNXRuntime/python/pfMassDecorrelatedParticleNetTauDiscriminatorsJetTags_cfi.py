import FWCore.ParameterSet.Config as cms

pfMassDecorrelatedParticleNetTauDiscriminatorsJetTags = cms.EDProducer(
   'BTagProbabilityToDiscriminator',
   discriminators = cms.VPSet(
      cms.PSet(
         name = cms.string('XbbvsQCD'),
         numerator = cms.VInputTag(
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probXbb'),
            ),
         denominator = cms.VInputTag(
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probXbb'),
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probQCD'),
            ),
         ),
      cms.PSet(
         name = cms.string('XccvsQCD'),
         numerator = cms.VInputTag(
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probXcc'),
            ),
         denominator = cms.VInputTag(
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probXcc'),
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probQCD'),
            ),
         ),
      cms.PSet(
         name = cms.string('XqqvsQCD'),
         numerator = cms.VInputTag(
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probXqq'),
            ),
         denominator = cms.VInputTag(
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probXqq'),
            cms.InputTag('pfMassDecorrelatedParticleNetTauJetTags', 'probQCD'),
            ),
         ),

      )
   )
