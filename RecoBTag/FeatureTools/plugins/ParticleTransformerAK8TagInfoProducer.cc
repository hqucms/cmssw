#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"

#include "Math/PtEtaPhiE4D.h"
#include "Math/LorentzVector.h"

using namespace btagbtvdeep;

class ParticleTransformerAK8TagInfoProducer : public edm::stream::EDProducer<> {
public:
  explicit ParticleTransformerAK8TagInfoProducer(const edm::ParameterSet &);
  ~ParticleTransformerAK8TagInfoProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  typedef std::vector<reco::DeepBoostedJetTagInfo> DeepBoostedJetTagInfoCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
  typedef reco::VertexCollection VertexCollection;
  typedef edm::View<reco::Candidate> CandidateView;

  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override {}

  void fillParticleFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet);
  void fillSVFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet);

  const double jet_radius_;
  const double min_jet_pt_;
  const double max_jet_eta_;
  const double min_pt_for_track_properties_;
  const double min_puppi_wgt_;
  const bool flip_ip_sign_;
  const double max_sip3dsig_;

  edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<SVCollection> sv_token_;
  edm::EDGetTokenT<CandidateView> pfcand_token_;
  edm::EDGetTokenT<CandidateView> lost_track_token_;

  bool use_puppi_value_map_;
  bool use_pvasq_value_map_;

  edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> pvasq_value_map_token_;
  edm::EDGetTokenT<edm::Association<VertexCollection>> pvas_token_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<SVCollection> svs_;
  edm::Handle<CandidateView> pfcands_;
  edm::Handle<CandidateView> lost_tracks_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;
  edm::Handle<edm::ValueMap<float>> puppi_value_map_;
  edm::Handle<edm::ValueMap<int>> pvasq_value_map_;
  edm::Handle<edm::Association<VertexCollection>> pvas_;

  const static std::vector<std::string> charged_particle_features_;
  const static std::vector<std::string> neutral_particle_features_;
  const static std::vector<std::string> sv_features_;
  const reco::Vertex *pv_ = nullptr;
};

const std::vector<std::string> ParticleTransformerAK8TagInfoProducer::charged_particle_features_{
    "cpfcandlt_VTX_ass",
    "cpfcandlt_lostInnerHits",
    "cpfcandlt_quality",
    "cpfcandlt_charge",
    "cpfcandlt_isEl",
    "cpfcandlt_isMu",
    "cpfcandlt_isChargedHad",
    "cpfcandlt_phirel",
    "cpfcandlt_etarel",
    "cpfcandlt_abseta",
    "cpfcandlt_normchi2",
    "cpfcandlt_dz",
    "cpfcandlt_dzsig",
    "cpfcandlt_dxy",
    "cpfcandlt_dxysig",
    "cpfcandlt_btagEtaRel",
    "cpfcandlt_btagPtRatio",
    "cpfcandlt_btagPParRatio",
    "cpfcandlt_btagSip3dVal",
    "cpfcandlt_btagSip3dSig",
    "cpfcandlt_btagJetDistVal",
    "cpfcandlt_mask",
    "cpfcandlt_pt_log_nopuppi",
    "cpfcandlt_e_log_nopuppi",
    "cpfcandlt_isLostTrack",
    "cpfcandlt_pixelBarrelLayersWithMeasurement",
    "cpfcandlt_pixelEndcapLayersWithMeasurement",
    "cpfcandlt_stripTECLayersWithMeasurement",
    "cpfcandlt_stripTIBLayersWithMeasurement",
    "cpfcandlt_stripTIDLayersWithMeasurement",
    "cpfcandlt_stripTOBLayersWithMeasurement",
    "cpfcandlt_px",
    "cpfcandlt_py",
    "cpfcandlt_pz",
    "cpfcandlt_energy"};

const std::vector<std::string> ParticleTransformerAK8TagInfoProducer::neutral_particle_features_{
    "npfcand_isGamma",
    "npfcand_isNeutralHad",
    "npfcand_phirel",
    "npfcand_etarel",
    "npfcand_abseta",
    "npfcand_mask",
    "npfcand_pt_log_nopuppi",
    "npfcand_e_log_nopuppi",
    "npfcand_px",
    "npfcand_py",
    "npfcand_pz",
    "npfcand_energy"};

const std::vector<std::string> ParticleTransformerAK8TagInfoProducer::sv_features_{"sv_mask",
                                                                                   "sv_phirel",
                                                                                   "sv_etarel",
                                                                                   "sv_abseta",
                                                                                   "sv_mass",
                                                                                   "sv_pt_log",
                                                                                   "sv_ntracks",
                                                                                   "sv_normchi2",
                                                                                   "sv_dxy",
                                                                                   "sv_dxysig",
                                                                                   "sv_d3d",
                                                                                   "sv_d3dsig",
                                                                                   "sv_px",
                                                                                   "sv_py",
                                                                                   "sv_pz",
                                                                                   "sv_energy"};

ParticleTransformerAK8TagInfoProducer::ParticleTransformerAK8TagInfoProducer(const edm::ParameterSet &iConfig)
    : jet_radius_(iConfig.getParameter<double>("jet_radius")),
      min_jet_pt_(iConfig.getParameter<double>("min_jet_pt")),
      max_jet_eta_(iConfig.getParameter<double>("max_jet_eta")),
      min_pt_for_track_properties_(iConfig.getParameter<double>("min_pt_for_track_properties")),
      min_puppi_wgt_(iConfig.getParameter<double>("min_puppi_wgt")),
      flip_ip_sign_(iConfig.getParameter<bool>("flip_ip_sign")),
      max_sip3dsig_(iConfig.getParameter<double>("sip3dSigMax")),
      jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
      pfcand_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pf_candidates"))),
      lost_track_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("lost_tracks"))),
      use_puppi_value_map_(false),
      use_pvasq_value_map_(false),
      track_builder_token_(
          esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))) {
  const auto &puppi_value_map_tag = iConfig.getParameter<edm::InputTag>("puppi_value_map");
  if (!puppi_value_map_tag.label().empty()) {
    puppi_value_map_token_ = consumes<edm::ValueMap<float>>(puppi_value_map_tag);
    use_puppi_value_map_ = true;
  }

  const auto &pvas_tag = iConfig.getParameter<edm::InputTag>("vertex_associator");
  if (!pvas_tag.label().empty()) {
    pvasq_value_map_token_ = consumes<edm::ValueMap<int>>(pvas_tag);
    pvas_token_ = consumes<edm::Association<VertexCollection>>(pvas_tag);
    use_pvasq_value_map_ = true;
  }

  produces<DeepBoostedJetTagInfoCollection>();
}

ParticleTransformerAK8TagInfoProducer::~ParticleTransformerAK8TagInfoProducer() {}

void ParticleTransformerAK8TagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // pfParticleTransformerAK8TagInfos
  edm::ParameterSetDescription desc;
  desc.add<double>("jet_radius", 0.8);
  desc.add<double>("min_jet_pt", 150);
  desc.add<double>("max_jet_eta", 99);
  desc.add<double>("min_pt_for_track_properties", -1);
  desc.add<double>("min_puppi_wgt", 0.01);
  desc.add<bool>("flip_ip_sign", false);
  desc.add<double>("sip3dSigMax", -1);
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
  desc.add<edm::InputTag>("pf_candidates", edm::InputTag("particleFlow"));
  desc.add<edm::InputTag>("lost_tracks", edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("jets", edm::InputTag("ak8PFJetsPuppi"));
  desc.add<edm::InputTag>("puppi_value_map", edm::InputTag("puppi"));
  desc.add<edm::InputTag>("vertex_associator", edm::InputTag("primaryVertexAssociation", "original"));
  descriptions.add("pfParticleTransformerAK8TagInfos", desc);
}

void ParticleTransformerAK8TagInfoProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  auto output_tag_infos = std::make_unique<DeepBoostedJetTagInfoCollection>();

  auto jets = iEvent.getHandle(jet_token_);

  iEvent.getByToken(vtx_token_, vtxs_);
  if (vtxs_->empty()) {
    // produce empty TagInfos in case no primary vertex
    iEvent.put(std::move(output_tag_infos));
    return;  // exit event
  }
  // primary vertex
  pv_ = &vtxs_->at(0);

  iEvent.getByToken(sv_token_, svs_);
  iEvent.getByToken(pfcand_token_, pfcands_);
  iEvent.getByToken(lost_track_token_, lost_tracks_);

  track_builder_ = iSetup.getHandle(track_builder_token_);

  if (use_puppi_value_map_) {
    iEvent.getByToken(puppi_value_map_token_, puppi_value_map_);
  }

  if (use_pvasq_value_map_) {
    iEvent.getByToken(pvasq_value_map_token_, pvasq_value_map_);
    iEvent.getByToken(pvas_token_, pvas_);
  }

  for (std::size_t jet_n = 0; jet_n < jets->size(); jet_n++) {
    const auto &jet = (*jets)[jet_n];
    edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);

    // create jet features
    DeepBoostedJetFeatures features;
    // declare all the feature variables (init as empty vector)
    for (const auto &name : charged_particle_features_) {
      features.add(name);
    }
    for (const auto &name : neutral_particle_features_) {
      features.add(name);
    }
    for (const auto &name : sv_features_) {
      features.add(name);
    }

    // fill values only if above pt threshold and has daughters, otherwise left empty
    bool fill_vars = true;
    if (jet.pt() < min_jet_pt_ || std::abs(jet.eta()) > max_jet_eta_)
      fill_vars = false;
    if (jet.numberOfDaughters() == 0)
      fill_vars = false;

    if (fill_vars) {
      fillParticleFeatures(features, jet);
      fillSVFeatures(features, jet);

      features.check_consistency(charged_particle_features_);
      features.check_consistency(neutral_particle_features_);
      features.check_consistency(sv_features_);
    }

    // this should always be done even if features are not filled
    output_tag_infos->emplace_back(features, jet_ref);
  }

  iEvent.put(std::move(output_tag_infos));
}

void ParticleTransformerAK8TagInfoProducer::fillParticleFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet) {
  // require the input to be a pat::Jet
  const auto *patJet = dynamic_cast<const pat::Jet *>(&jet);
  if (!patJet) {
    throw edm::Exception(edm::errors::InvalidReference) << "Input is not a pat::Jet.";
  }

  // do nothing if jet does not have constituents
  if (jet.numberOfDaughters() == 0)
    return;

  // some jet properties
  math::XYZVector jet_dir = jet.momentum().Unit();
  GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
  const float etasign = jet.eta() > 0 ? 1 : -1;

  TrackInfoBuilder trkinfo(track_builder_);

  auto puppiWgt = [&](const reco::CandidatePtr &cand) {
    const auto *pack_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));
    float wgt = 1.;
    if (pack_cand) {
      if (use_puppi_value_map_)
        wgt = (*puppi_value_map_)[cand];
    } else if (reco_cand) {
      if (use_puppi_value_map_) {
        wgt = (*puppi_value_map_)[cand];
      } else {
        throw edm::Exception(edm::errors::InvalidReference) << "Puppi value map is missing";
      }
    } else {
      throw edm::Exception(edm::errors::InvalidReference)
          << "Cannot convert to either pat::PackedCandidate or reco::PFCandidate";
    }
    return wgt;
  };

  std::vector<reco::CandidatePtr> cpfPtrs, npfPtrs;
  std::map<reco::CandidatePtr::key_type, bool> isLostTrackMap;

  for (const auto &dau : jet.daughterPtrVector()) {
    // get the original reco/packed candidate not scaled by the puppi weight
    auto cand = pfcands_->ptrAt(dau.key());
    // remove particles w/ extremely low puppi weights
    if ((puppiWgt(cand)) < min_puppi_wgt_)
      continue;
    // only when computing the nagative tagger: remove charged candidates with high sip3d
    if (flip_ip_sign_ && cand->charge()) {
      trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, *pv_);
      if (trkinfo.getTrackSip3dSig() > max_sip3dsig_)
        continue;
    }
    if (cand->charge() != 0) {
      cpfPtrs.push_back(cand);
      isLostTrackMap[cand.key()] = false;
    } else {
      npfPtrs.push_back(cand);
    }
  }
  // lost tracks: fill to cpfcands
  for (size_t i = 0; i < lost_tracks_->size(); ++i) {
    auto cand = lost_tracks_->ptrAt(i);
    if (reco::deltaR(*cand, jet) < jet_radius_) {
      cpfPtrs.push_back(cand);
      isLostTrackMap[cand.key()] = true;
    }
  }
  // sort by pt
  std::sort(cpfPtrs.begin(), cpfPtrs.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });
  std::sort(npfPtrs.begin(), npfPtrs.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });

  // reserve space
  for (const auto &name : charged_particle_features_) {
    fts.reserve(name, cpfPtrs.size());
  }

  auto useTrackProperties = [&](const reco::PFCandidate *reco_cand) {
    const auto *trk = reco_cand->bestTrack();
    return trk != nullptr && trk->pt() > min_pt_for_track_properties_;
  };

  for (const auto &cand : cpfPtrs) {
    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));

    const float ip_sign = flip_ip_sign_ ? -1 : 1;

    auto candP4 = cand->p4();
    if (packed_cand) {
      fts.fill("cpfcandlt_VTX_ass", packed_cand->pvAssociationQuality());
      fts.fill("cpfcandlt_lostInnerHits", packed_cand->lostInnerHits());
      fts.fill("cpfcandlt_quality", packed_cand->bestTrack() ? packed_cand->bestTrack()->qualityMask() : 0);

      fts.fill("cpfcandlt_charge", packed_cand->charge());
      fts.fill("cpfcandlt_isEl", std::abs(packed_cand->pdgId()) == 11);
      fts.fill("cpfcandlt_isMu", std::abs(packed_cand->pdgId()) == 13);
      fts.fill("cpfcandlt_isChargedHad", std::abs(packed_cand->pdgId()) == 211);
      fts.fill("cpfcandlt_isLostTrack", isLostTrackMap[cand.key()]);

      // impact parameters
      fts.fill("cpfcandlt_dz", ip_sign * packed_cand->dz());
      fts.fill("cpfcandlt_dxy", ip_sign * packed_cand->dxy());
      fts.fill("cpfcandlt_dzsig", packed_cand->bestTrack() ? ip_sign * packed_cand->dz() / packed_cand->dzError() : 0);
      fts.fill("cpfcandlt_dxysig",
               packed_cand->bestTrack() ? ip_sign * packed_cand->dxy() / packed_cand->dxyError() : 0);

    } else if (reco_cand) {
      // get vertex association quality
      int pv_ass_quality = 0;  // fallback value
      float vtx_ass = 0;
      if (use_pvasq_value_map_) {
        pv_ass_quality = (*pvasq_value_map_)[cand];
        const reco::VertexRef &PV_orig = (*pvas_)[cand];
        vtx_ass = vtx_ass_from_pfcand(*reco_cand, pv_ass_quality, PV_orig);
      } else {
        throw edm::Exception(edm::errors::InvalidReference) << "Vertex association missing";
      }

      fts.fill("cpfcandlt_VTX_ass", vtx_ass);
      fts.fill("cpfcandlt_lostInnerHits", useTrackProperties(reco_cand) ? lost_inner_hits_from_pfcand(*reco_cand) : 0);
      fts.fill("cpfcandlt_quality", useTrackProperties(reco_cand) ? quality_from_pfcand(*reco_cand) : 0);

      fts.fill("cpfcandlt_charge", reco_cand->charge());
      fts.fill("cpfcandlt_isEl", std::abs(reco_cand->pdgId()) == 11);
      fts.fill("cpfcandlt_isMu", std::abs(reco_cand->pdgId()) == 13);
      fts.fill("cpfcandlt_isChargedHad", std::abs(reco_cand->pdgId()) == 211);
      fts.fill("cpfcandlt_isLostTrack", isLostTrackMap[cand.key()]);

      // impact parameters
      const auto *trk = reco_cand->bestTrack();
      float dz = trk ? ip_sign * trk->dz(pv_->position()) : 0;
      float dxy = trk ? ip_sign * trk->dxy(pv_->position()) : 0;
      fts.fill("cpfcandlt_dz", dz);
      fts.fill("cpfcandlt_dzsig", trk ? dz / trk->dzError() : 0);
      fts.fill("cpfcandlt_dxy", dxy);
      fts.fill("cpfcandlt_dxysig", trk ? dxy / trk->dxyError() : 0);
    }

    // basic kinematics
    fts.fill("cpfcandlt_px", candP4.px());
    fts.fill("cpfcandlt_py", candP4.py());
    fts.fill("cpfcandlt_pz", candP4.pz());
    fts.fill("cpfcandlt_energy", candP4.energy());

    fts.fill("cpfcandlt_phirel", reco::deltaPhi(candP4, jet));
    fts.fill("cpfcandlt_etarel", etasign * (candP4.eta() - jet.eta()));
    fts.fill("cpfcandlt_abseta", std::abs(candP4.eta()));

    fts.fill("cpfcandlt_mask", 1);
    fts.fill("cpfcandlt_pt_log_nopuppi", std::log(cand->pt()));
    fts.fill("cpfcandlt_e_log_nopuppi", std::log(cand->energy()));

    const reco::Track *trk = nullptr;
    if (packed_cand) {
      trk = packed_cand->bestTrack();
    } else if (reco_cand && useTrackProperties(reco_cand)) {
      trk = reco_cand->bestTrack();
    }
    if (trk) {
      fts.fill("cpfcandlt_normchi2", std::floor(trk->normalizedChi2()));

      trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, *pv_);
      fts.fill("cpfcandlt_btagEtaRel", trkinfo.getTrackEtaRel());
      fts.fill("cpfcandlt_btagPtRatio", trkinfo.getTrackPtRatio());
      fts.fill("cpfcandlt_btagPParRatio", trkinfo.getTrackPParRatio());
      fts.fill("cpfcandlt_btagSip3dVal", ip_sign * trkinfo.getTrackSip3dVal());
      fts.fill("cpfcandlt_btagSip3dSig", ip_sign * trkinfo.getTrackSip3dSig());
      fts.fill("cpfcandlt_btagJetDistVal", trkinfo.getTrackJetDistVal());

      fts.fill("cpfcandlt_pixelBarrelLayersWithMeasurement", trk->hitPattern().pixelBarrelLayersWithMeasurement());
      fts.fill("cpfcandlt_pixelEndcapLayersWithMeasurement", trk->hitPattern().pixelEndcapLayersWithMeasurement());
      fts.fill("cpfcandlt_stripTIBLayersWithMeasurement", trk->hitPattern().stripTIBLayersWithMeasurement());
      fts.fill("cpfcandlt_stripTIDLayersWithMeasurement", trk->hitPattern().stripTIDLayersWithMeasurement());
      fts.fill("cpfcandlt_stripTOBLayersWithMeasurement", trk->hitPattern().stripTOBLayersWithMeasurement());
      fts.fill("cpfcandlt_stripTECLayersWithMeasurement", trk->hitPattern().stripTECLayersWithMeasurement());

    } else {
      fts.fill("cpfcandlt_normchi2", 999);

      fts.fill("cpfcandlt_btagEtaRel", 0);
      fts.fill("cpfcandlt_btagPtRatio", 0);
      fts.fill("cpfcandlt_btagPParRatio", 0);
      fts.fill("cpfcandlt_btagSip3dVal", 0);
      fts.fill("cpfcandlt_btagSip3dSig", 0);
      fts.fill("cpfcandlt_btagJetDistVal", 0);

      fts.fill("cpfcandlt_pixelBarrelLayersWithMeasurement", 0);
      fts.fill("cpfcandlt_pixelEndcapLayersWithMeasurement", 0);
      fts.fill("cpfcandlt_stripTIBLayersWithMeasurement", 0);
      fts.fill("cpfcandlt_stripTIDLayersWithMeasurement", 0);
      fts.fill("cpfcandlt_stripTOBLayersWithMeasurement", 0);
      fts.fill("cpfcandlt_stripTECLayersWithMeasurement", 0);
    }
  }

  // fill neutral candidate features
  for (const auto &cand : npfPtrs) {
    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));

    auto candP4 = cand->p4();
    if (packed_cand) {
      fts.fill("npfcand_isGamma", std::abs(packed_cand->pdgId()) == 22);
      fts.fill("npfcand_isNeutralHad", std::abs(packed_cand->pdgId()) == 130);
    } else if (reco_cand) {
      fts.fill("npfcand_isGamma", std::abs(reco_cand->pdgId()) == 22);
      fts.fill("npfcand_isNeutralHad", std::abs(reco_cand->pdgId()) == 130);
    }

    // basic kinematics
    fts.fill("npfcand_px", candP4.px());
    fts.fill("npfcand_py", candP4.py());
    fts.fill("npfcand_pz", candP4.pz());
    fts.fill("npfcand_energy", candP4.energy());

    fts.fill("npfcand_phirel", reco::deltaPhi(candP4, jet));
    fts.fill("npfcand_etarel", etasign * (candP4.eta() - jet.eta()));
    fts.fill("npfcand_abseta", std::abs(candP4.eta()));

    fts.fill("npfcand_mask", 1);
    fts.fill("npfcand_pt_log_nopuppi", std::log(cand->pt()));
    fts.fill("npfcand_e_log_nopuppi", std::log(cand->energy()));
  }
}

void ParticleTransformerAK8TagInfoProducer::fillSVFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet) {
  std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
  for (const auto &sv : *svs_) {
    if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
      jetSVs.push_back(&sv);
    }
  }
  // sort by dxy significance
  std::sort(jetSVs.begin(),
            jetSVs.end(),
            [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) {
              return sv_vertex_comparator(*sva, *svb, *pv_);
            });

  // reserve space
  for (const auto &name : sv_features_) {
    fts.reserve(name, jetSVs.size());
  }

  const float etasign = jet.eta() > 0 ? 1 : -1;

  for (const auto *sv : jetSVs) {
    // basic kinematics
    fts.fill("sv_mask", 1);

    fts.fill("sv_px", sv->px());
    fts.fill("sv_py", sv->py());
    fts.fill("sv_pz", sv->pz());
    fts.fill("sv_energy", sv->energy());

    fts.fill("sv_phirel", reco::deltaPhi(*sv, jet));
    fts.fill("sv_etarel", etasign * (sv->eta() - jet.eta()));
    fts.fill("sv_abseta", std::abs(sv->eta()));
    fts.fill("sv_mass", sv->mass());
    fts.fill("sv_pt_log", std::log(sv->pt()));

    // sv properties
    fts.fill("sv_ntracks", sv->numberOfDaughters());
    fts.fill("sv_normchi2", sv->vertexNormalizedChi2());

    const auto &dxy = vertexDxy(*sv, *pv_);
    fts.fill("sv_dxy", dxy.value());
    fts.fill("sv_dxysig", dxy.significance());

    const auto &d3d = vertexD3d(*sv, *pv_);
    fts.fill("sv_d3d", d3d.value());
    fts.fill("sv_d3dsig", d3d.significance());
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(ParticleTransformerAK8TagInfoProducer);
