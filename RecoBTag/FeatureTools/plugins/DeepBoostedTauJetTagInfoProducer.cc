#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/BTauReco/interface/DeepBoostedTauJetTagInfo.h"
#include <TVector3.h>

using namespace btagbtvdeep;

class DeepBoostedTauJetTagInfoProducer : public edm::stream::EDProducer<> {
public:
  explicit DeepBoostedTauJetTagInfoProducer(const edm::ParameterSet &);
  ~DeepBoostedTauJetTagInfoProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  typedef std::vector<reco::DeepBoostedTauJetTagInfo> DeepBoostedTauJetTagInfoCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
  typedef reco::VertexCollection VertexCollection;
  typedef edm::View<reco::Candidate> CandidateView;

  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override {}

  void fillParticleFeatures(DeepBoostedTauJetFeatures &fts, const reco::Jet &jet);
  void fillSVFeatures(DeepBoostedTauJetFeatures &fts, const reco::Jet &jet);

  const double jet_radius_;
  const double min_jet_pt_;
  const double max_jet_eta_;
  const double min_pt_for_track_properties_;
  const bool use_puppiP4_;
  const bool include_neutrals_;
  const bool sort_by_sip2dsig_;
  const double min_puppi_wgt_;
  const bool flip_ip_sign_;
  const double max_sip3dsig_;

  edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<SVCollection> sv_token_;
  edm::EDGetTokenT<CandidateView> pfcand_token_;

  bool use_puppi_value_map_;
  bool use_pvasq_value_map_;

  edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> pvasq_value_map_token_;
  edm::EDGetTokenT<edm::Association<VertexCollection>> pvas_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<SVCollection> svs_;
  edm::Handle<CandidateView> pfcands_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;
  edm::Handle<edm::ValueMap<float>> puppi_value_map_;
  edm::Handle<edm::ValueMap<int>> pvasq_value_map_;
  edm::Handle<edm::Association<VertexCollection>> pvas_;

  const static std::vector<std::string> particle_features_;
  const static std::vector<std::string> sv_features_;
  const reco::Vertex *pv_ = nullptr;
};

const std::vector<std::string> DeepBoostedTauJetTagInfoProducer::particle_features_{
  "pfcand_mask",
    "jet_pfcand_pt_log", "jet_pfcand_energy_log", "jet_pfcand_deta", "jet_pfcand_dphi", "jet_pfcand_eta",
    "jet_pfcand_charge", "jet_pfcand_frompv", "jet_pfcand_nlostinnerhits",
    "jet_pfcand_track_chi2",
    "jet_pfcand_track_qual",
    "jet_pfcand_dz",
    "jet_pfcand_dzsig",
    "jet_pfcand_dxy",
    "jet_pfcand_dxysig",
    "jet_pfcand_etarel",
    "jet_pfcand_pperp_ratio",
    "jet_pfcand_ppara_ratio",
    "jet_pfcand_trackjet_d3d",
    "jet_pfcand_trackjet_d3dsig",
    "jet_pfcand_trackjet_dist",
    "jet_pfcand_trackjet_decayL",
    "jet_pfcand_npixhits",
    "jet_pfcand_nstriphits",
    "jet_pfcand_id",
    "jet_pfcand_hcalfraction",
    "jet_pfcand_calofraction",
};


const std::vector<std::string> DeepBoostedTauJetTagInfoProducer::sv_features_{
  "sv_mask",
    "jet_sv_pt_log",
    "jet_sv_mass",
    "jet_sv_deta",
    "jet_sv_dphi",
    "jet_sv_eta",
    "jet_sv_ntrack",
    "jet_sv_chi2",
    "jet_sv_dxy",
    "jet_sv_dxysig",
    "jet_sv_d3d",
    "jet_sv_d3dsig"
};



DeepBoostedTauJetTagInfoProducer::DeepBoostedTauJetTagInfoProducer(const edm::ParameterSet &iConfig)
    : jet_radius_(iConfig.getParameter<double>("jet_radius")),
      min_jet_pt_(iConfig.getParameter<double>("min_jet_pt")),
      max_jet_eta_(iConfig.getParameter<double>("max_jet_eta")),
      min_pt_for_track_properties_(iConfig.getParameter<double>("min_pt_for_track_properties")),
      use_puppiP4_(iConfig.getParameter<bool>("use_puppiP4")),
      include_neutrals_(iConfig.getParameter<bool>("include_neutrals")),
      sort_by_sip2dsig_(iConfig.getParameter<bool>("sort_by_sip2dsig")),
      min_puppi_wgt_(iConfig.getParameter<double>("min_puppi_wgt")),
      flip_ip_sign_(iConfig.getParameter<bool>("flip_ip_sign")),
      max_sip3dsig_(iConfig.getParameter<double>("sip3dSigMax")),
      jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
      pfcand_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pf_candidates"))),
      use_puppi_value_map_(false),
      use_pvasq_value_map_(false) {
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

  produces<DeepBoostedTauJetTagInfoCollection>();
}

DeepBoostedTauJetTagInfoProducer::~DeepBoostedTauJetTagInfoProducer() {}

void DeepBoostedTauJetTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // pfDeepBoostedTauJetTagInfos
  edm::ParameterSetDescription desc;
  desc.add<double>("jet_radius", 0.8);
  desc.add<double>("min_jet_pt", 150);
  desc.add<double>("max_jet_eta", 99);
  desc.add<double>("min_pt_for_track_properties", -1);
  desc.add<bool>("use_puppiP4", true);
  desc.add<bool>("include_neutrals", true);
  desc.add<bool>("sort_by_sip2dsig", false);
  desc.add<double>("min_puppi_wgt", 0.01);
  desc.add<bool>("flip_ip_sign", false);
  desc.add<double>("sip3dSigMax", -1);
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
  desc.add<edm::InputTag>("pf_candidates", edm::InputTag("particleFlow"));
  desc.add<edm::InputTag>("jets", edm::InputTag("ak8PFJetsPuppi"));
  desc.add<edm::InputTag>("puppi_value_map", edm::InputTag("puppi"));
  desc.add<edm::InputTag>("vertex_associator", edm::InputTag("primaryVertexAssociation", "original"));
  descriptions.add("pfDeepBoostedTauJetTagInfos", desc);
}

void DeepBoostedTauJetTagInfoProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  auto output_tag_infos = std::make_unique<DeepBoostedTauJetTagInfoCollection>();

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

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);

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
    DeepBoostedTauJetFeatures features;
    // declare all the feature variables (init as empty vector)
    for (const auto &name : particle_features_) {
      features.add(name);
    }
    for (const auto &name : sv_features_) {
      features.add(name);
    }

    // fill values only if above pt threshold and has daughters, otherwise left
    // empty
    //bool fill_vars = true;
    //if (jet.pt() < min_jet_pt_ || std::abs(jet.eta()) > max_jet_eta_)
      //fill_vars = false;
    //if (jet.numberOfDaughters() == 0)
      //fill_vars = false;


    if (!( (jet.pt() < min_jet_pt_) || (std::abs(jet.eta()) > max_jet_eta_) || (jet.numberOfDaughters() == 0) )) {
      std::cout << " features \n";
      fillParticleFeatures(features, jet);
      fillSVFeatures(features, jet);
      
      std::cout << "check_consistency\n";
      features.check_consistency(particle_features_);
      features.check_consistency(sv_features_);
    }
    


    // this should always be done even if features are not filled
    output_tag_infos->emplace_back(features, jet_ref);
  }

  iEvent.put(std::move(output_tag_infos));
}

void DeepBoostedTauJetTagInfoProducer::fillParticleFeatures(DeepBoostedTauJetFeatures &fts, const reco::Jet &jet) {
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
  TVector3 jet_direction(jet.momentum().Unit().x(), jet.momentum().Unit().y(), jet.momentum().Unit().z());
  GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
  const float etasign = jet.eta() > 0 ? 1 : -1;

  std::map<reco::CandidatePtr::key_type, float> puppi_wgt_cache;
  auto puppiWgt = [&](const reco::CandidatePtr &cand) {
    const auto *pack_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));
    float wgt = 1.;
    if (pack_cand) {
      wgt = pack_cand->puppiWeight();
    } else if (reco_cand) {
      if (use_puppi_value_map_) {
        wgt = (*puppi_value_map_)[cand];
      } else {
        throw edm::Exception(edm::errors::InvalidReference) << "Puppi value map is missing";
      }
    } else {
      throw edm::Exception(edm::errors::InvalidReference) << "Cannot convert to either pat::PackedCandidate or "
                                                             "reco::PFCandidate";
    }
    puppi_wgt_cache[cand.key()] = wgt;
    return wgt;
  };

  std::vector<reco::CandidatePtr> daughters;
  for (const auto &dau : jet.daughterPtrVector()) {
    // remove particles w/ extremely low puppi weights
    // [Note] use jet daughters here to get the puppiWgt correctly
    if ((puppiWgt(dau)) < min_puppi_wgt_)
      continue;
    // from here: get the original reco/packed candidate not scaled by the puppi weight
    auto cand = pfcands_->ptrAt(dau.key());
    // base requirements on PF candidates
    //if (use_pnettau_features_ and cand->pt() < min_pt_for_pfcandidates_)
    //  continue;
    // charged candidate selection (for Higgs Interaction Net)
    if (!include_neutrals_ && (cand->charge() == 0 || cand->pt() < min_pt_for_track_properties_))
      continue;
    // only when computing the nagative tagger: remove charged candidates with high sip3d
    if (flip_ip_sign_ && cand->charge()) {
      TrackInfoBuilder trkinfo(track_builder_);
      trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, *pv_);
      if (trkinfo.getTrackSip3dSig() > max_sip3dsig_)
        continue;
    }
    daughters.push_back(cand);
  }

  std::vector<btagbtvdeep::SortingClass<reco::CandidatePtr>> c_sorted;
  if (sort_by_sip2dsig_) {
    // sort charged pf candidates by 2d impact parameter significance
    for (const auto &cand : daughters) {
      TrackInfoBuilder trkinfo(track_builder_);
      trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, *pv_);
      c_sorted.emplace_back(cand,
                            trkinfo.getTrackSip2dSig(),
                            -btagbtvdeep::mindrsvpfcand(*svs_, &(*cand), jet_radius_),
                            cand->pt() / jet.pt());
      std::sort(c_sorted.begin(), c_sorted.end(), btagbtvdeep::SortingClass<reco::CandidatePtr>::compareByABCInv);
    }
    for (unsigned int i = 0; i < c_sorted.size(); i++) {
      const auto &c = c_sorted.at(i);
      const auto &cand = c.get();
      daughters.at(i) = cand;
    }
  } else {
    if (use_puppiP4_) {
      // sort by Puppi-weighted pt
      std::sort(daughters.begin(),
                daughters.end(),
                [&puppi_wgt_cache](const reco::CandidatePtr &a, const reco::CandidatePtr &b) {
                  return puppi_wgt_cache.at(a.key()) * a->pt() > puppi_wgt_cache.at(b.key()) * b->pt();
                });
    } else {
      // sort by original pt (not Puppi-weighted)
      std::sort(daughters.begin(), daughters.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });
    }
  }

  // reserve space
  for (const auto &name : particle_features_)
    fts.reserve(name, daughters.size());


  auto useTrackProperties = [&](const reco::PFCandidate *reco_cand) {
    const auto *trk = reco_cand->bestTrack();
    return trk != nullptr && trk->pt() > min_pt_for_track_properties_;
  };

  for (const auto &cand : daughters) {
    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));

    if (!include_neutrals_ &&
        ((packed_cand && !packed_cand->hasTrackDetails()) || (reco_cand && !useTrackProperties(reco_cand))))
      continue;

    const float ip_sign = flip_ip_sign_ ? -1 : 1;

    auto candP4 = use_puppiP4_ ? puppi_wgt_cache.at(cand.key()) * cand->p4() : cand->p4();
    auto candP3 = use_puppiP4_ ? puppi_wgt_cache.at(cand.key()) * cand->momentum() : cand->momentum();

    if (packed_cand) {
      
      reco::VertexRef pv_ass = reco::VertexRef(vtxs_, 0);
      math::XYZPoint pv_ass_pos = pv_ass->position();

      TVector3 cand_direction(candP3.x(), candP3.y(), candP3.z());

      fts.fill("jet_pfcand_pt_log", std::log(candP4.pt()));
      fts.fill("jet_pfcand_energy_log", std::log(candP4.energy()));
      fts.fill("jet_pfcand_eta", candP4.eta());
      fts.fill("jet_pfcand_deta", jet_direction.Eta() - cand_direction.Eta());
      fts.fill("jet_pfcand_dphi", jet_direction.DeltaPhi(cand_direction));
      fts.fill("jet_pfcand_charge", packed_cand->charge());
      fts.fill("jet_pfcand_etarel", reco::btau::etaRel(jet_dir, candP3));
      fts.fill("jet_pfcand_pperp_ratio", jet_direction.Perp(cand_direction) / cand_direction.Mag());
      fts.fill("jet_pfcand_ppara_ratio", jet_direction.Dot(cand_direction) / cand_direction.Mag());
      fts.fill("jet_pfcand_frompv", packed_cand->fromPV());
      fts.fill("jet_pfcand_dz", packed_cand->dz(pv_ass_pos));
      fts.fill("jet_pfcand_dxy", packed_cand->dxy(pv_ass_pos));
      fts.fill("jet_pfcand_nlostinnerhits", packed_cand->lostInnerHits());
      fts.fill("jet_pfcand_npixhits", packed_cand->numberOfPixelHits());
      fts.fill("jet_pfcand_nstriphits", packed_cand->stripLayersWithMeasurement());
      
      if(abs(packed_cand->pdgId()) == 11 and packed_cand->charge() != 0)
	fts.fill("jet_pfcand_id",0);
      else if(abs(packed_cand->pdgId()) == 13 and packed_cand->charge() != 0)
	fts.fill("jet_pfcand_id",1);
      else if(abs(packed_cand->pdgId()) == 22 and packed_cand->charge() == 0)
	fts.fill("jet_pfcand_id",2);
      else if(abs(packed_cand->pdgId()) != 22 and packed_cand->charge() == 0)
	fts.fill("jet_pfcand_id",3);
      else if(abs(packed_cand->pdgId()) != 11 and abs(packed_cand->pdgId()) != 13 and packed_cand->charge() != 0)
	fts.fill("jet_pfcand_id",4);
      else
	fts.fill("jet_pfcand_id",-1);
      
      fts.fill("jet_pfcand_hcalfraction",packed_cand->hcalFraction());
      fts.fill("jet_pfcand_calofraction",packed_cand->caloFraction());
      fts.fill("pfcand_mask", 1);
      
      //fts.fill("jet_pfcand_dzsig", fabs(packed_cand->dz(pv_ass_pos)) / packed_cand->dzError());
      //fts.fill("jet_pfcand_dxysig", fabs(packed_cand->dxy(pv_ass_pos)) / packed_cand->dxyError());
      
      //fts.fill("jet_pfcand_dzsig", packed_cand->bestTrack() ? ip_sign * packed_cand->dz() / packed_cand->dzError() : 0);
      //fts.fill("jet_pfcand_dxysig", packed_cand->bestTrack() ? ip_sign * packed_cand->dxy() / packed_cand->dxyError() : 0);

      fts.fill("jet_pfcand_dzsig", packed_cand->bestTrack() ? fabs(packed_cand->dz(pv_ass_pos)) / packed_cand->dzError() : 0);
      fts.fill("jet_pfcand_dxysig", packed_cand->bestTrack() ? fabs(packed_cand->dxy(pv_ass_pos)) / packed_cand->dxyError() : 0);
    }

    const reco::Track *trk = nullptr;
    if (packed_cand) {
      trk = packed_cand->bestTrack();
    } else if (reco_cand && useTrackProperties(reco_cand)) {
      trk = reco_cand->bestTrack();
    }

    
    if (trk) {
      fts.fill("jet_pfcand_track_chi2", trk->normalizedChi2());
      fts.fill("jet_pfcand_track_qual", trk->qualityMask());
      /*
      auto cov = [&](unsigned i, unsigned j) { return trk->covariance(i, j); };
      fts.fill("jet_pfcand_dptdpt", cov(0, 0));
      fts.fill("jet_pfcand_detadeta", cov(1, 1));
      fts.fill("jet_pfcand_dphidphi", cov(2, 2));
      fts.fill("jet_pfcand_dxydxy", cov(3, 3));
      fts.fill("jet_pfcand_dzdz", cov(4, 4));
      fts.fill("jet_pfcand_dxydz", cov(3, 4));
      fts.fill("jet_pfcand_dphidxy", cov(2, 3));
      fts.fill("jet_pfcand_dlambdadz", cov(1, 4));
      */
      reco::TransientTrack transientTrack = track_builder_->build(*trk);
      Measurement1D meas_ip2d    = IPTools::signedTransverseImpactParameter(transientTrack, jet_ref_track_dir, *pv_).second;
      Measurement1D meas_ip3d    = IPTools::signedImpactParameter3D(transientTrack, jet_ref_track_dir, *pv_).second;
      Measurement1D meas_jetdist = IPTools::jetTrackDistance(transientTrack, jet_ref_track_dir, *pv_).second;
      Measurement1D meas_decayl  = IPTools::signedDecayLength3D(transientTrack, jet_ref_track_dir, *pv_).second;

      fts.fill("jet_pfcand_trackjet_d3d", meas_ip3d.value());
      fts.fill("jet_pfcand_trackjet_d3dsig", fabs(meas_ip3d.significance()));
      fts.fill("jet_pfcand_trackjet_dist", -meas_jetdist.value());
      fts.fill("jet_pfcand_trackjet_decayL", meas_decayl.value());
    }
    else {

      fts.fill("jet_pfcand_trackjet_d3d", 0.);
      fts.fill("jet_pfcand_trackjet_d3dsig", 0.);
      fts.fill("jet_pfcand_trackjet_dist", 0.);
      fts.fill("jet_pfcand_trackjet_decayL", 0.);
      
      fts.fill("jet_pfcand_track_chi2", 0.);
      fts.fill("jet_pfcand_track_qual", 0.);
      /*
      fts.fill("jet_pfcand_trackjet_d3d", 0.);
      fts.fill("jet_pfcand_trackjet_d3dsig", 0.);
      fts.fill("jet_pfcand_trackjet_dist", 0.);
      fts.fill("jet_pfcand_trackjet_decayL", 0.);

      fts.fill("jet_pfcand_normchi2", 999);

      fts.fill("jet_pfcand_dptdpt", 0);
      fts.fill("jet_pfcand_detadeta", 0);
      fts.fill("jet_pfcand_dphidphi", 0);
      fts.fill("jet_pfcand_dxydxy", 0);
      fts.fill("jet_pfcand_dzdz", 0);
      fts.fill("jet_pfcand_dxydz", 0);
      fts.fill("jet_pfcand_dphidxy", 0);
      fts.fill("jet_pfcand_dlambdadz", 0);

      fts.fill("jet_pfcand_btagEtaRel", 0);
      fts.fill("jet_pfcand_btagPtRatio", 0);
      fts.fill("jet_pfcand_btagPParRatio", 0);
      fts.fill("jet_pfcand_btagSip2dVal", 0);
      fts.fill("jet_pfcand_btagSip2dSig", 0);
      fts.fill("jet_pfcand_btagSip3dVal", 0);
      fts.fill("jet_pfcand_btagSip3dSig", 0);
      fts.fill("jet_pfcand_btagJetDistVal", 0);
      */
    }
    
    /*
    const reco::Track *trk = nullptr;
    if (packed_cand) {
      trk = packed_cand->bestTrack();
    } else if (reco_cand && useTrackProperties(reco_cand)) {
      trk = reco_cand->bestTrack();
    }
    if (trk) {
      fts.fill("pfcand_normchi2", std::floor(trk->normalizedChi2()));

      // track covariance
      TrackInfoBuilder trkinfo(track_builder_);
      trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, *pv_);
      fts.fill("pfcand_btagEtaRel", trkinfo.getTrackEtaRel());
      fts.fill("pfcand_btagPtRatio", trkinfo.getTrackPtRatio());
      fts.fill("pfcand_btagPParRatio", trkinfo.getTrackPParRatio());
      fts.fill("pfcand_btagSip2dVal", ip_sign * trkinfo.getTrackSip2dVal());
      fts.fill("pfcand_btagSip2dSig", ip_sign * trkinfo.getTrackSip2dSig());
      fts.fill("pfcand_btagSip3dVal", ip_sign * trkinfo.getTrackSip3dVal());
      fts.fill("pfcand_btagSip3dSig", ip_sign * trkinfo.getTrackSip3dSig());
      fts.fill("pfcand_btagJetDistVal", trkinfo.getTrackJetDistVal());
    }
    */

  }

}

void DeepBoostedTauJetTagInfoProducer::fillSVFeatures(DeepBoostedTauJetFeatures &fts, const reco::Jet &jet) {
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

  GlobalVector jet_global_vec(jet.px(), jet.py(), jet.pz());
  const float etasign = jet.eta() > 0 ? 1 : -1;

  for (const auto *sv : jetSVs) {

    fts.fill("sv_mask", 1);
    fts.fill("jet_sv_pt_log", log(sv->pt()));
    fts.fill("jet_sv_eta", sv->eta());
    fts.fill("jet_sv_mass", sv->mass());
    fts.fill("jet_sv_deta", sv->eta() - jet.eta());
    fts.fill("jet_sv_dphi", sv->phi() - jet.phi());
    fts.fill("jet_sv_ntrack", sv->numberOfDaughters());
    fts.fill("jet_sv_chi2", sv->vertexNormalizedChi2());

    reco::Vertex::CovarianceMatrix csv;
    sv->fillVertexCovariance(csv);
    reco::Vertex svtx(sv->vertex(), csv);

    VertexDistanceXY dxy;
    auto valxy = dxy.signedDistance(svtx, *pv_, jet_global_vec);
    fts.fill("jet_sv_dxy", valxy.value());
    fts.fill("jet_sv_dxysig", fabs(valxy.significance()));

    VertexDistance3D d3d;
    auto val3d = d3d.signedDistance(svtx, *pv_, jet_global_vec);
    fts.fill("jet_sv_d3d", val3d.value());
    fts.fill("jet_sv_d3dsig", fabs(val3d.significance()));

  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(DeepBoostedTauJetTagInfoProducer);
