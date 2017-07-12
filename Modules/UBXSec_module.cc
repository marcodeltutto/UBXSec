////////////////////////////////////////////////////////////////////////
// Class:       UBXSec
// Plugin Type: analyzer (art v2_05_00)
// File:        UBXSec_module.cc
//
// Generated at Fri Jan  27 09:44:39 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class UBXSec
 *
 * \ingroup UBXSec
 *
 * \brief Art analyzer module
 * 
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/03/10 $
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Friday, March 10, 2017 at 12:32:31
 *
 */

// Art include
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "uboone/UBXSec/DataTypes/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"

// LArSoft include
#include "uboone/UBFlashFinder/PECalib.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/MCBase/MCDataHolder.h"
#include "lardataobj/MCBase/MCHitCollection.h"

// Algorithms include
#include "uboone/UBXSec/Algorithms/UBXSecHelper.h"
#include "uboone/UBXSec/Algorithms/VertexCheck.h"
#include "uboone/UBXSec/Algorithms/McPfpMatch.h"
#include "uboone/UBXSec/Algorithms/FindDeadRegions.h"

// Root include
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"


namespace ubxsec {
  struct Hit3D_t {
    double x;
    double y;
    double z;
    double q;
  };
}


class UBXSec;


class UBXSec : public art::EDAnalyzer {
public:
  explicit UBXSec(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBXSec(UBXSec const &) = delete;
  UBXSec(UBXSec &&) = delete;
  UBXSec & operator = (UBXSec const &) = delete;
  UBXSec & operator = (UBXSec &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  FindDeadRegions deadRegionsFinder;
  ubxsec::McPfpMatch mcpfpMatcher;
  ::pmtana::PECalib _pecalib;

  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_tag_producer;
  std::string _neutrino_flash_match_producer;
  std::string _cosmic_flash_match_producer;
  std::string _opflash_producer_beam;
  std::string _acpt_producer;
  bool _recursiveMatching = false;
  bool _debug = true;
  int _minimumHitRequirement; ///< Minimum number of hits in at least a plane for a track
  bool _use_genie_info; ///< Turn this off if looking at cosmic only files
  double _beam_spill_start; 
  double _beam_spill_end;

  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;
  const simb::Origin_t COSMIC_ORIGIN   = simb::kCosmicRay;

  bool _is_data, _is_mc;

  /// Maps used for PFParticle truth matching
  typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
  typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
  typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
  typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;


  TTree* _tree1;
  int _run, _subrun, _event;
  int _muon_is_reco;
  double _muon_reco_pur = -9999;
  double _muon_reco_eff = -9999;
  double _true_muon_mom = -9999;
  double _true_muon_mom_matched = -9999;
  int _nPFPtagged, _muon_is_flash_tagged;
  double _muon_tag_score;
  double _fm_score;
  int _fv, _ccnc, _nupdg;
  double _nu_e;
  double _recon_muon_start_x, _recon_muon_start_y, _recon_muon_start_z;
  double _recon_muon_end_x, _recon_muon_end_y, _recon_muon_end_z;
  double _mc_muon_start_x, _mc_muon_start_y, _mc_muon_start_z;
  double _mc_muon_end_x, _mc_muon_end_y, _mc_muon_end_z;
  int _mc_muon_contained;
  double _vtx_resolution;

  int _nslices;
  std::vector<double> _slc_flsmatch_score, _slc_flsmatch_qllx, _slc_flsmatch_tpcx, _slc_flsmatch_t0, _slc_flsmatch_hypoz;
  std::vector<double> _slc_flsmatch_xfixed_chi2, _slc_flsmatch_xfixed_ll;
  std::vector<double> _slc_flsmatch_cosmic_score, _slc_flsmatch_cosmic_t0;
  std::vector<double> _slc_nuvtx_x, _slc_nuvtx_y, _slc_nuvtx_z;
  std::vector<int> _slc_nuvtx_fv;
  std::vector<double> _slc_vtxcheck_angle;
  std::vector<int> _slc_origin;
  std::vector<int> _slc_nhits_u, _slc_nhits_v, _slc_nhits_w;
  std::vector<double> _slc_longesttrack_length;
  std::vector<int> _slc_acpt_outoftime;
  std::vector<int> _slc_crosses_top_boundary;
  std::vector<int> _slc_nuvtx_closetodeadregion_u, _slc_nuvtx_closetodeadregion_v, _slc_nuvtx_closetodeadregion_w;
  std::vector<double> _slc_kalman_chi2;
  std::vector<int> _slc_kalman_ndof;
  std::vector<bool> _slc_passed_min_track_quality;
  std::vector<double> _slc_n_intime_pe_closestpmt;
  std::vector<double> _slc_maxdistance_vtxtrack;

  int _nbeamfls;
  std::vector<double> _beamfls_time, _beamfls_pe, _beamfls_z;
  bool _no_mcflash_but_op_activity; ///< is true if we don't have a neutrino MCFlash in the event, but there is a recon flash in the beam spill
  std::vector<std::vector<double>> _beamfls_spec, _slc_flshypo_spec, _slc_flshypo_xfixed_spec;
  std::vector<double> _numc_flash_spec;
  int _nsignal;
  int _is_swtriggered;

  std::vector<double> _mctrk_start_x, _mctrk_start_y, _mctrk_start_z;
  std::vector<double> _trk_start_x, _trk_start_y, _trk_start_z;
  std::vector<double> _vtx_x, _vtx_y, _vtx_z;
  std::vector<double> _tvtx_x, _tvtx_y, _tvtx_z;

  TTree* _tree2;
  int _total_matches, _nmatch;
  std::vector<double> _hypo_spec, _beam_spec, _fixx_spec;
  double _score;
  int _is_muon;

  TH2F * _deadRegion2P;
  TH2F * _deadRegion3P;
};


UBXSec::UBXSec(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p) 
{
  _pfp_producer                   = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel                 = p.get<std::string>("HitProducer");
  _geantModuleLabel               = p.get<std::string>("GeantModule");
  _spacepointLabel                = p.get<std::string>("SpacePointProducer");
  _neutrino_flash_match_producer  = p.get<std::string>("NeutrinoFlashMatchProducer");
  _cosmic_flash_match_producer    = p.get<std::string>("CosmicFlashMatchProducer");
  _opflash_producer_beam          = p.get<std::string>("OpFlashBeamProducer");
  _acpt_producer                  = p.get<std::string>("ACPTProducer");
    
  _use_genie_info                 = p.get<bool>("UseGENIEInfo", false);
  _minimumHitRequirement          = p.get<int>("MinimumHitRequirement", 3);

  _beam_spill_start               = p.get<double>("BeamSpillStart", 3.2);
  _beam_spill_end                 = p.get<double>("BeamSpillEnd",   4.8);

  _pecalib.Configure(p.get<fhicl::ParameterSet>("PECalib"));

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",                  &_run,                   "run/I");
  _tree1->Branch("subrun",               &_subrun,                "subrun/I");
  _tree1->Branch("event",                &_event,                 "event/I");
  _tree1->Branch("muon_is_reco",         &_muon_is_reco,          "muon_is_reco/I");
  _tree1->Branch("muon_reco_pur",        &_muon_reco_pur,         "muon_reco_pur/D");
  _tree1->Branch("muon_reco_eff",        &_muon_reco_eff,         "muon_reco_eff/D");
  _tree1->Branch("true_muon_mom",        &_true_muon_mom,         "true_muon_mom/D");
  _tree1->Branch("true_muon_mom_matched",&_true_muon_mom_matched, "true_muon_mom_matched/D");
  _tree1->Branch("nPFPtagged",           &_nPFPtagged,            "nPFPtagged/I");
  _tree1->Branch("muon_is_flash_tagged", &_muon_is_flash_tagged,  "muon_is_flash_tagged/I");
  _tree1->Branch("muon_tag_score",       &_muon_tag_score,        "muon_tag_score/D");
  _tree1->Branch("fm_score",             &_fm_score,              "fm_score/D");
  _tree1->Branch("fv",                   &_fv,                    "fv/I");
  _tree1->Branch("ccnc",                 &_ccnc,                  "ccnc/I");
  _tree1->Branch("nupdg",                &_nupdg,                 "nupdg/I");
  _tree1->Branch("nu_e",                 &_nu_e,                  "nu_e/D");
  _tree1->Branch("recon_muon_start_x",   &_recon_muon_start_x,    "recon_muon_start_x/D");
  _tree1->Branch("recon_muon_start_y",   &_recon_muon_start_y,    "recon_muon_start_y/D");
  _tree1->Branch("recon_muon_start_z",   &_recon_muon_start_z,    "recon_muon_start_z/D");
  _tree1->Branch("recon_muon_end_x",     &_recon_muon_end_x,      "recon_muon_end_x/D");
  _tree1->Branch("recon_muon_end_y",     &_recon_muon_end_y,      "recon_muon_end_y/D");
  _tree1->Branch("recon_muon_end_z",     &_recon_muon_end_z,      "recon_muon_end_z/D");
  _tree1->Branch("mc_muon_start_x",      &_mc_muon_start_x,       "mc_muon_start_x/D");
  _tree1->Branch("mc_muon_start_y",      &_mc_muon_start_y,       "mc_muon_start_y/D");
  _tree1->Branch("mc_muon_start_z",      &_mc_muon_start_z,       "mc_muon_start_z/D");
  _tree1->Branch("mc_muon_end_x",        &_mc_muon_end_x,         "mc_muon_end_x/D");
  _tree1->Branch("mc_muon_end_y",        &_mc_muon_end_y,         "mc_muon_end_y/D");
  _tree1->Branch("mc_muon_end_z",        &_mc_muon_end_z,         "mc_muon_end_z/D");
  _tree1->Branch("mc_muon_contained",    &_mc_muon_contained,     "mc_muon_contained/I");
  _tree1->Branch("is_swtriggered",       &_is_swtriggered,        "is_swtriggered/I");
  _tree1->Branch("vtx_resolution",       &_vtx_resolution,        "vtx_resolution/D");

  _tree1->Branch("nslices",                        &_nslices,            "nslices/I");
  _tree1->Branch("slc_flsmatch_score",             "std::vector<double>", &_slc_flsmatch_score);
  _tree1->Branch("slc_flsmatch_qllx",              "std::vector<double>", &_slc_flsmatch_qllx);
  _tree1->Branch("slc_flsmatch_tpcx",              "std::vector<double>", &_slc_flsmatch_tpcx);
  _tree1->Branch("slc_flsmatch_t0",                "std::vector<double>", &_slc_flsmatch_t0);
  _tree1->Branch("slc_flsmatch_hypoz",             "std::vector<double>", &_slc_flsmatch_hypoz);
  _tree1->Branch("slc_flsmatch_xfixed_chi2",       "std::vector<double>", &_slc_flsmatch_xfixed_chi2);
  _tree1->Branch("slc_flsmatch_xfixed_ll",         "std::vector<double>", &_slc_flsmatch_xfixed_ll);
  _tree1->Branch("slc_flsmatch_cosmic_score",      "std::vector<double>", &_slc_flsmatch_cosmic_score);
  _tree1->Branch("slc_flsmatch_cosmic_t0",         "std::vector<double>", &_slc_flsmatch_cosmic_t0);
  _tree1->Branch("slc_nuvtx_x",                    "std::vector<double>", &_slc_nuvtx_x);
  _tree1->Branch("slc_nuvtx_y",                    "std::vector<double>", &_slc_nuvtx_y);
  _tree1->Branch("slc_nuvtx_z",                    "std::vector<double>", &_slc_nuvtx_z);
  _tree1->Branch("slc_nuvtx_fv",                   "std::vector<int>",    &_slc_nuvtx_fv);
  _tree1->Branch("slc_vtxcheck_angle",             "std::vector<double>", &_slc_vtxcheck_angle);
  _tree1->Branch("slc_origin",                     "std::vector<int>",    &_slc_origin);
  _tree1->Branch("slc_nhits_u",                    "std::vector<int>",    &_slc_nhits_u);
  _tree1->Branch("slc_nhits_v",                    "std::vector<int>",    &_slc_nhits_v);
  _tree1->Branch("slc_nhits_w",                    "std::vector<int>",    &_slc_nhits_w);
  _tree1->Branch("slc_longesttrack_length",        "std::vector<double>", &_slc_longesttrack_length);
  _tree1->Branch("slc_acpt_outoftime",             "std::vector<int>",    &_slc_acpt_outoftime);
  _tree1->Branch("slc_crosses_top_boundary",       "std::vector<int>",    &_slc_crosses_top_boundary);
  _tree1->Branch("slc_nuvtx_closetodeadregion_u",  "std::vector<int>",    &_slc_nuvtx_closetodeadregion_u);
  _tree1->Branch("slc_nuvtx_closetodeadregion_v",  "std::vector<int>",    &_slc_nuvtx_closetodeadregion_v);
  _tree1->Branch("slc_nuvtx_closetodeadregion_w",  "std::vector<int>",    &_slc_nuvtx_closetodeadregion_w);
  _tree1->Branch("slc_kalman_chi2",                "std::vector<double>", &_slc_kalman_chi2);
  _tree1->Branch("slc_kalman_ndof",                "std::vector<int>",    &_slc_kalman_ndof);
  _tree1->Branch("slc_passed_min_track_quality",   "std::vector<bool>",   &_slc_passed_min_track_quality);
  _tree1->Branch("slc_n_intime_pe_closestpmt",     "std::vector<double>", &_slc_n_intime_pe_closestpmt);
  _tree1->Branch("slc_maxdistance_vtxtrack",       "std::vector<double>", &_slc_maxdistance_vtxtrack);

  _tree1->Branch("nbeamfls",                   &_nbeamfls,                         "nbeamfls/I");
  _tree1->Branch("beamfls_time",               "std::vector<double>",              &_beamfls_time);
  _tree1->Branch("beamfls_pe",                 "std::vector<double>",              &_beamfls_pe);
  _tree1->Branch("beamfls_z",                  "std::vector<double>",              &_beamfls_z);
  _tree1->Branch("no_mcflash_but_op_activity", &_no_mcflash_but_op_activity,       "no_mcflash_but_op_activity/O");
  _tree1->Branch("beamfls_spec",               "std::vector<std::vector<double>>", &_beamfls_spec);
  _tree1->Branch("numc_flash_spec",            "std::vector<double>",              &_numc_flash_spec);
  _tree1->Branch("slc_flshypo_xfixed_spec",    "std::vector<std::vector<double>>", &_slc_flshypo_xfixed_spec);
  _tree1->Branch("slc_flshypo_spec",           "std::vector<std::vector<double>>", &_slc_flshypo_spec);
  _tree1->Branch("nsignal",                    &_nsignal,                          "nsignal/I");

  _tree1->Branch("mctrk_start_x",        "std::vector<double>", &_mctrk_start_x);
  _tree1->Branch("mctrk_start_y",        "std::vector<double>", &_mctrk_start_y);
  _tree1->Branch("mctrk_start_z",        "std::vector<double>", &_mctrk_start_z);
  _tree1->Branch("trk_start_x",          "std::vector<double>", &_trk_start_x);
  _tree1->Branch("trk_start_y",          "std::vector<double>", &_trk_start_y);
  _tree1->Branch("trk_start_z",          "std::vector<double>", &_trk_start_z);
  _tree1->Branch("vtx_x",                "std::vector<double>", &_vtx_x);
  _tree1->Branch("vtx_y",                "std::vector<double>", &_vtx_y);
  _tree1->Branch("vtx_z",                "std::vector<double>", &_vtx_z);
  _tree1->Branch("tvtx_x",               "std::vector<double>", &_tvtx_x);
  _tree1->Branch("tvtx_y",               "std::vector<double>", &_tvtx_y);
  _tree1->Branch("tvtx_z",               "std::vector<double>", &_tvtx_z);

  _tree2 = fs->make<TTree>("matchtree","");
  _tree2->Branch("run",                &_run,                "run/I");
  _tree2->Branch("subrun",             &_subrun,             "subrun/I");
  _tree2->Branch("event",              &_event,              "event/I");
  _tree2->Branch("total_matches",      &_total_matches,      "total_matches/I");
  _tree2->Branch("nmatch",             &_nmatch,             "nmatch/I");
  _tree2->Branch("score",              &_score,              "score/D");
  _tree2->Branch("hypo_spec",          "std::vector<double>", &_hypo_spec);
  _tree2->Branch("beam_spec",          "std::vector<double>", &_beam_spec);
  _tree2->Branch("fixx_spec",          "std::vector<double>", &_fixx_spec);
  _tree2->Branch("is_muon",            &_is_muon,            "is_muon/I");
  _tree2->Branch("muon_is_reco",       &_muon_is_reco,       "muon_is_reco/I");

  _deadRegion2P = fs->make<TH2F>("deadRegion2P","deadRegion2P", 10350,0.0,1035.0,2300,-115.0,115.0);
  _deadRegion3P = fs->make<TH2F>("deadRegion3P","deadRegion3P", 10350,0.0,1035.0,2300,-115.0,115.0);
}

void UBXSec::analyze(art::Event const & e)
{

  if(_debug) std::cout << "********** UBXSec starts" << std::endl;
  if(_debug) std::cout << "event: " << e.id().event() << std::endl;

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  _is_data = e.isRealData();
  _is_mc   = !_is_data;

  if (_use_genie_info && _is_data) {
    std::cout << "[UBXSec] You have asked to use GENIE info but you are running on a data file.";
    std::cout << " _use_genie_info will be switched to false." << std::endl;
    _use_genie_info = false;
  }

  if (_is_data) {
    std::cout << "[UBXSec] Running on a real data file. No MC-PFP matching will be attempted." << std::endl;
  } else {
    mcpfpMatcher.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);
  }

  art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;

  // Collect tracks
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, pfParticleToTrackMap);

  // Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector  recoParticleVector;
  lar_pandora::PFParticleVector  recoNeutrinoVector;
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

  // Do the MCParticle to PFParticle matching
  lar_pandora::MCParticlesToPFParticles matchedParticles;    // This is a map: MCParticle to matched PFParticle
  lar_pandora::MCParticlesToHits        matchedParticleHits;

  if (_is_mc) {
    mcpfpMatcher.GetRecoToTrueMatches(matchedParticles, matchedParticleHits);
  }


  // *******************
  // MC-PFP Match Study
  // *******************

  std::vector<art::Ptr<recob::PFParticle>> taggedPFP;
  std::vector<double>                      taggedPFPscore;
  std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP;
  std::vector<art::Ptr<recob::PFParticle>> cosmicOriginPFP;
  art::Ptr<recob::PFParticle>              muonPFP;
  art::Ptr<simb::MCParticle>               muonMCParticle;
  _muon_is_reco = 0;
  _mc_muon_contained = 0;

  if (_is_data) goto doanalysis;

  neutrinoOriginPFP.clear(); 
  cosmicOriginPFP.clear();

  // Loop over true particle and find the cosmic related ones
  for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
      iter1 != iterEnd1; ++iter1) {

    art::Ptr<simb::MCParticle>  mc_par = iter1->first;   // The MCParticle 
    art::Ptr<recob::PFParticle> pf_par = iter1->second;  // The matched PFParticle

    const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());

    if (!mc_truth) {
      std::cerr << "[UBXSec] Problem with MCTruth pointer." << std::endl;
      continue;
    }

    if (mc_truth->Origin() == COSMIC_ORIGIN) {

      cosmicOriginPFP.emplace_back(pf_par);

      double end[3];
      end[0] = mc_par->EndX();
      end[1] = mc_par->EndY();
      end[2] = mc_par->EndZ();
      if ( (mc_par->PdgCode() == 13 || mc_par->PdgCode() == -13) && UBXSecHelper::InFV(end) ){
        std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Is stopping muon" << std::endl;

        lar_pandora::VertexVector          vertexVector;
        lar_pandora::PFParticlesToVertices particlesToVertices;
        lar_pandora::LArPandoraHelper::CollectVertices(e, _pfp_producer, vertexVector, particlesToVertices);

        lar_pandora::VertexVector vertex_v = particlesToVertices.find(pf_par)->second;
        double xyz[3];
        vertex_v[0]->XYZ(xyz);

        std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM The PFP has vtx x="<<xyz[0]<<" y="<<xyz[1]<<" z="<<xyz[2] << std::endl;
      }
         
    }
  }

  // Loop over true particle and find the neutrino related ones
  for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
             iter1 != iterEnd1; ++iter1) {

     art::Ptr<simb::MCParticle>  mc_par = iter1->first;   // The MCParticle 
     art::Ptr<recob::PFParticle> pf_par = iter1->second;  // The matched PFParticle

     const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());

     if (!mc_truth) {
       std::cerr << "[UBXSec] Problem with MCTruth pointer." << std::endl;
       continue;
     }

     if (mc_truth->Origin() == NEUTRINO_ORIGIN) {
       if (_debug) {
         std::cout << "Neutrino related track found." << std::endl;
         std::cout << "Process (0==CC, 1==NC) " << mc_truth->GetNeutrino().CCNC()         << std::endl;
         std::cout << "Neutrino PDG           " << mc_truth->GetNeutrino().Nu().PdgCode() << std::endl;
         std::cout << "PDG  " << mc_par->PdgCode() << std::endl;
         std::cout << "Mass " << mc_par->Mass()    << std::endl;
         std::cout << "Proc " << mc_par->Process() << std::endl;
         std::cout << "Vx   " << mc_par->Vx()      << std::endl;
         std::cout << "Vy   " << mc_par->Vy()      << std::endl;
         std::cout << "Vz   " << mc_par->Vz()      << std::endl;
         std::cout << "T    " << mc_par->T()       << std::endl;
         double timeCorrection = 343.75;
         std::cout << "Remeber a time correction of " << timeCorrection << std::endl;
         auto iter =  matchedParticleHits.find(mc_par);
         std::cout << "Related hits: " << (iter->second).size() << std::endl;
       }    
       if (_debug) {
         std::cout << "  The related PFP: " << std::endl;
         std::cout << "  has ID: " << pf_par->Self() << std::endl;
       }

       neutrinoOriginPFP.emplace_back(pf_par);

       // If we matched a muon
       if (mc_par->PdgCode() == 13 && mc_par->Mother() == 0) {
         muonMCParticle = mc_par;
         muonPFP = pf_par;
         _muon_is_reco = 1;


         // Muon track puritity and efficiency
         _muon_reco_pur = _muon_reco_eff = -9999;
         auto iter = recoParticlesToHits.find(pf_par);
         if (iter != recoParticlesToHits.end()) {
           UBXSecHelper::GetTrackPurityAndEfficiency((*iter).second, _muon_reco_pur, _muon_reco_eff);
         }
         _true_muon_mom_matched = mc_par->P();
         //std::cout << "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY efficiency: " << eff << "  purity "  << pur << std::endl;

         lar_pandora::PFParticlesToTracks::const_iterator it =  pfParticleToTrackMap.find(pf_par);
         if (it != pfParticleToTrackMap.end()) {

           lar_pandora::TrackVector trk_v = it->second;
           std::cout << "Track vector length is: " << trk_v.size() << std::endl;
           art::Ptr<recob::Track>   trk   = trk_v[0];

           std::cout << "Recon muon start x " << trk->Vertex().X() << std::endl;
           std::cout << "Recon muon start y " << trk->Vertex().Y() << std::endl;
           std::cout << "Recon muon start z " << trk->Vertex().Z() << std::endl;

           _recon_muon_start_x = trk->Vertex().X();
           _recon_muon_start_y = trk->Vertex().Y();
           _recon_muon_start_z = trk->Vertex().Z();

           std::cout << "Recon muon end x " << trk->End().X() << std::endl;
           std::cout << "Recon muon end y " << trk->End().Y() << std::endl;
           std::cout << "Recon muon end z " << trk->End().Z() << std::endl;

           _recon_muon_end_x = trk->End().X();
           _recon_muon_end_y = trk->End().Y();
           _recon_muon_end_z = trk->End().Z();

           std::cout << "MC muon start x " << mc_par->Vx() << std::endl;
           std::cout << "MC muon start y " << mc_par->Vy() << std::endl;
           std::cout << "MC muon start z " << mc_par->Vz() << std::endl;

           _mc_muon_start_x = mc_par->Vx();
           _mc_muon_start_y = mc_par->Vy();
           _mc_muon_start_z = mc_par->Vz();

           std::cout << "MC muon end x " << mc_par->EndX() << std::endl;
           std::cout << "MC muon end y " << mc_par->EndY() << std::endl;
           std::cout << "MC muon end z " << mc_par->EndZ() << std::endl;
       
           _mc_muon_end_x = mc_par->EndX();
           _mc_muon_end_y = mc_par->EndY();
           _mc_muon_end_z = mc_par->EndZ();
         }


         double start[3] = {_mc_muon_start_x, _mc_muon_start_y, _mc_muon_start_z};
         double stop[3]  = {_mc_muon_end_x,   _mc_muon_end_y,   _mc_muon_end_z};

         if (UBXSecHelper::InFV(start) && UBXSecHelper::InFV(stop))
           _mc_muon_contained = 1;
       }  
     }
  }
  if (_debug) std::cout << "Neutrino related PFPs in this event: " << neutrinoOriginPFP.size() << std::endl;




  // *******************
  // Analysis
  // *******************

  doanalysis:

  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(_pfp_producer, pfpHandle);
  art::FindManyP<ubana::FlashMatch> pfpToNeutrinoFlashMatchAssns(pfpHandle, e, _neutrino_flash_match_producer);  

  // ACPT
  art::Handle<std::vector<anab::T0> > t0_h;
  e.getByLabel(_acpt_producer,t0_h);
  if(!t0_h.isValid()) {
    std::cout << "[UBXSec] T0 product not found..." << std::endl;
    //throw std::exception();
  }
  if(t0_h->empty()) {
    std::cout << "[UBXSec] t0 is empty." << std::endl;
  }
   
  art::FindManyP<recob::Track> track_ptr_coll_v(t0_h, e, _acpt_producer); 

  art::Handle<std::vector<recob::Track>> track_h;
  e.getByLabel(_pfp_producer,track_h);
  if (!track_h.isValid() || track_h->empty()) {
    std::cout << "[UBXSec] Track handle is not valid or empty." << std::endl;
    //throw std::exception();
  }

  art::FindManyP<recob::OpFlash> opfls_ptr_coll_v(track_h, e, _acpt_producer);
  

  // Kalman Track
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);
  if(!pfp_h.isValid()){
    std::cout << "Track product " << _pfp_producer << " not found..." << std::endl;
    throw std::exception();
  }
  if(pfp_h->empty()) {
    std::cout << "PFP " << _pfp_producer << " is empty." << std::endl;
  }

  art::FindManyP<recob::Track> trk_kalman_v(pfp_h, e, "pandoraNuKalmanTrack");


  // Check if golden
  /*
  bool is_golden = false;
  for (auto const & mctrk : (*mctrk_h)) {
    if (mctrk.Origin() == NEUTRINO_ORIGIN && mctrk.PdgCode() == 14) {
      std::cout << "PPPPPPProcess is " << mctrk.Process() << std::endl;
      for (size_t pt = 0; pt < mctrk.NumberTrajectoryPoints(); pt++){
        float ptNearDeadRegion = 0;
        if (NearDeadReg2P( mctrk.Vy(pt), mctrk.Vz(pt), 0.6) {
          ptNearDeadRegion++;
        }
      } // loop over trj points
      if (ptNearDeadRegion/(float)mctrk.NumberTrajectoryPoints() > 0.05) {
        is_golden = false; 
        break;
      } 
      double start[3] = {mctrk.Vx(),   mctrk.Vy(),   mctrk.Vz()};
      double end[3]   = {mctrk.EndX(), mctrk.EndY(), mctrk.EndZ()};
      if (UBXSecHelper::InFV(start) && UBXSecHelper::InFV(end)) {
        is_golden = true;
        break;
      }
    }
  }
  std::cout << "       is good track? " << is_golden << std::endl;
  */

  // Check if truth nu is in FV
  // Collecting GENIE particles
  if(_use_genie_info) {
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (e.getByLabel("generator",mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

    int iList = 0; // 1 nu int per spill
    double truth_nu_vtx[3] = {mclist[iList]->GetNeutrino().Nu().Vx(),
                              mclist[iList]->GetNeutrino().Nu().Vy(),
                              mclist[iList]->GetNeutrino().Nu().Vz()};
    if (UBXSecHelper::InFV(truth_nu_vtx)) _fv = 1;
    else _fv = 0;
    _ccnc    = mclist[iList]->GetNeutrino().CCNC();
    _nupdg   = mclist[iList]->GetNeutrino().Nu().PdgCode();
    _nu_e    = mclist[iList]->GetNeutrino().Nu().E();

    _tvtx_x.clear(); _tvtx_y.clear(); _tvtx_z.clear();
    for(size_t n = 0; n < mclist.size(); n++ ) {
      _tvtx_x.emplace_back(mclist[n]->GetNeutrino().Nu().Vx());
      _tvtx_y.emplace_back(mclist[n]->GetNeutrino().Nu().Vy());
      _tvtx_z.emplace_back(mclist[n]->GetNeutrino().Nu().Vz());
    }

    _nsignal = 0;
    if(_nupdg==14 && _ccnc==0 && _fv==1) _nsignal=1; 

    // Also save muon momentum if is signal
    _true_muon_mom = -9999.;
    if (_nsignal == 1) {
      for (int p = 0; p < mclist[iList]->NParticles(); p++) {
        auto const & mcp = mclist[iList]->GetParticle(p);
        if (mcp.Mother() != 0) continue;
        if (mcp.PdgCode() != 13) continue;
        _true_muon_mom = mcp.P();
      }

    }
  }

  // OpHits related 
  lar_pandora::PFParticlesToSpacePoints pfp_to_spacept;
  lar_pandora::SpacePointsToHits spacept_to_hits;

  lar_pandora::PFParticleVector temp2;
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, temp2, pfp_to_spacept);

  lar_pandora::SpacePointVector temp3;
  lar_pandora::LArPandoraHelper::CollectSpacePoints (e, _pfp_producer, temp3, spacept_to_hits);

  art::Handle<std::vector<recob::OpHit>> ophit_h;
  e.getByLabel("ophitBeam", ophit_h);
  if(!ophit_h.isValid()) {
    std::cout << "[UBXSec] Cannot locate OpHits." << std::endl;
  }
  // Save the number of slices in this event
  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;
  UBXSecHelper::GetTPCObjects(e, _pfp_producer, pfp_v_v, track_v_v);

  _nslices = pfp_v_v.size();
  _slc_flsmatch_score.resize(_nslices, -9999);
  _slc_flsmatch_qllx.resize(_nslices, -9999);
  _slc_flsmatch_tpcx.resize(_nslices, -9999);
  _slc_flsmatch_t0.resize(_nslices, -9999);
  _slc_flsmatch_hypoz.resize(_nslices, -9999);
  _slc_flsmatch_xfixed_chi2.resize(_nslices, -9999);
  _slc_flsmatch_xfixed_ll.resize(_nslices, -9999);
  _slc_nuvtx_x.resize(_nslices);
  _slc_nuvtx_y.resize(_nslices);
  _slc_nuvtx_z.resize(_nslices);
  _slc_nuvtx_fv.resize(_nslices);
  _slc_vtxcheck_angle.resize(_nslices);
  _slc_origin.resize(_nslices);
  _slc_flshypo_xfixed_spec.resize(_nslices);
  _slc_flshypo_spec.resize(_nslices);
  _slc_nhits_u.resize(_nslices, -9999);
  _slc_nhits_v.resize(_nslices, -9999);
  _slc_nhits_w.resize(_nslices, -9999);
  _slc_flsmatch_cosmic_score.resize(_nslices, -9999);
  _slc_flsmatch_cosmic_t0.resize(_nslices, -9999);
  _slc_longesttrack_length.resize(_nslices, -9999);
  _slc_acpt_outoftime.resize(_nslices, -9999);
  _slc_crosses_top_boundary.resize(_nslices, -9999);
  _slc_nuvtx_closetodeadregion_u.resize(_nslices, -9999);
  _slc_nuvtx_closetodeadregion_v.resize(_nslices, -9999);
  _slc_nuvtx_closetodeadregion_w.resize(_nslices, -9999);  
  _slc_kalman_chi2.resize(_nslices, -9999);
  _slc_kalman_ndof.resize(_nslices, -9999);  
  _slc_passed_min_track_quality.resize(_nslices, -9999);
  _slc_n_intime_pe_closestpmt.resize(_nslices, -9999);
  _slc_maxdistance_vtxtrack.resize(_nslices, -9999);

  std::cout << "UBXSec - SAVING INFORMATION" << std::endl;
  _vtx_resolution = -9999;
  for (unsigned int slice = 0; slice < pfp_v_v.size(); slice++){
    std::cout << ">>> SLICE" << slice << std::endl;

    // Slice origin (0 is neutrino, 1 is cosmic)
    _slc_origin[slice] = UBXSecHelper::GetSliceOrigin(neutrinoOriginPFP, cosmicOriginPFP, pfp_v_v[slice]);

    // Reco vertex
    double reco_nu_vtx[3];
    UBXSecHelper::GetNuVertexFromTPCObject(e, _pfp_producer, pfp_v_v[slice], reco_nu_vtx);
    _slc_nuvtx_x[slice] = reco_nu_vtx[0];
    _slc_nuvtx_y[slice] = reco_nu_vtx[1];
    _slc_nuvtx_z[slice] = reco_nu_vtx[2];
    _slc_nuvtx_fv[slice] = (UBXSecHelper::InFV(reco_nu_vtx) ? 1 : 0);
    std::cout << "    Reco vertex saved" << std::endl;

    // Vertex resolution
    if (_slc_origin[slice] == 0) {
      _vtx_resolution = sqrt( pow(_slc_nuvtx_y[slice]-_tvtx_y[0], 2) + pow(_slc_nuvtx_z[slice]-_tvtx_z[0], 2) );
    } 

    // Neutrino Flash match
    _slc_flsmatch_score[slice] = -9999;
    art::Ptr<recob::PFParticle> NuPFP = UBXSecHelper::GetNuPFP(pfp_v_v[slice]);
    std::vector<art::Ptr<ubana::FlashMatch>> pfpToFlashMatch_v = pfpToNeutrinoFlashMatchAssns.at(NuPFP.key());
    if (pfpToFlashMatch_v.size() > 1) {
      std::cout << "    More than one flash match per nu pfp ?!" << std::endl;
      continue;
    } else if (pfpToFlashMatch_v.size() == 0){
      // do nothing
    } else {
      _slc_flsmatch_score[slice]       = pfpToFlashMatch_v[0]->GetScore(); 
      _slc_flsmatch_qllx[slice]        = pfpToFlashMatch_v[0]->GetEstimatedX();
      _slc_flsmatch_tpcx[slice]        = pfpToFlashMatch_v[0]->GetTPCX();
      _slc_flsmatch_t0[slice]          = pfpToFlashMatch_v[0]->GetT0();
      _slc_flsmatch_hypoz[slice]       = UBXSecHelper::GetFlashZCenter(pfpToFlashMatch_v[0]->GetHypoFlashSpec());
      _slc_flsmatch_xfixed_chi2[slice] = pfpToFlashMatch_v[0]->GetXFixedChi2();
      _slc_flsmatch_xfixed_ll[slice]   = pfpToFlashMatch_v[0]->GetXFixedLl();
      _slc_flshypo_xfixed_spec[slice]  = pfpToFlashMatch_v[0]->GetXFixedHypoFlashSpec();
      _slc_flshypo_spec[slice]         = pfpToFlashMatch_v[0]->GetHypoFlashSpec();
      for (auto v : _slc_flshypo_spec[slice]) std::cout << "PE: " << v << std::endl;

      std::cout << "    FM score: " << _slc_flsmatch_score[slice] << std::endl;
    }

    // Cosmic Flash Match
    _slc_flsmatch_cosmic_score[slice] = -9999;
    /*
    std::vector<art::Ptr<ubana::FlashMatch>> pfpToCosmicFlashMatch_v = pfpToCosmicFlashMatchAssns.at(NuPFP.key());
    if (pfpToCosmicFlashMatch_v.size() > 1) {
      std::cout << "    More than one flash match per nu pfp!" << std::endl;
      continue;
    } else if (pfpToCosmicFlashMatch_v.size() == 0){
      std::cout << "    PFP to flash match ass for cosmic is zero." << std::endl;
      //continue;
    } else if (pfpToCosmicFlashMatch_v.size() == 1){
      //std::cout << "pfpToCosmicFlashMatch_v[0]->GetScore() is " << pfpToCosmicFlashMatch_v[0]->GetScore() << std::endl;
      //std::cout << "pfpToCosmicFlashMatch_v[0]->GetT0() is " << pfpToCosmicFlashMatch_v[0]->GetT0() << std::endl;
      _slc_flsmatch_cosmic_score[slice] = pfpToCosmicFlashMatch_v[0]->GetScore();
      _slc_flsmatch_cosmic_t0[slice]    = pfpToCosmicFlashMatch_v[0]->GetT0();
    } else {
      std::cout << "    I don't know what fucking case this is." << std::endl;
    }
    */

    // Hits
    int nhits_u, nhits_v, nhits_w;
    UBXSecHelper::GetNumberOfHitsPerPlane(e, _pfp_producer, track_v_v[slice], nhits_u, nhits_v, nhits_w);
    _slc_nhits_u[slice] = nhits_u;
    _slc_nhits_v[slice] = nhits_v;
    _slc_nhits_w[slice] = nhits_w;

    // Longest track and check boundary
    recob::Track lt;
    if (UBXSecHelper::GetLongestTrackFromTPCObj(track_v_v[slice], lt)){
      _slc_longesttrack_length[slice] = lt.Length();
      int vtx_ok;
      _slc_crosses_top_boundary[slice] = (UBXSecHelper::IsCrossingTopBoundary(lt, vtx_ok) ? 1 : 0);
    } else {
      _slc_longesttrack_length[slice] = -9999;
    }

    // ACPT
    _slc_acpt_outoftime[slice] = 0;
    for (unsigned int t = 0; t < track_v_v[slice].size(); t++) {
      if(opfls_ptr_coll_v.at(track_v_v[slice][t].key()).size()>1) {
        std::cout << "[UBXSec] More than 1 association found (ACPT)!" << std::endl;
        throw std::exception();
      } else if (opfls_ptr_coll_v.at(track_v_v[slice][t].key()).size()==0){
        continue;
      } else {
        art::Ptr<recob::OpFlash> flash_ptr = opfls_ptr_coll_v.at(track_v_v[slice][t].key()).at(0);
        if (flash_ptr->Time() < _beam_spill_start || flash_ptr->Time() > _beam_spill_end) {
          _slc_acpt_outoftime[slice] = 1;
        }
      }
    }

    // Track quality
    _slc_kalman_chi2[slice] = -9999;
    for (unsigned int t = 0; t < pfp_v_v[slice].size(); t++) {
      if(trk_kalman_v.at(pfp_v_v[slice][t].key()).size()>1) {
        std::cout << "[UBXSec] TQ more than one track per PFP, ntracks " << trk_kalman_v.at(pfp_v_v[slice][t].key()).size() << std::endl;
      } else if (trk_kalman_v.at(pfp_v_v[slice][t].key()).size()==0){
        continue;
      } else {
        art::Ptr<recob::Track> trk_ptr = trk_kalman_v.at(pfp_v_v[slice][t].key()).at(0);
        _slc_kalman_chi2[slice] = trk_ptr->Chi2();
        _slc_kalman_ndof[slice] = trk_ptr->Ndof();
      }
    }
    bool goodTrack = false;
    for (auto trk : track_v_v[slice]) {
      if (!deadRegionsFinder.NearDeadReg2P( (trk->Vertex()).Y(), (trk->Vertex()).Z(), 0.6 )  &&
          !deadRegionsFinder.NearDeadReg2P( (trk->End()).Y(),    (trk->End()).Z(),    0.6 )  &&
          UBXSecHelper::TrackPassesHitRequirment(e, _pfp_producer, trk, _minimumHitRequirement) ) {
        goodTrack = true;
        continue;
      }
    }
    if (goodTrack) _slc_passed_min_track_quality[slice] = true;
    else _slc_passed_min_track_quality[slice] = false;

    // Channel status
    _slc_nuvtx_closetodeadregion_u[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 0) ? 1 : 0);
    _slc_nuvtx_closetodeadregion_v[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 1) ? 1 : 0);
    _slc_nuvtx_closetodeadregion_w[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 2) ? 1 : 0);

    // Vertex check
    recob::Vertex slice_vtx;
    UBXSecHelper::GetNuVertexFromTPCObject(e, _pfp_producer, pfp_v_v[slice], slice_vtx);
    ubxsec::VertexCheck vtxCheck(track_v_v[slice], slice_vtx);
    _slc_vtxcheck_angle[slice] = vtxCheck.AngleBetweenLongestTracks();
    
    // OpHits
    std::vector<ubxsec::Hit3D_t> hit3d_v;
    hit3d_v.clear();
    for (auto pfp : pfp_v_v[slice]) {
      auto iter = pfp_to_spacept.find(pfp);
      if (iter != pfp_to_spacept.end()) {
        //std::cout << "[UBXSec] Found related spacepoints, size is " << (iter->second).size() << std::endl;
      } else {
        std::cout << "[UBXSec] Can't find spacepoints for pfp with pdg " << pfp->PdgCode() << std::endl;
        continue;
      }
      // Loop through the hits associated 
      for (auto sp_pt : (iter->second)) {
        auto iter2 = spacept_to_hits.find(sp_pt);
        if (iter2 != spacept_to_hits.end()) {
          //std::cout << "[UBXSec] Founds hits associated to this sp_pt" << std::endl;
        } else {
          std::cout << "[UBXSec] Can't find hits ass to this sp_pt" << std::endl;
          continue;
        }   
        // Save sp_pt position and hit charge for all the sp_pt you have
        auto hit = iter2->second;
        ubxsec::Hit3D_t thishit;
        thishit.x = sp_pt->XYZ()[0];
        thishit.y = sp_pt->XYZ()[1];
        thishit.z = sp_pt->XYZ()[2];
        thishit.q = hit->Integral();
        hit3d_v.emplace_back(thishit);
      }
    }
    std::cout << "[UBXSec] For this TPC object we have " << hit3d_v.size() << " Hit3D_t hits." << std::endl;

    // Now construct average position
    double sumx = 0, sumy = 0, sumz = 0;
    double totq = 0;
    for (auto hit3d : hit3d_v) {
      sumx += hit3d.q * hit3d.x;
      sumy += hit3d.q * hit3d.y;
      sumz += hit3d.q * hit3d.z;
      
      totq += hit3d.q; 
    }
    double charge_center[3] = {sumx / totq, sumy / totq, sumz / totq};

    int this_opch = UBXSecHelper::GetClosestPMT(charge_center);

    // Look at the opHits from this pmt
    int n_intime_ophits = 0;
    double n_intime_pe = 0;
    for (size_t oh = 0; oh < ophit_h->size(); oh++) {
      auto const & ophit = (*ophit_h)[oh];
      if (ophit.OpChannel() != this_opch) continue;
      if (ophit.PeakTime() > _beam_spill_start && ophit.PeakTime() < _beam_spill_end) {
        n_intime_ophits ++;
        size_t opdet = geo->OpDetFromOpChannel(ophit.OpChannel());
        n_intime_pe += _pecalib.BeamPE(opdet,ophit.Area(),ophit.Amplitude());
      }     
    } // end loop ophit

    _slc_n_intime_pe_closestpmt[slice] = n_intime_pe;


    // Distance from recon nu vertex to thefar away track in TPCObject
    //_slc_maxdistance_vtxtrack = UBXSecHelper::GetMaxTrackVertexDistance();


    std::cout << "UBXSec - INFORMATION SAVED" << std::endl;
  } // slice loop



  /* Dead regions
  //art::ServiceHandle<FindDeadRegions> deadRegionsFinder;
  FindDeadRegions deadRegionsFinder;
  deadRegionsFinder.GetDeadRegionHisto2P(_deadRegion2P);
  deadRegionsFinder.GetDeadRegionHisto3P(_deadRegion3P);
  */

  // Flashes
  ::art::Handle<std::vector<recob::OpFlash>> beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);
  if( !beamflash_h.isValid() || beamflash_h->empty() ) {
    std::cerr << "Don't have good flashes." << std::endl;
  }

  _nbeamfls = beamflash_h->size();
  _beamfls_pe.resize(_nbeamfls);
  _beamfls_time.resize(_nbeamfls);
  _beamfls_z.resize(_nbeamfls);
  _beamfls_spec.resize(_nbeamfls);

  for (size_t n = 0; n < beamflash_h->size(); n++) {
    auto const& flash = (*beamflash_h)[n];
    _beamfls_pe[n]   = flash.TotalPE();
    _beamfls_time[n] = flash.Time();
    _beamfls_z[n]    = flash.ZCenter();

    _beamfls_spec[n].resize(32);
    for (unsigned int i = 0; i < 32; i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      _beamfls_spec[n][opdet] = flash.PE(i);
    }
  }
  



  if (_is_mc) {

    // SW Trigger
    art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
    e.getByLabel("swtrigger", softwareTriggerHandle);

    if (softwareTriggerHandle.isValid()){ 
      if (softwareTriggerHandle->getNumberOfAlgorithms() == 1) {
        std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
        std::cout << "SW trigger name: " << algoNames[0] << std::endl;
        _is_swtriggered = (softwareTriggerHandle->passedAlgo(algoNames[0]) ? 1 : 0);
      }
    }

    // MC Flash
    ::art::Handle<std::vector<recob::OpFlash> > nuMcflash_h;
    e.getByLabel("NeutrinoMCFlash",nuMcflash_h);
    if( !nuMcflash_h.isValid() || nuMcflash_h->empty() ) {
      std::cerr << "Don't have neutrino MC flashes." << std::endl;
      //return;
    } else {

      auto const& flash = (*nuMcflash_h)[0];
      _numc_flash_spec.resize(geo->NOpDets());
      for (unsigned int i = 0; i < geo->NOpDets(); i++) {
        unsigned int opdet = geo->OpDetFromOpChannel(i);
        _numc_flash_spec[opdet] = flash.PE(i);
      }
    }

    // MCFlash vs op activity
    bool opActivityInBeamSpill = false;
    // Check if there are recon beam flashed in the beam spill window
    for (auto reco_fls_time : _beamfls_time) {
      if (reco_fls_time > _beam_spill_start && reco_fls_time < _beam_spill_end) {
         opActivityInBeamSpill = true;
       }
    }
    if (nuMcflash_h->size() > 0) {
      std::cout << "We have a neutrino MCFlash, and its time is: " << (*nuMcflash_h)[0].Time() << std::endl;    
    } else if (nuMcflash_h->size() == 0) {
      if(opActivityInBeamSpill) {
        std::cout << "No MCFlash but optical activity in the beam spill." << std::endl;
        _no_mcflash_but_op_activity = true;
      }
    }

  }


  /* MCHits
  ::art::Handle< std::vector<sim::MCHitCollection> > mcHit_h;
  e.getByLabel("mchitfinder",mcHit_h);
  if( !mcHit_h.isValid() || mcHit_h->empty() ) {
    std::cerr << "Don't have MCHits." << std::endl;
    return;
  }

  std::map<int,int> geantTrackIDToNumberOfMCHit;

  std::cout << "Number of MCHits: " << mcHit_h->size() << std::endl;
  for (auto const & mchit_v : *mcHit_h){
    std::cout << "CHANNEL " << mchit_v.Channel() << std::endl;
    auto iter = geantTrackIDToNumberOfMCHit.find(1480857);
    if (iter != geantTrackIDToNumberOfMCHit.end())
      std::cout << "   For this channel the muon has # of MCHits = " << (*iter).second << std::endl;
    for (auto const & mchit : mchit_v){
      //std::cout << "Hit track id: " << mchit.PartTrackId() << std::endl;
      auto iter = geantTrackIDToNumberOfMCHit.find(mchit.PartTrackId());
      if (iter != geantTrackIDToNumberOfMCHit.end()) (*iter).second += 1;
      else geantTrackIDToNumberOfMCHit[mchit.PartTrackId()] = 1;
      if (mchit.PartTrackId() == 1480857)
        std::cout << "From mchit I get that the muon has energy: " << mchit.PartEnergy() << " and the hit time is " << mchit.PeakTime() << std::endl;
    }
  }

  //for(auto const& key_value : geantTrackIDToNumberOfMCHit){
    //std::cout << "Track id: " << key_value.first << "  number of hits: " << key_value.second << std::endl;
  //} 


  std::cout << "The muon MC particle track id is: " << muonMCParticle->TrackId() << std::endl;
  auto iter = geantTrackIDToNumberOfMCHit.find(muonMCParticle->TrackId());
  if (iter != geantTrackIDToNumberOfMCHit.end()) {
    const simb::MCParticle * mcpar_bt = bt->TrackIDToParticle((*iter).first);
    std::cout << "Number of MCHits for the muon: " << (*iter).second << std::endl;
    std::cout << "The backtracked particle has pdg: " << mcpar_bt->PdgCode() << std::endl;
  }

  auto titer = matchedParticleHits.find(muonMCParticle);
  if (titer != matchedParticleHits.end()){
    std::cout << "Number of muon recon hits: " << ((*titer).second).size() << std::endl;
  }
  */



 


  if(_debug) std::cout << "[UBXSec] Filling tree now." << std::endl;
  _tree1->Fill();

  if(_debug) std::cout << "********** UBXSec ends" << std::endl;
}



DEFINE_ART_MODULE(UBXSec)
