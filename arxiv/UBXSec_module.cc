////////////////////////////////////////////////////////////////////////
// Class:       UBXSec
// Plugin Type: producer 
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
 * \brief Art producer module
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
#include "art/Framework/Core/EDProducer.h"
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
#include "uboone/UBXSec/DataTypes/MCGhost.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"

#include "uboone/UBXSec/DataTypes/UBXSecEvent.h"

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
//#include "uboone/UBXSec/Algorithms/McPfpMatch.h"
#include "uboone/UBXSec/Algorithms/FindDeadRegions.h"
#include "uboone/UBXSec/Algorithms/MuonCandidateFinder.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"

// Root include
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"

#include <fstream>


namespace ubxsec {
  struct Hit3D_t {
    double x;
    double y;
    double z;
    double q;
  };
}


class UBXSec;


class UBXSec : public art::EDProducer {
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
  void produce(art::Event & e) override;
  void endSubRun(art::SubRun &sr) override;

private:

  /// Prints MC particles from GENIE on the screen
  void PrintMC(std::vector<art::Ptr<simb::MCTruth>> mclist);

  FindDeadRegions deadRegionsFinder;
  //ubxsec::McPfpMatch mcpfpMatcher;
  ::ubana::FiducialVolume _fiducial_volume;
  ::pmtana::PECalib _pecalib;

  // To be set via fcl parameters
  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_tag_producer;
  std::string _neutrino_flash_match_producer;
  std::string _cosmic_flash_match_producer;
  std::string _opflash_producer_beam;
  std::string _acpt_producer;
  std::string _tpcobject_producer;
  std::string _potsum_producer;
  std::string _potsum_instance;
  std::string _particle_id_producer;
  std::string _mc_ghost_producer;
  std::string _geocosmictag_producer;
  bool _debug = true;                   ///< Debug mode
  int _minimumHitRequirement;           ///< Minimum number of hits in at least a plane for a track
  bool _use_genie_info;                 ///< Turn this off if looking at cosmic only files
  double _beam_spill_start;             ///< Start time of the beam spill (us) 
  double _beam_spill_end;               ///< Start time of the beam spill (us) 
  double _total_pe_cut;                 ///< PE cut to be applied to beam flash
  double _geo_cosmic_score_cut;         ///< Cut on the score of the pandoraNu geo cosmic tagger
  double _tolerance_track_multiplicity; ///< Tolerance to consider a track coming from the nu reco vertex

  // Constants
  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;
  const simb::Origin_t COSMIC_ORIGIN   = simb::kCosmicRay;

  // To be filled within module
  bool _is_data, _is_mc;
  double _candidate_flash_time;
  double _drift_velocity;

  // Outputs trees
  TTree* _tree1;
  UBXSecEvent *ubxsec_event = new UBXSecEvent();
 
  /*
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
  bool _is_signal;
  double _mc_muon_contained;
  double _vtx_resolution;

  int _nslices;
  std::vector<double> _slc_flsmatch_score, _slc_flsmatch_qllx, _slc_flsmatch_tpcx, _slc_flsmatch_t0, _slc_flsmatch_hypoz;
  std::vector<double> _slc_flsmatch_xfixed_chi2, _slc_flsmatch_xfixed_ll;
  std::vector<double> _slc_flsmatch_cosmic_score, _slc_flsmatch_cosmic_t0;
  std::vector<double> _slc_nuvtx_x, _slc_nuvtx_y, _slc_nuvtx_z;
  std::vector<int> _slc_nuvtx_fv;
  std::vector<double> _slc_vtxcheck_angle;
  std::vector<int> _slc_origin, _slc_origin_extra;
  std::vector<int> _slc_nhits_u, _slc_nhits_v, _slc_nhits_w;
  std::vector<double> _slc_longesttrack_length;
  std::vector<double> _slc_longesttrack_phi, _slc_longesttrack_theta;
  std::vector<bool> _slc_longesttrack_iscontained;
  std::vector<bool> _slc_muoncandidate_exists;
  std::vector<double> _slc_muoncandidate_length;
  std::vector<double> _slc_muoncandidate_phi, _slc_muoncandidate_theta;
  std::vector<int> _slc_acpt_outoftime;
  std::vector<int> _slc_crosses_top_boundary;
  std::vector<int> _slc_nuvtx_closetodeadregion_u, _slc_nuvtx_closetodeadregion_v, _slc_nuvtx_closetodeadregion_w;
  std::vector<double> _slc_kalman_chi2;
  std::vector<int> _slc_kalman_ndof;
  std::vector<bool> _slc_passed_min_track_quality, _slc_passed_min_vertex_quality;
  std::vector<double> _slc_n_intime_pe_closestpmt;
  std::vector<double> _slc_maxdistance_vtxtrack;
  std::vector<int> _slc_npfp, _slc_ntrack, _slc_nshower;
  std::vector<bool> _slc_iscontained;
  std::vector<int> _slc_mult_pfp, _slc_mult_track, _slc_mult_shower, _slc_mult_track_tolerance;
  std::vector<bool> _slc_geocosmictag;

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

  double _pot;
  */

  /* 
  TTree* _tree2;
  int _total_matches, _nmatch;
  std::vector<double> _hypo_spec, _beam_spec, _fixx_spec;
  double _score;
  int _is_muon;
  */

  TH2F * _deadRegion2P;
  TH2F * _deadRegion3P;
  
  TH1D * _h_pida_proton,     * _h_pida_muon,     * _h_pida_pion,     * _h_pida_kaon;
  TH2D * _h_pida_len_proton, * _h_pida_len_muon, * _h_pida_len_pion, * _h_pida_len_kaon;

  TTree* _sr_tree;
  int _sr_run, _sr_subrun; 
  double _sr_begintime, _sr_endtime;
  double _sr_pot;

  std::ofstream _csvfile;
  std::ofstream _run_subrun_list_file;
};


UBXSec::UBXSec(fhicl::ParameterSet const & p) {

  ::art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;

  _pfp_producer                   = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel                 = p.get<std::string>("HitProducer");
  _geantModuleLabel               = p.get<std::string>("GeantModule");
  _spacepointLabel                = p.get<std::string>("SpacePointProducer");
  _neutrino_flash_match_producer  = p.get<std::string>("NeutrinoFlashMatchProducer");
  _cosmic_flash_match_producer    = p.get<std::string>("CosmicFlashMatchProducer");
  _opflash_producer_beam          = p.get<std::string>("OpFlashBeamProducer");
  _acpt_producer                  = p.get<std::string>("ACPTProducer");
  _tpcobject_producer             = p.get<std::string>("TPCObjectProducer");
  _potsum_producer                = p.get<std::string>("POTSummaryProducer");
  _potsum_instance                = p.get<std::string>("POTSummaryInstance");
  _particle_id_producer           = p.get<std::string>("ParticleIDProducer");
  _mc_ghost_producer              = p.get<std::string>("MCGhostProducer");
  _geocosmictag_producer          = p.get<std::string>("GeoCosmicTaggerProducer");

  _use_genie_info                 = p.get<bool>("UseGENIEInfo", false);
  _minimumHitRequirement          = p.get<int>("MinimumHitRequirement", 3);

  _beam_spill_start               = p.get<double>("BeamSpillStart", 3.2);
  _beam_spill_end                 = p.get<double>("BeamSpillEnd",   4.8);
  _total_pe_cut                   = p.get<double>("TotalPECut",     50);

  _geo_cosmic_score_cut           = p.get<double>("GeoCosmicScoreCut", 0.6);
  _tolerance_track_multiplicity   = p.get<double>("ToleranceTrackMultiplicity", 5.);

  _pecalib.Configure(p.get<fhicl::ParameterSet>("PECalib"));

  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"), 
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");

  _tree1->Branch("ubxsec_event_split", &ubxsec_event, 16000, 99);
 /*
  _tree1->Branch("run",                            &_run,                   "run/I");
  _tree1->Branch("subrun",                         &_subrun,                "subrun/I");
  _tree1->Branch("event",                          &_event,                 "event/I");
  _tree1->Branch("muon_is_reco",                   &_muon_is_reco,          "muon_is_reco/I");
  _tree1->Branch("muon_reco_pur",                  &_muon_reco_pur,         "muon_reco_pur/D");
  _tree1->Branch("muon_reco_eff",                  &_muon_reco_eff,         "muon_reco_eff/D");
  _tree1->Branch("true_muon_mom",                  &_true_muon_mom,         "true_muon_mom/D");
  _tree1->Branch("true_muon_mom_matched",          &_true_muon_mom_matched, "true_muon_mom_matched/D");
  _tree1->Branch("nPFPtagged",                     &_nPFPtagged,            "nPFPtagged/I");
  _tree1->Branch("muon_is_flash_tagged",           &_muon_is_flash_tagged,  "muon_is_flash_tagged/I");
  _tree1->Branch("muon_tag_score",                 &_muon_tag_score,        "muon_tag_score/D");
  _tree1->Branch("fm_score",                       &_fm_score,              "fm_score/D");
  _tree1->Branch("fv",                             &_fv,                    "fv/I");
  _tree1->Branch("ccnc",                           &_ccnc,                  "ccnc/I");
  _tree1->Branch("nupdg",                          &_nupdg,                 "nupdg/I");
  _tree1->Branch("is_signal",                      &_is_signal,             "is_signal/O");
  _tree1->Branch("nu_e",                           &_nu_e,                  "nu_e/D");
  _tree1->Branch("mc_muon_contained",              &_mc_muon_contained,     "mc_muon_contained/O");
  _tree1->Branch("is_swtriggered",                 &_is_swtriggered,        "is_swtriggered/I");
  _tree1->Branch("vtx_resolution",                 &_vtx_resolution,        "vtx_resolution/D");

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
  _tree1->Branch("slc_origin_extra",               "std::vector<int>",    &_slc_origin_extra);
  _tree1->Branch("slc_nhits_u",                    "std::vector<int>",    &_slc_nhits_u);
  _tree1->Branch("slc_nhits_v",                    "std::vector<int>",    &_slc_nhits_v);
  _tree1->Branch("slc_nhits_w",                    "std::vector<int>",    &_slc_nhits_w);
  _tree1->Branch("slc_longesttrack_length",        "std::vector<double>", &_slc_longesttrack_length);
  _tree1->Branch("slc_longesttrack_phi",           "std::vector<double>", &_slc_longesttrack_phi);
  _tree1->Branch("slc_longesttrack_theta",         "std::vector<double>", &_slc_longesttrack_theta);
  _tree1->Branch("slc_longesttrack_iscontained",   "std::vector<bool>",   &_slc_longesttrack_iscontained);
  _tree1->Branch("slc_muoncandidate_exists",       "std::vector<bool>",   &_slc_muoncandidate_exists);
  _tree1->Branch("slc_muoncandidate_length",       "std::vector<double>", &_slc_muoncandidate_length);
  _tree1->Branch("slc_muoncandidate_phi"   ,       "std::vector<double>", &_slc_muoncandidate_phi);
  _tree1->Branch("slc_muoncandidate_theta",        "std::vector<double>", &_slc_muoncandidate_theta);
  _tree1->Branch("slc_acpt_outoftime",             "std::vector<int>",    &_slc_acpt_outoftime);
  _tree1->Branch("slc_crosses_top_boundary",       "std::vector<int>",    &_slc_crosses_top_boundary);
  _tree1->Branch("slc_nuvtx_closetodeadregion_u",  "std::vector<int>",    &_slc_nuvtx_closetodeadregion_u);
  _tree1->Branch("slc_nuvtx_closetodeadregion_v",  "std::vector<int>",    &_slc_nuvtx_closetodeadregion_v);
  _tree1->Branch("slc_nuvtx_closetodeadregion_w",  "std::vector<int>",    &_slc_nuvtx_closetodeadregion_w);
  _tree1->Branch("slc_kalman_chi2",                "std::vector<double>", &_slc_kalman_chi2);
  _tree1->Branch("slc_kalman_ndof",                "std::vector<int>",    &_slc_kalman_ndof);
  _tree1->Branch("slc_passed_min_track_quality",   "std::vector<bool>",   &_slc_passed_min_track_quality);
  _tree1->Branch("slc_passed_min_vertex_quality",  "std::vector<bool>",   &_slc_passed_min_vertex_quality);
  _tree1->Branch("slc_n_intime_pe_closestpmt",     "std::vector<double>", &_slc_n_intime_pe_closestpmt);
  _tree1->Branch("slc_maxdistance_vtxtrack",       "std::vector<double>", &_slc_maxdistance_vtxtrack);
  _tree1->Branch("slc_npfp",                       "std::vector<int>",    &_slc_npfp);
  _tree1->Branch("slc_ntrack",                     "std::vector<int>",    &_slc_ntrack);
  _tree1->Branch("slc_nshower",                    "std::vector<int>",    &_slc_nshower);
  _tree1->Branch("slc_iscontained",                "std::vector<bool>",   &_slc_iscontained);
  _tree1->Branch("slc_mult_pfp",                   "std::vector<int>",    &_slc_mult_pfp);
  _tree1->Branch("slc_mult_track",                 "std::vector<int>",    &_slc_mult_track);
  _tree1->Branch("slc_mult_shower",                "std::vector<int>",    &_slc_mult_shower);
  _tree1->Branch("slc_mult_track_tolerance",       "std::vector<int>",    &_slc_mult_track_tolerance);
  _tree1->Branch("slc_geocosmictag",               "std::vector<bool>",   &_slc_geocosmictag);

  _tree1->Branch("nbeamfls",                       &_nbeamfls,                         "nbeamfls/I");
  _tree1->Branch("beamfls_time",                   "std::vector<double>",              &_beamfls_time);
  _tree1->Branch("beamfls_pe",                     "std::vector<double>",              &_beamfls_pe);
  _tree1->Branch("beamfls_z",                      "std::vector<double>",              &_beamfls_z);
  _tree1->Branch("no_mcflash_but_op_activity",     &_no_mcflash_but_op_activity,       "no_mcflash_but_op_activity/O");
  _tree1->Branch("beamfls_spec",                   "std::vector<std::vector<double>>", &_beamfls_spec);
  _tree1->Branch("numc_flash_spec",                "std::vector<double>",              &_numc_flash_spec);
  _tree1->Branch("slc_flshypo_xfixed_spec",        "std::vector<std::vector<double>>", &_slc_flshypo_xfixed_spec);
  _tree1->Branch("slc_flshypo_spec",               "std::vector<std::vector<double>>", &_slc_flshypo_spec);
  _tree1->Branch("nsignal",                        &_nsignal,                          "nsignal/I");

  _tree1->Branch("mctrk_start_x",                  "std::vector<double>", &_mctrk_start_x);
  _tree1->Branch("mctrk_start_y",                  "std::vector<double>", &_mctrk_start_y);
  _tree1->Branch("mctrk_start_z",                  "std::vector<double>", &_mctrk_start_z);
  _tree1->Branch("trk_start_x",                    "std::vector<double>", &_trk_start_x);
  _tree1->Branch("trk_start_y",                    "std::vector<double>", &_trk_start_y);
  _tree1->Branch("trk_start_z",                    "std::vector<double>", &_trk_start_z);
  _tree1->Branch("vtx_x",                          "std::vector<double>", &_vtx_x);
  _tree1->Branch("vtx_y",                          "std::vector<double>", &_vtx_y);
  _tree1->Branch("vtx_z",                          "std::vector<double>", &_vtx_z);
  _tree1->Branch("tvtx_x",                         "std::vector<double>", &_tvtx_x);
  _tree1->Branch("tvtx_y",                         "std::vector<double>", &_tvtx_y);
  _tree1->Branch("tvtx_z",                         "std::vector<double>", &_tvtx_z);

  _tree1->Branch("pot",                            &_pot,                "pot/D");
  */


  /*
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
  */

  _deadRegion2P = fs->make<TH2F>("deadRegion2P","deadRegion2P", 10350,0.0,1035.0,2300,-115.0,115.0);
  _deadRegion3P = fs->make<TH2F>("deadRegion3P","deadRegion3P", 10350,0.0,1035.0,2300,-115.0,115.0);

  _h_pida_muon = fs->make<TH1D>("h_pida_muon", "Muon tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_proton = fs->make<TH1D>("h_pida_proton", "Proton tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_pion = fs->make<TH1D>("h_pida_pion", "Pion tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_kaon = fs->make<TH1D>("h_pida_kaon", "Kaon tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);

  _h_pida_len_muon = fs->make<TH2D>("h_pida_len_muon", "Muon tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_proton = fs->make<TH2D>("h_pida_len_proton", "Proton tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_pion = fs->make<TH2D>("h_pida_len_pion", "Pion tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_kaon = fs->make<TH2D>("h_pida_len_kaon", "Kaon tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run",                &_sr_run,                "run/I");
  _sr_tree->Branch("subrun",             &_sr_subrun,             "subrun/I");
  _sr_tree->Branch("begintime",          &_sr_begintime,          "begintime/D");
  _sr_tree->Branch("endtime",            &_sr_endtime,            "endtime/D");
  _sr_tree->Branch("pot",                &_sr_pot,                "pot/D");

  _csvfile.open ("pida_trklen.csv", std::ofstream::out | std::ofstream::trunc);
  _csvfile << "pida,trklen,y" << std::endl;

  _run_subrun_list_file.open ("run_subrub_list.txt", std::ofstream::out | std::ofstream::trunc);
}



void UBXSec::produce(art::Event & e) {

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

  /*
  if (_is_data) {
    std::cout << "[UBXSec] Running on a real data file. No MC-PFP matching will be attempted." << std::endl;
  } else {
    mcpfpMatcher.Configure(e, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);
  }
  */

  ::art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;

  // Use '_detp' to find 'efield' and 'temp'
  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield = _detp -> Efield();
  double temp   = _detp -> Temperature();
  // Determine the drift velocity from 'efield' and 'temp'
  _drift_velocity = _detp -> DriftVelocity(efield,temp);
  if (_debug) std::cout << "Using drift velocity = " << _drift_velocity << " cm/us, with E = " << efield << ", and T = " << temp << std::endl;

  // Collect tracks
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, pfParticleToTrackMap);

  // Collect PFParticles and match Reco Particles to Hits
  //lar_pandora::PFParticleVector  recoParticleVector;
  //lar_pandora::PFParticleVector  recoNeutrinoVector;
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  //lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
  //lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    std::cout << "[UBXSec] Cannote locate ubana::TPCObject." << std::endl;
  }
  art::FindManyP<ubana::FlashMatch> tpcobjToFlashMatchAssns(tpcobj_h, e, _neutrino_flash_match_producer);
  art::FindManyP<recob::Track>      tpcobjToTrackAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::PFParticle> tpcobjToPFPAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToCosmicTagAssns(tpcobj_h, e, _geocosmictag_producer);

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

  // Get Tracks
  art::Handle<std::vector<recob::Track>> track_h;
  e.getByLabel(_pfp_producer,track_h);
  if (!track_h.isValid() || track_h->empty()) {
    std::cout << "[UBXSec] Track handle is not valid or empty." << std::endl;
    //throw std::exception();
  }
  std::vector<art::Ptr<recob::Track>> track_p_v;
  art::fill_ptr_vector(track_p_v, track_h);
  art::FindManyP<recob::OpFlash> opfls_ptr_coll_v(track_h, e, _acpt_producer);
  
  // Get PFP
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);
  if(!pfp_h.isValid()){
    std::cout << "[UBXSec] PFP product " << _pfp_producer << " not found..." << std::endl;
    //throw std::exception();
  }
  if(pfp_h->empty()) {
    std::cout << "[UBXSec] PFP " << _pfp_producer << " is empty." << std::endl;
  }
  art::FindManyP<recob::Track> tracks_from_pfp(pfp_h, e, _pfp_producer);

  // Get Giuseppe's Kalman Tracks
  art::FindManyP<recob::Track> trk_kalman_v(pfp_h, e, "pandoraNuKalmanTrack");

  // Get PID information
  art::FindMany<anab::ParticleID> particleids_from_track (track_h, e, _particle_id_producer);
  if (!particleids_from_track.isValid()) {
    std::cout << "[UBXSec] anab::ParticleID is not valid." << std::endl;
  }
  std::cout << "[UBXSec] Numeber of particleids_from_track " << particleids_from_track.size() << std::endl;

  // Fill a std::map Track->ParticleID
  std::map<recob::Track, anab::ParticleID> track_to_pid_map;
  for (auto track : track_p_v) {
    std::vector<const anab::ParticleID*> pids = particleids_from_track.at(track.key());

    if(pids.size() == 0) 
      continue;

    for (auto pid : pids) {
      if (!pid->PlaneID().isValid) continue;
      int planenum = pid->PlaneID().Plane;
      if (planenum != 2) continue;
      track_to_pid_map[(*track)] = (*pid);
      continue;
    }
  }

  // Get Ghosts
  art::Handle<std::vector<ubana::MCGhost> > ghost_h;
  e.getByLabel(_mc_ghost_producer,ghost_h);
  if(!ghost_h.isValid()){
    std::cout << "[UBXSec] MCGhost product " << _mc_ghost_producer << " not found..." << std::endl;
    //throw std::exception();
  }
  art::FindManyP<ubana::MCGhost>   mcghost_from_pfp   (pfp_h,   e, _mc_ghost_producer);
  art::FindManyP<simb::MCParticle> mcpar_from_mcghost (ghost_h, e, _mc_ghost_producer); 

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
    if (_debug) std::cout << "[UBXSec] Reco beam flash pe: " << std::endl;
    for (unsigned int i = 0; i < 32; i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      _beamfls_spec[n][opdet] = flash.PE(i);
      if (_beamfls_time[n] > _beam_spill_start && _beamfls_time[n] < _beam_spill_end) {
        if (flash.TotalPE() > _total_pe_cut) 
          _candidate_flash_time = flash.Time();
        if (_debug) std::cout << "\t PMT " << opdet << ": " << _beamfls_spec[n][opdet] << std::endl;
      }
    }
  }

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

    if (mclist[iList]->Origin() == NEUTRINO_ORIGIN) {

      if (_debug) this->PrintMC(mclist); 

      double truth_nu_vtx[3] = {mclist[iList]->GetNeutrino().Nu().Vx(),
                                mclist[iList]->GetNeutrino().Nu().Vy(),
                                mclist[iList]->GetNeutrino().Nu().Vz()};
      if (_fiducial_volume.InFV(truth_nu_vtx)) _fv = 1;
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
    } // neutrino origin
    else {
      _ccnc = -1;
      _nupdg = -1;
      _nu_e = -1;
      _true_muon_mom = -9999.;
    }
  } else {
    _ccnc = -1;
    _nupdg = -1;
    _nu_e = -1;
    _true_muon_mom = -9999.;
  }

  _is_signal = false;
  if (_ccnc == 0 && _nupdg == 14 && _fv == 1) {
    _is_signal = true;
  }

  lar_pandora::PFParticlesToSpacePoints pfp_to_spacept;
  lar_pandora::SpacePointsToHits spacept_to_hits;

  lar_pandora::PFParticleVector temp2;
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, temp2, pfp_to_spacept);

  lar_pandora::SpacePointVector temp3;
  lar_pandora::LArPandoraHelper::CollectSpacePoints(e, _pfp_producer, temp3, spacept_to_hits);

  art::Handle<std::vector<recob::OpHit>> ophit_h;
  e.getByLabel("ophitBeam", ophit_h);
  if(!ophit_h.isValid()) {
    std::cout << "[UBXSec] Cannot locate OpHits." << std::endl;
  }

 
  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;
  for (size_t slice = 0; slice < tpcobj_h->size(); slice++) {
    track_v_v.push_back(tpcobjToTrackAssns.at(slice));
    pfp_v_v.push_back(tpcobjToPFPAssns.at(slice));
  }


  _nslices = tpcobj_h->size();
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
  _slc_origin_extra.resize(_nslices);
  _slc_flshypo_xfixed_spec.resize(_nslices);
  _slc_flshypo_spec.resize(_nslices);
  _slc_nhits_u.resize(_nslices, -9999);
  _slc_nhits_v.resize(_nslices, -9999);
  _slc_nhits_w.resize(_nslices, -9999);
  _slc_flsmatch_cosmic_score.resize(_nslices, -9999);
  _slc_flsmatch_cosmic_t0.resize(_nslices, -9999);
  _slc_longesttrack_length.resize(_nslices, -9999);
  _slc_longesttrack_phi.resize(_nslices, -9999);
  _slc_longesttrack_theta.resize(_nslices, -9999);
  _slc_longesttrack_iscontained.resize(_nslices, -9999);
  _slc_muoncandidate_exists.resize(_nslices, -9999);
  _slc_muoncandidate_length.resize(_nslices, -9999);
  _slc_muoncandidate_phi.resize(_nslices, -9999);
  _slc_muoncandidate_theta.resize(_nslices, -9999);
  _slc_acpt_outoftime.resize(_nslices, -9999);
  _slc_crosses_top_boundary.resize(_nslices, -9999);
  _slc_nuvtx_closetodeadregion_u.resize(_nslices, -9999);
  _slc_nuvtx_closetodeadregion_v.resize(_nslices, -9999);
  _slc_nuvtx_closetodeadregion_w.resize(_nslices, -9999);  
  _slc_kalman_chi2.resize(_nslices, -9999);
  _slc_kalman_ndof.resize(_nslices, -9999);  
  _slc_passed_min_track_quality.resize(_nslices, -9999);
  _slc_passed_min_vertex_quality.resize(_nslices, -9999);
  _slc_n_intime_pe_closestpmt.resize(_nslices, -9999);
  _slc_maxdistance_vtxtrack.resize(_nslices, -9999);
  _slc_npfp.resize(_nslices, -9999);
  _slc_ntrack.resize(_nslices, -9999);
  _slc_nshower.resize(_nslices, -9999);
  _slc_iscontained.resize(_nslices, -9999);
  _slc_mult_pfp.resize(_nslices, -9999);
  _slc_mult_track.resize(_nslices, -9999);
  _slc_mult_shower.resize(_nslices, -9999);
  _slc_mult_track_tolerance.resize(_nslices, -9999);
  _slc_geocosmictag.resize(_nslices, false);

  std::cout << "[UBXSec] --- SAVING INFORMATION" << std::endl;
  _vtx_resolution = -9999;
 
  for (unsigned int slice = 0; slice < tpcobj_h->size(); slice++){
    std::cout << "[UBXSec] >>> SLICE " << slice << std::endl;

    ubana::TPCObject tpcobj = (*tpcobj_h)[slice];

    _slc_npfp[slice]    = tpcobj.GetNPFP();
    _slc_ntrack[slice]  = tpcobj.GetNTracks();
    _slc_nshower[slice] = tpcobj.GetNShowers();

    // Slice origin 
    _slc_origin[slice] = tpcobj.GetOrigin();
    std::cout << "[UBXSec] \t Origin is " << _slc_origin[slice] << std::endl;

    // Slice origin extra
    _slc_origin_extra[slice] = tpcobj.GetOriginExtra();
    std::cout << "[UBXSec] \t Origin extra is " << _slc_origin_extra[slice] << std::endl;

    // Containment
    _slc_iscontained[slice] = UBXSecHelper::TracksAreContained(tpcobj.GetTracks());

    // Reco vertex
    double reco_nu_vtx_raw[3];
    recob::Vertex tpcobj_nu_vtx = tpcobj.GetVertex();
    tpcobj_nu_vtx.XYZ(reco_nu_vtx_raw);

    // X position correction
    double reco_nu_vtx[3];
    UBXSecHelper::GetTimeCorrectedPoint(reco_nu_vtx_raw, reco_nu_vtx, _candidate_flash_time, _drift_velocity);
    _slc_nuvtx_x[slice] = reco_nu_vtx[0];
    _slc_nuvtx_y[slice] = reco_nu_vtx[1];
    _slc_nuvtx_z[slice] = reco_nu_vtx[2];
    _slc_nuvtx_fv[slice] = (_fiducial_volume.InFV(reco_nu_vtx) ? 1 : 0);
    std::cout << "[UBXSec] \t Reco vertex is " << _slc_nuvtx_x[slice] << ", " << _slc_nuvtx_y[slice] << ", " << _slc_nuvtx_z[slice] << std::endl; 
    std::cout << "[UBXSec] \t Reco vertex is " << (_slc_nuvtx_fv[slice]==1 ? "in" : "ouside") << " the FV." << std::endl;

    // Through-going?
    _slc_geocosmictag[slice] = false;
    std::vector<art::Ptr<anab::CosmicTag>> geo_cosmic_tags = tpcobjToCosmicTagAssns.at(slice);
    if(geo_cosmic_tags.size() == 0 || geo_cosmic_tags.size() > 1) {
      std::cout << "[UBXSec] \t More than one Geo Cosmic Tag match per tpcobj ?!" << std::endl;
    } else {
      auto ct = geo_cosmic_tags.at(0);
      if (ct->CosmicScore() > _geo_cosmic_score_cut) {
        std::cout << "[UBXSec] \t This slice has been tagged as through-going cosmic" << std::endl;
        _slc_geocosmictag[slice] = true;
      }
    }
    
    // Vertex resolution
    if (_slc_origin[slice] == ubana::kBeamNeutrino) {
      _vtx_resolution = sqrt( pow(_slc_nuvtx_y[slice]-_tvtx_y[0], 2) + pow(_slc_nuvtx_z[slice]-_tvtx_z[0], 2) );
    } 

    // Multiplicity
    int p, t, s;
    tpcobj.GetMultiplicity(p, t, s);
    _slc_mult_pfp[slice] = p;
    _slc_mult_track[slice] = t;
    _slc_mult_shower[slice] = s;
    _slc_mult_track_tolerance[slice] = tpcobj.GetNTracksCloseToVertex(_tolerance_track_multiplicity);

    // Neutrino Flash match
    _slc_flsmatch_score[slice] = -9999;
    std::vector<art::Ptr<ubana::FlashMatch>> pfpToFlashMatch_v = tpcobjToFlashMatchAssns.at(slice);
    if (pfpToFlashMatch_v.size() > 1) {
      std::cout << "[UBXSec] \t More than one flash match per nu pfp ?!" << std::endl;
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
      //for (auto v : _slc_flshypo_spec[slice]) std::cout << "Hypo PE: " << v << std::endl;

      std::cout << "[UBXSec] \t FM score:       " << _slc_flsmatch_score[slice] << std::endl;
      std::cout << "[UBXSec] \t qllx - tpcx is: " << _slc_flsmatch_qllx[slice] - _slc_flsmatch_tpcx[slice] << std::endl;
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
      _slc_longesttrack_phi[slice]   = UBXSecHelper::GetCorrectedPhi(lt, tpcobj_nu_vtx);
      _slc_longesttrack_theta[slice] = UBXSecHelper::GetCorrectedCosTheta(lt, tpcobj_nu_vtx);
      _slc_longesttrack_iscontained[slice] = UBXSecHelper::TrackIsContained(lt);
      int vtx_ok;
      _slc_crosses_top_boundary[slice] = (UBXSecHelper::IsCrossingTopBoundary(lt, vtx_ok) ? 1 : 0);
    } else {
      _slc_longesttrack_length[slice] = -9999;
    }

    // ACPT
    _slc_acpt_outoftime[slice] = 0;
    for (unsigned int t = 0; t < track_v_v[slice].size(); t++) {
      if(opfls_ptr_coll_v.at(track_v_v[slice][t].key()).size()>1) {
        std::cout << "[UBXSec] \t More than 1 association found (ACPT)!" << std::endl;
        //throw std::exception();
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
        std::cout << "[UBXSec] \t TQ more than one track per PFP, ntracks " << trk_kalman_v.at(pfp_v_v[slice][t].key()).size() << std::endl;
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

    // Vertex quality
    recob::Vertex slice_vtx = tpcobj.GetVertex();
    double slice_vtx_xyz[3];
    slice_vtx.XYZ(slice_vtx_xyz);
    _slc_passed_min_vertex_quality[slice] = true;
    if (deadRegionsFinder.NearDeadReg2P(slice_vtx_xyz[1], slice_vtx_xyz[2], 5.0))
      _slc_passed_min_vertex_quality[slice] = false;

    // Channel status
    _slc_nuvtx_closetodeadregion_u[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 0) ? 1 : 0);
    _slc_nuvtx_closetodeadregion_v[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 1) ? 1 : 0);
    _slc_nuvtx_closetodeadregion_w[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 2) ? 1 : 0);

    // Vertex check
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
        std::cout << "[UBXSec] \t Can't find spacepoints for pfp with pdg " << pfp->PdgCode() << std::endl;
        continue;
      }
      // Loop through the hits associated 
      for (auto sp_pt : (iter->second)) {
        auto iter2 = spacept_to_hits.find(sp_pt);
        if (iter2 != spacept_to_hits.end()) {
          //std::cout << "[UBXSec] Founds hits associated to this sp_pt" << std::endl;
        } else {
          std::cout << "[UBXSec] \t Can't find hits ass to this sp_pt" << std::endl;
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
    std::cout << "[UBXSec] \t For this TPC object we have " << hit3d_v.size() << " Hit3D_t hits." << std::endl;

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

    // Muon Candidate
    ubana::MuonCandidateFinder muon_finder;
    muon_finder.SetTPCObject(tpcobj);
    muon_finder.SetTrackToPIDMap(track_to_pid_map);
    recob::Track candidate_track;

    if (muon_finder.GetCandidateTrack(candidate_track)) {
      _slc_muoncandidate_exists[slice] = true;
      _slc_muoncandidate_length[slice] = candidate_track.Length();
      _slc_muoncandidate_phi[slice]    = UBXSecHelper::GetCorrectedPhi(candidate_track, tpcobj_nu_vtx); 
      _slc_muoncandidate_theta[slice]  = UBXSecHelper::GetCorrectedCosTheta(candidate_track, tpcobj_nu_vtx);
    } else {
      _slc_muoncandidate_exists[slice] = false;
      _slc_muoncandidate_length[slice] = -9999;
      _slc_muoncandidate_phi[slice]    = -9999;
      _slc_muoncandidate_theta[slice]  = -9999;
    }

    // Particle ID
    auto pfps_from_tpcobj = tpcobjToPFPAssns.at(slice);

    for (auto pfp : pfps_from_tpcobj){

      std::cout << "[UBXSec] \t This is PFP " << pfp->Self()  << std::endl;
      std::vector<art::Ptr<ubana::MCGhost>> mcghosts = mcghost_from_pfp.at(pfp.key());
      std::vector<art::Ptr<simb::MCParticle>> mcpars;
      int pdg = -1;
      if (mcghosts.size() == 0 || mcghosts.size() > 1 ) {
        std::cout << "[UBXSec] \t\t mcghosts is ether 0 or >1" << std::endl;
        continue;
      }

      mcpars = mcpar_from_mcghost.at(mcghosts[0].key());
      pdg = mcpars[0]->PdgCode();
      std::cout << "[UBXSec] \t\t MCPar has pdg " << pdg << std::endl;
      const auto mc_truth = bt->TrackIDToMCTruth(mcpars[0]->TrackId()); 
      if (!mc_truth) {
        std::cerr << "[UBXSec] Problem with MCTruth pointer." << std::endl;
        continue;
      }
      if (mc_truth->Origin() == simb::kBeamNeutrino &&
          mcpars[0]->PdgCode() == 13 && mcpars[0]->Mother() == 0) {

        _muon_reco_pur = _muon_reco_eff = -9999;
        auto iter = recoParticlesToHits.find(pfp);
        if (iter != recoParticlesToHits.end()) {
          UBXSecHelper::GetTrackPurityAndEfficiency((*iter).second, _muon_reco_pur, _muon_reco_eff);
        }
        _true_muon_mom_matched = mcpars[0]->P();

       if (_fiducial_volume.InFV(mcpars[0]->Vx(), mcpars[0]->Vy(), mcpars[0]->Vz())  && _fiducial_volume.InFV(mcpars[0]->EndX(), mcpars[0]->EndY(), mcpars[0]->EndZ())) {
         _mc_muon_contained = true;
       }
      } 
      

      std::vector<art::Ptr<recob::Track>> tracks = tracks_from_pfp.at(pfp.key());
      std::cout << "[UBXSec] \t\t n tracks ass to this pfp: " << tracks.size() << std::endl;
      for (auto track : tracks) {

        std::vector<const anab::ParticleID*> pids = particleids_from_track.at(track.key());
        if(pids.size() == 0) std::cout << "[UBXSec] \t\t Zero ParticleID" << std::endl;
        if(pids.size() > 1) {
          std::cout << "[UBXSec] \t\t ParticleID vector is bigger than 1. Only one saved." << std::endl;
        }
        for (auto pid : pids) {
          if (!pid->PlaneID().isValid) continue;
          int planenum = pid->PlaneID().Plane;
          if (planenum < 0 || planenum > 2) continue;
          std::cout << "[UBXSec] \t\t ParticleID PIDA is " << pid->PIDA() << ", plane is " << planenum << std::endl;
          if (/*_is_signal && (_slc_origin[slice] == 0 || _slc_origin[slice] == 2) &&*/ planenum == 2) {
            if (pdg == 13) {
              _h_pida_muon->Fill(pid->PIDA());
              _h_pida_len_muon->Fill(pid->PIDA(), track->Length());
              if( pid->PIDA() > 0 && pid->PIDA() < 50. ) _csvfile << pid->PIDA() << "," << track->Length() << "," << "1" << std::endl;
            } else if (pdg == 2212) {
              _h_pida_proton->Fill(pid->PIDA());
              _h_pida_len_proton->Fill(pid->PIDA(), track->Length());
              if( pid->PIDA() > 0 && pid->PIDA() < 50. ) _csvfile << pid->PIDA() << "," << track->Length() << "," << "0" << std::endl;
            } else if (pdg == 211) {
              _h_pida_pion->Fill(pid->PIDA());
              _h_pida_len_pion->Fill(pid->PIDA(), track->Length());
            } else if (pdg == 321) {
              _h_pida_kaon->Fill(pid->PIDA());
              _h_pida_len_kaon->Fill(pid->PIDA(), track->Length());
            }
          }
        }
      }
    }
  
    std::cout << "[UBXSec] --- SLICE INFORMATION SAVED" << std::endl;
  } // slice loop



  /* Dead regions
  //art::ServiceHandle<FindDeadRegions> deadRegionsFinder;
  FindDeadRegions deadRegionsFinder;
  deadRegionsFinder.GetDeadRegionHisto2P(_deadRegion2P);
  deadRegionsFinder.GetDeadRegionHisto3P(_deadRegion3P);
  */




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


  // POT

  art::Handle< sumdata::POTSummary > potsum_h;
  if(e.getByLabel(_potsum_producer, potsum_h))
    _pot = potsum_h->totpot;
  else
    _pot = 0.;

 


  if(_debug) std::cout << "[UBXSec] Filling tree now." << std::endl;
  _tree1->Fill();

  if(_debug) std::cout << "********** UBXSec ends" << std::endl;

  return;
}



void UBXSec::endSubRun(art::SubRun& sr) {

  if (_debug) std::cout << "[UBXSec::endSubRun] Starts" << std::endl;

  // Saving run and subrun number on file so that we can run Zarko's script easily
  _run_subrun_list_file << sr.run() << " " << sr.subRun() << std::endl;

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> potsum_h;

  // MC
  if (_is_mc) {
    if (_debug) std::cout << "[UBXSec::endSubRun] Getting POT for MC" << std::endl;
    if(sr.getByLabel(_potsum_producer, potsum_h)) {
      if (_debug) std::cout << "[UBXSec::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot;
    }
    else
      _sr_pot = 0.;
  }

  // Data
  if (_is_data) {
    if (_debug) std::cout << "[UBXSec::endSubRun] Getting POT for DATA, producer " << _potsum_producer << ", instance " << _potsum_instance << std::endl;
    if (sr.getByLabel(_potsum_producer, _potsum_instance, potsum_h)){
      if (_debug) std::cout << "[UBXSec::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot;
    }
    else
      _sr_pot = 0;
  }

  _sr_tree->Fill();

  if (_debug) std::cout << "[UBXSec::endSubRun] Ends" << std::endl;
}

void UBXSec::PrintMC(std::vector<art::Ptr<simb::MCTruth>> mclist) {

  std::cout << "[UBXSec] ================= MC Information ================= [UBXSec]" << std::endl;

  int iList = 0;
  std::cout << " NEUTRINO:" << std::endl;
  if (mclist[iList]->NeutrinoSet()) {
    std::cout << "\tPDG      " << mclist[iList]->GetNeutrino().Nu().PdgCode() << std::endl;
    std::cout << "\tCC/NC?   " << (mclist[iList]->GetNeutrino().CCNC() == 0 ? "CC" : "NC") << std::endl;
    std::cout << "\tMode     " << mclist[iList]->GetNeutrino().Mode() << std::endl;
    std::cout << "\tQSqr     " << mclist[iList]->GetNeutrino().QSqr() << std::endl;
    std::cout << "\tW        " << mclist[iList]->GetNeutrino().W() << std::endl;
    std::cout << "\tX        " << mclist[iList]->GetNeutrino().X() << std::endl;
    std::cout << "\tY        " << mclist[iList]->GetNeutrino().Y() << std::endl;
    std::cout << "\tHitNuc   " << mclist[iList]->GetNeutrino().HitNuc() << std::endl;
    std::cout << "\tE        " << mclist[iList]->GetNeutrino().Nu().E() << std::endl;
    std::cout << "\tVx       " << mclist[iList]->GetNeutrino().Nu().Vx() << std::endl;
    std::cout << "\tVy       " << mclist[iList]->GetNeutrino().Nu().Vy() << std::endl;
    std::cout << "\tVz       " << mclist[iList]->GetNeutrino().Nu().Vz() << std::endl;
  } else
    std::cout << "\t---No Neutrino information---" << std::endl;

  std::cout << std::endl;
  std::cout << " PRIMARIES (only with status code==1):" << std::endl;
  for (int p = 0; p < mclist[0]->NParticles(); p++) {
    const simb::MCParticle mc_par = mclist[0]->GetParticle(p);
    if (mc_par.StatusCode() != 1) continue;
    std::cout << "\tPDG           " << mc_par.PdgCode() << std::endl;
    std::cout << "\tStart process " << mc_par.Process() << std::endl;
    std::cout << "\tEnd process   " << mc_par.EndProcess() << std::endl;
    std::cout << "\tEnergy        " << mc_par.E() << std::endl;
    std::cout << "\tMomentum      " << mc_par.P() << std::endl;
    std::cout << "\tVertex        " << mc_par.Vx() << ", " << mc_par.Vy() << ", " << mc_par.Vz() << std::endl;
    std::cout << "\tStatus Code   " << mc_par.StatusCode() << std::endl << std::endl;
  }

  std::cout << "[UBXSec] ================= MC Information ================= [UBXSec]" << std::endl;
}

DEFINE_ART_MODULE(UBXSec)
