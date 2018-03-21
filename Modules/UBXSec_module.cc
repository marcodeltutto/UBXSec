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
 * \author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 *
 * \version producer 
 *
 * \date 2017/03/10
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
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "uboone/UBXSec/DataTypes/FlashMatch.h"
#include "uboone/UBXSec/DataTypes/MCGhost.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfo/DetectorClocks.h"
#include "uboone/EventWeight/MCEventWeight.h"

#include "uboone/UBXSec/DataTypes/UBXSecEvent.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"


#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


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
#include "uboone/UBXSec/Algorithms/NuMuCCEventSelection.h"
#include "uboone/UBXSec/Algorithms/TrackQuality.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// Root include
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

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
  ::ubana::MuonCandidateFinder _muon_finder;
  ::ubana::NuMuCCEventSelection _event_selection;
  ::pmtana::PECalib _pecalib;
  ::trkf::TrackMomentumCalculator _trk_mom_calculator;

  // Database to understand particle pdg
  const TDatabasePDG* _database_pdg = TDatabasePDG::Instance();

  // Detector info service
  ::detinfo::DetectorProperties const* fDetectorProperties;

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
  std::string _candidateconsistency_producer;
  std::string _mcsfitresult_mu_producer;
  std::string _mcsfitresult_pi_producer;
  std::string _calorimetry_producer;
  std::string _eventweight_producer;
  std::string _genie_eventweight_pm1_producer;
  std::string _genie_eventweight_multisim_producer;
  std::string _flux_eventweight_multisim_producer;
  bool _debug = true;                   ///< Debug mode
  int _minimumHitRequirement;           ///< Minimum number of hits in at least a plane for a track
  double _minimumDistDeadReg;           ///< Minimum distance the track end points can have to a dead region
  bool _use_genie_info;                 ///< Turn this off if looking at cosmic only files
  double _beam_spill_start;             ///< Start time of the beam spill (us) 
  double _beam_spill_end;               ///< Start time of the beam spill (us) 
  double _total_pe_cut;                 ///< PE cut to be applied to beam flash
  double _geo_cosmic_score_cut;         ///< Cut on the score of the pandoraNu geo cosmic tagger
  double _tolerance_track_multiplicity; ///< Tolerance to consider a track coming from the nu reco vertex
  double _min_track_len;                ///< Min track length for momentum calculation
  bool _make_ophit_csv;                 ///< If true makea a csv file with ophit info
  bool _make_pida_csv;                  ///< If true makea a csv file with pida/tracklength info

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
  TTree* _tree2;
  int _total_matches, _nmatch;
  std::vector<double> _hypo_spec, _beam_spec, _fixx_spec;
  double _score;
  int _is_muon;
  */

  TH2F * _deadRegion2P;
  TH2F * _deadRegion3P;
  
  // PIDA related variables
  TH1D * _h_pida_proton,     * _h_pida_muon,     * _h_pida_pion,     * _h_pida_kaon;
  TH2D * _h_pida_len_proton, * _h_pida_len_muon, * _h_pida_len_pion, * _h_pida_len_kaon;

  // Momentum Related Variables
  TH2D* _h_mom_true_mcs; ///< 2D histogram of true muon momentum VS reconstructed (using MCS)
  TH2D* _h_mom_true_mcs_contained; ///< 2D histogram of true muon momentum VS reconstructed (using MCS) (contained tracks)
  TH2D* _h_mom_true_mcs_uncontained; ///< 2D histogram of true muon momentum VS reconstructed (using MCS) (uncontained tracks)
  TH2D* _h_mom_true_range_contained; ///< 2D histogram of true muon momentum VS reconstructed (using Length) (contained tracks)
  TH2D* _h_mom_range_mcs_contained;  ///< 2D histogram of reconstructed (using MCS) muon momentum VS reconstructed (using Length) (contained tracks)
  TH1D* _h_mcs_cosmic_track_direction; ///< Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward)
  TH2D* _h_mcs_cosmic_track_direction_deltall; /// Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward) VS delta LL from MCS fit
  TH2D* _h_mcs_cosmic_track_direction_ratioll; /// Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward) VS ratio LL from MCS fit
  TTree *_mom_tree_contained, *_mom_tree_uncontained;
  int _run, _subrun, _event;
  double _mom_true_contained;
  double _mom_mcs_contained;
  double _mom_range_contained;
  double _mom_true_uncontained;
  double _mom_mcs_uncontained;
  TTree *_mcs_cosmic_track_direction_tree;
  double _mcs_cosmic_track_direction;
  double _mcs_cosmic_track_downll;
  double _mcs_cosmic_track_upll;
  TTree *_mom_cosmic_tree;
  double _mom_cosmic_true;
  double _mom_cosmic_mcs;
  double _mom_cosmic_mcs_downforced;
  double _mom_cosmic_range;
  bool _mom_cosmic_down;

  TTree* _sr_tree;
  int _sr_run, _sr_subrun; 
  double _sr_begintime, _sr_endtime;
  double _sr_pot;

  std::ofstream _csvfile, _csvfile2;
  std::ofstream _run_subrun_list_file;
};


UBXSec::UBXSec(fhicl::ParameterSet const & p) {

  //::art::ServiceHandle<cheat::BackTracker> bt;
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
  _candidateconsistency_producer  = p.get<std::string>("CandidateConsistencyProducer");
  _mcsfitresult_mu_producer       = p.get<std::string>("MCSFitResultMuProducer");
  _mcsfitresult_pi_producer       = p.get<std::string>("MCSFitResultPiProducer");
  _calorimetry_producer           = p.get<std::string>("CalorimetryProducer");
  _eventweight_producer           = p.get<std::string>("EventWeightProducer");
  _genie_eventweight_pm1_producer = p.get<std::string>("GenieEventWeightPMOneProducer");
  _genie_eventweight_multisim_producer = p.get<std::string>("GenieEventWeightMultisimProducer");
  _flux_eventweight_multisim_producer = p.get<std::string>("FluxEventWeightMultisimProducer");

  _use_genie_info                 = p.get<bool>("UseGENIEInfo", false);
  _minimumHitRequirement          = p.get<int>("MinimumHitRequirement", 3);
  _minimumDistDeadReg             = p.get<double>("MinimumDistanceToDeadRegion", 5.);

  _beam_spill_start               = p.get<double>("BeamSpillStart", 3.2);
  _beam_spill_end                 = p.get<double>("BeamSpillEnd",   4.8);
  _total_pe_cut                   = p.get<double>("TotalPECut",     50);

  _geo_cosmic_score_cut           = p.get<double>("GeoCosmicScoreCut", 0.6);
  _tolerance_track_multiplicity   = p.get<double>("ToleranceTrackMultiplicity", 5.);

  _min_track_len                  = p.get<double>("MinTrackLength", 0.1);

  _make_ophit_csv                 = p.get<bool>("MakeOpHitCSV", false);
  _make_pida_csv                  = p.get<bool>("MakePIDACSV", false);

  _pecalib.Configure(p.get<fhicl::ParameterSet>("PECalib"));

  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"), 
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();

  _muon_finder.Configure(p.get<fhicl::ParameterSet>("MuonCandidateFinderSettings"));

  _muon_finder.PrintConfig();

  _event_selection.Configure(p.get<fhicl::ParameterSet>("NuMuCCSelectionSettings"));

  _event_selection.PrintConfig();

  _trk_mom_calculator.SetMinLength(_min_track_len);

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
/*
  auto index = gROOT->GetListOfSpecials()->IndexOf(_database_pdg);
  if (index == -1)
    _database_pdg = new TDatabasePDG;
  else 
    _database_pdg = (TDatabasePDG*)gROOT->GetListOfSpecials()->At(index);

  
  if (!TDatabasePDG::fgInstance)
    _database_pdg = TDatabasePDG::fgInstance;
  else 
    _database_pdg = new TDatabasePDG;
*/

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");

  int bufsize    = 16000;
  int splitlevel = 99;
  _tree1->Branch("ubxsec_event_split", &ubxsec_event, bufsize, splitlevel);

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

  _h_mom_true_mcs = fs->make<TH2D>("h_mom_true_mcs", ";True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_mcs_contained = fs->make<TH2D>("h_mom_true_mcs_contained", "Contained;True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_mcs_uncontained = fs->make<TH2D>("h_mom_true_mcs_uncontained", "Uncontained;True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_range_contained = fs->make<TH2D>("h_mom_true_range_contained", "Contained;True Muon Momentum [GeV];Reconstructed (via Length) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_range_mcs_contained = fs->make<TH2D>("h_mom_range_mcs_contained", "Contained;Reconstructed (via Length) Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);

  _h_mcs_cosmic_track_direction = fs->make<TH1D>("h_mcs_cosmic_track_direction", "0: down, 1: up;;", 2, 0, 2);
  _h_mcs_cosmic_track_direction_deltall = fs->make<TH2D>("h_mcs_cosmic_track_direction_deltall", ";0: down, 1: up;(FWD - BWD) LL;", 2, 0, 2, 500, -1e-6, 1e-6);
  _h_mcs_cosmic_track_direction_ratioll = fs->make<TH2D>("h_mcs_cosmic_track_direction_ratioll", ";0: down, 1: up;(FWD / BWD) LL;", 2, 0, 2, 500, 0, 2);

  _mom_tree_contained = fs->make<TTree>("mom_tree_contained","");
  _mom_tree_contained->Branch("run",                 &_run,                 "run/I");
  _mom_tree_contained->Branch("subrun",              &_subrun,              "subrun/I");
  _mom_tree_contained->Branch("event",               &_event,               "event/I");
  _mom_tree_contained->Branch("mom_true_contained",  &_mom_true_contained,  "mom_true_contained/D");
  _mom_tree_contained->Branch("mom_mcs_contained",   &_mom_mcs_contained,   "mom_mcs_contained/D");
  _mom_tree_contained->Branch("mom_range_contained", &_mom_range_contained, "mom_range_contained/D");

  _mom_tree_uncontained = fs->make<TTree>("mom_tree_uncontained","");
  _mom_tree_uncontained->Branch("run",                  &_run,                   "run/I");
  _mom_tree_uncontained->Branch("subrun",               &_subrun,                "subrun/I");
  _mom_tree_uncontained->Branch("event",                &_event,                 "event/I");
  _mom_tree_uncontained->Branch("mom_true_uncontained", &_mom_true_uncontained,  "mom_true_uncontained/D");
  _mom_tree_uncontained->Branch("mom_mcs_uncontained",  &_mom_mcs_uncontained,   "mom_mcs_uncontained/D");

  _mcs_cosmic_track_direction_tree = fs->make<TTree>("mcs_cosmic_track_direction_tree","");
  _mcs_cosmic_track_direction_tree->Branch("run",                 &_run,                 "run/I");
  _mcs_cosmic_track_direction_tree->Branch("subrun",              &_subrun,              "subrun/I");
  _mcs_cosmic_track_direction_tree->Branch("event",               &_event,               "event/I");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_direction",  &_mcs_cosmic_track_direction,  "mcs_cosmic_track_direction/D");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_downll",     &_mcs_cosmic_track_downll,     "mcs_cosmic_track_downll/D");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_upll",       &_mcs_cosmic_track_upll,       "mcs_cosmic_track_upll/D");

  _mom_cosmic_tree = fs->make<TTree>("mom_cosmic_tree","");
  _mom_cosmic_tree->Branch("run",                       &_run,                       "run/I");
  _mom_cosmic_tree->Branch("subrun",                    &_subrun,                    "subrun/I");
  _mom_cosmic_tree->Branch("event",                     &_event,                     "event/I");
  _mom_cosmic_tree->Branch("mom_cosmic_true",           &_mom_cosmic_true,           "mom_cosmic_true/D");
  _mom_cosmic_tree->Branch("mom_cosmic_mcs",            &_mom_cosmic_mcs,            "mom_cosmic_mcs/D");
  _mom_cosmic_tree->Branch("mom_cosmic_mcs_downforced", &_mom_cosmic_mcs_downforced, "mom_cosmic_mcs_downforced/D");
  _mom_cosmic_tree->Branch("mom_cosmic_range",          &_mom_cosmic_range,          "mom_cosmic_range/D");
  _mom_cosmic_tree->Branch("mom_cosmic_down",           &_mom_cosmic_down,           "mom_cosmic_down/O");

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run",                &_sr_run,                "run/I");
  _sr_tree->Branch("subrun",             &_sr_subrun,             "subrun/I");
  _sr_tree->Branch("begintime",          &_sr_begintime,          "begintime/D");
  _sr_tree->Branch("endtime",            &_sr_endtime,            "endtime/D");
  _sr_tree->Branch("pot",                &_sr_pot,                "pot/D");

  if(_make_pida_csv) _csvfile.open ("pida_trklen.csv", std::ofstream::out | std::ofstream::trunc);
  if(_make_pida_csv) _csvfile << "pida,trklen,y" << std::endl;

  if(_make_ophit_csv) _csvfile2.open("ophit.csv", std::ofstream::out | std::ofstream::trunc);
  if(_make_ophit_csv) _csvfile2 << "ophit,opdet,time,pe" << std::endl;

  _run_subrun_list_file.open ("run_subrub_list.txt", std::ofstream::out | std::ofstream::trunc);




  produces<std::vector<ubana::SelectionResult>>();
  produces<art::Assns<ubana::SelectionResult, ubana::TPCObject>>();

  // For the neutrino id filter
  produces<art::Assns<recob::Vertex, recob::Track>>();
  produces< art::Assns<recob::Vertex, recob::PFParticle>>();

}



void UBXSec::produce(art::Event & e) {

  if(_debug) std::cout << "********** UBXSec starts" << std::endl;
  if(_debug) std::cout << "Run: " << e.id().run() << 
                          ", subRun: " << e.id().subRun() <<
                          ", event: " << e.id().event()  << std::endl;


  // Instantiate the output
  std::unique_ptr< std::vector<ubana::SelectionResult>>                   selectionResultVector           (new std::vector<ubana::SelectionResult>);
  std::unique_ptr< art::Assns<ubana::SelectionResult, ubana::TPCObject>>  assnOutSelectionResultTPCObject (new art::Assns<ubana::SelectionResult, ubana::TPCObject>);

  std::unique_ptr<art::Assns<recob::Vertex, recob::Track>>      vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
  std::unique_ptr<art::Assns<recob::Vertex, recob::PFParticle>> vertexPFParticleAssociations(new art::Assns<recob::Vertex, recob::PFParticle>);

  // Initialize the UBXSecEvent
  ubxsec_event->Init();

  _run = ubxsec_event->run = e.id().run();
  _subrun = ubxsec_event->subrun = e.id().subRun();
  _event = ubxsec_event->event  = e.id().event();

  _is_data = e.isRealData();
  _is_mc   = !_is_data;

  if (_use_genie_info && _is_data) {
    std::cout << "[UBXSec] You have asked to use GENIE info but you are running on a data file.";
    std::cout << " _use_genie_info will be switched to false." << std::endl;
    _use_genie_info = false;
  }

  //::art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;

  // Prepare the dead region finder
  std::cout << "Recreate ch status map" << std::endl;
  const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  for (unsigned int ch = 0; ch < 8256; ch++) {
    deadRegionsFinder.SetChannelStatus(ch, chanFilt.Status(ch));
  }
  std::cout << "Now force reload BWires" << std::endl;
  deadRegionsFinder.CreateBWires();

  // Use '_detp' to find 'efield' and 'temp'
  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield = _detp -> Efield();
  double temp   = _detp -> Temperature();
  // Determine the drift velocity from 'efield' and 'temp'
  _drift_velocity = _detp -> DriftVelocity(efield,temp);
  if (_debug) std::cout << "[UBXSec] Using drift velocity = " << _drift_velocity << " cm/us, with E = " << efield << ", and T = " << temp << std::endl;

  // Collect tracks
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::TracksToHits           trackToHitsMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, pfParticleToTrackMap);
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, trackToHitsMap);

  // Collect showers
  lar_pandora::ShowerVector           _shower_v;
  lar_pandora::PFParticlesToShowers   _pfp_to_shower_map;;
  lar_pandora::LArPandoraHelper::CollectShowers(e, _pfp_producer, _shower_v, _pfp_to_shower_map);

  // Collect PFParticles and match Reco Particles to Hits
  //lar_pandora::PFParticleVector  recoParticleVector;
  //lar_pandora::PFParticleVector  recoNeutrinoVector;
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  //lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
  //lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kUseDaughters, true);

  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    std::cout << "[UBXSec] Cannote locate ubana::TPCObject." << std::endl;
  }
  art::FindManyP<ubana::FlashMatch> tpcobjToFlashMatchAssns(tpcobj_h, e, _neutrino_flash_match_producer);
  art::FindManyP<recob::Track>      tpcobjToTrackAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Shower>     tpcobjToShowerAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::PFParticle> tpcobjToPFPAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Vertex>     tpcobjToVertexAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToCosmicTagAssns(tpcobj_h, e, _geocosmictag_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToConsistency(tpcobj_h, e, _candidateconsistency_producer);

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
  art::FindManyP<recob::PFParticle> pfp_from_track(track_h, e, _pfp_producer);
  art::FindManyP<anab::Calorimetry> calos_from_track(track_h, e, _calorimetry_producer);

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
  std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfp_h);
  art::FindManyP<recob::Track> tracks_from_pfp(pfp_h, e, _pfp_producer);

  ubxsec_event->n_pfp = ubxsec_event->n_pfp_primary = 0;
  for (size_t i = 0; i < pfp_h->size(); i++) {
    ubxsec_event->n_pfp++;
    if ((*pfp_h)[i].IsPrimary())
      ubxsec_event->n_pfp_primary++;
  }

  // Get Giuseppe's Kalman Tracks
  art::FindManyP<recob::Track> trk_kalman_v(pfp_h, e, "pandoraNuKalmanTrack");

  // Get PID information
  art::FindManyP<anab::ParticleID> particleids_from_track (track_h, e, _particle_id_producer);
  if (!particleids_from_track.isValid()) {
    std::cout << "[UBXSec] anab::ParticleID is not valid." << std::endl;
  }
  std::cout << "[UBXSec] Numeber of particleids_from_track " << particleids_from_track.size() << std::endl;

  // Fill a std::map Track->ParticleID
  std::map<art::Ptr<recob::Track>, art::Ptr<anab::ParticleID>> track_to_pid_map;
  for (auto track : track_p_v) {
    std::vector<art::Ptr<anab::ParticleID>> pids = particleids_from_track.at(track.key());

    if(pids.size() == 0) 
      continue;

    for (auto pid : pids) {
      if (!pid->PlaneID().isValid) continue;
      int planenum = pid->PlaneID().Plane;
      if (planenum != 2) continue;
      track_to_pid_map[track] = pid;
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

  // Get MCSFitResult - Muon
  art::Handle<std::vector<recob::MCSFitResult> > mcsfitresult_mu_h;
  e.getByLabel(_mcsfitresult_mu_producer,mcsfitresult_mu_h);
  if(!mcsfitresult_mu_h.isValid()){
    std::cout << "[UBXSec] MCSFitResult product " << _mcsfitresult_mu_producer << " not found..." << std::endl;
    //throw std::exception();
  } 
  std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_mu_v;
  art::fill_ptr_vector(mcsfitresult_mu_v, mcsfitresult_mu_h);

  // Get MCSFitResult - Pi
  art::Handle<std::vector<recob::MCSFitResult> > mcsfitresult_pi_h;
  e.getByLabel(_mcsfitresult_pi_producer,mcsfitresult_pi_h);
  if(!mcsfitresult_pi_h.isValid()){
    std::cout << "[UBXSec] MCSFitResult product " << _mcsfitresult_pi_producer << " not found..." << std::endl;
    //throw std::exception();
  }
  std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_pi_v;
  art::fill_ptr_vector(mcsfitresult_pi_v, mcsfitresult_pi_h);

  /*
  std::cout << "mcsfitresult_mu_v.at(0)->fwdLogLikelihood(): " << mcsfitresult_mu_v.at(0)->fwdLogLikelihood() << std::endl;
  std::cout << "mcsfitresult_mu_v.at(0)->fwdMomentum(): " << mcsfitresult_mu_v.at(0)->fwdMomentum() << std::endl;
  std::cout << "mcsfitresult_mu_v.at(0)->bwdLogLikelihood(): " << mcsfitresult_mu_v.at(0)->bwdLogLikelihood() << std::endl;
  std::cout << "mcsfitresult_mu_v.at(0)->bwdMomentum(): " << mcsfitresult_mu_v.at(0)->bwdMomentum() << std::endl;
  std::cout << "mcsfitresult_mu_v.at(0)->bestMomentum(): " << mcsfitresult_mu_v.at(0)->bestMomentum() << std::endl;
  std::cout << "mcsfitresult_mu_v.at(0)->bestLogLikelihood(): " << mcsfitresult_mu_v.at(0)->bestLogLikelihood() << std::endl;
  std::cout << "mcsfitresult_mu_v.at(0)->deltaLogLikelihood(): " << mcsfitresult_mu_v.at(0)->deltaLogLikelihood() << std::endl;  
  */

  // Get the BNB Correction Weights
  double bnb_weight = 1.;
  if (_is_mc) {
    art::Handle<std::vector<evwgh::MCEventWeight>> eventweight_h;
    e.getByLabel(_eventweight_producer, eventweight_h);
    if(!eventweight_h.isValid()){
      std::cout << "[UBXSec] MCEventWeight product " << _eventweight_producer << " not found..." << std::endl;
      //throw std::exception();
    }
    std::vector<art::Ptr<evwgh::MCEventWeight>> eventweight_v;
    art::fill_ptr_vector(eventweight_v, eventweight_h);
    if (eventweight_v.size() > 0) {
      art::Ptr<evwgh::MCEventWeight> evt_wgt = eventweight_v.at(0);
      for (auto entry : evt_wgt->fWeight) {
        if (entry.first.find("bnbcorrection") != std::string::npos) {
          bnb_weight *= entry.second.at(0);
          std::cout << "[UBXSec] BNB Correction Weight: " << bnb_weight << std::endl;
        }
      }
    }
  }

  ubxsec_event->bnb_weight = bnb_weight;

  // GENIE reweigthing (systematics - pm1sigma)
  ubxsec_event->ResetGenieEventWeightVectorsPM1();
  if (_is_mc) {
    art::Handle<std::vector<evwgh::MCEventWeight>> genieeventweight_h;
    e.getByLabel(_genie_eventweight_pm1_producer, genieeventweight_h);
    if(!genieeventweight_h.isValid()){
      std::cout << "[UBXSec] MCEventWeight for GENIE reweight pm1sigma, product " << _genie_eventweight_pm1_producer << " not found..." << std::endl;
      //throw std::exception();
    } else {
      std::vector<art::Ptr<evwgh::MCEventWeight>> genieeventweight_v;
      art::fill_ptr_vector(genieeventweight_v, genieeventweight_h);
      if (genieeventweight_v.size() > 0) {
        art::Ptr<evwgh::MCEventWeight> evt_wgt = genieeventweight_v.at(0); // Just for the first nu interaction
        std::map<std::string, std::vector<double>> evtwgt_map = evt_wgt->fWeight;
        int countFunc = 0;
        // loop over the map and save the name of the function and the vector of weights for each function
        for(auto it : evtwgt_map) {
          std::string func_name = it.first;
          std::vector<double> weight_v = it.second; 
          ubxsec_event->evtwgt_genie_pm1_funcname.push_back(func_name);
          ubxsec_event->evtwgt_genie_pm1_weight.push_back(weight_v);
          ubxsec_event->evtwgt_genie_pm1_nweight.push_back(weight_v.size());
          countFunc++;
        }
        ubxsec_event->evtwgt_genie_pm1_nfunc = countFunc;
      }
    }
  }

  // GENIE reweigthing (systematics - multisim)
  ubxsec_event->ResetGenieEventWeightVectorsMultisim();
  if (_is_mc) {
    art::Handle<std::vector<evwgh::MCEventWeight>> genieeventweight_h;
    e.getByLabel(_genie_eventweight_multisim_producer, genieeventweight_h);
    if(!genieeventweight_h.isValid()){
      std::cout << "[UBXSec] MCEventWeight for GENIE reweight multisim, product " << _genie_eventweight_multisim_producer << " not found..." << std::endl;
      //throw std::exception();
    } else {
      std::vector<art::Ptr<evwgh::MCEventWeight>> genieeventweight_v;
      art::fill_ptr_vector(genieeventweight_v, genieeventweight_h);
      if (genieeventweight_v.size() > 0) {
        art::Ptr<evwgh::MCEventWeight> evt_wgt = genieeventweight_v.at(0); // Just for the first nu interaction
        std::map<std::string, std::vector<double>> evtwgt_map = evt_wgt->fWeight;
        int countFunc = 0;
        // loop over the map and save the name of the function and the vector of weights for each function
        for(auto it : evtwgt_map) {
          std::string func_name = it.first;
          std::vector<double> weight_v = it.second; 
          ubxsec_event->evtwgt_genie_multisim_funcname.push_back(func_name);
          ubxsec_event->evtwgt_genie_multisim_weight.push_back(weight_v);
          ubxsec_event->evtwgt_genie_multisim_nweight.push_back(weight_v.size());
          countFunc++;
        }
        ubxsec_event->evtwgt_genie_multisim_nfunc = countFunc;
      }
    }
  }

  // FLUX reweigthing (systematics - multisim)
  ubxsec_event->ResetFluxEventWeightVectorsMultisim();
  if (_is_mc) {
    art::Handle<std::vector<evwgh::MCEventWeight>> fluxeventweight_h;
    e.getByLabel(_flux_eventweight_multisim_producer, fluxeventweight_h);
    if(!fluxeventweight_h.isValid()){
      std::cout << "[UBXSec] MCEventWeight for FLUX reweight multisim, product " << _flux_eventweight_multisim_producer << " not found..." << std::endl;
      //throw std::exception();
    } else {
      std::vector<art::Ptr<evwgh::MCEventWeight>> fluxeventweight_v;
      art::fill_ptr_vector(fluxeventweight_v, fluxeventweight_h);
      if (fluxeventweight_v.size() > 0) {
        art::Ptr<evwgh::MCEventWeight> evt_wgt = fluxeventweight_v.at(0); // Just for the first nu interaction
        std::map<std::string, std::vector<double>> evtwgt_map = evt_wgt->fWeight;
        int countFunc = 0;
        // loop over the map and save the name of the function and the vector of weights for each function
        for(auto it : evtwgt_map) {
          std::string func_name = it.first;
          std::vector<double> weight_v = it.second; 
          ubxsec_event->evtwgt_flux_multisim_funcname.push_back(func_name);
          ubxsec_event->evtwgt_flux_multisim_weight.push_back(weight_v);
          ubxsec_event->evtwgt_flux_multisim_nweight.push_back(weight_v.size());
          countFunc++;
        }
        ubxsec_event->evtwgt_flux_multisim_nfunc = countFunc;
      }
    }
  }


  
  
  // pandoraCosmic PFPs (for cosmic removal studies)
  art::Handle<std::vector<recob::PFParticle>> pfp_cosmic_h;
  e.getByLabel("pandoraCosmic",pfp_cosmic_h);
  if(pfp_cosmic_h.isValid()){
    std::vector<art::Ptr<recob::PFParticle>> pfp_cosmic_v;
    art::fill_ptr_vector(pfp_cosmic_v, pfp_cosmic_h);
    ubxsec_event->n_primary_cosmic_pfp = 0;
    for (auto p : pfp_cosmic_v) {
      if (!p->IsPrimary()) continue;
      ubxsec_event->n_primary_cosmic_pfp++;
    }
  } else {
    std::cout << "[UBXSec] pandoraCosmic PFP product not found..." << std::endl;
  }

  // Flashes
  ::art::Handle<std::vector<recob::OpFlash>> beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);
  if( !beamflash_h.isValid() || beamflash_h->empty() ) {
    std::cout << "[UBXSec] Don't have good flashes." << std::endl;
  }

  ubxsec_event->nbeamfls = beamflash_h->size();
  ubxsec_event->beamfls_pe.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_time.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_z.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_spec.resize(ubxsec_event->nbeamfls);

  for (size_t n = 0; n < beamflash_h->size(); n++) {
    auto const& flash = (*beamflash_h)[n];
    ubxsec_event->beamfls_pe[n]   = flash.TotalPE();
    ubxsec_event->beamfls_time[n] = flash.Time();
    ubxsec_event->beamfls_z[n]    = flash.ZCenter();

    ubxsec_event->beamfls_spec[n].resize(32);
    ubxsec_event->candidate_flash_time = 0.;
    ubxsec_event->candidate_flash_z = 0.;
    double min_pe = -1;
    //if (_debug) std::cout << "[UBXSec] Reco beam flash pe: " << std::endl;
    for (unsigned int i = 0; i < 32; i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      ubxsec_event->beamfls_spec[n][opdet] = flash.PE(i);
      if (ubxsec_event->beamfls_time[n] > _beam_spill_start && ubxsec_event->beamfls_time[n] < _beam_spill_end) {
        // Find largest flash above threshold
        if (flash.TotalPE() > _total_pe_cut && flash.TotalPE() > min_pe) { 
          ubxsec_event->candidate_flash_time = flash.Time();
          ubxsec_event->candidate_flash_z = flash.ZCenter();
          min_pe = flash.TotalPE();
        }
        //if (_debug) std::cout << "\t PMT " << opdet << ": " << ubxsec_event->beamfls_spec[n][opdet] << std::endl;
      }
    }
  }


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
      if (_fiducial_volume.InFV(truth_nu_vtx)) ubxsec_event->fv = 1;
      else ubxsec_event->fv = 0;

      int n_genie_particles = 0;
      int n_genie_particles_charged = 0;
      for (int p = 0; p < mclist[iList]->NParticles(); p++) {
        const simb::MCParticle mc_par = mclist[iList]->GetParticle(p);
        if (mc_par.StatusCode() != 1) continue;
        n_genie_particles ++;
        const TParticlePDG* par_pdg = _database_pdg->GetParticle(mc_par.PdgCode());
        if (!par_pdg) continue;
        if (par_pdg->Charge() == 0) continue;
        n_genie_particles_charged ++;
      }

      ubxsec_event->ccnc            = mclist[iList]->GetNeutrino().CCNC();
      ubxsec_event->mode            = mclist[iList]->GetNeutrino().Mode();
      ubxsec_event->nupdg           = mclist[iList]->GetNeutrino().Nu().PdgCode();
      ubxsec_event->nu_e            = mclist[iList]->GetNeutrino().Nu().E();
      ubxsec_event->lep_costheta    = mclist[iList]->GetNeutrino().Lepton().Pz() / mclist[iList]->GetNeutrino().Lepton().P();
      ubxsec_event->lep_phi         = UBXSecHelper::GetPhi(mclist[iList]->GetNeutrino().Lepton().Px(), 
                                                           mclist[iList]->GetNeutrino().Lepton().Py(),
                                                           mclist[iList]->GetNeutrino().Lepton().Pz()); 
      ubxsec_event->genie_mult      = n_genie_particles;
      ubxsec_event->genie_mult_ch   = n_genie_particles_charged;

      ubxsec_event->tvtx_x.clear(); ubxsec_event->tvtx_x.clear(); ubxsec_event->tvtx_z.clear();
      for(size_t n = 0; n < mclist.size(); n++ ) {
        ubxsec_event->tvtx_x.emplace_back(mclist[n]->GetNeutrino().Nu().Vx());
        ubxsec_event->tvtx_y.emplace_back(mclist[n]->GetNeutrino().Nu().Vy());
        ubxsec_event->tvtx_z.emplace_back(mclist[n]->GetNeutrino().Nu().Vz());
      }

      ubxsec_event->nsignal = 0;
      if(ubxsec_event->nupdg==14 && ubxsec_event->ccnc==0 && ubxsec_event->fv==1) ubxsec_event->nsignal=1; 

      // Also save muon momentum if is signal
      ubxsec_event->true_muon_mom = -9999.;
      if (ubxsec_event->nsignal == 1) {
        for (int p = 0; p < mclist[iList]->NParticles(); p++) {
          auto const & mcp = mclist[iList]->GetParticle(p);
          if (mcp.Mother() != 0) continue;
          if (mcp.PdgCode() != 13) continue;
          ubxsec_event->true_muon_mom = mcp.P();
        }
      }
    } // neutrino origin
    else {
      ubxsec_event->ccnc = -1;
      ubxsec_event->nupdg = -1;
      ubxsec_event->nu_e = -1;
      ubxsec_event->lep_costheta = -9999.;
      ubxsec_event->true_muon_mom = -9999.;
    }
  } else {
    ubxsec_event->ccnc = -1;
    ubxsec_event->nupdg = -1;
    ubxsec_event->nu_e = -1;
    ubxsec_event->lep_costheta = -9999.;
    ubxsec_event->true_muon_mom = -9999.;
  }

  ubxsec_event->is_signal = false;
  if (ubxsec_event->ccnc == 0 && ubxsec_event->nupdg == 14 && ubxsec_event->fv == 1) {
    ubxsec_event->is_signal = true;
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

  art::Handle<std::vector<recob::OpHit>> ophit_cosmic_h;
  e.getByLabel("ophitCosmic", ophit_cosmic_h);
  if(!ophit_cosmic_h.isValid()) {
    std::cout << "[UBXSec] Cannot locate OpHits from ophitCosmic." << std::endl;
  }

  /*art::Handle<std::vector<raw::OpDetWaveform> > waveform_h;
  e.getByLabel("saturation", "OpdetCosmicHighGain", waveform_h);
  if(!waveform_h.isValid()) {
    std::cout << "[UBXSec] Cannot locate OpDetWaveform from saturation." << std::endl;
  }
  std::vector<raw::OpDetWaveform> const& waveform_v(*waveform_h);*/

  auto const& detectorClocks (*lar::providerFrom< detinfo::DetectorClocksService >());
  std::cout << "Trigger Time: " << detectorClocks.TriggerTime() << std::endl;
  std::cout << "Tick Period:  " << detectorClocks.OpticalClock().TickPeriod() << std::endl;


  //std::cout << "Printing waveforms" << std::endl;
  //for (auto w : waveform_v) {
  //  std::cout << "timestamp: " << w.TimeStamp() 
  //            << ", relative time: " << w.TimeStamp() - detectorClocks.TriggerTime() 
  //            << ", channel: " << geo->OpDetFromOpChannel(w.ChannelNumber())<< std::endl;
  //}


  // Check if the muon is reconstructed
  for (auto p : pfp_v) {
    auto mcghosts = mcghost_from_pfp.at(p.key());
    if (mcghosts.size() > 0) {
      art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghosts.at(0).key()).at(0);
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
      //const auto mc_truth = bt->TrackIDToMCTruth(mcpar->TrackId());
      if (mc_truth) {
        if (mc_truth->Origin() == simb::kBeamNeutrino 
            && mcpar->PdgCode() == 13 && mcpar->Mother() == 0) {
          ubxsec_event->muon_is_reco = true;
        }
      }
    }
  }
 
  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::ShowerVector    > shower_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;
  for (size_t slice = 0; slice < tpcobj_h->size(); slice++) {
    track_v_v.push_back(tpcobjToTrackAssns.at(slice));
    shower_v_v.push_back(tpcobjToShowerAssns.at(slice));
    pfp_v_v.push_back(tpcobjToPFPAssns.at(slice));
  }


  ubxsec_event->nslices = tpcobj_h->size();
  ubxsec_event->ResizeVectors(tpcobj_h->size());


  std::cout << "[UBXSec] --- SAVING INFORMATION" << std::endl;
  //ubxsec_event->vtx_resolution = -9999;

  ubxsec_event->n_tpcobj_nu_origin = 0;
  ubxsec_event->n_tpcobj_cosmic_origin = 0;

  std::vector<art::Ptr<recob::Track>> muon_candidate_track_per_slice_v; 
  std::vector<art::Ptr<recob::PFParticle>> muon_candidate_pfparticle_per_slice_v;
  std::vector<art::Ptr<recob::Vertex>> neutrino_candidate_vertex_per_slice_v;
  muon_candidate_track_per_slice_v.resize(ubxsec_event->nslices);
  muon_candidate_pfparticle_per_slice_v.resize(ubxsec_event->nslices);
  neutrino_candidate_vertex_per_slice_v.resize(ubxsec_event->nslices);


  //
  // THIS IS THE MAIN LOOP OVER THE 
  // TPCOBJECTS CANDIDATES IN THIS EVENT
  //

  for (unsigned int slice = 0; slice < tpcobj_h->size(); slice++){

    std::cout << "[UBXSec] >>> SLICE " << slice << std::endl;

    ubana::TPCObject tpcobj = (*tpcobj_h)[slice];

    ubxsec_event->slc_npfp[slice]    = tpcobj.GetNPFP();
    ubxsec_event->slc_ntrack[slice]  = tpcobj.GetNTracks();
    ubxsec_event->slc_nshower[slice] = tpcobj.GetNShowers();

    // Slice origin 
    ubxsec_event->slc_origin[slice] = tpcobj.GetOrigin();
    std::cout << "[UBXSec] \t Origin is " << ubxsec_event->slc_origin[slice] << std::endl;

    if (tpcobj.GetOrigin() == ubana::kBeamNeutrino || tpcobj.GetOrigin() == ubana::kMixed)
      ubxsec_event->n_tpcobj_nu_origin ++;
    else 
      ubxsec_event->n_tpcobj_cosmic_origin ++;

    // Slice origin extra
    ubxsec_event->slc_origin_extra[slice] = tpcobj.GetOriginExtra();
    std::cout << "[UBXSec] \t Origin extra is " << ubxsec_event->slc_origin_extra[slice] << std::endl;

    // Containment
    ubxsec_event->slc_iscontained[slice] = UBXSecHelper::TracksAreContained(tpcobj.GetTracks());

    // Reco vertex
    double reco_nu_vtx_raw[3];
    recob::Vertex tpcobj_nu_vtx = tpcobj.GetVertex();
    //tpcobj_nu_vtx.XYZ(reco_nu_vtx_raw);
    std::vector<art::Ptr<recob::Vertex>> recob_vtx_v = tpcobjToVertexAssns.at(slice);
    if (recob_vtx_v.size() > 0) {
      recob_vtx_v.at(0)->XYZ(reco_nu_vtx_raw);
      neutrino_candidate_vertex_per_slice_v.at(slice) = recob_vtx_v.at(0);
    } else {
      reco_nu_vtx_raw[0] = reco_nu_vtx_raw[1] = reco_nu_vtx_raw[2] = -9999;
    }

    // X position correction (time offset)
    double reco_nu_vtx[3];
    UBXSecHelper::GetTimeCorrectedPoint(reco_nu_vtx_raw, reco_nu_vtx, ubxsec_event->candidate_flash_time, _drift_velocity);

    // Space Charge correction
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    std::vector<double> sce_corr = SCE->GetPosOffsets(reco_nu_vtx[0],
                                                      reco_nu_vtx[1],
                                                      reco_nu_vtx[2]);
    std::cout << "[UBXSec] \t SCE correction in x, y, z = " << sce_corr.at(0) 
                                                    << ", " << sce_corr.at(1) 
                                                    << ", " << sce_corr.at(2) << std::endl;
    //reco_nu_vtx[0] += sce_corr.at(0);
    //reco_nu_vtx[1] -= sce_corr.at(1);
    //reco_nu_vtx[2] -= sce_corr.at(2);

    // Emplacing reco vertex
    ubxsec_event->slc_nuvtx_x[slice] = reco_nu_vtx[0];
    ubxsec_event->slc_nuvtx_y[slice] = reco_nu_vtx[1];
    ubxsec_event->slc_nuvtx_z[slice] = reco_nu_vtx[2];
    ubxsec_event->slc_nuvtx_fv[slice] = (_fiducial_volume.InFV(reco_nu_vtx) ? 1 : 0);
    std::cout << "[UBXSec] \t Reco vertex is " << ubxsec_event->slc_nuvtx_x[slice] << ", " << ubxsec_event->slc_nuvtx_y[slice] << ", " << ubxsec_event->slc_nuvtx_z[slice] << std::endl; 
    std::cout << "[UBXSec] \t Reco vertex is " << (ubxsec_event->slc_nuvtx_fv[slice]==1 ? "in" : "ouside") << " the FV." << std::endl;

    

    // Through-going?
    ubxsec_event->slc_geocosmictag[slice] = false;
    std::vector<art::Ptr<anab::CosmicTag>> geo_cosmic_tags = tpcobjToCosmicTagAssns.at(slice);
    if(geo_cosmic_tags.size() == 0 || geo_cosmic_tags.size() > 1) {
      std::cout << "[UBXSec] \t More than one Geo Cosmic Tag match per tpcobj ?!" << std::endl;
    } else {
      auto ct = geo_cosmic_tags.at(0);
      if (ct->CosmicScore() > _geo_cosmic_score_cut) {
        std::cout << "[UBXSec] \t This slice has been tagged as through-going cosmic" << std::endl;
        ubxsec_event->slc_geocosmictag[slice] = true;
      }
    }
    
    // Vertex resolution
    if (ubxsec_event->slc_origin[slice] == ubana::kBeamNeutrino) {
      ubxsec_event->vtx_resolution = sqrt( pow(ubxsec_event->slc_nuvtx_y[slice]-ubxsec_event->tvtx_y[0], 2) + pow(ubxsec_event->slc_nuvtx_z[slice]-ubxsec_event->tvtx_z[0], 2) );
    } 

    // Multiplicity
    int p, t, s;
    tpcobj.GetMultiplicity(p, t, s);
    ubxsec_event->slc_mult_pfp[slice] = p;
    ubxsec_event->slc_mult_track[slice] = t;
    ubxsec_event->slc_mult_shower[slice] = s;
    ubxsec_event->slc_mult_track_tolerance[slice] = tpcobj.GetNTracksCloseToVertex(_tolerance_track_multiplicity);

    // Candidate Consistency
    ubxsec_event->slc_consistency[slice] = true;
    std::vector<art::Ptr<anab::CosmicTag>> consistency_tags = tpcobjToConsistency.at(slice);
    if(consistency_tags.size() == 0 || consistency_tags.size() > 1) {
      std::cout << "[UBXSec] \t More than one Consistency Tag match per tpcobj ?!" << std::endl;
    } else {
      auto ct = consistency_tags.at(0);
      std::cout << "[UBXSec] \t Candidate Consistency Score: " << ct->CosmicScore() << std::endl;
      if (ct->CosmicType() != anab::CosmicTagID_t::kNotTagged) {
        std::cout << "[UBXSec] \t This slice has been tagged as not consistent. Type " << ct->CosmicType() 
                  << ", score: " << ct->CosmicScore() << std::endl;
        ubxsec_event->slc_consistency[slice] = false;
      }
      ubxsec_event->slc_consistency_score[slice] = ct->CosmicScore();
    }

    // Neutrino Flash match
    ubxsec_event->slc_flsmatch_score[slice] = -9999;
    std::vector<art::Ptr<ubana::FlashMatch>> pfpToFlashMatch_v = tpcobjToFlashMatchAssns.at(slice);
    if (pfpToFlashMatch_v.size() > 1) {
      std::cout << "[UBXSec] \t More than one flash match per nu pfp ?!" << std::endl;
      continue;
    } else if (pfpToFlashMatch_v.size() == 0){
      std::cout << "[UBXSec] \t No Flash-Match for this TPCObject" << std::endl;
    } else {
      ubxsec_event->slc_flsmatch_score[slice]       = pfpToFlashMatch_v[0]->GetScore(); 
      ubxsec_event->slc_flsmatch_qllx[slice]        = pfpToFlashMatch_v[0]->GetEstimatedX();
      ubxsec_event->slc_flsmatch_tpcx[slice]        = pfpToFlashMatch_v[0]->GetTPCX();
      ubxsec_event->slc_flsmatch_t0[slice]          = pfpToFlashMatch_v[0]->GetT0();
      ubxsec_event->slc_flsmatch_hypoz[slice]       = UBXSecHelper::GetFlashZCenter(pfpToFlashMatch_v[0]->GetHypoFlashSpec());
      ubxsec_event->slc_flsmatch_xfixed_chi2[slice] = pfpToFlashMatch_v[0]->GetXFixedChi2();
      ubxsec_event->slc_flsmatch_xfixed_ll[slice]   = pfpToFlashMatch_v[0]->GetXFixedLl();
      ubxsec_event->slc_flshypo_xfixed_spec[slice]  = pfpToFlashMatch_v[0]->GetXFixedHypoFlashSpec();
      ubxsec_event->slc_flshypo_spec[slice]         = pfpToFlashMatch_v[0]->GetHypoFlashSpec();
      //for (auto v : _slc_flshypo_spec[slice]) std::cout << "Hypo PE: " << v << std::endl;

      std::cout << "[UBXSec] \t FM score:       " << ubxsec_event->slc_flsmatch_score[slice] << std::endl;
      std::cout << "[UBXSec] \t qllx - tpcx is: " << ubxsec_event->slc_flsmatch_qllx[slice] - ubxsec_event->slc_flsmatch_tpcx[slice] << std::endl;
    }

    // Hits
    int nhits_u, nhits_v, nhits_w;
    UBXSecHelper::GetNumberOfHitsPerPlane(e, _pfp_producer, track_v_v[slice], nhits_u, nhits_v, nhits_w);
    ubxsec_event->slc_nhits_u[slice] = nhits_u;
    ubxsec_event->slc_nhits_v[slice] = nhits_v;
    ubxsec_event->slc_nhits_w[slice] = nhits_w;

    // Longest track and check boundary
    recob::Track lt;
    if (UBXSecHelper::GetLongestTrackFromTPCObj(track_v_v[slice], lt)){
      ubxsec_event->slc_longesttrack_length[slice] = lt.Length();
      ubxsec_event->slc_longesttrack_phi[slice]   = UBXSecHelper::GetCorrectedPhi(lt, tpcobj_nu_vtx);
      ubxsec_event->slc_longesttrack_theta[slice] = UBXSecHelper::GetCorrectedCosTheta(lt, tpcobj_nu_vtx);
      ubxsec_event->slc_longesttrack_iscontained[slice] = UBXSecHelper::TrackIsContained(lt);
      int vtx_ok;
      ubxsec_event->slc_crosses_top_boundary[slice] = (UBXSecHelper::IsCrossingTopBoundary(lt, vtx_ok) ? 1 : 0);
    } else {
      ubxsec_event->slc_longesttrack_length[slice] = -9999;
    }

    // Longest shower
    recob::Shower ls;
    if (UBXSecHelper::GetLongestShowerFromTPCObj(shower_v_v[slice], ls)) {
      ubxsec_event->slc_longestshower_length[slice] = ls.Length();
      ubxsec_event->slc_longestshower_openangle[slice] = ls.OpenAngle();
      ubxsec_event->slc_longestshower_startx[slice] = ls.ShowerStart().X();
      ubxsec_event->slc_longestshower_starty[slice] = ls.ShowerStart().Y();
      ubxsec_event->slc_longestshower_startz[slice] = ls.ShowerStart().Z();
      ubxsec_event->slc_longestshower_phi[slice]   = UBXSecHelper::GetPhi(ls.Direction());
      ubxsec_event->slc_longestshower_theta[slice] = UBXSecHelper::GetCosTheta(ls.Direction());
    }

    // ACPT
    ubxsec_event->slc_acpt_outoftime[slice] = 0;
    for (unsigned int t = 0; t < track_v_v[slice].size(); t++) {
      if(opfls_ptr_coll_v.at(track_v_v[slice][t].key()).size()>1) {
        std::cout << "[UBXSec] \t More than 1 association found (ACPT)!" << std::endl;
        //throw std::exception();
      } else if (opfls_ptr_coll_v.at(track_v_v[slice][t].key()).size()==0){
        continue;
      } else {
        art::Ptr<recob::OpFlash> flash_ptr = opfls_ptr_coll_v.at(track_v_v[slice][t].key()).at(0);
        if (flash_ptr->Time() < _beam_spill_start || flash_ptr->Time() > _beam_spill_end) {
          ubxsec_event->slc_acpt_outoftime[slice] = 1;
        }
      }
    }

    // Track quality
    ubxsec_event->slc_kalman_chi2[slice] = -9999;
    for (unsigned int t = 0; t < pfp_v_v[slice].size(); t++) {
      if(trk_kalman_v.at(pfp_v_v[slice][t].key()).size()>1) {
        std::cout << "[UBXSec] \t TQ more than one track per PFP, ntracks " << trk_kalman_v.at(pfp_v_v[slice][t].key()).size() << std::endl;
      } else if (trk_kalman_v.at(pfp_v_v[slice][t].key()).size()==0){
        continue;
      } else {
        art::Ptr<recob::Track> trk_ptr = trk_kalman_v.at(pfp_v_v[slice][t].key()).at(0);
        ubxsec_event->slc_kalman_chi2[slice] = trk_ptr->Chi2();
        ubxsec_event->slc_kalman_ndof[slice] = trk_ptr->Ndof();
      }
    }
    bool goodTrack = true;
    for (auto trk : track_v_v[slice]) {
      if (deadRegionsFinder.NearDeadReg2P( (trk->Vertex()).Y(), (trk->Vertex()).Z(), _minimumDistDeadReg )  ||
          deadRegionsFinder.NearDeadReg2P( (trk->End()).Y(),    (trk->End()).Z(),    _minimumDistDeadReg )  ||
          deadRegionsFinder.NearDeadRegCollection(trk->Vertex().Z(), _minimumDistDeadReg) ||
          deadRegionsFinder.NearDeadRegCollection(trk->End().Z(),    _minimumDistDeadReg) ||
          !UBXSecHelper::TrackPassesHitRequirment(e, _pfp_producer, trk, _minimumHitRequirement) ) {
        goodTrack = false;
        break;
      }
    } 

    if (goodTrack) ubxsec_event->slc_passed_min_track_quality[slice] = true;
    else ubxsec_event->slc_passed_min_track_quality[slice] = false;

    std::cout << "[UBXSec] \t " << (goodTrack ? "Passed" : "Did not pass") 
              << " minimum track quality." << std::endl;

    // Vertex quality
    recob::Vertex slice_vtx = tpcobj.GetVertex();
    double slice_vtx_xyz[3];
    slice_vtx.XYZ(slice_vtx_xyz);
    ubxsec_event->slc_passed_min_vertex_quality[slice] = true;
    if (deadRegionsFinder.NearDeadReg2P(slice_vtx_xyz[1], slice_vtx_xyz[2], _minimumDistDeadReg))
      ubxsec_event->slc_passed_min_vertex_quality[slice] = false;

    // Channel status
    ubxsec_event->slc_nuvtx_closetodeadregion_u[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 0) ? 1 : 0);
    ubxsec_event->slc_nuvtx_closetodeadregion_v[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 1) ? 1 : 0);
    ubxsec_event->slc_nuvtx_closetodeadregion_w[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 2) ? 1 : 0);

    // Vertex check
    ubxsec::VertexCheck vtxCheck(track_v_v[slice], slice_vtx);
    ubxsec_event->slc_vtxcheck_angle[slice] = vtxCheck.AngleBetweenLongestTracks();
    
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
      size_t opdet = geo->OpDetFromOpChannel(ophit.OpChannel());
      //std::cout << "OpHit::  OpDet: " << opdet
      //          << ", PeakTime: " << ophit.PeakTime()
      //          << ", PE: " << _pecalib.BeamPE(opdet,ophit.Area(),ophit.Amplitude()) << std::endl;
      if(_make_ophit_csv) _csvfile2 << oh << "," << opdet << "," << ophit.PeakTime() << "," << _pecalib.BeamPE(opdet,ophit.Area(),ophit.Amplitude()) << std::endl;
      if (ophit.OpChannel() != this_opch) continue;
      if (ophit.PeakTime() > _beam_spill_start && ophit.PeakTime() < _beam_spill_end) {
        n_intime_ophits ++;
        n_intime_pe += _pecalib.BeamPE(opdet,ophit.Area(),ophit.Amplitude());
      }     
    } // end loop ophit

    ubxsec_event->slc_n_intime_pe_closestpmt[slice] = n_intime_pe;

    /*for (size_t oh = 0; oh < ophit_cosmic_h->size(); oh++) {
      auto const & ophit = (*ophit_cosmic_h)[oh];
      if (ophit.PeakTime() < -150 || ophit.PeakTime() > -50) continue;
      size_t opdet = geo->OpDetFromOpChannel(ophit.OpChannel());
      std::cout << "Cosmic Disc OpHit::  OpDet: " << opdet
                << ", PeakTime: " << ophit.PeakTime()
                << ", PE: " << _pecalib.CosmicPE(opdet,ophit.Area(),ophit.Amplitude()) << std::endl;
    }*/

    // Distance from recon nu vertex to thefar away track in TPCObject
    //_slc_maxdistance_vtxtrack = UBXSecHelper::GetMaxTrackVertexDistance();

    // Other showers in the event
    std::vector<art::Ptr<recob::Shower>> other_showers;
    bool ignore_shower = false;
    for (size_t s = 0; s < _shower_v.size(); s++) {

      ignore_shower = false;

      // Check this shower is not in this TPCObject
      for (size_t this_s = 0; this_s < shower_v_v.at(slice).size(); this_s++) {
        if (_shower_v.at(s).key() == shower_v_v.at(slice).at(this_s).key()) {
          ignore_shower = true;
          continue;
        }
      }

      if (ignore_shower) continue;
     
      other_showers.push_back(_shower_v.at(s));
    }

    double max_length = -1, index_max_length = -1;
    double max_costheta = -1e9, index_max_costheta = -1;
    double min_flashvtxdistance = 1e9, index_min_flashvtxdistance = -1;

    for (size_t s = 0; s < other_showers.size(); s++) {

      if (other_showers.at(s)->Length() > max_length) {
        max_length = other_showers.at(s)->Length();
        index_max_length = s;
      }

      double costheta = UBXSecHelper::GetCosTheta(other_showers.at(s)->Direction());
      if (costheta > max_costheta) {
        costheta = max_costheta;
        index_max_costheta = s;
      }

      double distance = std::abs(other_showers.at(s)->ShowerStart().Z() - ubxsec_event->candidate_flash_z);
      if (distance < min_flashvtxdistance) {
        min_flashvtxdistance = distance;
        index_min_flashvtxdistance = s;
      }
    }

    if (index_max_length != -1) {
      auto shower = other_showers.at(index_max_length);
      ubxsec_event->slc_othershowers_longest_length[slice] = shower->Length();
      ubxsec_event->slc_othershowers_longest_startx[slice] = shower->ShowerStart().X();
      ubxsec_event->slc_othershowers_longest_starty[slice] = shower->ShowerStart().Y();
      ubxsec_event->slc_othershowers_longest_startz[slice] = shower->ShowerStart().Z();
      ubxsec_event->slc_othershowers_longest_phi[slice] = UBXSecHelper::GetPhi(shower->Direction());
      ubxsec_event->slc_othershowers_longest_theta[slice] = UBXSecHelper::GetCosTheta(shower->Direction());
      ubxsec_event->slc_othershowers_longest_openangle[slice] = shower->OpenAngle();
    } else {
      ubxsec_event->slc_othershowers_longest_length[slice] = -1;
      ubxsec_event->slc_othershowers_longest_startx[slice] = -1;
      ubxsec_event->slc_othershowers_longest_starty[slice] = -1;
      ubxsec_event->slc_othershowers_longest_startz[slice] = -1;
      ubxsec_event->slc_othershowers_longest_phi[slice] = -1;
      ubxsec_event->slc_othershowers_longest_theta[slice] = -1;
      ubxsec_event->slc_othershowers_longest_openangle[slice] = -1;      
    }
    if (index_max_costheta != -1) {
      auto shower = other_showers.at(index_max_costheta);
      ubxsec_event->slc_othershowers_forward_length[slice] = shower->Length();
      ubxsec_event->slc_othershowers_forward_startx[slice] = shower->ShowerStart().X();
      ubxsec_event->slc_othershowers_forward_starty[slice] = shower->ShowerStart().Y();
      ubxsec_event->slc_othershowers_forward_startz[slice] = shower->ShowerStart().Z();
      ubxsec_event->slc_othershowers_forward_phi[slice] = UBXSecHelper::GetPhi(shower->Direction());
      ubxsec_event->slc_othershowers_forward_theta[slice] = UBXSecHelper::GetCosTheta(shower->Direction());
      ubxsec_event->slc_othershowers_forward_openangle[slice] = shower->OpenAngle();
    } else {
      ubxsec_event->slc_othershowers_forward_length[slice] = -1;
      ubxsec_event->slc_othershowers_forward_startx[slice] = -1;
      ubxsec_event->slc_othershowers_forward_starty[slice] = -1;
      ubxsec_event->slc_othershowers_forward_startz[slice] = -1;
      ubxsec_event->slc_othershowers_forward_phi[slice] = -1;
      ubxsec_event->slc_othershowers_forward_theta[slice] = -1;
      ubxsec_event->slc_othershowers_forward_openangle[slice] = -1;
    }
    if (index_min_flashvtxdistance != -1) {
      auto shower = other_showers.at(index_min_flashvtxdistance);
      ubxsec_event->slc_othershowers_flashmatch_length[slice] = shower->Length();
      ubxsec_event->slc_othershowers_flashmatch_startx[slice] = shower->ShowerStart().X();
      ubxsec_event->slc_othershowers_flashmatch_starty[slice] = shower->ShowerStart().Y();
      ubxsec_event->slc_othershowers_flashmatch_startz[slice] = shower->ShowerStart().Z();
      ubxsec_event->slc_othershowers_flashmatch_phi[slice] = UBXSecHelper::GetPhi(shower->Direction());
      ubxsec_event->slc_othershowers_flashmatch_theta[slice] = UBXSecHelper::GetCosTheta(shower->Direction());
      ubxsec_event->slc_othershowers_flashmatch_openangle[slice] = shower->OpenAngle();
    } else {
      ubxsec_event->slc_othershowers_flashmatch_length[slice] = -1;
      ubxsec_event->slc_othershowers_flashmatch_startx[slice] = -1;
      ubxsec_event->slc_othershowers_flashmatch_starty[slice] = -1;
      ubxsec_event->slc_othershowers_flashmatch_startz[slice] = -1;
      ubxsec_event->slc_othershowers_flashmatch_phi[slice] = -1;
      ubxsec_event->slc_othershowers_flashmatch_theta[slice] = -1;
      ubxsec_event->slc_othershowers_longest_openangle[slice] = -1;
    }

    std::cout << "[UBXSec] Shower info saved." << std::endl;


    // Muon Candidate
    _muon_finder.Reset();
    _muon_finder.SetTracks(track_v_v[slice]);
    _muon_finder.SetTrackToPIDMap(track_to_pid_map);
    art::Ptr<recob::Track> candidate_track;

    bool muon_cand_exists = _muon_finder.GetCandidateTrack(candidate_track);

    if (muon_cand_exists) {

      bool fully_contained = _fiducial_volume.InFV(candidate_track->Vertex(), candidate_track->End());

      ubxsec_event->slc_muoncandidate_exists[slice]    = true;
      ubxsec_event->slc_muoncandidate_contained[slice] = fully_contained;
      ubxsec_event->slc_muoncandidate_length[slice]    = candidate_track->Length();
      ubxsec_event->slc_muoncandidate_phi[slice]       = UBXSecHelper::GetCorrectedPhi((*candidate_track), tpcobj_nu_vtx); 
      ubxsec_event->slc_muoncandidate_theta[slice]     = UBXSecHelper::GetCorrectedCosTheta((*candidate_track), tpcobj_nu_vtx);
      ubxsec_event->slc_muoncandidate_mom_range[slice] = _trk_mom_calculator.GetTrackMomentum(candidate_track->Length(), 13);
      //ubxsec_event->slc_muoncandidate_mom_mcs[slice]   = _trk_mom_calculator.GetMomentumMultiScatterLLHD(candidate_track);

      // For MCS first check the track direction is rigth
      TVector3 temp(reco_nu_vtx[0], reco_nu_vtx[1], reco_nu_vtx[2]);
      bool track_direction_correct = (candidate_track->Vertex() - temp).Mag() < (candidate_track->End() - temp).Mag();
      if (track_direction_correct) {
        ubxsec_event->slc_muoncandidate_mom_mcs[slice] = mcsfitresult_mu_v.at(candidate_track.key())->fwdMomentum();
        ubxsec_event->slc_muoncandidate_mcs_ll[slice]  = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
        ubxsec_event->slc_muoncandidate_mom_mcs_pi[slice] = mcsfitresult_pi_v.at(candidate_track.key())->fwdMomentum();
        std::cout << "Muon MCS LL: " << mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood() << std::endl;
        std::cout << "Pion MCS LL: " << mcsfitresult_pi_v.at(candidate_track.key())->fwdLogLikelihood() << std::endl;
      } else {
        ubxsec_event->slc_muoncandidate_mom_mcs[slice] = mcsfitresult_mu_v.at(candidate_track.key())->bwdMomentum();
        ubxsec_event->slc_muoncandidate_mcs_ll[slice]  = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
        ubxsec_event->slc_muoncandidate_mom_mcs_pi[slice] = mcsfitresult_pi_v.at(candidate_track.key())->bwdMomentum();
        std::cout << "Muon MCS LL: " << mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood() << std::endl;
        std::cout << "Pion MCS LL: " << mcsfitresult_pi_v.at(candidate_track.key())->bwdLogLikelihood() << std::endl;
      }
      // Also see if the track is recon going downwards (for cosmic studies)
      bool track_going_down = candidate_track->Vertex().Y() > candidate_track->End().Y();

      // Look at calorimetry for the muon candidate
      std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_track.at(candidate_track.key());
      ubxsec_event->slc_muoncandidate_dqdx_v[slice] = UBXSecHelper::GetDqDxVector(calos);
      ubxsec_event->slc_muoncandidate_dqdx_trunc[slice] = UBXSecHelper::GetDqDxTruncatedMean(calos);
      ubxsec_event->slc_muoncandidate_dqdx_u_trunc[slice] = UBXSecHelper::GetDqDxTruncatedMean(calos, 0);
      ubxsec_event->slc_muoncandidate_dqdx_v_trunc[slice] = UBXSecHelper::GetDqDxTruncatedMean(calos, 1);
      ubxsec_event->slc_muoncandidate_mip_consistency[slice] = _muon_finder.MIPConsistency(ubxsec_event->slc_muoncandidate_dqdx_trunc[slice],
                                                                                           ubxsec_event->slc_muoncandidate_length[slice]);
      ubxsec_event->slc_muoncandidate_mip_consistency2[slice] = _muon_finder.SVMPredict(ubxsec_event->slc_muoncandidate_dqdx_trunc[slice],
                                                                                        ubxsec_event->slc_muoncandidate_length[slice]);
      std::cout << "[UBXSec] \t Truncated mean dQ/ds for candidate is (plane 0): " << ubxsec_event->slc_muoncandidate_dqdx_trunc[slice] << std::endl;
      std::cout << "[UBXSec] \t Truncated mean dQ/ds for candidate is (plane 1): " << ubxsec_event->slc_muoncandidate_dqdx_u_trunc[slice] << std::endl;
      std::cout << "[UBXSec] \t Truncated mean dQ/ds for candidate is (plane 2): " << ubxsec_event->slc_muoncandidate_dqdx_v_trunc[slice] << std::endl;
      std::cout << "[UBXSec] \t MIP consistent ? : " << (ubxsec_event->slc_muoncandidate_mip_consistency[slice] ? "YES" : "NO") << std::endl;

      // Get the related PFP
      art::Ptr<recob::PFParticle> candidate_pfp = pfp_from_track.at(candidate_track.key()).at(0);
      const auto mcghosts = mcghost_from_pfp.at(candidate_pfp.key());
      if (mcghosts.size() > 0) {
        art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghost_from_pfp.at(candidate_pfp.key()).at(0).key()).at(0);
        const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
        //const auto mc_truth = bt->TrackIDToMCTruth(mcpar->TrackId());
        ubxsec_event->slc_muoncandidate_truepdg[slice] = mcpar->PdgCode();
        if (mc_truth) {

          // Check the true origin of the candidate PFP
          if (mc_truth->Origin() == simb::kBeamNeutrino) {
            ubxsec_event->slc_muoncandidate_trueorigin[slice] = ubana::kBeamNeutrino;
          } else if (mc_truth->Origin() == simb::kCosmicRay) {
            ubxsec_event->slc_muoncandidate_trueorigin[slice] = ubana::kCosmicRay;
          }

          // Now make momentum distributions
          if (mc_truth->Origin() == simb::kBeamNeutrino && mcpar->PdgCode() == 13 && mcpar->Mother() == 0) {
            _h_mom_true_mcs->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            if (fully_contained) {
              _mom_true_contained = mcpar->P();
              _mom_mcs_contained = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
              _mom_range_contained = ubxsec_event->slc_muoncandidate_mom_range[slice];
              _mom_tree_contained->Fill();
              _h_mom_true_mcs_contained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
              _h_mom_true_range_contained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_range[slice]);
              _h_mom_range_mcs_contained->Fill(ubxsec_event->slc_muoncandidate_mom_range[slice],
                                               ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            } else {
              _mom_true_uncontained = mcpar->P();
              _mom_mcs_uncontained = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
              _mom_tree_uncontained->Fill();
              _h_mom_true_mcs_uncontained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            }
          }
          if (mc_truth->Origin() == simb::kCosmicRay && (mcpar->PdgCode() == 13 || mcpar->PdgCode() == -13)) {
            _mom_cosmic_true = mcpar->P();
            _mom_cosmic_mcs = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
            _mom_cosmic_mcs_downforced = track_going_down ?   mcsfitresult_mu_v.at(candidate_track.key())->fwdMomentum() 
                                                            : mcsfitresult_mu_v.at(candidate_track.key())->bwdMomentum(); 
            _mom_cosmic_range = ubxsec_event->slc_muoncandidate_mom_range[slice];
            _mom_cosmic_down = track_going_down;
            _mom_cosmic_tree->Fill();
          }
        }
        //std::cout << ">>>>>>>>>>>>>>>>>>>>> MCP has pdg " << mcpar->PdgCode() << std::endl;
      }
 

      if (_debug) std::cout << "[UBXSec] \t Muon Candidate Found" << std::endl;
      if (_debug) std::cout << "[UBXSec] \t \t Length:          " << ubxsec_event->slc_muoncandidate_length[slice] << std::endl;
      if (_debug) std::cout << "[UBXSec] \t \t Mom by Range:    " << ubxsec_event->slc_muoncandidate_mom_range[slice] << std::endl;
      if (_debug) std::cout << "[UBXSec] \t \t Mom by MCS:      " << ubxsec_event->slc_muoncandidate_mom_mcs[slice] << std::endl;
      if (_debug) std::cout << "[UBXSec] \t \t Fully Contained? " << (ubxsec_event->slc_muoncandidate_contained[slice] ? "YES" : "NO") << std::endl;

      // Try look at MCS for stopping muons
      bool down_track = candidate_track->Vertex().Y() > candidate_track->End().Y();
      double f_ll = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
      double b_ll = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
      ubxsec_event->slc_muoncandidate_mcs_delta_ll[slice] = f_ll - b_ll;
      if (!down_track) ubxsec_event->slc_muoncandidate_mcs_delta_ll[slice] = b_ll - f_ll;


      //
      // Look at residuals
      //
      std::vector<TVector3> hit_v; // a vec of hits from coll plane
      std::vector<TVector3> track_v; // a vec of hits from coll plane

      // Collect hits
      auto iter = trackToHitsMap.find(candidate_track);
      if (iter != trackToHitsMap.end()) {
        std::vector<art::Ptr<recob::Hit>> hits = iter->second;
        for (auto hit : hits) {
          if (hit->View() == 2) {
            TVector3 h (hit->WireID().Wire, hit->PeakTime(), 0);
            //std::cout << "emplacing hit with wire " << h.X() << ", and time " << h.Y() << std::endl;
            hit_v.emplace_back(h);
          }
        }
      }

      // Collect track points
      for (size_t i = 0; i < candidate_track->NumberTrajectoryPoints(); i++) {
        try {
          if (candidate_track->HasValidPoint(i)) {
            TVector3 trk_pt = candidate_track->LocationAtPoint(i);
            double wire = geo->NearestWire(trk_pt, 2);
            double time = fDetectorProperties->ConvertXToTicks(trk_pt.X(), geo::PlaneID(0,0,2));
            TVector3 p (wire, time, 0.);
            //std::cout << "emplacing track point on wire " << p.X() << ", and time " << p.Y() << std::endl;
            track_v.emplace_back(p);
          }
        } catch (...) {
          continue;
        }
      }
      
      if (_debug) std::cout << "[UBXSec] \t \t Hit points: " << hit_v.size() << ", track points: " << track_v.size() << std::endl;
      ubana::TrackQuality _track_quality;
      _track_quality.SetTrackPoints(track_v);
      _track_quality.SetHitCollection(hit_v);
      std::pair<double, double> residual_mean_std = _track_quality.GetResiduals();
      if (_debug) std::cout << "[UBXSec] \t \t Residuals, mean " << residual_mean_std.first << ", std " << residual_mean_std.second << std::endl;
      std::pair<double,int> dist_wire_pair = _track_quality.GetTrackGap();

      int start_wire = dist_wire_pair.second;
      int end_wire = dist_wire_pair.second + dist_wire_pair.first;

      if (start_wire > end_wire) 
        std::swap(start_wire, end_wire);

      // Create a channel to wire map
      std::map<int, int> wire_to_channel;
      for (unsigned int ch = 0; ch < 8256; ch++) {
        std::vector< geo::WireID > wire_v = geo->ChannelToWire(ch);
        wire_to_channel[wire_v[0].Wire] = ch;
      }

      int n_dead_wires = 0;

      for (int wire = start_wire; wire < end_wire; wire++) {

        int channel = wire_to_channel[wire];

        // Channel statuses: 1=dead, 3=noisy, 4=good
        if (chanFilt.Status(channel) < 4) {
          n_dead_wires++;
        }
      }

      std::cout << "Gap of " << end_wire-start_wire << ", n of dead wires in it: " << n_dead_wires << std::endl;
      double r = _track_quality.GetR();
      std::cout << "The r value is: " << r << std::endl;


      double n_hits_in_cluster = 0;
      auto it = recoParticlesToHits.find(candidate_pfp);
      if (it != recoParticlesToHits.end()) {
        for (auto h : it->second) {
          if (h->View() == 2) {
            n_hits_in_cluster++;
          }
        }
      }

      double ratio = (double)hit_v.size()/n_hits_in_cluster;

      std::cout << "Percentage of used hit in cluster = " << ratio << std::endl;

      // Also look at the scattering angle
      std::vector<TVector3> dir_v;
      for (size_t p = 0; p < candidate_track->NumberTrajectoryPoints(); p++) {
        if (!candidate_track->HasValidPoint(p)) continue;
        dir_v.push_back(candidate_track->DirectionAtPoint(p));
      }
      std::vector<double> angle_v;
      for (size_t p = 0; p < dir_v.size()-1; p++) {
        double angle = dir_v.at(p).Angle(dir_v.at(p+1));
        angle_v.push_back(angle);
      }
      std::sort(angle_v.begin(), angle_v.end());

      double max_angle = -1;
      if (angle_v.size() != 0)
        max_angle = angle_v.at(angle_v.size()-1) / TMath::Pi() * 180.;

      std::cout << "Max scattering angle = " << max_angle << std::endl;


      ubxsec_event->slc_muoncandidate_residuals_mean[slice] = residual_mean_std.first;
      ubxsec_event->slc_muoncandidate_residuals_std[slice] = residual_mean_std.second;
      ubxsec_event->slc_muoncandidate_wiregap[slice] = end_wire-start_wire;
      ubxsec_event->slc_muoncandidate_wiregap_dead[slice] = n_dead_wires;
      ubxsec_event->slc_muoncandidate_linearity[slice] = r;
      ubxsec_event->slc_muoncandidate_perc_used_hits_in_cluster[slice] = ratio;
      ubxsec_event->slc_muoncandidate_maxscatteringangle[slice] = max_angle;





      muon_candidate_track_per_slice_v.at(slice) = candidate_track;
      muon_candidate_pfparticle_per_slice_v.at(slice) = candidate_pfp;


    } else {

      ubxsec_event->slc_muoncandidate_exists[slice]    = false;
      ubxsec_event->slc_muoncandidate_length[slice]    = -9999;
      ubxsec_event->slc_muoncandidate_phi[slice]       = -9999;
      ubxsec_event->slc_muoncandidate_theta[slice]     = -9999;
      ubxsec_event->slc_muoncandidate_mom_range[slice] = -9999;
      ubxsec_event->slc_muoncandidate_mom_mcs[slice]   = -9999;

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
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpars[0]->TrackId());
      //const auto mc_truth = bt->TrackIDToMCTruth(mcpars[0]->TrackId()); 
      if (!mc_truth) {
        std::cerr << "[UBXSec] Problem with MCTruth pointer." << std::endl;
        continue;
      }
      if (mc_truth->Origin() == simb::kBeamNeutrino &&
          mcpars[0]->PdgCode() == 13 && mcpars[0]->Mother() == 0) {

        //ubxsec_event->muon_is_reco = true;
        //ubxsec_event->muon_reco_pur = ubxsec_event->muon_reco_eff = -9999;
        auto iter = recoParticlesToHits.find(pfp);
        if (iter != recoParticlesToHits.end()) {
          //CHECK! UBXSecHelper::GetTrackPurityAndEfficiency((*iter).second, ubxsec_event->muon_reco_pur, ubxsec_event->muon_reco_eff);
        }
        ubxsec_event->true_muon_mom_matched = mcpars[0]->P();

        if (_fiducial_volume.InFV(mcpars[0]->Vx(), mcpars[0]->Vy(), mcpars[0]->Vz()) && 
          _fiducial_volume.InFV(mcpars[0]->EndX(), mcpars[0]->EndY(), mcpars[0]->EndZ())) {
          ubxsec_event->mc_muon_contained = true;
        }
      } 
      

      std::vector<art::Ptr<recob::Track>> tracks = tracks_from_pfp.at(pfp.key());
      std::cout << "[UBXSec] \t\t n tracks ass to this pfp: " << tracks.size() << std::endl;
      for (auto track : tracks) {

        std::vector<art::Ptr<anab::ParticleID>> pids = particleids_from_track.at(track.key());
        if(pids.size() == 0) std::cout << "[UBXSec] \t\t Zero ParticleID" << std::endl;
        if(pids.size() > 1) {
          std::cout << "[UBXSec] \t\t ParticleID vector is bigger than 1. Only one saved." << std::endl;
        }
        for (auto pid : pids) {
          if (!pid->PlaneID().isValid) continue;
          int planenum = pid->PlaneID().Plane;
          if (planenum < 0 || planenum > 2) continue;
          std::cout << "[UBXSec] \t\t ParticleID PIDA is " << pid->PIDA() << ", plane is " << planenum << std::endl;
          if (/*_is_signal && (ubxsec_event->slc_origin[slice] == 0 || ubxsec_event->slc_origin[slice] == 2) &&*/ planenum == 2) {
            if (pdg == 13) {
              _h_pida_muon->Fill(pid->PIDA());
              _h_pida_len_muon->Fill(pid->PIDA(), track->Length());
              if( pid->PIDA() > 0 && pid->PIDA() < 50. && _make_pida_csv) _csvfile << pid->PIDA() << "," << track->Length() << "," << "1" << std::endl;
            } else if (pdg == 2212) {
              _h_pida_proton->Fill(pid->PIDA());
              _h_pida_len_proton->Fill(pid->PIDA(), track->Length());
              if( pid->PIDA() > 0 && pid->PIDA() < 50. && _make_pida_csv) _csvfile << pid->PIDA() << "," << track->Length() << "," << "0" << std::endl;
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
 


    // MCS - Track direction, study cosmic direction
    // Take the muon candidate track for a TPCObject
    // with cosmic origin, and check the track direction
    // (up/down) from MCS result
    if (muon_cand_exists && tpcobj.GetOrigin() == ubana::kCosmicRay) {
      bool best_fwd = mcsfitresult_mu_v.at(candidate_track.key())->isBestFwd();
      bool down_track = candidate_track->Vertex().Y() > candidate_track->End().Y();
      //double deltall = mcsfitresult_mu_v.at(candidate_track.key())->deltaLogLikelihood();
      //double ratioll = _mcs_cosmic_track_fwdll / _mcs_cosmic_track_bwdll;  
      if (down_track && best_fwd) {               // Track is reco going down and mcs agrees (true for cosmic)
        _h_mcs_cosmic_track_direction->Fill(0);
        //_h_mcs_cosmic_track_direction_deltall->Fill(0., deltall);
        //_h_mcs_cosmic_track_direction_ratioll->Fill(0., ratioll);
        _mcs_cosmic_track_direction = 0;
        _mcs_cosmic_track_downll = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
        _mcs_cosmic_track_upll = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
      } else if (!down_track && best_fwd) {       // Track is reco going up and mcs agrees
        _h_mcs_cosmic_track_direction->Fill(1);
        //_h_mcs_cosmic_track_direction_deltall->Fill(1., deltall);
        //_h_mcs_cosmic_track_direction_ratioll->Fill(1., ratioll);
        _mcs_cosmic_track_direction = 1;
        _mcs_cosmic_track_downll = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
        _mcs_cosmic_track_upll = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
      } else if (down_track && !best_fwd) {       // Track is reco going down and mcs disagrees
        _h_mcs_cosmic_track_direction->Fill(1);
        //_h_mcs_cosmic_track_direction_deltall->Fill(1., deltall);
        //_h_mcs_cosmic_track_direction_ratioll->Fill(1., ratioll);
        _mcs_cosmic_track_direction = 2;
        _mcs_cosmic_track_downll = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
        _mcs_cosmic_track_upll = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
      } else if (!down_track && !best_fwd) {      // Track is reco going up and mcs disagrees (true for cosmic)
        _h_mcs_cosmic_track_direction->Fill(0);
        //_h_mcs_cosmic_track_direction_deltall->Fill(0., deltall);
        //_h_mcs_cosmic_track_direction_ratioll->Fill(0., ratioll);
        _mcs_cosmic_track_direction = 3;
        _mcs_cosmic_track_downll = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
        _mcs_cosmic_track_upll = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
      }
      _mcs_cosmic_track_direction_tree->Fill();
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
        ubxsec_event->is_swtriggered = (softwareTriggerHandle->passedAlgo(algoNames[0]) ? 1 : 0);
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
      ubxsec_event->numc_flash_spec.resize(geo->NOpDets());
      for (unsigned int i = 0; i < geo->NOpDets(); i++) {
        unsigned int opdet = geo->OpDetFromOpChannel(i);
        ubxsec_event->numc_flash_spec[opdet] = flash.PE(i);
      }
    }

    // MCFlash vs op activity
    bool opActivityInBeamSpill = false;
    // Check if there are recon beam flashed in the beam spill window
    for (auto reco_fls_time : ubxsec_event->beamfls_time) {
      if (reco_fls_time > _beam_spill_start && reco_fls_time < _beam_spill_end) {
         opActivityInBeamSpill = true;
       }
    }
    if (nuMcflash_h->size() > 0) {
      std::cout << "We have a neutrino MCFlash, and its time is: " << (*nuMcflash_h)[0].Time() << std::endl;    
    } else if (nuMcflash_h->size() == 0) {
      if(opActivityInBeamSpill) {
        std::cout << "No MCFlash but optical activity in the beam spill." << std::endl;
        ubxsec_event->no_mcflash_but_op_activity = true;
      }
    }

  }


  // POT
  art::Handle< sumdata::POTSummary > potsum_h;
  if(e.getByLabel(_potsum_producer, potsum_h))
    ubxsec_event->pot = potsum_h->totpot;
  else
    ubxsec_event->pot = 0.;

 




  // *********************
  // Event Selection
  // *********************

  _event_selection.SetEvent(ubxsec_event);

  size_t slice_index;
  std::string reason = "no_failure";
  std::map<std::string,bool> failure_map;
  bool is_selected = _event_selection.IsSelected(slice_index, failure_map);
  if (_debug) std::cout << "[UBXSec] >>>>>>>>>>>>>>>>>>>>>> Is Selected? " << (is_selected ? "YES" : "NO") << std::endl;
  bool first = true;
  if (_debug) {
    for (auto iter : failure_map) {
      std::cout << "[UBXSec] Cut: " << iter.first << "  >>>  " << (iter.second ? "PASSED" : "NOT PASSED") << std::endl;
      if (first && !iter.second) {
        reason = "fail_" + iter.first;
        first = false;
      } 
    }
  }

  std::cout << "[UBXSec] Selection Failed at Cut: " << reason << std::endl;

  ::ubana::SelectionResult selection_result;
  selection_result.SetSelectionType("numu_cc_inclusive");
  selection_result.SetFailureReason(reason);
  selection_result.SetCutFlowStatus(failure_map);

  // *********************
  // Save Event Selection Output in the Event
  // *********************

  if (!is_selected) {

    selection_result.SetSelectionStatus(false);

    ubxsec_event->is_selected = false;

    selectionResultVector->emplace_back(std::move(selection_result));
    //util::CreateAssn(*this, e, *selectionResultVector, tpcobj_v, *assnOutSelectionResultTPCObject);

  } else {

    selection_result.SetSelectionStatus(true);

    ubxsec_event->is_selected = true;

    // Grab the selected TPCObject
    std::vector<art::Ptr<ubana::TPCObject>> tpcobj_v;
    art::fill_ptr_vector(tpcobj_v, tpcobj_h);
    art::Ptr<ubana::TPCObject> tpcobj = tpcobj_v.at(slice_index);

    if (_debug) std::cout << "[UBXSec] >>>>>>>>>>>>>>>>>>>>>> Selected TPCObject with index " << slice_index << std::endl;
 
    // Prepare the tpcobj output
    std::vector<art::Ptr<ubana::TPCObject>>  out_tpcobj_v;
    out_tpcobj_v.resize(1);
    out_tpcobj_v.at(0) = tpcobj;

    selectionResultVector->emplace_back(std::move(selection_result));
    util::CreateAssn(*this, e, *selectionResultVector, out_tpcobj_v, *assnOutSelectionResultTPCObject);

    // For the TPCNeutrinoID Filter
    util::CreateAssn(*this, e, muon_candidate_track_per_slice_v.at(slice_index), neutrino_candidate_vertex_per_slice_v.at(slice_index), *vertexTrackAssociations);
    util::CreateAssn(*this, e, muon_candidate_pfparticle_per_slice_v.at(slice_index), neutrino_candidate_vertex_per_slice_v.at(slice_index), *vertexPFParticleAssociations);

  }

  if(_debug) std::cout << "[UBXSec] Filling tree now." << std::endl;
  _tree1->Fill();

  e.put(std::move(selectionResultVector));
  e.put(std::move(assnOutSelectionResultTPCObject));

  e.put(std::move(vertexTrackAssociations));
  e.put(std::move(vertexPFParticleAssociations));


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

    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    std::vector<double> sce_corr = SCE->GetPosOffsets(mclist[iList]->GetNeutrino().Nu().Vx(),
                                                      mclist[iList]->GetNeutrino().Nu().Vy(),
                                                      mclist[iList]->GetNeutrino().Nu().Vz());
    std::cout << "\t\t SCE correction in x, y, z = " << sce_corr.at(0) 
                                                     << ", " << sce_corr.at(1) 
                                                     << ", " << sce_corr.at(2) << std::endl;
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
