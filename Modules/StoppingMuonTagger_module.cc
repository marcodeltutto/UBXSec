////////////////////////////////////////////////////////////////////////
// Class:       StoppingMuonTagger
// Plugin Type: producer (art v2_05_00)
// File:        StoppingMuonTagger_module.cc
//
// Generated at Fri Oct  6 15:55:20 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class StoppingMuonTagger
 *
 * \ingroup UBXSec
 *
 * \brief Art producer module that tags stopping muons
 * 
 *
 * \author Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
 *
 * \version producer (art v2_05_00)
 *
 * \date 2017/03/10
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Fri Oct  6 15:55:20 2017
 *
 */

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
#include "lardata/Utilities/AssociationUtil.h"

#include <memory>
#include <limits>

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "uboone/UBXSec/DataTypes/FlashMatch.h"
#include "uboone/UBXSec/DataTypes/MCGhost.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

// LArSoft include
#include "uboone/UBFlashFinder/PECalib.h"
#include "larsim/MCCheater/BackTracker.h"

//Algorithms include
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"
#include "uboone/UBXSec/HitCosmicTag/Base/DataTypes.h"
#include "uboone/UBXSec/HitCosmicTag/Base/CosmicTagManager.h"
#include "uboone/UBXSec/HitCosmicTag/Algorithms/StopMuMichel.h"
#include "uboone/UBXSec/HitCosmicTag/Algorithms/StopMuBragg.h"
#include "uboone/UBXSec/HitCosmicTag/Algorithms/CosmicSimpleMIP.h"

// Root include
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"

const std::vector<float> endPt1 = {-9999., -9999., -9999.};
const std::vector<float> endPt2 = {-9999., -9999., -9999.};

class StoppingMuonTagger;


class StoppingMuonTagger : public art::EDProducer {
public:
  explicit StoppingMuonTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StoppingMuonTagger(StoppingMuonTagger const &) = delete;
  StoppingMuonTagger(StoppingMuonTagger &&) = delete;
  StoppingMuonTagger & operator = (StoppingMuonTagger const &) = delete;
  StoppingMuonTagger & operator = (StoppingMuonTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  ::cosmictag::CosmicTagManager _ct_manager;

  ::art::ServiceHandle<geo::Geometry> geo;
  ::detinfo::DetectorProperties const* fDetectorProperties;

  ::ubana::FiducialVolume _fiducial_volume;

  ::trkf::TrajectoryMCSFitter _mcs_fitter;
  ::recob::MCSFitResult _result;

  /// Takes a pointer to a point, and, if outside the dector, puts it at the detector border
  void ContainPoint(double * point);

  /// Takes a recob::Track and returns true if is a stopping muon according to MCS algo
  bool IsStopMuMCS(art::Ptr<recob::Track> t, double & delta_ll);

  std::string _tpcobject_producer;
  std::string _pfp_producer;
  std::string _cluster_producer;
  std::string _track_producer;

  double _coplanar_cut = 5.; ///< If a track start X and end X difference is below this value, is consiered to be coplanar to a collection plane wire 

  bool _use_mcs; ///< If true, uses mcs fit to reject cosmics
  double _mcs_delta_ll_cut = -5.; ///< Cut on MCS delta loglikelihood (used if _use_mcs is true)

  bool _debug; ///< Debug flag
  bool _create_tree; ///< If true, creates a tree with mcs delta ll info, to make plots

  TH1D * _h_nstopmu;
  TH1D * _h_stopmu_type;

  TTree* _tree1;
  int _run, _subrun, _event;
  int _origin, _origin_extra;
  bool _fv;
  double _delta_ll;
  double _length;
};


StoppingMuonTagger::StoppingMuonTagger(fhicl::ParameterSet const & p) 
  : _mcs_fitter(p.get< fhicl::ParameterSet >("MCSFitter")) {


  _ct_manager.Configure(p.get<cosmictag::Config_t>("CosmicTagManager"));

  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  std::cout << "[StoppingMuonTagger] FV: " << std::endl;
  _fiducial_volume.PrintConfig();

  _tpcobject_producer = p.get<std::string>("TPCObjectProducer",  "TPCObjectMaker::UBXSec");
  _pfp_producer       = p.get<std::string>("PFParticleProducer", "pandoraNu::UBXSec");
  _cluster_producer   = p.get<std::string>("ClusterProducer",    "pandoraNu::UBXSec");
  _track_producer     = p.get<std::string>("TrackProducer",      "pandoraNu::UBXSec");
 
  _use_mcs            = p.get<bool>("UseMCS", false);
  _mcs_delta_ll_cut   = p.get<double>("MCSDeltaLLCut", -5.);

  _coplanar_cut       = p.get<double>("CoplanarCut",   5.);

  _debug = p.get<bool>("DebugMode", false);
  _create_tree = p.get<bool>("CreateTree", true);

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 

  art::ServiceHandle<art::TFileService> fs;
  _h_nstopmu = fs->make<TH1D>("h_nstopmu", ";Stopping Muons Per Event;", 10, 0, 10);
  _h_stopmu_type = fs->make<TH1D>("h_stopmu_type", "Stopping Muon Chategories;;", 10, 0, 10);


  if (_create_tree) {
    _tree1 = fs->make<TTree>("tree","");
    _tree1->Branch("run",           &_run,                 "run/I");
    _tree1->Branch("subrun",        &_subrun,              "subrun/I");
    _tree1->Branch("event",         &_event,               "event/I");
    _tree1->Branch("origin",        &_origin,              "origin/I");
    _tree1->Branch("origin_extra",  &_origin_extra,        "origin_extra/I");
    _tree1->Branch("fv",            &_fv,                  "fv/O");
    _tree1->Branch("delta_ll",      &_delta_ll,            "delta_ll/D");
    _tree1->Branch("length",        &_length,              "length/D");
  }

  produces< std::vector<anab::CosmicTag>>();
  produces< art::Assns<anab::CosmicTag,   recob::Track>>();
  produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();

}



void StoppingMuonTagger::produce(art::Event & e) {

  if (_debug) std::cout <<"[StoppingMuonTagger] Starts." << std::endl;

  if(_create_tree) {
    _run = e.id().run();
    _subrun = e.id().subRun();
    _event = e.id().event();
  }

  // Instantiate the output
  std::unique_ptr<std::vector<anab::CosmicTag>> cosmicTagVector (new std::vector<anab::CosmicTag>);
  std::unique_ptr<art::Assns<anab::CosmicTag, recob::Track>> assnOutCosmicTagTrack(new art::Assns<anab::CosmicTag, recob::Track>);
  std::unique_ptr<art::Assns<recob::PFParticle, anab::CosmicTag>> assnOutCosmicTagPFParticle(new art::Assns<recob::PFParticle, anab::CosmicTag>);


  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    throw cet::exception(__PRETTY_FUNCTION__) << "Cannote locate ubana::TPCObject." << std::endl;
  }

  art::FindManyP<recob::Track>      tpcobjToTrackAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::PFParticle> tpcobjToPFPAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Shower>     tpcobjToShowerAssns(tpcobj_h, e, _tpcobject_producer);

  std::vector<art::Ptr<ubana::TPCObject>> tpcobj_v;
  art::fill_ptr_vector(tpcobj_v, tpcobj_h);

  // Collect pfparticles and map pfp->spacepoint
  lar_pandora::PFParticleVector pfp_v;
  lar_pandora::PFParticlesToSpacePoints pfps_to_spacepoints;
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, pfp_v, pfps_to_spacepoints);

  // Collect pfparticles and map pfp->cluster
  lar_pandora::PFParticlesToClusters pfps_to_clusters;
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, pfp_v, pfps_to_clusters);

  // Collect clusters and map cluster->hit
  lar_pandora::ClusterVector cluster_v;
  lar_pandora::ClustersToHits clusters_to_hits;
  lar_pandora::LArPandoraHelper::CollectClusters(e, _cluster_producer, cluster_v, clusters_to_hits);

  // Collect tracks and map pfp->track
  lar_pandora::TrackVector track_v;
  lar_pandora::PFParticlesToTracks pfps_to_tracks;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _track_producer, track_v, pfps_to_tracks);


  int n_stopmu = 0;

  for (size_t i = 0; i < tpcobj_v.size(); i++) {

    if (_debug) std::cout << "[StoppingMuonTagger] >>>>> TPCObject " << i << std::endl;

    art::Ptr<ubana::TPCObject> tpcobj = tpcobj_v.at(i);

    std::vector<art::Ptr<recob::Track>>      tracks  = tpcobjToTrackAssns.at(tpcobj.key());
    std::vector<art::Ptr<recob::PFParticle>> pfps    = tpcobjToPFPAssns.at(tpcobj.key());
    std::vector<art::Ptr<recob::Shower>>     showers = tpcobjToShowerAssns.at(tpcobj.key());

    if (_debug) std::cout << "[StoppingMuonTagger] Origin is " << tpcobj->GetOrigin() << std::endl;
    if (_debug) std::cout << "[StoppingMuonTagger] Origin extra is " << tpcobj->GetOriginExtra() << std::endl;

    bool ignore_this = false;

    if (!e.isRealData()) {

      if (tpcobj->GetOrigin() != ubana::TPCObjectOrigin::kCosmicRay &&
          tpcobj->GetOriginExtra() != ubana::TPCObjectOriginExtra::kStoppingMuon) {

        if (_debug) std::cout << "[StoppingMuonTagger] Found TPCObject representing a stopping cosmic muon." << std::endl;
        n_stopmu++;

        if (tracks.size() == 1 && showers.size() == 0)          // One track - bin 0
          _h_stopmu_type->Fill(0);
        else if (tracks.size() == 0 && showers.size() == 1)     // One shower - bin 1
          _h_stopmu_type->Fill(1);
        else if (tracks.size() == 1 && showers.size() == 1)     // One track One shower - bin 2
          _h_stopmu_type->Fill(2);
        else if (tracks.size() == 2 && showers.size() == 0)     // Two track Zero shower - bin 3
          _h_stopmu_type->Fill(3);
        else if (tracks.size() == 0 && showers.size() == 2)     // Zero track Two shower - bin 4
          _h_stopmu_type->Fill(4);
        else                                                    // Others
          _h_stopmu_type->Fill(5);

      }

    }

    if (_debug) { 
      std::cout <<"[StoppingMuonTagger]\t Number of tracks for this TPCObject:  " << tracks.size()  << std::endl;
      std::cout <<"[StoppingMuonTagger]\t Number of showers for this TPCObject: " << showers.size() << std::endl;
      std::cout <<"[StoppingMuonTagger]\t Number of PFPs for this TPCObject:    " << pfps.size()    << std::endl;
   }


    // Get all the SpacePoints from the PFPs (this will be used to look for 
    // the highest point in the TPCObject). Also get all the clusters for this
    // TPCObject, and use the clusters to go to the hits. Collect all
    // these hits in a vector.
    // When we look at the primary pfp (a cosmic by fiat in pandoraCosmic), 
    // get the track and then check that is not coplanar with wires in the 
    // collection plane.
  
    // These two are needed to make the associations later
    art::Ptr<recob::PFParticle> primary_pfp;
    std::vector<art::Ptr<recob::Track>> primary_track_v;

    std::vector<art::Ptr<recob::SpacePoint>> sp_v;
    sp_v.clear();

    std::vector<art::Ptr<recob::Hit>> hit_v; 
    hit_v.clear();

    if (pfps.size() == 0) {
      ignore_this = true;    
      continue;
    }

    bool collection_coplanar = false;

    for (auto p : pfps) {
      auto iter = pfps_to_spacepoints.find(p);
      if (iter == pfps_to_spacepoints.end()) { 
        continue;
      }
      sp_v.reserve(sp_v.size() + iter->second.size());
      sp_v.insert(sp_v.end(), iter->second.begin(), iter->second.end());

      // Find clusters first ...
      auto iter2 = pfps_to_clusters.find(p);
      if (iter2 == pfps_to_clusters.end()) { 
        continue;
      }

      // ... then find hits
      for (auto c : iter2->second) {
        auto iter3 = clusters_to_hits.find(c);
        if (iter3 == clusters_to_hits.end()) {
          std::cout << "[StoppingMuonTagger] Cluster in TPCObject not found by pandora !?" << std::endl;
          throw std::exception(); 
        }

        hit_v.reserve(hit_v.size() + iter3->second.size());
        hit_v.insert(hit_v.end(), iter3->second.begin(), iter3->second.end());
      } 

      if (p->IsPrimary() && !lar_pandora::LArPandoraHelper::IsNeutrino(p)) {

        auto iter4 = pfps_to_tracks.find(p);
        if (iter4 == pfps_to_tracks.end()) {
          std::cout << "[StoppingMuonTagger] PFParticle in TPCObject not found by pandora !?" << std::endl;
          //throw std::exception();
          ignore_this = true;
          continue;
        }

        primary_pfp = p;
        primary_track_v = iter4->second;

        if (primary_track_v.size() == 0) { 
          ignore_this = true;
          continue; 
        } else {
          double deltax = primary_track_v.at(0)->Vertex().Z() - primary_track_v.at(0)->End().Z();
          if (_debug) std::cout << "[StoppingMuonTagger] Delta x is " << deltax << std::endl;
          collection_coplanar = std::abs(primary_track_v.at(0)->Vertex().Z() - primary_track_v.at(0)->End().Z()) 
                                < _coplanar_cut;
        }
      }
    }

    if (hit_v.size() == 0) {
      continue;
    }

    if (ignore_this || !primary_pfp) {
      continue;
    }

    if (_debug) std::cout << "[StoppingMuonTagger] Primary PFP is " << primary_pfp->Self() << std::endl;

    //
    // First exclude spacepoints outside the tpc
    //
    std::vector<art::Ptr<recob::SpacePoint>> temp;
    ::geoalgo::AABox tpcvol(0., (-1.)*geo->DetHalfHeight(), 
                            0., geo->DetHalfWidth()*2, 
                            geo->DetHalfHeight(), geo->DetLength());

    for (auto s : sp_v) {
      const double *xyz = s->XYZ();
      ::geoalgo::Vector point (xyz[0], xyz[1], xyz[2]);
      if (tpcvol.Contain(point)) {
        temp.push_back(s);
      }
    }
    sp_v = temp;

    //
    // Now get the highest point
    //
    std::sort(sp_v.begin(), sp_v.end(),
              [](art::Ptr<recob::SpacePoint> a, art::Ptr<recob::SpacePoint> b) -> bool
              {
                const double *xyz_a = a->XYZ();
                const double *xyz_b = b->XYZ();
                return xyz_a[1] > xyz_b[1];
              });

    if (sp_v.size() == 0) {
      if (_debug) std::cout << "[StoppingMuonTagger] Not enough spacepoints." << std::endl;
      continue;
    }

    const double *highest_point_c = sp_v.at(0)->XYZ();
    double highest_point[3] = {highest_point_c[0], highest_point_c[1], highest_point_c[2]};
    this->ContainPoint(highest_point);

    // Creating an approximate start hit on plane 2
    int highest_w = geo->NearestWire(highest_point, 2) ;//* geo->WirePitch(geo::PlaneID(0,0,2));
    double highest_t = fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0,0,2))/4.;//highest_point[0];
    if (_debug) std::cout << "[StoppingMuonTagger] Highest point: wire: " << geo->NearestWire(highest_point, 2) 
                       << ", time: " << fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0,0,2)) 
                       << std::endl;
    cosmictag::SimpleHit start_highest;
    start_highest.time = highest_t;
    start_highest.wire = highest_w;
    start_highest.plane = 2;

    // Creating an approximate start hit on plane 1 (used if collection coplanar)
    highest_w = geo->NearestWire(highest_point, 1) ;
    highest_t = fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0,0,1))/4.;
    if (_debug) std::cout << "[StoppingMuonTagger] Highest point: wire: " << geo->NearestWire(highest_point, 1) 
                       << ", time: " << fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0,0,1)) 
                       << std::endl;
    cosmictag::SimpleHit start_highest_plane1;
    start_highest_plane1.time = highest_t;
    start_highest_plane1.wire = highest_w;
    start_highest_plane1.plane = 1;


    //
    // Now get the point outfv
    //
    art::Ptr<recob::Track> trk = primary_track_v.at(0);
    bool start_fv = _fiducial_volume.InFV(trk->Vertex());
    bool end_fv   = _fiducial_volume.InFV(trk->End());
    double point_outfv[3] = {-999, -999, -999};
    if (!start_fv) {point_outfv[0] = trk->Vertex().X(); point_outfv[1] = trk->Vertex().Y(); point_outfv[2] = trk->Vertex().Z();}
    if (!end_fv) {point_outfv[0] = trk->End().X(); point_outfv[1] = trk->End().Y(); point_outfv[2] = trk->End().Z();}
    this->ContainPoint(point_outfv);
    int outfv_w = geo->NearestWire(point_outfv, 2);
    double outfv_t = fDetectorProperties->ConvertXToTicks(point_outfv[0], geo::PlaneID(0,0,2))/4.;
    if (_debug) std::cout << "[StoppingMuonTagger] OutFV point: wire: " << geo->NearestWire(point_outfv, 2) 
                       << ", time: " << fDetectorProperties->ConvertXToTicks(point_outfv[0], geo::PlaneID(0,0,2)) 
                       << std::endl;
    // Creating an approximate start hit
    cosmictag::SimpleHit start_outfv;
    start_outfv.time = outfv_t;
    start_outfv.wire = outfv_w;
    start_outfv.plane = 2;

    /*
    bool use_outfv_point = false;
     
    TVector3 pt_1 (highest_point[0], highest_point[1], highest_point[2]);
    TVector3 pt_2 (point_outfv[0], point_outfv[1], point_outfv[2]);
    if ((pt_1-pt_2).Mag() > 5 && point_outfv[0] != -999)
      use_outfv_point = true;
    */

    if (_debug) std::cout << "[StoppingMuonTagger] Now create simple hit vector, size " << hit_v.size() << std::endl;

    
    //std::cout << "Wire inclination for plane 1 (should give 60 degrees): " << geo->WireAngleToVertical(geo::View_t::kV) - 0.5*::util::pi<>()<< std::endl;
    //std::cout << "Wire pitch for plane 1 (should give 3mm): " << geo->WirePitch(geo::PlaneID(0,0,1)) << std::endl;

    //
    // Create SimpleHit vector with hits in collection plane only
    // 
    std::vector<cosmictag::SimpleHit> simple_hit_v;
    std::vector<cosmictag::SimpleHit> simple_hit_v_plane1;
    for (auto h : hit_v) {


      cosmictag::SimpleHit sh;

      sh.t = fDetectorProperties->ConvertTicksToX(h->PeakTime(), geo::PlaneID(0,0,h->View()));
      sh.w = h->WireID().Wire * geo->WirePitch(geo::PlaneID(0,0,h->View()));

      sh.plane = h->View();
      sh.integral = h->Integral();

      sh.time = h->PeakTime() / 4;
      sh.wire = h->WireID().Wire; 

      if (h->View() == 2) {
        simple_hit_v.emplace_back(sh);
      }

      if (h->View() == 1) {
        simple_hit_v_plane1.emplace_back(sh);
      }

    }


    if (_debug) std::cout << "[StoppingMuonTagger] Simple hit vector size " << simple_hit_v.size() << std::endl;


    bool is_cosmic = false;

    //
    // Running with highest point as start hit
    //

    _ct_manager.Reset();

    // Emplacing simple hits to the manager
    // Generally use plane 2, but use plane 1 
    // if trak is collection coplanar
    cosmictag::SimpleCluster sc(simple_hit_v);
    if (collection_coplanar) {
      if (_debug) std::cout << "[StoppingMuonTagger] Collection coplanar, using hits on plane 1" << std::endl;
      sc._s_hit_v = simple_hit_v_plane1;
    }
    _ct_manager.Emplace(std::move(sc));

    // Emplace the start hit
    if (!collection_coplanar) {
      _ct_manager.SetStartHit(std::move(start_highest));
    }
    else {
      if (_debug) std::cout << "[StoppingMuonTagger] Collection coplanar, setting start hit from plane 1" << std::endl;
      _ct_manager.SetStartHit(std::move(start_highest_plane1));
    }

    //if (collection_coplanar) {
      //if (_debug) std::cout << "[StoppingMuonTagger] This object is collection coplanar" << std::endl;
      //_ct_manager.CollectionCoplanar(true);
    //}

    // Running the cluster analyser
    bool passed = _ct_manager.Run();

    bool ct_result_michel = false;
    bool ct_result_bragg = false;
   // bool ct_result_simplemip = false;

    if (passed) {

      if (_debug) _ct_manager.PrintClusterStatus();

      //if (_debug) _ct_manager.PrintOnFile(i);


      cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();

      // Michel algo
      if (_debug) ((cosmictag::StopMuMichel*)(_ct_manager.GetCustomAlgo("StopMuMichel")))->PrintConfig();
      ct_result_michel = ((cosmictag::StopMuMichel*)(_ct_manager.GetCustomAlgo("StopMuMichel")))->IsStopMuMichel(processed_cluster);
      if(_debug) std::cout << "[StoppingMuonTagger] Is stopping muon (michel)? " << (ct_result_michel ? "YES" : "NO") << std::endl;



      // Bragg algo
      bool vtx_in_fv = _fiducial_volume.InFV(highest_point);
      if (_debug) ((cosmictag::StopMuBragg*)(_ct_manager.GetCustomAlgo("StopMuBragg")))->PrintConfig();
      ct_result_bragg = ((cosmictag::StopMuBragg*)(_ct_manager.GetCustomAlgo("StopMuBragg")))->IsStopMuBragg(processed_cluster) && !vtx_in_fv;
      if(_debug) std::cout << "[StoppingMuonTagger] Is stopping muon (bragg)? " << (ct_result_bragg ? "YES" : "NO") << std::endl;


      // CosmicSimpleMIP
      //((cosmictag::CosmicSimpleMIP*)(_ct_manager.GetCustomAlgo("CosmicSimpleMIP")))->PrintConfig();
      //ct_result_simplemip = ((cosmictag::CosmicSimpleMIP*)(_ct_manager.GetCustomAlgo("CosmicSimpleMIP")))->IsCosmicSimpleMIP(processed_cluster);
      //if(_debug) std::cout << "[StoppingMuonTagger] Is simple MIP? " << (ct_result_simplemip ? "YES" : "NO") << std::endl;

    }

    if (ct_result_michel || ct_result_bragg /*|| ct_result_simplemip*/)
      is_cosmic = true;




    //
    // Also try with the MCS fitter
    //

    bool result_mcs = false;
    double delta_ll;
    if (primary_track_v.size() != 0) { 
      result_mcs = this->IsStopMuMCS(primary_track_v.at(0), delta_ll);
    }

    if (_debug) std::cout << "[StoppingMuonTagger] MCS thinks " << (result_mcs ? "is" : "is not") << " a stopping muon." << std::endl;
    
    if(_create_tree) {
      _origin = tpcobj->GetOrigin();
      _origin_extra = tpcobj->GetOriginExtra();
      _delta_ll = delta_ll;
      _length = primary_track_v.at(0)->Length();

      _fv = false;

      if (!e.isRealData()) {

        art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
        std::vector<art::Ptr<simb::MCTruth> > mclist;
        if (e.getByLabel("generator",mctruthListHandle))
          art::fill_ptr_vector(mclist, mctruthListHandle);

        int iList = 0; // 1 nu int per spill

        if (mclist[iList]->Origin() == simb::kBeamNeutrino) {

          double truth_nu_vtx[3] = {mclist[iList]->GetNeutrino().Nu().Vx(),
                                    mclist[iList]->GetNeutrino().Nu().Vy(),
                                    mclist[iList]->GetNeutrino().Nu().Vz()};
          if (_fiducial_volume.InFV(truth_nu_vtx)) 
            _fv = true;
        }
      }

      _tree1->Fill();
    }


    if (result_mcs && _use_mcs) is_cosmic = true;




    //
    // Now Emplace CosmiTag
    //

    double cosmicScore = 0.;
    anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kNotTagged;
    if (is_cosmic) {
      cosmicScore = 1.;
      if (result_mcs && _use_mcs) cosmicScore = 1.77;
      tag_id = anab::CosmicTagID_t::kGeometry_Y; // ... don't have a proper id for these
    }
    if (_debug) std::cout << "[StoppingMuonTagger] Cosmic Score is " << cosmicScore << std::endl;
    cosmicTagVector->emplace_back(endPt1, endPt2, cosmicScore, tag_id); 
    util::CreateAssn(*this, e, *cosmicTagVector, primary_track_v, *assnOutCosmicTagTrack);
    util::CreateAssn(*this, e, *cosmicTagVector, primary_pfp, *assnOutCosmicTagPFParticle);
    
  } // TPCObject loop

  _h_nstopmu->Fill(n_stopmu);

  e.put(std::move(cosmicTagVector));
  e.put(std::move(assnOutCosmicTagTrack));
  e.put(std::move(assnOutCosmicTagPFParticle));

  if (_debug) std::cout <<"[StoppingMuonTagger] Ends." << std::endl;

}



bool StoppingMuonTagger::IsStopMuMCS(art::Ptr<recob::Track> t, double & delta_ll) {

  delta_ll = -9999;

  TVector3 track_start = t->Vertex();
  TVector3 track_end   = t->End();

  // Check it goes donwards
  if (track_end.Y() > track_start.Y())
    std::swap(track_start, track_end);

  // Check the upper point is outside the FV
  bool vtx_contained = _fiducial_volume.InFV(track_start);

  if (vtx_contained) {
    if (_debug) std::cout <<"[StoppingMuonTagger] Vertex is contained, MCS will quit." << std::endl;
    return false;
  }

  _result = _mcs_fitter.fitMcs(*t);

  double fwd_ll = _result.fwdLogLikelihood();
  double bwd_ll = _result.bwdLogLikelihood(); 

  if (_debug) std::cout <<"[StoppingMuonTagger] FWD " << fwd_ll << ", BWD " << bwd_ll << std::endl;

  if (track_start.Y() > track_end.Y())
    delta_ll = fwd_ll - bwd_ll;
  else
    delta_ll = bwd_ll - fwd_ll;

  if (delta_ll == -9999)
    return false;

  if (_debug) std::cout <<"[StoppingMuonTagger] DELTA " << delta_ll << ", cut at " << _mcs_delta_ll_cut << std::endl;
  if (delta_ll < _mcs_delta_ll_cut)
    return true;
  else
    return false;

}



void StoppingMuonTagger::ContainPoint(double * point) {

  double x = point[0];
  double y = point[1];
  double z = point[2];

  double e = std::numeric_limits<double>::epsilon();

  if (x < 0. + e)
    point[0] = 0. + e;
  if (x > 2.*geo->DetHalfWidth() - e)
    point[0] = 2.*geo->DetHalfWidth() - e;

  if (y < -geo->DetHalfWidth() + e)
    point[1] = -geo->DetHalfWidth() + e;
  if (y > geo->DetHalfWidth() - e)
    point[1] = geo->DetHalfWidth() - e;

  if (z < 0. + e) 
    point[2] = 0.+ e;
  if (z > geo->DetLength() - e)
    point[2] = geo->DetLength() - e;

}


DEFINE_ART_MODULE(StoppingMuonTagger)
