////////////////////////////////////////////////////////////////////////
// Class:       StoppingMuonTagger
// Plugin Type: producer (art v2_05_00)
// File:        StoppingMuonTagger_module.cc
//
// Generated at Fri Oct  6 15:55:20 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

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
#include "lardataobj/RecoBase/Shower.h"
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
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

// LArSoft include
#include "uboone/UBFlashFinder/PECalib.h"
#include "larsim/MCCheater/BackTracker.h"

//Algorithms include
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"
#include "uboone/UBXSec/Algorithms/StoppingMuonTaggerHelper.h"

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

  //::art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;
  //::art::ServiceHandle<detinfo::DetectorPropertiesService> det;
  ::detinfo::DetectorProperties const* fDetectorProperties;

  ::ubana::FiducialVolume _fiducial_volume;

  ::ubana::StoppingMuonTaggerHelper _helper;

  ::trkf::TrajectoryMCSFitter _mcs_fitter;

  /// Takes a pointer to a point, and if outside the dector puts it at the det border
  void ContainPoint(double * point);

  /// Takes a recob::Tracks and returns true if is a stopping muon
  bool IsStopMuMCS(art::Ptr<recob::Track> t);

  std::string _tpcobject_producer;
  std::string _pfp_producer;
  std::string _cluster_producer;
  std::string _track_producer;

  double _coplanar_cut = 5.;
  double _mcs_delta_ll_cut = -5.;

  bool _debug;

  TH1D * _h_nstopmu;
  TH1D * _h_stopmu_type;
};


StoppingMuonTagger::StoppingMuonTagger(fhicl::ParameterSet const & p) 
  : _mcs_fitter(p.get< fhicl::ParameterSet >("MCSFitter")) {


  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();

  _tpcobject_producer = p.get<std::string>("TPCObjectProducer",  "TPCObjectMaker::UBXSec");
  _pfp_producer       = p.get<std::string>("PFParticleProducer", "pandoraNu::UBXSec");
  _cluster_producer   = p.get<std::string>("ClusterProducer",    "pandoraNu::UBXSec");
  _track_producer     = p.get<std::string>("TrackProducer",    "pandoraNu::UBXSec");

  _coplanar_cut       = p.get<double>("CoplanarCut",   5.);
  _mcs_delta_ll_cut   = p.get<double>("MCSDeltaLLCut", -5.);

  _helper.Configure(p.get<fhicl::ParameterSet>("AlgorithmConfiguration"));

  _helper.PrintConfig();

  _debug = p.get<bool>("DebugMode", false);

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 

  art::ServiceHandle<art::TFileService> fs;
  _h_nstopmu = fs->make<TH1D>("h_nstopmu", ";Stopping Muons Per Event;", 10, 0, 10);
  _h_stopmu_type = fs->make<TH1D>("h_stopmu_type", "Stopping Muon Chategories;;", 10, 0, 10);


  produces< std::vector<anab::CosmicTag>>();
  produces< art::Assns<anab::CosmicTag,   recob::Track>>();
  produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();

}



void StoppingMuonTagger::produce(art::Event & e) {

  if (_debug) std::cout <<"[StoppingMuonTagger] Starts." << std::endl;

  // Instantiate the output
  std::unique_ptr<std::vector<anab::CosmicTag>> cosmicTagVector (new std::vector<anab::CosmicTag>);
  std::unique_ptr<art::Assns<anab::CosmicTag, recob::Track>> assnOutCosmicTagTrack(new art::Assns<anab::CosmicTag, recob::Track>);
  std::unique_ptr<art::Assns<recob::PFParticle, anab::CosmicTag>> assnOutCosmicTagPFParticle(new art::Assns<recob::PFParticle, anab::CosmicTag>);


  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    //std::cout << "[StoppingMuonTagger] Cannote locate ubana::TPCObject." << std::endl;
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
          if (_debug) std::cout << "delta x is " << deltax << std::endl;
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

    // First exclude spacepoints outside the tpc
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

    // Now get the highest point
    std::sort(sp_v.begin(), sp_v.end(),
              [](art::Ptr<recob::SpacePoint> a, art::Ptr<recob::SpacePoint> b) -> bool
              {
                const double *xyz_a = a->XYZ();
                const double *xyz_b = b->XYZ();
                return xyz_a[1] > xyz_b[1];
              });

    if (sp_v.size() == 0) {
      std::cout << "Not enough spacepoints." << std::endl;
      continue;
    }

    const double *highest_point_c = sp_v.at(0)->XYZ();
    double highest_point[3] = {highest_point_c[0], highest_point_c[1], highest_point_c[2]};
    //raw::ChannelID_t ch = geo->NearestChannel(highest_point, 2);
    this->ContainPoint(highest_point);
    int highest_w = geo->NearestWire(highest_point, 2) ;//* geo->WirePitch(geo::PlaneID(0,0,2));
    double highest_t = fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0,0,2))/4.;//highest_point[0];
    //size_t time = fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0,0,2));
    if (_debug) std::cout << "[StoppingMuonTagger] Highest point: wire: " << geo->NearestWire(highest_point, 2) 
                       << ", time: " << fDetectorProperties->ConvertXToTicks(highest_point[0], geo::PlaneID(0,0,2)) 
                       << std::endl;

    /*
    double time2cm = 0.1;
    double triggerOffset = 3200;
    double planeOrigin[3] = {0.0, -0.3, -0.6};
    double offset_collection = triggerOffset * time2cm - planeOrigin[2];
    size_t time = highest_point[0]/time2cm + offset_collection;  
    */

    if (_debug) std::cout << "[StoppingMuonTagger] Now create simple hit vector, size " << hit_v.size() << std::endl;

    // Create SimpleHit vector
    ubana::SimpleHitVector shit_v;
    for (auto h : hit_v) {
      ubana::SimpleHit shit;
      shit.t = fDetectorProperties->ConvertTicksToX(h->PeakTime(), geo::PlaneID(0,0,2));
      shit.w = h->WireID().Wire * geo->WirePitch(geo::PlaneID(0,0,2));

      shit.plane = h->View();
      shit.integral = h->Integral();

      shit.time = h->PeakTime() / 4;
      shit.wire = h->WireID().Wire;

      //if (_debug) std::cout << "Emplacing hit with time " << shit.time*4 << ", and wire " << shit.wire << ", on plane " << shit.plane << std::endl;
      shit_v.emplace_back(shit);
    }

    if (_debug) std::cout << "[StoppingMuonTagger] Simple hit vector size " << shit_v.size() << std::endl;

    // Clear the algorithm
    if (_debug) std::cout << "[StoppingMuonTagger] Clear algo" << std::endl;
    _helper.Clear();

    if (collection_coplanar) {
      _helper.SetMaxAllowedHitDistance(10);
    } else {
      _helper.SetMaxAllowedHitDistance(6); 
    }

    // Set the highest point and the FV
    _helper.SetFiducialVolume(_fiducial_volume);
    _helper.SetVertexPoint(highest_point);

    // Emplace to the helper
    if (_debug) std::cout << "[StoppingMuonTagger] Now emplace hits" << std::endl;
    _helper.Emplace(shit_v);

    // Only collection plane hits
    if (_debug) std::cout << "[StoppingMuonTagger] Now filter hits by plane" << std::endl;
    size_t n_hits = _helper.FilterByPlane(2);

    if (n_hits == 0) continue;

    // Set the start hit
    if (_debug) std::cout << "[StoppingMuonTagger] Now set start hit" << std::endl;
    _helper.SetStartHit(highest_t, highest_w, 2); 

    // Order hits
    if (_debug) std::cout << "[StoppingMuonTagger] Now order hits" << std::endl;
    n_hits = _helper.OrderHits();

    if (n_hits == 0) continue;

    // Filter hist is belonging to same wire
    n_hits = _helper.FilterOnSingleWire();

    if (n_hits == 0) continue;

    // dQds
    if (_debug) std::cout << "[StoppingMuonTagger] Now calculate dqds" << std::endl; 
    _helper.CalculatedQds();    

    // dQds Slider
    if (_debug) std::cout << "[StoppingMuonTagger] Now perform window slider" << std::endl;
    bool status = _helper.PerformdQdsSlider();

    if (!status)
      continue;

    // Linearity
    if (_debug) std::cout << "[StoppingMuonTagger] Now calculate local linearity" << std::endl;
    _helper.CalculateLocalLinearity();

     _helper.PrintOnFile(i);

    // Make a decision
    if (_debug) std::cout << "[StoppingMuonTagger] Now make a decision" << std::endl;
    bool result = _helper.MakeDecision(ubana::kAlgoMichel);

    if (_debug) std::cout << "[StoppingMuonTagger] Is stopping muon (michel)? " << (result ? "YES" : "NO") << std::endl;


    bool result_bragg = _helper.MakeDecision(ubana::kAlgoBragg);
    if (_debug) std::cout << "[StoppingMuonTagger] Is stopping muon (bragg)? " << (result_bragg ? "YES" : "NO") << std::endl;

    bool result_simplemip = _helper.MakeDecision(ubana::kAlgoSimpleMIP);
    if (_debug) std::cout << "[StoppingMuonTagger] Is simple MIP? " << (result_simplemip ? "YES" : "NO") << std::endl;
    //bool result_curvature = _helper.MakeDecision(ubana::kAlgoCurvature);
    //if (_debug) std::cout << "[StoppingMuonTagger] Is stopping muon (curvature)? " << (result_curvature ? "YES" : "NO") << std::endl;

    /* Also try with the MCS fitter
    bool result_mcs = false;
    if (primary_track_v.size() != 0) { 
      result_mcs = this->IsStopMuMCS(primary_track_v.at(0));
    }

    if (_debug) std::cout << "[StoppingMuonTagger] MCS thinks " << (result_mcs ? "is" : "is not") << " a stopping muon" << std::endl;
    */

    double cosmicScore = 0.;
    anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kNotTagged;
    if (result || result_bragg || result_simplemip/*|| result_mcs*/) {
      cosmicScore = 1.;
      tag_id = anab::CosmicTagID_t::kGeometry_Y;
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



bool StoppingMuonTagger::IsStopMuMCS(art::Ptr<recob::Track> t) {

  TVector3 track_start = t->Vertex();
  TVector3 track_end   = t->End();

  // Check it goes donwards
  if (track_end.Y() > track_start.Y())
    std::swap(track_start, track_end);

  // Check the upper point is outside the FV
  bool vtx_contained = _fiducial_volume.InFV(track_start);

  if (vtx_contained)
    return false;

  recob::MCSFitResult result = _mcs_fitter.fitMcs(*t);

  double fwd_ll = result.fwdLogLikelihood();
  double bwd_ll = result.bwdLogLikelihood(); 

  if (_debug) std::cout <<"[StoppingMuonTagger] FWD " << fwd_ll << ", BWD " << bwd_ll << std::endl;

  double delta_ll;

  if (track_start.Y() > track_end.Y())
    delta_ll = fwd_ll - bwd_ll;
  else
    delta_ll = bwd_ll - fwd_ll;

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
