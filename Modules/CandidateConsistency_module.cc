////////////////////////////////////////////////////////////////////////
// Class:       CandidateConsistency
// Plugin Type: producer (art v2_05_00)
// File:        CandidateConsistency_module.cc
//
// Generated at Sun Dec 31 07:59:34 2017 by Marco Del Tutto using cetskelgen
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

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <memory>

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

//Algorithms include
#include "uboone/UBXSec/Algorithms/UBXSecHelper.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"
#include "uboone/UBXSec/Algorithms/StoppingMuonTaggerHelper.h"
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

class CandidateConsistency;


class CandidateConsistency : public art::EDProducer {
public:
  explicit CandidateConsistency(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CandidateConsistency(CandidateConsistency const &) = delete;
  CandidateConsistency(CandidateConsistency &&) = delete;
  CandidateConsistency & operator = (CandidateConsistency const &) = delete;
  CandidateConsistency & operator = (CandidateConsistency &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::string _tpcobject_producer;
  std::string _shower_producer;
  std::string _track_producer;
  double      _tolerance;
  double      _dqds_threshold;
  double      _linearity_threshold;
  bool        _debug;

  void ContainPoint(double *);

  ::art::ServiceHandle<geo::Geometry> geo;
  ::detinfo::DetectorProperties const* fDetectorProperties;

  ::cosmictag::CosmicTagManager _ct_manager;

};


CandidateConsistency::CandidateConsistency(fhicl::ParameterSet const & p)
{

  _tpcobject_producer  = p.get<std::string>("TPCObjectProducer",  "TPCObjectMaker::UBXSec");
  _shower_producer     = p.get<std::string>("ShowerProducer",     "pandoraNu::UBXSec");
  _track_producer      = p.get<std::string>("TrackProducer",      "pandoraNu::UBXSec");

  _tolerance           = p.get<double>("Tolerance", 5.);
  _dqds_threshold      = p.get<double>("DqDsThreshold", 60000);
  _linearity_threshold = p.get<double>("LinearityThreshold", 0.7);

  _debug               = p.get<bool>("DebugMode", true);

  _ct_manager.Configure(p.get<cosmictag::Config_t>("CosmicTagManager"));

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

  produces< std::vector<anab::CosmicTag>>();
  produces< art::Assns<ubana::TPCObject, anab::CosmicTag>>();

}

void CandidateConsistency::produce(art::Event & e)
{

  if(_debug) std::cout << "[CandidateConsistency] Starts." << std::endl;

  // Instantiate the output
  std::unique_ptr<std::vector<anab::CosmicTag>> cosmicTagVector (new std::vector<anab::CosmicTag>);
  std::unique_ptr<art::Assns<ubana::TPCObject, anab::CosmicTag>> assnOutCosmicTagTPCObject(new art::Assns<ubana::TPCObject, anab::CosmicTag>);

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

  // Collect tracks and map tracks->hits
  lar_pandora::TrackVector track_v;
  lar_pandora::TracksToHits tracks_to_hits;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _track_producer, track_v, tracks_to_hits);

  // Collect showers and map showers->hits
  lar_pandora::ShowerVector shower_v;
  lar_pandora::ShowersToHits showers_to_hits;
  lar_pandora::LArPandoraHelper::CollectShowers(e, _shower_producer, shower_v, showers_to_hits);


  for (size_t i = 0; i < tpcobj_v.size(); i++) {

    if (_debug) std::cout << "[StoppingMuonTagger] >>>>> TPCObject " << i << std::endl;

    bool consistency_failed = false;

    art::Ptr<ubana::TPCObject> tpcobj = tpcobj_v.at(i);

    std::vector<art::Ptr<recob::Track>>      tracks  = tpcobjToTrackAssns.at(tpcobj.key());
    std::vector<art::Ptr<recob::PFParticle>> pfps    = tpcobjToPFPAssns.at(tpcobj.key());
    std::vector<art::Ptr<recob::Shower>>     showers = tpcobjToShowerAssns.at(tpcobj.key());

    recob::Vertex vertex = tpcobj->GetVertex();
    double xyz[3];
    vertex.XYZ(xyz);
    TVector3 vtx_xyz(xyz[0], xyz[1], xyz[2]);

    // Collect objects that are close to the vertex
    std::vector<art::Ptr<recob::Track>>  selected_tracks;
    std::vector<TVector3>                selected_tracks_vertex;
    std::vector<art::Ptr<recob::Shower>> selected_showers;
    std::vector<TVector3>                selected_showers_vertex;

    // 1) Tracks
    for (auto t : tracks) {

      if ( (t->Vertex() - vtx_xyz).Mag() < _tolerance ) {
        selected_tracks.push_back(t);
        selected_tracks_vertex.push_back(t->Vertex());
      }

      if ( (t->End() - vtx_xyz).Mag() < _tolerance ) {
        selected_tracks.push_back(t);
        selected_tracks_vertex.push_back(t->End());
      }


    }

    // 2) Showers
    for (auto s : showers) {

      if ( (s->ShowerStart() - vtx_xyz).Mag() < _tolerance ) {
        selected_showers.push_back(s);
        selected_showers_vertex.push_back(s->ShowerStart());
      }

    }



    //
    // Now analyse the selected objects
    //

    // 1) Tracks 

    if(_debug) std::cout << "[CandidateConsistency] Working of Tracks" << std::endl;

    for (size_t i = 0; i < selected_tracks.size(); i++) {

      if(_debug) std::cout << "[CandidateConsistency] Track number " << i << std::endl;

      auto iter = tracks_to_hits.find(selected_tracks.at(i));
      if (iter == tracks_to_hits.end()) {
        std::cout << "[CandidateConsistency] Track in TPCObject not found by Pandora?!" << std::endl;
        continue;
      }
      auto hits = iter->second;

      // Take only collection plane hits
      std::vector<cosmictag::SimpleHit> simple_hit_v;
      for (auto h : hits) {

        if (h->View() != 2) continue;
        
        cosmictag::SimpleHit sh;
        sh.t = fDetectorProperties->ConvertTicksToX(h->PeakTime(), geo::PlaneID(0,0,2));
        sh.w = h->WireID().Wire * geo->WirePitch(geo::PlaneID(0,0,2));

        sh.plane = h->View();
        sh.integral = h->Integral();

        sh.time = h->PeakTime() / 4;
        sh.wire = h->WireID().Wire;

        simple_hit_v.emplace_back(sh);

      }

      // Emplacing simple hits to the manager
      cosmictag::SimpleCluster sc(simple_hit_v);
      _ct_manager.Emplace(std::move(sc));

      // Creating an approximate start hit
      double vertex[3] = {selected_tracks_vertex.at(i).X(), selected_tracks_vertex.at(i).Y(), selected_tracks_vertex.at(i).Z()};
      this->ContainPoint(vertex);
      double vertex_t = fDetectorProperties->ConvertXToTicks(vertex[0], geo::PlaneID(0,0,2))/4.;
      int vertex_w    = geo->NearestWire(vertex, 2);

      cosmictag::SimpleHit start;
      start.time = vertex_t;
      start.wire = vertex_w;
      start.plane = 2;
      _ct_manager.SetStartHit(std::move(start));

      // Running the cluster analyser
      bool passed = _ct_manager.Run();

      if (passed) {
        _ct_manager.PrintClusterStatus();
      } else {
        std::cout << "not passed" << std::endl;
        continue;
      }

      cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();
      std::vector<double> dqds_v = processed_cluster._dqds_slider;

      // Calulated the trunc median
      double dqds_trunc = UBXSecHelper::GetDqDxTruncatedMean(dqds_v);
      if(_debug) std::cout << "[CandidateConsistency] dqds_trunc is " << dqds_trunc << std::endl;

    } // selected tracks loop



    // 1) Showers 

    if(_debug) std::cout << "[CandidateConsistency] Working of Showers" << std::endl;
    
    for (size_t i = 0; i < selected_showers.size(); i++) {

      if(_debug) std::cout << "[CandidateConsistency] Shower number " << i << std::endl;

      auto iter = showers_to_hits.find(selected_showers.at(i));
      if (iter == showers_to_hits.end()) {
        std::cout << "[CandidateConsistency] Shower in TPCObject not found by Pandora?!" << std::endl;
        continue;
      }
      auto hits = iter->second;

      // Take only collection plane hits
      std::vector<cosmictag::SimpleHit> simple_hit_v;
      for (auto h : hits) {

        if (h->View() != 2) continue;
        
        cosmictag::SimpleHit sh;
        sh.t = fDetectorProperties->ConvertTicksToX(h->PeakTime(), geo::PlaneID(0,0,2));
        sh.w = h->WireID().Wire * geo->WirePitch(geo::PlaneID(0,0,2));

        sh.plane = h->View();
        sh.integral = h->Integral();

        sh.time = h->PeakTime() / 4;
        sh.wire = h->WireID().Wire;

        simple_hit_v.emplace_back(sh);

      }

      // Emplacing simple hits to the manager
      cosmictag::SimpleCluster sc(simple_hit_v);
      _ct_manager.Emplace(std::move(sc));

      // Creating an approximate start hit
      double vertex[3] = {selected_showers_vertex.at(i).X(), selected_showers_vertex.at(i).Y(), selected_showers_vertex.at(i).Z()};
      double vertex_t = fDetectorProperties->ConvertXToTicks(selected_showers_vertex.at(i).X(), geo::PlaneID(0,0,2))/4.;
      int vertex_w    = geo->NearestWire(vertex, 2);

      cosmictag::SimpleHit start;
      start.time = vertex_t;
      start.wire = vertex_w;
      start.plane = 2;
      _ct_manager.SetStartHit(std::move(start));

      // Running the cluster analyser
      bool passed = _ct_manager.Run();

      if (passed) {
        _ct_manager.PrintClusterStatus();
      } else {
        std::cout << "Not passed" << std::endl;
        continue;
      }

      cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();
      std::vector<double> dqds_v = processed_cluster._dqds_slider;

      // Calulated the trunc median
      double dqds_trunc = UBXSecHelper::GetDqDxTruncatedMean(dqds_v);
      if(_debug) std::cout << "[CandidateConsistency] dqds_trunc is " << dqds_trunc << std::endl;

      // Verify linearity
      std::vector<double> linearity_v = processed_cluster._linearity_v;
      bool linearity_failed = false;
      for (auto l : linearity_v) {
        if (l < _linearity_threshold)
          linearity_failed = true;
      }

      // If we have a shower, with low dqds, and whose hits
      // are not linear, then this is likely an electron
      if (dqds_trunc < _dqds_threshold && linearity_failed) {
        consistency_failed = true;
        if(_debug) std::cout << "[CandidateConsistency] Consistency Failed." << std::endl;
      }


    } // selected showers loop

    double cosmicScore = 0.;
    anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kNotTagged;
    if (consistency_failed) {
      cosmicScore = 1.;
      tag_id = anab::CosmicTagID_t::kGeometry_Y; // ... don't have a proper id for these
    }
    if (_debug) std::cout << "[CandidateConsistency] Cosmic Score is " << cosmicScore << std::endl;
    cosmicTagVector->emplace_back(endPt1, endPt2, cosmicScore, tag_id); 
    util::CreateAssn(*this, e, *cosmicTagVector, tpcobj, *assnOutCosmicTagTPCObject);

  } // tpcobject loop

  e.put(std::move(cosmicTagVector));
  e.put(std::move(assnOutCosmicTagTPCObject));

  if(_debug) std::cout << "[CandidateConsistency] Ends." << std::endl;
}

void CandidateConsistency::ContainPoint(double * point) {

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


DEFINE_ART_MODULE(CandidateConsistency)
