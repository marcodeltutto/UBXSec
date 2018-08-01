////////////////////////////////////////////////////////////////////////
// Class:       GeoCosmicTagger
// Plugin Type: producer (art v2_05_00)
// File:        GeoCosmicTagger_module.cc
//
// Generated at Wed Oct  4 13:25:06 2017 by Marco Del Tutto using cetskelgen
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

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/geo.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include <memory>

class GeoCosmicTagger;


class GeoCosmicTagger : public art::EDProducer {
public:
  explicit GeoCosmicTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GeoCosmicTagger(GeoCosmicTagger const &) = delete;
  GeoCosmicTagger(GeoCosmicTagger &&) = delete;
  GeoCosmicTagger & operator = (GeoCosmicTagger const &) = delete;
  GeoCosmicTagger & operator = (GeoCosmicTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::string _track_producer;
  std::string _pfp_producer;
  std::string _tpcobject_producer;

  double fDetHalfHeight;
  double fDetWidth;
  double fDetLength;

  double fTPCXBoundary;
  double fTPCYBoundary;
  double fTPCZBoundary;

  bool fDebug;

  const std::vector<float> endPt1 = {-9999., -9999., -9999.};
  const std::vector<float> endPt2 = {-9999., -9999., -9999.};

  void GetCosmicTag(double, double, double, double, double, double, anab::CosmicTagID_t &, double &);
};


GeoCosmicTagger::GeoCosmicTagger(fhicl::ParameterSet const & p) {

  _track_producer                 = p.get<std::string>("TrackProducer");
  _pfp_producer                   = p.get<std::string>("PFParticleProducer");
  _tpcobject_producer             = p.get<std::string>("TPCObjectProducer");

  fTPCXBoundary = p.get< float >("TPCXBoundary", 5);
  fTPCYBoundary = p.get< float >("TPCYBoundary", 5);
  fTPCZBoundary = p.get< float >("TPCZBoundary", 5);

  fDebug = p.get< bool >("DebugMode", true);

  if (fDebug) {
    std::cout << "[GeoCosmicTagger] fTPCXBoundary " << fTPCXBoundary << std::endl;
    std::cout << "[GeoCosmicTagger] fTPCYBoundary " << fTPCYBoundary << std::endl;
    std::cout << "[GeoCosmicTagger] fTPCZBoundary " << fTPCZBoundary << std::endl;
  }

  produces< std::vector<anab::CosmicTag>>();
  produces< art::Assns<anab::CosmicTag,   ubana::TPCObject>>();
}

void GeoCosmicTagger::produce(art::Event & e) {

  // Instantiate the output
  std::unique_ptr< std::vector< anab::CosmicTag>>                  cosmicTagVector        (new std::vector<anab::CosmicTag>);
  std::unique_ptr< art::Assns<anab::CosmicTag, ubana::TPCObject>>  assnOutCosmicTagTPCObj (new art::Assns<anab::CosmicTag,ubana::TPCObject>);

  auto const* geo = lar::providerFrom<geo::Geometry>();

  fDetHalfHeight = geo->DetHalfHeight();
  fDetWidth      = 2.*geo->DetHalfWidth();
  fDetLength     = geo->DetLength();

  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    std::cout << "[GeoCosmicTagger] Cannote locate ubana::TPCObject." << std::endl;
    e.put( std::move(cosmicTagVector)       );
    e.put( std::move(assnOutCosmicTagTPCObj));
    return;
  }
  std::vector<art::Ptr<ubana::TPCObject>> tpcobj_v;
  art::fill_ptr_vector(tpcobj_v, tpcobj_h);
  art::FindManyP<recob::PFParticle> pfps_from_tpcobj(tpcobj_h, e, _tpcobject_producer);

  // Get PFParticles from the Event
  art::Handle<std::vector<recob::PFParticle>> pfp_h;
  e.getByLabel(_pfp_producer, pfp_h);
  if (!pfp_h.isValid()) {
    std::cout << "[GeoCosmicTagger] Cannote locate recob::PFParticle." << std::endl;
    e.put( std::move(cosmicTagVector)       );
    e.put( std::move(assnOutCosmicTagTPCObj));
    return;
  }
  art::FindManyP<recob::SpacePoint> spacepoints_from_pfp(pfp_h, e, _pfp_producer);


  // TPCObject loop
  for (size_t i = 0; i < tpcobj_v.size(); i++) {

    art::Ptr<ubana::TPCObject> tpcobj = tpcobj_v.at(i);
    std::vector<art::Ptr<recob::PFParticle>> pfps = pfps_from_tpcobj.at(tpcobj.key());

    std::vector<recob::SpacePoint> sp_v;
    for (auto pfp : pfps) {
      std::vector<art::Ptr<recob::SpacePoint>> spacepoints = spacepoints_from_pfp.at(pfp.key());
      for (auto sp : spacepoints) {
        sp_v.push_back(*sp);
      }
    }
    if (fDebug) std::cout << "[GeoCosmicTagger] This TPCObject (" << i << ") has " << sp_v.size() << " SpacePoints." << std::endl;

    anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kNotTagged;
    double cosmicScore = 0;

    std::vector<art::Ptr<ubana::TPCObject>>  tpcobj_v;
    tpcobj_v.resize(1);
    tpcobj_v.at(0) = tpcobj;

    if (sp_v.size() == 0) {
      cosmicTagVector->emplace_back(endPt1, endPt2, cosmicScore, tag_id);
      util::CreateAssn(*this, e, *cosmicTagVector, tpcobj_v, *assnOutCosmicTagTPCObj);
 
      if (fDebug)  std::cout << "[GeoCosmicTagger]\t 0 spacepoints for this tpcobj" << std::endl;
      if (fDebug)  std::cout << "[GeoCosmicTagger]\t cosmicScore assigned is " << cosmicScore << std::endl;
      continue;
    }

    // *************
    // Sort SpacePints by X position (anode-closer first)
    // *************
/*    std::sort(sp_v.begin(), sp_v.end(),
              [](recob::SpacePoint a, recob::SpacePoint b) -> bool
              {
                const double* xyz_a = a.XYZ();
                const double* xyz_b = b.XYZ();
                return xyz_a[0] < xyz_b[0];
              });

    double trackEndPt1_X = sp_v.at(0).XYZ()[0];
    double trackEndPt2_X = sp_v.at(sp_v.size()-1).XYZ()[0];

    if (fDebug) {
      std::cout << "[GeoCosmicTagger]\t End Point 1 X: " << trackEndPt1_X << std::endl;
      std::cout << "[GeoCosmicTagger]\t End Point 2 X: " << trackEndPt2_X << std::endl;
    }
*/
    // *************
    // Sort SpacePints by Y position (higher first)
    // *************
    std::sort(sp_v.begin(), sp_v.end(),
              [](recob::SpacePoint a, recob::SpacePoint b) -> bool
              {
                const double* xyz_a = a.XYZ();
                const double* xyz_b = b.XYZ();
                return xyz_a[1] > xyz_b[1];
              });

    double trackEndPt1_Y = sp_v.at(0).XYZ()[1];
    double trackEndPt2_Y = sp_v.at(sp_v.size()-1).XYZ()[1];

    //for (auto s : sp_v) {
      //std::cout << "[GeoCosmicTagger]\t\t y sp: " << s.XYZ()[1]  << std::endl;
    //}

    if (fDebug) {
      std::cout << "[GeoCosmicTagger]\t End Point 1 Y: " << trackEndPt1_Y << std::endl;
      std::cout << "[GeoCosmicTagger]\t End Point 2 Y: " << trackEndPt2_Y << std::endl;
    }

    double trackEndPt1_X = sp_v.at(0).XYZ()[0];
    double trackEndPt2_X = sp_v.at(sp_v.size()-1).XYZ()[0];

    if (fDebug) {
      std::cout << "[GeoCosmicTagger]\t End Point 1 X: " << trackEndPt1_X << std::endl;
      std::cout << "[GeoCosmicTagger]\t End Point 2 X: " << trackEndPt2_X << std::endl;
    }

    // *************
    // Sort SpacePints by Z position (upstream first)
    // *************
/*    std::sort(sp_v.begin(), sp_v.end(),
              [](recob::SpacePoint a, recob::SpacePoint b) -> bool
              {
                const double* xyz_a = a.XYZ();
                const double* xyz_b = b.XYZ();
                return xyz_a[2] < xyz_b[2];
              });
*/
    double trackEndPt1_Z = sp_v.at(0).XYZ()[2];
    double trackEndPt2_Z = sp_v.at(sp_v.size()-1).XYZ()[2];

    if(fDebug) {
      std::cout << "[GeoCosmicTagger]\t End Point 1 Z: " << trackEndPt1_Z << std::endl;
      std::cout << "[GeoCosmicTagger]\t End Point 2 Z: " << trackEndPt2_Z << std::endl;
    }

    this->GetCosmicTag(trackEndPt1_X, trackEndPt2_X,
                       trackEndPt1_Y, trackEndPt2_Y,
                       trackEndPt1_Z, trackEndPt2_Z,
                       tag_id, cosmicScore);

    if(fDebug) std::cout << "[GeoCosmicTagger]\t tag_id is " << tag_id << std::endl;


    cosmicTagVector->emplace_back(endPt1, endPt2, cosmicScore, tag_id);
    util::CreateAssn(*this, e, *cosmicTagVector, tpcobj_v, *assnOutCosmicTagTPCObj);

    if (fDebug)  std::cout << "[GeoCosmicTagger]\t cosmicScore assigned is " << cosmicScore << std::endl;

  } // end TPCObject loop

  e.put( std::move(cosmicTagVector)       );
  e.put( std::move(assnOutCosmicTagTPCObj));
  
}






//____________________________________________________________________________________________
void GeoCosmicTagger::GetCosmicTag(double trackEndPt1_X, double trackEndPt2_X, 
                                   double trackEndPt1_Y, double trackEndPt2_Y, 
                                   double trackEndPt1_Z, double trackEndPt2_Z,
                                   anab::CosmicTagID_t &tag_id,
                                   double &cosmicScore) {


    tag_id = anab::CosmicTagID_t::kNotTagged;

    int isCosmic = 0;

    // In below we check entry and exit points. Note that a special case of a particle entering
    // and exiting the same surface is considered to be running parallel to the surface and NOT
    // entering and exiting.
    // Also, in what follows we make no assumptions on which end point is the "start" or
    // "end" of the track being considered.
    unsigned boundaryMask[] = {0,0};

    // Check x extents - note that uboone coordinaes system has x=0 at edge
    // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
    // Also note that, in theory, any cosmic ray entering or exiting the X surfaces will have presumably
    // been removed already by the checking of "out of time" hits... but this will at least label
    // neutrino interaction tracks which exit through the X surfaces of the TPC
    if      (fDetWidth - trackEndPt1_X < fTPCXBoundary) boundaryMask[0] = 0x1;
    else if (            trackEndPt1_X < fTPCXBoundary) boundaryMask[0] = 0x2;

    if      (fDetWidth - trackEndPt2_X < fTPCXBoundary) boundaryMask[1] = 0x1;
    else if (            trackEndPt2_X < fTPCXBoundary) boundaryMask[1] = 0x2;

    // Check y extents (note coordinate system change)
    // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
    if      (fDetHalfHeight - trackEndPt1_Y < fTPCYBoundary) boundaryMask[0] = 0x10;
    else if (fDetHalfHeight + trackEndPt1_Y < fTPCYBoundary) boundaryMask[0] = 0x20;

    if      (fDetHalfHeight - trackEndPt2_Y < fTPCYBoundary) boundaryMask[1] = 0x10;
    else if (fDetHalfHeight + trackEndPt2_Y < fTPCYBoundary) boundaryMask[1] = 0x20;

    // Check z extents
    // Note this counts the case where the track enters and exits the same surface as a "1", not a "2"
    if      (fDetLength - trackEndPt1_Z < fTPCZBoundary) boundaryMask[0] = 0x100;
    else if (             trackEndPt1_Z < fTPCZBoundary) boundaryMask[0] = 0x200;

    if      (fDetLength - trackEndPt2_Z < fTPCZBoundary) boundaryMask[1] = 0x100;
    else if (             trackEndPt2_Z < fTPCZBoundary) boundaryMask[1] = 0x200;

    unsigned trackMask = boundaryMask[0] | boundaryMask[1];
    int      nBitsSet(0);

    for(int idx=0; idx<12; idx++) if (trackMask & (0x1<<idx)) nBitsSet++;

    // This should check for the case of a track which is both entering and exiting
    // but we consider entering and exiting the z boundaries to be a special case (should it be?)
    if(nBitsSet > 1)
    {
      if ((trackMask & 0x300) != 0x300)
      {
        isCosmic = 2;
        if(fDebug) std::cout << "[GeoCosmicTagger]\t Geometry x,y." << std::endl;
        if      ((trackMask &  0x3) ==  0x3)                tag_id = anab::CosmicTagID_t::kGeometry_XX;
        else if ((trackMask & 0x30) == 0x30)                tag_id = anab::CosmicTagID_t::kGeometry_YY;
        else if ((trackMask &  0x3) && (trackMask &  0x30)) tag_id = anab::CosmicTagID_t::kGeometry_XY;
        else if ((trackMask &  0x3) && (trackMask & 0x300)) tag_id = anab::CosmicTagID_t::kGeometry_XZ;
        else                                                tag_id = anab::CosmicTagID_t::kGeometry_YZ;
      }
      // This is the special case of track which appears to enter/exit z boundaries
      else
      {
        if(fDebug) std::cout << "[GeoCosmicTagger]\t Geometry z." << std::endl;
        isCosmic = 3;
        tag_id   = anab::CosmicTagID_t::kGeometry_ZZ;
      }
    }
    // This looks for track which enters/exits a boundary but has other endpoint in TPC
    else if (nBitsSet > 0)
    {
      if(fDebug) std::cout << "[GeoCosmicTagger]\t Geometry. Enter or Exit but not both." << std::endl;
      isCosmic = 4 ;
      if      (trackMask &   0x3) tag_id = anab::CosmicTagID_t::kGeometry_X;
      else if (trackMask &  0x30) tag_id = anab::CosmicTagID_t::kGeometry_Y;
      else if (trackMask & 0x300) tag_id = anab::CosmicTagID_t::kGeometry_Z;
    }

    cosmicScore = isCosmic > 0 ? 1. : 0.;

    // Handle special cases
    if      (isCosmic == 3) cosmicScore = 0.4;   // Enter/Exit at opposite Z boundaries
    else if (isCosmic == 4) cosmicScore = 0.5;   // Enter or Exit but not both

    return;
}

DEFINE_ART_MODULE(GeoCosmicTagger)
