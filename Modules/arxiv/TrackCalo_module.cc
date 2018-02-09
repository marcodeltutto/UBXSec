  //////////////////////////////////////////////////
  //
  // TrackCalo class
  //
  // maddalena.antonello@lngs.infn.it
  // ornella.palamara@lngs.infn.it
  // ART port echurch@fnal.gov
  //  This algorithm is designed to perform the calorimetric reconstruction 
  //  of the 3D reconstructed tracks
  ////////////////////////////////////////////////////////////////////////

extern "C" {
  #include <sys/types.h>
  #include <sys/stat.h>
}
  #include <vector>
  #include <string>
  #include <math.h>
  #include <algorithm>
  #include <iostream>
  #include <fstream>

  #include "lardata/DetectorInfoServices/LArPropertiesService.h"
  #include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
  #include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

  #include "lardataobj/RecoBase/Hit.h"
  #include "lardataobj/RecoBase/SpacePoint.h"
  #include "lardataobj/RecoBase/Track.h"
  #include "lardataobj/RecoBase/TrackHitMeta.h"
  #include "lardataobj/AnalysisBase/TrackCalo.h"
  #include "lardataobj/AnalysisBase/T0.h"
  #include "lardata/Utilities/AssociationUtil.h"
  #include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
  #include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
  #include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
  #include "larcorealg/Geometry/PlaneGeo.h"
  #include "larcorealg/Geometry/WireGeo.h"

  // ROOT includes
  #include <TROOT.h>
  #include <TFile.h>
  #include <TTree.h>
  #include <TBranch.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TMath.h>
  #include <TGraph.h>
  #include <TF1.h>
  #include <TVector3.h>

  // Framework includes
  #include "art/Framework/Core/EDProducer.h"
  #include "art/Framework/Core/ModuleMacros.h" 
  #include "canvas/Persistency/Common/FindManyP.h"
  #include "art/Framework/Principal/Event.h" 
  #include "fhiclcpp/ParameterSet.h" 
  #include "art/Framework/Principal/Handle.h" 
  #include "canvas/Persistency/Common/Ptr.h" 
  #include "canvas/Persistency/Common/PtrVector.h" 
  #include "art/Framework/Services/Registry/ServiceHandle.h" 
  #include "art/Framework/Services/Optional/TFileService.h" 
  #include "art/Framework/Services/Optional/TFileDirectory.h" 
  #include "messagefacility/MessageLogger/MessageLogger.h" 


class TrackCalo : public art::EDProducer {

public:

  explicit TrackCalo(fhicl::ParameterSet const& pset); 
  virtual ~TrackCalo();

  void beginJob(); 
      //    void endJob();

  void produce(art::Event& evt);

private:

  std::string fTrackModuleLabel;
  std::string fSpacePointModuleLabel;
  std::string fT0ModuleLabel;
  bool fUseArea;
  bool fFlipTrack_dQdx; //flip track direction if significant rise of dQ/dx at the track start

  int fnsps;
  std::vector<int>    fwire;
  std::vector<double> ftime;
  std::vector<double> fstime;
  std::vector<double> fetime;
  std::vector<double> fMIPs;
  std::vector<double> fdQdx;
  std::vector<double> fdEdx;
  std::vector<double> fResRng;
  std::vector<double> fpitch;
  std::vector<TVector3> fXYZ;

protected: 

}; // class TrackCalo



  //-------------------------------------------------
TrackCalo::TrackCalo(fhicl::ParameterSet const& pset)
{

  fTrackModuleLabel = pset.get< std::string >("TrackModuleLabel", "pandoraNu::UBXSec");
  fSpacePointModuleLabel = pset.get< std::string >("SpacePointModuleLabel", "pandoraNu::UBXSec");
  fT0ModuleLabel = pset.get< std::string >("T0ModuleLabel", "pandoraNu::UBXSec");
  fUseArea = pset.get< bool >("UseArea", true);
  fFlipTrack_dQdx = pset.get< bool >("FlipTrack_dQdx", false);

  produces< std::vector<anab::TrackCalo>              >();
  produces< art::Assns<recob::Track, anab::TrackCalo> >();
}

  //-------------------------------------------------
TrackCalo::~TrackCalo()
{

}

  //-------------------------------------------------
void TrackCalo::beginJob()
{
  return;
}

  //------------------------------------------------------------------------------------//
void TrackCalo::produce(art::Event& evt)
{ 
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

  // Get Geometry
  art::ServiceHandle<geo::Geometry> geom;

  // channel quality
  lariov::ChannelStatusProvider const& channelStatus
  = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

  size_t nplanes = geom->Nplanes();

  //create anab::TrackCalo objects and make association with recob::Track
  std::unique_ptr< std::vector<anab::TrackCalo> > calorimetrycol(new std::vector<anab::TrackCalo>);
  std::unique_ptr< art::Assns<recob::Track, anab::TrackCalo> > assn(new art::Assns<recob::Track, anab::TrackCalo>);

  //art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit>        fmht(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); //this has more information about hit-track association, only available in PMA for now

  for (size_t trkIter = 0; trkIter < tracklist.size(); ++trkIter){   

      decltype(auto) larEnd = tracklist[trkIter]->Trajectory().End();
      //store track directional cosines
      double trackCosStart[3]={0.,0.,0.};
      double trackCosEnd[3]={0.,0.,0.};
      tracklist[trkIter]->Direction(trackCosStart,trackCosEnd);

      // Some variables for the hit
      float time;          //hit time at maximum
      float stime;         //hit start time 
      float etime;         //hit end time 
      uint32_t     channel = 0;//channel number
      unsigned int cstat   = 0;    //hit cryostat number 
      unsigned int tpc     = 0;    //hit tpc number 
      unsigned int wire    = 0;   //hit wire number 
      unsigned int plane   = 0;  //hit plane number

      std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(trkIter);
      
      
      std::vector< std::vector<unsigned int> > hits(nplanes);

      art::FindManyP<recob::SpacePoint> fmspts(allHits, evt, fSpacePointModuleLabel);
      for (size_t ah = 0; ah< allHits.size(); ++ah){
        hits[allHits[ah]->WireID().Plane].push_back(ah);
      }
      //get hits in each plane
      for (size_t ipl = 0; ipl < nplanes; ++ipl){//loop over all wire planes

        geo::PlaneID planeID;//(cstat,tpc,ipl);

        fwire.clear();
        ftime.clear();
        fstime.clear();
        fetime.clear();
        fMIPs.clear();
        fdQdx.clear();
        fdEdx.clear();
        fpitch.clear();
        fResRng.clear();
        fXYZ.clear();

        double Kin_En = 0.;
        double Trk_Length = 0.;
        std::vector<double> vdEdx;
        std::vector<double> vresRange;
        std::vector<double> vdQdx;
        std::vector<double> deadwire; //residual range for dead wires
        std::vector<TVector3> vXYZ;

        //range of wire signals
        unsigned int wire0 = 100000;
        unsigned int wire1 = 0;
        double PIDA = 0;
        int nPIDA = 0;

        // determine track direction. Fill residual range array
        bool GoingDS = true;
        // find the track direction by comparing US and DS charge BB
        double USChg = 0;
        double DSChg = 0;
        // temp array holding distance betweeen space points
        std::vector<double> spdelta;
        //int nht = 0; //number of hits
        fnsps = 0; //number of space points
        std::vector<double> ChargeBeg;
        std::stack<double> ChargeEnd;     

        // find track pitch
        double fTrkPitch = 0;
        for (size_t itp = 0; itp < tracklist[trkIter]->NumberTrajectoryPoints(); ++itp){
          const TVector3& pos = tracklist[trkIter]->LocationAtPoint(itp);
          const double Position[3] = { pos.X(), pos.Y(), pos.Z() };
          geo::TPCID tpcid = geom->FindTPCAtPosition ( Position );
          if (tpcid.isValid) {
            try{
              fTrkPitch = lar::util::TrackPitchInView(*tracklist[trkIter], geom->Plane(ipl).View(), itp);
            }
            catch( cet::exception &e){
              mf::LogWarning("TrackCalo") << "caught exception " 
              << e << "\n setting pitch (C) to "
              << util::kBogusD;
              fTrkPitch = 0;
            }
            break;
          }
        }

        // find the separation between all space points
        double xx = 0.,yy = 0.,zz = 0.;

        //save track 3d points
        std::vector<double> trkx;
        std::vector<double> trky;
        std::vector<double> trkz;
        std::vector<double> trkw;
        std::vector<double> trkx0;
        for (size_t i = 0; i<hits[ipl].size(); ++i){
  	//Get space points associated with the hit
         std::vector< art::Ptr<recob::SpacePoint> > sptv = fmspts.at(hits[ipl][i]);
         for (size_t j = 0; j < sptv.size(); ++j){

  	  double t = allHits[hits[ipl][i]]->PeakTime() - TickT0; // Want T0 here? Otherwise ticks to x is wrong?
  	  double x = detprop->ConvertTicksToX(t, allHits[hits[ipl][i]]->WireID().Plane, allHits[hits[ipl][i]]->WireID().TPC, allHits[hits[ipl][i]]->WireID().Cryostat);
  	  double w = allHits[hits[ipl][i]]->WireID().Wire;
  	  if (TickT0){
       trkx.push_back(sptv[j]->XYZ()[0]-detprop->ConvertTicksToX(TickT0, allHits[hits[ipl][i]]->WireID().Plane, allHits[hits[ipl][i]]->WireID().TPC, allHits[hits[ipl][i]]->WireID().Cryostat));
     }
     else{
       trkx.push_back(sptv[j]->XYZ()[0]);
     }
     trky.push_back(sptv[j]->XYZ()[1]);
     trkz.push_back(sptv[j]->XYZ()[2]);
     trkw.push_back(w);
     trkx0.push_back(x);
   }
 }
        for (size_t ihit = 0; ihit < hits[ipl].size(); ++ihit){//loop over all hits on each wire plane

  	//std::cout<<ihit<<std::endl;

         if (!planeID.isValid){
           plane = allHits[hits[ipl][ihit]]->WireID().Plane;
           tpc   = allHits[hits[ipl][ihit]]->WireID().TPC;
           cstat = allHits[hits[ipl][ihit]]->WireID().Cryostat;
           planeID.Cryostat = cstat;
           planeID.TPC = tpc;
           planeID.Plane = plane;
           planeID.isValid = true;
         }

         wire = allHits[hits[ipl][ihit]]->WireID().Wire;
  	time = allHits[hits[ipl][ihit]]->PeakTime(); // What about here? T0 
  	stime = allHits[hits[ipl][ihit]]->PeakTimeMinusRMS();
  	etime = allHits[hits[ipl][ihit]]->PeakTimePlusRMS();
  	
  	double charge = allHits[hits[ipl][ihit]]->PeakAmplitude();
  	if (fUseArea) charge = allHits[hits[ipl][ihit]]->Integral();




  	//get 3d coordinate and track pitch for the current hit
  	//not all hits are associated with space points, the method uses neighboring spacepts to interpolate
  	double xyz3d[3];
  	double pitch;
    bool fBadhit = false;
    if (fmthm.isValid()){
      auto vhit = fmthm.at(trkIter);
      auto vmeta = fmthm.data(trkIter);
      for (size_t ii = 0; ii<vhit.size(); ++ii){
        if (vhit[ii].key() == allHits[hits[ipl][ihit]].key()){
          if (!tracklist[trkIter]->HasValidPoint(vmeta[ii]->Index())){
            fBadhit = true;
            continue;
          }
          double angleToVert = geom->WireAngleToVertical(vhit[ii]->View(), vhit[ii]->WireID().TPC, vhit[ii]->WireID().Cryostat) - 0.5*::util::pi<>();
          const TVector3& dir = tracklist[trkIter]->DirectionAtPoint(vmeta[ii]->Index());
          double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() + std::cos(angleToVert)*dir.Z());
          if (cosgamma){
            pitch = geom->WirePitch(0,1,0)/cosgamma;
          }
          else{
            pitch = 0;
          }
          TVector3 loc = tracklist[trkIter]->LocationAtPoint(vmeta[ii]->Index());
          xyz3d[0] = loc.X();
          xyz3d[1] = loc.Y();
          xyz3d[2] = loc.Z();
          break;
        }
      }
    }
    else
      GetPitch(allHits[hits[ipl][ihit]], trkx, trky, trkz, trkw, trkx0, xyz3d, pitch, TickT0);

    if (fBadhit) continue;
  	if (xyz3d[2]<-100) continue; //hit not on track
  	if (pitch<=0) pitch = fTrkPitch;
  	if (!pitch) continue;

    if(fnsps == 0) {
      xx = xyz3d[0];
      yy = xyz3d[1];
      zz = xyz3d[2];
      spdelta.push_back(0);
    } else {
      double dx = xyz3d[0] - xx;
      double dy = xyz3d[1] - yy;
      double dz = xyz3d[2] - zz;
      spdelta.push_back(sqrt(dx*dx + dy*dy + dz*dz));
      Trk_Length += spdelta.back();
      xx = xyz3d[0];
      yy = xyz3d[1];
      zz = xyz3d[2];
    }

    ChargeBeg.push_back(charge);
    ChargeEnd.push(charge);

    double MIPs = charge;
    double dQdx = MIPs/pitch;

  
    ++fnsps;
  }
  if (!fnsps){
  	//std::cout << "Adding the aforementioned positions..." << std::endl;
  	calorimetrycol->push_back(anab::TrackCalo(util::kBogusD,
      vdEdx,
      vdQdx,
      vresRange,
      deadwire,
      util::kBogusD,
      fpitch,
      vXYZ,
      planeID));
  	util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);
  	continue;
  }
  for (int isp = 0; isp<fnsps; ++isp){
  	if (isp>3) break;
  	USChg += ChargeBeg[isp];
  }
  int countsp = 0;
  while (!ChargeEnd.empty()){
  	if (countsp>3) break;
  	DSChg += ChargeEnd.top();
  	ChargeEnd.pop();
  	++countsp;
  }
        // Going DS if charge is higher at the end
  GoingDS = (DSChg > USChg) || (!fFlipTrack_dQdx);
        // determine the starting residual range and fill the array
  fResRng.resize(fnsps);
  if(GoingDS) {
    fResRng[fnsps - 1] = spdelta[fnsps - 1] / 2;
    for(int isp = fnsps - 2; isp > -1; isp--) {
      fResRng[isp] = fResRng[isp+1] + spdelta[isp+1];
    }
  } else {
    fResRng[0] = spdelta[1] / 2;
    for(int isp = 1; isp < fnsps; isp++) {
      fResRng[isp] = fResRng[isp-1] + spdelta[isp];
    }
  }


  double Ai = -1;
  for (int i = 0; i < fnsps; ++i){//loop over all 3D points
    vresRange.push_back(fResRng[i]);
    vdEdx.push_back(fdEdx[i]);
    vdQdx.push_back(fdQdx[i]);
    vXYZ.push_back(fXYZ[i]);
  	if (i!=0 && i!= fnsps-1){//ignore the first and last point
  	  // Calculate PIDA 
      Ai = fdEdx[i] * pow(fResRng[i],0.42);
      nPIDA++;
      PIDA += Ai;
    }
  }//end looping over 3D points
  if(nPIDA > 0) {
    PIDA = PIDA / (double)nPIDA;
  } 
  else {
    PIDA = -1;
  }

  
        //std::cout << "Adding at the end but still same fXYZ" << std::endl;
  calorimetrycol->push_back(anab::TrackCalo(Kin_En,
    vdEdx,
    vdQdx,
    vresRange,
    deadwire,
    Trk_Length,
    fpitch,
    vXYZ,
    planeID));
  util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);

    }//end looping over planes
  }//end looping over tracks

  evt.put(std::move(calorimetrycol));
  evt.put(std::move(assn));

  return;
}

  

DEFINE_ART_MODULE(TrackCalo)


