////////////////////////////////////////////////////////////////////////
// Class:       FindDeadRegions
// Plugin Type: service (art v2_05_00)
// File:        FindDeadRegions_service.cc
//
// Generated at Sat Apr 15 10:25:29 2017 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "TH2F.h"

#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

struct BoundaryWire {
  unsigned int wire_num;
  float y_start;
  float z_start;
  float y_end;
  float z_end;
  bool isLowWire;
};

class FindDeadRegions;

class FindDeadRegions {
public:
  explicit FindDeadRegions(fhicl::ParameterSet const & p, art::ActivityRegistry & areg);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  /// Returns true is the passed point is close to a dead region given a tolerance considering two planes only
  virtual bool NearDeadReg2P(float yVal, float zVal, float tolerance);

  /// Returns true is the passed point is close to a dead region given a tolerance considering all three planes
  virtual bool NearDeadReg3P(float yVal, float zVal, float tolerance);

  /// Return a root 2D histogram (y v.s. z) containing the detector dead regions considering two planes only
  virtual TH2F* GetDeadRegionHisto2P();

  /// Return a root 2D histogram (y v.s. z) containing the detector dead regions considering all three planes
  virtual TH2F* GetDeadRegionHisto3P();

private:

  void LoadBWires();

  std::vector<BoundaryWire> BWires_U; ///< Contains list of wires marking the boundaries of dead regions (U plane)
  std::vector<BoundaryWire> BWires_V; ///< Contains list of wires marking the boundaries of dead regions (V plane)
  std::vector<BoundaryWire> BWires_Y; ///< Contains list of wires marking the boundaries of dead regions (Y plane)

};


FindDeadRegions::FindDeadRegions(fhicl::ParameterSet const & p, art::ActivityRegistry & areg)
// :
// Initialize member data here.
{

  /*
  ::art::ServiceHandle<geo::Geometry> geo;
  double reco_nu_vtx[3];
  raw::ChannelID_t ch = geo->NearestChannel(reco_nu_vtx, 2);
  std::cout << ch << std::endl;
*/
}


void FindDeadRegions::LoadBWires() {


  //::art::ServiceHandle<geo::Geometry> geo;
 // geo::GeometryCore const* geo = lar::providerFrom<geo::Geometry>();
  //raw::ChannelID_t ch = 0;

  //double reco_nu_vtx[3]; 
  //raw::ChannelID_t ch = geo->NearestChannel(reco_nu_vtx, 2);
  //std::cout << ch << std::endl;
  //  std::vector< geo::WireID > wire_v = geo->ChannelToWire(ch);

/*  std::cout << "wire_v has size: " << wire_v.size() << std::endl;

  geo::WireGeo wire_g = geo->Wire (wire_v[0]);
  double xyz[3];
  wire_g.GetStart (xyz);
  std::cout << "channel start, x="<<xyz[0]<<", y="<<xyz[1]<<", z="<<xyz[2] << std::endl;
  wire_g.GetEnd (xyz);
  std::cout << "channel end, x="<<xyz[0]<<", y="<<xyz[1]<<", z="<<xyz[2] << std::endl;
*/

  std::ifstream geofile;
  geofile.open("ChannelWireGeometry_v2.txt");

  std::string string_channel;
  std::string string_plane;
  std::string string_wire;
  std::string string_sx;
  std::string string_sy;
  std::string string_sz;
  std::string string_ex;
  std::string string_ey;
  std::string string_ez;

  unsigned int channel;
  unsigned int plane;
  unsigned int wire;
  float sx;
  float sy;
  float sz;
  float ex;
  float ey;
  float ez;

  std::vector<unsigned int> channelVec;
  std::vector<unsigned int> planeVec;
  std::vector<unsigned int> wireVec;
  std::vector<float> sxVec;
  std::vector<float> syVec;
  std::vector<float> szVec;
  std::vector<float> exVec;
  std::vector<float> eyVec;
  std::vector<float> ezVec;

  std::string dummy;

  getline(geofile,dummy);
  getline(geofile,dummy);
  getline(geofile,dummy);
  getline(geofile,dummy);
  getline(geofile,dummy);
  getline(geofile,dummy);
  getline(geofile,dummy);
  getline(geofile,dummy);
  getline(geofile,dummy);

  while(geofile >> string_channel)
  {
    geofile >> string_plane;
    geofile >> string_wire;
    geofile >> string_sx;
    geofile >> string_sy;
    geofile >> string_sz;
    geofile >> string_ex;
    geofile >> string_ey;
    geofile >> string_ez;

    channel = atoi(string_channel.c_str());
    plane = atoi(string_plane.c_str());
    wire = atoi(string_wire.c_str());
    sx = atof(string_sx.c_str());
    sy = atof(string_sy.c_str());
    sz = atof(string_sz.c_str());
    ex = atof(string_ex.c_str());
    ey = atof(string_ey.c_str());
    ez = atof(string_ez.c_str());

    channelVec.push_back(channel);
    planeVec.push_back(plane);
    wireVec.push_back(wire);
    sxVec.push_back(sx);
    syVec.push_back(sy);
    szVec.push_back(sz);
    exVec.push_back(ex);
    eyVec.push_back(ey);
    ezVec.push_back(ez);
  }

  geofile.close();

  std::ifstream chanstatfile;
  chanstatfile.open("ChanStatus.txt");

  std::string string_CSchannel;
  std::string string_CSstatus;

  unsigned int CSchannel;
  bool CSstatus;

  std::vector<unsigned int> CSchannelVec;
  std::vector<bool> CSstatusVec;

  while(chanstatfile >> string_CSchannel)
  {
    chanstatfile >> string_CSstatus;

    CSchannel = atoi(string_CSchannel.c_str());
    CSstatus = atoi(string_CSstatus.c_str());

    CSchannelVec.push_back(CSchannel);
    CSstatusVec.push_back(CSstatus);
  }

  chanstatfile.close();

  bool isGoodChannel;

  isGoodChannel = true;
  for(int i = 0; i < 2400; i++) {
    if((CSstatusVec.at(i) == false) && (isGoodChannel == true)) {
      isGoodChannel = false;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = true;

      BWires_U.push_back(BWire);
    }
    else if((CSstatusVec.at(i) == true) && (isGoodChannel == false)) {
      isGoodChannel = true;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i-1);
      BWire.y_start = syVec.at(i-1);
      BWire.z_start = szVec.at(i-1);
      BWire.y_end = eyVec.at(i-1);
      BWire.z_end = ezVec.at(i-1);
      BWire.isLowWire = false;

      BWires_U.push_back(BWire);
    }
    else if((i == 2399) && (CSstatusVec.at(i) == false) && (isGoodChannel == false)) {
      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = false;

      BWires_U.push_back(BWire);

    }
  }

  isGoodChannel = true;
  for(int i = 2400; i < 4800; i++) {
    if((CSstatusVec.at(i) == false) && (isGoodChannel == true)) {
      isGoodChannel = false;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = true;

      BWires_V.push_back(BWire);
    }
    else if((CSstatusVec.at(i) == true) && (isGoodChannel == false)) {
      isGoodChannel = true;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i-1);
      BWire.y_start = syVec.at(i-1);
      BWire.z_start = szVec.at(i-1);
      BWire.y_end = eyVec.at(i-1);
      BWire.z_end = ezVec.at(i-1);
      BWire.isLowWire = false;

      BWires_V.push_back(BWire);
    }
    else if((i == 4799) && (CSstatusVec.at(i) == false) && (isGoodChannel == false)) {
      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = false;

      BWires_V.push_back(BWire);
    }
  }

  isGoodChannel = true;
  for(int i = 4800; i < 8256; i++) {
    if((CSstatusVec.at(i) == false) && (isGoodChannel == true)) {
      isGoodChannel = false;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = true;

      BWires_Y.push_back(BWire);
    }
    else if((CSstatusVec.at(i) == true) && (isGoodChannel == false)) {
      isGoodChannel = true;

      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i-1);
      BWire.y_start = syVec.at(i-1);
      BWire.z_start = szVec.at(i-1);
      BWire.y_end = eyVec.at(i-1);
      BWire.z_end = ezVec.at(i-1);
      BWire.isLowWire = false;

      BWires_Y.push_back(BWire);
    }
    else if((i == 8255) && (CSstatusVec.at(i) == false) && (isGoodChannel == false)) {
      BoundaryWire BWire;
      BWire.wire_num = wireVec.at(i);
      BWire.y_start = syVec.at(i);
      BWire.z_start = szVec.at(i);
      BWire.y_end = eyVec.at(i);
      BWire.z_end = ezVec.at(i);
      BWire.isLowWire = false;

      BWires_Y.push_back(BWire);
    }
  }

  return;
}


//___________________________________________________________________________________________________
bool FindDeadRegions::NearDeadReg2P(float yVal, float zVal, float tolerance) {

  float minDist_U = 100000.0;
  float minDist_V = 100000.0;
  float minDist_Y = 100000.0;

  for (unsigned int i = 0; i < BWires_U.size(); i++) {
    float m = (BWires_U[i].y_end-BWires_U[i].y_start)/(BWires_U[i].z_end-BWires_U[i].z_start);
    float b = BWires_U[i].y_start - m*BWires_U[i].z_start;
    float dist = fabs(yVal-m*zVal-b)/sqrt(pow(m,2)+1.0);

    if (dist < minDist_U) {
      minDist_U = dist;
    }

    if (BWires_U[i].isLowWire == true) {
      float m_next = (BWires_U[i+1].y_end-BWires_U[i+1].y_start)/(BWires_U[i+1].z_end-BWires_U[i+1].z_start);
      float b_next = BWires_U[i+1].y_start - m*BWires_U[i+1].z_start;
      float y1 = m*zVal + b;
      float y2 = m_next*zVal + b_next;

      if ((yVal <= y1) && (yVal >= y2)) {
        minDist_U = 0.0;
      }
    }
  }

  for (unsigned int i = 0; i < BWires_V.size(); i++) {
    float m = (BWires_V[i].y_end-BWires_V[i].y_start)/(BWires_V[i].z_end-BWires_V[i].z_start);
    float b = BWires_V[i].y_start - m*BWires_V[i].z_start;
    float dist = fabs(yVal-m*zVal-b)/sqrt(pow(m,2)+1.0);

    if (dist < minDist_V) {
      minDist_V = dist;
    }
 
    if (BWires_V[i].isLowWire == true) {
      float m_next = (BWires_V[i+1].y_end-BWires_V[i+1].y_start)/(BWires_V[i+1].z_end-BWires_V[i+1].z_start);
      float b_next = BWires_V[i+1].y_start - m*BWires_V[i+1].z_start;
      float y1 = m*zVal + b;
      float y2 = m_next*zVal + b_next;

      if ((yVal >= y1) && (yVal <= y2)) {
        minDist_V = 0.0;
      }
    }
  }

  for (unsigned int i = 0; i < BWires_Y.size(); i++) {
    float dist = fabs(zVal-BWires_Y[i].z_start);

    if (dist < minDist_Y) {
      minDist_Y = dist;
    }

    if (BWires_Y[i].isLowWire == true) {
      float z1 = BWires_Y[i].z_start;
      float z2 = BWires_Y[i+1].z_start;

      if ((zVal >= z1) && (zVal <= z2)) {
        minDist_Y = 0.0;
      }
    }
  }

  if ((minDist_U < tolerance) && (minDist_V < tolerance)) {
    return true;
  }
  else if ((minDist_U < tolerance) && (minDist_Y < tolerance)) {
    return true;
  }
  else if ((minDist_V < tolerance) && (minDist_Y < tolerance)) {
    return true;
  }
  else {
    return false;
  }
}




//____________________________________________________________________________________________________
bool FindDeadRegions::NearDeadReg3P(float yVal, float zVal, float tolerance) {

  float minDist = 100000.0;

  for (unsigned int i = 0; i < BWires_U.size(); i++) {
    float m = (BWires_U[i].y_end-BWires_U[i].y_start)/(BWires_U[i].z_end-BWires_U[i].z_start);
    float b = BWires_U[i].y_start - m*BWires_U[i].z_start;
    float dist = fabs(yVal-m*zVal-b)/sqrt(pow(m,2)+1.0);

    if (dist < minDist) {
      minDist = dist;
    }

    if (BWires_U[i].isLowWire == true) {
      float m_next = (BWires_U[i+1].y_end-BWires_U[i+1].y_start)/(BWires_U[i+1].z_end-BWires_U[i+1].z_start);
      float b_next = BWires_U[i+1].y_start - m*BWires_U[i+1].z_start;
      float y1 = m*zVal + b;
      float y2 = m_next*zVal + b_next;

      if ((yVal <= y1) && (yVal >= y2)) {
        minDist = 0.0;
      }
    }
  }

  for (unsigned int i = 0; i < BWires_V.size(); i++) {
    float m = (BWires_V[i].y_end-BWires_V[i].y_start)/(BWires_V[i].z_end-BWires_V[i].z_start);
    float b = BWires_V[i].y_start - m*BWires_V[i].z_start;
    float dist = fabs(yVal-m*zVal-b)/sqrt(pow(m,2)+1.0);

    if (dist < minDist) {
      minDist = dist;
    }

    if (BWires_V[i].isLowWire == true) {
      float m_next = (BWires_V[i+1].y_end-BWires_V[i+1].y_start)/(BWires_V[i+1].z_end-BWires_V[i+1].z_start);
      float b_next = BWires_V[i+1].y_start - m*BWires_V[i+1].z_start;
      float y1 = m*zVal + b;
      float y2 = m_next*zVal + b_next;

      if ((yVal >= y1) && (yVal <= y2)) {
        minDist = 0.0;
      }
    }
  }

  for (unsigned int i = 0; i < BWires_Y.size(); i++) {
    float dist = fabs(zVal-BWires_Y[i].z_start);

    if (dist < minDist) {
      minDist = dist;
    }

    if (BWires_Y[i].isLowWire == true) {
      float z1 = BWires_Y[i].z_start;
      float z2 = BWires_Y[i+1].z_start;

      if ((zVal >= z1) && (zVal <= z2)) {
        minDist = 0.0;
      }
    }
  }

  if (minDist < tolerance) {
    return true;
  }
  else {
    return false;
  }
}



//__________________________________________________________________________________________________
TH2F* FindDeadRegions::GetDeadRegionHisto2P(){

  TH2F *deadReg2P = new TH2F("deadReg2P","",10350,0.0,1035.0,2300,-115.0,115.0);

  float tolVal = 0.5;

  for (int i = 1; i <= deadReg2P->GetNbinsX(); i++) {
    for (int j = 1; j <= deadReg2P->GetNbinsY(); j++) {
      if (NearDeadReg2P(deadReg2P->GetYaxis()->GetBinCenter(j),deadReg2P->GetXaxis()->GetBinCenter(i),tolVal) == true) {
        deadReg2P->SetBinContent(i,j,1.0);
      }
    }
  }

  return deadReg2P;
}



//__________________________________________________________________________________________________
TH2F* FindDeadRegions::GetDeadRegionHisto3P(){

  TH2F *deadReg3P = new TH2F("deadReg3P","",10350,0.0,1035.0,2300,-115.0,115.0);

  float tolVal = 0.5;

  for (int i = 1; i <= deadReg3P->GetNbinsX(); i++) {
    for (int j = 1; j <= deadReg3P->GetNbinsY(); j++) {
      if (NearDeadReg3P(deadReg3P->GetYaxis()->GetBinCenter(j),deadReg3P->GetXaxis()->GetBinCenter(i),tolVal) == true) {
        deadReg3P->SetBinContent(i,j,1.0);
      }
    }
  }

  return deadReg3P;
}





DECLARE_ART_SERVICE(FindDeadRegions, LEGACY)
DEFINE_ART_SERVICE(FindDeadRegions)
