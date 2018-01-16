#ifndef FIDUCIALVOLUME_CXX
#define FIDUCIALVOLUME_CXX

#include "FiducialVolume.h"
#include <iostream>

namespace ubana {

  FiducialVolume::FiducialVolume()
  {
  }

  void FiducialVolume::Configure(fhicl::ParameterSet const& pset, double det_half_height, double det_width, double det_length)
  {
    _border_x_low    = pset.get<double>("BorderXLow",  10);
    _border_x_high   = pset.get<double>("BorderXHigh", 10);
    _border_y_low    = pset.get<double>("BorderYLow",  20);
    _border_y_high   = pset.get<double>("BorderYHigh", 20);
    _border_z_low    = pset.get<double>("BorderZLow",  10);
    _border_z_high   = pset.get<double>("BorderZHigh", 10);

    _det_half_height = det_half_height;
    _det_width       = det_width;
    _det_length      = det_length;

    _configured = true;
  }

  void FiducialVolume::PrintConfig() {

    std::cout << "--- FiducialVolume configuration:" << std::endl;
    std::cout << "---   _border_x_low    = " << _border_x_low << std::endl;
    std::cout << "---   _border_x_high   = " << _border_x_high << std::endl;
    std::cout << "---   _border_y_low    = " << _border_y_low << std::endl;
    std::cout << "---   _border_y_high   = " << _border_y_high << std::endl;
    std::cout << "---   _border_z_low    = " << _border_z_low << std::endl;
    std::cout << "---   _border_z_high   = " << _border_z_high << std::endl;
    std::cout << "---   _det_half_height = " << _det_half_height << std::endl;
    std::cout << "---   _det_width       = " << _det_width << std::endl;
    std::cout << "---   _det_length      = " << _det_length << std::endl;

  }

  bool FiducialVolume::InFV(double* x) {

    return this->InFV(x[0], x[1], x[2]);

  }

  bool FiducialVolume::InFV(TVector3 x) {

    return this->InFV(x.X(), x.Y(), x.Z());

  }

  bool FiducialVolume::InFV(TVector3 x1, TVector3 x2) {

    return (this->InFV(x1.X(), x1.Y(), x1.Z()) && this->InFV(x2.X(), x2.Y(), x2.Z()));

  }

  bool FiducialVolume::InFV(double x, double y, double z) {

    // Construct a vector
    ::geoalgo::Vector the_point(x, y, z);

    // Construct the fiducial volume
    ::geoalgo::AABox fidvol(_border_x_low, 
                            -1.*_det_half_height + _border_y_low, 
                            _border_z_low,
                            _det_width - _border_x_high, 
                            _det_half_height - _border_y_high, 
                            _det_length - _border_z_high);

    // Check if the vector is in the FV
    if(fidvol.Contain(the_point)) {
      return true;
    } else {
      return false;
    }


  }

}


#endif
