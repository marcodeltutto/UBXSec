/**
 * \file FiducialVolume.h
 *
 * \ingroup UBXSec
 * 
 * \brief Class def header for a class FiducialVolume
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

    @{*/
#ifndef FIDUCIALVOLUME_H
#define FIDUCIALVOLUME_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoVector.h"

namespace ubana{
  
  /**
   \class FiducialVolume
   User defined class FiducialVolume ... these comments are used to generate
   doxygen documentation!
 */

  class FiducialVolume {
    
  public:
    
    /// Default constructor
    FiducialVolume();

    /// Default destructor
    ~FiducialVolume(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p, double det_half_height, double det_width, double det_length);

    /// Printd the current configuration
    void PrintConfig();

    /// Returns true if the point is in the FV
    bool InFV(double x, double y, double z);

    /// Returns true if the point is in the FV 
    bool InFV(double* x);

  protected:

    double _det_half_height;
    double _det_width;
    double _det_length;

    double _border_x_low;
    double _border_x_high;
    double _border_y_low;
    double _border_y_high;
    double _border_z_low;        
    double _border_z_high;

    bool _configured = false;
  };
}

#endif
/** @} */ // end of doxygen group 

