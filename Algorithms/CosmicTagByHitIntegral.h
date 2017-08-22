/**
 * \class CosmicTagByHitIntegral
 *
 * \ingroup UBXSec
 *
 * \brief CosmicTag
 * 
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/08/18 $
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Friday, August 18, 2017 at 08:53:50
 *
 */


#ifndef COSMICTAGBYHITINTEGRAL
#define COSMICTAGBYHITINTEGRAL

#include <iostream>
#include <fstream>

//#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h"
#include "uboone/UBXSec/Algorithms/CosmicTagToolInterface.h"


namespace ubana {

  class CosmicTagByHitIntegral : public CosmicTagToolInterface {
  
  public:

    /// Constructor
    explicit CosmicTagByHitIntegral(fhicl::ParameterSet const& ps) :
      _test {ps.get<double>("test")},
      _counter(0)
    {
      _csvfile.open ("cosmictagbyhitintegral.csv", std::ofstream::out | std::ofstream::trunc);
      _csvfile << "n,wire,integral" << std::endl;
    }

    /// Constructor
    explicit CosmicTagByHitIntegral() : 
      _test (0), _counter(0)
    {
      _csvfile.open ("cosmictagbyhitintegral.csv", std::ofstream::out | std::ofstream::trunc);
      _csvfile << "n,wire,integral" << std::endl;
    }

    /// Description
    bool IsCosmic(SimpleHitVector) override;


  private:

    double _test;
    int _counter; 
    std::ofstream _csvfile;
  };
}

//DEFINE_ART_CLASS_TOOL(ubana::CosmicTagByHitIntegral)

#endif //  COSMICTAGBYHITINTEGRAL
