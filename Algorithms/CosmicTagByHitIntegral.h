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
//#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h"
#include "uboone/UBXSec/Algorithms/CosmicTagToolInterface.h"


namespace ubana {

  class CosmicTagByHitIntegral : public CosmicTagToolInterface {
  
  public:

    /// Constructor
    explicit CosmicTagByHitIntegral(fhicl::ParameterSet const& ps) :
      _test {ps.get<double>("test")}
    {}

    /// Description
    bool IsCosmic(SimpleHitVector) override;


  private:
    double _test;
    
  };
}

//DEFINE_ART_CLASS_TOOL(ubana::CosmicTagByHitIntegral)

#endif //  COSMICTAGBYHITINTEGRAL
