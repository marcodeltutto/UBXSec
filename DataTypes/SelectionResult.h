/**
 * \class ubana::SelectionResult
 *
 * \ingroup UBXSec
 *
 * \brief Data product to store a flash matching results
 * 
 *
 * \author $Author: Marco Del Tutto<marco.deltutto@physics.ox.ac.uk> $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2017/10/09 $
 *
 * Contact: marco.deltutto@physics.ox.ac.uk
 *
 * Created on: Monday, October 09, 2017 at 10:21:43
 *
 */

#ifndef SelectionResult_h
#define SelectionResult_h

#include <vector>
#include <string>

namespace ubana {

  class SelectionResult {

  public:

    SelectionResult();
    virtual ~SelectionResult();

    void SetSelectionStatus(bool status)    {_selection_status = status;}
    void SetSelectionType(std::string type) {_selection_type = type;}

    std::string GetSelectionType()   {return _selection_type;}
    bool        GetSelectionStatus() {return _selection_status;}

  private:

    std::string _selection_type;    ///< Event selection type (ccnumu, nue, ...)
    bool        _selection_status;  ///< True if event selected

 };
}

#endif /* SelectionResult_h */
