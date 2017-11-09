/**
 * \class ubana::SelectionResult
 *
 * \ingroup UBXSec
 *
 * \brief Data product to store the event selection results
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
#include <map>

namespace ubana {

  class SelectionResult {

  public:

    SelectionResult();
    virtual ~SelectionResult();

    void SetSelectionStatus(bool status)                    {_selection_status = status;}
    void SetSelectionType(std::string type)                 {_selection_type = type;}
    void SetFailureReason(std::string reason)               {_failure_reason = reason;}
    void SetCutFlowStatus(std::map<std::string,bool> input) {_cut_flow_status = input;}

    /// Returns the type of event selection run (numu cc inclusive / ccpi0 / others...)
    const std::string & GetSelectionType() const;

    /// Returns the status of the selection (true is event passed, flase otherwise)
    const bool & GetSelectionStatus() const;

    /// If GetSelectionStatus() returns false, it returns a string containing the reason for the failure
    const std::string & GetFailureReason() const;

    /// Returns the effect of each cut (all true if event is selected)
    const std::map<std::string,bool> & GetCutFlowStatus() const;

  private:

    std::string                _selection_type;    ///< Event selection type (ccnumu, nue, ...)
    std::string                _failure_reason;    ///< Reason for selection failure
    bool                       _selection_status;  ///< True if event selected
    std::map<std::string,bool> _cut_flow_status;   ///< Stores the effect of each cut

 };
}

#endif /* SelectionResult_h */
