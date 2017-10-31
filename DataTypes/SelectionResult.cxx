#include "SelectionResult.h"

namespace ubana {

  SelectionResult::SelectionResult() {
  }

  SelectionResult::~SelectionResult(){
  }

  const std::string & SelectionResult::GetSelectionType()   const {return _selection_type;   }

  const bool        & SelectionResult::GetSelectionStatus() const {return _selection_status; }

  const std::string & SelectionResult::GetFailureReason()   const {return _failure_reason;   }
}


