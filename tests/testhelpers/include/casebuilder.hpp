#ifndef CAVC_CASEBUILDER_HPP
#define CAVC_CASEBUILDER_HPP

#include <cavc/plinesegment.hpp>
#include <vector>

using PlineVertex = cavc::PlineVertex<double>;
class CaseBuilder {
public:
  static std::vector<PlineVertex> simpleRectangle();

  static std::vector<PlineVertex> twoBadArcCase();

  static std::vector<PlineVertex> figureEightCase();

  static std::vector<PlineVertex> complexSelfInterectCase();

  static std::vector<PlineVertex> closedLineArcCase();

  static std::vector<PlineVertex> offsetCase();

  static std::vector<PlineVertex> quarterArcCase();

  static std::vector<PlineVertex> positiveCircle();

  static std::vector<PlineVertex> negativeCircle();

  static std::vector<std::vector<PlineVertex>> simpleBoolCase();
};
#endif // CAVC_CASEBUILDER_HPP
