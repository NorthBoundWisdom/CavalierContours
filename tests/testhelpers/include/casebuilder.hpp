#ifndef CAVC_CASEBUILDER_HPP
#define CAVC_CASEBUILDER_HPP

#include "cavaliercontours.h"
#include <vector>

class CaseBuilder {
  static std::vector<cavc_vertex> simpleRectangle();

  static std::vector<cavc_vertex> twoBadArcCase();

  static std::vector<cavc_vertex> figureEightCase();

  static std::vector<cavc_vertex> complexSelfInterectCase();

  static std::vector<cavc_vertex> closedLineArcCase();
  static std::vector<cavc_vertex> offsetCase();
  static std::vector<cavc_vertex> quarterArcCase();
  static std::vector<cavc_vertex> positiveCircle();
  static std::vector<cavc_vertex> negativeCircle();

  static std::vector<std::vector<cavc_vertex>> simpleBoolCase();
};
#endif // CAVC_CASEBUILDER_HPP
