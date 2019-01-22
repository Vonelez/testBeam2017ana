#ifndef REALDATATESTS_RESOLUTION_H
#define REALDATATESTS_RESOLUTION_H

#include "TGraphErrors.h"
#include <iostream>

class Resolution {
public:
  Resolution(TGraphErrors *resolGraph);

  virtual ~Resolution();

private:
  void init();

  TGraphErrors *resolGraph;

  void weightedMean();
};


#endif //REALDATATESTS_RESOLUTION_H
