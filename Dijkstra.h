#ifndef DIJKSTRA_H
#define DIJKSTRA_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern double Dijkstra(const coder::array<float, 2U> &AorV, const coder::array<
        float, 2U> &xyCorE, double SID, double FID, std::vector<double> &pathIndexes);

#endif