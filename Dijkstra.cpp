//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  Dijkstra.cpp
//
//  Code generation for function 'Dijkstra'
//


// Include files
#include "Dijkstra.h"
#include "find.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "rt_nonfinite.h"

// Function Definitions
double Dijkstra(const coder::array<float, 2U> &AorV, const coder::array<float,
                2U> &xyCorE, double SID, double FID)
{
  coder::array<double, 1U> iTable;
  coder::array<double, 1U> jTable;
  coder::array<double, 1U> minCost;
  coder::array<float, 2U> A;
  coder::array<float, 1U> d;
  coder::array<int, 2U> E;
  coder::array<int, 1U> j;
  coder::array<int, 1U> nodeIndex;
  coder::array<boolean_T, 1U> b_E;
  coder::array<boolean_T, 1U> isSettled;
  double b_I;
  int idx;
  int jj;
  int nv;
  int nx;
  boolean_T exitg1;

  //  Process inputs
  //  -------------------------------------------------------------------
  if ((AorV.size(0) == 1) && (AorV.size(1) == 1)) {
    d.set_size(1);
    d[0] = AorV[0];
  } else {
    nx = AorV.size(0);
    nv = AorV.size(1);
    if (nx < nv) {
      nv = nx;
    }

    if (0 < AorV.size(1)) {
      nx = nv;
    } else {
      nx = 0;
    }

    d.set_size(nx);
    idx = nx - 1;
    for (jj = 0; jj <= idx; jj++) {
      d[jj] = AorV[jj + AorV.size(0) * jj];
    }
  }

  nv = d.size(0);
  A.set_size(d.size(0), d.size(0));
  nx = d.size(0) * d.size(0);
  for (idx = 0; idx < nx; idx++) {
    A[idx] = 0.0F;
  }

  for (nx = 0; nx < nv; nx++) {
    A[nx + A.size(0) * nx] = d[nx];
  }

  nx = AorV.size(0) * AorV.size(1);
  A.set_size(AorV.size(0), AorV.size(1));
  for (idx = 0; idx < nx; idx++) {
    A[idx] = AorV[idx] - A[idx];
  }

  //  Inputs = (A,cost)
  //  Convert adjacency matrix to edge list
  nx = A.size(0) * A.size(1);
  if (nx == 0) {
    nodeIndex.set_size(0);
    j.set_size(0);
  } else {
    idx = 0;
    nodeIndex.set_size(nx);
    j.set_size(nx);
    nv = 1;
    jj = 1;
    exitg1 = false;
    while ((!exitg1) && (jj <= A.size(1))) {
      boolean_T guard1 = false;
      guard1 = false;
      if (A[(nv + A.size(0) * (jj - 1)) - 1] != 0.0F) {
        idx++;
        nodeIndex[idx - 1] = nv;
        j[idx - 1] = jj;
        if (idx >= nx) {
          exitg1 = true;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        nv++;
        if (nv > A.size(0)) {
          nv = 1;
          jj++;
        }
      }
    }

    if (nx == 1) {
      if (idx == 0) {
        nodeIndex.set_size(0);
        j.set_size(0);
      }
    } else {
      if (1 > idx) {
        nx = 0;
      } else {
        nx = idx;
      }

      nodeIndex.set_size(nx);
      if (1 > idx) {
        idx = 0;
      }

      j.set_size(idx);
    }
  }

  E.set_size(nodeIndex.size(0), 2);
  nx = nodeIndex.size(0);
  for (idx = 0; idx < nx; idx++) {
    E[idx] = nodeIndex[idx];
  }

  nx = j.size(0);
  for (idx = 0; idx < nx; idx++) {
    E[idx + E.size(0)] = j[idx];
  }

  //  Find the minimum costs and paths using Dijkstra's Algorithm
  //  Initializations
  iTable.set_size(AorV.size(0));
  nx = AorV.size(0);
  for (idx = 0; idx < nx; idx++) {
    iTable[idx] = rtNaN;
  }

  minCost.set_size(AorV.size(0));
  nx = AorV.size(0);
  for (idx = 0; idx < nx; idx++) {
    minCost[idx] = rtInf;
  }

  isSettled.set_size(AorV.size(0));
  nx = AorV.size(0);
  for (idx = 0; idx < nx; idx++) {
    isSettled[idx] = false;
  }

  b_I = SID;
  minCost[static_cast<int>(SID) - 1] = 0.0;
  iTable[static_cast<int>(SID) - 1] = 0.0;
  isSettled[static_cast<int>(SID) - 1] = true;

  //  Execute Dijkstra's Algorithm for this vertex
  exitg1 = false;
  while (!(exitg1 || isSettled[static_cast<int>(FID) - 1])) {
    //  Update the table
    jTable.set_size(iTable.size(0));
    nx = iTable.size(0);
    for (idx = 0; idx < nx; idx++) {
      jTable[idx] = iTable[idx];
    }

    iTable[static_cast<int>(b_I) - 1] = rtNaN;
    nx = E.size(0);
    b_E.set_size(E.size(0));
    for (idx = 0; idx < nx; idx++) {
      b_E[idx] = (E[idx] == b_I);
    }

    coder::eml_find(b_E, nodeIndex);

    //  Calculate the costs to the neighbor nodes and record paths
    idx = nodeIndex.size(0);
    for (nv = 0; nv < idx; nv++) {
      nx = E[(nodeIndex[nv] + E.size(0)) - 1] - 1;
      if (!isSettled[nx]) {
        float c;
        c = xyCorE[(static_cast<int>(b_I) + xyCorE.size(0) * nx) - 1];
        if (rtIsNaN(jTable[nx]) || (jTable[nx] > static_cast<float>(jTable[
              static_cast<int>(b_I) - 1]) + c)) {
          iTable[nx] = static_cast<float>(jTable[static_cast<int>(b_I) - 1]) + c;
        } else {
          iTable[nx] = jTable[nx];
        }
      }
    }

    //  Find values in the table
    b_E.set_size(iTable.size(0));
    nx = iTable.size(0);
    for (idx = 0; idx < nx; idx++) {
      b_E[idx] = rtIsNaN(iTable[idx]);
    }

    nx = b_E.size(0);
    for (idx = 0; idx < nx; idx++) {
      b_E[idx] = !b_E[idx];
    }

    coder::eml_find(b_E, nodeIndex);
    jTable.set_size(nodeIndex.size(0));
    nx = nodeIndex.size(0);
    for (idx = 0; idx < nx; idx++) {
      jTable[idx] = nodeIndex[idx];
    }

    if (jTable.size(0) == 0) {
      exitg1 = true;
    } else {
      //  Settle the minimum value in the table
      nv = jTable.size(0);
      if (jTable.size(0) <= 2) {
        if (jTable.size(0) == 1) {
          nx = 0;
        } else {
          double b_d;
          b_d = iTable[static_cast<int>(jTable[1]) - 1];
          if ((iTable[static_cast<int>(jTable[0]) - 1] > b_d) || (rtIsNaN
               (iTable[static_cast<int>(jTable[0]) - 1]) && (!rtIsNaN(b_d)))) {
            nx = 1;
          } else {
            nx = 0;
          }
        }
      } else {
        if (!rtIsNaN(iTable[static_cast<int>(jTable[0]) - 1])) {
          idx = 1;
        } else {
          boolean_T exitg2;
          idx = 0;
          jj = 2;
          exitg2 = false;
          while ((!exitg2) && (jj <= jTable.size(0))) {
            if (!rtIsNaN(iTable[static_cast<int>(jTable[jj - 1]) - 1])) {
              idx = jj;
              exitg2 = true;
            } else {
              jj++;
            }
          }
        }

        if (idx == 0) {
          nx = 0;
        } else {
          b_I = iTable[static_cast<int>(jTable[idx - 1]) - 1];
          nx = idx - 1;
          idx++;
          for (jj = idx; jj <= nv; jj++) {
            double b_d;
            b_d = iTable[static_cast<int>(jTable[jj - 1]) - 1];
            if (b_I > b_d) {
              b_I = b_d;
              nx = jj - 1;
            }
          }
        }
      }

      b_I = jTable[nx];
      nx = static_cast<int>(jTable[nx]) - 1;
      minCost[nx] = iTable[nx];
      isSettled[nx] = true;
    }
  }

  //  Store costs and paths
  return minCost[static_cast<int>(FID) - 1];
}

// End of code generation (Dijkstra.cpp)
