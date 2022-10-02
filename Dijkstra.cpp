#include <iostream>
#include "Dijkstra.h"
#include "find.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <vector>

// Function Definitions
double Dijkstra(const coder::array<float, 2U> &AdjacMatrix,
                const coder::array<float, 2U> &CostMatrix, double SID, double FID, std::vector<double> &pathIndexes) {

    coder::array<double, 1U> iTable;
    coder::array<double, 1U> jTable;
    coder::array<double, 1U> minCost;
    coder::array<int, 2U> E;
    coder::array<int, 1U> j;
    coder::array<int, 1U> nodeIndex;
    coder::array<boolean_T, 1U> auxArray1;
    coder::array<boolean_T, 1U> auxArray2;
    coder::array<boolean_T, 1U> isSettled;
    double b_I;
    int b_J;
    int idx;
    int nx;

    //  Process inputs
    nx = AdjacMatrix.size(0) * AdjacMatrix.size(1);
    nodeIndex.set_size(nx);
    j.set_size(nx);
    int counter = 0;
    int iIndex = 1;
    int jIndex = 1;
    while (jIndex <= AdjacMatrix.size(1)) {
        if (AdjacMatrix[(iIndex + AdjacMatrix.size(0) * (jIndex - 1)) - 1] != 0.0F) {
            nodeIndex[counter] = iIndex;
            j[counter] = jIndex;
            counter++;
        }
        iIndex++;
        if (iIndex > AdjacMatrix.size(0)) {
            iIndex = 1;
            jIndex++;
        }
    }

    nodeIndex.set_size(counter);
    j.set_size(counter);

    E.set_size(counter, 2);
    auxArray1.set_size(counter);
    for (int index = 0; index < counter; index++) {
        E[index] = nodeIndex[index];
        E[index + counter] = j[index];
    }

    //  Find the minimum costs and paths using Dijkstra's Algorithm
    //  Initializations
    int rows = AdjacMatrix.size(0);
    iTable.set_size(rows);
    jTable.set_size(rows);
    minCost.set_size(rows);
    isSettled.set_size(rows);
    auxArray2.set_size(rows);
    for (int index = 0; index < rows; index++) {
        iTable[index] = rtNaN;
        minCost[index] = rtInf;
        isSettled[index] = false;
    }

    b_I = SID;
    minCost[static_cast<int>(SID) - 1] = 0.0;
    iTable[static_cast<int>(SID) - 1] = 0.0;
    isSettled[static_cast<int>(SID) - 1] = true;

    std::vector<std::vector<double>> paths(rows, std::vector<double>(0));
    paths[static_cast<int>(SID) - 1].push_back(SID);

    //  Execute Dijkstra's Algorithm for this vertex
    while (!isSettled[static_cast<int>(FID) - 1]) {
        //  Update the table
        for (int index = 0; index < rows; index++) {
            jTable[index] = iTable[index];
        }

        iTable[static_cast<int>(b_I) - 1] = rtNaN;
        for (int index = 0; index < counter; index++) {
            auxArray1[index] = (E[index] == b_I);
        }

        coder::eml_find(auxArray1, nodeIndex);

        //  Calculate the costs to the neighbor nodes and record paths
        nx = nodeIndex.size(0);
        for (int index = 0; index < nx; index++) {
            b_J = E[(nodeIndex[index] + counter) - 1] - 1;
            if (!isSettled[b_J]) {
                float c = CostMatrix[(static_cast<int>(b_I) + rows * b_J) - 1];
                if (rtIsNaN(jTable[b_J]) ||
                    (jTable[b_J] > static_cast<float>(jTable[static_cast<int>(b_I) - 1]) + c)) {
                    iTable[b_J] = static_cast<float>(jTable[static_cast<int>(b_I) - 1]) + c;
                    paths[static_cast<int>(b_J) - 1] = (paths[static_cast<int>(b_I) - 1]);
                    paths[static_cast<int>(b_J) - 1].push_back(b_J);
                } else {
                    iTable[b_J] = jTable[b_J];
                }
            }
        }

        //  Find values in the table
        for (int index = 0; index < rows; index++) {
            auxArray2[index] = !rtIsNaN(iTable[index]);
        }

        coder::eml_find(auxArray2, nodeIndex);

        //  Settle the minimum value in the table
        nx = nodeIndex.size(0);
        idx = 1;

        b_I = iTable[static_cast<int>(nodeIndex[idx - 1]) - 1];
        int aux = idx - 1;
        idx++;
        for (int index = idx; index <= nx; index++) {
            double d = iTable[static_cast<int>(nodeIndex[index - 1]) - 1];
            if (b_I > d) {
                b_I = d;
                aux = index - 1;
            }
        }

        b_I = nodeIndex[aux];
        nx = static_cast<int>(b_I) - 1;
        minCost[nx] = iTable[nx];
        isSettled[nx] = true;
    }

    //  Store costs and paths
    pathIndexes = paths[static_cast<int>(FID) - 1];
    return minCost[static_cast<int>(FID) - 1];
}
