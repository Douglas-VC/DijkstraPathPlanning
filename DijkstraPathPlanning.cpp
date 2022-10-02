#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/common/common.h>
#include <pcl/octree/octree_search.h>
#include <cmath>
#include <chrono>

using namespace std::chrono;

#include "coder_array.h"
#include "Dijkstra.h"

typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointCloud;
typedef pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> Octree;
typedef coder::array<float, 2> Matrix;
typedef std::vector<double> Vector;

void
printExecutionTime(std::chrono::high_resolution_clock::time_point start,
                   std::chrono::high_resolution_clock::time_point stop,
                   const std::string &mode) {

    if (mode == "setup") {
        std::chrono::duration<double, std::milli> duration = stop - start;
        std::cout << "Setup levou "
                  << duration.count() / 1000
                  << " segundos para executar."
                  << std::endl;
    } else if (mode == "dijkstra") {
        std::chrono::duration<double, std::milli> duration = stop - start;
        std::cout << "Dijkstra levou "
                  << duration.count() / 1000
                  << " segundos para executar."
                  << std::endl;
    }

}

PointCloud importPointCloud() {

    PointCloud cloud(new pcl::PointCloud<pcl::PointXYZ>);

    if (pcl::io::loadPCDFile<pcl::PointXYZ>("PointClouds/campinhoremaster3.pcd", *cloud) == -1) {
        PCL_ERROR("Couldn't read point cloud file\n");
        exit(-1);
    }

    std::cout << "Loaded "
              << cloud->width * cloud->height
              << " data points from point cloud. "
              << std::endl;

    return cloud;
}

Octree transformPointCloudToOctree(const PointCloud &cloud, float resolution) {
    Octree octree(resolution);

    octree.setInputCloud(cloud);
    octree.addPointsFromInputCloud();
    return octree;
}

void getMatrixDimensions(const PointCloud &cloud, int &rows, int &columns, float graphResolution) {
    pcl::PointXYZ minPt, maxPt;
    pcl::getMinMax3D(*cloud, minPt, maxPt);

    rows = floor(fabs((double) (maxPt.x - minPt.x)) / graphResolution) + 1;
    columns = floor(fabs((double) (maxPt.y - minPt.y)) / graphResolution) + 1;
}

void populateAdjacencyMatrix(Matrix &AdjacencyMatrix, int rows, int columns) {
    int N = rows * columns;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j > i + 2 * columns) {
                break;
            } else if ((j == i + 1 && (i % (columns - 1) != 0 || i == 0))
                       || (j == i + columns)
                       || (j == i + columns - 1 && j % (columns - 1) != 0)
                       || (j == i + columns + 1 && (i % (columns - 1) != 0 || i == 0))) {
                AdjacencyMatrix.at(i, j) = 1;
                AdjacencyMatrix.at(j, i) = 1;
            }
        }
    }
}

std::vector<double> makeStepVector(double beg, double step, double end) {
    std::vector<double> vec;
    vec.reserve(fabs((end - beg)) / step + 1);
    while (beg <= end) {
        vec.push_back(beg);
        beg += step;
    }
    return vec;
}

void getGraphCoordinates(Octree octree, const PointCloud &cloud, Matrix &graphCoordinates, int rows,
                         int columns,
                         float graphResolution) {
    pcl::PointXYZ minPt, maxPt;
    std::vector<int> pointIdxNKNSearch;
    std::vector<float> pointNKNSquaredDistance;

    pcl::getMinMax3D(*cloud, minPt, maxPt);

    auto iValues = makeStepVector(minPt.x, graphResolution, maxPt.x);
    auto jValues = makeStepVector(minPt.y, graphResolution, maxPt.y);

    pcl::PointXYZ searchPoint;
    searchPoint.x = 0.0f;
    searchPoint.y = 0.0f;
    searchPoint.z = 0.0f;

    int index = 0;

    for (int j = columns - 1; j >= 0; j--) {
        for (int i = 0; i < rows; i++) {
            searchPoint.x = (float) iValues[i];
            searchPoint.y = (float) jValues[j];
            octree.nearestKSearch(searchPoint, 1, pointIdxNKNSearch, pointNKNSquaredDistance);
            graphCoordinates.at(index, 1) = (*cloud)[pointIdxNKNSearch[0]].x;
            graphCoordinates.at(index, 2) = (*cloud)[pointIdxNKNSearch[0]].y;
            graphCoordinates.at(index, 3) = (*cloud)[pointIdxNKNSearch[0]].z;
            index += 1;
        }
    }
}

float dist3D(float p1_x, float p1_y, float p1_z, float p2_x, float p2_y, float p2_z) {
    return sqrt(pow(p1_x - p2_x, 2) + pow(p1_y - p2_y, 2) + pow(p1_z - p2_z, 2) * 1.0);
}

void calculateCostMatrix(Matrix &AdjacencyMatrix, Matrix &CostMatrix, Matrix &graphCoordinates, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (AdjacencyMatrix.at(i, j) == 1) {
                CostMatrix.at(i, j) = dist3D(
                        graphCoordinates.at(i, 1), graphCoordinates.at(i, 2), graphCoordinates.at(i, 3),
                        graphCoordinates.at(j, 1), graphCoordinates.at(j, 2), graphCoordinates.at(j, 3));
            }
        }
    }
}

int getCorrespondentGraphNodes(Matrix &graphCoordinates, int N, float pointX, float pointY) {
    int index = 0;
    float distance = std::numeric_limits<float>::infinity();
    float temp_distance = 0;

    for (int i = 0; i < N; i++) {
        temp_distance = sqrt(pow(pointX - graphCoordinates.at(i, 1), 2) + pow(pointY - graphCoordinates.at(i, 2), 2));
        if (temp_distance < distance) {
            distance = temp_distance;
            index = i;
        }
    }

    return index;
}

int main() {
    auto start = high_resolution_clock::now();

    float octreeResolution = 0.1f; // Resolução da octree
    float graphResolution = 0.5f; // Resolução do grafo em metros

    float startPositionX = 13.8f;
    float startPositionY = -4.1f;
    float goalPositionX = -15.0f;
    float goalPositionY = 1.25f;

    int startNodeID = 0;
    int goalNodeID = 0;

    PointCloud cloud = importPointCloud();
    Octree octree = transformPointCloudToOctree(cloud, octreeResolution);

    Matrix AdjacencyMatrix;
    Matrix CostMatrix;
    Matrix graphCoordinates;
    Vector pathIndexes;

    int rows = 0;
    int columns = 0;
    getMatrixDimensions(cloud, rows, columns, graphResolution);
    int N = rows * columns;

    AdjacencyMatrix.set_size(N, N);
    CostMatrix.set_size(N, N);
    graphCoordinates.set_size(N * 3, 3);

    populateAdjacencyMatrix(AdjacencyMatrix, rows, columns);
    getGraphCoordinates(octree, cloud, graphCoordinates, rows, columns, graphResolution);
    calculateCostMatrix(AdjacencyMatrix, CostMatrix, graphCoordinates, rows * columns);
    startNodeID = getCorrespondentGraphNodes(graphCoordinates, N, startPositionX, startPositionY);
    goalNodeID = getCorrespondentGraphNodes(graphCoordinates, N, goalPositionX, goalPositionY);

    auto stop = high_resolution_clock::now();
    printExecutionTime(start, stop, "setup");

    /*------------Dijkstra------------*/

    start = high_resolution_clock::now();

    double totalCost = Dijkstra(AdjacencyMatrix, CostMatrix, startNodeID, goalNodeID, pathIndexes);

    stop = high_resolution_clock::now();
    printExecutionTime(start, stop, "dijkstra");

    std::cout << "Custo total: " << totalCost << std::endl;

    return 0;
}
