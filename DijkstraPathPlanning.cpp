#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/common/common.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/octree/octree_search.h>
#include <cmath>
#include <chrono>

using namespace std::chrono;

#include "coder_array.h"
#include "Dijkstra.h"

typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointCloud;
typedef pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> Octree;

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

void getMatrixDimensions(const PointCloud &cloud, int &rows, int &columns, float resolution) {
    pcl::PointXYZ minPt, maxPt;
    pcl::getMinMax3D(*cloud, minPt, maxPt);

    rows = floor(fabs((double) (maxPt.x - minPt.x)) / resolution);
    columns = floor(fabs((double) (maxPt.y - minPt.y)) / resolution);
}

void populateAdjacencyMatrix(coder::array<float, 2> &AdjacencyMatrix, int rows, int columns) {
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

void getKNearestNeighboursFromOctree(Octree octree, const PointCloud &cloud, int K, pcl::PointXYZ searchPoint) {

    std::vector<int> pointIdxNKNSearch;
    std::vector<float> pointNKNSquaredDistance;

    std::cout << "K nearest neighbor search at (" << searchPoint.x
              << " " << searchPoint.y
              << " " << searchPoint.z
              << ") with K=" << K << std::endl;

    if (octree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0) {
        for (std::size_t i = 0; i < pointIdxNKNSearch.size(); ++i)
            std::cout << "    " << (*cloud)[pointIdxNKNSearch[i]].x
                      << " " << (*cloud)[pointIdxNKNSearch[i]].y
                      << " " << (*cloud)[pointIdxNKNSearch[i]].z
                      << " (squared distance: " << pointNKNSquaredDistance[i] << ")" << std::endl;
    }
}

int main() {
    auto start = high_resolution_clock::now();

    float octreeResolution = 0.1f; // Resolução da octree
    float graphResolution = 0.6f; // Resolução do grafo em metros

    PointCloud cloud = importPointCloud();
    Octree octree = transformPointCloudToOctree(cloud, octreeResolution);

    coder::array<float, 2> AdjacencyMatrix;
    coder::array<float, 2> CostMatrix;

    int rows = 0;
    int columns = 0;
    getMatrixDimensions(cloud, rows, columns, graphResolution);
    int N = rows * columns;

    AdjacencyMatrix.set_size(N, N);
    CostMatrix.set_size(N, N);

    populateAdjacencyMatrix(AdjacencyMatrix, rows, columns);

    double resultado = Dijkstra(AdjacencyMatrix, CostMatrix, 1200, 1500);

    std::cout << "Custo total: " << resultado << std::endl;


//
//
//
//
//
//
//
//
//
//

    auto stop = high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = stop - start;
    std::cout << "Código levou "
              << duration.count() / 1000
              << " segundos para executar."
              << std::endl;

//    int K = 10;
//    pcl::PointXYZ searchPoint;
//    searchPoint.x = 6.25f;
//    searchPoint.y = -14.15f;
//    searchPoint.z = 0.35f;
//    getKNearestNeighboursFromOctree(octree, cloud, K, searchPoint);

    return 0;
}
