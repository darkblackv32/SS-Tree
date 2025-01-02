#include <iostream>
#include <vector>
#include <unordered_set>
#include <random>
#include "Point.h"
#include "Data.h"
#include "SSTree.h"
#include <chrono> 

constexpr size_t NUM_POINTS = 10000;
constexpr size_t MAX_POINTS_PER_NODE = 20;

/*
 * Helper functions
 */

std::vector<Data*> generateRandomData(size_t numPoints) {
    std::vector<Data*> data;
    for (size_t i = 0; i < numPoints; ++i) {
        Point embedding = Point::random();
        std::string imagePath = "eda_" + std::to_string(i) + ".jpg";
        Data* dataPoint = new Data(embedding, imagePath);
        data.push_back(dataPoint);
    }
    return data;
}

void collectDataDFS(SSNode* node, std::unordered_set<Data*>& treeData) {
    if (node->getIsLeaf()) {
        for (const auto& d : node->getData()) {
            treeData.insert(d);
        }
    } else {
        for (const auto& child : node->getChildren()) {
            collectDataDFS(child, treeData);
        }
    }
}

/*
 * Testing functions
 */

// Test 1: Check if all data is present in the tree
bool allDataPresent(const SSTree& tree, const std::vector<Data*>& data) {
    std::unordered_set<Data*> dataSet(data.begin(), data.end());
    std::unordered_set<Data*> treeData;

    collectDataDFS(tree.getRoot(), treeData);
    for (const auto& d : dataSet) {
        if (treeData.find(d) == treeData.end()) {
            return false;
        }
    }
    for (const auto& d : treeData) {
        if (dataSet.find(d) == dataSet.end()) {
            return false;
        }
    }
    return true;
}

// Test 2: Check if all leaves are at the same level
bool leavesAtSameLevelDFS(SSNode* node, int level, int& leafLevel) {
    if (node->getIsLeaf()) {
        if (leafLevel == -1) leafLevel = level;
        return leafLevel == level;
    }
    for (const auto& child : node->getChildren()) {
        if (!leavesAtSameLevelDFS(child, level + 1, leafLevel)) return false;
    }
    return true;
}

bool leavesAtSameLevel(SSNode* root) {
    int leafLevel = -1;
    return leavesAtSameLevelDFS(root, 0, leafLevel);
}

// Test 3: Check if no node exceeds the maximum number of children
bool noNodeExceedsMaxChildrenDFS(SSNode* node, size_t maxPointsPerNode) {
    if (node->getChildren().size() > maxPointsPerNode) return false;
    for (const auto& child : node->getChildren()) {
        if (!noNodeExceedsMaxChildrenDFS(child, maxPointsPerNode)) return false;
    }
    return true;
}

bool noNodeExceedsMaxChildren(SSNode* root, size_t maxPointsPerNode) {
    return noNodeExceedsMaxChildrenDFS(root, maxPointsPerNode);
}

// Test 4: Check if all points are inside the bounding sphere of their respective nodes
bool sphereCoversAllPointsDFS(SSNode* node) {
    if (!node->getIsLeaf()) return true;
    const Point& centroid = node->getCentroid();
    float radius = node->getRadius();
    for (const auto& data : node->getData()) {
        if (Point::distance(centroid, data->getEmbedding()) > radius) return false;
    }
    return true;
}

bool dfsSphereCoversAllPoints(SSNode* node) {
    if (node->getIsLeaf()) {
        return sphereCoversAllPointsDFS(node);
    } else {
        for (const auto& child : node->getChildren()) {
            if (!dfsSphereCoversAllPoints(child)) return false;
        }
    }
    return true;
}

bool sphereCoversAllPoints(SSNode* root) {
    return dfsSphereCoversAllPoints(root);
}

// Test 5: Check if all children are inside the bounding sphere of their parent node
bool sphereCoversAllChildrenSpheresDFS(SSNode* node) {
    if (node->getIsLeaf()) return true;
    const Point& centroid = node->getCentroid();
    float radius = node->getRadius();
    for (const auto& child : node->getChildren()) {
        const Point& childCentroid = child->getCentroid();
        float childRadius = child->getRadius();
        if (Point::distance(centroid, childCentroid) + childRadius > radius) return false;
    }
    return true;
}

bool dfsSphereCoversAllChildrenSpheres(SSNode* node) {
    if (!sphereCoversAllChildrenSpheresDFS(node)) return false;
    for (const auto& child : node->getChildren()) {
        if (!dfsSphereCoversAllChildrenSpheres(child)) return false;
    }
    return true;
}

bool sphereCoversAllChildrenSpheres(SSNode* root) {
    return dfsSphereCoversAllChildrenSpheres(root);
}

// Test 6: Verify KNN search consistency by comparing tree results with manually sorted neighbors.
bool correctKnnSearch(const SSTree &tree, std::vector<Data *> &data) {
    Point query = Point::random();
    size_t k = 1;
    auto resultUsingTree = tree.knn(query, k);
    std::sort(data.begin(), data.end(), [&query](Data *a, Data *b) {
        return a->getEmbedding().distance(query) < b->getEmbedding().distance(query);
    });
    data.resize(k);
    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i] != resultUsingTree[i]) {
            return false;
        }
    }
    return true;
}


int main() {

    auto start = std::chrono::high_resolution_clock::now();

    auto data = generateRandomData(NUM_POINTS);
    SSTree tree(MAX_POINTS_PER_NODE);
    for (const auto &d: data) {
        tree.insert(d);
    }

    bool allPresent = allDataPresent(tree, data);
    bool sameLevel = leavesAtSameLevel(tree.getRoot());
    bool noExceed = noNodeExceedsMaxChildren(tree.getRoot(), MAX_POINTS_PER_NODE);
    bool spherePoints = sphereCoversAllPoints(tree.getRoot());
    bool sphereChildren = sphereCoversAllChildrenSpheres(tree.getRoot());
    bool testKnn = correctKnnSearch(tree, data);

    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed = end - start; 

    std::cout << "All data present: " << (allPresent ? "Yes" : "No") << std::endl;
    std::cout << "Leaf nodes at the same level: " << (sameLevel ? "Yes" : "No") << std::endl;
    std::cout << "No exceeding the child limit per node: " << (noExceed ? "Yes" : "No") << std::endl;
    std::cout << "Hypersphere covers all points in leaf nodes: " << (spherePoints ? "Yes" : "No") << std::endl;
    std::cout << "Hypersphere covers all internal node hyperspheres: "
            << (sphereChildren ? "Yes" : "No") << std::endl;
    std::cout << "Performs KNN search: " << (testKnn ? "Yes" : "No") << std::endl;

    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    std::cout << "Happy ending! :D" << std::endl;

    return 0;
}

