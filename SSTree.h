#ifndef SSTREE_H
#define SSTREE_H

#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>
#include <queue>
#include "Point.h"
#include "Data.h"

class SSNode {
    size_t maxPointsPerNode;
    Point centroid;
    float radius;
    bool isLeaf;

    SSNode* parent;
    std::vector<SSNode*> children;
    std::vector<Data*> _data;

    // For searching
    SSNode* findClosestChild(const Point& target);

    // For insertion
    void updateBoundingEnvelope();
    size_t directionOfMaxVariance();
    std::pair<SSNode*, SSNode*> split();
    size_t findSplitIndex(size_t coordinateIndex);
    std::vector<Point> getEntriesCentroids() const;
    size_t minVarianceSplit(const std::vector<float>& values);

public:
    explicit SSNode(const Point& centroid, float radius=0.0f, bool isLeaf=true, SSNode* parent=nullptr , size_t M = 4)
        : centroid(centroid), radius(radius), isLeaf(isLeaf), parent(parent) , maxPointsPerNode(M){}

    // Checks if a point is inside the bounding sphere
    bool intersectsPoint(const Point& point) const;

    // Getters
    const Point& getCentroid() const { return centroid; }
    float getRadius() const { return radius; }
    const std::vector<SSNode*>& getChildren() const { return children; }
    const std::vector<Data*>& getData () const { return    _data; }
    bool getIsLeaf() const { return isLeaf; }
    SSNode* getParent() const { return parent; }

    // Insertion
    SSNode* searchParentLeaf(SSNode* node, const Point& target);
    std::pair<SSNode*, SSNode*> insert(SSNode*& node, Data* data);

    // Search
    SSNode* search(SSNode* node, Data* _data);

    friend class SSTree;
};

class SSTree {
    SSNode* root;
    size_t maxPointsPerNode;

public:
    SSTree(size_t maxPointsPerNode) : maxPointsPerNode(maxPointsPerNode), root(nullptr) {}

    void insert(Data* _data);
    SSNode* search(Data* _data);

    SSNode * getRoot() const {
        return root;
    };

    std::vector<Data*> knn(const Point& query, size_t k) const;
};

#endif // SSTREE_H