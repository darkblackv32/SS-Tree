#include "SSTree.h"

float calculateMean(const std::vector<Point>& centroids, size_t dimension) {
    float sum = 0.0f;
    for (const auto& point : centroids) {
        sum += point[dimension];
    }
    return sum / centroids.size();
}

float calculateVariance(const std::vector<Point>& centroids, size_t dimension) {
    
    float meanValue = calculateMean(centroids, dimension);

    float varianceSum = 0.0f;
    for (const auto& point : centroids) {
        float deviation = point[dimension] - meanValue;
        varianceSum += deviation * deviation; 
    }

    return varianceSum / static_cast<float>(centroids.size());
}

/**
 * intersectsPoint
 * Checks if a point is inside the bounding sphere of the node.
 * @param point: Point to verify.
 * @return bool: Returns true if the point is inside the sphere; otherwise, false.
 */

bool SSNode::intersectsPoint(const Point& point) const {
    return Point::distance(centroid, point) <= radius;
}

/**
 * findClosestChild
 * Finds the closest child to a given point.
 * @param target: The target point to find the closest child.
 * @return SSNode*: Returns a pointer to the closest child.
 */

SSNode* SSNode::findClosestChild(const Point& target) {
    return *std::min_element(children.begin(), children.end(), [&target](SSNode* a, SSNode* b) {
        return a->getCentroid().distance(target) < b->getCentroid().distance(target);
    });
}


/**
 * updateBoundingEnvelope
 * Updates the centroid and radius of the node based on internal nodes or data.
 */

void SSNode::updateBoundingEnvelope() {
    std::vector<Point> entryCentroids = getEntriesCentroids();

    for (size_t dim = 0; dim < DIM; dim++) {
        this->centroid[dim] = calculateMean(entryCentroids, dim);
    }

    float maxRadius = 0.0f;

    if (this->isLeaf) {
        for (const auto& entry : this->_data) {
            float distanceToCentroid = Point::distance(this->centroid, entry->getEmbedding());
            maxRadius = std::max(maxRadius, distanceToCentroid);
        }
    } else {
        for (const auto& child : this->children) {
            float distanceToChild = Point::distance(this->centroid, child->centroid) + child->radius;
            maxRadius = std::max(maxRadius, distanceToChild);
        }
    }

    this->radius = maxRadius;
}

/**
 * directionOfMaxVariance
 * Calculates and returns the index of the direction of maximum variance.
 * @return size_t: Index of the direction of maximum variance.
 */

size_t SSNode::directionOfMaxVariance() {
    float highestVariance = 0.0f;
    size_t maxVarianceDirection = 0;

    const auto centroids = this->getEntriesCentroids();

    for (size_t dim = 0; dim < DIM; ++dim) {
        float currentVariance = calculateVariance(centroids, dim);

        if (currentVariance > highestVariance) {
            highestVariance = currentVariance;
            maxVarianceDirection = dim;
        }
    }

    return maxVarianceDirection;
}

/**
 * split
 * Splits the node and returns the newly created node.
 * Implementation similar to an R-tree.
 * @return SSNode*: Pointer to the new node created by the split.
 */

std::pair<SSNode*, SSNode*> SSNode::split() {
    size_t splitDimension = directionOfMaxVariance(); 
    size_t splitIndex = findSplitIndex(splitDimension);

    SSNode* leftNode = new SSNode(centroid, radius, isLeaf, this, maxPointsPerNode);
    SSNode* rightNode = new SSNode(centroid, radius, isLeaf, this, maxPointsPerNode);

    if (isLeaf) {
        leftNode->_data.assign(_data.begin(), _data.begin() + splitIndex);
        rightNode->_data.assign(_data.begin() + splitIndex, _data.end());
    } else {
        leftNode->children.assign(children.begin(), children.begin() + splitIndex);
        rightNode->children.assign(children.begin() + splitIndex, children.end());

        for (auto* child : leftNode->children) {
            child->parent = leftNode;
        }
        for (auto* child : rightNode->children) {
            child->parent = rightNode;
        }
    }

    leftNode->updateBoundingEnvelope();
    rightNode->updateBoundingEnvelope();

    this->updateBoundingEnvelope();

    return {leftNode, rightNode};
}


/**
 * findSplitIndex
 * Finds the split index on a specific coordinate.
 * @param coordinateIndex: Index of the coordinate to find the split index.
 * @return size_t: Split index.
 */

size_t SSNode::findSplitIndex(size_t coordinateIndex) {
    std::vector<float> coordinateValues;

    if (isLeaf) {
        std::sort(_data.begin(), _data.end(),
            [coordinateIndex](const Data* lhs, const Data* rhs) {
                return lhs->getEmbedding()[coordinateIndex] < rhs->getEmbedding()[coordinateIndex];
            });

        for (const auto& entry : _data) {
            coordinateValues.push_back(entry->getEmbedding()[coordinateIndex]);
        }
    } else {

        std::sort(children.begin(), children.end(),
            [coordinateIndex](const SSNode* lhs, const SSNode* rhs) {
                return lhs->getCentroid()[coordinateIndex] < rhs->getCentroid()[coordinateIndex];
            });

        for (const auto& child : children) {
            coordinateValues.push_back(child->getCentroid()[coordinateIndex]);
        }
    }

    return minVarianceSplit(coordinateValues);
}

/**
 * getEntriesCentroids
 * Returns the centroids of the entries.
 * These centroids can be points stored in the leaves or the centroids of child nodes in internal nodes.
 * @return std::vector<Point>: Vector of entry centroids.
 */

std::vector<Point> SSNode::getEntriesCentroids() const {
    std::vector<Point> centroids;
    
    if (isLeaf) {
        for (const auto& data : _data) {
            centroids.push_back(data->getEmbedding());
        }
    } else {
        for (const auto& child : children) {
            centroids.push_back(child->getCentroid());
        }
    }
    
    return centroids;
}


/**
 * minVarianceSplit
 * Finds the optimal split index for a list of values such that the sum of variances of the two resulting partitions is minimized.
 * @param values: Vector of values to find the minimum variance index.
 * @return size_t: Index of minimum variance.
 */

size_t SSNode::minVarianceSplit(const std::vector<float>& values) {
    int M = maxPointsPerNode, m = 1;
    float min_s = std::numeric_limits<float>::max();
    size_t idx_min = 0;

    for (int i = m; i <= M - m; i++) {
        float meanLeft = std::accumulate(values.begin(), values.begin() + i, 0.0f) / i;
        float meanRight = std::accumulate(values.begin() + i, values.end(), 0.0f) / (values.size() - i);

        float varLeft = 0.0f, varRight = 0.0f;
        for (size_t j = 0; j < i; j++) {
            varLeft += std::pow(values[j] - meanLeft, 2);
        }
        for (size_t j = i; j < values.size(); j++) {
            varRight += std::pow(values[j] - meanRight, 2);
        }

        float sumVariance = varLeft + varRight;
        if (sumVariance < min_s) {
            min_s = sumVariance;
            idx_min = i;
        }
    }

    return idx_min;
}

/**
 * searchParentLeaf
 * Searches for the appropriate leaf node to insert a point.
 * @param node: Node from which to start the search.
 * @param target: Target point for the search.
 * @return SSNode*: Appropriate leaf node for insertion.
 */

SSNode* SSNode::searchParentLeaf(SSNode* node, const Point& target) {
    if(node->isLeaf ) return node;
    return searchParentLeaf(node->findClosestChild(target) , target);
}

/**
 * insert
 * Inserts data into the node, splitting it if necessary.
 * @param node: Node where the insertion will take place.
 * @param _data: Data to insert.
 * @return SSNode*: New root node if split occurred, otherwise nullptr.
 */

std::pair<SSNode*, SSNode*> SSNode::insert(SSNode*& node, Data* data) {
    if (node->isLeaf) {
        if (std::find(node->_data.begin(), node->_data.end(), data) != node->_data.end()) {
            return {nullptr, nullptr};
        }

        node->_data.push_back(data);
        node->updateBoundingEnvelope();

        if (node->_data.size() <= node->maxPointsPerNode) {
            return {nullptr, nullptr};
        }

        return node->split();
    }

    SSNode* closestChild = node->findClosestChild(data->getEmbedding());

    auto [leftSplit, rightSplit] = insert(closestChild, data);

    if (!leftSplit && !rightSplit) {
        node->updateBoundingEnvelope();
        return {nullptr, nullptr};
    }

    auto it = std::find(node->children.begin(), node->children.end(), closestChild);
    
    if (it != node->children.end()) {
        node->children.erase(it);
    }

    node->children.push_back(leftSplit);
    node->children.push_back(rightSplit);

    node->updateBoundingEnvelope();

    if (node->children.size() <= node->maxPointsPerNode) {
        return {nullptr, nullptr};
    }

    return node->split();
}


/**
 * search
 * Searches for a specific data in the tree.
 * @param node: Node from which to start the search.
 * @param _data: Data to search for.
 * @return SSNode*: Node containing the data (or nullptr if not found).
 */

SSNode* SSNode::search(SSNode* node, Data* _data) {
    if(node->isLeaf) {
        for(auto & point : node->_data) {
            if(point == _data) {
                return node;
            }
        }
    }

    else {
        for(auto & child : node->children) {
            if(child != nullptr && child->intersectsPoint(_data->getEmbedding())) {
                SSNode* ans = search(child , _data);
                if(ans != nullptr) 
                    return ans;
            }
        }
    }
    return nullptr;
}



/**
 * insert
 * Inserts data into the tree.
 * @param _data: Data to insert.
 */

void SSTree::insert(Data* _data) {
    if (root == nullptr) root = new SSNode(_data->getEmbedding(),0,true,nullptr,maxPointsPerNode);
    auto p = root->insert(root , _data);
    SSNode* n1 = p.first; SSNode*n2 = p.second;
    if (n1 != nullptr) {
        root = new SSNode(_data->getEmbedding(), 0 , true , nullptr , maxPointsPerNode);
        root->children.push_back(n1);
        root->children.push_back(n2);
        root->isLeaf = false;
    }
}


/**
 * search
 * Searches for a specific data in the tree.
 * @param _data: Data to search for.
 * @return SSNode*: Node containing the data (or nullptr if not found).
 */

SSNode* SSTree::search(Data* _data) {
    return root->search(root,_data);
}


/**
 * knn-search
 * Returns the k nearest neighbors.
 * @param k: number of neighbors
 * @param query: point from which to find the k nearest neighbors
 * @return std::vector<Data*>: List containing the k nearest neighbors
 */

std::vector<Data*> SSTree::knn(const Point& query, size_t k) const {
    if (!root) {
        return {}; 
    }

    auto compare = [](const std::pair<const SSNode*, float>& a, const std::pair<const SSNode*, float>& b) {
        return a.second > b.second;
    };

    std::priority_queue<std::pair<const SSNode*, float>, std::vector<std::pair<const SSNode*, float>>, decltype(compare)> nodeQueue(compare);

    auto dataCompare = [&query](const Data* a, const Data* b) {
        return a->getEmbedding().distance(query) < b->getEmbedding().distance(query);
    };

    std::priority_queue<Data*, std::vector<Data*>, decltype(dataCompare)> nearestNeighbors(dataCompare);

    nodeQueue.emplace(root, query.distance(root->getCentroid()) - root->getRadius());

    while (!nodeQueue.empty()) {
        auto [currentNode, nodeDistance] = nodeQueue.top();
        nodeQueue.pop();

        if (!nearestNeighbors.empty() && nodeDistance > nearestNeighbors.top()->getEmbedding().distance(query)) {
            continue;
        }

        if (currentNode->getIsLeaf()) {
            for (const auto& data : currentNode->getData()) {
                float dataDistance = data->getEmbedding().distance(query);
                if (nearestNeighbors.size() < k) {
                    nearestNeighbors.push(data);
                } else if (dataDistance < nearestNeighbors.top()->getEmbedding().distance(query)) {
                    nearestNeighbors.pop();
                    nearestNeighbors.push(data);
                }
            }
        } else {
            for (const auto& child : currentNode->getChildren()) {
                float childDistance = query.distance(child->getCentroid()) - child->getRadius();
                if (!nearestNeighbors.empty() && childDistance > nearestNeighbors.top()->getEmbedding().distance(query)) {
                    continue;
                }
                nodeQueue.emplace(child, childDistance);
            }
        }
    }

    std::vector<Data*> ans;
    while (!nearestNeighbors.empty()) {
        ans.push_back(nearestNeighbors.top());
        nearestNeighbors.pop();
    }

    std::reverse(ans.begin(), ans.end());

    return ans;
}