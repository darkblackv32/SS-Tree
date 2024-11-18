#include "SSTree.h"


/**
 * intersectsPoint
 * Verifica si un punto está dentro de la esfera delimitadora del nodo.
 * @param point: Punto a verificar.
 * @return bool: Retorna true si el punto está dentro de la esfera, de lo contrario false.
 */
bool SSNode::intersectsPoint(const Point& point) const {
    return Point::distance(centroid, point) <= radius;
}

/**
 * findClosestChild
 * Encuentra el hijo más cercano a un punto dado.
 * @param target: El punto objetivo para encontrar el hijo más cercano.
 * @return SSNode*: Retorna un puntero al hijo más cercano.
 */


SSNode* SSNode::findClosestChild(const Point& target) {
    return *std::min_element(children.begin(), children.end(), [&target](SSNode* a, SSNode* b) {
        return a->getCentroid().distance(target) < b->getCentroid().distance(target);
    });
}

/**
 * updateBoundingEnvelope
 * Actualiza el centroide y el radio del nodo basándose en los nodos internos o datos.
 */
void SSNode::updateBoundingEnvelope() {
    Point n_cent = Point::Zero(); 
    float n_rad = 0.0f;

    if (isLeaf) {
        for (const auto& data : _data) {
            n_cent += data->getEmbedding();
        }
        n_cent /= _data.size();

        for (const auto& data : _data) {
            n_rad = std::max(n_rad, n_cent.distance(data->getEmbedding()));
        }
    } else {
        for (const auto& child : children) {
            n_cent += child->getCentroid();
        }
        n_cent /= children.size();

        for (const auto& child : children) {
            n_rad = std::max(n_rad, 
                                 n_cent.distance(child->getCentroid()) + child->getRadius());
        }
    }

    centroid = n_cent;
    radius = n_rad;
}


/**
 * directionOfMaxVariance
 * Calcula y retorna el índice de la dirección de máxima varianza.
 * @return size_t: Índice de la dirección de máxima varianza.
 */
size_t SSNode::directionOfMaxVariance() {
    
    std::vector<float> med_vals = std::vector<float>(DIM, 0.0f);

    if (isLeaf) {
        for (auto data : _data) {
            for (size_t i = 0; i < DIM; i++) {
                med_vals[i] += data->getEmbedding()[i];
            }
        }

        for (size_t i = 0; i < DIM; i++) {
            med_vals[i] /= _data.size();
        }
    } else {
        for (auto child : children) {
            for (size_t i = 0; i < DIM; i++) {
                med_vals[i] += child->getCentroid()[i];
            }
        }

        for (size_t i = 0; i < DIM; i++) {
            med_vals[i] /= children.size();
        }
    }

    std::vector<float> variances = std::vector<float>(DIM, 0.0f);

    if (isLeaf) {
        for (auto data : _data) {
            for (size_t i = 0; i < DIM; i++) {
                variances[i] += std::pow(data->getEmbedding()[i] - med_vals[i], 2);
            }
        }
    } else {
        for (auto child : children) {
            for (size_t i = 0; i < DIM; i++) {
                variances[i] += std::pow(child->getCentroid()[i] - med_vals[i], 2);
            }
        }
    }

    size_t min_idxvar = 0;

    for (size_t i = 1; i < DIM; i++) {
        if (variances[i] > variances[min_idxvar]) {
            min_idxvar = i;
        }
    }

    return min_idxvar;
}




/**
 * split
 * Divide el nodo y retorna el nuevo nodo creado.
 * Implementación similar a R-tree.
 * @return SSNode*: Puntero al nuevo nodo creado por la división.
 */

SSNode* SSNode::split() {
    size_t dir_split = directionOfMaxVariance();

    // Determine the split index
    size_t idx_split = findSplitIndex(dir_split);

    SSNode* n_node = new SSNode(Point::Zero(), 0.0f, isLeaf, parent);

    if (isLeaf) {
        auto s_mid = _data.begin() + idx_split;
        std::nth_element(_data.begin(), s_mid, _data.end(),
            [dir_split](const Data* a, const Data* b) {
                return a->getEmbedding()[dir_split] < b->getEmbedding()[dir_split];
            });

        // Assign and erase the second half
        n_node->_data.assign(s_mid, _data.end());
        _data.erase(s_mid, _data.end());
    } else {
        auto s_mid = children.begin() + idx_split;
        std::nth_element(children.begin(), s_mid, children.end(),
            [dir_split](const SSNode* a, const SSNode* b) {
                return a->getCentroid()[dir_split] < b->getCentroid()[dir_split];
            });

        n_node->children.assign(s_mid, children.end());
        children.erase(s_mid, children.end());
    }

    updateBoundingEnvelope();
    n_node->updateBoundingEnvelope();

    return n_node;
}



/**
 * findidx_split
 * Encuentra el índice de división en una coordenada específica.
 * @param coordinateIndex: Índice de la coordenada para encontrar el índice de división.
 * @return size_t: Índice de la división.
 */
size_t SSNode::findSplitIndex(size_t coordinateIndex) {


    std::vector<float> vals;

    if (isLeaf) {
        for (const auto& data : _data) {
            vals.push_back(data->getEmbedding()[coordinateIndex]);
        }
    } else {
        for (const auto& child : children) {
            vals.push_back(child->getCentroid()[coordinateIndex]);
        }
    }

    return minVarianceSplit(vals);
}


/**
 * getEntriesCentroids
 * Devuelve los centroides de las entradas.
 * Estos centroides pueden ser puntos almacenados en las hojas o los centroides de los nodos hijos en los nodos internos.
 * @return std::vector<Point>: Vector de centroides de las entradas.
 */
std::vector<Point> SSNode::getEntriesCentroids() {
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
 * Encuentra el índice de división óptimo para una lista de valores, de tal manera que la suma de las varianzas de las dos particiones resultantes sea mínima.
 * @param vals: Vector de valores para encontrar el índice de mínima varianza.
 * @return size_t: Índice de mínima varianza.
 */

size_t SSNode::minVarianceSplit(const std::vector<float>& vals) {
    int M = maxPointsPerNode, m = 1;
    float min_s = std::numeric_limits<float>::max();
    size_t idx_min = 0;

    for (int i = m; i <= M - m; i++) {
        float meanLeft = std::accumulate(vals.begin(), vals.begin() + i, 0.0f) / i;
        float meanRight = std::accumulate(vals.begin() + i, vals.end(), 0.0f) / (vals.size() - i);

        float varLeft = 0.0f, varRight = 0.0f;
        for (size_t j = 0; j < i; j++) {
            varLeft += std::pow(vals[j] - meanLeft, 2);
        }
        for (size_t j = i; j < vals.size(); j++) {
            varRight += std::pow(vals[j] - meanRight, 2);
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
 * Busca el nodo hoja adecuado para insertar un punto.
 * @param node: Nodo desde el cual comenzar la búsqueda.
 * @param target: Punto objetivo para la búsqueda.
 * @return SSNode*: Nodo hoja adecuado para la inserción.
 */
SSNode* SSNode::searchParentLeaf(SSNode* node, const Point& target) {
    if (!node) return nullptr;
    
    if (node->isLeaf) return node;
    
    SSNode* closestChild = node->findClosestChild(target);
    return searchParentLeaf(closestChild, target);
}

/**
 * insert
 * Inserta un dato en el nodo, dividiéndolo si es necesario.
 * @param node: Nodo donde se realizará la inserción.
 * @param _data: Dato a insertar.
 * @return SSNode*: Nuevo nodo raíz si se dividió, de lo contrario nullptr.
 */
SSNode* SSNode::insert(SSNode* node, Data* data) {
    if(node->getIsLeaf()){ 
        if (node->maxPointsPerNode > node->_data.size()) { 
            
            node->_data.push_back(data);

            node->updateBoundingEnvelope();

            return nullptr;

        } else { 

            SSNode* s_auxpr = node->split(); 

            return s_auxpr;
        }

    } else { 
        SSNode* s_auxpr = node->findClosestChild(data->getEmbedding());

        SSNode* s_ret = s_auxpr->insert(s_auxpr, data);

        if (s_ret == nullptr) { 

            node->updateBoundingEnvelope();

            return nullptr;

        } else { 
            if (node->children.size() < node->maxPointsPerNode){ 

                node->children.push_back(s_ret);

                node->updateBoundingEnvelope();

                return nullptr;
            } else { 
                node->children.push_back(s_ret);

                s_auxpr = node->split();

                return s_auxpr;
            }
        }
    }  
}


/**
 * search
 * Busca un dato específico en el árbol.
 * @param node: Nodo desde el cual comenzar la búsqueda.
 * @param _data: Dato a buscar.
 * @return SSNode*: Nodo que contiene el dato (o nullptr si no se encuentra).
 */
SSNode* SSNode::search(SSNode* node, Data* data) {
    if (node == nullptr) return nullptr;

    if (node->getIsLeaf()){ 

        for (int i = 0; i < node->children.size(); i++){
            SSNode* ans = node->search(node->children[i], data);

            if (ans != nullptr){
                return ans;
            }
        }
        
        return nullptr;
    } else { 
        Data* ans = nullptr;

        for (int i = 0; i < node->_data.size(); i++) {
            if (data->getPath() == node->_data[i]->getPath()){
                ans = node->_data[i];
                break;
            }
        }

        if (ans != nullptr) return node;
    
        return nullptr;
    }
}


/**
 * insert
 * Inserta un dato en el árbol.
 * @param _data: Dato a insertar.
 */
void SSTree::insert(Data* data) {
    if (this->root == nullptr) {
        this->root = new SSNode(data->getEmbedding());
        this->root->insert(this->root, data);
    } else {
        this->root->insert(this->root, data);
    }
}


/**
 * search
 * Busca un dato específico en el árbol.
 * @param _data: Dato a buscar.
 * @return SSNode*: Nodo que contiene el dato (o nullptr si no se encuentra).
 */
SSNode* SSTree::search(Data* data) {
    return root->search(root, data);
}