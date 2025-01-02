# SS-TREE

### Overview

SS-TREE is a high-performance spatial data structure designed for efficient querying and indexing of multi-dimensional data. It provides a tree-based structure that allows fast nearest-neighbor searches, efficient data insertion, and range queries. The tree's design is based on spatial partitioning principles and optimized for high-dimensional spaces, leveraging techniques such as bounding volumes and hierarchical indexing.

This repository contains the implementation of the SS-TREE data structure, with core features including:

- **Dynamic Insertion**: Supports insertion of data points while maintaining balanced tree structure and minimizing search complexity.
- **Bounding Volume Hierarchies**: Uses bounding spheres and other spatial heuristics to optimize search and insertion operations.
- **Efficient KNN Search**: Allows for fast retrieval of the k-nearest neighbors to a given query point.

This repository is intended for research purposes and can be used in various applications involving high-dimensional spatial data, such as machine learning, computer vision, robotics, and more.

### Prerequisites

To build and run the SS-TREE implementation, the following dependencies must be installed:

- **Eigen3**: A C++ template library for linear algebra. It is required for vector and matrix operations within the SS-TREE structure.

You can install **Eigen3** using the following commands:

#### On Ubuntu (or other Debian-based systems):
```bash
sudo apt-get install libeigen3-dev
```
#### On MacOS:
```bash
brew install eigen
```
### Building the Project

The project includes a `Makefile` to simplify the build process. To compile and run the project, execute the following commands:

1. Clone the repository:
    ```bash
    git clone [https://github.com/darkblackv32/ss-tree/]
    cd ss-tree
    ```

2. Install the required dependency, **Eigen3**, by following the instructions in the [Prerequisites](#prerequisites) section.

3. Run the `make` command to compile the project:
    ```bash
    make run
    ```
        
After execution, the test results will be displayed in the terminal, indicating whether the ss-tree implementation passed or failed the various tests.
