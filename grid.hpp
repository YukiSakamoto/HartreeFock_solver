#pragma once

#include <cmath>
#include "Eigen/Core"
#include "common.hpp"
#include "system.hpp"
#include "orbitals.hpp"
#include <cstdio>
#include <string>
#include <vector>

#include "lebedev_194.h"

struct GridPoint {
    GridPoint(const Vector3Real &x, REAL weight = 1.): pos(x), weight1(weight), weight2(1.)
    {}

    REAL weight1;
    REAL weight2;   // Becke scheme
    Vector3Real pos;
    REAL radial_;
};

struct AtomicGrid {
    AtomicGrid(const Vector3Real &pos, REAL weight = 1.) :pos(pos), weight(weight)
    {}

    std::vector<GridPoint> grid_points;
    Vector3Real pos;
    REAL weight;
    
    void push_back(const GridPoint& gp) {
        this->grid_points.push_back(gp);
    }
    void append(const GridPoint& gp) {
        this->grid_points.push_back(gp);
    }
    size_t size() const {   
        return grid_points.size();  
    }
};

std::vector<AtomicGrid>
GenerateGrid(const System &system);
void
PrintGridWeight(const AtomicGrid &grid);
