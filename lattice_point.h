/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   lattice_point.h
 * Author: Marshall
 *
 * Created on March 4, 2022, 10:17 AM
 */

#ifndef LATTICE_POINT_H
#define LATTICE_POINT_H

#include <algorithm>
#include "structures.h"

class LatticePoint {
public:
    int x = 0;
    int y = 0;
    int z = 0;
    LatticePoint(int x = 0, int y = 0, int z = 0) : x{x}, y{y}, z{z} {}
    int norm(){
        return (x * x + y * y + z * z);
    }
    int norm() const{
        return (x * x + y * y + z * z);
    }
    bool operator < (const LatticePoint& lattice) const {
        return (norm() < lattice.norm());
    }
    bool operator > (const LatticePoint& lattice) const {
        return (norm() > lattice.norm());
    }
    
    void print() {
        printf("(%d, %d, %d), norm = %d\n", x, y, z, norm());
    }
};

class Grid {
public:
    LatticePoint* lattice;
    int numLattice;
    int symDimSum;
    vec3<int> symmetryAxes;
    Grid(int boxLength, vec3<int> symmetryAxes) : symmetryAxes{symmetryAxes} {
        assert(boxLength % 2 == 1);
        int n = boxLength / 2;
        int ii, i, j, k, latticeIndex;
        vec3<int> tempSymDim = symmetryAxes;
        
        symDimSum = symmetryAxes.x + symmetryAxes.y + symmetryAxes.z;
        
        if (symDimSum == 1){
            numLattice = n + 1;
        }
        else if (symDimSum == 2){
            numLattice = n * (n + 1) + 1;
        }
        else if (symDimSum == 3){
            numLattice = n * (n * n + n + 1) + 1;
        }
        else {
            numLattice = 1;
        }

        lattice = new LatticePoint[numLattice];
        
        latticeIndex = 0;
        
        for (ii = 0; ii < symDimSum; ii++){
            i = tempSymDim.x;
            do {
                j = tempSymDim.y;
                do {
                    k = tempSymDim.z;
                    do {
                        lattice[latticeIndex] = LatticePoint(i, j, k);
                        latticeIndex++;
                        k++;
                    } while (tempSymDim.z != 0 && k <= n);
                    j++;
                } while (tempSymDim.y != 0 && j <= n);
                i++;
            } while (tempSymDim.x != 0 && i <= n);

            if (tempSymDim.x == 0){
                if (tempSymDim.y == 0){
                    tempSymDim.z = 0;
                }
                else {
                    tempSymDim.y = 0;
                }
            }
            else {
                tempSymDim.x = 0;
            }
        }
        lattice[latticeIndex] = LatticePoint();

        std::sort(lattice, lattice + numLattice, std::greater<LatticePoint>());
    }
    ~Grid(){
        delete[] lattice;
    }
};

#endif /* LATTICE_POINT_H */

