#include "orbit.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>


int main(void)
{
    double stdGrvPrm = 3.986004418E+14, seed = 521 * 1000;
    Vector3D pos(seed * 12, seed * 4, seed * 3),
        vel(-2435, 7305, 1200);

    Orbit integrated(stdGrvPrm, pos, vel, 10.0, false);
    Orbit predicted(stdGrvPrm, pos, vel, (std::size_t)1048576);
    predicted.saveElements("predicted");

    double meanDist = predicted.getMeanDist(integrated);
    std::cout << "\nMean Distance: " << std::setprecision(14)
              << meanDist << std::endl;

    Vector3D newPos = predicted.getPosVector(10000);
    std::cout << "\nnewpos = Vector3D(" << newPos.x << ", "
              << newPos.y << ",  "
              << newPos.z << ")" << std::endl;

    Vector3D newVel = predicted.getVelVector(10000);
    std::cout << "newvel = Vector3D(" << newVel.x << ", "
              << newVel.y << ",  "
              << newVel.z << ")" << std::endl;

    std::system("pause");

    return 0;
}