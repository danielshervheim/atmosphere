// Daniel Shervheim, 2019
// danielshervheim.com

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "vec2.h"

#include <string>

class Utilities
{
public:
    // Returns whether or not a ray intersects a circle. If the ray intersects
    // the circle, then the `t_min` and `t_max` parameters contain the hit points
    // for the parametric ray equation.
    static bool RayCircleIntersection(vec2 ray_origin, vec2 ray_direction,
        vec2 circle_center, double circle_radius,
        double& t_min, double& t_max);

    // Writes the float array to the file path. Returns whether or not the write
    // was successful.
    static bool WriteFloatArrayToFile(std::string file_path, float* array, int array_length);

    // Allocates memory for the given double array and converts it into a float
    // array. You must manually delete the memory once you are done with it by
    // calling `delete [] ARRAYNAME`.
    static float* ConvertDoubleArrayToFloatArray(double* array, int array_length);

    static void ConvertXyzToRgb(double x, double y, double z, double& r, double& g, double& b);
};

#endif  // UTILITIES_H_
