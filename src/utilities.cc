// Daniel Shervheim, 2019
// danielshervheim.com

#include "utilities.h"

#include "coefficients.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

bool Utilities::DirectoryExists(const char* path)
{
    // Source:
    // https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c

    struct stat info;
    if(stat(path, &info) != 0)
    {
        // Cannot access.
        return false;
    }
    else if (info.st_mode & S_IFDIR)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Utilities::RayCircleIntersection(vec2 ray_origin, vec2 ray_direction,
    vec2 circle_center, double circle_radius,
    double& t_min, double& t_max)
{
    t_min = 0.0;
    t_max = 0.0;

    vec2 circle_to_ray_origin = ray_origin - circle_center;

    double a = vec2::dot(ray_direction, ray_direction);
    double b = 2.0 * vec2::dot(ray_direction, circle_to_ray_origin);
    double c = vec2::dot(circle_to_ray_origin, circle_to_ray_origin) - circle_radius*circle_radius;

    double discriminant = b*b - (4.0 * a * c);

    if (discriminant < 0)
    {
        return false;
    }

    discriminant = sqrt(discriminant);

    double t0 = (-b + discriminant) / (2.0 * a);
    double t1 = (-b - discriminant) / (2.0 * a);

    t_min = fmin(t0, t1);
    t_max = fmax(t0, t1);

    // If t_max is < 0, then the ray points away from the circle.
    return t_max > 0.0;
};

bool Utilities::WriteFloatArrayToFile(std::string file_path, float* array, int array_length)
{
    std::ofstream out(file_path, std::ios::out | std::ios::binary);

    if(!out)
    {
        return false;
    }

    try
    {
        out.write((char*)array, sizeof array[0] * array_length);
        out.close();
    }
    catch(const std::exception &exc)
    {
        (void)exc;
        return false;
    }
    return true;
}

float* Utilities::ConvertDoubleArrayToFloatArray(double* array, int array_length)
{
    float* tmp = new float[array_length];

    for (int i = 0; i < array_length; i++)
    {
        tmp[i] = (float)array[i];
    }

    return tmp;
}

void Utilities::ConvertXyzToRgb(double x, double y, double z,
    double& r, double& g, double& b)
{
    r = x* 3.2406 + y*-1.5372 + z*-0.4986;
    g = x*-0.9689 + y* 1.8758 + z* 0.0415;
    b = x* 0.0557 + y*-0.2040 + z* 1.0570;
}
