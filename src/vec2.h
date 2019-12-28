// Daniel Shervheim, 2019
// danielshervheim.com

#ifndef VEC2_H_
#define VEC2_H_


#include <cmath>

class vec2
{
public:
    double x;
    double y;

    vec2()
    {
        x = 0.0;
        y = 0.0;
    };

    vec2(double x0, double y0)
    {
        x = x0;
        y = y0;
    };

    double magnitude()
    {
        return sqrt(x*x + y*y);
    };

    vec2 normalized()
    {
        return vec2(x, y) / magnitude();
    };

    static double dot(vec2 a, vec2 b)
    {
        return (a.x * b.x) + (a.y * b.y);
    }

    static double magnitude(vec2 a)
    {
        return a.magnitude();
    }

    static vec2 lerp(vec2 a, vec2 b, double t)
    {
        return a + (b-a)*t;
    }

    static vec2 normalized(vec2 a)
    {
        return a.normalized();
    }

    // Addition operator overload.
    vec2 operator + (const vec2& obj)
    {
        vec2 res;
        res.x = x + obj.x;
        res.y = y + obj.y;
        return res;
    };

    // Subtraction operator overload.
    vec2 operator - (const vec2& obj)
    {
        vec2 res;
        res.x = x - obj.x;
        res.y = y - obj.y;
        return res;
    };

    // Multiplication (scalar) operator overload.
    vec2 operator * (const double val)
    {
        vec2 res;
        res.x = x * val;
        res.y = y * val;
        return res;
    };

    // Division (scalar) operator overload.
    vec2 operator / (const double val)
    {
        vec2 res;
        res.x = x / val;
        res.y = y / val;
        return res;
    };
};

#endif  // VEC2_H_
