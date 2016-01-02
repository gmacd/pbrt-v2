#include "vector.h"

#include "pbrt.h"
#include "geometry.h"


// Define DISABLE_SIMD to force disable any SIMD code
// Undefined at bottom of file
#if defined __AVX2__ && !defined DISABLE_SIMD
#define USE_SIMD
#endif


Vector::Vector(const Point &p):
    x(p.x), y(p.y), z(p.z)
{
    Assert(!HasNaNs());
}

Vector::Vector(const Normal &n):
    x(n.x), y(n.y), z(n.z)
{
    Assert(!n.HasNaNs());
}


#undef USE_SIMD
