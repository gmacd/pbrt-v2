#pragma once

#include <xmmintrin.h>

#include "pbrt.h"


// Define DISABLE_SIMD to force disable any SIMD code
// Undefined at bottom of file
#if defined __AVX2__ && !defined DISABLE_SIMD
#define USE_SIMD
#endif

typedef __m128 float4_t;


class
__attribute__ ((aligned(16)))
Vector
{
public:
    union
    {
        struct { float x, y, z; };
        float4_t _vec;
    };


    Vector()
    {
        x = y = z = 0.f;
    }
    
    Vector(float xx, float yy, float zz):
        x(xx), y(yy), z(zz)
    {
        Assert(!HasNaNs());
    }
    
    Vector(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y; z = v.z;
    }

    Vector(float4_t v):
        _vec(v)
    {
    }

    
    explicit Vector(const Point &p);
    explicit Vector(const Normal &n);

    
    Vector &operator=(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x = v.x; y = v.y; z = v.z;
        return *this;
    }

    bool operator==(const Vector &v) const
    {
        return x == v.x && y == v.y && z == v.z;
    }
    
    bool operator!=(const Vector &v) const
    {
        return x != v.x || y != v.y || z != v.z;
    }
    
    Vector operator+(const Vector &v) const
    {
        Assert(!v.HasNaNs());
#ifdef USE_SIMD
        return Vector(_mm_add_ps(_vec, v._vec));
#else
        return Vector(x + v.x, y + v.y, z + v.z);
#endif
    }
    
    Vector& operator+=(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    
    Vector operator-(const Vector &v) const
    {
        Assert(!v.HasNaNs());
        return Vector(x - v.x, y - v.y, z - v.z);
    }
    
    Vector& operator-=(const Vector &v)
    {
        Assert(!v.HasNaNs());
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    
    Vector operator*(float f) const
    {
        return Vector(f*x, f*y, f*z);
    }
    
    Vector &operator*=(float f)
    {
        Assert(!isnan(f));
        x *= f; y *= f; z *= f;
        return *this;
    }
    
    Vector operator/(float f) const
    {
        Assert(f != 0);
        float inv = 1.f / f;
        return Vector(x * inv, y * inv, z * inv);
    }
    
    Vector &operator/=(float f)
    {
        Assert(f != 0);
        float inv = 1.f / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    
    Vector operator-() const
    {
        return Vector(-x, -y, -z);
    }
    
    float operator[](int i) const
    {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    
    float &operator[](int i)
    {
        Assert(i >= 0 && i <= 2);
        return (&x)[i];
    }
    

    bool HasNaNs() const
    {
        return isnan(x) || isnan(y) || isnan(z);
    }
        
    float LengthSquared() const
    {
#ifdef USE_SIMD
        return _mm_cvtss_f32(_mm_dp_ps(_vec, _vec, 0x71));
#else
        return x*x + y*y + z*z;
#endif
    }
    
    float Length() const
    {
#ifdef USE_SIMD
        return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(_vec, _vec, 0x71)));
#else
        return sqrtf(LengthSquared());
#endif
    }
};


#undef USE_SIMD
