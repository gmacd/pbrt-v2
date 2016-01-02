#include "catch.hpp"

#define DISABLE_SIMD
#include "vector.h"


TEST_CASE("BasicOperations", "[vector]")
{
    Vector v1(2, 2, 2);
    Vector v2(3, 3, 3);

    auto vRes = v1+v2;
    CHECK(vRes == Vector(5, 5, 5));

    vRes = v1-v2;
    CHECK(vRes == Vector(-1, -1, -1));

    vRes = -v1;
    CHECK(vRes == Vector(-2, -2, -2));
}

TEST_CASE("Length", "[vector]")
{
    Vector v(0, 0, 0);
    
    auto len = v.Length();
    auto lenSq = v.LengthSquared();

    CHECK(len == 0);
    CHECK(lenSq == 0);

    v = Vector(1, 1, 1);
    len = v.Length();
    lenSq = v.LengthSquared();

    CHECK(len == Approx(1.73205));
    CHECK(lenSq == Approx(3));

    v = Vector(-1, -1, -1);
    len = v.Length();
    lenSq = v.LengthSquared();

    CHECK(len == Approx(1.73205));
    CHECK(lenSq == Approx(3));
}
