#include "catch.hpp"
#include "geometry.h"
#include "transform.h"

TEST_CASE("Alignment", "[alignment]")
{
    CHECK(alignof(Vector) == 16);
    CHECK(sizeof(Vector) == 16);

    CHECK(alignof(Point) == 16);
    CHECK(sizeof(Point) == 16);
    
    CHECK(alignof(Matrix4x4) == 16);
    CHECK(sizeof(Matrix4x4) == 64);

    CHECK(alignof(Transform) == 16);
}
