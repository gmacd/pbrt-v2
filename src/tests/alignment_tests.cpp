#include "catch.hpp"
#include "geometry.h"
#include "transform.h"

TEST_CASE("Alignment", "[alignment]")
{
    CHECK(__alignof__(Vector) == 16);
    CHECK(sizeof(Vector) == 16);

    CHECK(__alignof__(Point) == 16);
    CHECK(sizeof(Point) == 16);
    
    CHECK(__alignof__(Matrix4x4) == 16);
    CHECK(sizeof(Matrix4x4) == 64);

    CHECK(__alignof__(Transform) == 16);
}
