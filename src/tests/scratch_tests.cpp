#include "catch.hpp"
#include "geometry.h"

TEST_CASE("Scratch", "[scratch]")
{
    CHECK(__alignof__(Vector) == 16);
    CHECK(sizeof(Vector) == 16);

    CHECK(__alignof__(Point) == 16);
    CHECK(sizeof(Point) == 16);
}
