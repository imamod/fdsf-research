#include "Common.h"
#include "ZetaFunction.h"

TEST_CASE("bernulli") {
    SECTION("odd") {
        CHECK(0 == bernulli(1));
        CHECK(0 == bernulli(3));
        CHECK(0 == bernulli(5));
        CHECK(0 == bernulli(7));
        CHECK(0 == bernulli(9));
        CHECK(0 == bernulli(11));
        CHECK(0 == bernulli(13));
        CHECK(0 == bernulli(15));
        CHECK(0 == bernulli(17));
        CHECK(0 == bernulli(19));
        CHECK(0 == bernulli(21));
        CHECK(0 == bernulli(23));
    }
    SECTION("even") {
        CHECK(_1 / 6 == bernulli(2));
        CHECK(-_1 / 30 == bernulli(4));
        CHECK(_1 / 42 == bernulli(6));
        CHECK(-_1 / 30 == bernulli(8));
        CHECK(_5 / 66 == bernulli(10));
        CHECK(-_691 / 2730 == bernulli(12));
        CHECK(_7 / 6 == bernulli(14));
        CHECK(-_3617 / 510 == bernulli(16));
        CHECK(_43867 / 798 == bernulli(18));
        CHECK(-_174611 / 330 == bernulli(20));
        CHECK(_854513 / 138 == bernulli(22));
        CHECK(-_236364091 / 2730 == bernulli(24));
    }
}

TEST_CASE("dzeta") {
    SECTION("odd") {
        CHECK(0 == zetaFunction(1));
        CHECK(0 == zetaFunction(3));
        CHECK(0 == zetaFunction(5));
        CHECK(0 == zetaFunction(7));
        CHECK(0 == zetaFunction(9));
        CHECK(0 == zetaFunction(11));
        CHECK(0 == zetaFunction(13));
        CHECK(0 == zetaFunction(15));
        CHECK(0 == zetaFunction(17));
        CHECK(0 == zetaFunction(19));
        CHECK(0 == zetaFunction(21));
        CHECK(0 == zetaFunction(23));
    }
    SECTION("even") {
        CHECK(_1 / 6 == zetaFunction(2) / pow(pi(), 2));
        CHECK(_1 / 90 == zetaFunction(4) / pow(pi(), 4));
        CHECK(_1 / 945 == zetaFunction(6) / pow(pi(), 6));
        CHECK(_1 / 9450 == zetaFunction(8) / pow(pi(), 8));
        CHECK(_1 / 93550 == zetaFunction(10) / pow(pi(), 10));
        CHECK(_691 / 638512875 == zetaFunction(12) / pow(pi(), 12));
        CHECK(_2 / 18243225 == zetaFunction(14) / pow(pi(), 14));
        CHECK(_3617 / 325641566250 == zetaFunction(16) / pow(pi(), 16));
        CHECK(_43867 / 38979295480125 == zetaFunction(18) / pow(pi(), 18));
        CHECK(_174611 / 1531329465290625 == zetaFunction(20) / pow(pi(), 20));
        CHECK(_77683 / 13447856940643124 == zetaFunction(22) / pow(pi(), 22));
        CHECK(_236364091 / BmpReal(201919571963756511232.0) == zetaFunction(24) / pow(pi(), 24));
    }
}
