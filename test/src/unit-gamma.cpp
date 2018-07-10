#include "Common.h"
#include "Gamma.h"

TEST_CASE("gamma") {
    SECTION("integer") {
        CHECK(factorial(0) == 1);
        CHECK(factorial(1) == 1);
        CHECK(factorial(3) == 6);
        CHECK(factorial(4) == 24);
        CHECK(factorial(5) == 120);
    }
    SECTION("noninteger") {
        const BmpReal PI = pi();
        CHECK(factorial(-1.5) == -2 * sqrt(PI));
        CHECK(factorial(-0.5) == sqrt(PI));
        CHECK(factorial(0.5) == sqrt(PI) / 2);
        CHECK(factorial(1.5) == 3 * sqrt(PI) / 4);
        CHECK(factorial(2.5) == 15 * sqrt(PI) / 8);
        CHECK(factorial(3.5) == 105 * sqrt(PI) / 16);
    }
}
