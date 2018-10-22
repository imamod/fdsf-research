#include "Common.h"

const inline BmpReal fdJmhalf2X30() {
    return 1790.3529161501788;
}

TEST_CASE("jRestore") {
    setPreciseOutput();
    const BmpReal x = 30.0;
    BmpReal delta = _1 - fdJmhalf2X30() / (2*x*x); // 0.005359...
    BmpReal jRestored = log(x) - 6*x*x*delta / (pi()*pi());
    std::cout << x << ") : " << jRestored << std::endl;
}
