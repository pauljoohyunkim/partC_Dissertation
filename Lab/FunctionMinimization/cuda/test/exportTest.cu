#include "../src/export.hpp"
#include <vector>

int main()
{
    JsonExporter exporter("hello.json");

    std::vector<double> v1 {1.0, 2.0, 3.1};
    std::vector<double> v2 {1.2, 2.3, 3.8};

    exporter << v1;
    exporter << v2;

    return 0;
}
