#include "Solver.h"
#include "GraphProcessor.h"
#include <string>
#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <json_file> <param1> <param2> <param3> <param4>" << std::endl;
        return 1;
    }

    std::string json_name = argv[1];
    double param1 = std::stod(argv[2]);
    double param2 = std::stod(argv[3]);
    double param3 = std::stod(argv[4]);
    double param4 = std::stod(argv[5]);

    Solver solver;
    solver.hyperparameters = {param1, param2, param3, param4};
    solver.readSceneGraph(json_name, 0.0f);
    solver.solve();

    return 0;
}