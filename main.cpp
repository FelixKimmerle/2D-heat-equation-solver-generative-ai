#include "HeatSolver.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " params.txt\n";
        return 1;
    }

    try
    {
        HeatParams params = read_parameters(argv[1]);
        auto u = solve_heat_equation(params);
        write_vtu(params.output_file, params, u);
        std::cout << "Computation finished. Result written to " << params.output_file << "\n";
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
