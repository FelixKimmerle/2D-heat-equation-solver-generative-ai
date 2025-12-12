#include "HeatSolver.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>

void test_read_parameters()
{
    // Create a temporary parameter file
    const char* fname = "test_params.txt";
    {
        std::ofstream out(fname);
        out << "nx 5\nny 7\nlx 1.0\nly 2.0\nalpha 0.1\n"
            << "dt 0.01\nt_final 0.1\nic_type 0\noutput_file test.vtu\n";
    }
    HeatParams p = read_parameters(fname);
    assert(p.nx == 5);
    assert(p.ny == 7);
    assert(p.lx == 1.0);
    assert(p.ly == 2.0);
    assert(p.alpha == 0.1);
    assert(p.dt == 0.01);
    assert(p.t_final == 0.1);
    assert(p.ic_type == 0);
    assert(p.output_file == "test.vtu");
    
    // Clean up temporary file (errors ignored as cleanup is best-effort)
    std::remove(fname);
}

void test_solver_dimensions()
{
    HeatParams p;
    p.nx = 10;
    p.ny = 8;
    p.lx = 1.0;
    p.ly = 1.0;
    p.alpha = 0.01;
    p.dt = 0.0001;
    p.t_final = 0.001;
    p.ic_type = 1;
    p.output_file = "";

    auto u = solve_heat_equation(p);
    assert(static_cast<int>(u.size()) == p.nx * p.ny);
}

int main()
{
    try
    {
        test_read_parameters();
        test_solver_dimensions();
        std::cout << "All tests passed.\n";
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Test failed with exception: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
