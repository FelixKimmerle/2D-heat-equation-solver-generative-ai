#ifndef HEAT_SOLVER_HPP
#define HEAT_SOLVER_HPP

#include <string>
#include <vector>

/**
 * @brief Parameter set for the 2D heat equation solver.
 */
struct HeatParams
{
    int nx = 0;              ///< Number of grid points in x
    int ny = 0;              ///< Number of grid points in y
    double lx = 1.0;         ///< Domain length in x
    double ly = 1.0;         ///< Domain length in y
    double alpha = 1.0;      ///< Diffusion coefficient
    double dt = 1e-3;        ///< Time step size
    double t_final = 0.1;    ///< Final time
    int ic_type = 0;         ///< Initial condition type (0: zero, 1: Gaussian)
    std::string output_file; ///< VTU output file name
};

/**
 * @brief Read solver parameters from a text file.
 *
 * The file format is "key value" pairs, '#' starts a comment.
 * Mandatory keys: nx, ny, lx, ly, alpha, dt, t_final, ic_type, output_file
 *
 * @param filename Path to the parameter file
 * @return HeatParams filled parameter structure
 * @throws std::runtime_error on missing keys or invalid values.
 */
HeatParams read_parameters(const std::string& filename);

/**
 * @brief Solve the 2D heat equation using an explicit finite difference scheme.
 *
 * PDE: u_t = alpha * (u_xx + u_yy) on [0,lx] x [0,ly]
 * with homogeneous Dirichlet boundary conditions (u=0 on boundary).
 *
 * @param params Configuration parameters
 * @return std::vector<double> Flattened solution of size nx*ny (row-major: i + nx*j)
 */
std::vector<double> solve_heat_equation(const HeatParams& params);

/**
 * @brief Write a VTU file containing the 2D scalar field.
 *
 * The mesh is a uniform Cartesian grid of (nx x ny) points with
 * spacing dx = lx/(nx-1), dy = ly/(ny-1).
 *
 * @param filename Path to VTU file
 * @param params Problem parameters (for grid definition)
 * @param u Solution vector of size nx*ny (row-major)
 * @throws std::runtime_error on I/O errors.
 */
void write_vtu(const std::string& filename, const HeatParams& params, const std::vector<double>& u);

#endif // HEAT_SOLVER_HPP
