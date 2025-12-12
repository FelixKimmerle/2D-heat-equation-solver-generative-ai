
#include "HeatSolver.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cmath>

static void check_required(const HeatParams& p)
{
    if (p.nx <= 1 || p.ny <= 1)
        throw std::runtime_error("nx and ny must be > 1");
    if (p.lx <= 0.0 || p.ly <= 0.0)
        throw std::runtime_error("lx and ly must be > 0");
    if (p.alpha <= 0.0)
        throw std::runtime_error("alpha must be > 0");
    if (p.dt <= 0.0 || p.t_final <= 0.0)
        throw std::runtime_error("dt and t_final must be > 0");
    if (p.output_file.empty())
        throw std::runtime_error("output_file must not be empty");
}

HeatParams read_parameters(const std::string& filename)
{
    HeatParams p;
    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("Cannot open parameter file: " + filename);

    std::string line;
    while (std::getline(in, line))
    {
        // Remove comments
        auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
            line = line.substr(0, comment_pos);
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key))
            continue;
        if (key == "nx")
            iss >> p.nx;
        else if (key == "ny")
            iss >> p.ny;
        else if (key == "lx")
            iss >> p.lx;
        else if (key == "ly")
            iss >> p.ly;
        else if (key == "alpha")
            iss >> p.alpha;
        else if (key == "dt")
            iss >> p.dt;
        else if (key == "t_final")
            iss >> p.t_final;
        else if (key == "ic_type")
            iss >> p.ic_type;
        else if (key == "output_file")
            iss >> p.output_file;
        // unknown keys are ignored
    }

    check_required(p);
    return p;
}

std::vector<double> solve_heat_equation(const HeatParams& params)
{
    const int nx = params.nx;
    const int ny = params.ny;
    const double lx = params.lx;
    const double ly = params.ly;
    const double alpha = params.alpha;
    const double dt = params.dt;
    const double t_final = params.t_final;

    const double dx = lx / (nx - 1);
    const double dy = ly / (ny - 1);

    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    // Stability condition for explicit scheme (approximate)
    double coeff = alpha * dt * (1.0 / dx2 + 1.0 / dy2);
    if (coeff >= 0.5)
    {
        std::cerr << "Warning: explicit scheme may be unstable; "
                  << "alpha*dt*(1/dx^2+1/dy^2) = " << coeff << " >= 0.5\n";
    }

    std::vector<double> u(nx * ny, 0.0);
    std::vector<double> u_new(nx * ny, 0.0);

    // Initial condition
    if (params.ic_type == 1)
    {
        // Gaussian in the center
        double cx = 0.5 * lx;
        double cy = 0.5 * ly;
        double sigma2 = 0.01; // variance
        for (int j = 0; j < ny; ++j)
        {
            double y = j * dy;
            for (int i = 0; i < nx; ++i)
            {
                double x = i * dx;
                double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                u[i + nx * j] = std::exp(-r2 / (2.0 * sigma2));
            }
        }
    }
    // else ic_type==0: already zero

    double t = 0.0;
    while (t < t_final)
    {
        // ensure last step hits (or slightly exceeds) t_final
        double dt_step = std::min(dt, t_final - t);

        // interior points update (Dirichlet u=0 at boundary)
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                int idx = i + nx * j;
                double uij = u[idx];
                double uxx = (u[idx - 1] - 2.0 * uij + u[idx + 1]) / dx2;
                double uyy = (u[idx - nx] - 2.0 * uij + u[idx + nx]) / dy2;
                u_new[idx] = uij + dt_step * alpha * (uxx + uyy);
            }
        }

        // boundaries: homogeneous Dirichlet
        for (int i = 0; i < nx; ++i)
        {
            u_new[i] = 0.0;                 // bottom
            u_new[i + nx * (ny - 1)] = 0.0; // top
        }
        for (int j = 0; j < ny; ++j)
        {
            u_new[0 + nx * j] = 0.0;        // left
            u_new[(nx - 1) + nx * j] = 0.0; // right
        }

        u.swap(u_new);
        t += dt_step;
    }

    return u;
}

void write_vtu(const std::string& filename, const HeatParams& params, const std::vector<double>& u)
{
    const int nx = params.nx;
    const int ny = params.ny;
    const double lx = params.lx;
    const double ly = params.ly;

    const double dx = lx / (nx - 1);
    const double dy = ly / (ny - 1);

    const int num_points = nx * ny;
    const int num_cells = (nx - 1) * (ny - 1);

    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("Cannot open output file: " + filename);

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">\n";

    // Points
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int j = 0; j < ny; ++j)
    {
        double y = j * dy;
        for (int i = 0; i < nx; ++i)
        {
            double x = i * dx;
            out << "          " << x << " " << y << " " << 0.0 << "\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // Cells (quads)
    out << "      <Cells>\n";

    // connectivity
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int j = 0; j < ny - 1; ++j)
    {
        for (int i = 0; i < nx - 1; ++i)
        {
            int p0 = i + nx * j;
            int p1 = (i + 1) + nx * j;
            int p2 = (i + 1) + nx * (j + 1);
            int p3 = i + nx * (j + 1);
            out << "          " << p0 << " " << p1 << " " << p2 << " " << p3 << "\n";
        }
    }
    out << "        </DataArray>\n";

    // offsets
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (int c = 0; c < num_cells; ++c)
    {
        offset += 4;
        out << "          " << offset << "\n";
    }
    out << "        </DataArray>\n";

    // types (VTK_QUAD = 9)
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int c = 0; c < num_cells; ++c)
    {
        out << "          9\n";
    }
    out << "        </DataArray>\n";

    out << "      </Cells>\n";

    // Point data (scalar field u)
    if (static_cast<int>(u.size()) != num_points)
        throw std::runtime_error("Size of solution vector does not match grid.");

    out << "      <PointData Scalars=\"temperature\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"temperature\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int i = 0; i < num_points; ++i)
    {
        out << "          " << u[i] << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </PointData>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}
