// error.cpp
#include <error.hpp>
#include <mesh.hpp>
#include <geomap.hpp>
#include <basis.hpp>
#include <quadrature.hpp>
#include <Eigen/Dense>

double computeL2Error(
    const FunctionSpace& V,
    const Eigen::VectorXd& sol,
    std::function<double(double,double)> exact)
{
    const Mesh& M = V.mesh();
    GeometricTransform G;
    double err2 = 0.0;

    // use order 4
    const QuadratureRule& qr = getRule(4);
    const std::vector<double>& qw = qr.weights;

    // Loop over elements
    for (int e = 0; e < int(M.numElements()); ++e) {
        G.update(M, e);
        double detJ = G.detJ();

        // Quadrature on this element
        for (int q = 0; q < 3; ++q) {
            double r = qr.points[q][0];
            double s = qr.points[q][1];

            Eigen::Vector2d p = G.map(r, s);

            // Evaluate FE solution at (r,s)
            double uh = 0.0;
            for (int i = 0; i < 3; ++i) {
                int gdof = V.localToGlobal(e, i);
                double phi_i = Basis().phi(i, r, s);
                uh += sol(gdof) * phi_i;
            }

            // Exact solution
            double uex = exact(p.x(), p.y());

            double diff = uh - uex;
            err2 += qw[q] * diff * diff * detJ;
        }
    }

    return std::sqrt(err2);
}