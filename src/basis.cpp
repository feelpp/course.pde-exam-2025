#include <basis.hpp>
#include <stdexcept>
#include <Eigen/Dense>

//-- raw (r,s) versions -----------------------------------

double Basis::phi(int i, double r, double s) const
{
    switch(i)
    {
        case 0: return 1.0 - r - s;
        case 1: return r;
        case 2: return s;
        default: throw std::out_of_range("Basis::phi: invalid index");
    }
}

void Basis::dphi(int i, double& dr, double& ds) const
{
    // ∂/∂r [1−r−s] = −1, ∂/∂s = −1
    // ∂/∂r [ r ]    =  1, ∂/∂s =  0
    // ∂/∂r [ s ]    =  0, ∂/∂s =  1
    switch(i)
    {
        case 0: dr = -1.0; ds = -1.0; break;
        case 1: dr =  1.0; ds =  0.0; break;
        case 2: dr =  0.0; ds =  1.0; break;
        default: throw std::out_of_range("Basis::dphi: invalid index");
    }
}

//-- Eigen‐vector overloads --------------------------------

// evaluate φ_i at a point p = (r,s)
double Basis::phi(int i, const Eigen::Vector2d& p) const
{
    return phi(i, p.x(), p.y());
}

// evaluate derivative component c (0 for ∂/∂r, 1 for ∂/∂s)
double Basis::dphi(int i, int c, const Eigen::Vector2d& p) const
{
    double dr, ds;
    dphi(i, dr, ds);
    return (c == 0 ? dr : ds);
}