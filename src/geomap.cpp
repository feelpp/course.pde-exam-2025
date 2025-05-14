// Geometry.cpp
#include <geomap.hpp>
#include <stdexcept>
#include <Eigen/Dense>


GeometricTransform::GeometricTransform()
  : coords_(Eigen::Matrix<double,2,3>::Zero()),
    J_(Eigen::Matrix2d::Zero()),
    detJ_(0.0)
{}

void GeometricTransform::update(const Mesh& mesh, int elemIndex)
{
    if(elemIndex < 0 || elemIndex >= int(mesh.numElements()))
        throw std::out_of_range("GeometricTransform::update: bad elemIndex");

    const auto& tri = mesh.elements().at(elemIndex);
    // load physical vertex coords using Eigen::Vector2d
    for(int i = 0; i < 3; ++i) {
        coords_.col(i) = mesh.nodes().at(tri.v[i]);
    }

    // reference derivatives for P1 are constant
    Eigen::Matrix<double,3,2> dPhi;
    for(int i = 0; i < 3; ++i) {
        double dr, ds;
        basis_.dphi(i, dr, ds);
        dPhi(i,0) = dr;
        dPhi(i,1) = ds;
    }

    // compute constant Jacobian and its determinant
    J_    = coords_ * dPhi;       
    detJ_ = J_.determinant();
}

Eigen::Vector2d GeometricTransform::map(double r, double s) const
{
    Eigen::Vector3d phi;
    for(int i = 0; i < 3; ++i)
        phi(i) = basis_.phi(i, r, s);
    return coords_ * phi;
}

const Eigen::Matrix2d& GeometricTransform::jacobian() const
{
    return J_;
}

double GeometricTransform::detJ() const
{
    return detJ_;
}