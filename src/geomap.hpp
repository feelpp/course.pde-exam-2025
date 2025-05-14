#pragma once

#include <mesh.hpp>
#include <basis.hpp>
#include <Eigen/Dense>

class GeometricTransform {
    public:
        GeometricTransform();
    
        // prepare mapping for element `elemIndex` of `mesh`
        void update(const Mesh& mesh, int elemIndex);
    
        // map (r,s) in reference triangle → (x,y) in physical element
        Eigen::Vector2d map(double r, double s) const;
    
        // constant Jacobian matrix ∂(x,y)/∂(r,s)
        const Eigen::Matrix2d& jacobian() const;
    
        // constant determinant of the Jacobian
        double detJ() const;
    
    private:
        Eigen::Matrix<double,2,3> coords_;    // [ x0 x1 x2 ; y0 y1 y2 ]
        Eigen::Matrix2d             J_;       // constant Jacobian
        double                      detJ_;    // |J_|
        Basis                       basis_;
    };