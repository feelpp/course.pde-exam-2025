#pragma once
#include <Eigen/Dense>


// P1 Lagrange basis on reference triangle
class Basis 
{
    public:
        // Evaluate shape function i (0,1,2) at (r,s) in ref. tri
        double phi(int i, double r, double s) const;
        // Derivatives in reference coords
        void dphi(int i, double& dr, double& ds) const;
    
        // Overloads using Eigen::Vector2d
        double phi(int i, const Eigen::Vector2d& p) const;
        // Component c = 0 => ∂/∂r, c = 1 => ∂/∂s
        double dphi(int i, int c, const Eigen::Vector2d& p) const;
};