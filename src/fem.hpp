// FEM.hpp
#pragma once

#include <mesh.hpp>
#include <basis.hpp>
#include <geomap.hpp>
#include <functionspace.hpp>
#include <quadrature.hpp>
#include <Eigen/Sparse>
#include <functional>
#include <vector>
#include <cmath>
#include <iostream>

// Dirichlet boundary detection by coordinate (unit‐square or triangle [0,1]²/[0,1]×[0,1])
inline bool isBoundaryNode(const Eigen::Vector2d& p, double tol = 1e-8) {
    return (std::abs(p.x()) < tol || std::abs(p.x() - 1.0) < tol ||
            std::abs(p.y()) < tol || std::abs(p.y() - 1.0) < tol);
}

using SpMat = Eigen::SparseMatrix<double>;
using Vec   = Eigen::VectorXd;

// Linear form: f(x,y)
class LinearForm {
public:
    LinearForm(const FunctionSpace& V,
               std::function<double(double,double)> f)
      : V_(V), f_(f), F_(V.numDofs())
    {
        assemble();
    }

    const Vec& vector() const { return F_; }

private:
    const FunctionSpace& V_;
    std::function<double(double,double)> f_;
    Vec F_;

    void assemble() {
        const Mesh& M = V_.mesh();
        GeometricTransform G;
        // P1 3‐point quadrature on reference triangle
        // use quadrature rule for degree 2
        const QuadratureRule& rule = getRule(2);


        F_.setZero();
        for(int e = 0; e < int(M.numElements()); ++e) {
            G.update(M, e);
            // local load vector
            Eigen::Vector3d Fe = Eigen::Vector3d::Zero();
            for(int q = 0; q < 3; ++q) {
                double r = rule.points[q].x();
                double s = rule.points[q].y();
                double w = rule.weights[q];
                auto p = G.map(r,s);
                double fx = f_(p.x(), p.y());
                for(int i = 0; i < 3; ++i) {
                    double phi_i = Basis().phi(i, r, s);
                    Fe(i) += w * fx * phi_i * G.detJ();
                }
            }
            // scatter to global
            for(int i = 0; i < 3; ++i)
                F_( V_.localToGlobal(e,i) ) += Fe(i);

            
        }
        // sum all terms in F and print
        //std::cout << "F=" << F_.sum() << "\n";
    }
};

// Bilinear form: stiffness for -Δu
class BilinearForm {
public:
    BilinearForm(const FunctionSpace& V)
      : V_(V), K_(V.numDofs(), V.numDofs())
    {
        assemble();
    }

    // Solve K u = F with Dirichlet g on boundary
// inside FEM.hpp, replace solve(...) with:
    Vec solve(const LinearForm& L,
        const std::map<std::string,std::function<double(double,double)>>& bc )
    {
        int n = int(V_.numDofs());
        Vec F = L.vector();

        // 1) Build isDir[] & gval[] from named BCs
        std::vector<bool>       isDir(n,false);
        Vec                     gval = Vec::Zero(n);

        for (auto const& [name, gfunc] : bc) {
        // only consider names containing "Dirichlet"
        if (name.find("Dirichlet") == std::string::npos) continue;

        int tag = V_.mesh().boundaryTag(name);
        if (tag == 0) {
            std::cerr << "Warning: no physical boundary named '" 
                        << name << "'\n";
            continue;
        }

        // mark all nodes on that tag
        for (int i = 0; i < n; ++i) {
            if (V_.mesh().nodeMarker(i) == tag) {
                isDir[i]   = true;
                auto& P    = V_.mesh().nodes()[i];
                gval(i)    = gfunc(P.x(), P.y());
            }
        }
        }

        // 2) Modify RHS: for Dirichlet rows F(i)=gval; for interior subtract K(i,j)*gval(j)
        for (int i = 0; i < n; ++i) {
            if (isDir[i]) {
                F(i) = gval(i);
            } 
        }

        // 3) Enforce Dirichlet in K_: zero rows+cols, then diag=1
        for (int i = 0; i < n; ++i) {
        if (!isDir[i]) continue;

        // zero out row i
        for (int c = 0; c < n; ++c) {
            for (auto it = Eigen::SparseMatrix<double>::InnerIterator(K_, c); it; ++it) {
                if (it.row() == i) it.valueRef() = 0.0;
            }
        }
        // place 1.0 on the diagonal
        K_.coeffRef(i,i) = 1.0;
        }

        K_.prune(0.0);

        // 4) Solve with SparseLU (or your solver of choice)
        Eigen::SparseLU<SpMat> solver;
        solver.analyzePattern(K_);
        solver.factorize(K_);
        //std::cout << "K_=" << K_ << "\n";
        //std::cout << "F=" << F.transpose() << "\n";
        //std::cout << "sol=" << solver.solve(F).transpose() << "\n";
        return solver.solve(F);
    }

private:
    const FunctionSpace& V_;
    SpMat                 K_;

    void assemble() {
        const Mesh& M = V_.mesh();
        GeometricTransform G;

        // use quadrature rule for degree 1
        const QuadratureRule& rule = getRule(1);
        // reference triangle vertices

        // reference derivatives for P1
        Eigen::Matrix<double,3,2> dPhi;
        for(int i=0;i<3;++i){
            double dr, ds;
            Basis().dphi(i, dr, ds);
            dPhi(i,0)=dr; dPhi(i,1)=ds;
        }

        std::vector<Eigen::Triplet<double>> trips;
        for(int e=0; e<int(M.numElements()); ++e) {
            G.update(M,e);
            Eigen::Matrix3d Ke = Eigen::Matrix3d::Zero();
            Eigen::Matrix2d Jinv = G.jacobian().inverse();
            double detJ = G.detJ();
            // element stiffness
            for(int i=0;i<3;++i) 
                for(int j=0;j<3;++j){
                    Eigen::Vector2d gi = Jinv * dPhi.row(i).transpose();
                    Eigen::Vector2d gj = Jinv * dPhi.row(j).transpose();
                    Ke(i,j) = gi.dot(gj) * detJ;
            }
            // scatter
            
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    trips.emplace_back(V_.localToGlobal(e,i), V_.localToGlobal(e,j), Ke(i,j));
        }
        K_.setFromTriplets(trips.begin(), trips.end());
        //std::cout << "K=" << K_ << "\n";
        checkSum();
    }
    // check that sum of rows is 0
    void checkSum() {
        for(int i=0;i<int(K_.outerSize());++i) {
            double sum = 0.0;
            for(SpMat::InnerIterator it(K_, i); it; ++it) {
                sum += it.value();
            }
            //if (std::abs(sum) > 1e-12)
            //    std::cerr << "Row " << i << " sum = " << sum << "\n";
        }
    }
};