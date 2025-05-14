// quadrature.hpp
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>

// A simple holder for points & weights on the reference triangle
struct QuadratureRule {
    std::vector<Eigen::Vector2d> points;
    std::vector<double>          weights;
};

// Returns a rule exact for all polynomials up to total degree `deg`
inline const QuadratureRule& getRule(int deg) {
    // static map from max‐degree → rule
    static const std::map<int, QuadratureRule> rules = {
        { 1, /* 1‐point rule exact deg=1 */
            { { {1./3,1./3} },
              {        0.5 } }  // area of ref triangle is 1/2
        },
        { 2, /* 3‐point rule exact deg=2 */
            { { {1./6,1./6}, {2./3,1./6}, {1./6,2./3} },
              { 1./6,      1./6,        1./6      } }
        },
        { 4, { 
            { 
              Eigen::Vector2d(0.445948490915965, 0.445948490915965),
              Eigen::Vector2d(0.445948490915965, 0.108103018168070),
              Eigen::Vector2d(0.108103018168070, 0.445948490915965),
              Eigen::Vector2d(0.091576213509771, 0.091576213509771),
              Eigen::Vector2d(0.091576213509771, 0.816847572980459),
              Eigen::Vector2d(0.816847572980459, 0.091576213509771)
            },
            {
              0.111690794839005, 
              0.111690794839005, 
              0.111690794839005,
              0.054975871827661,
              0.054975871827661,
              0.054975871827661
            }
          }
        }
        // add more if you like
    };

    auto it = rules.lower_bound(deg);
    if (it == rules.end()) {
        throw std::out_of_range("No quadrature rule for degree " + std::to_string(deg));
    }
    return it->second;
}

// Integrate an arbitrary f(r,s) using the rule for degree `deg`:
inline double integrateTriangle(
    std::function<double(double,double)> f,
    int exactDegree)
{
    const auto& qr = getRule(exactDegree);
    double sum = 0;
    for (size_t q = 0; q < qr.points.size(); ++q) {
        auto c = qr.points[q].array();
        sum += qr.weights[q] * f(c[0], c[1]);
    }
    return sum;
}