#pragma once

#include <functionspace.hpp>
#include <functional>

double computeL2Error(
    const FunctionSpace& V,
    const Eigen::VectorXd& sol,
    std::function<double(double,double)> exact);
