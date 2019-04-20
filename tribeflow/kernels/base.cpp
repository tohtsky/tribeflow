#include "base.hpp"

NoopKernel::NoopKernel(){}

double NoopKernel::pdf(double x, int z, const StampLists & stamps) { return 1.0;}

void NoopKernel::build(size_t n_trace, size_t nz, vector<double> priors) {
    P = Eigen::MatrixXd(0, 0);
}

Eigen::MatrixXd NoopKernel::get_state() { return P ; }
void NoopKernel::update_state(const Eigen::MatrixXd & P) {}

