#include "base.hpp"

NoopKernel::NoopKernel(){}

double NoopKernel::pdf(double x, int z, const StampLists & stamps) { return 1.0;}

void NoopKernel::build(size_t n_trace, size_t nz, vector<double> priors) {
    P = Eigen::MatrixXd(0, 0);
}

Eigen::MatrixXd NoopKernel::get_state() { return P ; }
void NoopKernel::update_state(const Eigen::MatrixXd & P) {}

void NoopKernel::m_step(const StampLists & stamps) {}

std::unique_ptr<KernelBase> kernel_factory(const std::string & name) {
    if (name=="noop") {
        auto ptr = new NoopKernel();
        return std::unique_ptr<KernelBase>(ptr);
    }
    throw std::runtime_error("unknown kernel name");
}

