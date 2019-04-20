#ifndef KERNEL_BASE_HPP
#define KERNEL_BASE_HPP
#include "../stamp_lists.hpp"
#include <Eigen/Eigen>

struct KernelBase {
    virtual double pdf(double x, int z, const StampLists & stamps) = 0;
    virtual void build(size_t n_trace, size_t nz, vector<double> priors) = 0;
    virtual Eigen::MatrixXd get_state() = 0;
    virtual void update_state(const Eigen::MatrixXd & P) = 0;
};

struct NoopKernel: KernelBase {
    NoopKernel();
    virtual double pdf(double x, int z, const StampLists & stamps) override;
    virtual void build(size_t n_trace, size_t nz, vector<double> priors) override;
    virtual Eigen::MatrixXd get_state() override;
    virtual void update_state(const Eigen::MatrixXd & P) override;
    private:
    Eigen::MatrixXd P;
    vector<double> priors;
};

#endif

