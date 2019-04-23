#ifndef KERNEL_BASE_HPP
#define KERNEL_BASE_HPP
#include "../defs.hpp"
#include "../stamp_lists.hpp"
#include <vector> 
#include <string>
#include <Eigen/Eigen>
#include <memory>

struct KernelBase {
    virtual double pdf(double x, int z, const StampLists & stamps) = 0;
    virtual void m_step(const StampLists & stamps) = 0; 
    virtual void build(size_t n_trace, size_t nz, vector<double> priors) = 0;
    virtual DoubleMatrix get_state() = 0;
    virtual void update_state(const DoubleMatrix & P) = 0;
    inline virtual ~KernelBase() {}
};

struct NoopKernel: KernelBase {
    NoopKernel();
    virtual double pdf(double x, int z, const StampLists & stamps) override;
    virtual void m_step(const StampLists & stamps) override; 
    virtual void build(size_t n_trace, size_t nz, vector<double> priors) override;
    virtual DoubleMatrix get_state() override;
    virtual void update_state(const DoubleMatrix & P) override;
    private:
    DoubleMatrix P;
    vector<double> priors;
};

std::unique_ptr<KernelBase> kernel_factory(const std::string & name);

#endif

