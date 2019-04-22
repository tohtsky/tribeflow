#include "plearn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using std::size_t;

PYBIND11_MODULE(tribeflowpp, m) {
    py::class_<OutPutData> (m, "OutPutData")
        .def(py::init<>())
        .def_readwrite("n_topics", &OutPutData::n_topics)
        .def_readwrite("alpha_zh", &OutPutData::alpha_zh)
        .def_readwrite("alpha_zh", &OutPutData::beta_zs)
        .def_readwrite("residency_priors", &OutPutData::residency_priors)
        .def_readwrite("n_iter", &OutPutData::n_iter)
        .def_readwrite("kernel_name", &OutPutData::kernel_name)
        .def_readwrite("burn_in", &OutPutData::burn_in)
        .def_readwrite("Theta_zh", &OutPutData::Theta_zh)
        .def_readwrite("Psi_sz", &OutPutData::Psi_sz)
        .def_readwrite("Dts", &OutPutData::Dts)
        .def_readwrite("Count_zh", &OutPutData::Count_zh)
        .def_readwrite("Count_sz", &OutPutData::Count_sz)
        .def_readwrite("count_h", &OutPutData::count_h)
        .def_readwrite("count_z", &OutPutData::count_z) 
        .def_readwrite("assign", &OutPutData::assign) 
        .def_readwrite("hyper2id", &OutPutData::hyper2id) 
        .def_readwrite("site2id", &OutPutData::site2id) 
        ;

    m.doc() = R"pbdoc(
        C++ Re implementation of Tribeflow
        -----------------------
        .. currentmodule:: tribeflowpp
        .. autosummary::
           :toctree: _generate
           learn
    )pbdoc";

    m.def("learn", &plearn, R"pbdoc(
        pybind parallel learn function
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
