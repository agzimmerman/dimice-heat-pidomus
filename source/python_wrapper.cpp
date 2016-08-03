#include <pybind11/pybind11.h>
#include <heat_problem.hpp>
namespace py = pybind11;
PYBIND11_PLUGIN(dealii_heat_problem) {
    py::module m("dealii_heat_problem", "pybind11 dealii_heat_problem plugin");
    py::class_<HeatProblem>(m, "HeatProblem")
        .def(py::init<std::string>())
        .def("solve", &HeatProblem::solve);
    return m.ptr();
}