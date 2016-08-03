#include <pybind11/pybind11.h>
#include <heat_problem.hpp>
namespace py = pybind11;
PYBIND11_PLUGIN(dealii_heat_problem) {
    py::module m("dealii_heat_problem", "pybind11 dealii_heat_problem plugin");
    py::class_<HeatProblem>(m, "HeatProblem")
        .def(py::init<const std::string, const int>(),
             py::arg("parameter_file_name")="/../inputs/sphere-cylinder.prm",
             py::arg("mpi_size")=1
        )
        .def("solve", &HeatProblem::solve);
    return m.ptr();
}