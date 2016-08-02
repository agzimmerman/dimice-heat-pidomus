#include <pybind11/pybind11.h>
namespace py = pybind11;
PYBIND11_PLUGIN(heat_problem) {
    py::module m("heat_problem", "pybind11 heat_problem plugin");
    py::class_<HeatProblem>(m, "HeatProblem")
        .def(py::init<>())
        .def("solve", &Pet::solve, py::arg("parameter_file_path")="");
    return m.ptr();
}