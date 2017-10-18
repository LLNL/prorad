
#include <Python.h>
#include <numpy/arrayobject.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "pmover.h"

static char module_docstring[] =
	"This module provides an interface for moving protons in parallel using C.";

static char pmover_docstring[] = "Calculate the trajectories of a set of protons through a plasma";

static PyObject *move_protons(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
	{"pmover", move_protons, METH_VARARGS, pmover_docstring},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_pmover(void) {
	PyObject *m = Py_InitModule3("_pmover", module_methods, module_docstring);
	if (m == NULL)
		return;

	/* Load 'numpy' functionality */
	import_array();
}

static PyObject *move_protons(PyObject *self, PyObject *args) {
	double qm, ds, dx, dy, dz, gx, gy, gz;
	int NP, ns, ngridx, ngridy, ngridz, ntraces, cyl_coords;
	PyObject *x_obj, *y_obj, *z_obj, *vx_obj, *vy_obj, *vz_obj, *prop_dir_obj, *gridvals_obj, *traces_obj;
	
	/* Parse the input tuple */
	if (!PyArg_ParseTuple(args, "ddddddddiiiiiiiOOOOOOOOO", &qm, &ds, 
					&dx, &dy, &dz, &gx, &gy, &gz,
					&NP, &ns, &ngridx, &ngridy, &ngridz, &ntraces, &cyl_coords,
					&x_obj, &y_obj, &z_obj, 
					&vx_obj, &vy_obj, &vz_obj, &prop_dir_obj, &gridvals_obj, &traces_obj))
		return NULL;
	
	/* Interpret the input objects as numpy arrays. */
	PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *z_array = PyArray_FROM_OTF(z_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *vx_array = PyArray_FROM_OTF(vx_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *vy_array = PyArray_FROM_OTF(vy_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *vz_array = PyArray_FROM_OTF(vz_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *prop_dir_array = PyArray_FROM_OTF(prop_dir_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *gridvals_array = PyArray_FROM_OTF(gridvals_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *traces_array = PyArray_FROM_OTF(traces_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	
	/* If that didn't work, throw an exception. */
	if (x_array == NULL || y_array == NULL || z_array == NULL || 
	    vx_array == NULL || vy_array == NULL || vz_array == NULL ||
	    gridvals_array == NULL || traces_array == NULL || prop_dir_array == NULL) {
		Py_XDECREF(x_array);
		Py_XDECREF(y_array);
		Py_XDECREF(z_array);
		Py_XDECREF(vx_array);
		Py_XDECREF(vy_array);
		Py_XDECREF(vz_array);
		Py_XDECREF(prop_dir_array);
		Py_XDECREF(gridvals_array);
		Py_XDECREF(traces_array);
		return NULL;
	}
	
	/* Get pointers to the data as C-types. */
	double *x  = (double*)PyArray_DATA(x_array);
	double *y  = (double*)PyArray_DATA(y_array);
	double *z  = (double*)PyArray_DATA(z_array);
	double *vx = (double*)PyArray_DATA(vx_array);
	double *vy = (double*)PyArray_DATA(vy_array);
	double *vz = (double*)PyArray_DATA(vz_array);
	double *prop_dir = (double*)PyArray_DATA(prop_dir_array);
	double *gridvals_flat = (double*) PyArray_DATA(gridvals_array);
	double *traces_flat = (double*) PyArray_DATA(traces_array);
	
	double grid_spacing[] = {dx, dy, dz};
	double grid_offset[] = {gx, gy, gz};
    
	/* Call the external C function to move the protons. */
	pmover(NP, ns, qm, ds, x, y, z, vx, vy, vz, prop_dir, ngridx, ngridy, ngridz, 
            gridvals_flat, grid_offset, grid_spacing, ntraces, traces_flat, cyl_coords);


	/* Clean up. */
	Py_DECREF(x_array);
	Py_DECREF(y_array);
	Py_DECREF(z_array);
	Py_DECREF(vx_array);
	Py_DECREF(vy_array);
	Py_DECREF(vz_array);
	Py_DECREF(prop_dir_array);
	Py_DECREF(gridvals_array);
	Py_DECREF(traces_array);

	return Py_BuildValue("");
}
