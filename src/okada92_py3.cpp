#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dc3d.h"
#include "disloc3d.h"
#include "Python.h"
#include "arrayobject.h"

static PyObject *py_Okada(PyObject *self, PyObject *args)
{
    // Arguments for disloc3d
    double *models;         // Dislocations [nmod x 10]: [length, width, depth, dip, strike, xc, yc, ss, ds, ts]
    int nmod;               // Number of dislocations
    double *stations;       // Stations [nstat x 3]: [xs, ys, zs]
    int nstat;              // Number of stations
    double mu;              // Shear modulus    (Used for stress computation)
    double nu;              // Poisson's ration
    double *uout;           // Displacements [nstat x 3]: [Ux, Uy, Uz]
    double *dout;           // Strains [nstat x 9]: [Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz]
    double *sout;           // Stresses [nstat x 6]: [Sxx, Sxy, Sxz, Syy, Syz, Szz]
    double *flagout;        // A flag per station [nstat]
    double *flagout2;       // Another flag per station and per dislocation [nstat x nmod]

    // Equivalent in python
    PyArrayObject *Py_xs, *Py_ys, *Py_zs;   // Station positions
    PyArrayObject *Py_xc, *Py_yc, *Py_depth;// Dislocation position
    PyArrayObject *Py_length, *Py_width;    // Dislocation size
    PyArrayObject *Py_dip, *Py_strike;      // Dislocation orientation
    PyArrayObject *Py_ss, *Py_ds, *Py_ts;   // Strike-, Dip- and Tensile slip 

    // import python arguments and check
    import_array1(NULL);
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!O!O!O!O!O!dd",
                &PyArray_Type, &Py_xs, 
                &PyArray_Type, &Py_ys, 
                &PyArray_Type, &Py_zs, 
                &PyArray_Type, &Py_xc, 
                &PyArray_Type, &Py_yc, 
                &PyArray_Type, &Py_depth,
                &PyArray_Type, &Py_length, 
                &PyArray_Type, &Py_width,
                &PyArray_Type, &Py_dip, 
                &PyArray_Type, &Py_strike, 
                &PyArray_Type, &Py_ss, 
                &PyArray_Type, &Py_ds, 
                &PyArray_Type, &Py_ts, 
                &mu, &nu)) return NULL;

    // Get problem size
    nmod = (int)Py_xc->dimensions[0];
    nstat = (int)Py_xs->dimensions[0];
    /*
    fprintf(stderr, "Problem size: \n");
    fprintf(stderr, "       %d dislocations \n", nmod);
    fprintf(stderr, "       %d stations \n", nstat);
    */

    // Allocate python objects for output 
    PyArrayObject *Py_u;
    PyArrayObject *Py_d;
    PyArrayObject *Py_s;
    npy_intp dims1 = (npy_intp) nstat*3;
    Py_u = (PyArrayObject *) PyArray_ZEROS(1, &dims1, NPY_DOUBLE, 0);  // displacements
    npy_intp dims2 = (npy_intp) nstat*9;
    Py_d = (PyArrayObject *) PyArray_ZEROS(1, &dims2, NPY_DOUBLE, 0);  // strains
    npy_intp dims3 = (npy_intp) nstat*6;
    Py_s = (PyArrayObject *) PyArray_ZEROS(1, &dims3, NPY_DOUBLE, 0);  // stress

    // Link c arrays to the data part of the pythons
    uout = (double *) Py_u->data;
    dout = (double *) Py_d->data;
    sout = (double *) Py_s->data;

    // Allocate python objects for flags
    PyArrayObject *Py_flag;
    PyArrayObject *Py_flag2;
    npy_intp dims4 = (npy_intp) nstat;
    Py_flag = (PyArrayObject *) PyArray_ZEROS(1, &dims4, NPY_DOUBLE, 0); // A flag
    npy_intp dims5 = (npy_intp) nstat*nmod;
    Py_flag2 = (PyArrayObject *) PyArray_ZEROS(1, &dims5, NPY_DOUBLE, 0); // Another flag

    // link the c flags to the python flags
    flagout = (double *) Py_flag->data;
    flagout2 = (double *) Py_flag2->data;

    // Allocate inputs for disloc3d
    models = (double *)calloc(nmod*10, sizeof(double));
    stations = (double *)calloc(nstat*3, sizeof(double));

    // Fill in models
    for (size_t i=0; i<nmod; ++i)
    {
       models[10*i] = ((double *) Py_length->data)[i];
       models[10*i+1] = ((double *) Py_width->data)[i];
       models[10*i+2] = ((double *) Py_depth->data)[i];
       models[10*i+3] = ((double *) Py_dip->data)[i];
       models[10*i+4] = ((double *) Py_strike->data)[i];
       models[10*i+5] = ((double *) Py_xc->data)[i];
       models[10*i+6] = ((double *) Py_yc->data)[i];
       models[10*i+7] = ((double *) Py_ss->data)[i];
       models[10*i+8] = ((double *) Py_ds->data)[i];
       models[10*i+9] = ((double *) Py_ts->data)[i];
       /*
       fprintf(stderr, "Models %d \n", i);
       fprintf(stderr, "    [ %f %f %f %f %f %f %f %f %f %f ] \n", models[10*i], models[10*i+1],
               models[10*i+2], models[10*i+3], models[10*i+4], models[10*i+5], models[10*i+6], models[10*i+7], 
               models[10*i+8], models[10*i+9]);
        */
    }

    // Fill in stations
    for (size_t i=0; i<nstat; ++i)
    {
        stations[3*i] = ((double *) Py_xs->data)[i];
        stations[3*i+1] = ((double *) Py_ys->data)[i];
        stations[3*i+2] = ((double *) Py_zs->data)[i];
    }

    // Call disloc3d
    disloc3d(models, nmod, stations, nstat, mu, nu, uout, dout, sout, flagout, flagout2);

    // Build the tuple for return
    PyObject *results = PyTuple_New(5);
    PyTuple_SetItem(results, 0, PyArray_Return(Py_u));
    PyTuple_SetItem(results, 1, PyArray_Return(Py_d));
    PyTuple_SetItem(results, 2, PyArray_Return(Py_s));
    PyTuple_SetItem(results, 3, PyArray_Return(Py_flag));
    PyTuple_SetItem(results, 4, PyArray_Return(Py_flag2));

    // Free everything
    free(models);
    free(stations);

    // All done
    return results;
}

PyDoc_STRVAR(
    pyOkada_doc,
    "Okada 92 function for computation of displacements, strain and stresses in 3D\n"
    "\n"
    "       Arguments: (xs, ys, zs, xc, yc, depth, length, width, dip, strike, StrikeSlip, DipSlip, Tensile, mu, nu)\n"
    "                   All arguments are arrays, except mu (Shear modulus) and nu (Poisson's ration).\n"
    "                   xs : Stations x coordinates dims=[nstat] \n"
    "                   ys : Stations y coordinates dims=[nstat] \n"
    "                   zs : Stations z coordinates dims=[nstat] \n"
    "                   xc : Dislocation Center x coordinates dims=[nmods] \n"
    "                   yc : Dislocation Center y coordinates dims=[nmods]\n"
    "                depth : Dislocation Center Depth dims=[nmods] \n"
    "               length : Dislocation Length dims=[nmods] \n"
    "                width : Dislocation Width  dims=[nmods] \n"
    "                  dip : Dislocation Dip dims=[nmods] \n"
    "               strike : Dislocation Strike dims=[nmods] \n"
    "           StrikeSlip : Strike-Slip displacement dims=[nmods] \n"
    "              DipSlip : Dip-Slip displacement dims=[nmods] \n"
    "              Tensile : Tensile displacement dims=[nmods] \n"
    "                   Mu : Shear Modulus dims=[1] \n"
    "                   Nu : Poisson's ratio dims=[1] \n"
    "\n"
    "       Returns: \n"
    "                    u : Displacements dims=[nstat*3] [Ux, Uy, Uz] \n"
    "                    d : Strain dims=[nstat*9] [Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz] \n"
    "                    s : Stress dims=[nstat*6] [Sxx, Sxy, Sxz, Syy, Syz, Szz] \n"
    "                 flag : A flag \n"
    "                flag2 : Another flag \n"
    );

static PyMethodDef _okada92_methods[]={
            {"okada92",py_Okada,METH_VARARGS, pyOkada_doc},
            {NULL,NULL}
};

static struct PyModuleDef _okada92module = {
    PyModuleDef_HEAD_INIT,
    "_okada92",
    NULL, 
    -1, 
    _okada92_methods
};

PyMODINIT_FUNC PyInit__okada92(void){
    return PyModule_Create(&_okada92module);
}

} /* Closing brace for extern "C" */

