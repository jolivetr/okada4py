#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dc3d.h"
#include "disloc3d.h"

int main()
{
    double *models;
    int nmod=1;
    double *stations;
    int nstat=3;
    double mu=3e10;
    double nu=0.25;
    double *uout;
    double *dout;
    double *sout;
    double *flagout;
    double *flagout2;

    // Input
    models = (double *)calloc(10*nmod, sizeof(double));
    stations = (double *)calloc(3*nstat, sizeof(double));

    // Output
    uout = (double *)calloc(3*nstat, sizeof(double));
    dout = (double *)calloc(9*nstat, sizeof(double));
    sout = (double *)calloc(6*nstat, sizeof(double));

    // Flags
    flagout = (double *)calloc(nstat, sizeof(double));
    flagout2 = (double *)calloc(nstat*nmod, sizeof(double));

    // Call disloc3d
    disloc3d(models, nmod, stations, nstat, mu, nu, uout, dout, sout, flagout, flagout2);

    // All done
    return 0;
}

