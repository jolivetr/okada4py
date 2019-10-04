#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dc3d.h"
#include "disloc3d.h"

//
// TODO: change flagout,flagout2 to int arrays
// 
// Note that flagout2 was previosly disabled (re-enabled since we're using only one station)
//

void disloc3d(double *models, int nmod,
              double *stations, int nstat,
              double mu, double nu,
              double *uout, double *dout, double *sout,
              double *flagout, double *flagout2)
{
    // input arguments
    //  models      [nmod x 10]
    //  stations    [nstat x 3]
    //  mu          [1 x 1]
    //  nu          [1 x 1]
    
    // matrices for return arguments
    //  uout        [nstat x 3]
    //  dout        [nstat x 9]
    //  sout        [nstat x 6]
    //  flagout     [nstat x 1]
    //  flagout2    [nstat x nmod]
    
    double lambda;
    double theta;

    double *model;
    double *stat;
    double *u, *d, *s;
    double *flag;
    double *flag2;
    double iret;
    const bool using_flag2 = true;

    double ux,  uy,  uz;
    double uxx, uxy, uxz;
    double uyx, uyy, uyz;
    double uzx, uzy, uzz;

    double uxt,  uyt,  uzt;
    double uxxt, uxyt, uxzt;
    double uyxt, uyyt, uyzt;
    double uzxt, uzyt, uzzt;

    double alpha;
    double x, y, z;

    double strike;
    double cs,ss;

    double dip;
    double cd,sd;

    double depth;
    double disl1, disl2, disl3;
    double al1, al2;
    double aw1, aw2;

    int i, j;


    lambda = 2*mu*nu/(1 - 2*nu);
    alpha  = (lambda + mu)/(lambda + 2*mu);

    for (i = 0; i < nstat; i++)
    {
        stat  = &stations[3*i];
        flag  = &flagout[i];
        if (using_flag2) { flag2 = &flagout2[nmod*i]; }

        *flag = 0;

        if (stat[2] > 0)
        {
            *flag += 1;
            fprintf(stderr, "Warning: Positive depth given. (station %d)\n", i);
        }
        else
        {
            uxt  = uyt  = uzt  = 0;
            uxxt = uxyt = uxzt = 0;
            uyxt = uyyt = uyzt = 0;
            uzxt = uzyt = uzzt = 0;
            
            for (j = 0; j < nmod; j++)
            {
                model = &models[10*j];

                strike = (model[4] - 90)*(M_PI/180.0);
                cs = cos(strike);
                ss = sin(strike);
                dip = model[3];
                cd = cos(dip * M_PI/180.0);
                sd = sin(dip * M_PI/180.0);
                disl1 = model[7];
                disl2 = model[8];
                disl3 = model[9];
                depth = model[2];
                al1 = model[0]/2;
                al2 = al1;
                aw1 = model[1]/2;
                aw2 = aw1;
                x =  cs * (-model[5] + stat[0]) - ss * (-model[6] + stat[1]);
                y =  ss * (-model[5] + stat[0]) + cs * (-model[6] + stat[1]);
                z = stat[2];

                // Small hack for patches touching the surface
                if (depth - sd*model[1]*0.5 < 0.)
                {
                    // Add 0.01 meters 
                    depth += 0.00001;
                }

                // generate warnings for unphysical models
                if ((model[0] <= 0.) ||
                    (model[1] <= 0.) ||
                    (depth - sd*model[1]*0.5 < 0.))
                {
                    *flag += 1;
                    fprintf(stderr, "%f %f %f \n", model[0], model[1], depth - sd*model[1]*0.5);
                    fprintf(stderr, "Warning: Unphysical model (station %d, subfault %d)\n",i,j);
                }
                else
                {
                    dc3d_(&alpha, &x, &y, &z, &depth, &dip,
                          &al1, &al2, &aw1, &aw2, &disl1, &disl2, &disl3,
                          &ux,  &uy,  &uz, 
                          &uxx, &uyx, &uzx,
                          &uxy, &uyy, &uzy,
                          &uxz, &uyz, &uzz,
                          &iret);

                    // fill flags - the rows contain flags for each station, one row per model/fault
                    *flag += iret;
                    if (using_flag2) { flag2[j] = iret; }

                    //
                    // rotate then add
                    //
                    uxt +=  cs*ux + ss*uy;
                    uyt += -ss*ux + cs*uy;
                    uzt += uz;

                    uxxt += cs*cs*uxx + cs*ss*(uxy + uyx) + ss*ss*uyy; // aaa
                    uxyt += cs*cs*uxy - ss*ss*uyx + cs*ss*(-uxx + uyy); // bbb
                    uxzt += cs*uxz + ss*uyz; // ccc
                    
                    uyxt += -(ss*(cs*uxx + ss*uxy)) + cs*(cs*uyx + ss*uyy); // ddd
                    uyyt += ss*ss*uxx - cs*ss*(uxy + uyx) + cs*cs*uyy; // eee
                    uyzt += -(ss*uxz) + cs*uyz; // fff
                    
                    uzxt += cs*uzx + ss*uzy; // ggg
                    uzyt += -(ss*uzx) + cs*uzy; // hhh
                    uzzt += uzz; // iii
                }
            }

            // jump to appropriate offset
            u = &uout[3*i];
            d = &dout[9*i];
            s = &sout[6*i];

            // assign outputs
            u[0] = uxt;
            u[1] = uyt;
            u[2] = uzt;
            
            d[0] = uxxt;
            d[1] = uxyt;
            d[2] = uxzt;

            d[3] = uyxt;
            d[4] = uyyt;
            d[5] = uyzt;

            d[6] = uzxt;
            d[7] = uzyt;
            d[8] = uzzt;

            // calculate stresses
            theta = d[0] + d[4] + d[8];
            s[0] = lambda*theta + 2*mu*d[0];
            s[1] = mu*(d[1] + d[3]);
            s[2] = mu*(d[2] + d[6]);
            s[3] = lambda*theta + 2*mu*d[4];
            s[4] = mu*(d[5] + d[7]);
            s[5] = lambda*theta + 2*mu*d[8];

        } // loop over models (subfaults)

    } // loop over stations
}
