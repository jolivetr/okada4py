/* disloc3d_mod2.f -- translated by f2c (version 20050501).
 * Link this file with libm (add linker flag -lm)
 */

#include "dc3d.h"
// #include "f2c.h"
typedef double doublereal;
typedef long int integer;
//#define abs(x) ((x) >= 0 ? (x) : -(x))

#include <stdio.h>
#include <math.h>
#include <cmath>
using std::abs;

/* Common Block Declarations */

union {
    struct {
        doublereal dummy[5], sd, cd, dummy2[5];
    } _1;
    struct {
        doublereal alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, 
                s2d, c2d;
    } _2;
} c0_;

#define c0_1 (c0_._1)
#define c0_2 (c0_._2)

union {
    struct {
        doublereal xi2, et2, q2, r__, dummy3[20];
    } _1;
    struct {
        doublereal xi2, et2, q2, r__, r2, r3, r5, y, d__, tt, alx, ale, x11, 
                y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz;
    } _2;
} c2_;

#define c2_1 (c2_._1)
#define c2_2 (c2_._2)

/* Subroutine */ int dc3d_(doublereal *alpha, doublereal *x, doublereal *y, 
        doublereal *z__, doublereal *depth, doublereal *dip, doublereal *al1, 
        doublereal *al2, doublereal *aw1, doublereal *aw2, doublereal *disl1, 
        doublereal *disl2, doublereal *disl3, doublereal *ux, doublereal *uy, 
        doublereal *uz, doublereal *uxx, doublereal *uyx, doublereal *uzx, 
        doublereal *uxy, doublereal *uyy, doublereal *uzy, doublereal *uxz, 
        doublereal *uyz, doublereal *uzz, doublereal *iret)
{
    /* Initialized data */

    static doublereal f0 = 0.;

    /* Builtin functions */
    //integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;
    static doublereal p, q, u[12];
    extern /* Subroutine */ int ua_(doublereal *, doublereal *, doublereal *, 
            doublereal *, doublereal *, doublereal *, doublereal *), ub_(
            doublereal *, doublereal *, doublereal *, doublereal *, 
            doublereal *, doublereal *, doublereal *);
    static doublereal du[12], et;
    extern /* Subroutine */ int uc_(doublereal *, doublereal *, doublereal *, 
            doublereal *, doublereal *, doublereal *, doublereal *, 
            doublereal *);
    static doublereal xi, zz, dd1, dd2, dd3, dua[12], dub[12], duc[12];
    static integer jet, jxi;
    static doublereal ddip;
    extern /* Subroutine */ int dccon0_(doublereal *, doublereal *), dccon2_(
            doublereal *, doublereal *, doublereal *, doublereal *, 
            doublereal *);
    static doublereal aalpha;

    /* Fortran I/O blocks */
    //static cilist io___2 = { 0, 6, 0, "(' ** POSITIVE Z WAS GIVEN IN SUB-DC3"
    //        "D')", 0 };


/*                                                                       04670000 */
/* ********************************************************************   04680000 */
/* *****                                                          *****   04690000 */
/* *****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   04700000 */
/* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   04710000 */
/* *****                         CODED BY  Y.OKADA ... SEP 1991   *****   04720002 */
/* *****                         REVISED   Y.OKADA ... NOV 1991   *****   04730002 */
/* *****                                                          *****   04740000 */
/* ********************************************************************   04750000 */
/*                                                                       04760000 */
/* ***** INPUT                                                            04770000 */
/* *****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           04780000 */
/* *****   X,Y,Z : COORDINATE OF OBSERVING POINT                          04790000 */
/* *****   DEPTH : SOURCE DEPTH                                           04800000 */
/* *****   DIP   : DIP-ANGLE (DEGREE)                                     04810000 */
/* *****   AL1,AL2   : FAULT LENGTH (-STRIKE,+STRIKE)                     04820000 */
/* *****   AW1,AW2   : FAULT WIDTH  ( DOWNDIP, UPDIP)                     04830000 */
/* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              04840000 */
/*                                                                       04850000 */
/* ***** OUTPUT                                                           04860000 */
/* *****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               04870000 */
/* *****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /             04880000 */
/* *****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )04890000 */
/* *****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     04900000 */
/* *****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  04910002 */
/*                                                                       04920000 */
/* -----                                                                  04970000 */
    if (*z__ > 0.f) {
        fprintf(stderr, "(' ** POSITIVE Z WAS GIVEN IN SUB-DC3D')");
        //s_wsfe(&io___2);
        //e_wsfe();
    }
    for (i__ = 1; i__ <= 12; ++i__) {
        u[i__ - 1] = f0;
        dua[i__ - 1] = f0;
        dub[i__ - 1] = f0;
        duc[i__ - 1] = f0;
/* L111: */
    }
    aalpha = *alpha;
    ddip = *dip;
    dccon0_(&aalpha, &ddip);
/* ======================================                                 05080000 */
/* =====  REAL-SOURCE CONTRIBUTION  =====                                 05090000 */
/* ======================================                                 05100000 */
    d__ = *depth + *z__;
    p = *y * c0_1.cd + d__ * c0_1.sd;
    q = *y * c0_1.sd - d__ * c0_1.cd;
    jxi = 0;
    jet = 0;
    if ((*x + *al1) * (*x - *al2) <= 0.f) {
        jxi = 1;
    }
    if ((p + *aw1) * (p - *aw2) <= 0.f) {
        jet = 1;
    }
    dd1 = *disl1;
    dd2 = *disl2;
    dd3 = *disl3;
/* -----                                                                  05210000 */
    for (k = 1; k <= 2; ++k) {
        if (k == 1) {
            et = p + *aw1;
        }
        if (k == 2) {
            et = p - *aw2;
        }
        for (j = 1; j <= 2; ++j) {
            if (j == 1) {
                xi = *x + *al1;
            }
            if (j == 2) {
                xi = *x - *al2;
            }
            dccon2_(&xi, &et, &q, &c0_1.sd, &c0_1.cd);
            if (jxi == 1 && q == f0 && et == f0) {
                goto L99;
            }
            if (jet == 1 && q == f0 && xi == f0) {
                goto L99;
            }
            ua_(&xi, &et, &q, &dd1, &dd2, &dd3, dua);
/* -----                                                                  05320000 */
            for (i__ = 1; i__ <= 10; i__ += 3) {
                du[i__ - 1] = -dua[i__ - 1];
                du[i__] = -dua[i__] * c0_1.cd + dua[i__ + 1] * c0_1.sd;
                du[i__ + 1] = -dua[i__] * c0_1.sd - dua[i__ + 1] * c0_1.cd;
                if (i__ < 10) {
                    goto L220;
                }
                du[i__ - 1] = -du[i__ - 1];
                du[i__] = -du[i__];
                du[i__ + 1] = -du[i__ + 1];
L220:
                ;
            }
            for (i__ = 1; i__ <= 12; ++i__) {
                if (j + k != 3) {
                    u[i__ - 1] += du[i__ - 1];
                }
                if (j + k == 3) {
                    u[i__ - 1] -= du[i__ - 1];
                }
/* L221: */
            }
/* -----                                                                  05460000 */
/* L222: */
        }
/* L223: */
    }
/* =======================================                                05490000 */
/* =====  IMAGE-SOURCE CONTRIBUTION  =====                                05500000 */
/* =======================================                                05510000 */
    zz = *z__;
    d__ = *depth - *z__;
    p = *y * c0_1.cd + d__ * c0_1.sd;
    q = *y * c0_1.sd - d__ * c0_1.cd;
    jet = 0;
    if ((p + *aw1) * (p - *aw2) <= 0.f) {
        jet = 1;
    }
/* -----                                                                  05580000 */
    for (k = 1; k <= 2; ++k) {
        if (k == 1) {
            et = p + *aw1;
        }
        if (k == 2) {
            et = p - *aw2;
        }
        for (j = 1; j <= 2; ++j) {
            if (j == 1) {
                xi = *x + *al1;
            }
            if (j == 2) {
                xi = *x - *al2;
            }
            dccon2_(&xi, &et, &q, &c0_1.sd, &c0_1.cd);
            ua_(&xi, &et, &q, &dd1, &dd2, &dd3, dua);
            ub_(&xi, &et, &q, &dd1, &dd2, &dd3, dub);
            uc_(&xi, &et, &q, &zz, &dd1, &dd2, &dd3, duc);
/* -----                                                                  05690000 */
            for (i__ = 1; i__ <= 10; i__ += 3) {
                du[i__ - 1] = dua[i__ - 1] + dub[i__ - 1] + *z__ * duc[i__ - 
                        1];
                du[i__] = (dua[i__] + dub[i__] + *z__ * duc[i__]) * c0_1.cd - 
                        (dua[i__ + 1] + dub[i__ + 1] + *z__ * duc[i__ + 1]) * 
                        c0_1.sd;
                du[i__ + 1] = (dua[i__] + dub[i__] - *z__ * duc[i__]) * 
                        c0_1.sd + (dua[i__ + 1] + dub[i__ + 1] - *z__ * duc[
                        i__ + 1]) * c0_1.cd;
                if (i__ < 10) {
                    goto L330;
                }
                du[9] += duc[0];
                du[10] = du[10] + duc[1] * c0_1.cd - duc[2] * c0_1.sd;
                du[11] = du[11] - duc[1] * c0_1.sd - duc[2] * c0_1.cd;
L330:
                ;
            }
            for (i__ = 1; i__ <= 12; ++i__) {
                if (j + k != 3) {
                    u[i__ - 1] += du[i__ - 1];
                }
                if (j + k == 3) {
                    u[i__ - 1] -= du[i__ - 1];
                }
/* L331: */
            }
/* -----                                                                  05850000 */
/* L333: */
        }
/* L334: */
    }
/* =====                                                                  05880000 */
    *ux = u[0];
    *uy = u[1];
    *uz = u[2];
    *uxx = u[3];
    *uyx = u[4];
    *uzx = u[5];
    *uxy = u[6];
    *uyy = u[7];
    *uzy = u[8];
    *uxz = u[9];
    *uyz = u[10];
    *uzz = u[11];
    *iret = 0.;
    return 0;
/* =======================================                                06030000 */
/* =====  IN CASE OF SINGULAR (R=0)  =====                                06040000 */
/* =======================================                                06050000 */
L99:
    *ux = f0;
    *uy = f0;
    *uz = f0;
    *uxx = f0;
    *uyx = f0;
    *uzx = f0;
    *uxy = f0;
    *uyy = f0;
    *uzy = f0;
    *uxz = f0;
    *uyz = f0;
    *uzz = f0;
    *iret = 1.;
    return 0;
} /* dc3d_ */

/* Subroutine */ int dccon0_(doublereal *alpha, doublereal *dip)
{
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal pi2 = 6.283185307179586;
    static doublereal eps = 1e-6;

    /* Local variables */
    static doublereal p18;

/*                                                                       09320000 */
/* *******************************************************************    09330000 */
/* *****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****    09340000 */
/* *******************************************************************    09350000 */
/*                                                                       09360000 */
/* ***** INPUT                                                            09370000 */
/* *****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           09380000 */
/* *****   DIP   : DIP-ANGLE (DEGREE)                                     09390000 */
/* ### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO   09400000 */
/*                                                                       09410000 */
/* everything in CO is used in this routine */
/* -----                                                                  09450000 */
    c0_2.alp1 = (f1 - *alpha) / f2;
    c0_2.alp2 = *alpha / f2;
    c0_2.alp3 = (f1 - *alpha) / *alpha;
    c0_2.alp4 = f1 - *alpha;
    c0_2.alp5 = *alpha;
/* -----                                                                  09510000 */
    p18 = pi2 / 360.;
    c0_2.sd = sin(*dip * p18);
    c0_2.cd = cos(*dip * p18);
    if (abs(c0_2.cd) < eps) {
        c0_2.cd = f0;
        if (c0_2.sd > f0) {
            c0_2.sd = f1;
        }
        if (c0_2.sd < f0) {
            c0_2.sd = -f1;
        }
    }
    c0_2.sdsd = c0_2.sd * c0_2.sd;
    c0_2.cdcd = c0_2.cd * c0_2.cd;
    c0_2.sdcd = c0_2.sd * c0_2.cd;
    c0_2.s2d = f2 * c0_2.sdcd;
    c0_2.c2d = c0_2.cdcd - c0_2.sdsd;
    return 0;
} /* dccon0_ */

/* Subroutine */ int ua_(doublereal *xi, doublereal *et, doublereal *q, 
        doublereal *disl1, doublereal *disl2, doublereal *disl3, doublereal *
        u)
{
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f2 = 2.;
    static doublereal pi2 = 6.283185307179586;

    static integer i__;
    static doublereal du[12], qx, qy, xy;

/*                                                                       06240000 */
/* ********************************************************************   06250000 */
/* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   06260000 */
/* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   06270000 */
/* ********************************************************************   06280000 */
/*                                                                       06290000 */
/* ***** INPUT                                                            06300000 */
/* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  06310000 */
/* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              06320000 */
/* ***** OUTPUT                                                           06330000 */
/* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     06340000 */
/*                                                                       06350000 */
    /* Parameter adjustments */
    --u;

    /* Function Body */
/* -----                                                                  06400000 */
    for (i__ = 1; i__ <= 12; ++i__) {
/* L111: */
        u[i__] = f0;
    }
    xy = *xi * c2_2.y11;
    qx = *q * c2_2.x11;
    qy = *q * c2_2.y11;
/* ======================================                                 06460000 */
/* =====  STRIKE-SLIP CONTRIBUTION  =====                                 06470000 */
/* ======================================                                 06480000 */
    if (*disl1 != f0) {
        du[0] = c2_2.tt / f2 + c0_2.alp2 * *xi * qy;
        du[1] = c0_2.alp2 * *q / c2_2.r__;
        du[2] = c0_2.alp1 * c2_2.ale - c0_2.alp2 * *q * qy;
        du[3] = -c0_2.alp1 * qy - c0_2.alp2 * c2_2.xi2 * *q * c2_2.y32;
        du[4] = -c0_2.alp2 * *xi * *q / c2_2.r3;
        du[5] = c0_2.alp1 * xy + c0_2.alp2 * *xi * c2_2.q2 * c2_2.y32;
        du[6] = c0_2.alp1 * xy * c0_2.sd + c0_2.alp2 * *xi * c2_2.fy + 
                c2_2.d__ / f2 * c2_2.x11;
        du[7] = c0_2.alp2 * c2_2.ey;
        du[8] = c0_2.alp1 * (c0_2.cd / c2_2.r__ + qy * c0_2.sd) - c0_2.alp2 * 
                *q * c2_2.fy;
        du[9] = c0_2.alp1 * xy * c0_2.cd + c0_2.alp2 * *xi * c2_2.fz + c2_2.y 
                / f2 * c2_2.x11;
        du[10] = c0_2.alp2 * c2_2.ez;
        du[11] = -c0_2.alp1 * (c0_2.sd / c2_2.r__ - qy * c0_2.cd) - c0_2.alp2 
                * *q * c2_2.fz;
        for (i__ = 1; i__ <= 12; ++i__) {
/* L222: */
            u[i__] += *disl1 / pi2 * du[i__ - 1];
        }
    }
/* ======================================                                 06650000 */
/* =====    DIP-SLIP CONTRIBUTION   =====                                 06660000 */
/* ======================================                                 06670000 */
    if (*disl2 != f0) {
        du[0] = c0_2.alp2 * *q / c2_2.r__;
        du[1] = c2_2.tt / f2 + c0_2.alp2 * *et * qx;
        du[2] = c0_2.alp1 * c2_2.alx - c0_2.alp2 * *q * qx;
        du[3] = -c0_2.alp2 * *xi * *q / c2_2.r3;
        du[4] = -qy / f2 - c0_2.alp2 * *et * *q / c2_2.r3;
        du[5] = c0_2.alp1 / c2_2.r__ + c0_2.alp2 * c2_2.q2 / c2_2.r3;
        du[6] = c0_2.alp2 * c2_2.ey;
        du[7] = c0_2.alp1 * c2_2.d__ * c2_2.x11 + xy / f2 * c0_2.sd + 
                c0_2.alp2 * *et * c2_2.gy;
        du[8] = c0_2.alp1 * c2_2.y * c2_2.x11 - c0_2.alp2 * *q * c2_2.gy;
        du[9] = c0_2.alp2 * c2_2.ez;
        du[10] = c0_2.alp1 * c2_2.y * c2_2.x11 + xy / f2 * c0_2.cd + 
                c0_2.alp2 * *et * c2_2.gz;
        du[11] = -c0_2.alp1 * c2_2.d__ * c2_2.x11 - c0_2.alp2 * *q * c2_2.gz;
        for (i__ = 1; i__ <= 12; ++i__) {
/* L333: */
            u[i__] += *disl2 / pi2 * du[i__ - 1];
        }
    }
/* ========================================                               06840000 */
/* =====  TENSILE-FAULT CONTRIBUTION  =====                               06850000 */
/* ========================================                               06860000 */
    if (*disl3 != f0) {
        du[0] = -c0_2.alp1 * c2_2.ale - c0_2.alp2 * *q * qy;
        du[1] = -c0_2.alp1 * c2_2.alx - c0_2.alp2 * *q * qx;
        du[2] = c2_2.tt / f2 - c0_2.alp2 * (*et * qx + *xi * qy);
        du[3] = -c0_2.alp1 * xy + c0_2.alp2 * *xi * c2_2.q2 * c2_2.y32;
        du[4] = -c0_2.alp1 / c2_2.r__ + c0_2.alp2 * c2_2.q2 / c2_2.r3;
        du[5] = -c0_2.alp1 * qy - c0_2.alp2 * *q * c2_2.q2 * c2_2.y32;
        du[6] = -c0_2.alp1 * (c0_2.cd / c2_2.r__ + qy * c0_2.sd) - c0_2.alp2 *
                 *q * c2_2.fy;
        du[7] = -c0_2.alp1 * c2_2.y * c2_2.x11 - c0_2.alp2 * *q * c2_2.gy;
        du[8] = c0_2.alp1 * (c2_2.d__ * c2_2.x11 + xy * c0_2.sd) + c0_2.alp2 *
                 *q * c2_2.hy;
        du[9] = c0_2.alp1 * (c0_2.sd / c2_2.r__ - qy * c0_2.cd) - c0_2.alp2 * 
                *q * c2_2.fz;
        du[10] = c0_2.alp1 * c2_2.d__ * c2_2.x11 - c0_2.alp2 * *q * c2_2.gz;
        du[11] = c0_2.alp1 * (c2_2.y * c2_2.x11 + xy * c0_2.cd) + c0_2.alp2 * 
                *q * c2_2.hz;
        for (i__ = 1; i__ <= 12; ++i__) {
/* L444: */
            u[i__] += *disl3 / pi2 * du[i__ - 1];
        }
    }
    return 0;
} /* ua_ */

/* Subroutine */ int ub_(doublereal *xi, doublereal *et, doublereal *q, 
        doublereal *disl1, doublereal *disl2, doublereal *disl3, doublereal *
        u)
{
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal pi2 = 6.283185307179586;

    /* Local variables */
    static integer i__;
    static doublereal x, d11, rd, du[12], qx, qy, xy, ai1, ai2, aj2, ai4, ai3,
             aj5, ak1, ak3, aj3, aj6, ak2, ak4, aj1, rd2, aj4;

/*                                                                       07080000 */
/* ********************************************************************   07090000 */
/* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   07100000 */
/* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   07110000 */
/* ********************************************************************   07120000 */
/*                                                                       07130000 */
/* ***** INPUT                                                            07140000 */
/* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  07150000 */
/* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              07160000 */
/* ***** OUTPUT                                                           07170000 */
/* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     07180000 */
/*                                                                       07190000 */
    /* Parameter adjustments */
    --u;

    /* Function Body */
/* -----                                                                  07240000 */
    rd = c2_2.r__ + c2_2.d__;
    d11 = f1 / (c2_2.r__ * rd);
    aj2 = *xi * c2_2.y / rd * d11;
    aj5 = -(c2_2.d__ + c2_2.y * c2_2.y / rd) * d11;
    if (c0_2.cd != f0) {
        if (*xi == f0) {
            ai4 = f0;
        } else {
            x = sqrt(c2_2.xi2 + c2_2.q2);
            ai4 = f1 / c0_2.cdcd * (*xi / rd * c0_2.sdcd + f2 * atan((*et * (
                    x + *q * c0_2.cd) + x * (c2_2.r__ + x) * c0_2.sd) / (*xi *
                     (c2_2.r__ + x) * c0_2.cd)));
        }
        ai3 = (c2_2.y * c0_2.cd / rd - c2_2.ale + c0_2.sd * log(rd)) / 
                c0_2.cdcd;
        ak1 = *xi * (d11 - c2_2.y11 * c0_2.sd) / c0_2.cd;
        ak3 = (*q * c2_2.y11 - c2_2.y * d11) / c0_2.cd;
        aj3 = (ak1 - aj2 * c0_2.sd) / c0_2.cd;
        aj6 = (ak3 - aj5 * c0_2.sd) / c0_2.cd;
    } else {
        rd2 = rd * rd;
        ai3 = (*et / rd + c2_2.y * *q / rd2 - c2_2.ale) / f2;
        ai4 = *xi * c2_2.y / rd2 / f2;
        ak1 = *xi * *q / rd * d11;
        ak3 = c0_2.sd / rd * (c2_2.xi2 * d11 - f1);
        aj3 = -(*xi) / rd2 * (c2_2.q2 * d11 - f1 / f2);
        aj6 = -c2_2.y / rd2 * (c2_2.xi2 * d11 - f1 / f2);
    }
/* -----                                                                  07510000 */
    xy = *xi * c2_2.y11;
    ai1 = -(*xi) / rd * c0_2.cd - ai4 * c0_2.sd;
    ai2 = log(rd) + ai3 * c0_2.sd;
    ak2 = f1 / c2_2.r__ + ak3 * c0_2.sd;
    ak4 = xy * c0_2.cd - ak1 * c0_2.sd;
    aj1 = aj5 * c0_2.cd - aj6 * c0_2.sd;
    aj4 = -xy - aj2 * c0_2.cd + aj3 * c0_2.sd;
/* =====                                                                  07590000 */
    for (i__ = 1; i__ <= 12; ++i__) {
/* L111: */
        u[i__] = f0;
    }
    qx = *q * c2_2.x11;
    qy = *q * c2_2.y11;
/* ======================================                                 07640000 */
/* =====  STRIKE-SLIP CONTRIBUTION  =====                                 07650000 */
/* ======================================                                 07660000 */
    if (*disl1 != f0) {
        du[0] = -(*xi) * qy - c2_2.tt - c0_2.alp3 * ai1 * c0_2.sd;
        du[1] = -(*q) / c2_2.r__ + c0_2.alp3 * c2_2.y / rd * c0_2.sd;
        du[2] = *q * qy - c0_2.alp3 * ai2 * c0_2.sd;
        du[3] = c2_2.xi2 * *q * c2_2.y32 - c0_2.alp3 * aj1 * c0_2.sd;
        du[4] = *xi * *q / c2_2.r3 - c0_2.alp3 * aj2 * c0_2.sd;
        du[5] = -(*xi) * c2_2.q2 * c2_2.y32 - c0_2.alp3 * aj3 * c0_2.sd;
        du[6] = -(*xi) * c2_2.fy - c2_2.d__ * c2_2.x11 + c0_2.alp3 * (xy + 
                aj4) * c0_2.sd;
        du[7] = -c2_2.ey + c0_2.alp3 * (f1 / c2_2.r__ + aj5) * c0_2.sd;
        du[8] = *q * c2_2.fy - c0_2.alp3 * (qy - aj6) * c0_2.sd;
        du[9] = -(*xi) * c2_2.fz - c2_2.y * c2_2.x11 + c0_2.alp3 * ak1 * 
                c0_2.sd;
        du[10] = -c2_2.ez + c0_2.alp3 * c2_2.y * d11 * c0_2.sd;
        du[11] = *q * c2_2.fz + c0_2.alp3 * ak2 * c0_2.sd;
        for (i__ = 1; i__ <= 12; ++i__) {
/* L222: */
            u[i__] += *disl1 / pi2 * du[i__ - 1];
        }
    }
/* ======================================                                 07830000 */
/* =====    DIP-SLIP CONTRIBUTION   =====                                 07840000 */
/* ======================================                                 07850000 */
    if (*disl2 != f0) {
        du[0] = -(*q) / c2_2.r__ + c0_2.alp3 * ai3 * c0_2.sdcd;
        du[1] = -(*et) * qx - c2_2.tt - c0_2.alp3 * *xi / rd * c0_2.sdcd;
        du[2] = *q * qx + c0_2.alp3 * ai4 * c0_2.sdcd;
        du[3] = *xi * *q / c2_2.r3 + c0_2.alp3 * aj4 * c0_2.sdcd;
        du[4] = *et * *q / c2_2.r3 + qy + c0_2.alp3 * aj5 * c0_2.sdcd;
        du[5] = -c2_2.q2 / c2_2.r3 + c0_2.alp3 * aj6 * c0_2.sdcd;
        du[6] = -c2_2.ey + c0_2.alp3 * aj1 * c0_2.sdcd;
        du[7] = -(*et) * c2_2.gy - xy * c0_2.sd + c0_2.alp3 * aj2 * c0_2.sdcd;
        du[8] = *q * c2_2.gy + c0_2.alp3 * aj3 * c0_2.sdcd;
        du[9] = -c2_2.ez - c0_2.alp3 * ak3 * c0_2.sdcd;
        du[10] = -(*et) * c2_2.gz - xy * c0_2.cd - c0_2.alp3 * *xi * d11 * 
                c0_2.sdcd;
        du[11] = *q * c2_2.gz - c0_2.alp3 * ak4 * c0_2.sdcd;
        for (i__ = 1; i__ <= 12; ++i__) {
/* L333: */
            u[i__] += *disl2 / pi2 * du[i__ - 1];
        }
    }
/* ========================================                               08020000 */
/* =====  TENSILE-FAULT CONTRIBUTION  =====                               08030000 */
/* ========================================                               08040000 */
    if (*disl3 != f0) {
        du[0] = *q * qy - c0_2.alp3 * ai3 * c0_2.sdsd;
        du[1] = *q * qx + c0_2.alp3 * *xi / rd * c0_2.sdsd;
        du[2] = *et * qx + *xi * qy - c2_2.tt - c0_2.alp3 * ai4 * c0_2.sdsd;
        du[3] = -(*xi) * c2_2.q2 * c2_2.y32 - c0_2.alp3 * aj4 * c0_2.sdsd;
        du[4] = -c2_2.q2 / c2_2.r3 - c0_2.alp3 * aj5 * c0_2.sdsd;
        du[5] = *q * c2_2.q2 * c2_2.y32 - c0_2.alp3 * aj6 * c0_2.sdsd;
        du[6] = *q * c2_2.fy - c0_2.alp3 * aj1 * c0_2.sdsd;
        du[7] = *q * c2_2.gy - c0_2.alp3 * aj2 * c0_2.sdsd;
        du[8] = -(*q) * c2_2.hy - c0_2.alp3 * aj3 * c0_2.sdsd;
        du[9] = *q * c2_2.fz + c0_2.alp3 * ak3 * c0_2.sdsd;
        du[10] = *q * c2_2.gz + c0_2.alp3 * *xi * d11 * c0_2.sdsd;
        du[11] = -(*q) * c2_2.hz + c0_2.alp3 * ak4 * c0_2.sdsd;
        for (i__ = 1; i__ <= 12; ++i__) {
/* L444: */
            u[i__] += *disl3 / pi2 * du[i__ - 1];
        }
    }
    return 0;
} /* ub_ */

/* Subroutine */ int uc_(doublereal *xi, doublereal *et, doublereal *q, 
        doublereal *z__, doublereal *disl1, doublereal *disl2, doublereal *
        disl3, doublereal *u)
{
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal f3 = 3.;
    static doublereal pi2 = 6.283185307179586;

    static doublereal c__, h__;
    static integer i__;
    static doublereal y0, z0, du[12], x53, y53, z32, z53, qq, qx, qy, qr, xy, 
            yy0, cdr, cqx, ppy, ppz, qqy, qqz;

/*                                                                       08260000 */
/* ********************************************************************   08270000 */
/* *****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****   08280000 */
/* *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   08290000 */
/* ********************************************************************   08300000 */
/*                                                                       08310000 */
/* ***** INPUT                                                            08320000 */
/* *****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM              08330000 */
/* *****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              08340000 */
/* ***** OUTPUT                                                           08350000 */
/* *****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     08360000 */
/*                                                                       08370000 */
    /* Parameter adjustments */
    --u;

    /* Function Body */
/* -----                                                                  08420000 */
    c__ = c2_2.d__ + *z__;
    x53 = (c2_2.r2 * 8. + c2_2.r__ * 9. * *xi + f3 * c2_2.xi2) * c2_2.x11 * 
            c2_2.x11 * c2_2.x11 / c2_2.r2;
    y53 = (c2_2.r2 * 8. + c2_2.r__ * 9. * *et + f3 * c2_2.et2) * c2_2.y11 * 
            c2_2.y11 * c2_2.y11 / c2_2.r2;
    h__ = *q * c0_2.cd - *z__;
    z32 = c0_2.sd / c2_2.r3 - h__ * c2_2.y32;
    z53 = f3 * c0_2.sd / c2_2.r5 - h__ * y53;
    y0 = c2_2.y11 - c2_2.xi2 * c2_2.y32;
    z0 = z32 - c2_2.xi2 * z53;
    ppy = c0_2.cd / c2_2.r3 + *q * c2_2.y32 * c0_2.sd;
    ppz = c0_2.sd / c2_2.r3 - *q * c2_2.y32 * c0_2.cd;
    qq = *z__ * c2_2.y32 + z32 + z0;
    qqy = f3 * c__ * c2_2.d__ / c2_2.r5 - qq * c0_2.sd;
    qqz = f3 * c__ * c2_2.y / c2_2.r5 - qq * c0_2.cd + *q * c2_2.y32;
    xy = *xi * c2_2.y11;
    qx = *q * c2_2.x11;
    qy = *q * c2_2.y11;
    qr = f3 * *q / c2_2.r5;
    cqx = c__ * *q * x53;
    cdr = (c__ + c2_2.d__) / c2_2.r3;
    yy0 = c2_2.y / c2_2.r3 - y0 * c0_2.cd;
/* =====                                                                  08630000 */
    for (i__ = 1; i__ <= 12; ++i__) {
/* L111: */
        u[i__] = f0;
    }
/* ======================================                                 08660000 */
/* =====  STRIKE-SLIP CONTRIBUTION  =====                                 08670000 */
/* ======================================                                 08680000 */
    if (*disl1 != f0) {
        du[0] = c0_2.alp4 * xy * c0_2.cd - c0_2.alp5 * *xi * *q * z32;
        du[1] = c0_2.alp4 * (c0_2.cd / c2_2.r__ + f2 * qy * c0_2.sd) - 
                c0_2.alp5 * c__ * *q / c2_2.r3;
        du[2] = c0_2.alp4 * qy * c0_2.cd - c0_2.alp5 * (c__ * *et / c2_2.r3 - 
                *z__ * c2_2.y11 + c2_2.xi2 * z32);
        du[3] = c0_2.alp4 * y0 * c0_2.cd - c0_2.alp5 * *q * z0;
        du[4] = -c0_2.alp4 * *xi * (c0_2.cd / c2_2.r3 + f2 * *q * c2_2.y32 * 
                c0_2.sd) + c0_2.alp5 * c__ * *xi * qr;
        du[5] = -c0_2.alp4 * *xi * *q * c2_2.y32 * c0_2.cd + c0_2.alp5 * *xi *
                 (f3 * c__ * *et / c2_2.r5 - qq);
        du[6] = -c0_2.alp4 * *xi * ppy * c0_2.cd - c0_2.alp5 * *xi * qqy;
        du[7] = c0_2.alp4 * f2 * (c2_2.d__ / c2_2.r3 - y0 * c0_2.sd) * 
                c0_2.sd - c2_2.y / c2_2.r3 * c0_2.cd - c0_2.alp5 * (cdr * 
                c0_2.sd - *et / c2_2.r3 - c__ * c2_2.y * qr);
        du[8] = -c0_2.alp4 * *q / c2_2.r3 + yy0 * c0_2.sd + c0_2.alp5 * (cdr *
                 c0_2.cd + c__ * c2_2.d__ * qr - (y0 * c0_2.cd + *q * z0) * 
                c0_2.sd);
        du[9] = c0_2.alp4 * *xi * ppz * c0_2.cd - c0_2.alp5 * *xi * qqz;
        du[10] = c0_2.alp4 * f2 * (c2_2.y / c2_2.r3 - y0 * c0_2.cd) * c0_2.sd 
                + c2_2.d__ / c2_2.r3 * c0_2.cd - c0_2.alp5 * (cdr * c0_2.cd + 
                c__ * c2_2.d__ * qr);
        du[11] = yy0 * c0_2.cd - c0_2.alp5 * (cdr * c0_2.sd - c__ * c2_2.y * 
                qr - y0 * c0_2.sdsd + *q * z0 * c0_2.cd);
        for (i__ = 1; i__ <= 12; ++i__) {
/* L222: */
            u[i__] += *disl1 / pi2 * du[i__ - 1];
        }
    }
/* ======================================                                 08860000 */
/* =====    DIP-SLIP CONTRIBUTION   =====                                 08870000 */
/* ======================================                                 08880000 */
    if (*disl2 != f0) {
        du[0] = c0_2.alp4 * c0_2.cd / c2_2.r__ - qy * c0_2.sd - c0_2.alp5 * 
                c__ * *q / c2_2.r3;
        du[1] = c0_2.alp4 * c2_2.y * c2_2.x11 - c0_2.alp5 * c__ * *et * *q * 
                c2_2.x32;
        du[2] = -c2_2.d__ * c2_2.x11 - xy * c0_2.sd - c0_2.alp5 * c__ * (
                c2_2.x11 - c2_2.q2 * c2_2.x32);
        du[3] = -c0_2.alp4 * *xi / c2_2.r3 * c0_2.cd + c0_2.alp5 * c__ * *xi *
                 qr + *xi * *q * c2_2.y32 * c0_2.sd;
        du[4] = -c0_2.alp4 * c2_2.y / c2_2.r3 + c0_2.alp5 * c__ * *et * qr;
        du[5] = c2_2.d__ / c2_2.r3 - y0 * c0_2.sd + c0_2.alp5 * c__ / c2_2.r3 
                * (f1 - f3 * c2_2.q2 / c2_2.r2);
        du[6] = -c0_2.alp4 * *et / c2_2.r3 + y0 * c0_2.sdsd - c0_2.alp5 * (
                cdr * c0_2.sd - c__ * c2_2.y * qr);
        du[7] = c0_2.alp4 * (c2_2.x11 - c2_2.y * c2_2.y * c2_2.x32) - 
                c0_2.alp5 * c__ * ((c2_2.d__ + f2 * *q * c0_2.cd) * c2_2.x32 
                - c2_2.y * *et * *q * x53);
        du[8] = *xi * ppy * c0_2.sd + c2_2.y * c2_2.d__ * c2_2.x32 + 
                c0_2.alp5 * c__ * ((c2_2.y + f2 * *q * c0_2.sd) * c2_2.x32 - 
                c2_2.y * c2_2.q2 * x53);
        du[9] = -(*q) / c2_2.r3 + y0 * c0_2.sdcd - c0_2.alp5 * (cdr * c0_2.cd 
                + c__ * c2_2.d__ * qr);
        du[10] = c0_2.alp4 * c2_2.y * c2_2.d__ * c2_2.x32 - c0_2.alp5 * c__ * 
                ((c2_2.y - f2 * *q * c0_2.sd) * c2_2.x32 + c2_2.d__ * *et * *
                q * x53);
        du[11] = -(*xi) * ppz * c0_2.sd + c2_2.x11 - c2_2.d__ * c2_2.d__ * 
                c2_2.x32 - c0_2.alp5 * c__ * ((c2_2.d__ - f2 * *q * c0_2.cd) *
                 c2_2.x32 - c2_2.d__ * c2_2.q2 * x53);
        for (i__ = 1; i__ <= 12; ++i__) {
/* L333: */
            u[i__] += *disl2 / pi2 * du[i__ - 1];
        }
    }
/* ========================================                               09050000 */
/* =====  TENSILE-FAULT CONTRIBUTION  =====                               09060000 */
/* ========================================                               09070000 */
    if (*disl3 != f0) {
        du[0] = -c0_2.alp4 * (c0_2.sd / c2_2.r__ + qy * c0_2.cd) - c0_2.alp5 *
                 (*z__ * c2_2.y11 - c2_2.q2 * z32);
        du[1] = c0_2.alp4 * f2 * xy * c0_2.sd + c2_2.d__ * c2_2.x11 - 
                c0_2.alp5 * c__ * (c2_2.x11 - c2_2.q2 * c2_2.x32);
        du[2] = c0_2.alp4 * (c2_2.y * c2_2.x11 + xy * c0_2.cd) + c0_2.alp5 * *
                q * (c__ * *et * c2_2.x32 + *xi * z32);
        du[3] = c0_2.alp4 * *xi / c2_2.r3 * c0_2.sd + *xi * *q * c2_2.y32 * 
                c0_2.cd + c0_2.alp5 * *xi * (f3 * c__ * *et / c2_2.r5 - f2 * 
                z32 - z0);
        du[4] = c0_2.alp4 * f2 * y0 * c0_2.sd - c2_2.d__ / c2_2.r3 + 
                c0_2.alp5 * c__ / c2_2.r3 * (f1 - f3 * c2_2.q2 / c2_2.r2);
        du[5] = -c0_2.alp4 * yy0 - c0_2.alp5 * (c__ * *et * qr - *q * z0);
        du[6] = c0_2.alp4 * (*q / c2_2.r3 + y0 * c0_2.sdcd) + c0_2.alp5 * (*
                z__ / c2_2.r3 * c0_2.cd + c__ * c2_2.d__ * qr - *q * z0 * 
                c0_2.sd);
        du[7] = -c0_2.alp4 * f2 * *xi * ppy * c0_2.sd - c2_2.y * c2_2.d__ * 
                c2_2.x32 + c0_2.alp5 * c__ * ((c2_2.y + f2 * *q * c0_2.sd) * 
                c2_2.x32 - c2_2.y * c2_2.q2 * x53);
        du[8] = -c0_2.alp4 * (*xi * ppy * c0_2.cd - c2_2.x11 + c2_2.y * 
                c2_2.y * c2_2.x32) + c0_2.alp5 * (c__ * ((c2_2.d__ + f2 * *q *
                 c0_2.cd) * c2_2.x32 - c2_2.y * *et * *q * x53) + *xi * qqy);
        du[9] = -(*et) / c2_2.r3 + y0 * c0_2.cdcd - c0_2.alp5 * (*z__ / 
                c2_2.r3 * c0_2.sd - c__ * c2_2.y * qr - y0 * c0_2.sdsd + *q * 
                z0 * c0_2.cd);
        du[10] = c0_2.alp4 * f2 * *xi * ppz * c0_2.sd - c2_2.x11 + c2_2.d__ * 
                c2_2.d__ * c2_2.x32 - c0_2.alp5 * c__ * ((c2_2.d__ - f2 * *q *
                 c0_2.cd) * c2_2.x32 - c2_2.d__ * c2_2.q2 * x53);
        du[11] = c0_2.alp4 * (*xi * ppz * c0_2.cd + c2_2.y * c2_2.d__ * 
                c2_2.x32) + c0_2.alp5 * (c__ * ((c2_2.y - f2 * *q * c0_2.sd) *
                 c2_2.x32 + c2_2.d__ * *et * *q * x53) + *xi * qqz);
        for (i__ = 1; i__ <= 12; ++i__) {
/* L444: */
            u[i__] += *disl3 / pi2 * du[i__ - 1];
        }
    }
    return 0;
} /* uc_ */

/* Subroutine */ int dccon2_(doublereal *xi, doublereal *et, doublereal *q, 
        doublereal *sd, doublereal *cd)
{
    /* Initialized data */

    static doublereal f0 = 0.;
    static doublereal f1 = 1.;
    static doublereal f2 = 2.;
    static doublereal eps = 1e-6;

    /* Local variables */
    static doublereal ret, rxi;

/*                                                                       10190000 */
/* ********************************************************************** 10200000 */
/* *****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** 10210000 */
/* ********************************************************************** 10220000 */
/*                                                                       10230000 */
/* ***** INPUT                                                            10240000 */
/* *****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  10250000 */
/* *****   SD,CD   : SIN, COS OF DIP-ANGLE                                10260000 */
/*                                                                       10270000 */
/* ### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER010280000 */
/*                                                                       10290000 */
/* -----                                                                  10330000 */
    if (abs(*xi) < eps) {
        *xi = f0;
    }
    if (abs(*et) < eps) {
        *et = f0;
    }
    if (abs(*q) < eps) {
        *q = f0;
    }
    c2_2.xi2 = *xi * *xi;
    c2_2.et2 = *et * *et;
    c2_2.q2 = *q * *q;
    c2_2.r2 = c2_2.xi2 + c2_2.et2 + c2_2.q2;
    c2_2.r__ = sqrt(c2_2.r2);
    if (c2_2.r__ == f0) {
        return 0;
    }
    c2_2.r3 = c2_2.r__ * c2_2.r2;
    c2_2.r5 = c2_2.r3 * c2_2.r2;
    c2_2.y = *et * *cd + *q * *sd;
    c2_2.d__ = *et * *sd - *q * *cd;
/* -----                                                                  10470000 */
    if (*q == f0) {
        c2_2.tt = f0;
    } else {
        c2_2.tt = atan(*xi * *et / (*q * c2_2.r__));
    }
/* -----                                                                  10530000 */
    if (*xi < f0 && *q == f0 && *et == f0) {
        c2_2.alx = -log(c2_2.r__ - *xi);
        c2_2.x11 = f0;
        c2_2.x32 = f0;
    } else {
        rxi = c2_2.r__ + *xi;
        c2_2.alx = log(rxi);
        c2_2.x11 = f1 / (c2_2.r__ * rxi);
        c2_2.x32 = (c2_2.r__ + rxi) * c2_2.x11 * c2_2.x11 / c2_2.r__;
    }
/* -----                                                                  10640000 */
    if (*et < f0 && *q == f0 && *xi == f0) {
        c2_2.ale = -log(c2_2.r__ - *et);
        c2_2.y11 = f0;
        c2_2.y32 = f0;
    } else {
        ret = c2_2.r__ + *et;
        c2_2.ale = log(ret);
        c2_2.y11 = f1 / (c2_2.r__ * ret);
        c2_2.y32 = (c2_2.r__ + ret) * c2_2.y11 * c2_2.y11 / c2_2.r__;
    }
/* -----                                                                  10750000 */
    c2_2.ey = *sd / c2_2.r__ - c2_2.y * *q / c2_2.r3;
    c2_2.ez = *cd / c2_2.r__ + c2_2.d__ * *q / c2_2.r3;
    c2_2.fy = c2_2.d__ / c2_2.r3 + c2_2.xi2 * c2_2.y32 * *sd;
    c2_2.fz = c2_2.y / c2_2.r3 + c2_2.xi2 * c2_2.y32 * *cd;
    c2_2.gy = f2 * c2_2.x11 * *sd - c2_2.y * *q * c2_2.x32;
    c2_2.gz = f2 * c2_2.x11 * *cd + c2_2.d__ * *q * c2_2.x32;
    c2_2.hy = c2_2.d__ * *q * c2_2.x32 + *xi * *q * c2_2.y32 * *sd;
    c2_2.hz = c2_2.y * *q * c2_2.x32 + *xi * *q * c2_2.y32 * *cd;
    return 0;
} /* dccon2_ */

