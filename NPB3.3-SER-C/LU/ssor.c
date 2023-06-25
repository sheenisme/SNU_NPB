//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB LU code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

#include "applu.incl"
#include "timers.h"
#include <math.h>
#include <stdio.h>

//---------------------------------------------------------------------
// to perform pseudo-time stepping SSOR iterations
// for five nonlinear pde's.
//---------------------------------------------------------------------
void ssor(int niter) {
    //---------------------------------------------------------------------
    // local variables
    //---------------------------------------------------------------------
    int i, j, k, m, n;
    int istep;
    double tmp_ssor, tv[ISIZ2][ISIZ1][5];
    double delunm[5];

    //---------------------------------------------------------------------
    // begin pseudo-time stepping iterations
    //---------------------------------------------------------------------
    tmp_ssor = 1.0 / (omega * (2.0 - omega));

    //---------------------------------------------------------------------
    // initialize a,b,c,d to zero (guarantees that page tables have been
    // formed, if applicable on given architecture, before timestepping).
    //---------------------------------------------------------------------
    for (j = 0; j < ISIZ2; j++) {
        for (i = 0; i < ISIZ1; i++) {
            for (n = 0; n < 5; n++) {
                for (m = 0; m < 5; m++) {
                    a[j][i][n][m] = 0.0;
                    b[j][i][n][m] = 0.0;
                    c[j][i][n][m] = 0.0;
                    d[j][i][n][m] = 0.0;
                }
            }
        }
    }
    for (i = 1; i <= t_last; i++) {
        timer_clear(i);
    }

    //---------------------------------------------------------------------
    // compute the steady-state residuals
    //---------------------------------------------------------------------
    rhs();

    //---------------------------------------------------------------------
    // compute the L2 norms of newton iteration residuals
    //---------------------------------------------------------------------
    l2norm(ISIZ1, ISIZ2, ISIZ3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);

    /*
  if ( ipr == 1 ) {
    printf("           Initial residual norms\n");
    printf("\n");
    printf(" \n RMS-norm of steady-state residual for "
           "first pde  = %12.5E\n"
           " RMS-norm of steady-state residual for "
           "second pde = %12.5E\n"
           " RMS-norm of steady-state residual for "
           "third pde  = %12.5E\n"
           " RMS-norm of steady-state residual for "
           "fourth pde = %12.5E\n"
           " RMS-norm of steady-state residual for "
           "fifth pde  = %12.5E\n",
           rsdnm[0], rsdnm[1], rsdnm[2], rsdnm[3], rsdnm[4]);
    printf("\nIteration RMS-residual of 5th PDE\n");
  }
  */

    for (i = 1; i <= t_last; i++) {
        timer_clear(i);
    }
    timer_start(1);

    //---------------------------------------------------------------------
    // the timestep loop
    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    // local variables
    //---------------------------------------------------------------------
    double r43;
    double c1345;
    double c34;
    double tmp, tmp1, tmp2, tmp3;
    double tmat[5][5], tv_tmp[5];
    double zero = 0.0;
    double one = 1.0;
    double two = 2.0;
    double three = 3.0;
    double four = 4.0;
    double cc1 = C1;
    double cc2 = C2;
    double cc3 = C3;
    double cc4 = C4;
    double cc5 = C5;

#pragma scop
    // // __Pencil set for Scale-W
    // __pencil_assume(timeron >= 1);
    // __pencil_assume(jst >= 1);
    // __pencil_assume(jend >= 32);
    // __pencil_assume(ist >= 1);
    // __pencil_assume(iend >= 32);
    // __pencil_assume(nz >= 33);
    // __pencil_assume(timeron <= 1);
    // __pencil_assume(jst <= 1);
    // __pencil_assume(jend <= 32);
    // __pencil_assume(ist <= 1);
    // __pencil_assume(iend <= 32);
    // __pencil_assume(nz <= 33);

    // __Pencil set for Scale-A
    __pencil_assume(timeron >= 1);
    __pencil_assume(jst >= 1);
    __pencil_assume(jend >= 63);
    __pencil_assume(ist >= 1);
    __pencil_assume(iend >= 63);
    __pencil_assume(nz >= 64);
    __pencil_assume(timeron <= 1);
    __pencil_assume(jst <= 1);
    __pencil_assume(jend <= 63);
    __pencil_assume(ist <= 1);
    __pencil_assume(iend <= 63);
    __pencil_assume(nz <= 64);

    // // __Pencil set for Scale-B
    // __pencil_assume(timeron >= 1);
    // __pencil_assume(jst >= 1);
    // __pencil_assume(jend >= 101);
    // __pencil_assume(ist >= 1);
    // __pencil_assume(iend >= 101);
    // __pencil_assume(nz >= 102);
    // __pencil_assume(timeron <= 1);
    // __pencil_assume(jst <= 1);
    // __pencil_assume(jend <= 101);
    // __pencil_assume(ist <= 1);
    // __pencil_assume(iend <= 101);
    // __pencil_assume(nz <= 102);

    // // __Pencil set for Scale-C
    // __pencil_assume(timeron >= 1);
    // __pencil_assume(jst >= 1);
    // __pencil_assume(jend >= 161);
    // __pencil_assume(ist >= 1);
    // __pencil_assume(iend >= 161);
    // __pencil_assume(nz >= 162);
    // __pencil_assume(timeron <= 1);
    // __pencil_assume(jst <= 1);
    // __pencil_assume(jend <= 161);
    // __pencil_assume(ist <= 1);
    // __pencil_assume(iend <= 161);
    // __pencil_assume(nz <= 162);

    __pencil_assume(niter >= 1);
    
    // printf("%d, %d,%d, %d,%d, %d, %d\n", timeron, jst, jend, ist,iend, niter, nz);

    for (int istep = 1; istep <= niter; istep++) {
        // if ( ( (istep % inorm) == 0 ) && ipr == 1 ) {
        //   printf(" \n     pseudo-time SSOR iteration no.=%4d\n\n", istep);
        // }
        //  if ((istep % 20) == 0 || istep == itmax || istep == 1) {
        //    if (niter > 1) printf(" Time step %4d\n", istep);
        //  }

        //---------------------------------------------------------------------
        // perform SSOR iteration
        //---------------------------------------------------------------------
        if (timeron)
            timer_start(t_rhs);

        // #pragma scop
        for (k = 1; k < nz - 1; k++) {
            for (j = jst; j < jend; j++) {
                for (i = ist; i < iend; i++) {
                    for (m = 0; m < 5; m++) {
                        rsd[k][j][i][m] = dt * rsd[k][j][i][m];
                    }
                }
            }
        }
        // #pragma endscop

        if (timeron)
            timer_stop(t_rhs);

        for (k = 1; k < nz - 1; k++) {
            //---------------------------------------------------------------------
            // form the lower triangular part of the jacobian matrix
            //---------------------------------------------------------------------
            if (timeron)
                timer_start(t_jacld);

            // jacld(k);
            {
                //---------------------------------------------------------------------
                // local variables
                //---------------------------------------------------------------------
                // int i, j;
                // double r43;
                // double c1345;
                // double c34;
                // double tmp1, tmp2, tmp3;

                r43 = (4.0 / 3.0);
                c1345 = cc1 * cc3 * cc4 * cc5;
                c34 = cc3 * cc4;

                // #pragma scop
                for (j = jst; j < jend; j++) {
                    for (i = ist; i < iend; i++) {
                        //---------------------------------------------------------------------
                        // form the block daigonal
                        //---------------------------------------------------------------------
                        tmp1 = rho_i[k][j][i];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        d[j][i][0][0] = one + dt * two * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
                        d[j][i][1][0] = zero;
                        d[j][i][2][0] = zero;
                        d[j][i][3][0] = zero;
                        d[j][i][4][0] = zero;

                        d[j][i][0][1] = -dt * two * (tx1 * r43 + ty1 + tz1) * c34 * tmp2 * u[k][j][i][1];
                        d[j][i][1][1] = one + dt * two * c34 * tmp1 * (tx1 * r43 + ty1 + tz1) + dt * two * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
                        d[j][i][2][1] = zero;
                        d[j][i][3][1] = zero;
                        d[j][i][4][1] = zero;

                        d[j][i][0][2] = -dt * two * (tx1 + ty1 * r43 + tz1) * c34 * tmp2 * u[k][j][i][2];
                        d[j][i][1][2] = zero;
                        d[j][i][2][2] = one + dt * two * c34 * tmp1 * (tx1 + ty1 * r43 + tz1) + dt * two * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
                        d[j][i][3][2] = zero;
                        d[j][i][4][2] = zero;

                        d[j][i][0][3] = -dt * two * (tx1 + ty1 + tz1 * r43) * c34 * tmp2 * u[k][j][i][3];
                        d[j][i][1][3] = zero;
                        d[j][i][2][3] = zero;
                        d[j][i][3][3] = one + dt * two * c34 * tmp1 * (tx1 + ty1 + tz1 * r43) + dt * two * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
                        d[j][i][4][3] = zero;

                        d[j][i][0][4] = -dt * two * (((tx1 * (r43 * c34 - c1345) + ty1 * (c34 - c1345) + tz1 * (c34 - c1345)) * (u[k][j][i][1] * u[k][j][i][1]) + (tx1 * (c34 - c1345) + ty1 * (r43 * c34 - c1345) + tz1 * (c34 - c1345)) * (u[k][j][i][2] * u[k][j][i][2]) + (tx1 * (c34 - c1345) + ty1 * (c34 - c1345) + tz1 * (r43 * c34 - c1345)) * (u[k][j][i][3] * u[k][j][i][3])) * tmp3 + (tx1 + ty1 + tz1) * c1345 * tmp2 * u[k][j][i][4]);
                        d[j][i][1][4] = dt * two * tmp2 * u[k][j][i][1] * (tx1 * (r43 * c34 - c1345) + ty1 * (c34 - c1345) + tz1 * (c34 - c1345));
                        d[j][i][2][4] = dt * two * tmp2 * u[k][j][i][2] * (tx1 * (c34 - c1345) + ty1 * (r43 * c34 - c1345) + tz1 * (c34 - c1345));
                        d[j][i][3][4] = dt * two * tmp2 * u[k][j][i][3] * (tx1 * (c34 - c1345) + ty1 * (c34 - c1345) + tz1 * (r43 * c34 - c1345));
                        d[j][i][4][4] = one + dt * two * (tx1 + ty1 + tz1) * c1345 * tmp1 + dt * two * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);

                        //---------------------------------------------------------------------
                        // form the first block sub-diagonal
                        //---------------------------------------------------------------------
                        tmp1 = rho_i[k - 1][j][i];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        a[j][i][0][0] = -dt * tz1 * dz1;
                        a[j][i][1][0] = zero;
                        a[j][i][2][0] = zero;
                        a[j][i][3][0] = -dt * tz2;
                        a[j][i][4][0] = zero;

                        a[j][i][0][1] = -dt * tz2 * (-(u[k - 1][j][i][1] * u[k - 1][j][i][3]) * tmp2) - dt * tz1 * (-c34 * tmp2 * u[k - 1][j][i][1]);
                        a[j][i][1][1] = -dt * tz2 * (u[k - 1][j][i][3] * tmp1) - dt * tz1 * c34 * tmp1 - dt * tz1 * dz2;
                        a[j][i][2][1] = zero;
                        a[j][i][3][1] = -dt * tz2 * (u[k - 1][j][i][1] * tmp1);
                        a[j][i][4][1] = zero;

                        a[j][i][0][2] = -dt * tz2 * (-(u[k - 1][j][i][2] * u[k - 1][j][i][3]) * tmp2) - dt * tz1 * (-c34 * tmp2 * u[k - 1][j][i][2]);
                        a[j][i][1][2] = zero;
                        a[j][i][2][2] = -dt * tz2 * (u[k - 1][j][i][3] * tmp1) - dt * tz1 * (c34 * tmp1) - dt * tz1 * dz3;
                        a[j][i][3][2] = -dt * tz2 * (u[k - 1][j][i][2] * tmp1);
                        a[j][i][4][2] = zero;

                        a[j][i][0][3] = -dt * tz2 * (-(u[k - 1][j][i][3] * tmp1) * (u[k - 1][j][i][3] * tmp1) + cc2 * qs[k - 1][j][i] * tmp1) - dt * tz1 * (-r43 * c34 * tmp2 * u[k - 1][j][i][3]);
                        a[j][i][1][3] = -dt * tz2 * (-cc2 * (u[k - 1][j][i][1] * tmp1));
                        a[j][i][2][3] = -dt * tz2 * (-cc2 * (u[k - 1][j][i][2] * tmp1));
                        a[j][i][3][3] = -dt * tz2 * (two - cc2) * (u[k - 1][j][i][3] * tmp1) - dt * tz1 * (r43 * c34 * tmp1) - dt * tz1 * dz4;
                        a[j][i][4][3] = -dt * tz2 * cc2;

                        a[j][i][0][4] = -dt * tz2 * ((cc2 * two * qs[k - 1][j][i] - cc1 * u[k - 1][j][i][4]) * u[k - 1][j][i][3] * tmp2) - dt * tz1 * (-(c34 - c1345) * tmp3 * (u[k - 1][j][i][1] * u[k - 1][j][i][1]) - (c34 - c1345) * tmp3 * (u[k - 1][j][i][2] * u[k - 1][j][i][2]) - (r43 * c34 - c1345) * tmp3 * (u[k - 1][j][i][3] * u[k - 1][j][i][3]) - c1345 * tmp2 * u[k - 1][j][i][4]);
                        a[j][i][1][4] = -dt * tz2 * (-cc2 * (u[k - 1][j][i][1] * u[k - 1][j][i][3]) * tmp2) - dt * tz1 * (c34 - c1345) * tmp2 * u[k - 1][j][i][1];
                        a[j][i][2][4] = -dt * tz2 * (-cc2 * (u[k - 1][j][i][2] * u[k - 1][j][i][3]) * tmp2) - dt * tz1 * (c34 - c1345) * tmp2 * u[k - 1][j][i][2];
                        a[j][i][3][4] = -dt * tz2 * (cc1 * (u[k - 1][j][i][4] * tmp1) - cc2 * (qs[k - 1][j][i] * tmp1 + u[k - 1][j][i][3] * u[k - 1][j][i][3] * tmp2)) - dt * tz1 * (r43 * c34 - c1345) * tmp2 * u[k - 1][j][i][3];
                        a[j][i][4][4] = -dt * tz2 * (cc1 * (u[k - 1][j][i][3] * tmp1)) - dt * tz1 * c1345 * tmp1 - dt * tz1 * dz5;

                        //---------------------------------------------------------------------
                        // form the second block sub-diagonal
                        //---------------------------------------------------------------------
                        tmp1 = rho_i[k][j - 1][i];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        b[j][i][0][0] = -dt * ty1 * dy1;
                        b[j][i][1][0] = zero;
                        b[j][i][2][0] = -dt * ty2;
                        b[j][i][3][0] = zero;
                        b[j][i][4][0] = zero;

                        b[j][i][0][1] = -dt * ty2 * (-(u[k][j - 1][i][1] * u[k][j - 1][i][2]) * tmp2) - dt * ty1 * (-c34 * tmp2 * u[k][j - 1][i][1]);
                        b[j][i][1][1] = -dt * ty2 * (u[k][j - 1][i][2] * tmp1) - dt * ty1 * (c34 * tmp1) - dt * ty1 * dy2;
                        b[j][i][2][1] = -dt * ty2 * (u[k][j - 1][i][1] * tmp1);
                        b[j][i][3][1] = zero;
                        b[j][i][4][1] = zero;

                        b[j][i][0][2] = -dt * ty2 * (-(u[k][j - 1][i][2] * tmp1) * (u[k][j - 1][i][2] * tmp1) + cc2 * (qs[k][j - 1][i] * tmp1)) - dt * ty1 * (-r43 * c34 * tmp2 * u[k][j - 1][i][2]);
                        b[j][i][1][2] = -dt * ty2 * (-cc2 * (u[k][j - 1][i][1] * tmp1));
                        b[j][i][2][2] = -dt * ty2 * ((two - cc2) * (u[k][j - 1][i][2] * tmp1)) - dt * ty1 * (r43 * c34 * tmp1) - dt * ty1 * dy3;
                        b[j][i][3][2] = -dt * ty2 * (-cc2 * (u[k][j - 1][i][3] * tmp1));
                        b[j][i][4][2] = -dt * ty2 * cc2;

                        b[j][i][0][3] = -dt * ty2 * (-(u[k][j - 1][i][2] * u[k][j - 1][i][3]) * tmp2) - dt * ty1 * (-c34 * tmp2 * u[k][j - 1][i][3]);
                        b[j][i][1][3] = zero;
                        b[j][i][2][3] = -dt * ty2 * (u[k][j - 1][i][3] * tmp1);
                        b[j][i][3][3] = -dt * ty2 * (u[k][j - 1][i][2] * tmp1) - dt * ty1 * (c34 * tmp1) - dt * ty1 * dy4;
                        b[j][i][4][3] = zero;

                        b[j][i][0][4] = -dt * ty2 * ((cc2 * two * qs[k][j - 1][i] - cc1 * u[k][j - 1][i][4]) * (u[k][j - 1][i][2] * tmp2)) - dt * ty1 * (-(c34 - c1345) * tmp3 * (u[k][j - 1][i][1] * u[k][j - 1][i][1]) - (r43 * c34 - c1345) * tmp3 * (u[k][j - 1][i][2] * u[k][j - 1][i][2]) - (c34 - c1345) * tmp3 * (u[k][j - 1][i][3] * u[k][j - 1][i][3]) - c1345 * tmp2 * u[k][j - 1][i][4]);
                        b[j][i][1][4] = -dt * ty2 * (-cc2 * (u[k][j - 1][i][1] * u[k][j - 1][i][2]) * tmp2) - dt * ty1 * (c34 - c1345) * tmp2 * u[k][j - 1][i][1];
                        b[j][i][2][4] = -dt * ty2 * (cc1 * (u[k][j - 1][i][4] * tmp1) - cc2 * (qs[k][j - 1][i] * tmp1 + u[k][j - 1][i][2] * u[k][j - 1][i][2] * tmp2)) - dt * ty1 * (r43 * c34 - c1345) * tmp2 * u[k][j - 1][i][2];
                        b[j][i][3][4] = -dt * ty2 * (-cc2 * (u[k][j - 1][i][2] * u[k][j - 1][i][3]) * tmp2) - dt * ty1 * (c34 - c1345) * tmp2 * u[k][j - 1][i][3];
                        b[j][i][4][4] = -dt * ty2 * (cc1 * (u[k][j - 1][i][2] * tmp1)) - dt * ty1 * c1345 * tmp1 - dt * ty1 * dy5;

                        //---------------------------------------------------------------------
                        // form the third block sub-diagonal
                        //---------------------------------------------------------------------
                        tmp1 = rho_i[k][j][i - 1];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        c[j][i][0][0] = -dt * tx1 * dx1;
                        c[j][i][1][0] = -dt * tx2;
                        c[j][i][2][0] = zero;
                        c[j][i][3][0] = zero;
                        c[j][i][4][0] = zero;

                        c[j][i][0][1] = -dt * tx2 * (-(u[k][j][i - 1][1] * tmp1) * (u[k][j][i - 1][1] * tmp1) + cc2 * qs[k][j][i - 1] * tmp1) - dt * tx1 * (-r43 * c34 * tmp2 * u[k][j][i - 1][1]);
                        c[j][i][1][1] = -dt * tx2 * ((two - cc2) * (u[k][j][i - 1][1] * tmp1)) - dt * tx1 * (r43 * c34 * tmp1) - dt * tx1 * dx2;
                        c[j][i][2][1] = -dt * tx2 * (-cc2 * (u[k][j][i - 1][2] * tmp1));
                        c[j][i][3][1] = -dt * tx2 * (-cc2 * (u[k][j][i - 1][3] * tmp1));
                        c[j][i][4][1] = -dt * tx2 * cc2;

                        c[j][i][0][2] = -dt * tx2 * (-(u[k][j][i - 1][1] * u[k][j][i - 1][2]) * tmp2) - dt * tx1 * (-c34 * tmp2 * u[k][j][i - 1][2]);
                        c[j][i][1][2] = -dt * tx2 * (u[k][j][i - 1][2] * tmp1);
                        c[j][i][2][2] = -dt * tx2 * (u[k][j][i - 1][1] * tmp1) - dt * tx1 * (c34 * tmp1) - dt * tx1 * dx3;
                        c[j][i][3][2] = zero;
                        c[j][i][4][2] = zero;

                        c[j][i][0][3] = -dt * tx2 * (-(u[k][j][i - 1][1] * u[k][j][i - 1][3]) * tmp2) - dt * tx1 * (-c34 * tmp2 * u[k][j][i - 1][3]);
                        c[j][i][1][3] = -dt * tx2 * (u[k][j][i - 1][3] * tmp1);
                        c[j][i][2][3] = zero;
                        c[j][i][3][3] = -dt * tx2 * (u[k][j][i - 1][1] * tmp1) - dt * tx1 * (c34 * tmp1) - dt * tx1 * dx4;
                        c[j][i][4][3] = zero;

                        c[j][i][0][4] = -dt * tx2 * ((cc2 * two * qs[k][j][i - 1] - cc1 * u[k][j][i - 1][4]) * u[k][j][i - 1][1] * tmp2) - dt * tx1 * (-(r43 * c34 - c1345) * tmp3 * (u[k][j][i - 1][1] * u[k][j][i - 1][1]) - (c34 - c1345) * tmp3 * (u[k][j][i - 1][2] * u[k][j][i - 1][2]) - (c34 - c1345) * tmp3 * (u[k][j][i - 1][3] * u[k][j][i - 1][3]) - c1345 * tmp2 * u[k][j][i - 1][4]);
                        c[j][i][1][4] = -dt * tx2 * (cc1 * (u[k][j][i - 1][4] * tmp1) - cc2 * (u[k][j][i - 1][1] * u[k][j][i - 1][1] * tmp2 + qs[k][j][i - 1] * tmp1)) - dt * tx1 * (r43 * c34 - c1345) * tmp2 * u[k][j][i - 1][1];
                        c[j][i][2][4] = -dt * tx2 * (-cc2 * (u[k][j][i - 1][2] * u[k][j][i - 1][1]) * tmp2) - dt * tx1 * (c34 - c1345) * tmp2 * u[k][j][i - 1][2];
                        c[j][i][3][4] = -dt * tx2 * (-cc2 * (u[k][j][i - 1][3] * u[k][j][i - 1][1]) * tmp2) - dt * tx1 * (c34 - c1345) * tmp2 * u[k][j][i - 1][3];
                        c[j][i][4][4] = -dt * tx2 * (cc1 * (u[k][j][i - 1][1] * tmp1)) - dt * tx1 * c1345 * tmp1 - dt * tx1 * dx5;
                    }
                }
                // #pragma endscop
            }

            if (timeron)
                timer_stop(t_jacld);

            //---------------------------------------------------------------------
            // perform the lower triangular solution
            //---------------------------------------------------------------------
            if (timeron)
                timer_start(t_blts);

            /* blts( ISIZ1, ISIZ2, ISIZ3,
            nx, ny, nz, k,
            omega,
            rsd,
            a, b, c, d,
            ist, iend, jst, jend,
            nx0, ny0 );
      */
            /* void blts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, int k,
                double omega, double v[][ldmy / 2 * 2 + 1][ldmx / 2 * 2 + 1][5],
                double ldz[ldmy][ldmx / 2 * 2 + 1][5][5],
                double ldy[ldmy][ldmx / 2 * 2 + 1][5][5],
                double ldx[ldmy][ldmx / 2 * 2 + 1][5][5],rsd[k]
                double d[ldmy][ldmx / 2 * 2 + 1][5][5], int ist, int iend,
                int jst, int jend, int nx0, int ny0)
      */
            {
                //---------------------------------------------------------------------
                // local variables
                //---------------------------------------------------------------------
                // int i, j, m;
                // double tmp, tmp1;
                // double tmat[5][5], tv[5];

                // #pragma scop
                for (j = jst; j < jend; j++) {
                    for (i = ist; i < iend; i++) {
                        for (m = 0; m < 5; m++) {
                            rsd[k][j][i][m] = rsd[k][j][i][m] - omega * (a[j][i][0][m] * rsd[k - 1][j][i][0] + a[j][i][1][m] * rsd[k - 1][j][i][1] + a[j][i][2][m] * rsd[k - 1][j][i][2] + a[j][i][3][m] * rsd[k - 1][j][i][3] + a[j][i][4][m] * rsd[k - 1][j][i][4]);
                        }
                    }
                }
                // #pragma endscop

                for (j = jst; j < jend; j++) {
                    for (i = ist; i < iend; i++) {
                        for (m = 0; m < 5; m++) {
                            tv_tmp[m] = rsd[k][j][i][m] - omega * (b[j][i][0][m] * rsd[k][j - 1][i][0] + c[j][i][0][m] * rsd[k][j][i - 1][0] + b[j][i][1][m] * rsd[k][j - 1][i][1] + c[j][i][1][m] * rsd[k][j][i - 1][1] + b[j][i][2][m] * rsd[k][j - 1][i][2] + c[j][i][2][m] * rsd[k][j][i - 1][2] + b[j][i][3][m] * rsd[k][j - 1][i][3] + c[j][i][3][m] * rsd[k][j][i - 1][3] + b[j][i][4][m] * rsd[k][j - 1][i][4] + c[j][i][4][m] * rsd[k][j][i - 1][4]);
                        }

                        //---------------------------------------------------------------------
                        // diagonal block inversion
                        //
                        // forward elimination
                        //---------------------------------------------------------------------
                        for (m = 0; m < 5; m++) {
                            tmat[m][0] = d[j][i][0][m];
                            tmat[m][1] = d[j][i][1][m];
                            tmat[m][2] = d[j][i][2][m];
                            tmat[m][3] = d[j][i][3][m];
                            tmat[m][4] = d[j][i][4][m];
                        }

                        tmp1 = one / tmat[0][0];
                        tmp = tmp1 * tmat[1][0];
                        tmat[1][1] = tmat[1][1] - tmp * tmat[0][1];
                        tmat[1][2] = tmat[1][2] - tmp * tmat[0][2];
                        tmat[1][3] = tmat[1][3] - tmp * tmat[0][3];
                        tmat[1][4] = tmat[1][4] - tmp * tmat[0][4];
                        tv_tmp[1] = tv_tmp[1] - tv_tmp[0] * tmp;

                        tmp = tmp1 * tmat[2][0];
                        tmat[2][1] = tmat[2][1] - tmp * tmat[0][1];
                        tmat[2][2] = tmat[2][2] - tmp * tmat[0][2];
                        tmat[2][3] = tmat[2][3] - tmp * tmat[0][3];
                        tmat[2][4] = tmat[2][4] - tmp * tmat[0][4];
                        tv_tmp[2] = tv_tmp[2] - tv_tmp[0] * tmp;

                        tmp = tmp1 * tmat[3][0];
                        tmat[3][1] = tmat[3][1] - tmp * tmat[0][1];
                        tmat[3][2] = tmat[3][2] - tmp * tmat[0][2];
                        tmat[3][3] = tmat[3][3] - tmp * tmat[0][3];
                        tmat[3][4] = tmat[3][4] - tmp * tmat[0][4];
                        tv_tmp[3] = tv_tmp[3] - tv_tmp[0] * tmp;

                        tmp = tmp1 * tmat[4][0];
                        tmat[4][1] = tmat[4][1] - tmp * tmat[0][1];
                        tmat[4][2] = tmat[4][2] - tmp * tmat[0][2];
                        tmat[4][3] = tmat[4][3] - tmp * tmat[0][3];
                        tmat[4][4] = tmat[4][4] - tmp * tmat[0][4];
                        tv_tmp[4] = tv_tmp[4] - tv_tmp[0] * tmp;

                        tmp1 = one / tmat[1][1];
                        tmp = tmp1 * tmat[2][1];
                        tmat[2][2] = tmat[2][2] - tmp * tmat[1][2];
                        tmat[2][3] = tmat[2][3] - tmp * tmat[1][3];
                        tmat[2][4] = tmat[2][4] - tmp * tmat[1][4];
                        tv_tmp[2] = tv_tmp[2] - tv_tmp[1] * tmp;

                        tmp = tmp1 * tmat[3][1];
                        tmat[3][2] = tmat[3][2] - tmp * tmat[1][2];
                        tmat[3][3] = tmat[3][3] - tmp * tmat[1][3];
                        tmat[3][4] = tmat[3][4] - tmp * tmat[1][4];
                        tv_tmp[3] = tv_tmp[3] - tv_tmp[1] * tmp;

                        tmp = tmp1 * tmat[4][1];
                        tmat[4][2] = tmat[4][2] - tmp * tmat[1][2];
                        tmat[4][3] = tmat[4][3] - tmp * tmat[1][3];
                        tmat[4][4] = tmat[4][4] - tmp * tmat[1][4];
                        tv_tmp[4] = tv_tmp[4] - tv_tmp[1] * tmp;

                        tmp1 = one / tmat[2][2];
                        tmp = tmp1 * tmat[3][2];
                        tmat[3][3] = tmat[3][3] - tmp * tmat[2][3];
                        tmat[3][4] = tmat[3][4] - tmp * tmat[2][4];
                        tv_tmp[3] = tv_tmp[3] - tv_tmp[2] * tmp;

                        tmp = tmp1 * tmat[4][2];
                        tmat[4][3] = tmat[4][3] - tmp * tmat[2][3];
                        tmat[4][4] = tmat[4][4] - tmp * tmat[2][4];
                        tv_tmp[4] = tv_tmp[4] - tv_tmp[2] * tmp;

                        tmp1 = one / tmat[3][3];
                        tmp = tmp1 * tmat[4][3];
                        tmat[4][4] = tmat[4][4] - tmp * tmat[3][4];
                        tv_tmp[4] = tv_tmp[4] - tv_tmp[3] * tmp;

                        //---------------------------------------------------------------------
                        // back substitution
                        //---------------------------------------------------------------------
                        rsd[k][j][i][4] = tv_tmp[4] / tmat[4][4];

                        tv_tmp[3] = tv_tmp[3] - tmat[3][4] * rsd[k][j][i][4];
                        rsd[k][j][i][3] = tv_tmp[3] / tmat[3][3];

                        tv_tmp[2] = tv_tmp[2] - tmat[2][3] * rsd[k][j][i][3] - tmat[2][4] * rsd[k][j][i][4];
                        rsd[k][j][i][2] = tv_tmp[2] / tmat[2][2];

                        tv_tmp[1] = tv_tmp[1] - tmat[1][2] * rsd[k][j][i][2] - tmat[1][3] * rsd[k][j][i][3] - tmat[1][4] * rsd[k][j][i][4];
                        rsd[k][j][i][1] = tv_tmp[1] / tmat[1][1];

                        tv_tmp[0] = tv_tmp[0] - tmat[0][1] * rsd[k][j][i][1] - tmat[0][2] * rsd[k][j][i][2] - tmat[0][3] * rsd[k][j][i][3] - tmat[0][4] * rsd[k][j][i][4];
                        rsd[k][j][i][0] = tv_tmp[0] / tmat[0][0];
                    }
                }
            }

            if (timeron)
                timer_stop(t_blts);
        }

        for (k = nz - 2; k > 0; k--) {
            //---------------------------------------------------------------------
            // form the strictly upper triangular part of the jacobian matrix
            //---------------------------------------------------------------------
            if (timeron)
                timer_start(t_jacu);

            // jacu(k);
            // void jacu(int k)
            {
                //---------------------------------------------------------------------
                // local variables
                //---------------------------------------------------------------------
                // int i, j;
                // double r43;
                // double c1345;
                // double c34;
                // double tmp1, tmp2, tmp3;

                r43 = (4.0 / 3.0);
                c1345 = cc1 * cc3 * cc4 * cc5;
                c34 = cc3 * cc4;
                for (j = jst; j < jend; j++) {
                    for (i = ist; i < iend; i++) {
                        //---------------------------------------------------------------------
                        // form the block daigonal
                        //---------------------------------------------------------------------
                        tmp1 = rho_i[k][j][i];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        d[j][i][0][0] = one + dt * two * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
                        d[j][i][1][0] = zero;
                        d[j][i][2][0] = zero;
                        d[j][i][3][0] = zero;
                        d[j][i][4][0] = zero;

                        d[j][i][0][1] = dt * two * (-tx1 * r43 - ty1 - tz1) * (c34 * tmp2 * u[k][j][i][1]);
                        d[j][i][1][1] = one + dt * two * c34 * tmp1 * (tx1 * r43 + ty1 + tz1) + dt * two * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
                        d[j][i][2][1] = zero;
                        d[j][i][3][1] = zero;
                        d[j][i][4][1] = zero;

                        d[j][i][0][2] = dt * two * (-tx1 - ty1 * r43 - tz1) * (c34 * tmp2 * u[k][j][i][2]);
                        d[j][i][1][2] = zero;
                        d[j][i][2][2] = one + dt * two * c34 * tmp1 * (tx1 + ty1 * r43 + tz1) + dt * two * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
                        d[j][i][3][2] = zero;
                        d[j][i][4][2] = zero;

                        d[j][i][0][3] = dt * two * (-tx1 - ty1 - tz1 * r43) * (c34 * tmp2 * u[k][j][i][3]);
                        d[j][i][1][3] = zero;
                        d[j][i][2][3] = zero;
                        d[j][i][3][3] = one + dt * two * c34 * tmp1 * (tx1 + ty1 + tz1 * r43) + dt * two * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
                        d[j][i][4][3] = zero;

                        d[j][i][0][4] = -dt * two * (((tx1 * (r43 * c34 - c1345) + ty1 * (c34 - c1345) + tz1 * (c34 - c1345)) * (u[k][j][i][1] * u[k][j][i][1]) + (tx1 * (c34 - c1345) + ty1 * (r43 * c34 - c1345) + tz1 * (c34 - c1345)) * (u[k][j][i][2] * u[k][j][i][2]) + (tx1 * (c34 - c1345) + ty1 * (c34 - c1345) + tz1 * (r43 * c34 - c1345)) * (u[k][j][i][3] * u[k][j][i][3])) * tmp3 + (tx1 + ty1 + tz1) * c1345 * tmp2 * u[k][j][i][4]);
                        d[j][i][1][4] = dt * two * (tx1 * (r43 * c34 - c1345) + ty1 * (c34 - c1345) + tz1 * (c34 - c1345)) * tmp2 * u[k][j][i][1];
                        d[j][i][2][4] = dt * two * (tx1 * (c34 - c1345) + ty1 * (r43 * c34 - c1345) + tz1 * (c34 - c1345)) * tmp2 * u[k][j][i][2];
                        d[j][i][3][4] = dt * two * (tx1 * (c34 - c1345) + ty1 * (c34 - c1345) + tz1 * (r43 * c34 - c1345)) * tmp2 * u[k][j][i][3];
                        d[j][i][4][4] = one + dt * two * (tx1 + ty1 + tz1) * c1345 * tmp1 + dt * two * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);

                        //---------------------------------------------------------------------
                        // form the first block sub-diagonal
                        //---------------------------------------------------------------------
                        tmp1 = rho_i[k][j][i + 1];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        a[j][i][0][0] = -dt * tx1 * dx1;
                        a[j][i][1][0] = dt * tx2;
                        a[j][i][2][0] = zero;
                        a[j][i][3][0] = zero;
                        a[j][i][4][0] = zero;

                        a[j][i][0][1] = dt * tx2 * (-(u[k][j][i + 1][1] * tmp1) * (u[k][j][i + 1][1] * tmp1) + cc2 * qs[k][j][i + 1] * tmp1) - dt * tx1 * (-r43 * c34 * tmp2 * u[k][j][i + 1][1]);
                        a[j][i][1][1] = dt * tx2 * ((two - cc2) * (u[k][j][i + 1][1] * tmp1)) - dt * tx1 * (r43 * c34 * tmp1) - dt * tx1 * dx2;
                        a[j][i][2][1] = dt * tx2 * (-cc2 * (u[k][j][i + 1][2] * tmp1));
                        a[j][i][3][1] = dt * tx2 * (-cc2 * (u[k][j][i + 1][3] * tmp1));
                        a[j][i][4][1] = dt * tx2 * cc2;

                        a[j][i][0][2] = dt * tx2 * (-(u[k][j][i + 1][1] * u[k][j][i + 1][2]) * tmp2) - dt * tx1 * (-c34 * tmp2 * u[k][j][i + 1][2]);
                        a[j][i][1][2] = dt * tx2 * (u[k][j][i + 1][2] * tmp1);
                        a[j][i][2][2] = dt * tx2 * (u[k][j][i + 1][1] * tmp1) - dt * tx1 * (c34 * tmp1) - dt * tx1 * dx3;
                        a[j][i][3][2] = zero;
                        a[j][i][4][2] = zero;

                        a[j][i][0][3] = dt * tx2 * (-(u[k][j][i + 1][1] * u[k][j][i + 1][3]) * tmp2) - dt * tx1 * (-c34 * tmp2 * u[k][j][i + 1][3]);
                        a[j][i][1][3] = dt * tx2 * (u[k][j][i + 1][3] * tmp1);
                        a[j][i][2][3] = zero;
                        a[j][i][3][3] = dt * tx2 * (u[k][j][i + 1][1] * tmp1) - dt * tx1 * (c34 * tmp1) - dt * tx1 * dx4;
                        a[j][i][4][3] = zero;

                        a[j][i][0][4] = dt * tx2 * ((cc2 * two * qs[k][j][i + 1] - cc1 * u[k][j][i + 1][4]) * (u[k][j][i + 1][1] * tmp2)) - dt * tx1 * (-(r43 * c34 - c1345) * tmp3 * (u[k][j][i + 1][1] * u[k][j][i + 1][1]) - (c34 - c1345) * tmp3 * (u[k][j][i + 1][2] * u[k][j][i + 1][2]) - (c34 - c1345) * tmp3 * (u[k][j][i + 1][3] * u[k][j][i + 1][3]) - c1345 * tmp2 * u[k][j][i + 1][4]);
                        a[j][i][1][4] = dt * tx2 * (cc1 * (u[k][j][i + 1][4] * tmp1) - cc2 * (u[k][j][i + 1][1] * u[k][j][i + 1][1] * tmp2 + qs[k][j][i + 1] * tmp1)) - dt * tx1 * (r43 * c34 - c1345) * tmp2 * u[k][j][i + 1][1];
                        a[j][i][2][4] = dt * tx2 * (-cc2 * (u[k][j][i + 1][2] * u[k][j][i + 1][1]) * tmp2) - dt * tx1 * (c34 - c1345) * tmp2 * u[k][j][i + 1][2];
                        a[j][i][3][4] = dt * tx2 * (-cc2 * (u[k][j][i + 1][3] * u[k][j][i + 1][1]) * tmp2) - dt * tx1 * (c34 - c1345) * tmp2 * u[k][j][i + 1][3];
                        a[j][i][4][4] = dt * tx2 * (cc1 * (u[k][j][i + 1][1] * tmp1)) - dt * tx1 * c1345 * tmp1 - dt * tx1 * dx5;

                        //---------------------------------------------------------------------
                        // form the second block sub-diagonal
                        //---------------------------------------------------------------------
                        tmp1 = rho_i[k][j + 1][i];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        b[j][i][0][0] = -dt * ty1 * dy1;
                        b[j][i][1][0] = zero;
                        b[j][i][2][0] = dt * ty2;
                        b[j][i][3][0] = zero;
                        b[j][i][4][0] = zero;

                        b[j][i][0][1] = dt * ty2 * (-(u[k][j + 1][i][1] * u[k][j + 1][i][2]) * tmp2) - dt * ty1 * (-c34 * tmp2 * u[k][j + 1][i][1]);
                        b[j][i][1][1] = dt * ty2 * (u[k][j + 1][i][2] * tmp1) - dt * ty1 * (c34 * tmp1) - dt * ty1 * dy2;
                        b[j][i][2][1] = dt * ty2 * (u[k][j + 1][i][1] * tmp1);
                        b[j][i][3][1] = zero;
                        b[j][i][4][1] = zero;

                        b[j][i][0][2] = dt * ty2 * (-(u[k][j + 1][i][2] * tmp1) * (u[k][j + 1][i][2] * tmp1) + cc2 * (qs[k][j + 1][i] * tmp1)) - dt * ty1 * (-r43 * c34 * tmp2 * u[k][j + 1][i][2]);
                        b[j][i][1][2] = dt * ty2 * (-cc2 * (u[k][j + 1][i][1] * tmp1));
                        b[j][i][2][2] = dt * ty2 * ((two - cc2) * (u[k][j + 1][i][2] * tmp1)) - dt * ty1 * (r43 * c34 * tmp1) - dt * ty1 * dy3;
                        b[j][i][3][2] = dt * ty2 * (-cc2 * (u[k][j + 1][i][3] * tmp1));
                        b[j][i][4][2] = dt * ty2 * cc2;

                        b[j][i][0][3] = dt * ty2 * (-(u[k][j + 1][i][2] * u[k][j + 1][i][3]) * tmp2) - dt * ty1 * (-c34 * tmp2 * u[k][j + 1][i][3]);
                        b[j][i][1][3] = zero;
                        b[j][i][2][3] = dt * ty2 * (u[k][j + 1][i][3] * tmp1);
                        b[j][i][3][3] = dt * ty2 * (u[k][j + 1][i][2] * tmp1) - dt * ty1 * (c34 * tmp1) - dt * ty1 * dy4;
                        b[j][i][4][3] = zero;

                        b[j][i][0][4] = dt * ty2 * ((cc2 * two * qs[k][j + 1][i] - cc1 * u[k][j + 1][i][4]) * (u[k][j + 1][i][2] * tmp2)) - dt * ty1 * (-(c34 - c1345) * tmp3 * (u[k][j + 1][i][1] * u[k][j + 1][i][1]) - (r43 * c34 - c1345) * tmp3 * (u[k][j + 1][i][2] * u[k][j + 1][i][2]) - (c34 - c1345) * tmp3 * (u[k][j + 1][i][3] * u[k][j + 1][i][3]) - c1345 * tmp2 * u[k][j + 1][i][4]);
                        b[j][i][1][4] = dt * ty2 * (-cc2 * (u[k][j + 1][i][1] * u[k][j + 1][i][2]) * tmp2) - dt * ty1 * (c34 - c1345) * tmp2 * u[k][j + 1][i][1];
                        b[j][i][2][4] = dt * ty2 * (cc1 * (u[k][j + 1][i][4] * tmp1) - cc2 * (qs[k][j + 1][i] * tmp1 + u[k][j + 1][i][2] * u[k][j + 1][i][2] * tmp2)) - dt * ty1 * (r43 * c34 - c1345) * tmp2 * u[k][j + 1][i][2];
                        b[j][i][3][4] = dt * ty2 * (-cc2 * (u[k][j + 1][i][2] * u[k][j + 1][i][3]) * tmp2) - dt * ty1 * (c34 - c1345) * tmp2 * u[k][j + 1][i][3];
                        b[j][i][4][4] = dt * ty2 * (cc1 * (u[k][j + 1][i][2] * tmp1)) - dt * ty1 * c1345 * tmp1 - dt * ty1 * dy5;

                        //---------------------------------------------------------------------
                        // form the third block sub-diagonal
                        //---------------------------------------------------------------------
                        tmp1 = rho_i[k + 1][j][i];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        c[j][i][0][0] = -dt * tz1 * dz1;
                        c[j][i][1][0] = zero;
                        c[j][i][2][0] = zero;
                        c[j][i][3][0] = dt * tz2;
                        c[j][i][4][0] = zero;

                        c[j][i][0][1] = dt * tz2 * (-(u[k + 1][j][i][1] * u[k + 1][j][i][3]) * tmp2) - dt * tz1 * (-c34 * tmp2 * u[k + 1][j][i][1]);
                        c[j][i][1][1] = dt * tz2 * (u[k + 1][j][i][3] * tmp1) - dt * tz1 * c34 * tmp1 - dt * tz1 * dz2;
                        c[j][i][2][1] = zero;
                        c[j][i][3][1] = dt * tz2 * (u[k + 1][j][i][1] * tmp1);
                        c[j][i][4][1] = zero;

                        c[j][i][0][2] = dt * tz2 * (-(u[k + 1][j][i][2] * u[k + 1][j][i][3]) * tmp2) - dt * tz1 * (-c34 * tmp2 * u[k + 1][j][i][2]);
                        c[j][i][1][2] = zero;
                        c[j][i][2][2] = dt * tz2 * (u[k + 1][j][i][3] * tmp1) - dt * tz1 * (c34 * tmp1) - dt * tz1 * dz3;
                        c[j][i][3][2] = dt * tz2 * (u[k + 1][j][i][2] * tmp1);
                        c[j][i][4][2] = zero;

                        c[j][i][0][3] = dt * tz2 * (-(u[k + 1][j][i][3] * tmp1) * (u[k + 1][j][i][3] * tmp1) + cc2 * (qs[k + 1][j][i] * tmp1)) - dt * tz1 * (-r43 * c34 * tmp2 * u[k + 1][j][i][3]);
                        c[j][i][1][3] = dt * tz2 * (-cc2 * (u[k + 1][j][i][1] * tmp1));
                        c[j][i][2][3] = dt * tz2 * (-cc2 * (u[k + 1][j][i][2] * tmp1));
                        c[j][i][3][3] = dt * tz2 * (two - cc2) * (u[k + 1][j][i][3] * tmp1) - dt * tz1 * (r43 * c34 * tmp1) - dt * tz1 * dz4;
                        c[j][i][4][3] = dt * tz2 * cc2;

                        c[j][i][0][4] = dt * tz2 * ((cc2 * two * qs[k + 1][j][i] - cc1 * u[k + 1][j][i][4]) * (u[k + 1][j][i][3] * tmp2)) - dt * tz1 * (-(c34 - c1345) * tmp3 * (u[k + 1][j][i][1] * u[k + 1][j][i][1]) - (c34 - c1345) * tmp3 * (u[k + 1][j][i][2] * u[k + 1][j][i][2]) - (r43 * c34 - c1345) * tmp3 * (u[k + 1][j][i][3] * u[k + 1][j][i][3]) - c1345 * tmp2 * u[k + 1][j][i][4]);
                        c[j][i][1][4] = dt * tz2 * (-cc2 * (u[k + 1][j][i][1] * u[k + 1][j][i][3]) * tmp2) - dt * tz1 * (c34 - c1345) * tmp2 * u[k + 1][j][i][1];
                        c[j][i][2][4] = dt * tz2 * (-cc2 * (u[k + 1][j][i][2] * u[k + 1][j][i][3]) * tmp2) - dt * tz1 * (c34 - c1345) * tmp2 * u[k + 1][j][i][2];
                        c[j][i][3][4] = dt * tz2 * (cc1 * (u[k + 1][j][i][4] * tmp1) - cc2 * (qs[k + 1][j][i] * tmp1 + u[k + 1][j][i][3] * u[k + 1][j][i][3] * tmp2)) - dt * tz1 * (r43 * c34 - c1345) * tmp2 * u[k + 1][j][i][3];
                        c[j][i][4][4] = dt * tz2 * (cc1 * (u[k + 1][j][i][3] * tmp1)) - dt * tz1 * c1345 * tmp1 - dt * tz1 * dz5;
                    }
                }
            }

            if (timeron)
                timer_stop(t_jacu);

            //---------------------------------------------------------------------
            // perform the upper triangular solution
            //---------------------------------------------------------------------
            if (timeron)
                timer_start(t_buts);

            /* buts( ISIZ1, ISIZ2, ISIZ3,
            nx, ny, nz, k,
            omega,
            rsd, tv,
            d, a, b, c,
            ist, iend, jst, jend,
            nx0, ny0 );
      */
            /*void buts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, int k,
                double omega, double v[][ldmy / 2 * 2 + 1][ldmx / 2 * 2 + 1][5],
                double tv[ldmy][ldmx / 2 * 2 + 1][5],
                double d[ldmy][ldmx / 2 * 2 + 1][5][5],
                double udx[ldmy][ldmx / 2 * 2 + 1][5][5],
                double udy[ldmy][ldmx / 2 * 2 + 1][5][5],
                double udz[ldmy][ldmx / 2 * 2 + 1][5][5], int ist, int iend,
                int jst, int jend, int nx0, int ny0)
      */
            {
                //---------------------------------------------------------------------
                // local variables
                //---------------------------------------------------------------------
                // int i, j, m;
                // double tmp, tmp1;
                // double tmat[5][5];

                // #pragma scop
                for (j = jend - 1; j >= jst; j--) {
                    for (i = iend - 1; i >= ist; i--) {
                        for (m = 0; m < 5; m++) {
                            tv[j][i][m] = omega * (c[j][i][0][m] * rsd[k + 1][j][i][0] + c[j][i][1][m] * rsd[k + 1][j][i][1] + c[j][i][2][m] * rsd[k + 1][j][i][2] + c[j][i][3][m] * rsd[k + 1][j][i][3] + c[j][i][4][m] * rsd[k + 1][j][i][4]);
                        }
                    }
                }

                for (j = jend - 1; j >= jst; j--) {
                    for (i = iend - 1; i >= ist; i--) {
                        for (m = 0; m < 5; m++) {
                            tv[j][i][m] =
                                tv[j][i][m] + omega * (b[j][i][0][m] * rsd[k][j + 1][i][0] + a[j][i][0][m] * rsd[k][j][i + 1][0] + b[j][i][1][m] * rsd[k][j + 1][i][1] + a[j][i][1][m] * rsd[k][j][i + 1][1] + b[j][i][2][m] * rsd[k][j + 1][i][2] + a[j][i][2][m] * rsd[k][j][i + 1][2] + b[j][i][3][m] * rsd[k][j + 1][i][3] + a[j][i][3][m] * rsd[k][j][i + 1][3] + b[j][i][4][m] * rsd[k][j + 1][i][4] + a[j][i][4][m] * rsd[k][j][i + 1][4]);
                        }

                        //---------------------------------------------------------------------
                        // diagonal block inversion
                        //---------------------------------------------------------------------
                        for (m = 0; m < 5; m++) {
                            tmat[m][0] = d[j][i][0][m];
                            tmat[m][1] = d[j][i][1][m];
                            tmat[m][2] = d[j][i][2][m];
                            tmat[m][3] = d[j][i][3][m];
                            tmat[m][4] = d[j][i][4][m];
                        }

                        tmp1 = one / tmat[0][0];
                        tmp = tmp1 * tmat[1][0];
                        tmat[1][1] = tmat[1][1] - tmp * tmat[0][1];
                        tmat[1][2] = tmat[1][2] - tmp * tmat[0][2];
                        tmat[1][3] = tmat[1][3] - tmp * tmat[0][3];
                        tmat[1][4] = tmat[1][4] - tmp * tmat[0][4];
                        tv[j][i][1] = tv[j][i][1] - tv[j][i][0] * tmp;

                        tmp = tmp1 * tmat[2][0];
                        tmat[2][1] = tmat[2][1] - tmp * tmat[0][1];
                        tmat[2][2] = tmat[2][2] - tmp * tmat[0][2];
                        tmat[2][3] = tmat[2][3] - tmp * tmat[0][3];
                        tmat[2][4] = tmat[2][4] - tmp * tmat[0][4];
                        tv[j][i][2] = tv[j][i][2] - tv[j][i][0] * tmp;

                        tmp = tmp1 * tmat[3][0];
                        tmat[3][1] = tmat[3][1] - tmp * tmat[0][1];
                        tmat[3][2] = tmat[3][2] - tmp * tmat[0][2];
                        tmat[3][3] = tmat[3][3] - tmp * tmat[0][3];
                        tmat[3][4] = tmat[3][4] - tmp * tmat[0][4];
                        tv[j][i][3] = tv[j][i][3] - tv[j][i][0] * tmp;

                        tmp = tmp1 * tmat[4][0];
                        tmat[4][1] = tmat[4][1] - tmp * tmat[0][1];
                        tmat[4][2] = tmat[4][2] - tmp * tmat[0][2];
                        tmat[4][3] = tmat[4][3] - tmp * tmat[0][3];
                        tmat[4][4] = tmat[4][4] - tmp * tmat[0][4];
                        tv[j][i][4] = tv[j][i][4] - tv[j][i][0] * tmp;

                        tmp1 = one / tmat[1][1];
                        tmp = tmp1 * tmat[2][1];
                        tmat[2][2] = tmat[2][2] - tmp * tmat[1][2];
                        tmat[2][3] = tmat[2][3] - tmp * tmat[1][3];
                        tmat[2][4] = tmat[2][4] - tmp * tmat[1][4];
                        tv[j][i][2] = tv[j][i][2] - tv[j][i][1] * tmp;

                        tmp = tmp1 * tmat[3][1];
                        tmat[3][2] = tmat[3][2] - tmp * tmat[1][2];
                        tmat[3][3] = tmat[3][3] - tmp * tmat[1][3];
                        tmat[3][4] = tmat[3][4] - tmp * tmat[1][4];
                        tv[j][i][3] = tv[j][i][3] - tv[j][i][1] * tmp;

                        tmp = tmp1 * tmat[4][1];
                        tmat[4][2] = tmat[4][2] - tmp * tmat[1][2];
                        tmat[4][3] = tmat[4][3] - tmp * tmat[1][3];
                        tmat[4][4] = tmat[4][4] - tmp * tmat[1][4];
                        tv[j][i][4] = tv[j][i][4] - tv[j][i][1] * tmp;

                        tmp1 = one / tmat[2][2];
                        tmp = tmp1 * tmat[3][2];
                        tmat[3][3] = tmat[3][3] - tmp * tmat[2][3];
                        tmat[3][4] = tmat[3][4] - tmp * tmat[2][4];
                        tv[j][i][3] = tv[j][i][3] - tv[j][i][2] * tmp;

                        tmp = tmp1 * tmat[4][2];
                        tmat[4][3] = tmat[4][3] - tmp * tmat[2][3];
                        tmat[4][4] = tmat[4][4] - tmp * tmat[2][4];
                        tv[j][i][4] = tv[j][i][4] - tv[j][i][2] * tmp;

                        tmp1 = one / tmat[3][3];
                        tmp = tmp1 * tmat[4][3];
                        tmat[4][4] = tmat[4][4] - tmp * tmat[3][4];
                        tv[j][i][4] = tv[j][i][4] - tv[j][i][3] * tmp;

                        //---------------------------------------------------------------------
                        // back substitution
                        //---------------------------------------------------------------------
                        tv[j][i][4] = tv[j][i][4] / tmat[4][4];

                        tv[j][i][3] = tv[j][i][3] - tmat[3][4] * tv[j][i][4];
                        tv[j][i][3] = tv[j][i][3] / tmat[3][3];

                        tv[j][i][2] = tv[j][i][2] - tmat[2][3] * tv[j][i][3] - tmat[2][4] * tv[j][i][4];
                        tv[j][i][2] = tv[j][i][2] / tmat[2][2];

                        tv[j][i][1] = tv[j][i][1] - tmat[1][2] * tv[j][i][2] - tmat[1][3] * tv[j][i][3] - tmat[1][4] * tv[j][i][4];
                        tv[j][i][1] = tv[j][i][1] / tmat[1][1];

                        tv[j][i][0] = tv[j][i][0] - tmat[0][1] * tv[j][i][1] - tmat[0][2] * tv[j][i][2] - tmat[0][3] * tv[j][i][3] - tmat[0][4] * tv[j][i][4];
                        tv[j][i][0] = tv[j][i][0] / tmat[0][0];

                        rsd[k][j][i][0] = rsd[k][j][i][0] - tv[j][i][0];
                        rsd[k][j][i][1] = rsd[k][j][i][1] - tv[j][i][1];
                        rsd[k][j][i][2] = rsd[k][j][i][2] - tv[j][i][2];
                        rsd[k][j][i][3] = rsd[k][j][i][3] - tv[j][i][3];
                        rsd[k][j][i][4] = rsd[k][j][i][4] - tv[j][i][4];
                    }
                }
                // #pragma endscop
            }

            if (timeron)
                timer_stop(t_buts);
        }

        //---------------------------------------------------------------------
        // update the variables
        //---------------------------------------------------------------------
        if (timeron)
            timer_start(t_add);

        // #pragma scop
        for (k = 1; k < nz - 1; k++) {
            for (j = jst; j < jend; j++) {
                for (i = ist; i < iend; i++) {
                    for (m = 0; m < 5; m++) {
                        u[k][j][i][m] = u[k][j][i][m] + tmp_ssor * rsd[k][j][i][m];
                    }
                }
            }
        }
        // #pragma endscop

        if (timeron)
            timer_stop(t_add);

        //---------------------------------------------------------------------
        // compute the max-norms of newton iteration corrections
        //---------------------------------------------------------------------
        if ((istep % inorm) == 0) {
            if (timeron) timer_start(t_l2norm);

            l2norm(ISIZ1, ISIZ2, ISIZ3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, delunm);
            /*void l2norm(int ldx, int ldy, int ldz, int nx0, int ny0, int nz0, int
         ist, int iend, int jst, int jend, double v[][ldy / 2 * 2 + 1][ldx / 2 *
         2 + 1][5], double sum[5]) */
            /*
            {
                //---------------------------------------------------------------------
                // local variables
                //---------------------------------------------------------------------
                // int i, j, k, m;

                for (m = 0; m < 5; m++) {
                delunm[m] = zero;
                }

                for (k = 1; k < nz0 - 1; k++) {
                for (j = jst; j < jend; j++) {
                    for (i = ist; i < iend; i++) {
                    for (m = 0; m < 5; m++) {
                        delunm[m] = delunm[m] + rsd[k][j][i][m] * rsd[k][j][i][m];
                    }
                    }
                }
                }

                for (m = 0; m < 5; m++) {
                delunm[m] = sqrt(delunm[m] / ((nx0 - 2) * (ny0 - 2) * (nz0 - 2)));
                }
            }
            */

            if (timeron) timer_stop(t_l2norm);
            /*
            if ( ipr == 1 ) {
                printf(" \n RMS-norm of SSOR-iteration correction "
                    "for first pde  = %12.5E\n"
                    " RMS-norm of SSOR-iteration correction "
                    "for second pde = %12.5E\n"
                    " RMS-norm of SSOR-iteration correction "
                    "for third pde  = %12.5E\n"
                    " RMS-norm of SSOR-iteration correction "
                    "for fourth pde = %12.5E\n",
                    " RMS-norm of SSOR-iteration correction "
                    "for fifth pde  = %12.5E\n",
                    delunm[0], delunm[1], delunm[2], delunm[3], delunm[4]);
            } else if ( ipr == 2 ) {
                printf("(%5d,%15.6f)\n", istep, delunm[4]);
            }
            */
        }

        //---------------------------------------------------------------------
        // compute the steady-state residuals
        //---------------------------------------------------------------------
        rhs();
        // void rhs()
        /*
        {
            //---------------------------------------------------------------------
            // local variables
            //---------------------------------------------------------------------
            // int i, j, k, m;
            double q;
            double utmp[ISIZ3][6], rtmp[ISIZ3][5];
            double u21, u31, u41;
            double u21i, u31i, u41i, u51i;
            double u21j, u31j, u41j, u51j;
            double u21k, u31k, u41k, u51k;
            double u21im1, u31im1, u41im1, u51im1;
            double u21jm1, u31jm1, u41jm1, u51jm1;
            double u21km1, u31km1, u41km1, u51km1;

            if (timeron)
                timer_start(t_rhs);
            for (k = 0; k < nz; k++) {
                for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    for (m = 0; m < 5; m++) {
                    rsd[k][j][i][m] = -frct[k][j][i][m];
                    }
                    tmp = one / u[k][j][i][0];
                    rho_i[k][j][i] = tmp;
                    qs[k][j][i] =
                        0.50 *
                        (u[k][j][i][1] * u[k][j][i][1] + u[k][j][i][2] * u[k][j][i][2] +
                        u[k][j][i][3] * u[k][j][i][3]) *
                        tmp;
                }
                }
            }

            if (timeron)
                timer_start(t_rhsx);
            //---------------------------------------------------------------------
            // xi-direction flux differences
            //---------------------------------------------------------------------
            for (k = 1; k < nz - 1; k++) {
                for (j = jst; j < jend; j++) {
                for (i = 0; i < nx; i++) {
                    flux[i][0] = u[k][j][i][1];
                    u21 = u[k][j][i][1] * rho_i[k][j][i];
                    q = qs[k][j][i];
                    flux[i][1] = u[k][j][i][1] * u21 + cc2 * (u[k][j][i][4] - q);
                    flux[i][2] = u[k][j][i][2] * u21;
                    flux[i][3] = u[k][j][i][3] * u21;
                    flux[i][4] = (cc1 * u[k][j][i][4] - cc2 * q) * u21;
                }

                for (i = ist; i < iend; i++) {
                    for (m = 0; m < 5; m++) {
                    rsd[k][j][i][m] =
                        rsd[k][j][i][m] - tx2 * (flux[i + 1][m] - flux[i - 1][m]);
                    }
                }

                for (i = ist; i < nx; i++) {
                    tmp = rho_i[k][j][i];

                    u21i = tmp * u[k][j][i][1];
                    u31i = tmp * u[k][j][i][2];
                    u41i = tmp * u[k][j][i][3];
                    u51i = tmp * u[k][j][i][4];

                    tmp = rho_i[k][j][i - 1];

                    u21im1 = tmp * u[k][j][i - 1][1];
                    u31im1 = tmp * u[k][j][i - 1][2];
                    u41im1 = tmp * u[k][j][i - 1][3];
                    u51im1 = tmp * u[k][j][i - 1][4];

                    flux[i][1] = (4.0 / 3.0) * tx3 * (u21i - u21im1);
                    flux[i][2] = tx3 * (u31i - u31im1);
                    flux[i][3] = tx3 * (u41i - u41im1);
                    flux[i][4] =
                        0.50 * (one - cc1 * cc5) * tx3 *
                            ((u21i * u21i + u31i * u31i + u41i * u41i) -
                            (u21im1 * u21im1 + u31im1 * u31im1 + u41im1 * u41im1)) +
                        (one / 6.0) * tx3 * (u21i * u21i - u21im1 * u21im1) +
                        cc1 * cc5 * tx3 * (u51i - u51im1);
                }

                for (i = ist; i < iend; i++) {
                    rsd[k][j][i][0] =
                        rsd[k][j][i][0] + dx1 * tx1 *
                                            (u[k][j][i - 1][0] - two * u[k][j][i][0] +
                                            u[k][j][i + 1][0]);
                    rsd[k][j][i][1] = rsd[k][j][i][1] +
                                    tx3 * cc3 * cc4 * (flux[i + 1][1] - flux[i][1]) +
                                    dx2 * tx1 *
                                        (u[k][j][i - 1][1] - two * u[k][j][i][1] +
                                        u[k][j][i + 1][1]);
                    rsd[k][j][i][2] = rsd[k][j][i][2] +
                                    tx3 * cc3 * cc4 * (flux[i + 1][2] - flux[i][2]) +
                                    dx3 * tx1 *
                                        (u[k][j][i - 1][2] - two * u[k][j][i][2] +
                                        u[k][j][i + 1][2]);
                    rsd[k][j][i][3] = rsd[k][j][i][3] +
                                    tx3 * cc3 * cc4 * (flux[i + 1][3] - flux[i][3]) +
                                    dx4 * tx1 *
                                        (u[k][j][i - 1][3] - two * u[k][j][i][3] +
                                        u[k][j][i + 1][3]);
                    rsd[k][j][i][4] = rsd[k][j][i][4] +
                                    tx3 * cc3 * cc4 * (flux[i + 1][4] - flux[i][4]) +
                                    dx5 * tx1 *
                                        (u[k][j][i - 1][4] - two * u[k][j][i][4] +
                                        u[k][j][i + 1][4]);
                }

                //---------------------------------------------------------------------
                // Fourth-order dissipation
                //---------------------------------------------------------------------
                for (m = 0; m < 5; m++) {
                    rsd[k][j][1][m] =
                        rsd[k][j][1][m] - dssp * (5.0 * u[k][j][1][m] -
                                                4.0 * u[k][j][2][m] + u[k][j][3][m]);
                    rsd[k][j][2][m] =
                        rsd[k][j][2][m] -
                        dssp * (-4.0 * u[k][j][1][m] + 6.0 * u[k][j][2][m] -
                                4.0 * u[k][j][3][m] + u[k][j][4][m]);
                }

                for (i = 3; i < nx - 3; i++) {
                    for (m = 0; m < 5; m++) {
                    rsd[k][j][i][m] =
                        rsd[k][j][i][m] -
                        dssp * (u[k][j][i - 2][m] - 4.0 * u[k][j][i - 1][m] +
                                6.0 * u[k][j][i][m] - 4.0 * u[k][j][i + 1][m] +
                                u[k][j][i + 2][m]);
                    }
                }

                for (m = 0; m < 5; m++) {
                    rsd[k][j][nx - 3][m] =
                        rsd[k][j][nx - 3][m] -
                        dssp * (u[k][j][nx - 5][m] - 4.0 * u[k][j][nx - 4][m] +
                                6.0 * u[k][j][nx - 3][m] - 4.0 * u[k][j][nx - 2][m]);
                    rsd[k][j][nx - 2][m] =
                        rsd[k][j][nx - 2][m] -
                        dssp * (u[k][j][nx - 4][m] - 4.0 * u[k][j][nx - 3][m] +
                                5.0 * u[k][j][nx - 2][m]);
                }
                }
            }
            if (timeron)
                timer_stop(t_rhsx);

            if (timeron)
                timer_start(t_rhsy);
            //---------------------------------------------------------------------
            // eta-direction flux differences
            //---------------------------------------------------------------------
            for (k = 1; k < nz - 1; k++) {
                for (i = ist; i < iend; i++) {
                for (j = 0; j < ny; j++) {
                    flux[j][0] = u[k][j][i][2];
                    u31 = u[k][j][i][2] * rho_i[k][j][i];

                    q = qs[k][j][i];

                    flux[j][1] = u[k][j][i][1] * u31;
                    flux[j][2] = u[k][j][i][2] * u31 + cc2 * (u[k][j][i][4] - q);
                    flux[j][3] = u[k][j][i][3] * u31;
                    flux[j][4] = (cc1 * u[k][j][i][4] - cc2 * q) * u31;
                }

                for (j = jst; j < jend; j++) {
                    for (m = 0; m < 5; m++) {
                    rsd[k][j][i][m] =
                        rsd[k][j][i][m] - ty2 * (flux[j + 1][m] - flux[j - 1][m]);
                    }
                }

                for (j = jst; j < ny; j++) {
                    tmp = rho_i[k][j][i];

                    u21j = tmp * u[k][j][i][1];
                    u31j = tmp * u[k][j][i][2];
                    u41j = tmp * u[k][j][i][3];
                    u51j = tmp * u[k][j][i][4];

                    tmp = rho_i[k][j - 1][i];
                    u21jm1 = tmp * u[k][j - 1][i][1];
                    u31jm1 = tmp * u[k][j - 1][i][2];
                    u41jm1 = tmp * u[k][j - 1][i][3];
                    u51jm1 = tmp * u[k][j - 1][i][4];

                    flux[j][1] = ty3 * (u21j - u21jm1);
                    flux[j][2] = (4.0 / 3.0) * ty3 * (u31j - u31jm1);
                    flux[j][3] = ty3 * (u41j - u41jm1);
                    flux[j][4] =
                        0.50 * (one - cc1 * cc5) * ty3 *
                            ((u21j * u21j + u31j * u31j + u41j * u41j) -
                            (u21jm1 * u21jm1 + u31jm1 * u31jm1 + u41jm1 * u41jm1)) +
                        (one / 6.0) * ty3 * (u31j * u31j - u31jm1 * u31jm1) +
                        cc1 * cc5 * ty3 * (u51j - u51jm1);
                }

                for (j = jst; j < jend; j++) {
                    rsd[k][j][i][0] =
                        rsd[k][j][i][0] + dy1 * ty1 *
                                            (u[k][j - 1][i][0] - two * u[k][j][i][0] +
                                            u[k][j + 1][i][0]);

                    rsd[k][j][i][1] = rsd[k][j][i][1] +
                                    ty3 * cc3 * cc4 * (flux[j + 1][1] - flux[j][1]) +
                                    dy2 * ty1 *
                                        (u[k][j - 1][i][1] - two * u[k][j][i][1] +
                                        u[k][j + 1][i][1]);

                    rsd[k][j][i][2] = rsd[k][j][i][2] +
                                    ty3 * cc3 * cc4 * (flux[j + 1][2] - flux[j][2]) +
                                    dy3 * ty1 *
                                        (u[k][j - 1][i][2] - two * u[k][j][i][2] +
                                        u[k][j + 1][i][2]);

                    rsd[k][j][i][3] = rsd[k][j][i][3] +
                                    ty3 * cc3 * cc4 * (flux[j + 1][3] - flux[j][3]) +
                                    dy4 * ty1 *
                                        (u[k][j - 1][i][3] - two * u[k][j][i][3] +
                                        u[k][j + 1][i][3]);

                    rsd[k][j][i][4] = rsd[k][j][i][4] +
                                    ty3 * cc3 * cc4 * (flux[j + 1][4] - flux[j][4]) +
                                    dy5 * ty1 *
                                        (u[k][j - 1][i][4] - two * u[k][j][i][4] +
                                        u[k][j + 1][i][4]);
                }
                }

                //---------------------------------------------------------------------
                // fourth-order dissipation
                //---------------------------------------------------------------------
                for (i = ist; i < iend; i++) {
                for (m = 0; m < 5; m++) {
                    rsd[k][1][i][m] =
                        rsd[k][1][i][m] - dssp * (5.0 * u[k][1][i][m] -
                                                4.0 * u[k][2][i][m] + u[k][3][i][m]);
                    rsd[k][2][i][m] =
                        rsd[k][2][i][m] -
                        dssp * (-4.0 * u[k][1][i][m] + 6.0 * u[k][2][i][m] -
                                4.0 * u[k][3][i][m] + u[k][4][i][m]);
                }
                }

                for (j = 3; j < ny - 3; j++) {
                for (i = ist; i < iend; i++) {
                    for (m = 0; m < 5; m++) {
                    rsd[k][j][i][m] =
                        rsd[k][j][i][m] -
                        dssp * (u[k][j - 2][i][m] - 4.0 * u[k][j - 1][i][m] +
                                6.0 * u[k][j][i][m] - 4.0 * u[k][j + 1][i][m] +
                                u[k][j + 2][i][m]);
                    }
                }
                }

                for (i = ist; i < iend; i++) {
                for (m = 0; m < 5; m++) {
                    rsd[k][ny - 3][i][m] =
                        rsd[k][ny - 3][i][m] -
                        dssp * (u[k][ny - 5][i][m] - 4.0 * u[k][ny - 4][i][m] +
                                6.0 * u[k][ny - 3][i][m] - 4.0 * u[k][ny - 2][i][m]);
                    rsd[k][ny - 2][i][m] =
                        rsd[k][ny - 2][i][m] -
                        dssp * (u[k][ny - 4][i][m] - 4.0 * u[k][ny - 3][i][m] +
                                5.0 * u[k][ny - 2][i][m]);
                }
                }
            }
            if (timeron)
                timer_stop(t_rhsy);

            if (timeron)
                timer_start(t_rhsz);
            //---------------------------------------------------------------------
            // zeta-direction flux differences
            //---------------------------------------------------------------------
            for (j = jst; j < jend; j++) {
                for (i = ist; i < iend; i++) {
                for (k = 0; k < nz; k++) {
                    utmp[k][0] = u[k][j][i][0];
                    utmp[k][1] = u[k][j][i][1];
                    utmp[k][2] = u[k][j][i][2];
                    utmp[k][3] = u[k][j][i][3];
                    utmp[k][4] = u[k][j][i][4];
                    utmp[k][5] = rho_i[k][j][i];
                }
                for (k = 0; k < nz; k++) {
                    flux[k][0] = utmp[k][3];
                    u41 = utmp[k][3] * utmp[k][5];

                    q = qs[k][j][i];

                    flux[k][1] = utmp[k][1] * u41;
                    flux[k][2] = utmp[k][2] * u41;
                    flux[k][3] = utmp[k][3] * u41 + cc2 * (utmp[k][4] - q);
                    flux[k][4] = (cc1 * utmp[k][4] - cc2 * q) * u41;
                }

                for (k = 1; k < nz - 1; k++) {
                    for (m = 0; m < 5; m++) {
                    rtmp[k][m] =
                        rsd[k][j][i][m] - tz2 * (flux[k + 1][m] - flux[k - 1][m]);
                    }
                }

                for (k = 1; k < nz; k++) {
                    tmp = utmp[k][5];

                    u21k = tmp * utmp[k][1];
                    u31k = tmp * utmp[k][2];
                    u41k = tmp * utmp[k][3];
                    u51k = tmp * utmp[k][4];

                    tmp = utmp[k - 1][5];

                    u21km1 = tmp * utmp[k - 1][1];
                    u31km1 = tmp * utmp[k - 1][2];
                    u41km1 = tmp * utmp[k - 1][3];
                    u51km1 = tmp * utmp[k - 1][4];

                    flux[k][1] = tz3 * (u21k - u21km1);
                    flux[k][2] = tz3 * (u31k - u31km1);
                    flux[k][3] = (4.0 / 3.0) * tz3 * (u41k - u41km1);
                    flux[k][4] =
                        0.50 * (one - cc1 * cc5) * tz3 *
                            ((u21k * u21k + u31k * u31k + u41k * u41k) -
                            (u21km1 * u21km1 + u31km1 * u31km1 + u41km1 * u41km1)) +
                        (one / 6.0) * tz3 * (u41k * u41k - u41km1 * u41km1) +
                        cc1 * cc5 * tz3 * (u51k - u51km1);
                }

                for (k = 1; k < nz - 1; k++) {
                    rtmp[k][0] = rtmp[k][0] + dz1 * tz1 *
                                                (utmp[k - 1][0] - two * utmp[k][0] +
                                                utmp[k + 1][0]);
                    rtmp[k][1] =
                        rtmp[k][1] + tz3 * cc3 * cc4 * (flux[k + 1][1] - flux[k][1]) +
                        dz2 * tz1 *
                            (utmp[k - 1][1] - two * utmp[k][1] + utmp[k + 1][1]);
                    rtmp[k][2] =
                        rtmp[k][2] + tz3 * cc3 * cc4 * (flux[k + 1][2] - flux[k][2]) +
                        dz3 * tz1 *
                            (utmp[k - 1][2] - two * utmp[k][2] + utmp[k + 1][2]);
                    rtmp[k][3] =
                        rtmp[k][3] + tz3 * cc3 * cc4 * (flux[k + 1][3] - flux[k][3]) +
                        dz4 * tz1 *
                            (utmp[k - 1][3] - two * utmp[k][3] + utmp[k + 1][3]);
                    rtmp[k][4] =
                        rtmp[k][4] + tz3 * cc3 * cc4 * (flux[k + 1][4] - flux[k][4]) +
                        dz5 * tz1 *
                            (utmp[k - 1][4] - two * utmp[k][4] + utmp[k + 1][4]);
                }

                //---------------------------------------------------------------------
                // fourth-order dissipation
                //---------------------------------------------------------------------
                for (m = 0; m < 5; m++) {
                    rsd[1][j][i][m] =
                        rtmp[1][m] -
                        dssp * (5.0 * utmp[1][m] - 4.0 * utmp[2][m] + utmp[3][m]);
                    rsd[2][j][i][m] =
                        rtmp[2][m] - dssp * (-4.0 * utmp[1][m] + 6.0 * utmp[2][m] -
                                            4.0 * utmp[3][m] + utmp[4][m]);
                }

                for (k = 3; k < nz - 3; k++) {
                    for (m = 0; m < 5; m++) {
                    rsd[k][j][i][m] =
                        rtmp[k][m] - dssp * (utmp[k - 2][m] - 4.0 * utmp[k - 1][m] +
                                            6.0 * utmp[k][m] - 4.0 * utmp[k + 1][m] +
                                            utmp[k + 2][m]);
                    }
                }

                for (m = 0; m < 5; m++) {
                    rsd[nz - 3][j][i][m] =
                        rtmp[nz - 3][m] -
                        dssp * (utmp[nz - 5][m] - 4.0 * utmp[nz - 4][m] +
                                6.0 * utmp[nz - 3][m] - 4.0 * utmp[nz - 2][m]);
                    rsd[nz - 2][j][i][m] =
                        rtmp[nz - 2][m] -
                        dssp * (utmp[nz - 4][m] - 4.0 * utmp[nz - 3][m] +
                                5.0 * utmp[nz - 2][m]);
                }
                }
            }
            if (timeron)
                timer_stop(t_rhsz);
            if (timeron)
                timer_stop(t_rhs);
        }
        */

        //---------------------------------------------------------------------
        // compute the max-norms of newton iteration residuals
        //---------------------------------------------------------------------
        if (((istep % inorm) == 0) || (istep == itmax)) {
            if (timeron)
                timer_start(t_l2norm);

            l2norm(ISIZ1, ISIZ2, ISIZ3, nx0, ny0, nz0, ist, iend, jst, jend, rsd,
                   rsdnm);
            /*
            void l2norm(int ldx, int ldy, int ldz, int nx0, int ny0, int nz0, int ist,
                        int iend, int jst, int jend,
                        double v[][ldy / 2 * 2 + 1][ldx / 2 * 2 + 1][5],
                        double sum[5])
            */
            /*
            {
                //---------------------------------------------------------------------
                // local variables
                //---------------------------------------------------------------------
                // int i, j, k, m;

                for (m = 0; m < 5; m++) {
                rsdnm[m] = zero;
                }

                for (k = 1; k < nz0 - 1; k++) {
                for (j = jst; j < jend; j++) {
                    for (i = ist; i < iend; i++) {
                    for (m = 0; m < 5; m++) {
                        rsdnm[m] = rsdnm[m] + rsd[k][j][i][m] * rsd[k][j][i][m];
                    }
                    }
                }
                }

                for (m = 0; m < 5; m++) {
                rsdnm[m] = sqrt(rsdnm[m] / ((nx0 - 2) * (ny0 - 2) * (nz0 - 2)));
                }
            }
            */

            if (timeron)
                timer_stop(t_l2norm);
            //
            // if ( ipr == 1 ) {
            //   printf(" \n RMS-norm of steady-state residual for "
            //          "first pde  = %12.5E\n"
            //          " RMS-norm of steady-state residual for "
            //          "second pde = %12.5E\n"
            //          " RMS-norm of steady-state residual for "
            //          "third pde  = %12.5E\n"
            //          " RMS-norm of steady-state residual for "
            //          "fourth pde = %12.5E\n"
            //          " RMS-norm of steady-state residual for "
            //          "fifth pde  = %12.5E\n",
            //          rsdnm[0], rsdnm[1], rsdnm[2], rsdnm[3], rsdnm[4]);
            // }
            //
        }

        //---------------------------------------------------------------------
        // check the newton-iteration residuals against the tolerance levels
        //---------------------------------------------------------------------
        // if ( ( rsdnm[0] < tolrsd[0] ) && ( rsdnm[1] < tolrsd[1] ) &&
        //      ( rsdnm[2] < tolrsd[2] ) && ( rsdnm[3] < tolrsd[3] ) &&
        //      ( rsdnm[4] < tolrsd[4] ) ) {
        //   //if (ipr == 1 ) {
        //   printf(" \n convergence was achieved after %4d pseudo-time steps\n",
        //       istep);
        //   //}
        //   break;
        // }
    }
#pragma endscop

    timer_stop(1);
    maxtime = timer_read(1);
}
