#include <assert.h>
#include <stdio.h>
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

    /* ppcg generated CPU code with AMP */
    #pragma scop
    
    #define ppcg_min(x,y)    (x < y ? x : y)
    #define ppcg_max(x,y)    (x > y ? x : y)
    float amp_lower_a_0[63][65][5][5];
    float amp_lower_a_1[63][65][5][5];
    float amp_lower_a_2[63][65][5][5];
    float amp_lower_a_3[63][65][5][5];
    float amp_lower_a_4[63][65][5][5];
    float amp_lower_b_0[63][65][5][5];
    float amp_lower_b_1[63][65][5][5];
    float amp_lower_b_2[63][65][5][5];
    float amp_lower_b_3[63][65][5][5];
    float amp_lower_b_4[63][65][5][5];
    float amp_lower_c_0[63][65][5][5];
    float amp_lower_c_1[63][65][5][5];
    float amp_lower_c_2[63][65][5][5];
    float amp_lower_c_3[63][65][5][5];
    float amp_lower_c_4[63][65][5][5];
    float amp_lower_c1345;
    float amp_lower_c34;
    float amp_lower_cc1;
    float amp_lower_cc2;
    float amp_lower_cc3;
    float amp_lower_cc4;
    float amp_lower_cc5;
    float amp_lower_d_0[63][65][5][5];
    float amp_lower_d_1[63][65][5][5];
    float amp_lower_d_2[63][65][5][5];
    float amp_lower_d_3[63][65][5][5];
    float amp_lower_d_4[63][65][5][5];
    float amp_lower_dt;
    float amp_lower_dx1;
    float amp_lower_dx2;
    float amp_lower_dx3;
    float amp_lower_dx4;
    float amp_lower_dx5;
    float amp_lower_dy1;
    float amp_lower_dy2;
    float amp_lower_dy3;
    float amp_lower_dy4;
    float amp_lower_dy5;
    float amp_lower_dz1;
    float amp_lower_dz2;
    float amp_lower_dz3;
    float amp_lower_dz4;
    float amp_lower_dz5;
    float amp_lower_omega;
    float amp_lower_one;
    float amp_lower_qs[64][65][65];
    float amp_lower_r43;
    float amp_lower_rho_i[64][65][65];
    float amp_lower_tmat_0[5][5];
    float amp_lower_tmat_1[5][5];
    float amp_lower_tmat_2[5][5];
    float amp_lower_tmat_3[5][5];
    float amp_lower_tmat_4[5][5];
    float amp_lower_tmp;
    float amp_lower_tmp1;
    float amp_lower_tmp2;
    float amp_lower_tmp3;
    float amp_lower_tmp_ssor;
    float amp_lower_tv[63][64][5];
    float amp_lower_tv_tmp[5];
    float amp_lower_two;
    float amp_lower_tx1;
    float amp_lower_tx2;
    float amp_lower_ty1;
    float amp_lower_ty2;
    float amp_lower_tz1;
    float amp_lower_tz2;
    float amp_lower_u[64][65][65][5];
    float amp_lower_zero;
    {
      for (int c0 = 1; c0 <= niter - (niter + 8) / 10; c0 += 1) {
        timer_start(5);
        for (int c1 = 1; c1 <= 62; c1 += 1)
          for (int c2 = 1; c2 <= 62; c2 += 1)
            for (int c3 = 1; c3 <= 62; c3 += 1)
              for (int c4 = 0; c4 <= 4; c4 += 1)
                rsd[c1][c2][c3][c4] = (dt * rsd[c1][c2][c3][c4]);
        timer_stop(5);
        for (int c1 = 1; c1 <= 62; c1 += 1) {
          timer_start(6);
          r43 = (4.0 / 3.0);
          c1345 = (((cc1 * cc3) * cc4) * cc5);
          c34 = (cc3 * cc4);
          for (int c2 = 1; c2 <= 62; c2 += 1)
            for (int c3 = 1; c3 <= 62; c3 += 1) {
              tmp1 = rho_i[c1][c2][c3];
              tmp2 = (tmp1 * tmp1);
              tmp3 = (tmp1 * tmp2);
              d[c2][c3][0][0] = (one + ((dt * two) * (((tx1 * dx1) + (ty1 * dy1)) + (tz1 * dz1))));
              d[c2][c3][1][0] = zero;
              d[c2][c3][2][0] = zero;
              d[c2][c3][3][0] = zero;
              d[c2][c3][4][0] = zero;
              d[c2][c3][0][1] = ((((((-dt) * two) * (((tx1 * r43) + ty1) + tz1)) * c34) * tmp2) * u[c1][c2][c3][1]);
              d[c2][c3][1][1] = ((one + ((((dt * two) * c34) * tmp1) * (((tx1 * r43) + ty1) + tz1))) + ((dt * two) * (((tx1 * dx2) + (ty1 * dy2)) + (tz1 * dz2))));
              d[c2][c3][2][1] = zero;
              d[c2][c3][3][1] = zero;
              d[c2][c3][4][1] = zero;
              d[c2][c3][0][2] = ((((((-dt) * two) * ((tx1 + (ty1 * r43)) + tz1)) * c34) * tmp2) * u[c1][c2][c3][2]);
              d[c2][c3][1][2] = zero;
              d[c2][c3][2][2] = ((one + ((((dt * two) * c34) * tmp1) * ((tx1 + (ty1 * r43)) + tz1))) + ((dt * two) * (((tx1 * dx3) + (ty1 * dy3)) + (tz1 * dz3))));
              d[c2][c3][3][2] = zero;
              d[c2][c3][4][2] = zero;
              d[c2][c3][0][3] = ((((((-dt) * two) * ((tx1 + ty1) + (tz1 * r43))) * c34) * tmp2) * u[c1][c2][c3][3]);
              d[c2][c3][1][3] = zero;
              d[c2][c3][2][3] = zero;
              d[c2][c3][3][3] = ((one + ((((dt * two) * c34) * tmp1) * ((tx1 + ty1) + (tz1 * r43)))) + ((dt * two) * (((tx1 * dx4) + (ty1 * dy4)) + (tz1 * dz4))));
              d[c2][c3][4][3] = zero;
              d[c2][c3][0][4] = (((-dt) * two) * ((((((((tx1 * ((r43 * c34) - c1345)) + (ty1 * (c34 - c1345))) + (tz1 * (c34 - c1345))) * (u[c1][c2][c3][1] * u[c1][c2][c3][1])) + ((((tx1 * (c34 - c1345)) + (ty1 * ((r43 * c34) - c1345))) + (tz1 * (c34 - c1345))) * (u[c1][c2][c3][2] * u[c1][c2][c3][2]))) + ((((tx1 * (c34 - c1345)) + (ty1 * (c34 - c1345))) + (tz1 * ((r43 * c34) - c1345))) * (u[c1][c2][c3][3] * u[c1][c2][c3][3]))) * tmp3) + (((((tx1 + ty1) + tz1) * c1345) * tmp2) * u[c1][c2][c3][4])));
              d[c2][c3][1][4] = ((((dt * two) * tmp2) * u[c1][c2][c3][1]) * (((tx1 * ((r43 * c34) - c1345)) + (ty1 * (c34 - c1345))) + (tz1 * (c34 - c1345))));
              d[c2][c3][2][4] = ((((dt * two) * tmp2) * u[c1][c2][c3][2]) * (((tx1 * (c34 - c1345)) + (ty1 * ((r43 * c34) - c1345))) + (tz1 * (c34 - c1345))));
              d[c2][c3][3][4] = ((((dt * two) * tmp2) * u[c1][c2][c3][3]) * (((tx1 * (c34 - c1345)) + (ty1 * (c34 - c1345))) + (tz1 * ((r43 * c34) - c1345))));
              d[c2][c3][4][4] = ((one + ((((dt * two) * ((tx1 + ty1) + tz1)) * c1345) * tmp1)) + ((dt * two) * (((tx1 * dx5) + (ty1 * dy5)) + (tz1 * dz5))));
              tmp1 = rho_i[c1 - 1][c2][c3];
              tmp2 = (tmp1 * tmp1);
              tmp3 = (tmp1 * tmp2);
              a[c2][c3][0][0] = (((-dt) * tz1) * dz1);
              a[c2][c3][1][0] = zero;
              a[c2][c3][2][0] = zero;
              a[c2][c3][3][0] = ((-dt) * tz2);
              a[c2][c3][4][0] = zero;
              a[c2][c3][0][1] = ((((-dt) * tz2) * ((-(u[c1 - 1][c2][c3][1] * u[c1 - 1][c2][c3][3])) * tmp2)) - ((dt * tz1) * (((-c34) * tmp2) * u[c1 - 1][c2][c3][1])));
              a[c2][c3][1][1] = (((((-dt) * tz2) * (u[c1 - 1][c2][c3][3] * tmp1)) - (((dt * tz1) * c34) * tmp1)) - ((dt * tz1) * dz2));
              a[c2][c3][2][1] = zero;
              a[c2][c3][3][1] = (((-dt) * tz2) * (u[c1 - 1][c2][c3][1] * tmp1));
              a[c2][c3][4][1] = zero;
              a[c2][c3][0][2] = ((((-dt) * tz2) * ((-(u[c1 - 1][c2][c3][2] * u[c1 - 1][c2][c3][3])) * tmp2)) - ((dt * tz1) * (((-c34) * tmp2) * u[c1 - 1][c2][c3][2])));
              a[c2][c3][1][2] = zero;
              a[c2][c3][2][2] = (((((-dt) * tz2) * (u[c1 - 1][c2][c3][3] * tmp1)) - ((dt * tz1) * (c34 * tmp1))) - ((dt * tz1) * dz3));
              a[c2][c3][3][2] = (((-dt) * tz2) * (u[c1 - 1][c2][c3][2] * tmp1));
              a[c2][c3][4][2] = zero;
              a[c2][c3][0][3] = ((((-dt) * tz2) * (((-(u[c1 - 1][c2][c3][3] * tmp1)) * (u[c1 - 1][c2][c3][3] * tmp1)) + ((cc2 * qs[c1 - 1][c2][c3]) * tmp1))) - ((dt * tz1) * ((((-r43) * c34) * tmp2) * u[c1 - 1][c2][c3][3])));
              a[c2][c3][1][3] = (((-dt) * tz2) * ((-cc2) * (u[c1 - 1][c2][c3][1] * tmp1)));
              a[c2][c3][2][3] = (((-dt) * tz2) * ((-cc2) * (u[c1 - 1][c2][c3][2] * tmp1)));
              a[c2][c3][3][3] = ((((((-dt) * tz2) * (two - cc2)) * (u[c1 - 1][c2][c3][3] * tmp1)) - ((dt * tz1) * ((r43 * c34) * tmp1))) - ((dt * tz1) * dz4));
              a[c2][c3][4][3] = (((-dt) * tz2) * cc2);
              a[c2][c3][0][4] = ((((-dt) * tz2) * (((((cc2 * two) * qs[c1 - 1][c2][c3]) - (cc1 * u[c1 - 1][c2][c3][4])) * u[c1 - 1][c2][c3][3]) * tmp2)) - ((dt * tz1) * ((((((-(c34 - c1345)) * tmp3) * (u[c1 - 1][c2][c3][1] * u[c1 - 1][c2][c3][1])) - (((c34 - c1345) * tmp3) * (u[c1 - 1][c2][c3][2] * u[c1 - 1][c2][c3][2]))) - ((((r43 * c34) - c1345) * tmp3) * (u[c1 - 1][c2][c3][3] * u[c1 - 1][c2][c3][3]))) - ((c1345 * tmp2) * u[c1 - 1][c2][c3][4]))));
              a[c2][c3][1][4] = ((((-dt) * tz2) * (((-cc2) * (u[c1 - 1][c2][c3][1] * u[c1 - 1][c2][c3][3])) * tmp2)) - ((((dt * tz1) * (c34 - c1345)) * tmp2) * u[c1 - 1][c2][c3][1]));
              a[c2][c3][2][4] = ((((-dt) * tz2) * (((-cc2) * (u[c1 - 1][c2][c3][2] * u[c1 - 1][c2][c3][3])) * tmp2)) - ((((dt * tz1) * (c34 - c1345)) * tmp2) * u[c1 - 1][c2][c3][2]));
              a[c2][c3][3][4] = ((((-dt) * tz2) * ((cc1 * (u[c1 - 1][c2][c3][4] * tmp1)) - (cc2 * ((qs[c1 - 1][c2][c3] * tmp1) + ((u[c1 - 1][c2][c3][3] * u[c1 - 1][c2][c3][3]) * tmp2))))) - ((((dt * tz1) * ((r43 * c34) - c1345)) * tmp2) * u[c1 - 1][c2][c3][3]));
              a[c2][c3][4][4] = (((((-dt) * tz2) * (cc1 * (u[c1 - 1][c2][c3][3] * tmp1))) - (((dt * tz1) * c1345) * tmp1)) - ((dt * tz1) * dz5));
              tmp1 = rho_i[c1][c2 - 1][c3];
              tmp2 = (tmp1 * tmp1);
              tmp3 = (tmp1 * tmp2);
              b[c2][c3][0][0] = (((-dt) * ty1) * dy1);
              b[c2][c3][1][0] = zero;
              b[c2][c3][2][0] = ((-dt) * ty2);
              b[c2][c3][3][0] = zero;
              b[c2][c3][4][0] = zero;
              b[c2][c3][0][1] = ((((-dt) * ty2) * ((-(u[c1][c2 - 1][c3][1] * u[c1][c2 - 1][c3][2])) * tmp2)) - ((dt * ty1) * (((-c34) * tmp2) * u[c1][c2 - 1][c3][1])));
              b[c2][c3][1][1] = (((((-dt) * ty2) * (u[c1][c2 - 1][c3][2] * tmp1)) - ((dt * ty1) * (c34 * tmp1))) - ((dt * ty1) * dy2));
              b[c2][c3][2][1] = (((-dt) * ty2) * (u[c1][c2 - 1][c3][1] * tmp1));
              b[c2][c3][3][1] = zero;
              b[c2][c3][4][1] = zero;
              b[c2][c3][0][2] = ((((-dt) * ty2) * (((-(u[c1][c2 - 1][c3][2] * tmp1)) * (u[c1][c2 - 1][c3][2] * tmp1)) + (cc2 * (qs[c1][c2 - 1][c3] * tmp1)))) - ((dt * ty1) * ((((-r43) * c34) * tmp2) * u[c1][c2 - 1][c3][2])));
              b[c2][c3][1][2] = (((-dt) * ty2) * ((-cc2) * (u[c1][c2 - 1][c3][1] * tmp1)));
              b[c2][c3][2][2] = (((((-dt) * ty2) * ((two - cc2) * (u[c1][c2 - 1][c3][2] * tmp1))) - ((dt * ty1) * ((r43 * c34) * tmp1))) - ((dt * ty1) * dy3));
              b[c2][c3][3][2] = (((-dt) * ty2) * ((-cc2) * (u[c1][c2 - 1][c3][3] * tmp1)));
              b[c2][c3][4][2] = (((-dt) * ty2) * cc2);
              b[c2][c3][0][3] = ((((-dt) * ty2) * ((-(u[c1][c2 - 1][c3][2] * u[c1][c2 - 1][c3][3])) * tmp2)) - ((dt * ty1) * (((-c34) * tmp2) * u[c1][c2 - 1][c3][3])));
              b[c2][c3][1][3] = zero;
              b[c2][c3][2][3] = (((-dt) * ty2) * (u[c1][c2 - 1][c3][3] * tmp1));
              b[c2][c3][3][3] = (((((-dt) * ty2) * (u[c1][c2 - 1][c3][2] * tmp1)) - ((dt * ty1) * (c34 * tmp1))) - ((dt * ty1) * dy4));
              b[c2][c3][4][3] = zero;
              b[c2][c3][0][4] = ((((-dt) * ty2) * ((((cc2 * two) * qs[c1][c2 - 1][c3]) - (cc1 * u[c1][c2 - 1][c3][4])) * (u[c1][c2 - 1][c3][2] * tmp2))) - ((dt * ty1) * ((((((-(c34 - c1345)) * tmp3) * (u[c1][c2 - 1][c3][1] * u[c1][c2 - 1][c3][1])) - ((((r43 * c34) - c1345) * tmp3) * (u[c1][c2 - 1][c3][2] * u[c1][c2 - 1][c3][2]))) - (((c34 - c1345) * tmp3) * (u[c1][c2 - 1][c3][3] * u[c1][c2 - 1][c3][3]))) - ((c1345 * tmp2) * u[c1][c2 - 1][c3][4]))));
              b[c2][c3][1][4] = ((((-dt) * ty2) * (((-cc2) * (u[c1][c2 - 1][c3][1] * u[c1][c2 - 1][c3][2])) * tmp2)) - ((((dt * ty1) * (c34 - c1345)) * tmp2) * u[c1][c2 - 1][c3][1]));
              b[c2][c3][2][4] = ((((-dt) * ty2) * ((cc1 * (u[c1][c2 - 1][c3][4] * tmp1)) - (cc2 * ((qs[c1][c2 - 1][c3] * tmp1) + ((u[c1][c2 - 1][c3][2] * u[c1][c2 - 1][c3][2]) * tmp2))))) - ((((dt * ty1) * ((r43 * c34) - c1345)) * tmp2) * u[c1][c2 - 1][c3][2]));
              b[c2][c3][3][4] = ((((-dt) * ty2) * (((-cc2) * (u[c1][c2 - 1][c3][2] * u[c1][c2 - 1][c3][3])) * tmp2)) - ((((dt * ty1) * (c34 - c1345)) * tmp2) * u[c1][c2 - 1][c3][3]));
              b[c2][c3][4][4] = (((((-dt) * ty2) * (cc1 * (u[c1][c2 - 1][c3][2] * tmp1))) - (((dt * ty1) * c1345) * tmp1)) - ((dt * ty1) * dy5));
              tmp1 = rho_i[c1][c2][c3 - 1];
              tmp2 = (tmp1 * tmp1);
              tmp3 = (tmp1 * tmp2);
              c[c2][c3][0][0] = (((-dt) * tx1) * dx1);
              c[c2][c3][1][0] = ((-dt) * tx2);
              c[c2][c3][2][0] = zero;
              c[c2][c3][3][0] = zero;
              c[c2][c3][4][0] = zero;
              c[c2][c3][0][1] = ((((-dt) * tx2) * (((-(u[c1][c2][c3 - 1][1] * tmp1)) * (u[c1][c2][c3 - 1][1] * tmp1)) + ((cc2 * qs[c1][c2][c3 - 1]) * tmp1))) - ((dt * tx1) * ((((-r43) * c34) * tmp2) * u[c1][c2][c3 - 1][1])));
              c[c2][c3][1][1] = (((((-dt) * tx2) * ((two - cc2) * (u[c1][c2][c3 - 1][1] * tmp1))) - ((dt * tx1) * ((r43 * c34) * tmp1))) - ((dt * tx1) * dx2));
              c[c2][c3][2][1] = (((-dt) * tx2) * ((-cc2) * (u[c1][c2][c3 - 1][2] * tmp1)));
              c[c2][c3][3][1] = (((-dt) * tx2) * ((-cc2) * (u[c1][c2][c3 - 1][3] * tmp1)));
              c[c2][c3][4][1] = (((-dt) * tx2) * cc2);
              c[c2][c3][0][2] = ((((-dt) * tx2) * ((-(u[c1][c2][c3 - 1][1] * u[c1][c2][c3 - 1][2])) * tmp2)) - ((dt * tx1) * (((-c34) * tmp2) * u[c1][c2][c3 - 1][2])));
              c[c2][c3][1][2] = (((-dt) * tx2) * (u[c1][c2][c3 - 1][2] * tmp1));
              c[c2][c3][2][2] = (((((-dt) * tx2) * (u[c1][c2][c3 - 1][1] * tmp1)) - ((dt * tx1) * (c34 * tmp1))) - ((dt * tx1) * dx3));
              c[c2][c3][3][2] = zero;
              c[c2][c3][4][2] = zero;
              c[c2][c3][0][3] = ((((-dt) * tx2) * ((-(u[c1][c2][c3 - 1][1] * u[c1][c2][c3 - 1][3])) * tmp2)) - ((dt * tx1) * (((-c34) * tmp2) * u[c1][c2][c3 - 1][3])));
              c[c2][c3][1][3] = (((-dt) * tx2) * (u[c1][c2][c3 - 1][3] * tmp1));
              c[c2][c3][2][3] = zero;
              c[c2][c3][3][3] = (((((-dt) * tx2) * (u[c1][c2][c3 - 1][1] * tmp1)) - ((dt * tx1) * (c34 * tmp1))) - ((dt * tx1) * dx4));
              c[c2][c3][4][3] = zero;
              c[c2][c3][0][4] = ((((-dt) * tx2) * (((((cc2 * two) * qs[c1][c2][c3 - 1]) - (cc1 * u[c1][c2][c3 - 1][4])) * u[c1][c2][c3 - 1][1]) * tmp2)) - ((dt * tx1) * ((((((-((r43 * c34) - c1345)) * tmp3) * (u[c1][c2][c3 - 1][1] * u[c1][c2][c3 - 1][1])) - (((c34 - c1345) * tmp3) * (u[c1][c2][c3 - 1][2] * u[c1][c2][c3 - 1][2]))) - (((c34 - c1345) * tmp3) * (u[c1][c2][c3 - 1][3] * u[c1][c2][c3 - 1][3]))) - ((c1345 * tmp2) * u[c1][c2][c3 - 1][4]))));
              c[c2][c3][1][4] = ((((-dt) * tx2) * ((cc1 * (u[c1][c2][c3 - 1][4] * tmp1)) - (cc2 * (((u[c1][c2][c3 - 1][1] * u[c1][c2][c3 - 1][1]) * tmp2) + (qs[c1][c2][c3 - 1] * tmp1))))) - ((((dt * tx1) * ((r43 * c34) - c1345)) * tmp2) * u[c1][c2][c3 - 1][1]));
              c[c2][c3][2][4] = ((((-dt) * tx2) * (((-cc2) * (u[c1][c2][c3 - 1][2] * u[c1][c2][c3 - 1][1])) * tmp2)) - ((((dt * tx1) * (c34 - c1345)) * tmp2) * u[c1][c2][c3 - 1][2]));
              c[c2][c3][3][4] = ((((-dt) * tx2) * (((-cc2) * (u[c1][c2][c3 - 1][3] * u[c1][c2][c3 - 1][1])) * tmp2)) - ((((dt * tx1) * (c34 - c1345)) * tmp2) * u[c1][c2][c3 - 1][3]));
              c[c2][c3][4][4] = (((((-dt) * tx2) * (cc1 * (u[c1][c2][c3 - 1][1] * tmp1))) - (((dt * tx1) * c1345) * tmp1)) - ((dt * tx1) * dx5));
            }
          timer_stop(6);
          timer_start(7);
          for (int c2 = 1; c2 <= 62; c2 += 1)
            for (int c3 = 1; c3 <= 62; c3 += 1)
              for (int c4 = 0; c4 <= 4; c4 += 1)
                rsd[c1][c2][c3][c4] = (rsd[c1][c2][c3][c4] - (omega * (((((a[c2][c3][0][c4] * rsd[c1 - 1][c2][c3][0]) + (a[c2][c3][1][c4] * rsd[c1 - 1][c2][c3][1])) + (a[c2][c3][2][c4] * rsd[c1 - 1][c2][c3][2])) + (a[c2][c3][3][c4] * rsd[c1 - 1][c2][c3][3])) + (a[c2][c3][4][c4] * rsd[c1 - 1][c2][c3][4]))));
          for (int c2 = 1; c2 <= 62; c2 += 1)
            for (int c3 = 1; c3 <= 62; c3 += 1) {
              for (int c4 = 0; c4 <= 4; c4 += 1)
                tv_tmp[c4] = (rsd[c1][c2][c3][c4] - (omega * ((((((((((b[c2][c3][0][c4] * rsd[c1][c2 - 1][c3][0]) + (c[c2][c3][0][c4] * rsd[c1][c2][c3 - 1][0])) + (b[c2][c3][1][c4] * rsd[c1][c2 - 1][c3][1])) + (c[c2][c3][1][c4] * rsd[c1][c2][c3 - 1][1])) + (b[c2][c3][2][c4] * rsd[c1][c2 - 1][c3][2])) + (c[c2][c3][2][c4] * rsd[c1][c2][c3 - 1][2])) + (b[c2][c3][3][c4] * rsd[c1][c2 - 1][c3][3])) + (c[c2][c3][3][c4] * rsd[c1][c2][c3 - 1][3])) + (b[c2][c3][4][c4] * rsd[c1][c2 - 1][c3][4])) + (c[c2][c3][4][c4] * rsd[c1][c2][c3 - 1][4]))));
              for (int c4 = 0; c4 <= 4; c4 += 1) {
                tmat[c4][0] = d[c2][c3][0][c4];
                tmat[c4][1] = d[c2][c3][1][c4];
                tmat[c4][2] = d[c2][c3][2][c4];
                tmat[c4][3] = d[c2][c3][3][c4];
                tmat[c4][4] = d[c2][c3][4][c4];
              }
              tmp1 = (one / tmat[0][0]);
              tmp = (tmp1 * tmat[1][0]);
              tmat[1][1] = (tmat[1][1] - (tmp * tmat[0][1]));
              tmat[1][2] = (tmat[1][2] - (tmp * tmat[0][2]));
              tmat[1][3] = (tmat[1][3] - (tmp * tmat[0][3]));
              tmat[1][4] = (tmat[1][4] - (tmp * tmat[0][4]));
              tv_tmp[1] = (tv_tmp[1] - (tv_tmp[0] * tmp));
              tmp = (tmp1 * tmat[2][0]);
              tmat[2][1] = (tmat[2][1] - (tmp * tmat[0][1]));
              tmat[2][2] = (tmat[2][2] - (tmp * tmat[0][2]));
              tmat[2][3] = (tmat[2][3] - (tmp * tmat[0][3]));
              tmat[2][4] = (tmat[2][4] - (tmp * tmat[0][4]));
              tv_tmp[2] = (tv_tmp[2] - (tv_tmp[0] * tmp));
              tmp = (tmp1 * tmat[3][0]);
              tmat[3][1] = (tmat[3][1] - (tmp * tmat[0][1]));
              tmat[3][2] = (tmat[3][2] - (tmp * tmat[0][2]));
              tmat[3][3] = (tmat[3][3] - (tmp * tmat[0][3]));
              tmat[3][4] = (tmat[3][4] - (tmp * tmat[0][4]));
              tv_tmp[3] = (tv_tmp[3] - (tv_tmp[0] * tmp));
              tmp = (tmp1 * tmat[4][0]);
              tmat[4][1] = (tmat[4][1] - (tmp * tmat[0][1]));
              tmat[4][2] = (tmat[4][2] - (tmp * tmat[0][2]));
              tmat[4][3] = (tmat[4][3] - (tmp * tmat[0][3]));
              tmat[4][4] = (tmat[4][4] - (tmp * tmat[0][4]));
              tv_tmp[4] = (tv_tmp[4] - (tv_tmp[0] * tmp));
              tmp1 = (one / tmat[1][1]);
              tmp = (tmp1 * tmat[2][1]);
              tmat[2][2] = (tmat[2][2] - (tmp * tmat[1][2]));
              tmat[2][3] = (tmat[2][3] - (tmp * tmat[1][3]));
              tmat[2][4] = (tmat[2][4] - (tmp * tmat[1][4]));
              tv_tmp[2] = (tv_tmp[2] - (tv_tmp[1] * tmp));
              tmp = (tmp1 * tmat[3][1]);
              tmat[3][2] = (tmat[3][2] - (tmp * tmat[1][2]));
              tmat[3][3] = (tmat[3][3] - (tmp * tmat[1][3]));
              tmat[3][4] = (tmat[3][4] - (tmp * tmat[1][4]));
              tv_tmp[3] = (tv_tmp[3] - (tv_tmp[1] * tmp));
              tmp = (tmp1 * tmat[4][1]);
              tmat[4][2] = (tmat[4][2] - (tmp * tmat[1][2]));
              tmat[4][3] = (tmat[4][3] - (tmp * tmat[1][3]));
              tmat[4][4] = (tmat[4][4] - (tmp * tmat[1][4]));
              tv_tmp[4] = (tv_tmp[4] - (tv_tmp[1] * tmp));
              tmp1 = (one / tmat[2][2]);
              tmp = (tmp1 * tmat[3][2]);
              tmat[3][3] = (tmat[3][3] - (tmp * tmat[2][3]));
              tmat[3][4] = (tmat[3][4] - (tmp * tmat[2][4]));
              tv_tmp[3] = (tv_tmp[3] - (tv_tmp[2] * tmp));
              tmp = (tmp1 * tmat[4][2]);
              tmat[4][3] = (tmat[4][3] - (tmp * tmat[2][3]));
              tmat[4][4] = (tmat[4][4] - (tmp * tmat[2][4]));
              tv_tmp[4] = (tv_tmp[4] - (tv_tmp[2] * tmp));
              tmp1 = (one / tmat[3][3]);
              tmp = (tmp1 * tmat[4][3]);
              tmat[4][4] = (tmat[4][4] - (tmp * tmat[3][4]));
              tv_tmp[4] = (tv_tmp[4] - (tv_tmp[3] * tmp));
              rsd[c1][c2][c3][4] = (tv_tmp[4] / tmat[4][4]);
              tv_tmp[3] = (tv_tmp[3] - (tmat[3][4] * rsd[c1][c2][c3][4]));
              rsd[c1][c2][c3][3] = (tv_tmp[3] / tmat[3][3]);
              tv_tmp[2] = ((tv_tmp[2] - (tmat[2][3] * rsd[c1][c2][c3][3])) - (tmat[2][4] * rsd[c1][c2][c3][4]));
              rsd[c1][c2][c3][2] = (tv_tmp[2] / tmat[2][2]);
              tv_tmp[1] = (((tv_tmp[1] - (tmat[1][2] * rsd[c1][c2][c3][2])) - (tmat[1][3] * rsd[c1][c2][c3][3])) - (tmat[1][4] * rsd[c1][c2][c3][4]));
              rsd[c1][c2][c3][1] = (tv_tmp[1] / tmat[1][1]);
              tv_tmp[0] = ((((tv_tmp[0] - (tmat[0][1] * rsd[c1][c2][c3][1])) - (tmat[0][2] * rsd[c1][c2][c3][2])) - (tmat[0][3] * rsd[c1][c2][c3][3])) - (tmat[0][4] * rsd[c1][c2][c3][4]));
              rsd[c1][c2][c3][0] = (tv_tmp[0] / tmat[0][0]);
            }
          timer_stop(7);
        }
        for (int c1 = -62; c1 < 0; c1 += 1) {
          timer_start(8);
          r43 = (4.0 / 3.0);
          c1345 = (((cc1 * cc3) * cc4) * cc5);
          c34 = (cc3 * cc4);
          for (int c2 = 1; c2 <= 62; c2 += 1)
            for (int c3 = 1; c3 <= 62; c3 += 1) {
              tmp1 = rho_i[-c1][c2][c3];
              tmp2 = (tmp1 * tmp1);
              tmp3 = (tmp1 * tmp2);
              d[c2][c3][0][0] = (one + ((dt * two) * (((tx1 * dx1) + (ty1 * dy1)) + (tz1 * dz1))));
              d[c2][c3][1][0] = zero;
              d[c2][c3][2][0] = zero;
              d[c2][c3][3][0] = zero;
              d[c2][c3][4][0] = zero;
              d[c2][c3][0][1] = (((dt * two) * ((((-tx1) * r43) - ty1) - tz1)) * ((c34 * tmp2) * u[-c1][c2][c3][1]));
              d[c2][c3][1][1] = ((one + ((((dt * two) * c34) * tmp1) * (((tx1 * r43) + ty1) + tz1))) + ((dt * two) * (((tx1 * dx2) + (ty1 * dy2)) + (tz1 * dz2))));
              d[c2][c3][2][1] = zero;
              d[c2][c3][3][1] = zero;
              d[c2][c3][4][1] = zero;
              d[c2][c3][0][2] = (((dt * two) * (((-tx1) - (ty1 * r43)) - tz1)) * ((c34 * tmp2) * u[-c1][c2][c3][2]));
              d[c2][c3][1][2] = zero;
              d[c2][c3][2][2] = ((one + ((((dt * two) * c34) * tmp1) * ((tx1 + (ty1 * r43)) + tz1))) + ((dt * two) * (((tx1 * dx3) + (ty1 * dy3)) + (tz1 * dz3))));
              d[c2][c3][3][2] = zero;
              d[c2][c3][4][2] = zero;
              d[c2][c3][0][3] = (((dt * two) * (((-tx1) - ty1) - (tz1 * r43))) * ((c34 * tmp2) * u[-c1][c2][c3][3]));
              d[c2][c3][1][3] = zero;
              d[c2][c3][2][3] = zero;
              d[c2][c3][3][3] = ((one + ((((dt * two) * c34) * tmp1) * ((tx1 + ty1) + (tz1 * r43)))) + ((dt * two) * (((tx1 * dx4) + (ty1 * dy4)) + (tz1 * dz4))));
              d[c2][c3][4][3] = zero;
              d[c2][c3][0][4] = (((-dt) * two) * ((((((((tx1 * ((r43 * c34) - c1345)) + (ty1 * (c34 - c1345))) + (tz1 * (c34 - c1345))) * (u[-c1][c2][c3][1] * u[-c1][c2][c3][1])) + ((((tx1 * (c34 - c1345)) + (ty1 * ((r43 * c34) - c1345))) + (tz1 * (c34 - c1345))) * (u[-c1][c2][c3][2] * u[-c1][c2][c3][2]))) + ((((tx1 * (c34 - c1345)) + (ty1 * (c34 - c1345))) + (tz1 * ((r43 * c34) - c1345))) * (u[-c1][c2][c3][3] * u[-c1][c2][c3][3]))) * tmp3) + (((((tx1 + ty1) + tz1) * c1345) * tmp2) * u[-c1][c2][c3][4])));
              d[c2][c3][1][4] = ((((dt * two) * (((tx1 * ((r43 * c34) - c1345)) + (ty1 * (c34 - c1345))) + (tz1 * (c34 - c1345)))) * tmp2) * u[-c1][c2][c3][1]);
              d[c2][c3][2][4] = ((((dt * two) * (((tx1 * (c34 - c1345)) + (ty1 * ((r43 * c34) - c1345))) + (tz1 * (c34 - c1345)))) * tmp2) * u[-c1][c2][c3][2]);
              d[c2][c3][3][4] = ((((dt * two) * (((tx1 * (c34 - c1345)) + (ty1 * (c34 - c1345))) + (tz1 * ((r43 * c34) - c1345)))) * tmp2) * u[-c1][c2][c3][3]);
              d[c2][c3][4][4] = ((one + ((((dt * two) * ((tx1 + ty1) + tz1)) * c1345) * tmp1)) + ((dt * two) * (((tx1 * dx5) + (ty1 * dy5)) + (tz1 * dz5))));
              tmp1 = rho_i[-c1][c2][c3 + 1];
              tmp2 = (tmp1 * tmp1);
              tmp3 = (tmp1 * tmp2);
              a[c2][c3][0][0] = (((-dt) * tx1) * dx1);
              a[c2][c3][1][0] = (dt * tx2);
              a[c2][c3][2][0] = zero;
              a[c2][c3][3][0] = zero;
              a[c2][c3][4][0] = zero;
              a[c2][c3][0][1] = (((dt * tx2) * (((-(u[-c1][c2][c3 + 1][1] * tmp1)) * (u[-c1][c2][c3 + 1][1] * tmp1)) + ((cc2 * qs[-c1][c2][c3 + 1]) * tmp1))) - ((dt * tx1) * ((((-r43) * c34) * tmp2) * u[-c1][c2][c3 + 1][1])));
              a[c2][c3][1][1] = ((((dt * tx2) * ((two - cc2) * (u[-c1][c2][c3 + 1][1] * tmp1))) - ((dt * tx1) * ((r43 * c34) * tmp1))) - ((dt * tx1) * dx2));
              a[c2][c3][2][1] = ((dt * tx2) * ((-cc2) * (u[-c1][c2][c3 + 1][2] * tmp1)));
              a[c2][c3][3][1] = ((dt * tx2) * ((-cc2) * (u[-c1][c2][c3 + 1][3] * tmp1)));
              a[c2][c3][4][1] = ((dt * tx2) * cc2);
              a[c2][c3][0][2] = (((dt * tx2) * ((-(u[-c1][c2][c3 + 1][1] * u[-c1][c2][c3 + 1][2])) * tmp2)) - ((dt * tx1) * (((-c34) * tmp2) * u[-c1][c2][c3 + 1][2])));
              a[c2][c3][1][2] = ((dt * tx2) * (u[-c1][c2][c3 + 1][2] * tmp1));
              a[c2][c3][2][2] = ((((dt * tx2) * (u[-c1][c2][c3 + 1][1] * tmp1)) - ((dt * tx1) * (c34 * tmp1))) - ((dt * tx1) * dx3));
              a[c2][c3][3][2] = zero;
              a[c2][c3][4][2] = zero;
              a[c2][c3][0][3] = (((dt * tx2) * ((-(u[-c1][c2][c3 + 1][1] * u[-c1][c2][c3 + 1][3])) * tmp2)) - ((dt * tx1) * (((-c34) * tmp2) * u[-c1][c2][c3 + 1][3])));
              a[c2][c3][1][3] = ((dt * tx2) * (u[-c1][c2][c3 + 1][3] * tmp1));
              a[c2][c3][2][3] = zero;
              a[c2][c3][3][3] = ((((dt * tx2) * (u[-c1][c2][c3 + 1][1] * tmp1)) - ((dt * tx1) * (c34 * tmp1))) - ((dt * tx1) * dx4));
              a[c2][c3][4][3] = zero;
              a[c2][c3][0][4] = (((dt * tx2) * ((((cc2 * two) * qs[-c1][c2][c3 + 1]) - (cc1 * u[-c1][c2][c3 + 1][4])) * (u[-c1][c2][c3 + 1][1] * tmp2))) - ((dt * tx1) * ((((((-((r43 * c34) - c1345)) * tmp3) * (u[-c1][c2][c3 + 1][1] * u[-c1][c2][c3 + 1][1])) - (((c34 - c1345) * tmp3) * (u[-c1][c2][c3 + 1][2] * u[-c1][c2][c3 + 1][2]))) - (((c34 - c1345) * tmp3) * (u[-c1][c2][c3 + 1][3] * u[-c1][c2][c3 + 1][3]))) - ((c1345 * tmp2) * u[-c1][c2][c3 + 1][4]))));
              a[c2][c3][1][4] = (((dt * tx2) * ((cc1 * (u[-c1][c2][c3 + 1][4] * tmp1)) - (cc2 * (((u[-c1][c2][c3 + 1][1] * u[-c1][c2][c3 + 1][1]) * tmp2) + (qs[-c1][c2][c3 + 1] * tmp1))))) - ((((dt * tx1) * ((r43 * c34) - c1345)) * tmp2) * u[-c1][c2][c3 + 1][1]));
              a[c2][c3][2][4] = (((dt * tx2) * (((-cc2) * (u[-c1][c2][c3 + 1][2] * u[-c1][c2][c3 + 1][1])) * tmp2)) - ((((dt * tx1) * (c34 - c1345)) * tmp2) * u[-c1][c2][c3 + 1][2]));
              a[c2][c3][3][4] = (((dt * tx2) * (((-cc2) * (u[-c1][c2][c3 + 1][3] * u[-c1][c2][c3 + 1][1])) * tmp2)) - ((((dt * tx1) * (c34 - c1345)) * tmp2) * u[-c1][c2][c3 + 1][3]));
              a[c2][c3][4][4] = ((((dt * tx2) * (cc1 * (u[-c1][c2][c3 + 1][1] * tmp1))) - (((dt * tx1) * c1345) * tmp1)) - ((dt * tx1) * dx5));
              tmp1 = rho_i[-c1][c2 + 1][c3];
              tmp2 = (tmp1 * tmp1);
              tmp3 = (tmp1 * tmp2);
              b[c2][c3][0][0] = (((-dt) * ty1) * dy1);
              b[c2][c3][1][0] = zero;
              b[c2][c3][2][0] = (dt * ty2);
              b[c2][c3][3][0] = zero;
              b[c2][c3][4][0] = zero;
              b[c2][c3][0][1] = (((dt * ty2) * ((-(u[-c1][c2 + 1][c3][1] * u[-c1][c2 + 1][c3][2])) * tmp2)) - ((dt * ty1) * (((-c34) * tmp2) * u[-c1][c2 + 1][c3][1])));
              b[c2][c3][1][1] = ((((dt * ty2) * (u[-c1][c2 + 1][c3][2] * tmp1)) - ((dt * ty1) * (c34 * tmp1))) - ((dt * ty1) * dy2));
              b[c2][c3][2][1] = ((dt * ty2) * (u[-c1][c2 + 1][c3][1] * tmp1));
              b[c2][c3][3][1] = zero;
              b[c2][c3][4][1] = zero;
              b[c2][c3][0][2] = (((dt * ty2) * (((-(u[-c1][c2 + 1][c3][2] * tmp1)) * (u[-c1][c2 + 1][c3][2] * tmp1)) + (cc2 * (qs[-c1][c2 + 1][c3] * tmp1)))) - ((dt * ty1) * ((((-r43) * c34) * tmp2) * u[-c1][c2 + 1][c3][2])));
              b[c2][c3][1][2] = ((dt * ty2) * ((-cc2) * (u[-c1][c2 + 1][c3][1] * tmp1)));
              b[c2][c3][2][2] = ((((dt * ty2) * ((two - cc2) * (u[-c1][c2 + 1][c3][2] * tmp1))) - ((dt * ty1) * ((r43 * c34) * tmp1))) - ((dt * ty1) * dy3));
              b[c2][c3][3][2] = ((dt * ty2) * ((-cc2) * (u[-c1][c2 + 1][c3][3] * tmp1)));
              b[c2][c3][4][2] = ((dt * ty2) * cc2);
              b[c2][c3][0][3] = (((dt * ty2) * ((-(u[-c1][c2 + 1][c3][2] * u[-c1][c2 + 1][c3][3])) * tmp2)) - ((dt * ty1) * (((-c34) * tmp2) * u[-c1][c2 + 1][c3][3])));
              b[c2][c3][1][3] = zero;
              b[c2][c3][2][3] = ((dt * ty2) * (u[-c1][c2 + 1][c3][3] * tmp1));
              b[c2][c3][3][3] = ((((dt * ty2) * (u[-c1][c2 + 1][c3][2] * tmp1)) - ((dt * ty1) * (c34 * tmp1))) - ((dt * ty1) * dy4));
              b[c2][c3][4][3] = zero;
              b[c2][c3][0][4] = (((dt * ty2) * ((((cc2 * two) * qs[-c1][c2 + 1][c3]) - (cc1 * u[-c1][c2 + 1][c3][4])) * (u[-c1][c2 + 1][c3][2] * tmp2))) - ((dt * ty1) * ((((((-(c34 - c1345)) * tmp3) * (u[-c1][c2 + 1][c3][1] * u[-c1][c2 + 1][c3][1])) - ((((r43 * c34) - c1345) * tmp3) * (u[-c1][c2 + 1][c3][2] * u[-c1][c2 + 1][c3][2]))) - (((c34 - c1345) * tmp3) * (u[-c1][c2 + 1][c3][3] * u[-c1][c2 + 1][c3][3]))) - ((c1345 * tmp2) * u[-c1][c2 + 1][c3][4]))));
              b[c2][c3][1][4] = (((dt * ty2) * (((-cc2) * (u[-c1][c2 + 1][c3][1] * u[-c1][c2 + 1][c3][2])) * tmp2)) - ((((dt * ty1) * (c34 - c1345)) * tmp2) * u[-c1][c2 + 1][c3][1]));
              b[c2][c3][2][4] = (((dt * ty2) * ((cc1 * (u[-c1][c2 + 1][c3][4] * tmp1)) - (cc2 * ((qs[-c1][c2 + 1][c3] * tmp1) + ((u[-c1][c2 + 1][c3][2] * u[-c1][c2 + 1][c3][2]) * tmp2))))) - ((((dt * ty1) * ((r43 * c34) - c1345)) * tmp2) * u[-c1][c2 + 1][c3][2]));
              b[c2][c3][3][4] = (((dt * ty2) * (((-cc2) * (u[-c1][c2 + 1][c3][2] * u[-c1][c2 + 1][c3][3])) * tmp2)) - ((((dt * ty1) * (c34 - c1345)) * tmp2) * u[-c1][c2 + 1][c3][3]));
              b[c2][c3][4][4] = ((((dt * ty2) * (cc1 * (u[-c1][c2 + 1][c3][2] * tmp1))) - (((dt * ty1) * c1345) * tmp1)) - ((dt * ty1) * dy5));
              tmp1 = rho_i[-c1 + 1][c2][c3];
              tmp2 = (tmp1 * tmp1);
              tmp3 = (tmp1 * tmp2);
              c[c2][c3][0][0] = (((-dt) * tz1) * dz1);
              c[c2][c3][1][0] = zero;
              c[c2][c3][2][0] = zero;
              c[c2][c3][3][0] = (dt * tz2);
              c[c2][c3][4][0] = zero;
              c[c2][c3][0][1] = (((dt * tz2) * ((-(u[-c1 + 1][c2][c3][1] * u[-c1 + 1][c2][c3][3])) * tmp2)) - ((dt * tz1) * (((-c34) * tmp2) * u[-c1 + 1][c2][c3][1])));
              c[c2][c3][1][1] = ((((dt * tz2) * (u[-c1 + 1][c2][c3][3] * tmp1)) - (((dt * tz1) * c34) * tmp1)) - ((dt * tz1) * dz2));
              c[c2][c3][2][1] = zero;
              c[c2][c3][3][1] = ((dt * tz2) * (u[-c1 + 1][c2][c3][1] * tmp1));
              c[c2][c3][4][1] = zero;
              c[c2][c3][0][2] = (((dt * tz2) * ((-(u[-c1 + 1][c2][c3][2] * u[-c1 + 1][c2][c3][3])) * tmp2)) - ((dt * tz1) * (((-c34) * tmp2) * u[-c1 + 1][c2][c3][2])));
              c[c2][c3][1][2] = zero;
              c[c2][c3][2][2] = ((((dt * tz2) * (u[-c1 + 1][c2][c3][3] * tmp1)) - ((dt * tz1) * (c34 * tmp1))) - ((dt * tz1) * dz3));
              c[c2][c3][3][2] = ((dt * tz2) * (u[-c1 + 1][c2][c3][2] * tmp1));
              c[c2][c3][4][2] = zero;
              c[c2][c3][0][3] = (((dt * tz2) * (((-(u[-c1 + 1][c2][c3][3] * tmp1)) * (u[-c1 + 1][c2][c3][3] * tmp1)) + (cc2 * (qs[-c1 + 1][c2][c3] * tmp1)))) - ((dt * tz1) * ((((-r43) * c34) * tmp2) * u[-c1 + 1][c2][c3][3])));
              c[c2][c3][1][3] = ((dt * tz2) * ((-cc2) * (u[-c1 + 1][c2][c3][1] * tmp1)));
              c[c2][c3][2][3] = ((dt * tz2) * ((-cc2) * (u[-c1 + 1][c2][c3][2] * tmp1)));
              c[c2][c3][3][3] = (((((dt * tz2) * (two - cc2)) * (u[-c1 + 1][c2][c3][3] * tmp1)) - ((dt * tz1) * ((r43 * c34) * tmp1))) - ((dt * tz1) * dz4));
              c[c2][c3][4][3] = ((dt * tz2) * cc2);
              c[c2][c3][0][4] = (((dt * tz2) * ((((cc2 * two) * qs[-c1 + 1][c2][c3]) - (cc1 * u[-c1 + 1][c2][c3][4])) * (u[-c1 + 1][c2][c3][3] * tmp2))) - ((dt * tz1) * ((((((-(c34 - c1345)) * tmp3) * (u[-c1 + 1][c2][c3][1] * u[-c1 + 1][c2][c3][1])) - (((c34 - c1345) * tmp3) * (u[-c1 + 1][c2][c3][2] * u[-c1 + 1][c2][c3][2]))) - ((((r43 * c34) - c1345) * tmp3) * (u[-c1 + 1][c2][c3][3] * u[-c1 + 1][c2][c3][3]))) - ((c1345 * tmp2) * u[-c1 + 1][c2][c3][4]))));
              c[c2][c3][1][4] = (((dt * tz2) * (((-cc2) * (u[-c1 + 1][c2][c3][1] * u[-c1 + 1][c2][c3][3])) * tmp2)) - ((((dt * tz1) * (c34 - c1345)) * tmp2) * u[-c1 + 1][c2][c3][1]));
              c[c2][c3][2][4] = (((dt * tz2) * (((-cc2) * (u[-c1 + 1][c2][c3][2] * u[-c1 + 1][c2][c3][3])) * tmp2)) - ((((dt * tz1) * (c34 - c1345)) * tmp2) * u[-c1 + 1][c2][c3][2]));
              c[c2][c3][3][4] = (((dt * tz2) * ((cc1 * (u[-c1 + 1][c2][c3][4] * tmp1)) - (cc2 * ((qs[-c1 + 1][c2][c3] * tmp1) + ((u[-c1 + 1][c2][c3][3] * u[-c1 + 1][c2][c3][3]) * tmp2))))) - ((((dt * tz1) * ((r43 * c34) - c1345)) * tmp2) * u[-c1 + 1][c2][c3][3]));
              c[c2][c3][4][4] = ((((dt * tz2) * (cc1 * (u[-c1 + 1][c2][c3][3] * tmp1))) - (((dt * tz1) * c1345) * tmp1)) - ((dt * tz1) * dz5));
            }
          timer_stop(8);
          timer_start(9);
          for (int c2 = -62; c2 < 0; c2 += 1)
            for (int c3 = -62; c3 < 0; c3 += 1)
              for (int c4 = 0; c4 <= 4; c4 += 1)
                tv[-c2][-c3][c4] = (omega * (((((c[-c2][-c3][0][c4] * rsd[-c1 + 1][-c2][-c3][0]) + (c[-c2][-c3][1][c4] * rsd[-c1 + 1][-c2][-c3][1])) + (c[-c2][-c3][2][c4] * rsd[-c1 + 1][-c2][-c3][2])) + (c[-c2][-c3][3][c4] * rsd[-c1 + 1][-c2][-c3][3])) + (c[-c2][-c3][4][c4] * rsd[-c1 + 1][-c2][-c3][4])));
          for (int c2 = -62; c2 < 0; c2 += 1)
            for (int c3 = -62; c3 < 0; c3 += 1) {
              for (int c4 = 0; c4 <= 4; c4 += 1)
                tv[-c2][-c3][c4] = (tv[-c2][-c3][c4] + (omega * ((((((((((b[-c2][-c3][0][c4] * rsd[-c1][-c2 + 1][-c3][0]) + (a[-c2][-c3][0][c4] * rsd[-c1][-c2][-c3 + 1][0])) + (b[-c2][-c3][1][c4] * rsd[-c1][-c2 + 1][-c3][1])) + (a[-c2][-c3][1][c4] * rsd[-c1][-c2][-c3 + 1][1])) + (b[-c2][-c3][2][c4] * rsd[-c1][-c2 + 1][-c3][2])) + (a[-c2][-c3][2][c4] * rsd[-c1][-c2][-c3 + 1][2])) + (b[-c2][-c3][3][c4] * rsd[-c1][-c2 + 1][-c3][3])) + (a[-c2][-c3][3][c4] * rsd[-c1][-c2][-c3 + 1][3])) + (b[-c2][-c3][4][c4] * rsd[-c1][-c2 + 1][-c3][4])) + (a[-c2][-c3][4][c4] * rsd[-c1][-c2][-c3 + 1][4]))));
              for (int c4 = 0; c4 <= 4; c4 += 1) {
                tmat[c4][0] = d[-c2][-c3][0][c4];
                tmat[c4][1] = d[-c2][-c3][1][c4];
                tmat[c4][2] = d[-c2][-c3][2][c4];
                tmat[c4][3] = d[-c2][-c3][3][c4];
                tmat[c4][4] = d[-c2][-c3][4][c4];
              }
              tmp1 = (one / tmat[0][0]);
              tmp = (tmp1 * tmat[1][0]);
              tmat[1][1] = (tmat[1][1] - (tmp * tmat[0][1]));
              tmat[1][2] = (tmat[1][2] - (tmp * tmat[0][2]));
              tmat[1][3] = (tmat[1][3] - (tmp * tmat[0][3]));
              tmat[1][4] = (tmat[1][4] - (tmp * tmat[0][4]));
              tv[-c2][-c3][1] = (tv[-c2][-c3][1] - (tv[-c2][-c3][0] * tmp));
              tmp = (tmp1 * tmat[2][0]);
              tmat[2][1] = (tmat[2][1] - (tmp * tmat[0][1]));
              tmat[2][2] = (tmat[2][2] - (tmp * tmat[0][2]));
              tmat[2][3] = (tmat[2][3] - (tmp * tmat[0][3]));
              tmat[2][4] = (tmat[2][4] - (tmp * tmat[0][4]));
              tv[-c2][-c3][2] = (tv[-c2][-c3][2] - (tv[-c2][-c3][0] * tmp));
              tmp = (tmp1 * tmat[3][0]);
              tmat[3][1] = (tmat[3][1] - (tmp * tmat[0][1]));
              tmat[3][2] = (tmat[3][2] - (tmp * tmat[0][2]));
              tmat[3][3] = (tmat[3][3] - (tmp * tmat[0][3]));
              tmat[3][4] = (tmat[3][4] - (tmp * tmat[0][4]));
              tv[-c2][-c3][3] = (tv[-c2][-c3][3] - (tv[-c2][-c3][0] * tmp));
              tmp = (tmp1 * tmat[4][0]);
              tmat[4][1] = (tmat[4][1] - (tmp * tmat[0][1]));
              tmat[4][2] = (tmat[4][2] - (tmp * tmat[0][2]));
              tmat[4][3] = (tmat[4][3] - (tmp * tmat[0][3]));
              tmat[4][4] = (tmat[4][4] - (tmp * tmat[0][4]));
              tv[-c2][-c3][4] = (tv[-c2][-c3][4] - (tv[-c2][-c3][0] * tmp));
              tmp1 = (one / tmat[1][1]);
              tmp = (tmp1 * tmat[2][1]);
              tmat[2][2] = (tmat[2][2] - (tmp * tmat[1][2]));
              tmat[2][3] = (tmat[2][3] - (tmp * tmat[1][3]));
              tmat[2][4] = (tmat[2][4] - (tmp * tmat[1][4]));
              tv[-c2][-c3][2] = (tv[-c2][-c3][2] - (tv[-c2][-c3][1] * tmp));
              tmp = (tmp1 * tmat[3][1]);
              tmat[3][2] = (tmat[3][2] - (tmp * tmat[1][2]));
              tmat[3][3] = (tmat[3][3] - (tmp * tmat[1][3]));
              tmat[3][4] = (tmat[3][4] - (tmp * tmat[1][4]));
              tv[-c2][-c3][3] = (tv[-c2][-c3][3] - (tv[-c2][-c3][1] * tmp));
              tmp = (tmp1 * tmat[4][1]);
              tmat[4][2] = (tmat[4][2] - (tmp * tmat[1][2]));
              tmat[4][3] = (tmat[4][3] - (tmp * tmat[1][3]));
              tmat[4][4] = (tmat[4][4] - (tmp * tmat[1][4]));
              tv[-c2][-c3][4] = (tv[-c2][-c3][4] - (tv[-c2][-c3][1] * tmp));
              tmp1 = (one / tmat[2][2]);
              tmp = (tmp1 * tmat[3][2]);
              tmat[3][3] = (tmat[3][3] - (tmp * tmat[2][3]));
              tmat[3][4] = (tmat[3][4] - (tmp * tmat[2][4]));
              tv[-c2][-c3][3] = (tv[-c2][-c3][3] - (tv[-c2][-c3][2] * tmp));
              tmp = (tmp1 * tmat[4][2]);
              tmat[4][3] = (tmat[4][3] - (tmp * tmat[2][3]));
              tmat[4][4] = (tmat[4][4] - (tmp * tmat[2][4]));
              tv[-c2][-c3][4] = (tv[-c2][-c3][4] - (tv[-c2][-c3][2] * tmp));
              tmp1 = (one / tmat[3][3]);
              tmp = (tmp1 * tmat[4][3]);
              tmat[4][4] = (tmat[4][4] - (tmp * tmat[3][4]));
              tv[-c2][-c3][4] = (tv[-c2][-c3][4] - (tv[-c2][-c3][3] * tmp));
              tv[-c2][-c3][4] = (tv[-c2][-c3][4] / tmat[4][4]);
              tv[-c2][-c3][3] = (tv[-c2][-c3][3] - (tmat[3][4] * tv[-c2][-c3][4]));
              tv[-c2][-c3][3] = (tv[-c2][-c3][3] / tmat[3][3]);
              tv[-c2][-c3][2] = ((tv[-c2][-c3][2] - (tmat[2][3] * tv[-c2][-c3][3])) - (tmat[2][4] * tv[-c2][-c3][4]));
              tv[-c2][-c3][2] = (tv[-c2][-c3][2] / tmat[2][2]);
              tv[-c2][-c3][1] = (((tv[-c2][-c3][1] - (tmat[1][2] * tv[-c2][-c3][2])) - (tmat[1][3] * tv[-c2][-c3][3])) - (tmat[1][4] * tv[-c2][-c3][4]));
              tv[-c2][-c3][1] = (tv[-c2][-c3][1] / tmat[1][1]);
              tv[-c2][-c3][0] = ((((tv[-c2][-c3][0] - (tmat[0][1] * tv[-c2][-c3][1])) - (tmat[0][2] * tv[-c2][-c3][2])) - (tmat[0][3] * tv[-c2][-c3][3])) - (tmat[0][4] * tv[-c2][-c3][4]));
              tv[-c2][-c3][0] = (tv[-c2][-c3][0] / tmat[0][0]);
              rsd[-c1][-c2][-c3][0] = (rsd[-c1][-c2][-c3][0] - tv[-c2][-c3][0]);
              rsd[-c1][-c2][-c3][1] = (rsd[-c1][-c2][-c3][1] - tv[-c2][-c3][1]);
              rsd[-c1][-c2][-c3][2] = (rsd[-c1][-c2][-c3][2] - tv[-c2][-c3][2]);
              rsd[-c1][-c2][-c3][3] = (rsd[-c1][-c2][-c3][3] - tv[-c2][-c3][3]);
              rsd[-c1][-c2][-c3][4] = (rsd[-c1][-c2][-c3][4] - tv[-c2][-c3][4]);
            }
          timer_stop(9);
        }
        timer_start(10);
        for (int c1 = 1; c1 <= 62; c1 += 1)
          for (int c2 = 1; c2 <= 62; c2 += 1)
            for (int c3 = 1; c3 <= 62; c3 += 1)
              for (int c4 = 0; c4 <= 4; c4 += 1)
                u[c1][c2][c3][c4] = (u[c1][c2][c3][c4] + (tmp_ssor * rsd[c1][c2][c3][c4]));
        timer_stop(10);
        if (((c0) % (inorm)) == 0) {
          if ((1)) {
            timer_start(11);
          }
          l2norm(64, 64, 64, (nx0), (ny0), (nz0), (1), (63), (1), (63), rsd, delunm);
          if ((1)) {
            timer_stop(11);
          }
        }
        rhs();
        if ((((c0) % (inorm)) == 0) || (c0 == itmax ? 1 : 0)) {
          if ((1)) {
            timer_start(11);
          }
          l2norm(64, 64, 64, (nx0), (ny0), (nz0), (1), (63), (1), (63), rsd, rsdnm);
          if ((1)) {
            timer_stop(11);
          }
        }
      }
      if (niter >= 2) {
        // amp_kernel
        // amp_lower
        {
          amp_lower_cc1 = (float)cc1;
          amp_lower_cc2 = (float)cc2;
          amp_lower_cc3 = (float)cc3;
          amp_lower_cc4 = (float)cc4;
          amp_lower_cc5 = (float)cc5;
          amp_lower_dt = (float)dt;
          amp_lower_dx1 = (float)dx1;
          amp_lower_dx2 = (float)dx2;
          amp_lower_dx3 = (float)dx3;
          amp_lower_dx4 = (float)dx4;
          amp_lower_dx5 = (float)dx5;
          amp_lower_dy1 = (float)dy1;
          amp_lower_dy2 = (float)dy2;
          amp_lower_dy3 = (float)dy3;
          amp_lower_dy4 = (float)dy4;
          amp_lower_dy5 = (float)dy5;
          amp_lower_dz1 = (float)dz1;
          amp_lower_dz2 = (float)dz2;
          amp_lower_dz3 = (float)dz3;
          amp_lower_dz4 = (float)dz4;
          amp_lower_dz5 = (float)dz5;
          amp_lower_omega = (float)omega;
          amp_lower_one = (float)one;
          for (int c0 = 0; c0 <= 63; c0 += 1)
            for (int c1 = ppcg_max(ppcg_max(0, c0 - 62), -c0 + 1); c1 <= ppcg_min(ppcg_min(63, c0 + 62), -c0 + 125); c1 += 1) {
              if (c0 >= 1 && c1 >= 1) {
                if (c0 <= 62 && c1 <= 62) {
                  amp_lower_qs[c0][c1][0] = (float)qs[c0][c1][0];
                  if (c0 == 1 && c1 == 1)
                    amp_lower_qs[1][1][1] = (float)qs[1][1][1];
                }
                for (int c2 = ppcg_max(1, -c0 - c1 + 4); c2 <= ppcg_min(ppcg_min(63, -c0 + 125), -c1 + 125); c2 += 1)
                  amp_lower_qs[c0][c1][c2] = (float)qs[c0][c1][c2];
              } else if (c1 == 0) {
                for (int c2 = 1; c2 <= 62; c2 += 1)
                  amp_lower_qs[c0][0][c2] = (float)qs[c0][0][c2];
              } else {
                for (int c2 = 1; c2 <= 62; c2 += 1)
                  amp_lower_qs[0][c1][c2] = (float)qs[0][c1][c2];
              }
            }
          for (int c0 = 0; c0 <= 63; c0 += 1)
            for (int c1 = ppcg_max(ppcg_max(0, c0 - 62), -c0 + 1); c1 <= ppcg_min(ppcg_min(63, c0 + 62), -c0 + 125); c1 += 1)
              for (int c2 = ppcg_max(ppcg_max(ppcg_max(ppcg_max(0, c0 - 62), -c0 + 1), c1 - 62), -c1 + 1); c2 <= ppcg_min(ppcg_min(ppcg_min(ppcg_min(63, c0 + 62), -c0 + 125), c1 + 62), -c1 + 125); c2 += 1)
                amp_lower_rho_i[c0][c1][c2] = (float)rho_i[c0][c1][c2];
          amp_lower_tmp_ssor = (float)tmp_ssor;
          amp_lower_two = (float)two;
          amp_lower_tx1 = (float)tx1;
          amp_lower_tx2 = (float)tx2;
          amp_lower_ty1 = (float)ty1;
          amp_lower_ty2 = (float)ty2;
          amp_lower_tz1 = (float)tz1;
          amp_lower_tz2 = (float)tz2;
          for (int c0 = 0; c0 <= 63; c0 += 1)
            for (int c1 = ppcg_max(ppcg_max(0, c0 - 62), -c0 + 1); c1 <= ppcg_min(ppcg_min(63, c0 + 62), -c0 + 125); c1 += 1)
              for (int c2 = ppcg_max(ppcg_max(ppcg_max(ppcg_max(0, c0 - 62), -c0 + 1), c1 - 62), -c1 + 1); c2 <= ppcg_min(ppcg_min(ppcg_min(ppcg_min(63, c0 + 62), -c0 + 125), c1 + 62), -c1 + 125); c2 += 1)
                for (int c3 = ppcg_max(ppcg_max(ppcg_max(ppcg_max(ppcg_max(ppcg_max(0, c0 - 62), -c0 + 1), c1 - 62), -c1 + 1), c2 - 62), -c2 + 1); c3 <= 4; c3 += 1)
                  amp_lower_u[c0][c1][c2][c3] = (float)u[c0][c1][c2][c3];
          amp_lower_zero = (float)zero;
          for (int c0 = niter - (niter + 8) / 10 + 1; c0 <= niter; c0 += 1) {
            timer_start(5);
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c2 = 1; c2 <= 62; c2 += 1)
                for (int c3 = 1; c3 <= 62; c3 += 1)
                  for (int c4 = 0; c4 <= 4; c4 += 1)
                    rsd[c1][c2][c3][c4] = (amp_lower_dt * rsd[c1][c2][c3][c4]);
            timer_stop(5);
            for (int c1 = 1; c1 <= 62; c1 += 1) {
              timer_start(6);
              amp_lower_r43 = (4.0 / 3.0);
              amp_lower_c1345 = (((amp_lower_cc1 * amp_lower_cc3) * amp_lower_cc4) * amp_lower_cc5);
              amp_lower_c34 = (amp_lower_cc3 * amp_lower_cc4);
              for (int c2 = 1; c2 <= 62; c2 += 1)
                for (int c3 = 1; c3 <= 62; c3 += 1) {
                  amp_lower_tmp1 = amp_lower_rho_i[c1][c2][c3];
                  amp_lower_tmp2 = (amp_lower_tmp1 * amp_lower_tmp1);
                  amp_lower_tmp3 = (amp_lower_tmp1 * amp_lower_tmp2);
                  amp_lower_d_0[c2][c3][0][0] = (amp_lower_one + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx1) + (amp_lower_ty1 * amp_lower_dy1)) + (amp_lower_tz1 * amp_lower_dz1))));
                  amp_lower_d_1[c2][c3][1][0] = amp_lower_zero;
                  amp_lower_d_2[c2][c3][2][0] = amp_lower_zero;
                  amp_lower_d_3[c2][c3][3][0] = amp_lower_zero;
                  amp_lower_d_4[c2][c3][4][0] = amp_lower_zero;
                  amp_lower_d_0[c2][c3][0][1] = ((((((-amp_lower_dt) * amp_lower_two) * (((amp_lower_tx1 * amp_lower_r43) + amp_lower_ty1) + amp_lower_tz1)) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3][1]);
                  amp_lower_d_1[c2][c3][1][1] = ((amp_lower_one + ((((amp_lower_dt * amp_lower_two) * amp_lower_c34) * amp_lower_tmp1) * (((amp_lower_tx1 * amp_lower_r43) + amp_lower_ty1) + amp_lower_tz1))) + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx2) + (amp_lower_ty1 * amp_lower_dy2)) + (amp_lower_tz1 * amp_lower_dz2))));
                  amp_lower_d_2[c2][c3][2][1] = amp_lower_zero;
                  amp_lower_d_3[c2][c3][3][1] = amp_lower_zero;
                  amp_lower_d_4[c2][c3][4][1] = amp_lower_zero;
                  amp_lower_d_0[c2][c3][0][2] = ((((((-amp_lower_dt) * amp_lower_two) * ((amp_lower_tx1 + (amp_lower_ty1 * amp_lower_r43)) + amp_lower_tz1)) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3][2]);
                  amp_lower_d_1[c2][c3][1][2] = amp_lower_zero;
                  amp_lower_d_2[c2][c3][2][2] = ((amp_lower_one + ((((amp_lower_dt * amp_lower_two) * amp_lower_c34) * amp_lower_tmp1) * ((amp_lower_tx1 + (amp_lower_ty1 * amp_lower_r43)) + amp_lower_tz1))) + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx3) + (amp_lower_ty1 * amp_lower_dy3)) + (amp_lower_tz1 * amp_lower_dz3))));
                  amp_lower_d_3[c2][c3][3][2] = amp_lower_zero;
                  amp_lower_d_4[c2][c3][4][2] = amp_lower_zero;
                  amp_lower_d_0[c2][c3][0][3] = ((((((-amp_lower_dt) * amp_lower_two) * ((amp_lower_tx1 + amp_lower_ty1) + (amp_lower_tz1 * amp_lower_r43))) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3][3]);
                  amp_lower_d_1[c2][c3][1][3] = amp_lower_zero;
                  amp_lower_d_2[c2][c3][2][3] = amp_lower_zero;
                  amp_lower_d_3[c2][c3][3][3] = ((amp_lower_one + ((((amp_lower_dt * amp_lower_two) * amp_lower_c34) * amp_lower_tmp1) * ((amp_lower_tx1 + amp_lower_ty1) + (amp_lower_tz1 * amp_lower_r43)))) + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx4) + (amp_lower_ty1 * amp_lower_dy4)) + (amp_lower_tz1 * amp_lower_dz4))));
                  amp_lower_d_4[c2][c3][4][3] = amp_lower_zero;
                  amp_lower_d_0[c2][c3][0][4] = (((-amp_lower_dt) * amp_lower_two) * ((((((((amp_lower_tx1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) + (amp_lower_ty1 * (amp_lower_c34 - amp_lower_c1345))) + (amp_lower_tz1 * (amp_lower_c34 - amp_lower_c1345))) * (amp_lower_u[c1][c2][c3][1] * amp_lower_u[c1][c2][c3][1])) + ((((amp_lower_tx1 * (amp_lower_c34 - amp_lower_c1345)) + (amp_lower_ty1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345))) + (amp_lower_tz1 * (amp_lower_c34 - amp_lower_c1345))) * (amp_lower_u[c1][c2][c3][2] * amp_lower_u[c1][c2][c3][2]))) + ((((amp_lower_tx1 * (amp_lower_c34 - amp_lower_c1345)) + (amp_lower_ty1 * (amp_lower_c34 - amp_lower_c1345))) + (amp_lower_tz1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345))) * (amp_lower_u[c1][c2][c3][3] * amp_lower_u[c1][c2][c3][3]))) * amp_lower_tmp3) + (((((amp_lower_tx1 + amp_lower_ty1) + amp_lower_tz1) * amp_lower_c1345) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3][4])));
                  amp_lower_d_1[c2][c3][1][4] = ((((amp_lower_dt * amp_lower_two) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3][1]) * (((amp_lower_tx1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) + (amp_lower_ty1 * (amp_lower_c34 - amp_lower_c1345))) + (amp_lower_tz1 * (amp_lower_c34 - amp_lower_c1345))));
                  amp_lower_d_2[c2][c3][2][4] = ((((amp_lower_dt * amp_lower_two) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3][2]) * (((amp_lower_tx1 * (amp_lower_c34 - amp_lower_c1345)) + (amp_lower_ty1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345))) + (amp_lower_tz1 * (amp_lower_c34 - amp_lower_c1345))));
                  amp_lower_d_3[c2][c3][3][4] = ((((amp_lower_dt * amp_lower_two) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3][3]) * (((amp_lower_tx1 * (amp_lower_c34 - amp_lower_c1345)) + (amp_lower_ty1 * (amp_lower_c34 - amp_lower_c1345))) + (amp_lower_tz1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345))));
                  amp_lower_d_4[c2][c3][4][4] = ((amp_lower_one + ((((amp_lower_dt * amp_lower_two) * ((amp_lower_tx1 + amp_lower_ty1) + amp_lower_tz1)) * amp_lower_c1345) * amp_lower_tmp1)) + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx5) + (amp_lower_ty1 * amp_lower_dy5)) + (amp_lower_tz1 * amp_lower_dz5))));
                  amp_lower_tmp1 = amp_lower_rho_i[c1 - 1][c2][c3];
                  amp_lower_tmp2 = (amp_lower_tmp1 * amp_lower_tmp1);
                  amp_lower_tmp3 = (amp_lower_tmp1 * amp_lower_tmp2);
                  amp_lower_a_0[c2][c3][0][0] = (((-amp_lower_dt) * amp_lower_tz1) * amp_lower_dz1);
                  amp_lower_a_1[c2][c3][1][0] = amp_lower_zero;
                  amp_lower_a_2[c2][c3][2][0] = amp_lower_zero;
                  amp_lower_a_3[c2][c3][3][0] = ((-amp_lower_dt) * amp_lower_tz2);
                  amp_lower_a_4[c2][c3][4][0] = amp_lower_zero;
                  amp_lower_a_0[c2][c3][0][1] = ((((-amp_lower_dt) * amp_lower_tz2) * ((-(amp_lower_u[c1 - 1][c2][c3][1] * amp_lower_u[c1 - 1][c2][c3][3])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tz1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1 - 1][c2][c3][1])));
                  amp_lower_a_1[c2][c3][1][1] = (((((-amp_lower_dt) * amp_lower_tz2) * (amp_lower_u[c1 - 1][c2][c3][3] * amp_lower_tmp1)) - (((amp_lower_dt * amp_lower_tz1) * amp_lower_c34) * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tz1) * amp_lower_dz2));
                  amp_lower_a_2[c2][c3][2][1] = amp_lower_zero;
                  amp_lower_a_3[c2][c3][3][1] = (((-amp_lower_dt) * amp_lower_tz2) * (amp_lower_u[c1 - 1][c2][c3][1] * amp_lower_tmp1));
                  amp_lower_a_4[c2][c3][4][1] = amp_lower_zero;
                  amp_lower_a_0[c2][c3][0][2] = ((((-amp_lower_dt) * amp_lower_tz2) * ((-(amp_lower_u[c1 - 1][c2][c3][2] * amp_lower_u[c1 - 1][c2][c3][3])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tz1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1 - 1][c2][c3][2])));
                  amp_lower_a_1[c2][c3][1][2] = amp_lower_zero;
                  amp_lower_a_2[c2][c3][2][2] = (((((-amp_lower_dt) * amp_lower_tz2) * (amp_lower_u[c1 - 1][c2][c3][3] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tz1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tz1) * amp_lower_dz3));
                  amp_lower_a_3[c2][c3][3][2] = (((-amp_lower_dt) * amp_lower_tz2) * (amp_lower_u[c1 - 1][c2][c3][2] * amp_lower_tmp1));
                  amp_lower_a_4[c2][c3][4][2] = amp_lower_zero;
                  amp_lower_a_0[c2][c3][0][3] = ((((-amp_lower_dt) * amp_lower_tz2) * (((-(amp_lower_u[c1 - 1][c2][c3][3] * amp_lower_tmp1)) * (amp_lower_u[c1 - 1][c2][c3][3] * amp_lower_tmp1)) + ((amp_lower_cc2 * amp_lower_qs[c1 - 1][c2][c3]) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tz1) * ((((-amp_lower_r43) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1 - 1][c2][c3][3])));
                  amp_lower_a_1[c2][c3][1][3] = (((-amp_lower_dt) * amp_lower_tz2) * ((-amp_lower_cc2) * (amp_lower_u[c1 - 1][c2][c3][1] * amp_lower_tmp1)));
                  amp_lower_a_2[c2][c3][2][3] = (((-amp_lower_dt) * amp_lower_tz2) * ((-amp_lower_cc2) * (amp_lower_u[c1 - 1][c2][c3][2] * amp_lower_tmp1)));
                  amp_lower_a_3[c2][c3][3][3] = ((((((-amp_lower_dt) * amp_lower_tz2) * (amp_lower_two - amp_lower_cc2)) * (amp_lower_u[c1 - 1][c2][c3][3] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tz1) * ((amp_lower_r43 * amp_lower_c34) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tz1) * amp_lower_dz4));
                  amp_lower_a_4[c2][c3][4][3] = (((-amp_lower_dt) * amp_lower_tz2) * amp_lower_cc2);
                  amp_lower_a_0[c2][c3][0][4] = ((((-amp_lower_dt) * amp_lower_tz2) * (((((amp_lower_cc2 * amp_lower_two) * amp_lower_qs[c1 - 1][c2][c3]) - (amp_lower_cc1 * amp_lower_u[c1 - 1][c2][c3][4])) * amp_lower_u[c1 - 1][c2][c3][3]) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tz1) * ((((((-(amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp3) * (amp_lower_u[c1 - 1][c2][c3][1] * amp_lower_u[c1 - 1][c2][c3][1])) - (((amp_lower_c34 - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[c1 - 1][c2][c3][2] * amp_lower_u[c1 - 1][c2][c3][2]))) - ((((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[c1 - 1][c2][c3][3] * amp_lower_u[c1 - 1][c2][c3][3]))) - ((amp_lower_c1345 * amp_lower_tmp2) * amp_lower_u[c1 - 1][c2][c3][4]))));
                  amp_lower_a_1[c2][c3][1][4] = ((((-amp_lower_dt) * amp_lower_tz2) * (((-amp_lower_cc2) * (amp_lower_u[c1 - 1][c2][c3][1] * amp_lower_u[c1 - 1][c2][c3][3])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_tz1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1 - 1][c2][c3][1]));
                  amp_lower_a_2[c2][c3][2][4] = ((((-amp_lower_dt) * amp_lower_tz2) * (((-amp_lower_cc2) * (amp_lower_u[c1 - 1][c2][c3][2] * amp_lower_u[c1 - 1][c2][c3][3])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_tz1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1 - 1][c2][c3][2]));
                  amp_lower_a_3[c2][c3][3][4] = ((((-amp_lower_dt) * amp_lower_tz2) * ((amp_lower_cc1 * (amp_lower_u[c1 - 1][c2][c3][4] * amp_lower_tmp1)) - (amp_lower_cc2 * ((amp_lower_qs[c1 - 1][c2][c3] * amp_lower_tmp1) + ((amp_lower_u[c1 - 1][c2][c3][3] * amp_lower_u[c1 - 1][c2][c3][3]) * amp_lower_tmp2))))) - ((((amp_lower_dt * amp_lower_tz1) * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1 - 1][c2][c3][3]));
                  amp_lower_a_4[c2][c3][4][4] = (((((-amp_lower_dt) * amp_lower_tz2) * (amp_lower_cc1 * (amp_lower_u[c1 - 1][c2][c3][3] * amp_lower_tmp1))) - (((amp_lower_dt * amp_lower_tz1) * amp_lower_c1345) * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tz1) * amp_lower_dz5));
                  amp_lower_tmp1 = amp_lower_rho_i[c1][c2 - 1][c3];
                  amp_lower_tmp2 = (amp_lower_tmp1 * amp_lower_tmp1);
                  amp_lower_tmp3 = (amp_lower_tmp1 * amp_lower_tmp2);
                  amp_lower_b_0[c2][c3][0][0] = (((-amp_lower_dt) * amp_lower_ty1) * amp_lower_dy1);
                  amp_lower_b_1[c2][c3][1][0] = amp_lower_zero;
                  amp_lower_b_2[c2][c3][2][0] = ((-amp_lower_dt) * amp_lower_ty2);
                  amp_lower_b_3[c2][c3][3][0] = amp_lower_zero;
                  amp_lower_b_4[c2][c3][4][0] = amp_lower_zero;
                  amp_lower_b_0[c2][c3][0][1] = ((((-amp_lower_dt) * amp_lower_ty2) * ((-(amp_lower_u[c1][c2 - 1][c3][1] * amp_lower_u[c1][c2 - 1][c3][2])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_ty1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2 - 1][c3][1])));
                  amp_lower_b_1[c2][c3][1][1] = (((((-amp_lower_dt) * amp_lower_ty2) * (amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_ty1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_ty1) * amp_lower_dy2));
                  amp_lower_b_2[c2][c3][2][1] = (((-amp_lower_dt) * amp_lower_ty2) * (amp_lower_u[c1][c2 - 1][c3][1] * amp_lower_tmp1));
                  amp_lower_b_3[c2][c3][3][1] = amp_lower_zero;
                  amp_lower_b_4[c2][c3][4][1] = amp_lower_zero;
                  amp_lower_b_0[c2][c3][0][2] = ((((-amp_lower_dt) * amp_lower_ty2) * (((-(amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_tmp1)) * (amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_tmp1)) + (amp_lower_cc2 * (amp_lower_qs[c1][c2 - 1][c3] * amp_lower_tmp1)))) - ((amp_lower_dt * amp_lower_ty1) * ((((-amp_lower_r43) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2 - 1][c3][2])));
                  amp_lower_b_1[c2][c3][1][2] = (((-amp_lower_dt) * amp_lower_ty2) * ((-amp_lower_cc2) * (amp_lower_u[c1][c2 - 1][c3][1] * amp_lower_tmp1)));
                  amp_lower_b_2[c2][c3][2][2] = (((((-amp_lower_dt) * amp_lower_ty2) * ((amp_lower_two - amp_lower_cc2) * (amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_ty1) * ((amp_lower_r43 * amp_lower_c34) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_ty1) * amp_lower_dy3));
                  amp_lower_b_3[c2][c3][3][2] = (((-amp_lower_dt) * amp_lower_ty2) * ((-amp_lower_cc2) * (amp_lower_u[c1][c2 - 1][c3][3] * amp_lower_tmp1)));
                  amp_lower_b_4[c2][c3][4][2] = (((-amp_lower_dt) * amp_lower_ty2) * amp_lower_cc2);
                  amp_lower_b_0[c2][c3][0][3] = ((((-amp_lower_dt) * amp_lower_ty2) * ((-(amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_u[c1][c2 - 1][c3][3])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_ty1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2 - 1][c3][3])));
                  amp_lower_b_1[c2][c3][1][3] = amp_lower_zero;
                  amp_lower_b_2[c2][c3][2][3] = (((-amp_lower_dt) * amp_lower_ty2) * (amp_lower_u[c1][c2 - 1][c3][3] * amp_lower_tmp1));
                  amp_lower_b_3[c2][c3][3][3] = (((((-amp_lower_dt) * amp_lower_ty2) * (amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_ty1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_ty1) * amp_lower_dy4));
                  amp_lower_b_4[c2][c3][4][3] = amp_lower_zero;
                  amp_lower_b_0[c2][c3][0][4] = ((((-amp_lower_dt) * amp_lower_ty2) * ((((amp_lower_cc2 * amp_lower_two) * amp_lower_qs[c1][c2 - 1][c3]) - (amp_lower_cc1 * amp_lower_u[c1][c2 - 1][c3][4])) * (amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_tmp2))) - ((amp_lower_dt * amp_lower_ty1) * ((((((-(amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp3) * (amp_lower_u[c1][c2 - 1][c3][1] * amp_lower_u[c1][c2 - 1][c3][1])) - ((((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_u[c1][c2 - 1][c3][2]))) - (((amp_lower_c34 - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[c1][c2 - 1][c3][3] * amp_lower_u[c1][c2 - 1][c3][3]))) - ((amp_lower_c1345 * amp_lower_tmp2) * amp_lower_u[c1][c2 - 1][c3][4]))));
                  amp_lower_b_1[c2][c3][1][4] = ((((-amp_lower_dt) * amp_lower_ty2) * (((-amp_lower_cc2) * (amp_lower_u[c1][c2 - 1][c3][1] * amp_lower_u[c1][c2 - 1][c3][2])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_ty1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1][c2 - 1][c3][1]));
                  amp_lower_b_2[c2][c3][2][4] = ((((-amp_lower_dt) * amp_lower_ty2) * ((amp_lower_cc1 * (amp_lower_u[c1][c2 - 1][c3][4] * amp_lower_tmp1)) - (amp_lower_cc2 * ((amp_lower_qs[c1][c2 - 1][c3] * amp_lower_tmp1) + ((amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_u[c1][c2 - 1][c3][2]) * amp_lower_tmp2))))) - ((((amp_lower_dt * amp_lower_ty1) * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1][c2 - 1][c3][2]));
                  amp_lower_b_3[c2][c3][3][4] = ((((-amp_lower_dt) * amp_lower_ty2) * (((-amp_lower_cc2) * (amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_u[c1][c2 - 1][c3][3])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_ty1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1][c2 - 1][c3][3]));
                  amp_lower_b_4[c2][c3][4][4] = (((((-amp_lower_dt) * amp_lower_ty2) * (amp_lower_cc1 * (amp_lower_u[c1][c2 - 1][c3][2] * amp_lower_tmp1))) - (((amp_lower_dt * amp_lower_ty1) * amp_lower_c1345) * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_ty1) * amp_lower_dy5));
                  amp_lower_tmp1 = amp_lower_rho_i[c1][c2][c3 - 1];
                  amp_lower_tmp2 = (amp_lower_tmp1 * amp_lower_tmp1);
                  amp_lower_tmp3 = (amp_lower_tmp1 * amp_lower_tmp2);
                  amp_lower_c_0[c2][c3][0][0] = (((-amp_lower_dt) * amp_lower_tx1) * amp_lower_dx1);
                  amp_lower_c_1[c2][c3][1][0] = ((-amp_lower_dt) * amp_lower_tx2);
                  amp_lower_c_2[c2][c3][2][0] = amp_lower_zero;
                  amp_lower_c_3[c2][c3][3][0] = amp_lower_zero;
                  amp_lower_c_4[c2][c3][4][0] = amp_lower_zero;
                  amp_lower_c_0[c2][c3][0][1] = ((((-amp_lower_dt) * amp_lower_tx2) * (((-(amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_tmp1)) * (amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_tmp1)) + ((amp_lower_cc2 * amp_lower_qs[c1][c2][c3 - 1]) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * ((((-amp_lower_r43) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3 - 1][1])));
                  amp_lower_c_1[c2][c3][1][1] = (((((-amp_lower_dt) * amp_lower_tx2) * ((amp_lower_two - amp_lower_cc2) * (amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * ((amp_lower_r43 * amp_lower_c34) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * amp_lower_dx2));
                  amp_lower_c_2[c2][c3][2][1] = (((-amp_lower_dt) * amp_lower_tx2) * ((-amp_lower_cc2) * (amp_lower_u[c1][c2][c3 - 1][2] * amp_lower_tmp1)));
                  amp_lower_c_3[c2][c3][3][1] = (((-amp_lower_dt) * amp_lower_tx2) * ((-amp_lower_cc2) * (amp_lower_u[c1][c2][c3 - 1][3] * amp_lower_tmp1)));
                  amp_lower_c_4[c2][c3][4][1] = (((-amp_lower_dt) * amp_lower_tx2) * amp_lower_cc2);
                  amp_lower_c_0[c2][c3][0][2] = ((((-amp_lower_dt) * amp_lower_tx2) * ((-(amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_u[c1][c2][c3 - 1][2])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tx1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3 - 1][2])));
                  amp_lower_c_1[c2][c3][1][2] = (((-amp_lower_dt) * amp_lower_tx2) * (amp_lower_u[c1][c2][c3 - 1][2] * amp_lower_tmp1));
                  amp_lower_c_2[c2][c3][2][2] = (((((-amp_lower_dt) * amp_lower_tx2) * (amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tx1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * amp_lower_dx3));
                  amp_lower_c_3[c2][c3][3][2] = amp_lower_zero;
                  amp_lower_c_4[c2][c3][4][2] = amp_lower_zero;
                  amp_lower_c_0[c2][c3][0][3] = ((((-amp_lower_dt) * amp_lower_tx2) * ((-(amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_u[c1][c2][c3 - 1][3])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tx1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3 - 1][3])));
                  amp_lower_c_1[c2][c3][1][3] = (((-amp_lower_dt) * amp_lower_tx2) * (amp_lower_u[c1][c2][c3 - 1][3] * amp_lower_tmp1));
                  amp_lower_c_2[c2][c3][2][3] = amp_lower_zero;
                  amp_lower_c_3[c2][c3][3][3] = (((((-amp_lower_dt) * amp_lower_tx2) * (amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tx1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * amp_lower_dx4));
                  amp_lower_c_4[c2][c3][4][3] = amp_lower_zero;
                  amp_lower_c_0[c2][c3][0][4] = ((((-amp_lower_dt) * amp_lower_tx2) * (((((amp_lower_cc2 * amp_lower_two) * amp_lower_qs[c1][c2][c3 - 1]) - (amp_lower_cc1 * amp_lower_u[c1][c2][c3 - 1][4])) * amp_lower_u[c1][c2][c3 - 1][1]) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tx1) * ((((((-((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) * amp_lower_tmp3) * (amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_u[c1][c2][c3 - 1][1])) - (((amp_lower_c34 - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[c1][c2][c3 - 1][2] * amp_lower_u[c1][c2][c3 - 1][2]))) - (((amp_lower_c34 - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[c1][c2][c3 - 1][3] * amp_lower_u[c1][c2][c3 - 1][3]))) - ((amp_lower_c1345 * amp_lower_tmp2) * amp_lower_u[c1][c2][c3 - 1][4]))));
                  amp_lower_c_1[c2][c3][1][4] = ((((-amp_lower_dt) * amp_lower_tx2) * ((amp_lower_cc1 * (amp_lower_u[c1][c2][c3 - 1][4] * amp_lower_tmp1)) - (amp_lower_cc2 * (((amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_u[c1][c2][c3 - 1][1]) * amp_lower_tmp2) + (amp_lower_qs[c1][c2][c3 - 1] * amp_lower_tmp1))))) - ((((amp_lower_dt * amp_lower_tx1) * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3 - 1][1]));
                  amp_lower_c_2[c2][c3][2][4] = ((((-amp_lower_dt) * amp_lower_tx2) * (((-amp_lower_cc2) * (amp_lower_u[c1][c2][c3 - 1][2] * amp_lower_u[c1][c2][c3 - 1][1])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_tx1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3 - 1][2]));
                  amp_lower_c_3[c2][c3][3][4] = ((((-amp_lower_dt) * amp_lower_tx2) * (((-amp_lower_cc2) * (amp_lower_u[c1][c2][c3 - 1][3] * amp_lower_u[c1][c2][c3 - 1][1])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_tx1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[c1][c2][c3 - 1][3]));
                  amp_lower_c_4[c2][c3][4][4] = (((((-amp_lower_dt) * amp_lower_tx2) * (amp_lower_cc1 * (amp_lower_u[c1][c2][c3 - 1][1] * amp_lower_tmp1))) - (((amp_lower_dt * amp_lower_tx1) * amp_lower_c1345) * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tx1) * amp_lower_dx5));
                }
              timer_stop(6);
              timer_start(7);
              for (int c2 = 1; c2 <= 62; c2 += 1)
                for (int c3 = 1; c3 <= 62; c3 += 1)
                  for (int c4 = 0; c4 <= 4; c4 += 1)
                    rsd[c1][c2][c3][c4] = (rsd[c1][c2][c3][c4] - (amp_lower_omega * (((((amp_lower_a_0[c2][c3][0][c4] * rsd[c1 - 1][c2][c3][0]) + (amp_lower_a_1[c2][c3][1][c4] * rsd[c1 - 1][c2][c3][1])) + (amp_lower_a_2[c2][c3][2][c4] * rsd[c1 - 1][c2][c3][2])) + (amp_lower_a_3[c2][c3][3][c4] * rsd[c1 - 1][c2][c3][3])) + (amp_lower_a_4[c2][c3][4][c4] * rsd[c1 - 1][c2][c3][4]))));
              for (int c2 = 1; c2 <= 62; c2 += 1)
                for (int c3 = 1; c3 <= 62; c3 += 1) {
                  for (int c4 = 0; c4 <= 4; c4 += 1)
                    amp_lower_tv_tmp[c4] = (rsd[c1][c2][c3][c4] - (amp_lower_omega * ((((((((((amp_lower_b_0[c2][c3][0][c4] * rsd[c1][c2 - 1][c3][0]) + (amp_lower_c_0[c2][c3][0][c4] * rsd[c1][c2][c3 - 1][0])) + (amp_lower_b_1[c2][c3][1][c4] * rsd[c1][c2 - 1][c3][1])) + (amp_lower_c_1[c2][c3][1][c4] * rsd[c1][c2][c3 - 1][1])) + (amp_lower_b_2[c2][c3][2][c4] * rsd[c1][c2 - 1][c3][2])) + (amp_lower_c_2[c2][c3][2][c4] * rsd[c1][c2][c3 - 1][2])) + (amp_lower_b_3[c2][c3][3][c4] * rsd[c1][c2 - 1][c3][3])) + (amp_lower_c_3[c2][c3][3][c4] * rsd[c1][c2][c3 - 1][3])) + (amp_lower_b_4[c2][c3][4][c4] * rsd[c1][c2 - 1][c3][4])) + (amp_lower_c_4[c2][c3][4][c4] * rsd[c1][c2][c3 - 1][4]))));
                  for (int c4 = 0; c4 <= 4; c4 += 1) {
                    amp_lower_tmat_0[c4][0] = amp_lower_d_0[c2][c3][0][c4];
                    amp_lower_tmat_1[c4][1] = amp_lower_d_1[c2][c3][1][c4];
                    amp_lower_tmat_2[c4][2] = amp_lower_d_2[c2][c3][2][c4];
                    amp_lower_tmat_3[c4][3] = amp_lower_d_3[c2][c3][3][c4];
                    amp_lower_tmat_4[c4][4] = amp_lower_d_4[c2][c3][4][c4];
                  }
                  amp_lower_tmp1 = (amp_lower_one / amp_lower_tmat_0[0][0]);
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_0[1][0]);
                  amp_lower_tmat_1[1][1] = (amp_lower_tmat_1[1][1] - (amp_lower_tmp * amp_lower_tmat_1[0][1]));
                  amp_lower_tmat_2[1][2] = (amp_lower_tmat_2[1][2] - (amp_lower_tmp * amp_lower_tmat_2[0][2]));
                  amp_lower_tmat_3[1][3] = (amp_lower_tmat_3[1][3] - (amp_lower_tmp * amp_lower_tmat_3[0][3]));
                  amp_lower_tmat_4[1][4] = (amp_lower_tmat_4[1][4] - (amp_lower_tmp * amp_lower_tmat_4[0][4]));
                  amp_lower_tv_tmp[1] = (amp_lower_tv_tmp[1] - (amp_lower_tv_tmp[0] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_0[2][0]);
                  amp_lower_tmat_1[2][1] = (amp_lower_tmat_1[2][1] - (amp_lower_tmp * amp_lower_tmat_1[0][1]));
                  amp_lower_tmat_2[2][2] = (amp_lower_tmat_2[2][2] - (amp_lower_tmp * amp_lower_tmat_2[0][2]));
                  amp_lower_tmat_3[2][3] = (amp_lower_tmat_3[2][3] - (amp_lower_tmp * amp_lower_tmat_3[0][3]));
                  amp_lower_tmat_4[2][4] = (amp_lower_tmat_4[2][4] - (amp_lower_tmp * amp_lower_tmat_4[0][4]));
                  amp_lower_tv_tmp[2] = (amp_lower_tv_tmp[2] - (amp_lower_tv_tmp[0] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_0[3][0]);
                  amp_lower_tmat_1[3][1] = (amp_lower_tmat_1[3][1] - (amp_lower_tmp * amp_lower_tmat_1[0][1]));
                  amp_lower_tmat_2[3][2] = (amp_lower_tmat_2[3][2] - (amp_lower_tmp * amp_lower_tmat_2[0][2]));
                  amp_lower_tmat_3[3][3] = (amp_lower_tmat_3[3][3] - (amp_lower_tmp * amp_lower_tmat_3[0][3]));
                  amp_lower_tmat_4[3][4] = (amp_lower_tmat_4[3][4] - (amp_lower_tmp * amp_lower_tmat_4[0][4]));
                  amp_lower_tv_tmp[3] = (amp_lower_tv_tmp[3] - (amp_lower_tv_tmp[0] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_0[4][0]);
                  amp_lower_tmat_1[4][1] = (amp_lower_tmat_1[4][1] - (amp_lower_tmp * amp_lower_tmat_1[0][1]));
                  amp_lower_tmat_2[4][2] = (amp_lower_tmat_2[4][2] - (amp_lower_tmp * amp_lower_tmat_2[0][2]));
                  amp_lower_tmat_3[4][3] = (amp_lower_tmat_3[4][3] - (amp_lower_tmp * amp_lower_tmat_3[0][3]));
                  amp_lower_tmat_4[4][4] = (amp_lower_tmat_4[4][4] - (amp_lower_tmp * amp_lower_tmat_4[0][4]));
                  amp_lower_tv_tmp[4] = (amp_lower_tv_tmp[4] - (amp_lower_tv_tmp[0] * amp_lower_tmp));
                  amp_lower_tmp1 = (amp_lower_one / amp_lower_tmat_1[1][1]);
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_1[2][1]);
                  amp_lower_tmat_2[2][2] = (amp_lower_tmat_2[2][2] - (amp_lower_tmp * amp_lower_tmat_2[1][2]));
                  amp_lower_tmat_3[2][3] = (amp_lower_tmat_3[2][3] - (amp_lower_tmp * amp_lower_tmat_3[1][3]));
                  amp_lower_tmat_4[2][4] = (amp_lower_tmat_4[2][4] - (amp_lower_tmp * amp_lower_tmat_4[1][4]));
                  amp_lower_tv_tmp[2] = (amp_lower_tv_tmp[2] - (amp_lower_tv_tmp[1] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_1[3][1]);
                  amp_lower_tmat_2[3][2] = (amp_lower_tmat_2[3][2] - (amp_lower_tmp * amp_lower_tmat_2[1][2]));
                  amp_lower_tmat_3[3][3] = (amp_lower_tmat_3[3][3] - (amp_lower_tmp * amp_lower_tmat_3[1][3]));
                  amp_lower_tmat_4[3][4] = (amp_lower_tmat_4[3][4] - (amp_lower_tmp * amp_lower_tmat_4[1][4]));
                  amp_lower_tv_tmp[3] = (amp_lower_tv_tmp[3] - (amp_lower_tv_tmp[1] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_1[4][1]);
                  amp_lower_tmat_2[4][2] = (amp_lower_tmat_2[4][2] - (amp_lower_tmp * amp_lower_tmat_2[1][2]));
                  amp_lower_tmat_3[4][3] = (amp_lower_tmat_3[4][3] - (amp_lower_tmp * amp_lower_tmat_3[1][3]));
                  amp_lower_tmat_4[4][4] = (amp_lower_tmat_4[4][4] - (amp_lower_tmp * amp_lower_tmat_4[1][4]));
                  amp_lower_tv_tmp[4] = (amp_lower_tv_tmp[4] - (amp_lower_tv_tmp[1] * amp_lower_tmp));
                  amp_lower_tmp1 = (amp_lower_one / amp_lower_tmat_2[2][2]);
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_2[3][2]);
                  amp_lower_tmat_3[3][3] = (amp_lower_tmat_3[3][3] - (amp_lower_tmp * amp_lower_tmat_3[2][3]));
                  amp_lower_tmat_4[3][4] = (amp_lower_tmat_4[3][4] - (amp_lower_tmp * amp_lower_tmat_4[2][4]));
                  amp_lower_tv_tmp[3] = (amp_lower_tv_tmp[3] - (amp_lower_tv_tmp[2] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_2[4][2]);
                  amp_lower_tmat_3[4][3] = (amp_lower_tmat_3[4][3] - (amp_lower_tmp * amp_lower_tmat_3[2][3]));
                  amp_lower_tmat_4[4][4] = (amp_lower_tmat_4[4][4] - (amp_lower_tmp * amp_lower_tmat_4[2][4]));
                  amp_lower_tv_tmp[4] = (amp_lower_tv_tmp[4] - (amp_lower_tv_tmp[2] * amp_lower_tmp));
                  amp_lower_tmp1 = (amp_lower_one / amp_lower_tmat_3[3][3]);
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_3[4][3]);
                  amp_lower_tmat_4[4][4] = (amp_lower_tmat_4[4][4] - (amp_lower_tmp * amp_lower_tmat_4[3][4]));
                  amp_lower_tv_tmp[4] = (amp_lower_tv_tmp[4] - (amp_lower_tv_tmp[3] * amp_lower_tmp));
                  rsd[c1][c2][c3][4] = (amp_lower_tv_tmp[4] / amp_lower_tmat_4[4][4]);
                  amp_lower_tv_tmp[3] = (amp_lower_tv_tmp[3] - (amp_lower_tmat_4[3][4] * rsd[c1][c2][c3][4]));
                  rsd[c1][c2][c3][3] = (amp_lower_tv_tmp[3] / amp_lower_tmat_3[3][3]);
                  amp_lower_tv_tmp[2] = ((amp_lower_tv_tmp[2] - (amp_lower_tmat_3[2][3] * rsd[c1][c2][c3][3])) - (amp_lower_tmat_4[2][4] * rsd[c1][c2][c3][4]));
                  rsd[c1][c2][c3][2] = (amp_lower_tv_tmp[2] / amp_lower_tmat_2[2][2]);
                  amp_lower_tv_tmp[1] = (((amp_lower_tv_tmp[1] - (amp_lower_tmat_2[1][2] * rsd[c1][c2][c3][2])) - (amp_lower_tmat_3[1][3] * rsd[c1][c2][c3][3])) - (amp_lower_tmat_4[1][4] * rsd[c1][c2][c3][4]));
                  rsd[c1][c2][c3][1] = (amp_lower_tv_tmp[1] / amp_lower_tmat_1[1][1]);
                  amp_lower_tv_tmp[0] = ((((amp_lower_tv_tmp[0] - (amp_lower_tmat_1[0][1] * rsd[c1][c2][c3][1])) - (amp_lower_tmat_2[0][2] * rsd[c1][c2][c3][2])) - (amp_lower_tmat_3[0][3] * rsd[c1][c2][c3][3])) - (amp_lower_tmat_4[0][4] * rsd[c1][c2][c3][4]));
                  rsd[c1][c2][c3][0] = (amp_lower_tv_tmp[0] / amp_lower_tmat_0[0][0]);
                }
              timer_stop(7);
            }
            for (int c1 = -62; c1 < 0; c1 += 1) {
              timer_start(8);
              amp_lower_r43 = (4.0 / 3.0);
              amp_lower_c1345 = (((amp_lower_cc1 * amp_lower_cc3) * amp_lower_cc4) * amp_lower_cc5);
              amp_lower_c34 = (amp_lower_cc3 * amp_lower_cc4);
              for (int c2 = 1; c2 <= 62; c2 += 1)
                for (int c3 = 1; c3 <= 62; c3 += 1) {
                  amp_lower_tmp1 = amp_lower_rho_i[-c1][c2][c3];
                  amp_lower_tmp2 = (amp_lower_tmp1 * amp_lower_tmp1);
                  amp_lower_tmp3 = (amp_lower_tmp1 * amp_lower_tmp2);
                  amp_lower_d_0[c2][c3][0][0] = (amp_lower_one + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx1) + (amp_lower_ty1 * amp_lower_dy1)) + (amp_lower_tz1 * amp_lower_dz1))));
                  amp_lower_d_1[c2][c3][1][0] = amp_lower_zero;
                  amp_lower_d_2[c2][c3][2][0] = amp_lower_zero;
                  amp_lower_d_3[c2][c3][3][0] = amp_lower_zero;
                  amp_lower_d_4[c2][c3][4][0] = amp_lower_zero;
                  amp_lower_d_0[c2][c3][0][1] = (((amp_lower_dt * amp_lower_two) * ((((-amp_lower_tx1) * amp_lower_r43) - amp_lower_ty1) - amp_lower_tz1)) * ((amp_lower_c34 * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3][1]));
                  amp_lower_d_1[c2][c3][1][1] = ((amp_lower_one + ((((amp_lower_dt * amp_lower_two) * amp_lower_c34) * amp_lower_tmp1) * (((amp_lower_tx1 * amp_lower_r43) + amp_lower_ty1) + amp_lower_tz1))) + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx2) + (amp_lower_ty1 * amp_lower_dy2)) + (amp_lower_tz1 * amp_lower_dz2))));
                  amp_lower_d_2[c2][c3][2][1] = amp_lower_zero;
                  amp_lower_d_3[c2][c3][3][1] = amp_lower_zero;
                  amp_lower_d_4[c2][c3][4][1] = amp_lower_zero;
                  amp_lower_d_0[c2][c3][0][2] = (((amp_lower_dt * amp_lower_two) * (((-amp_lower_tx1) - (amp_lower_ty1 * amp_lower_r43)) - amp_lower_tz1)) * ((amp_lower_c34 * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3][2]));
                  amp_lower_d_1[c2][c3][1][2] = amp_lower_zero;
                  amp_lower_d_2[c2][c3][2][2] = ((amp_lower_one + ((((amp_lower_dt * amp_lower_two) * amp_lower_c34) * amp_lower_tmp1) * ((amp_lower_tx1 + (amp_lower_ty1 * amp_lower_r43)) + amp_lower_tz1))) + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx3) + (amp_lower_ty1 * amp_lower_dy3)) + (amp_lower_tz1 * amp_lower_dz3))));
                  amp_lower_d_3[c2][c3][3][2] = amp_lower_zero;
                  amp_lower_d_4[c2][c3][4][2] = amp_lower_zero;
                  amp_lower_d_0[c2][c3][0][3] = (((amp_lower_dt * amp_lower_two) * (((-amp_lower_tx1) - amp_lower_ty1) - (amp_lower_tz1 * amp_lower_r43))) * ((amp_lower_c34 * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3][3]));
                  amp_lower_d_1[c2][c3][1][3] = amp_lower_zero;
                  amp_lower_d_2[c2][c3][2][3] = amp_lower_zero;
                  amp_lower_d_3[c2][c3][3][3] = ((amp_lower_one + ((((amp_lower_dt * amp_lower_two) * amp_lower_c34) * amp_lower_tmp1) * ((amp_lower_tx1 + amp_lower_ty1) + (amp_lower_tz1 * amp_lower_r43)))) + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx4) + (amp_lower_ty1 * amp_lower_dy4)) + (amp_lower_tz1 * amp_lower_dz4))));
                  amp_lower_d_4[c2][c3][4][3] = amp_lower_zero;
                  amp_lower_d_0[c2][c3][0][4] = (((-amp_lower_dt) * amp_lower_two) * ((((((((amp_lower_tx1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) + (amp_lower_ty1 * (amp_lower_c34 - amp_lower_c1345))) + (amp_lower_tz1 * (amp_lower_c34 - amp_lower_c1345))) * (amp_lower_u[-c1][c2][c3][1] * amp_lower_u[-c1][c2][c3][1])) + ((((amp_lower_tx1 * (amp_lower_c34 - amp_lower_c1345)) + (amp_lower_ty1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345))) + (amp_lower_tz1 * (amp_lower_c34 - amp_lower_c1345))) * (amp_lower_u[-c1][c2][c3][2] * amp_lower_u[-c1][c2][c3][2]))) + ((((amp_lower_tx1 * (amp_lower_c34 - amp_lower_c1345)) + (amp_lower_ty1 * (amp_lower_c34 - amp_lower_c1345))) + (amp_lower_tz1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345))) * (amp_lower_u[-c1][c2][c3][3] * amp_lower_u[-c1][c2][c3][3]))) * amp_lower_tmp3) + (((((amp_lower_tx1 + amp_lower_ty1) + amp_lower_tz1) * amp_lower_c1345) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3][4])));
                  amp_lower_d_1[c2][c3][1][4] = ((((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) + (amp_lower_ty1 * (amp_lower_c34 - amp_lower_c1345))) + (amp_lower_tz1 * (amp_lower_c34 - amp_lower_c1345)))) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3][1]);
                  amp_lower_d_2[c2][c3][2][4] = ((((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * (amp_lower_c34 - amp_lower_c1345)) + (amp_lower_ty1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345))) + (amp_lower_tz1 * (amp_lower_c34 - amp_lower_c1345)))) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3][2]);
                  amp_lower_d_3[c2][c3][3][4] = ((((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * (amp_lower_c34 - amp_lower_c1345)) + (amp_lower_ty1 * (amp_lower_c34 - amp_lower_c1345))) + (amp_lower_tz1 * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)))) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3][3]);
                  amp_lower_d_4[c2][c3][4][4] = ((amp_lower_one + ((((amp_lower_dt * amp_lower_two) * ((amp_lower_tx1 + amp_lower_ty1) + amp_lower_tz1)) * amp_lower_c1345) * amp_lower_tmp1)) + ((amp_lower_dt * amp_lower_two) * (((amp_lower_tx1 * amp_lower_dx5) + (amp_lower_ty1 * amp_lower_dy5)) + (amp_lower_tz1 * amp_lower_dz5))));
                  amp_lower_tmp1 = amp_lower_rho_i[-c1][c2][c3 + 1];
                  amp_lower_tmp2 = (amp_lower_tmp1 * amp_lower_tmp1);
                  amp_lower_tmp3 = (amp_lower_tmp1 * amp_lower_tmp2);
                  amp_lower_a_0[c2][c3][0][0] = (((-amp_lower_dt) * amp_lower_tx1) * amp_lower_dx1);
                  amp_lower_a_1[c2][c3][1][0] = (amp_lower_dt * amp_lower_tx2);
                  amp_lower_a_2[c2][c3][2][0] = amp_lower_zero;
                  amp_lower_a_3[c2][c3][3][0] = amp_lower_zero;
                  amp_lower_a_4[c2][c3][4][0] = amp_lower_zero;
                  amp_lower_a_0[c2][c3][0][1] = (((amp_lower_dt * amp_lower_tx2) * (((-(amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_tmp1)) * (amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_tmp1)) + ((amp_lower_cc2 * amp_lower_qs[-c1][c2][c3 + 1]) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * ((((-amp_lower_r43) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3 + 1][1])));
                  amp_lower_a_1[c2][c3][1][1] = ((((amp_lower_dt * amp_lower_tx2) * ((amp_lower_two - amp_lower_cc2) * (amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * ((amp_lower_r43 * amp_lower_c34) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * amp_lower_dx2));
                  amp_lower_a_2[c2][c3][2][1] = ((amp_lower_dt * amp_lower_tx2) * ((-amp_lower_cc2) * (amp_lower_u[-c1][c2][c3 + 1][2] * amp_lower_tmp1)));
                  amp_lower_a_3[c2][c3][3][1] = ((amp_lower_dt * amp_lower_tx2) * ((-amp_lower_cc2) * (amp_lower_u[-c1][c2][c3 + 1][3] * amp_lower_tmp1)));
                  amp_lower_a_4[c2][c3][4][1] = ((amp_lower_dt * amp_lower_tx2) * amp_lower_cc2);
                  amp_lower_a_0[c2][c3][0][2] = (((amp_lower_dt * amp_lower_tx2) * ((-(amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_u[-c1][c2][c3 + 1][2])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tx1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3 + 1][2])));
                  amp_lower_a_1[c2][c3][1][2] = ((amp_lower_dt * amp_lower_tx2) * (amp_lower_u[-c1][c2][c3 + 1][2] * amp_lower_tmp1));
                  amp_lower_a_2[c2][c3][2][2] = ((((amp_lower_dt * amp_lower_tx2) * (amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tx1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * amp_lower_dx3));
                  amp_lower_a_3[c2][c3][3][2] = amp_lower_zero;
                  amp_lower_a_4[c2][c3][4][2] = amp_lower_zero;
                  amp_lower_a_0[c2][c3][0][3] = (((amp_lower_dt * amp_lower_tx2) * ((-(amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_u[-c1][c2][c3 + 1][3])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tx1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3 + 1][3])));
                  amp_lower_a_1[c2][c3][1][3] = ((amp_lower_dt * amp_lower_tx2) * (amp_lower_u[-c1][c2][c3 + 1][3] * amp_lower_tmp1));
                  amp_lower_a_2[c2][c3][2][3] = amp_lower_zero;
                  amp_lower_a_3[c2][c3][3][3] = ((((amp_lower_dt * amp_lower_tx2) * (amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tx1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tx1) * amp_lower_dx4));
                  amp_lower_a_4[c2][c3][4][3] = amp_lower_zero;
                  amp_lower_a_0[c2][c3][0][4] = (((amp_lower_dt * amp_lower_tx2) * ((((amp_lower_cc2 * amp_lower_two) * amp_lower_qs[-c1][c2][c3 + 1]) - (amp_lower_cc1 * amp_lower_u[-c1][c2][c3 + 1][4])) * (amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_tmp2))) - ((amp_lower_dt * amp_lower_tx1) * ((((((-((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) * amp_lower_tmp3) * (amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_u[-c1][c2][c3 + 1][1])) - (((amp_lower_c34 - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[-c1][c2][c3 + 1][2] * amp_lower_u[-c1][c2][c3 + 1][2]))) - (((amp_lower_c34 - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[-c1][c2][c3 + 1][3] * amp_lower_u[-c1][c2][c3 + 1][3]))) - ((amp_lower_c1345 * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3 + 1][4]))));
                  amp_lower_a_1[c2][c3][1][4] = (((amp_lower_dt * amp_lower_tx2) * ((amp_lower_cc1 * (amp_lower_u[-c1][c2][c3 + 1][4] * amp_lower_tmp1)) - (amp_lower_cc2 * (((amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_u[-c1][c2][c3 + 1][1]) * amp_lower_tmp2) + (amp_lower_qs[-c1][c2][c3 + 1] * amp_lower_tmp1))))) - ((((amp_lower_dt * amp_lower_tx1) * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3 + 1][1]));
                  amp_lower_a_2[c2][c3][2][4] = (((amp_lower_dt * amp_lower_tx2) * (((-amp_lower_cc2) * (amp_lower_u[-c1][c2][c3 + 1][2] * amp_lower_u[-c1][c2][c3 + 1][1])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_tx1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3 + 1][2]));
                  amp_lower_a_3[c2][c3][3][4] = (((amp_lower_dt * amp_lower_tx2) * (((-amp_lower_cc2) * (amp_lower_u[-c1][c2][c3 + 1][3] * amp_lower_u[-c1][c2][c3 + 1][1])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_tx1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1][c2][c3 + 1][3]));
                  amp_lower_a_4[c2][c3][4][4] = ((((amp_lower_dt * amp_lower_tx2) * (amp_lower_cc1 * (amp_lower_u[-c1][c2][c3 + 1][1] * amp_lower_tmp1))) - (((amp_lower_dt * amp_lower_tx1) * amp_lower_c1345) * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tx1) * amp_lower_dx5));
                  amp_lower_tmp1 = amp_lower_rho_i[-c1][c2 + 1][c3];
                  amp_lower_tmp2 = (amp_lower_tmp1 * amp_lower_tmp1);
                  amp_lower_tmp3 = (amp_lower_tmp1 * amp_lower_tmp2);
                  amp_lower_b_0[c2][c3][0][0] = (((-amp_lower_dt) * amp_lower_ty1) * amp_lower_dy1);
                  amp_lower_b_1[c2][c3][1][0] = amp_lower_zero;
                  amp_lower_b_2[c2][c3][2][0] = (amp_lower_dt * amp_lower_ty2);
                  amp_lower_b_3[c2][c3][3][0] = amp_lower_zero;
                  amp_lower_b_4[c2][c3][4][0] = amp_lower_zero;
                  amp_lower_b_0[c2][c3][0][1] = (((amp_lower_dt * amp_lower_ty2) * ((-(amp_lower_u[-c1][c2 + 1][c3][1] * amp_lower_u[-c1][c2 + 1][c3][2])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_ty1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1][c2 + 1][c3][1])));
                  amp_lower_b_1[c2][c3][1][1] = ((((amp_lower_dt * amp_lower_ty2) * (amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_ty1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_ty1) * amp_lower_dy2));
                  amp_lower_b_2[c2][c3][2][1] = ((amp_lower_dt * amp_lower_ty2) * (amp_lower_u[-c1][c2 + 1][c3][1] * amp_lower_tmp1));
                  amp_lower_b_3[c2][c3][3][1] = amp_lower_zero;
                  amp_lower_b_4[c2][c3][4][1] = amp_lower_zero;
                  amp_lower_b_0[c2][c3][0][2] = (((amp_lower_dt * amp_lower_ty2) * (((-(amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_tmp1)) * (amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_tmp1)) + (amp_lower_cc2 * (amp_lower_qs[-c1][c2 + 1][c3] * amp_lower_tmp1)))) - ((amp_lower_dt * amp_lower_ty1) * ((((-amp_lower_r43) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1][c2 + 1][c3][2])));
                  amp_lower_b_1[c2][c3][1][2] = ((amp_lower_dt * amp_lower_ty2) * ((-amp_lower_cc2) * (amp_lower_u[-c1][c2 + 1][c3][1] * amp_lower_tmp1)));
                  amp_lower_b_2[c2][c3][2][2] = ((((amp_lower_dt * amp_lower_ty2) * ((amp_lower_two - amp_lower_cc2) * (amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_ty1) * ((amp_lower_r43 * amp_lower_c34) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_ty1) * amp_lower_dy3));
                  amp_lower_b_3[c2][c3][3][2] = ((amp_lower_dt * amp_lower_ty2) * ((-amp_lower_cc2) * (amp_lower_u[-c1][c2 + 1][c3][3] * amp_lower_tmp1)));
                  amp_lower_b_4[c2][c3][4][2] = ((amp_lower_dt * amp_lower_ty2) * amp_lower_cc2);
                  amp_lower_b_0[c2][c3][0][3] = (((amp_lower_dt * amp_lower_ty2) * ((-(amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_u[-c1][c2 + 1][c3][3])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_ty1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1][c2 + 1][c3][3])));
                  amp_lower_b_1[c2][c3][1][3] = amp_lower_zero;
                  amp_lower_b_2[c2][c3][2][3] = ((amp_lower_dt * amp_lower_ty2) * (amp_lower_u[-c1][c2 + 1][c3][3] * amp_lower_tmp1));
                  amp_lower_b_3[c2][c3][3][3] = ((((amp_lower_dt * amp_lower_ty2) * (amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_ty1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_ty1) * amp_lower_dy4));
                  amp_lower_b_4[c2][c3][4][3] = amp_lower_zero;
                  amp_lower_b_0[c2][c3][0][4] = (((amp_lower_dt * amp_lower_ty2) * ((((amp_lower_cc2 * amp_lower_two) * amp_lower_qs[-c1][c2 + 1][c3]) - (amp_lower_cc1 * amp_lower_u[-c1][c2 + 1][c3][4])) * (amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_tmp2))) - ((amp_lower_dt * amp_lower_ty1) * ((((((-(amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp3) * (amp_lower_u[-c1][c2 + 1][c3][1] * amp_lower_u[-c1][c2 + 1][c3][1])) - ((((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_u[-c1][c2 + 1][c3][2]))) - (((amp_lower_c34 - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[-c1][c2 + 1][c3][3] * amp_lower_u[-c1][c2 + 1][c3][3]))) - ((amp_lower_c1345 * amp_lower_tmp2) * amp_lower_u[-c1][c2 + 1][c3][4]))));
                  amp_lower_b_1[c2][c3][1][4] = (((amp_lower_dt * amp_lower_ty2) * (((-amp_lower_cc2) * (amp_lower_u[-c1][c2 + 1][c3][1] * amp_lower_u[-c1][c2 + 1][c3][2])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_ty1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1][c2 + 1][c3][1]));
                  amp_lower_b_2[c2][c3][2][4] = (((amp_lower_dt * amp_lower_ty2) * ((amp_lower_cc1 * (amp_lower_u[-c1][c2 + 1][c3][4] * amp_lower_tmp1)) - (amp_lower_cc2 * ((amp_lower_qs[-c1][c2 + 1][c3] * amp_lower_tmp1) + ((amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_u[-c1][c2 + 1][c3][2]) * amp_lower_tmp2))))) - ((((amp_lower_dt * amp_lower_ty1) * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1][c2 + 1][c3][2]));
                  amp_lower_b_3[c2][c3][3][4] = (((amp_lower_dt * amp_lower_ty2) * (((-amp_lower_cc2) * (amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_u[-c1][c2 + 1][c3][3])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_ty1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1][c2 + 1][c3][3]));
                  amp_lower_b_4[c2][c3][4][4] = ((((amp_lower_dt * amp_lower_ty2) * (amp_lower_cc1 * (amp_lower_u[-c1][c2 + 1][c3][2] * amp_lower_tmp1))) - (((amp_lower_dt * amp_lower_ty1) * amp_lower_c1345) * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_ty1) * amp_lower_dy5));
                  amp_lower_tmp1 = amp_lower_rho_i[-c1 + 1][c2][c3];
                  amp_lower_tmp2 = (amp_lower_tmp1 * amp_lower_tmp1);
                  amp_lower_tmp3 = (amp_lower_tmp1 * amp_lower_tmp2);
                  amp_lower_c_0[c2][c3][0][0] = (((-amp_lower_dt) * amp_lower_tz1) * amp_lower_dz1);
                  amp_lower_c_1[c2][c3][1][0] = amp_lower_zero;
                  amp_lower_c_2[c2][c3][2][0] = amp_lower_zero;
                  amp_lower_c_3[c2][c3][3][0] = (amp_lower_dt * amp_lower_tz2);
                  amp_lower_c_4[c2][c3][4][0] = amp_lower_zero;
                  amp_lower_c_0[c2][c3][0][1] = (((amp_lower_dt * amp_lower_tz2) * ((-(amp_lower_u[-c1 + 1][c2][c3][1] * amp_lower_u[-c1 + 1][c2][c3][3])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tz1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1 + 1][c2][c3][1])));
                  amp_lower_c_1[c2][c3][1][1] = ((((amp_lower_dt * amp_lower_tz2) * (amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_tmp1)) - (((amp_lower_dt * amp_lower_tz1) * amp_lower_c34) * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tz1) * amp_lower_dz2));
                  amp_lower_c_2[c2][c3][2][1] = amp_lower_zero;
                  amp_lower_c_3[c2][c3][3][1] = ((amp_lower_dt * amp_lower_tz2) * (amp_lower_u[-c1 + 1][c2][c3][1] * amp_lower_tmp1));
                  amp_lower_c_4[c2][c3][4][1] = amp_lower_zero;
                  amp_lower_c_0[c2][c3][0][2] = (((amp_lower_dt * amp_lower_tz2) * ((-(amp_lower_u[-c1 + 1][c2][c3][2] * amp_lower_u[-c1 + 1][c2][c3][3])) * amp_lower_tmp2)) - ((amp_lower_dt * amp_lower_tz1) * (((-amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1 + 1][c2][c3][2])));
                  amp_lower_c_1[c2][c3][1][2] = amp_lower_zero;
                  amp_lower_c_2[c2][c3][2][2] = ((((amp_lower_dt * amp_lower_tz2) * (amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tz1) * (amp_lower_c34 * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tz1) * amp_lower_dz3));
                  amp_lower_c_3[c2][c3][3][2] = ((amp_lower_dt * amp_lower_tz2) * (amp_lower_u[-c1 + 1][c2][c3][2] * amp_lower_tmp1));
                  amp_lower_c_4[c2][c3][4][2] = amp_lower_zero;
                  amp_lower_c_0[c2][c3][0][3] = (((amp_lower_dt * amp_lower_tz2) * (((-(amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_tmp1)) * (amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_tmp1)) + (amp_lower_cc2 * (amp_lower_qs[-c1 + 1][c2][c3] * amp_lower_tmp1)))) - ((amp_lower_dt * amp_lower_tz1) * ((((-amp_lower_r43) * amp_lower_c34) * amp_lower_tmp2) * amp_lower_u[-c1 + 1][c2][c3][3])));
                  amp_lower_c_1[c2][c3][1][3] = ((amp_lower_dt * amp_lower_tz2) * ((-amp_lower_cc2) * (amp_lower_u[-c1 + 1][c2][c3][1] * amp_lower_tmp1)));
                  amp_lower_c_2[c2][c3][2][3] = ((amp_lower_dt * amp_lower_tz2) * ((-amp_lower_cc2) * (amp_lower_u[-c1 + 1][c2][c3][2] * amp_lower_tmp1)));
                  amp_lower_c_3[c2][c3][3][3] = (((((amp_lower_dt * amp_lower_tz2) * (amp_lower_two - amp_lower_cc2)) * (amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tz1) * ((amp_lower_r43 * amp_lower_c34) * amp_lower_tmp1))) - ((amp_lower_dt * amp_lower_tz1) * amp_lower_dz4));
                  amp_lower_c_4[c2][c3][4][3] = ((amp_lower_dt * amp_lower_tz2) * amp_lower_cc2);
                  amp_lower_c_0[c2][c3][0][4] = (((amp_lower_dt * amp_lower_tz2) * ((((amp_lower_cc2 * amp_lower_two) * amp_lower_qs[-c1 + 1][c2][c3]) - (amp_lower_cc1 * amp_lower_u[-c1 + 1][c2][c3][4])) * (amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_tmp2))) - ((amp_lower_dt * amp_lower_tz1) * ((((((-(amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp3) * (amp_lower_u[-c1 + 1][c2][c3][1] * amp_lower_u[-c1 + 1][c2][c3][1])) - (((amp_lower_c34 - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[-c1 + 1][c2][c3][2] * amp_lower_u[-c1 + 1][c2][c3][2]))) - ((((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345) * amp_lower_tmp3) * (amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_u[-c1 + 1][c2][c3][3]))) - ((amp_lower_c1345 * amp_lower_tmp2) * amp_lower_u[-c1 + 1][c2][c3][4]))));
                  amp_lower_c_1[c2][c3][1][4] = (((amp_lower_dt * amp_lower_tz2) * (((-amp_lower_cc2) * (amp_lower_u[-c1 + 1][c2][c3][1] * amp_lower_u[-c1 + 1][c2][c3][3])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_tz1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1 + 1][c2][c3][1]));
                  amp_lower_c_2[c2][c3][2][4] = (((amp_lower_dt * amp_lower_tz2) * (((-amp_lower_cc2) * (amp_lower_u[-c1 + 1][c2][c3][2] * amp_lower_u[-c1 + 1][c2][c3][3])) * amp_lower_tmp2)) - ((((amp_lower_dt * amp_lower_tz1) * (amp_lower_c34 - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1 + 1][c2][c3][2]));
                  amp_lower_c_3[c2][c3][3][4] = (((amp_lower_dt * amp_lower_tz2) * ((amp_lower_cc1 * (amp_lower_u[-c1 + 1][c2][c3][4] * amp_lower_tmp1)) - (amp_lower_cc2 * ((amp_lower_qs[-c1 + 1][c2][c3] * amp_lower_tmp1) + ((amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_u[-c1 + 1][c2][c3][3]) * amp_lower_tmp2))))) - ((((amp_lower_dt * amp_lower_tz1) * ((amp_lower_r43 * amp_lower_c34) - amp_lower_c1345)) * amp_lower_tmp2) * amp_lower_u[-c1 + 1][c2][c3][3]));
                  amp_lower_c_4[c2][c3][4][4] = ((((amp_lower_dt * amp_lower_tz2) * (amp_lower_cc1 * (amp_lower_u[-c1 + 1][c2][c3][3] * amp_lower_tmp1))) - (((amp_lower_dt * amp_lower_tz1) * amp_lower_c1345) * amp_lower_tmp1)) - ((amp_lower_dt * amp_lower_tz1) * amp_lower_dz5));
                }
              timer_stop(8);
              timer_start(9);
              for (int c2 = -62; c2 < 0; c2 += 1)
                for (int c3 = -62; c3 < 0; c3 += 1)
                  for (int c4 = 0; c4 <= 4; c4 += 1)
                    amp_lower_tv[-c2][-c3][c4] = (amp_lower_omega * (((((amp_lower_c_0[-c2][-c3][0][c4] * rsd[-c1 + 1][-c2][-c3][0]) + (amp_lower_c_1[-c2][-c3][1][c4] * rsd[-c1 + 1][-c2][-c3][1])) + (amp_lower_c_2[-c2][-c3][2][c4] * rsd[-c1 + 1][-c2][-c3][2])) + (amp_lower_c_3[-c2][-c3][3][c4] * rsd[-c1 + 1][-c2][-c3][3])) + (amp_lower_c_4[-c2][-c3][4][c4] * rsd[-c1 + 1][-c2][-c3][4])));
              for (int c2 = -62; c2 < 0; c2 += 1)
                for (int c3 = -62; c3 < 0; c3 += 1) {
                  for (int c4 = 0; c4 <= 4; c4 += 1)
                    amp_lower_tv[-c2][-c3][c4] = (amp_lower_tv[-c2][-c3][c4] + (amp_lower_omega * ((((((((((amp_lower_b_0[-c2][-c3][0][c4] * rsd[-c1][-c2 + 1][-c3][0]) + (amp_lower_a_0[-c2][-c3][0][c4] * rsd[-c1][-c2][-c3 + 1][0])) + (amp_lower_b_1[-c2][-c3][1][c4] * rsd[-c1][-c2 + 1][-c3][1])) + (amp_lower_a_1[-c2][-c3][1][c4] * rsd[-c1][-c2][-c3 + 1][1])) + (amp_lower_b_2[-c2][-c3][2][c4] * rsd[-c1][-c2 + 1][-c3][2])) + (amp_lower_a_2[-c2][-c3][2][c4] * rsd[-c1][-c2][-c3 + 1][2])) + (amp_lower_b_3[-c2][-c3][3][c4] * rsd[-c1][-c2 + 1][-c3][3])) + (amp_lower_a_3[-c2][-c3][3][c4] * rsd[-c1][-c2][-c3 + 1][3])) + (amp_lower_b_4[-c2][-c3][4][c4] * rsd[-c1][-c2 + 1][-c3][4])) + (amp_lower_a_4[-c2][-c3][4][c4] * rsd[-c1][-c2][-c3 + 1][4]))));
                  for (int c4 = 0; c4 <= 4; c4 += 1) {
                    amp_lower_tmat_0[c4][0] = amp_lower_d_0[-c2][-c3][0][c4];
                    amp_lower_tmat_1[c4][1] = amp_lower_d_1[-c2][-c3][1][c4];
                    amp_lower_tmat_2[c4][2] = amp_lower_d_2[-c2][-c3][2][c4];
                    amp_lower_tmat_3[c4][3] = amp_lower_d_3[-c2][-c3][3][c4];
                    amp_lower_tmat_4[c4][4] = amp_lower_d_4[-c2][-c3][4][c4];
                  }
                  amp_lower_tmp1 = (amp_lower_one / amp_lower_tmat_0[0][0]);
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_0[1][0]);
                  amp_lower_tmat_1[1][1] = (amp_lower_tmat_1[1][1] - (amp_lower_tmp * amp_lower_tmat_1[0][1]));
                  amp_lower_tmat_2[1][2] = (amp_lower_tmat_2[1][2] - (amp_lower_tmp * amp_lower_tmat_2[0][2]));
                  amp_lower_tmat_3[1][3] = (amp_lower_tmat_3[1][3] - (amp_lower_tmp * amp_lower_tmat_3[0][3]));
                  amp_lower_tmat_4[1][4] = (amp_lower_tmat_4[1][4] - (amp_lower_tmp * amp_lower_tmat_4[0][4]));
                  amp_lower_tv[-c2][-c3][1] = (amp_lower_tv[-c2][-c3][1] - (amp_lower_tv[-c2][-c3][0] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_0[2][0]);
                  amp_lower_tmat_1[2][1] = (amp_lower_tmat_1[2][1] - (amp_lower_tmp * amp_lower_tmat_1[0][1]));
                  amp_lower_tmat_2[2][2] = (amp_lower_tmat_2[2][2] - (amp_lower_tmp * amp_lower_tmat_2[0][2]));
                  amp_lower_tmat_3[2][3] = (amp_lower_tmat_3[2][3] - (amp_lower_tmp * amp_lower_tmat_3[0][3]));
                  amp_lower_tmat_4[2][4] = (amp_lower_tmat_4[2][4] - (amp_lower_tmp * amp_lower_tmat_4[0][4]));
                  amp_lower_tv[-c2][-c3][2] = (amp_lower_tv[-c2][-c3][2] - (amp_lower_tv[-c2][-c3][0] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_0[3][0]);
                  amp_lower_tmat_1[3][1] = (amp_lower_tmat_1[3][1] - (amp_lower_tmp * amp_lower_tmat_1[0][1]));
                  amp_lower_tmat_2[3][2] = (amp_lower_tmat_2[3][2] - (amp_lower_tmp * amp_lower_tmat_2[0][2]));
                  amp_lower_tmat_3[3][3] = (amp_lower_tmat_3[3][3] - (amp_lower_tmp * amp_lower_tmat_3[0][3]));
                  amp_lower_tmat_4[3][4] = (amp_lower_tmat_4[3][4] - (amp_lower_tmp * amp_lower_tmat_4[0][4]));
                  amp_lower_tv[-c2][-c3][3] = (amp_lower_tv[-c2][-c3][3] - (amp_lower_tv[-c2][-c3][0] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_0[4][0]);
                  amp_lower_tmat_1[4][1] = (amp_lower_tmat_1[4][1] - (amp_lower_tmp * amp_lower_tmat_1[0][1]));
                  amp_lower_tmat_2[4][2] = (amp_lower_tmat_2[4][2] - (amp_lower_tmp * amp_lower_tmat_2[0][2]));
                  amp_lower_tmat_3[4][3] = (amp_lower_tmat_3[4][3] - (amp_lower_tmp * amp_lower_tmat_3[0][3]));
                  amp_lower_tmat_4[4][4] = (amp_lower_tmat_4[4][4] - (amp_lower_tmp * amp_lower_tmat_4[0][4]));
                  amp_lower_tv[-c2][-c3][4] = (amp_lower_tv[-c2][-c3][4] - (amp_lower_tv[-c2][-c3][0] * amp_lower_tmp));
                  amp_lower_tmp1 = (amp_lower_one / amp_lower_tmat_1[1][1]);
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_1[2][1]);
                  amp_lower_tmat_2[2][2] = (amp_lower_tmat_2[2][2] - (amp_lower_tmp * amp_lower_tmat_2[1][2]));
                  amp_lower_tmat_3[2][3] = (amp_lower_tmat_3[2][3] - (amp_lower_tmp * amp_lower_tmat_3[1][3]));
                  amp_lower_tmat_4[2][4] = (amp_lower_tmat_4[2][4] - (amp_lower_tmp * amp_lower_tmat_4[1][4]));
                  amp_lower_tv[-c2][-c3][2] = (amp_lower_tv[-c2][-c3][2] - (amp_lower_tv[-c2][-c3][1] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_1[3][1]);
                  amp_lower_tmat_2[3][2] = (amp_lower_tmat_2[3][2] - (amp_lower_tmp * amp_lower_tmat_2[1][2]));
                  amp_lower_tmat_3[3][3] = (amp_lower_tmat_3[3][3] - (amp_lower_tmp * amp_lower_tmat_3[1][3]));
                  amp_lower_tmat_4[3][4] = (amp_lower_tmat_4[3][4] - (amp_lower_tmp * amp_lower_tmat_4[1][4]));
                  amp_lower_tv[-c2][-c3][3] = (amp_lower_tv[-c2][-c3][3] - (amp_lower_tv[-c2][-c3][1] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_1[4][1]);
                  amp_lower_tmat_2[4][2] = (amp_lower_tmat_2[4][2] - (amp_lower_tmp * amp_lower_tmat_2[1][2]));
                  amp_lower_tmat_3[4][3] = (amp_lower_tmat_3[4][3] - (amp_lower_tmp * amp_lower_tmat_3[1][3]));
                  amp_lower_tmat_4[4][4] = (amp_lower_tmat_4[4][4] - (amp_lower_tmp * amp_lower_tmat_4[1][4]));
                  amp_lower_tv[-c2][-c3][4] = (amp_lower_tv[-c2][-c3][4] - (amp_lower_tv[-c2][-c3][1] * amp_lower_tmp));
                  amp_lower_tmp1 = (amp_lower_one / amp_lower_tmat_2[2][2]);
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_2[3][2]);
                  amp_lower_tmat_3[3][3] = (amp_lower_tmat_3[3][3] - (amp_lower_tmp * amp_lower_tmat_3[2][3]));
                  amp_lower_tmat_4[3][4] = (amp_lower_tmat_4[3][4] - (amp_lower_tmp * amp_lower_tmat_4[2][4]));
                  amp_lower_tv[-c2][-c3][3] = (amp_lower_tv[-c2][-c3][3] - (amp_lower_tv[-c2][-c3][2] * amp_lower_tmp));
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_2[4][2]);
                  amp_lower_tmat_3[4][3] = (amp_lower_tmat_3[4][3] - (amp_lower_tmp * amp_lower_tmat_3[2][3]));
                  amp_lower_tmat_4[4][4] = (amp_lower_tmat_4[4][4] - (amp_lower_tmp * amp_lower_tmat_4[2][4]));
                  amp_lower_tv[-c2][-c3][4] = (amp_lower_tv[-c2][-c3][4] - (amp_lower_tv[-c2][-c3][2] * amp_lower_tmp));
                  amp_lower_tmp1 = (amp_lower_one / amp_lower_tmat_3[3][3]);
                  amp_lower_tmp = (amp_lower_tmp1 * amp_lower_tmat_3[4][3]);
                  amp_lower_tmat_4[4][4] = (amp_lower_tmat_4[4][4] - (amp_lower_tmp * amp_lower_tmat_4[3][4]));
                  amp_lower_tv[-c2][-c3][4] = (amp_lower_tv[-c2][-c3][4] - (amp_lower_tv[-c2][-c3][3] * amp_lower_tmp));
                  amp_lower_tv[-c2][-c3][4] = (amp_lower_tv[-c2][-c3][4] / amp_lower_tmat_4[4][4]);
                  amp_lower_tv[-c2][-c3][3] = (amp_lower_tv[-c2][-c3][3] - (amp_lower_tmat_4[3][4] * amp_lower_tv[-c2][-c3][4]));
                  amp_lower_tv[-c2][-c3][3] = (amp_lower_tv[-c2][-c3][3] / amp_lower_tmat_3[3][3]);
                  amp_lower_tv[-c2][-c3][2] = ((amp_lower_tv[-c2][-c3][2] - (amp_lower_tmat_3[2][3] * amp_lower_tv[-c2][-c3][3])) - (amp_lower_tmat_4[2][4] * amp_lower_tv[-c2][-c3][4]));
                  amp_lower_tv[-c2][-c3][2] = (amp_lower_tv[-c2][-c3][2] / amp_lower_tmat_2[2][2]);
                  amp_lower_tv[-c2][-c3][1] = (((amp_lower_tv[-c2][-c3][1] - (amp_lower_tmat_2[1][2] * amp_lower_tv[-c2][-c3][2])) - (amp_lower_tmat_3[1][3] * amp_lower_tv[-c2][-c3][3])) - (amp_lower_tmat_4[1][4] * amp_lower_tv[-c2][-c3][4]));
                  amp_lower_tv[-c2][-c3][1] = (amp_lower_tv[-c2][-c3][1] / amp_lower_tmat_1[1][1]);
                  amp_lower_tv[-c2][-c3][0] = ((((amp_lower_tv[-c2][-c3][0] - (amp_lower_tmat_1[0][1] * amp_lower_tv[-c2][-c3][1])) - (amp_lower_tmat_2[0][2] * amp_lower_tv[-c2][-c3][2])) - (amp_lower_tmat_3[0][3] * amp_lower_tv[-c2][-c3][3])) - (amp_lower_tmat_4[0][4] * amp_lower_tv[-c2][-c3][4]));
                  amp_lower_tv[-c2][-c3][0] = (amp_lower_tv[-c2][-c3][0] / amp_lower_tmat_0[0][0]);
                  rsd[-c1][-c2][-c3][0] = (rsd[-c1][-c2][-c3][0] - amp_lower_tv[-c2][-c3][0]);
                  rsd[-c1][-c2][-c3][1] = (rsd[-c1][-c2][-c3][1] - amp_lower_tv[-c2][-c3][1]);
                  rsd[-c1][-c2][-c3][2] = (rsd[-c1][-c2][-c3][2] - amp_lower_tv[-c2][-c3][2]);
                  rsd[-c1][-c2][-c3][3] = (rsd[-c1][-c2][-c3][3] - amp_lower_tv[-c2][-c3][3]);
                  rsd[-c1][-c2][-c3][4] = (rsd[-c1][-c2][-c3][4] - amp_lower_tv[-c2][-c3][4]);
                }
              timer_stop(9);
            }
            timer_start(10);
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c2 = 1; c2 <= 62; c2 += 1)
                for (int c3 = 1; c3 <= 62; c3 += 1)
                  for (int c4 = 0; c4 <= 4; c4 += 1)
                    amp_lower_u[c1][c2][c3][c4] = (amp_lower_u[c1][c2][c3][c4] + (amp_lower_tmp_ssor * rsd[c1][c2][c3][c4]));
            timer_stop(10);
            if (((c0) % (inorm)) == 0) {
              if ((1)) {
                timer_start(11);
              }
              l2norm(64, 64, 64, (nx0), (ny0), (nz0), (1), (63), (1), (63), rsd, delunm);
              if ((1)) {
                timer_stop(11);
              }
            }
            rhs();
            if ((((c0) % (inorm)) == 0) || (c0 == itmax ? 1 : 0)) {
              if ((1)) {
                timer_start(11);
              }
              l2norm(64, 64, 64, (nx0), (ny0), (nz0), (1), (63), (1), (63), rsd, rsdnm);
              if ((1)) {
                timer_stop(11);
              }
            }
          }
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c2 = 1; c2 <= 62; c2 += 1)
                for (int c3 = 0; c3 <= 4; c3 += 1)
                  u[c0][c1][c2][c3] = (double)amp_lower_u[c0][c1][c2][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                d[c0][c1][4][c3] = (double)amp_lower_d_4[c0][c1][4][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                d[c0][c1][3][c3] = (double)amp_lower_d_3[c0][c1][3][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                d[c0][c1][2][c3] = (double)amp_lower_d_2[c0][c1][2][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                d[c0][c1][1][c3] = (double)amp_lower_d_1[c0][c1][1][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                d[c0][c1][0][c3] = (double)amp_lower_d_0[c0][c1][0][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                c[c0][c1][4][c3] = (double)amp_lower_c_4[c0][c1][4][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                c[c0][c1][3][c3] = (double)amp_lower_c_3[c0][c1][3][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                c[c0][c1][2][c3] = (double)amp_lower_c_2[c0][c1][2][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                c[c0][c1][1][c3] = (double)amp_lower_c_1[c0][c1][1][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                c[c0][c1][0][c3] = (double)amp_lower_c_0[c0][c1][0][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                b[c0][c1][4][c3] = (double)amp_lower_b_4[c0][c1][4][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                b[c0][c1][3][c3] = (double)amp_lower_b_3[c0][c1][3][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                b[c0][c1][2][c3] = (double)amp_lower_b_2[c0][c1][2][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                b[c0][c1][1][c3] = (double)amp_lower_b_1[c0][c1][1][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                b[c0][c1][0][c3] = (double)amp_lower_b_0[c0][c1][0][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                a[c0][c1][4][c3] = (double)amp_lower_a_4[c0][c1][4][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                a[c0][c1][3][c3] = (double)amp_lower_a_3[c0][c1][3][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                a[c0][c1][2][c3] = (double)amp_lower_a_2[c0][c1][2][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                a[c0][c1][1][c3] = (double)amp_lower_a_1[c0][c1][1][c3];
          for (int c0 = 1; c0 <= 62; c0 += 1)
            for (int c1 = 1; c1 <= 62; c1 += 1)
              for (int c3 = 0; c3 <= 4; c3 += 1)
                a[c0][c1][0][c3] = (double)amp_lower_a_0[c0][c1][0][c3];
        }
      }
    }
    #pragma endscop

    timer_stop(1);
    maxtime = timer_read(1);
}
