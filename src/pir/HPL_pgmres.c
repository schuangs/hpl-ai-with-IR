/*
 *  by Junkang Huang, Dec. 2020
 * 
 *  based on the implementation of parallel GMRES:
 *      Parallelization Of The GMRES ---by Morgan Görtz, Lund University.
 *  
 *  and the pioneer work implementing Householder Transformations into GMRES:
 *      Implementation Of The GMRES Method Using Householder Transformations Method
 *      ---by Homer F. Walker, 1988
 */
# include "hpl.h"
# include "mpi.h"

# include <math.h>


double sign(double x)
{
    return (x >= 0) ? 1 : -1;
}

/* 
 * givens_rotation():
 * 
 * 1. perform:
 *   v <- Jk-1Jk-2...J1J0v
 * 
 * 2. solve for Jk:
 *   s.t. Jkv = (v0,v1,...,vk-2,somevalue,0,...,0)
 * 
 * 3. perform:
 *   v <- Jkv
 *   w <- Jkw
 * 
 * 4. append v to R, that is:
 *   R = [R, v]
 */
void givens_rotations
(
    HPL_T_grid *                    GRID,       /* processes grid information */
    HPL_T_pdmat *                   A,          /* local A */
    double *                        v,          /* kth column of H */
    double *                        w,          /* rhs */
    double *                        R,          /* R matrix */
    double *                        sinus,      /* sin(theta) */
    double *                        cosus,      /* cos(theta) */
    const int                       k,          /* offset */
    const int                       MM          /* restart size */
)
{
    /* local variables */
    int pi, pi1, ii, ii1, i;
    double tmp;

    /* update v */
    for (i = 0; i < k; ++i)
    {
        /* calculate the process row which possess v[i] and v[i+1] */
        HPL_indxg2lp(&ii,  &pi,  i,   A->nb, A->nb, 0, GRID->nprow);
        HPL_indxg2lp(&ii1, &pi1, i+1, A->nb, A->nb, 0, GRID->nprow);

        if (pi == pi1)
        {/* if two elememts in one process row, just perform local update */
            if (GRID->myrow == pi)
            {
                /* update v */
                tmp    = cosus[i]*v[ii] - sinus[i]*v[ii1];
                v[ii1] = sinus[i]*v[ii] + cosus[i]*v[ii1];
                v[ii]  = tmp;
            }
        }
        else
        {/* if two elememts in different process row, some communication required */
            if (GRID->myrow == pi)
            {
                /* update v */
                HPL_dsend(&v[ii], 1, pi1, 0, GRID->col_comm);
                HPL_drecv(&tmp,   1, pi1, 0, GRID->col_comm);
                v[ii] = cosus[i]*v[ii] - sinus[i]*tmp;
            }
            if (GRID->myrow == pi1)
            {   
                /* update v */
                HPL_dsend(&v[ii1], 1, pi, 0, GRID->col_comm);
                HPL_drecv(&tmp,    1, pi, 0, GRID->col_comm);
                v[ii1] = sinus[i]*tmp + cosus[i]*v[ii1];
            }
        }
    }

    /* solve for Jk */

    /* calculate the process row which possess v[k] and v[k+1] */
    HPL_indxg2lp(&ii,  &pi,  k,   A->nb, A->nb, 0, GRID->nprow);
    HPL_indxg2lp(&ii1, &pi1, k+1, A->nb, A->nb, 0, GRID->nprow);

    if (pi == pi1)
    {/* if two elements in the same process row */ 
        if (GRID->myrow == pi)
        {
            /* calculate sin and cos for Jk */
            cosus[k] = v[ii]   / sqrt(v[ii]*v[ii] + v[ii1]*v[ii1]);
            sinus[k] = -v[ii1] / sqrt(v[ii]*v[ii] + v[ii1]*v[ii1]);

            /* update v */
            tmp    = cosus[k]*v[ii] - sinus[k]*v[ii1];
            v[ii1] = sinus[k]*v[ii] + cosus[k]*v[ii1];
            v[ii]  = tmp;
        }
        /* broadcast sin and cos */
        HPL_broadcast(&cosus[k], 1, HPL_DOUBLE, pi, GRID->col_comm);
        HPL_broadcast(&sinus[k], 1, HPL_DOUBLE, pi, GRID->col_comm);

        /* update w */
        tmp    = cosus[k]*w[k] - sinus[k]*w[k+1];
        w[k+1] = sinus[k]*w[k] + cosus[k]*w[k+1];
        w[k]  = tmp;
    }
    else
    {/* if two elements not in the same process row */

        /* calculate sin and cos for Jk */
        if (GRID->myrow == pi)
        {
            HPL_drecv(&tmp,   1, pi1, 1, GRID->col_comm);
            cosus[k]  = v[ii] / sqrt(v[ii]*v[ii] + tmp*tmp);
            sinus[k]  = -tmp  / sqrt(v[ii]*v[ii] + tmp*tmp);
        }
        else if (GRID->myrow == pi1)
        {
            HPL_dsend(&v[ii1], 1, pi, 1, GRID->col_comm);
        }
        /* broadcast sin and cos */
        HPL_broadcast(&cosus[k], 1, HPL_DOUBLE, pi, GRID->col_comm);
        HPL_broadcast(&sinus[k], 1, HPL_DOUBLE, pi, GRID->col_comm);

        /* update v */
        if (GRID->myrow == pi)
        {
            HPL_dsend(&v[ii], 1, pi1, 2, GRID->col_comm);
            HPL_drecv(&tmp,   1, pi1, 2, GRID->col_comm);
            v[ii] = cosus[k]*v[ii] - sinus[k]*tmp;
        }
        if (GRID->myrow == pi1)
        {   
            HPL_dsend(&v[ii1], 1, pi, 2, GRID->col_comm);
            HPL_drecv(&tmp,    1, pi, 2, GRID->col_comm);
            v[ii1] = sinus[k]*tmp + cosus[k]*v[ii1];
        }
        /* update w */
        tmp    = cosus[k]*w[k] - sinus[k]*w[k+1];
        w[k+1] = sinus[k]*w[k] + cosus[k]*w[k+1];
        w[k]   = tmp;
    }
    /* update R */
    for (i = 0; i < MM; ++i)
    {
        HPL_indxg2lp(&ii, &pi, i, A->nb, A->nb, 0, GRID->nprow);
        if (pi == 0)
        {/* if v[i] already in process row 0 */
            if (GRID->myrow == 0)
            {
                /* just perform local R update on process row 0 */
                *Mptr(R, i, k, MM) = v[ii];
            }
        }
        else
        {/* if v[i] is in another process row */
            if (GRID->myrow == pi)
            {   
                /* send v[i] to process row 0, i+3 is just a tag 
                    in case of message mismatch */
                HPL_dsend(&v[ii], 1, 0, i+3, GRID->col_comm);
            }
            else if (GRID->myrow == 0)
            {
                /* process row 0 receive v[i], and update local R */
                HPL_drecv(Mptr(R, i, k, MM), 1, pi, i+3, GRID->col_comm);
            }
        }
    }
    /* broadcast R in process row 0 to all */
    HPL_broadcast(Mptr(R, 0, k, MM), MM, HPL_DOUBLE, 0, GRID->col_comm);

    /* end of givens_rotations() */
}

/* 
 * generateHouseholder():
 *  
 * solve for Householder vector u:
 *   s.t. Pkx = (I-2uuT)x = [x0,x1,..,xk-2,alpha,0,..0], alpha != 0
 */
void generateHouseholder
(
    HPL_T_grid *                    GRID,       /* processes grid information */
    HPL_T_pdmat *                   A,          /* local A */
    const double *                  x,          /* local object vector pointer */
    double *                        u,          /* local result Householder Vector */
    const int                       k,          /* order of the Householder */
    double *                        alpha       /* result variable */
)
{
    /* local variables */
    const int myrow = GRID->myrow;
    int i, ig, mp = A->mp, pi;
    double r = 0;

    for(i = 0; i < mp; i++)
    {
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);
        if(ig >= k)
        {
            /* load u[k:] with x[k:] for further operation */
            u[i] = x[i];
            /* calculate (xk*xk) + (xk+1*xk+1) + ...*/
            r += x[i]*x[i];
        }
        else
        {
            /* u[:k] should be 0 */
            u[i] = 0;
        }
    }
    HPL_indxg2lp(&i, &pi, k, A->nb, A->nb, 0, GRID->nprow);
    /* Get the total r on process row which possess u[k] */
    HPL_reduce(&r, 1, HPL_DOUBLE, HPL_sum, pi, GRID->col_comm);

    if(myrow == pi)
    {/* perform computation on process pi */
        /* calculate alpha and r */
        *alpha = -sign(x[i]) * sqrt(r);
        r = sqrt( 0.5*((*alpha)*(*alpha)-x[i]*(*alpha)));
        /* compute the first nonzero value of the transformation vector u */
        u[i]=x[i]-*alpha;
    }
    /* send r and alpha to all process */
    HPL_broadcast(&r, 1, HPL_DOUBLE, pi, GRID->col_comm);
    HPL_broadcast(alpha, 1, HPL_DOUBLE, pi, GRID->col_comm);
    /* apply 1/2r on u for all processes */
    for(i = 0; i < mp; i++){
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);
        if(ig >= k){
            u[i] *= (1./(2.*r));
        }
    }

    /* end of generateHouseholder() */
}

/*
 * applyHouseholder()
 * 
 * perform:
 *    y = Pkx = (I-2uuT)x = x - 2u(x, u)
 */
void applyHouseholder
(
    HPL_T_grid *                    GRID,       /* processes grid information */
    HPL_T_pdmat *                   A,          /* local A */
    const double *                  x,          /* local object vector pointer */
    double *                        u,          /* local result Householder Vector */
    const int                       k,          /* order of the Householder */
    double *                        y          /* target vector */
)
{
    /* local variables */
    double segsum = 0;
    int i, ig;

    /* calculate (x, u) */
    for(i = 0; i < A->mp; i++)
    {
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);
        if(ig >= k)
        {
            segsum += u[i] * x[i];
        }
    }
    /* sum (x, u) */
    HPL_all_reduce(&segsum, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm);
    
    /* calculate y = x - 2u(x, u) */
    for(i = 0; i < A->mp; i++)
    {
        ig = HPL_indxl2g(i, A->nb, A->nb, GRID->myrow, 0, GRID->nprow);
        if(ig >= k)
        {
            y[i] = x[i] - 2 * segsum * u[i];
        }
        else
        {
            y[i] = x[i];
        }
    }

    /* end of applyHouseholder() */
}

/*
 *  HPL_pgmres():
 * 
 */
int HPL_pgmres
(
   HPL_T_grid *                     GRID,
   HPL_T_palg *                     ALGO,
   HPL_T_pdmat *                    A,          /* local A */
   const double *                   preL,       /* preconditioning matrix L-1*/
   const double *                   preU,       /* preconditioning matrix U-1*/
   const int                        rmp,        /* local lda of preU and preL */
   const double *                   b,          /* local rhs */
   double *                         x,          /* local solution vector */
   double                           TOL,        /* tolerance of residual */
   const int                        MM,         /* restart size */
   const int                        MAXIT       /* maximum # of total iteration */
)
{
    /* local variables */
    int i, j, k = 0, start, ready = 0, index, pindex;
    double norm, currenterror, tmp;
    int mp = A->mp, nq = A->nq-1;


    /* distributed storages: each process row stores a part of data */
    double * v   = (double*)malloc(mp*sizeof(double));
    double * u   = (double*)malloc(mp*sizeof(double));
    double * H   = (double*)malloc(mp*MM*sizeof(double));
    double * rhs = (double*)malloc(mp*sizeof(double));

    /* replicated storage: all processes store the whole data */
    double * cosus = (double*)malloc((MM+1)*sizeof(double));
    double * sinus = (double*)malloc((MM+1)*sizeof(double));
    double * w     = (double*)malloc((MM+1)*sizeof(double));
    double * R     = (double*)malloc(MM*MM*sizeof(double));

    /* precondition b into rhs, that is: rhs = U-1L-1b*/
    HPL_dgemv(HplColumnMajor, HplNoTrans, mp, nq, 1, preL, rmp, b, 1, 0, rhs, 1);
    HPL_all_reduce(rhs, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);
    HPL_dgemv(HplColumnMajor, HplNoTrans, mp, nq, 1, preU, rmp, rhs, 1, 0, v, 1);
    HPL_all_reduce(v, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);

    /* no initial guess so the first residual r0 is just b */
    memcpy(rhs, v, mp*sizeof(double));

    norm = 0;
    /* calculate the norm of r0 */
    for (i = 0; i < mp; ++i)
    {
        norm += v[i] * v[i];
    }
    HPL_all_reduce(&norm, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm);
    norm = sqrt(norm);

    /* store the current error */
    currenterror = norm;

    /* creating w */
    memset(w, 0, (MM+1)*sizeof(double));
    /* generate and apply the first Housholder transformation
        Householder vector stored in u */
    generateHouseholder(GRID, A, v, u, 0, &w[0]);

    /* stop if approximation is good enough, zero solution returned. */
    if(currenterror < TOL)
    {
        ready = 1;
        memset(x, 0, mp*sizeof(double));
    }

    /* ------------------------------------------------- */
    /*  Restart Iterations                               */
    /* ------------------------------------------------- */
    for(start = 0; start <= MAXIT && !ready; start++)
    {
        /* do the same as above to start the method */
        if(start)
        {
            /* there is initial guess stored in x here from last iteration */
            /* calculate v = Ax */
            HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, HPL_rone,
                  A->A, A->ld, x, 1, 0, v, 1 );
            HPL_all_reduce(v, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);

            /* preconditioning A, using u as auxiliary storage */
            HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, HPL_rone,
                  preL, rmp, v, 1, 0, u, 1 );
            HPL_all_reduce(u, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);
            memcpy(v, u, mp*sizeof(double));
            HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, HPL_rone,
                  preU, rmp, v, 1, 0, u, 1 );
            HPL_all_reduce(v, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);
            memcpy(v, u, mp*sizeof(double));
            /* v = rhs - v = rhs - Ax */
            for (i = 0; i < mp; ++i)
            {
                v[i] = rhs[i] - v[i];
            }
            /* generate P0v = [alpha,0,..,0] */
            memset(w, 0, (MM+1)*sizeof(double));
            generateHouseholder(GRID, A, v, u, 0, &w[0]);
        }
        /* ------------------------------------------------ */
        /* Householder transformations and Givens rotations */
        /* ------------------------------------------------ */
        for(k = 0; k < MM; k++)
        {
            /* store the current trasformation in H */
            memcpy(Mptr(H, 0, k, mp), u, mp*sizeof(double));

            /* load ek into v */
            memset(v, 0, mp*sizeof(double));
            HPL_indxg2lp(&index, &pindex, k, A->nb, A->nb, 0, GRID->nprow);
            if (GRID->myrow == pindex)
            {
                v[index] = 1;
            }
            /* apply the last k + 1 Householder transformations in reverse order:
                that is : v = P0P1..Pkv */
            for(i = k; i >= 0; i--)
            {
                applyHouseholder(GRID, A, v, Mptr(H, 0, i, mp), i, v);
            }

            /* calculate u = AP0P1..Pkv */
            HPL_dgemv(HplColumnMajor, HplNoTrans, mp, nq, HPL_rone, 
                A->A, A->ld, v, 1, 0, u, 1);
            HPL_all_reduce(u, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);
            memcpy(v, u, mp*sizeof(double));

            /* preconditioning A, using u as auxiliary storage */
            HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, HPL_rone,
                  preL, rmp, v, 1, 0, u, 1 );
            HPL_all_reduce(u, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);
            memcpy(v, u, mp*sizeof(double));
            HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, HPL_rone,
                  preU, rmp, v, 1, 0, u, 1 );
            HPL_all_reduce(v, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);
            memcpy(v, u, mp*sizeof(double));

            /* apply last k + 1 Householder transformations: 
                that is : v = PkPk-1...P0AP0P1...Pkek*/
            for(i = 0; i <= k; i++)
            {
                applyHouseholder(GRID, A, v, Mptr(H, 0, i, mp), i, v);
            }
            /* generate and apply the last transformation */
            if(k < A->n - 1)
            {
                generateHouseholder(GRID, A, v, u, k+1, &tmp);

                /* apply this transformation: v<-Pk+1v */
                HPL_indxg2lp(&index, &pindex, k+1, A->nb, A->nb, 0, GRID->nprow);
                if (GRID->myrow == pindex)
                {
                    v[index] = tmp;
                }
                for(int i = k+2; i < A->n; i++)
                {
                    HPL_indxg2lp(&index, &pindex, i, A->nb, A->nb, 0, GRID->nprow);
                    if (GRID->myrow == pindex)
                    {
                        v[index] = 0;
                    }
                }
            }

            /* now v is the kth column of Hm. Apply Givens rotations on v */

            tmp = 0;
            /* generate and apply Givens rotations on v and w */
            /* sinus and cosus store the previous rotation parameters */
            givens_rotations(GRID, A, v, w, R, sinus, cosus, k, MM);
            /* tmp stored the last element of w, which is the current residual */
            tmp = w[k+1];
            /* store the current error */
            currenterror = fabs(tmp);
            // printf("Err: %.10f\n", currenterror);
            /* check if the solution is good enough */
            if(fabs(tmp) < TOL)
            {
                ready = 1;
                break;
            }
        }
        /* ------------------------------------ */
        /*  Solve for new solution x            */
        /* ------------------------------------ */
        if(k == MM)
        {
            k--;
        }

        /* solve Ry = w, R is upper-tri, and w will be overwritten by solution y */
        HPL_dtrsv(HplColumnMajor, HplUpper, HplNoTrans, HplNonUnit, MM, R, MM, w, 1);

        /* calculate the new solution */
        for(i = 0; i <= k; i++)
        {
            /* set v into a unit vector, with v[i] = 1 */
            memset(v, 0, mp*sizeof(double));
            HPL_indxg2lp(&index, &pindex, i, A->nb, A->nb, 0, GRID->nprow);
            if (GRID->myrow == pindex)
            {
                v[index] = 1;
            }
            /* apply last i + 1 householder transformations in reverse order */
            for(j = i; j >= 0; j--)
            {
                applyHouseholder(GRID, A, v, Mptr(H, 0, j, mp), j, v);
            }

            /* update x: perform x += yi*vi */
            for (j = 0; j < mp; ++j)
            {
                x[j] += v[j]*w[i];
            }
        }

        /* if the error is small enough, stop. 
            otherwise another iteration will be initiated. */
        if(currenterror < TOL)
        {
            ready = 1;
        }
    }
    /* check if we have done maximum number of starts */
    if(start > MAXIT)
    {
        start = MAXIT;
    }

    if (v)      free(v);
    if (u)      free(u);
    if (H)      free(H);
    if (rhs)    free(rhs);
    if (cosus)  free(cosus);
    if (sinus)  free(sinus);
    if (w)      free(w);
    if (R)      free(R);

    /* return total number of iterations performed */
    return (start * MM + k + 1);

    /* end of pgmres() */
}