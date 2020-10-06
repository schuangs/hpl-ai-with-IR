
# include "hpl.h"
# include "mpi.h"


double sign(double x)
{
    return (x >= 0) ? 1 : -1;
}

/* 
 * givens_rotation():
 * 
 * performing:
 *   v <- Jk-1Jk-2...J1J0v
 * 
 * solve:
 *   Jk,  s.t. Jkv = (v0,v1,...,vk-2,somevalue,0,...,0)
 * 
 * perform:
 *   v <- Jkv
 *   w <- Jkw
 * 
 * append v to R, that is:
 *   R = [R, v]
 */
void givens_rotations
(
    HPL_T_grid *                    GRID,
    HPL_T_pdmat *                   A,          /* local A */
    double *                        v,          /* kth column of H */
    double *                        w,          /* rhs */
    double *                        R,          /* R matrix */
    double *                        sinus,      /* sin(theta) */
    double *                        cosinus,    /* cos(theta) */
    const int                       k           /* offset */
)
{
    int pi, pi1, ii, ii1;
    double tmp;

    /* update v */
    for (int i = 0; i < k; ++i)
    {
        /* calculate the process row which possess v[i] and v[i+1] */
        HPL_indxg2lp(&ii,  &pi,  i,   A->nb, A->nb, 0, GRID->nprow);
        HPL_indxg2lp(&ii1, &pi1, i+1, A->nb, A->nb, 0, GRID->nprow);

        if (pi == pi1)
        {/* if two elememts in one process row, just perform local update */
            if (GRID->myrow == pi)
            {
                /* update v */
                tmp    = cosinus[i]*v[ii] - sinus[i]*v[ii1];
                v[ii1] = sinus[i]*v[ii]   + cosinus[i]*v[ii1];
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
                v[ii] = cosinus[i]*v[ii] - sinus[i]*tmp;
            }
            if (GRID->myrow == pi1)
            {   
                /* update v */
                HPL_dsend(&v[ii1], 1, pi, 0, GRID->col_comm);
                HPL_drecv(&tmp,    1, pi, 0, GRID->col_comm);
                v[ii1] = sinus[i]*tmp + cosinus[i]*v[ii1];
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
            cosinus[k] = v[ii]   / sqrt(v[ii]*v[ii] + v[ii1]*v[ii1]);
            sinus  [k] = -v[ii1] / sqrt(v[ii]*v[ii] + v[ii1]*v[ii1]);

            /* update v */
            tmp    = cosinus[k]*v[ii] - sinus[k]*v[ii1];
            v[ii1] = sinus[k]*v[ii]   + cosinus[k]*v[ii1];
            v[ii]  = tmp;
            /* update w */
            tmp    = cosinus[k]*w[ii] - sinus[k]*w[ii1];
            w[ii1] = sinus[k]*w[ii]   + cosinus[k]*w[ii1];
            w[ii]  = tmp;
        }
    }
    else
    {/* if two elements not in the same process row */

        /* calculate sin and cos for Jk */
        if (GRID->myrow == pi)
        {
            HPL_drecv(&tmp,   1, pi1, 1, GRID->col_comm);
            cosinus[k]  = v[ii]   / sqrt(v[ii]*v[ii] + tmp*tmp);
            sinus[k]    = -tmp    / sqrt(v[ii]*v[ii] + tmp*tmp);
        }
        else if (GRID->myrow == pi1)
        {
            HPL_dsend(&v[ii1], 1, pi, 1, GRID->col_comm);
        }
        HPL_broadcast(&cosinus[k], 1, HPL_DOUBLE, pi, GRID->col_comm);
        HPL_broadcast(&sinus[k],   1, HPL_DOUBLE, pi, GRID->col_comm);

        /* update v and w */
        if (GRID->myrow == pi)
        {
            /* update v */
            HPL_dsend(&v[ii], 1, pi1, 2, GRID->col_comm);
            HPL_drecv(&tmp,   1, pi1, 2, GRID->col_comm);
            v[ii] = cosinus[k]*v[ii] - sinus[k]*tmp;
            /* update w */
            HPL_dsend(&w[ii], 1, pi1, 2, GRID->col_comm);
            HPL_drecv(&tmp,   1, pi1, 2, GRID->col_comm);
            w[ii] = cosinus[k]*w[ii] - sinus[k]*tmp;
        }
        if (GRID->myrow == pi1)
        {   
            /* update v */
            HPL_dsend(&v[ii1], 1, pi, 2, GRID->col_comm);
            HPL_drecv(&tmp,    1, pi, 2, GRID->col_comm);
            v[ii1] = sinus[k]*tmp + cosinus[k]*v[ii1];
            /* update w */
            HPL_dsend(&w[ii1], 1, pi, 2, GRID->col_comm);
            HPL_drecv(&tmp,    1, pi, 2, GRID->col_comm);
            w[ii1] = sinus[k]*tmp + cosinus[k]*w[ii1];
        }

        /* update R */
        memcpy(Mptr(R, 0, k, A->n), v, A->mp);
    }
}

void generateHouseholder
(
    HPL_T_grid *                    GRID,
    const double *                  x,          /* global target vector pointer */
    double *                        u,          /* global result Householder Vector */
    const int                       k,          /* order of the Householder */
    double *                        alpha,      /* result variable */
    const int                       li,         /* local start index of x */
    const int                       ui,         /* local end index of x */
)
{
    const int myrow = GRID->myrow;
    double segsum = 0;

    for(int i = li; i < ui; i++)
    {
        if(i >= k)
        {
            //transfer the input vector to the output
            u[i] = x[i];
            //Calculate the sum part of alpha
            segsum += x[i]*x[i];
        }
    }
    double r=0;
    //Get the result on process row 0
    HPL_reduce(&segsum, &r, 1, HPL_DOUBLE, HPL_sum, 0, GRID->col_comm);
    //do the final calculations
    if(myrow == 0)
    {
        //Generate alpha and r
        *alpha = -sign(x[k]) * sqrt(r);
        r = sqrt( 0.5*((*alpha)*(*alpha)-x[k]*(*alpha)));
        //get the first value of the transformation vector
        u[k]=x[k]-*alpha;
    }else
    {
        *alpha=0;
    }
    //Send r to all process:
    HPL_broadcast(&r, 1, HPL_DOUBLE, 0, GRID->col_comm);
    //and apply 1/2r:
    for(int i = li; i < ui; i++){
        if(i >= k){
            u[i] *= (1./(2.*r));
        }
    }
}

void applyHouseholder
(
    HPL_T_grid *                    GRID,
    const double *                  u,
    const double *                  x,
    double *                        y,
    const int                       k,
    const int                       li,         /* local start index of x */
    const int                       ui,         /* local end index of x */
)
{
    double segsum = 0;
    for(int i = li; i < ui; i++)
    {
        if(i >= k)
        {
            segsum += u[i] * x[i];
        }
    }
    double tc = 0;
    //Broadcast the sum
    HPL_all_reduce(&segsum, &tc, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm);
    //The changed part:
    for(int i = li; i < x.ui; i++)
    {
        if(i >= k)
        {
            y[i]=x[i]-2*tc*u[i];
        }
    }
}

int pgmres
(
   HPL_T_grid *                     GRID,
   HPL_T_palg *                     ALGO,
   HPL_T_pdmat *                    A,          /* local A */
   const double *                   b,          /* local rhs */
   double *                         x,          /* local solution vector */
   const double *                   FACT,       /* local LU factors of A */
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
    double * v       = (double*)malloc(mp*sizeof(double));
    double * u       = (double*)malloc(mp*sizeof(double));
    double * H       = (double*)malloc(mp*MM*sizeof(double));
    double * R       = (double*)malloc(mp*MM*sizeof(double));

    /* replicated storage: all processes store the whole data */
    double * cosinus = (double*)malloc((MM+1)*sizeof(double));
    double * sinus   = (double*)malloc((MM+1)*sizeof(double));
    double * w       = (double*)malloc((MM+1)*sizeof(double));

    // doubleMatrix R(MM,MM);
    // mpi_doubleMatrix H(MM,n);
    // doubleVector cosinus(MM+1);
    // doubleVector sinus(MM+1);
    //Create a serial linear solver instance
    // serialsolver <doubleVector,doubleMatrix> Solver;
    //Do the first multiplication
    // v = (*Matrix) * (*approx);
    // v = (*rhs) - v;

    /* no initial guess so the first residual r0 is just b */
    memcpy(v, b, mp*sizeof(double));

    norm = 0;
    /* calculate the norm of r0 */
    for (i = 0; i < mp; ++i)
    {
        norm += v[i] * v[i];
    }
    HPL_all_reduce(&norm, 1, HPL_DOUBLE, HPL_sum, GRID->col_comm);
    norm = sqrt(norm);

    //store the current error
    currenterror = norm;

    /* Generate and apply the first Housholder transformation
        Housedolder vector stored in u */
    generateHouseholder(GRID, v, u, 0, &w[0]);

    /* creating w */
    for(i = 0; i < MM + 1; i++)
    {
        if(i > 0)
        {
            w[i] = 0;
        }
    }
    //Stop if approximation is good from the start
    if(currenterror < TOL)
    {
        ready = 1;
    }
    //-------------------------------------------------
    //Number of starts
    //-------------------------------------------------
    for(start = 0; start <= MAXIT && !ready; start++)
    {
        //Do the same as above to start the method
        if(start)
        {
            // v = (*Matrix)*(*approx);
            // v = (*rhs) - v;

            /* there is initial guess here from last iteration */
            memcpy(v, b, mp*sizeof(double));
            HPL_dgemv( HplColumnMajor, HplNoTrans, mp, nq, -HPL_rone,
                  A->A, A->ld, x, 1, HPL_rone, v, 1 );
            HPL_all_reduce(v, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);
            generateHouseholder(GRID, v, u, 0, &w[0]);
            for(int i = 0; i < MM + 1; i++)
            {
                if(i > 0)
                {
                    w[i] = 0;
                }
            }
        }
        //------------------------------------
        //The GMRES Method
        //------------------------------------
        for(k = 0; k < MM; k++)
        {
            //store the current trasformation in H
            memcpy(Mptr(H, 0, k, mp), u, mp*sizeof(double));
            // H[k] = u;

            /* calculating Pk..P1AP1...Pkek */
            /* load ek into v */
            memset(v, 0, mp*sizeof(double));
            // for(j = v.li();j < v.ui(); j++)
            // {
            //     v[j] = 0;
            // }
            HPL_indxg2lp(&index, &pindex, k, A->nb, A->nb, 0, GRID->nprow);
            if (GRID->myrow == pindex)
            {
                v[index] = 1;
            }
            /* Apply the last k + 1 Householder transformations in reverse order: */
            for(i = k; i >= 0; i--)
            {
                applyHouseholder(Mptr(H, 0, i, mp), v, v, i);
            }

            HPL_dgemv(HplColumnMajor, HplNoTrans, mp, nq, HPL_rone, 
                A->A, A->ld, v, 1, 0, u, 1);
            HPL_all_reduce(u, mp, HPL_DOUBLE, HPL_sum, GRID->row_comm);
            memcpy(v, u, mp*sizeof(double));
            //Apply last k + 1 Householder transformations:
            for(i = 0; i <= k; i++)
            {
                applyHouseholder(Mptr(H, 0, i, mp), v, v, i);
            }
            //Generate and apply the last transformation
            if(k < A->n - 1)
            {
                //Let u be the zero vector
                memset(u, 0, mp*sizeof(double));
                generateHouseholder(v, u, k+1, &tmp);

                /* Apply this transformation: v<-Pk+1v */
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

            /* Givens Rotation */

            tmp = 0;
            //Generate and apply the givens rotations on v and w
            /* sinus and cosinus store the previous rotation parameters */
            givens_rotations(v, w, R, sinus, cosinus, k);
            /* tmp stored the last element of w, which is the current residual */
            tmp = w[k+1];
            //store the current error
            currenterror = fabs(tmp);
            //Check if the solution is good enough
            if(fabs(tmp) < TOL)
            {
                ready = 1;
                break;
            }
        }
        //------------------------------------
        //The Solver
        //------------------------------------
        if(k == MM)
        {
            k--;
        }
        //Solve the triangular system and transfer it to u
        if(id == 0)
        {
            doubleVector x = Solver.solveUpperTriangular(R, w, k+1);
            //transfer the solution to u
            for(i = 0; i<k+1; i++){
                u[i] = x[i];
            }
        }
        //Calculate the new approximation
        for(i = 0; i <= k; i++)
        {
            //Unit vector
            for(j = v.li();j < v.ui();j++)
            {
                v[j] = 0;
            }
            v[i] = 1;
            //Apply last i + 1 householder transformations in reverse order:
            for(j = i; j >= 0; j--)
            {
                applyHouseholder(H[j], v, v, j);
            }
            //Get the coeficiant we wish to multiply the vectors with
            double c = 0;
            c = u[i];
            //make sure all parts has it
            MPI_Bcast(&c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            v *= c; //multiply with scalar
            (*approx) += v; //Get the new approximation
        }
        //If the error is small enough stop, else continue
        if(currenterror < stopping_res)
        {
            ready = 1;
        }
    }
    //Check if we have done maximum number of starts
    if(start > MAXIT)
    {
        start = MAXIT;
    }
    return (start * MM + k + 1);
}