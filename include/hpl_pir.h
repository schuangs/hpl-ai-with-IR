/*
 * By Junkang Huang, August, 2020
 * 
 */
#ifndef HPL_PIR_H
#define HPL_PIR_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_misc.h"
#include "hpl_panel.h"
#include "hpl_pgesv.h"
#include "hpl_grid.h"
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_pir
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_T_pdmat *,
   HPL_T_pmat *
) );

int                              HPL_pgmres
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_pdmat *,
   const double *,
   const double *,
   const int,
   const double *,
   double *,
   double,
   const int,
   const int
) );

void                             cal_pre
STDC_ARGS( (
   HPL_T_grid *,  
   const double *,
   double *,      
   double *,      
   const int,     
   const int,
   const int,    
   const int,     
   const int     
) );

#endif