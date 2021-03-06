/* Tarang-4.0
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-4.0 .
 *
 * Tarang-4.0 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-4.0 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-4.0; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */

/*! \file	basis_basicfn.h
 * 
 * @brief basic functions declarations that are common to all basis functions.
 *
 *	shell(n) = \f$  K' \in (R^{sh}(n-1), R^{sh}(n) ] \f$.  
 *	with typical \f$ R^{sh} = 0,2,4,8,..., R_{max},\infty \f$.
 *
 *	ring(n,m) = \f$ K' = (R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ] \f$.
 *	with typical \f$ R^{sh} = 0,4,8,..., R_{max},\infty \f$ and
 *		\f$ \Theta = 0:\Theta_{max}\f$ divided in m intervals.  \f$\Theta_{max} \f$ could be
 *		\f$ \pi \f$ or \f$ \pi/2 \f$ (see function Get_max_anis_theta).
 *
 * @version 4.0 Parallel version
 * @author  M. K. Verma
 * @date	Sept 2008
 * @bug		Precision=8 for double. Put condition for float.
 */ 

#ifndef _EXTERN_DEF_VARS_H
#define _EXTERN_DEF_VARS_H

#include <fftw3-mpi.h>
 
#include "blitz_random.h"			// for random variables - based on blitz lib

using namespace blitz ;
 
//*********************************************************************************************
// DEF vars  
#ifndef DP
#define DP double
#endif

#ifndef complx
#define complx  complex<DP>
#endif

   
//*********************************************************************************************
// Constant declarations
/// I = sqrt(-1).
const complx  I = complx(0,1.0);		

/// minusI = -sqrt(-1).			
const complx  minusI = complx(0,-1.0);

/// minus2I = -2*sqrt(-1).			
const complx  minus2I = complx(0,-2.0);	

// Number of digits for output files
// NOTE:  for double only.. For float put MY_PRECISION = 6
const int MY_PRECISION = 12;


const DP  MYEPS = 1E-15;
const DP  MYEPS2 = 1E-5;

const int MY_MAX_INT = 5000;
// cut off while reading diagnostic_procedure() array from input file and similar ops

/// Infinite radius.. All the modes outside -- for flux and shelltr calc
const DP INF_RADIUS = 10000; 

   
//*********************************************************************************************    
// extern vars
// MPI extern variables 

extern int my_id;						/**<	External: My process id	*/
extern int numprocs;					/**<	External: No of processors	*/
extern const int master_id;				/**<	External: Id of master proc	*/

extern ptrdiff_t local_N1;				/**<	External: local_N1 = N1/numprocs in the 
																	currentproc */
extern ptrdiff_t local_N1_start;		/**<	External: start of i1 in the currentproc */

extern ptrdiff_t local_N2;				/**<	External: local_N2 = N2/numprocs in the 
														currentproc; for transposed array */
extern ptrdiff_t local_N2_start;		/**<	External: start of i2 in the currentproc; 
														for transposed array */

extern MPI_Status status;

// FFT extern vars
// for fftw_original
extern fftw_plan r2c_plan_FOUR, c2r_plan_FOUR;

// for split fftw
extern fftw_plan r2c_2d_plan_FOUR, c2r_2d_plan_FOUR;
extern fftw_plan c2c_1d_forward_plan_FOUR, c2c_1d_inverse_plan_FOUR;

extern fftw_plan r2c_plan_SCFT, c2r_plan_SCFT, sintr_plan_SCFT, costr_plan_SCFT, 
					isintr_plan_SCFT, icostr_plan_SCFT; 

extern fftw_plan r2c_1d_plan_SCFT, c2r_1d_plan_SCFT;

extern int  globalvar_fftw_original_switch;
extern int	globalvar_anisotropy_switch;			// 1,2,3 for x,y,z directions
extern int	globalvar_waveno_switch;				// 0 for actual (default), 1 for grid

extern string	globalvar_prog_kind;				// program kind like INC_FLUID

// for RB
extern string	globalvar_Pr_switch;				// Prandtl number switch (PRLARGNE..) for RB
extern string	globalvar_RB_Uscaling;				// UBscaling (ULARGE... ) for RB
extern DP		globalvar_Ra;						// Rayleigh number
extern DP		globalvar_r;						// normalized Rayleigh number
extern DP		globalvar_Pr;						// Prandtl number
extern DP		globalvar_temperature_grad;			// +1 for convection; -1 for stratification; 
													// factor for u3 in temperature eqn
extern DP		globalvar_alpha_DT;	


// random extern vars
					
extern Uniform<DP> SPECrand;


#endif

//***********************************  End of extern_vars.h  ********************************




