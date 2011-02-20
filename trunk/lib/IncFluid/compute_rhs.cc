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


/*! \file  compute_rhs.cc
 * 
 * @brief Compute rhs.
 * @sa void IncFluid::Compute_rhs() 
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "IncFluid.h"

//*********************************************************************************************

/** @brief Compute rhs for fluid simulation
 * 
 * @note Status before entering the fn
 *		nlin[i] contains Dj T[Vr[j]*Vr[i]] - f[i];
 *		F.p contains pressure;
 * 
 * @return (rhs) = nlin[i] =  -Dj T[Vr[j]*Vr[i]] + f[i] - Di [p(k)] (-grad(pressure)) 	
 */
void IncFluid::Compute_rhs()       
{	

	Add_pressure_gradient();				// nlin = nlin + grad(p)

	*nlin1 = -*nlin1;  
	*nlin2 = -*nlin2;
	*nlin3 = -*nlin3;
}



//*********************************************************************************************

/** @brief Compute rhs for fluid simulation with scalar
 * 
 * @note Status before entering the fn
 *		nlin[i] contains Dj T[Vr[j]*Vr[i]] - f[i];
 *		F.p contains pressure;
 *		T.nlin = T[U.grad T]	
 * 
 * @return (rhs) = nlin[i] =  -Dj T[Vr[j]*Vr[i]] + f[i] - Di [p(k)] (-grad(pressure)) 
 *			(T.nlin) = T.nlin = -T.nlin	
 */

void IncFluid::Compute_rhs(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Compute_rhs_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG")
			 || (globalvar_prog_kind == "NonBoussinesq"))
		Compute_rhs_RB(T);
	
}


void IncFluid::Compute_rhs_scalar(IncSF& T)       
{

	Compute_rhs();				// Fluid:  U.nlin[i] = -U.nlin[i] -FT[grad(p)]

	*T.nlin = -*T.nlin;			// For scalar: rhs = -T.nlin  
}


void IncFluid::Compute_rhs_RB(IncSF& T)       
{
	
	Compute_rhs();
	
	if (globalvar_Pr_switch != "PRZERO")  	
		*T.nlin = -*T.nlin;	
}

//*********************************************************************************************

/** @brief Compute rhs for fluid simulation with a vector
 * 
 * @note Status before entering the fn
 *		nlin[i] contains Dj T[Vr[j]*Vr[i]] - f[i];
 *		F.p contains pressure;
 *		W.nlin = T[U.grad B - B.grad U]	
 * 
 * @return (rhs) = nlin[i] =  -Dj T[Vr[j]*Vr[i]] + f[i] - Di [p(k)] (-grad(pressure)) 
 *			(W.nlin) = W.nlin[i] = -W.nlin[i]
 */
void IncFluid::Compute_rhs(IncVF& W)       
{

	Compute_rhs();					// Fluid:  U.nlin[i] = -U.nlin[i] -FT[grad(p)]
	   
	*W.nlin1 = -*W.nlin1;			// VF-W rhs[i] = -nlin[i]
	*W.nlin2 = -*W.nlin2;
	*W.nlin3 = -*W.nlin3;    
}


//*********************************************************************************************

/** @brief Compute rhs for fluid simulation with a vector and a scalar
 * 
 * @note Status before entering the fn
 *		nlin[i] contains Dj T[Vr[j]*Vr[i]] - f[i];
 *		F.p contains pressure;
 *		W.nlin = T[U.grad B - B.grad U]	
 *		T.nlin = T[U.grad T]
 * 
 * @return (rhs) = nlin[i] =  -Dj T[Vr[j]*Vr[i]] + f[i] - Di [p(k)] (-grad(pressure)) 
 *			(W.nlin) = W.nlin[i] = -W.nlin[i]
 *			(T.nlin) = T.nlin = -T.nlin	
 */
void IncFluid::Compute_rhs(IncVF& W, IncSF& T)       
{
	Compute_rhs(W);
	
	*T.nlin = -*T.nlin;			// For scalar: rhs = -T.nlin     
}


//**********************************   End of compute_rhs.cc  *********************************

