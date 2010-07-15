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

/*! \file  compute_dt.cc
 * 
 * @brief Compute dt using CFL criteria.
 *			Returns min(dx/Urms, Tdt_fixed) to avoid higher dt.
 *
 * @sa DP IncFluid::Get_dt()
 * @note	Reference Pope for CFL condition
 *
 * @author  M. K. Verma
 * @version 4.0   MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "IncFluid.h"


//*********************************************************************************************

DP IncFluid::Get_dt()
{

	CV_Compute_totalenergy_diss();
	
	DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
	DP dx = (2*M_PI)/kmax;
	DP dt = dx / (sqrt(CV_total_energy) *20);

	MPI_Bcast( &dt, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
	
	return min(dt, Tdt_fixed);
}	


//
// Scalar: Same as fluid
//

DP IncFluid::Get_dt(IncSF& T)
{
	return Get_dt();
}	

//
// Vector
//

DP IncFluid::Get_dt(IncVF& W)
{
	
	CV_Compute_totalenergy_diss();
	W.CV_Compute_totalenergy_diss();
	
	DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
	DP dx = (2*M_PI)/kmax;
	DP dt = dx / (sqrt(CV_total_energy+W.CV_total_energy) *20);
	
	if (my_id == master_id)
	{
		TinyVector<DP,3>  B0;	
		B0 = real((*W.V1)(0,0,0)), real((*W.V2)(0,0,0)), real((*W.V3)(0,0,0));
		DP B0sqr = dot(B0, B0);
		
		if (B0sqr > MYEPS)
			dt = min(dt, dx/(sqrt(B0sqr) * 20));
	}		
		
	MPI_Bcast( &dt, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
	
	return min(dt, Tdt_fixed);
}	

//
// Vector + scalar
//

DP IncFluid::Get_dt(IncVF& W, IncSF& T)
{

	return Get_dt(W);
}



//**********************************   End of compute_dt.cc  **********************************


