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

/*! \file energy_tr_B0.cc
 * 
 * @brief  MHD: Computes \f$ \Im( (\vec{B0} \cdot \vec{K'}) \vec{V}.\vec{W}^{*} ) \f$ over 
 *			a shell or ring. Effect of B0 to energy transfer..
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */


#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************


/// Compute for ET shells \f$ \Im( (\vec{B0} \cdot \vec{K'}) \vec{V}.\vec{W}^{*} ) \f$.
/// Done to compare them with shell-to-shell transfers.
void IncVF::Compute_shell_ET_B0(IncVF& W)
{
	
	TinyVector<DP,3> B0;
	
	if (my_id == master_id)
		B0 = real((*W.V1)(0,0,0)), real((*W.V2)(0,0,0)), real((*W.V3)(0,0,0));
	
	int data_size = 3;
	MPI_Bcast(reinterpret_cast<double*>(B0.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
	
	Shell_mult_all_imagVW_B0(basis_type, alias_switch, N, B0, *V1, *V2, *V3, 
								*W.V1, *W.V2, *W.V3, 
								*shell_radius, *energy_tr_shell_B0, kfactor);	
												
}


//*********************************************************************************************

/// Compute for ET rings \f$ \Im( (\vec{B0} \cdot \vec{K'}) \vec{V}.\vec{W}^{*} ) \f$.
/// Done to compare them with shell-to-shell transfers.
void IncVF::Compute_ring_ET_B0(IncVF& W)
{
	TinyVector<DP,3> B0;
	
	if (my_id == master_id)
		B0 = real((*W.V1)(0,0,0)), real((*W.V2)(0,0,0)), real((*W.V3)(0,0,0));
	
	int data_size = 3;
	MPI_Bcast(reinterpret_cast<double*>(B0.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
	
	Ring_mult_all_imagVW_B0(basis_type, alias_switch, N, B0, *V1, *V2, *V3, 
								*W.V1, *W.V2, *W.V3,
								*ring_shell_radius, *sector_angle_ring_tr,
								*energy_tr_ring_B0, kfactor);	
												
}



//*********************************************************************************************

/// Compute for ET cylindrical rings 
///		\f$ \Im( (\vec{B0} \cdot \vec{K'}) \vec{V}.\vec{W}^{*} ) \f$.
/// Done to compare them with shell-to-shell transfers.
void IncVF::Compute_cylinder_ring_ET_B0(IncVF& W)
{

	TinyVector<DP,3> B0;
	
	if (my_id == master_id)
		B0 = real((*W.V1)(0,0,0)), real((*W.V2)(0,0,0)), real((*W.V3)(0,0,0));
	
	int data_size = 3;
	MPI_Bcast(reinterpret_cast<double*>(B0.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
	
	
	Cyl_ring_mult_all_imagVW_B0(basis_type, alias_switch, N, B0, *V1, *V2, *V3, 
								*W.V1, *W.V2, *W.V3,
								*cylinder_shell_radius, *cylinder_kpll_array_tr,
								*energy_tr_cylinder_ring_B0, kfactor);	
													
}

//******************************** End of energy_tr_B0.cc  ************************************



