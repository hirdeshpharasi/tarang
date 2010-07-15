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
/*! \file  compute_shell_tr.cc
 * 
 * @brief  Computes shell-to-shell transfer between V/W to V/W.
 *
 *	The Giver field is filled in shell  (region A1).  
 *	We do the same for the receiver field (region A2). <BR>
 *  
 *  The shell-from-index = 1:no-shells-1.
 *   shell-mult-all() provides us shell-to-shell energy transfer from these shells to
 *	 1:no_shells.  We output shell-to-shell(1:no-shells-1, 1:no-shells).
 * 
 *	Definitions given in EnergyTr.h.  We compute the energy transfer for these regions.
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */


#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************

// Fluid
//
void IncVF::Compute_kinetic_helicity_shell_tr()
{

	(*shelltoshell_hk) = 0.0;
	
	
	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{
		Fill_shell(shell_from_index);	
		
		EnergyTr_Compute_nlin();								
		// nlin = U.grad U<	

		
		Shell_mult_vorticity_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
									*V1, *V2, *V3, *shell_radius, *temp_shell_tr, kfactor);
									// results(shell_index) in (*temp_shell_tr)(index)
								
		(*shelltoshell_hk)(shell_from_index, Range::all()) = - (*temp_shell_tr)/2;	
		
		
		
		Fill_shell_vorticity(shell_from_index);	
		
		EnergyTr_Compute_nlin();										
		// nlin = U.grad omega	
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
					   *V1, *V2, *V3, *shell_radius, *temp_shell_tr, kfactor);
		
		
		(*shelltoshell_hk)(shell_from_index, Range::all()) 
			= (*shelltoshell_hk)(shell_from_index, Range::all()) - (*temp_shell_tr)/2;
		
		
		
		
		Fill_shell(shell_from_index);	
		
		EnergyTr_Compute_nlin_vorticity_helper();										
		// nlin = omega.grad U<
		
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
								 *V1, *V2, *V3, *shell_radius, *temp_shell_tr, kfactor);
		
		(*shelltoshell_hk)(shell_from_index, Range::all()) 
			= (*shelltoshell_hk)(shell_from_index, Range::all()) + (*temp_shell_tr)/2;
			
	}

}
	


//*********************************************************************************************
//	 MHD   
//
void IncVF::Compute_magnetic_helicity_shell_tr(IncVF& W)
{
		
	(*W.shelltoshell_hk) = 0.0;
	
	
	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{
		
		Fill_shell(shell_from_index);	
		// G = U<
		
		EnergyTr_Compute_nlin(W);								
		// nlin = W.grad U<	
	
		Shell_mult_vector_potential_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
					   *W.V1, *W.V2, *W.V3, *shell_radius, *temp_shell_tr, kfactor);
		
		(*W.shelltoshell_hk)(shell_from_index, Range::all()) =  (*temp_shell_tr)/2;
		// flux_hm = W.grad U<. a>
		
		
		
		
		Fill_shell(shell_from_index, W);	
		// G = W<
		
		EnergyTr_Compute_nlin();	
		// nlin = U.grad W<
		
		Shell_mult_vector_potential_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
								*W.V1, *W.V2, *W.V3, *shell_radius, *temp_shell_tr, kfactor);
		
		(*W.shelltoshell_hk)(shell_from_index, Range::all()) 
			= (*W.shelltoshell_hk)(shell_from_index, Range::all()) - (*temp_shell_tr)/2;
		// flux_hm = -U.grad W<. a>
			
		
		
		Fill_shell(shell_from_index, W);	
		// G = W<
		
		EnergyTr_Compute_nlin_UcrossB();
		// nlin = u x W<
		
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
						*W.V1, *W.V2, *W.V3, *shell_radius, *temp_shell_tr, kfactor);
								
		(*W.shelltoshell_hk)(shell_from_index, Range::all()) 
			= (*W.shelltoshell_hk)(shell_from_index, Range::all()) + (*temp_shell_tr)/2;	
		
		

	}

}



//*********************************************************************************************
// Fluid: 2D
//
void IncVF::Compute_enstrophy_shell_tr()
{
	
	if (N[2] == 1)
	{
		(*shelltoshell_hk) = 0.0;
		
		for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
		{
			
			Fill_shell_vorticity(shell_from_index);	
			// G = i k x u(k): 2D
			
			EnergyTr_Compute_nlin();										
			// nlin = U.grad omega	
			
			Shell_mult_vorticity_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
									 *V1, *V2, *V3, *shell_radius, *temp_shell_tr, kfactor);
									// results(shell_index) in (*temp_shell_tr)(index)
			
			
			(*shelltoshell_hk)(shell_from_index, Range::all()) = - (*temp_shell_tr);
		}
	}
	
	
}

//*********************************************************************************************
//	 MHD   
//
void IncVF::Compute_magnetic_enstrophy_shell_tr(IncVF& W)
{	
	
	if (N[2] == 1)
	{
		(*W.shelltoshell_hk) = 0.0;
		
		for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
		{
			
			Fill_shell(shell_from_index, W);	
			// G = W
			
			EnergyTr_Compute_nlin_UcrossB();
			// nlin = u x W< : 2D
			
			
			Shell_mult_vector_potential_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
								*W.V1, *W.V2, *W.V3, *shell_radius, *temp_shell_tr, kfactor);
			
			(*W.shelltoshell_hk)(shell_from_index, Range::all()) =  (*temp_shell_tr);
		}	
	}
	
}




//******************************  End of Compute_shell_tr.cc  *********************************




