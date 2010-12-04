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


/*! \file  compute_ring_tr.cc
 * 
 * @brief  Computes spherical ring-to-ring transfer between V/W to V/W.
 *
 *	The Giver field is filled in a ring  (region A1).  
 *	We do the same for the receiver a ring (region A2). <BR>
 *
 *	We exclude the last shell (Rmax,infty) for both Giver and Receiver rings. This is to
 *	save computer time.
 * 
 *	Definitions given in EnergyTr.h.  We compute the energy transfer for these regions.
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug	 No known bugs
 */

/**********************************************************************************************
 *
 *	No transfer for the last ring (max-possible inner radius, INF_RADIUS) 				
 * 
 *********************************************************************************************/  

#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************

void IncVF::Compute_ring_tr()
{
	
	(*ring_to_ring_self) = 0.0;
	
	
	// skip the last shell -- outer rad = infty			
	for (int ring_shell_from_i = 1; ring_shell_from_i < no_ring_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= no_sectors_ring_tr; sector_from_i++)
		{	
				
			Fill_ring(ring_shell_from_i, sector_from_i);
		
			EnergyTr_Compute_nlin();												
			// nlin = U.grad Um	

		
			Ring_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, *V1, *V2, *V3, 
							*ring_shell_radius, *sector_angle_ring_tr, *temp_ring_tr, kfactor);
			
															
			(*ring_to_ring_self)(ring_shell_from_i, sector_from_i, Range::all(), Range::all()) 
							= -*temp_ring_tr;	
										
		}
}



//*********************************************************************************************
//	SCALAR   
//

void IncVF::Compute_ring_tr(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Compute_ring_tr_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Compute_ring_tr_RB(T);
}


void IncVF::Compute_ring_tr_scalar(IncSF& T)
{
	// U to U
	Compute_ring_tr();
	
	// T to T
	(*ring_to_ring_SF) = 0.0;
	
	
	// skip the last shell -- outer rad = infty				
	for (int ring_shell_from_i = 1; ring_shell_from_i < no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= no_sectors_ring_tr; sector_from_i++)
		{
						
			Fill_ring(ring_shell_from_i, sector_from_i, T);
		
			EnergyTr_Compute_nlin(T);												
			// nlin = U.grad Tm	
	
			Ring_mult_all(basis_type, alias_switch, N, *nlin1, *T.F, *ring_shell_radius, 
										*sector_angle_ring_tr, *temp_ring_tr, kfactor);	
											
			(*ring_to_ring_SF)(ring_shell_from_i, sector_from_i, 
										Range::all(), Range::all()) = -*temp_ring_tr;	
						
		}	
}


// RBC
void IncVF::Compute_ring_tr_RB(IncSF& T)
{
	
	(*ring_to_ring_self) = 0.0;	
	(*ring_to_ring_SF) = 0.0;
	
	if (globalvar_Pr_switch == "PRZERO") 
		Compute_flux();
	
	else if (globalvar_Pr_switch == "PRINFTY")		// fill only Temperature flux
	{
		
		// skip the last shell -- outer rad = infty				
		for (int ring_shell_from_i = 1; ring_shell_from_i < no_shells; ring_shell_from_i++) 
			for (int sector_from_i = 1; sector_from_i <= no_sectors_ring_tr; sector_from_i++)
			{
				
				Fill_ring(ring_shell_from_i, sector_from_i, T);	
				
				EnergyTr_Compute_nlin(T);									
				// nlin = U.grad Tm	
				
				
				Ring_mult_all(basis_type, alias_switch, N, *nlin1, *T.F, *ring_shell_radius, 
							  *sector_angle_ring_tr, *temp_ring_tr, kfactor);	
				
				(*ring_to_ring_SF)(ring_shell_from_i, sector_from_i, 
								   Range::all(), Range::all()) = -*temp_ring_tr;	
				
			}
	}
	
	else
		Compute_ring_tr_scalar(T);
}


//*********************************************************************************************
//	 MHD   
//
void IncVF::Compute_ring_tr(IncVF& W)
{
	Compute_ring_tr();
	
	
	// W to W
	(*W.ring_to_ring_self) = 0.0;

	
	for (int ring_shell_from_i = 1; ring_shell_from_i < no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= no_sectors_ring_tr; sector_from_i++)
		{	
		
			Fill_ring(ring_shell_from_i, sector_from_i, W);	
			
			EnergyTr_Compute_nlin();											
			// nlin = U.grad Wm	
			
			Ring_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3,
										*W.V1, *W.V2, *W.V3, 
										*ring_shell_radius, *sector_angle_ring_tr, 
										*temp_ring_tr, kfactor);	
											
			(*W.ring_to_ring_self)(ring_shell_from_i, sector_from_i, 
										Range::all(), Range::all()) = -*temp_ring_tr;	
			
		}						

	
	
	//
	// U to W
	(*ring_to_ring_VF) = 0.0;
	
	
	for (int ring_shell_from_i = 1; ring_shell_from_i < no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= no_sectors_ring_tr; sector_from_i++)
		{	
				
			Fill_ring(ring_shell_from_i, sector_from_i);
	
			EnergyTr_Compute_nlin(W);											// nlin = W.grad Um	
			
			
			Ring_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
										*W.V1, *W.V2, *W.V3, 
										*ring_shell_radius, *sector_angle_ring_tr, 
										*temp_ring_tr, kfactor);	
											
			(*ring_to_ring_VF)(ring_shell_from_i, sector_from_i, 
										Range::all(), Range::all()) = *temp_ring_tr;	
			
		}
	
	
	// Shell_to_shell transfers for Elsasser vars
	//
	
	(*ring_to_ring_Elsasser) = 0.0;
	(*W.ring_to_ring_Elsasser) = 0.0;

	
	IncVF::UB_to_Elsasser_field(W);											
	// U=Zp=(U+B); B=Zm=(U-B);
	
	// ring_to_ring: Zp to Zp
	for (int ring_shell_from_i = 1; ring_shell_from_i < no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= no_sectors_ring_tr; sector_from_i++)
		{	
 	
			Fill_ring(ring_shell_from_i, sector_from_i);				
			// G = Zp<									
		
			EnergyTr_Compute_nlin(W);											
			// nlin = Zm.grad Zp<	

			Ring_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, *V1, *V2, *V3, 
										*ring_shell_radius, *sector_angle_ring_tr, 
										*temp_ring_tr, kfactor);		
										
			(*ring_to_ring_Elsasser)(ring_shell_from_i, sector_from_i, 
										Range::all(), Range::all()) = -*temp_ring_tr;	

		}
	
	// ring_to_ring: Zm to Zm	
	for (int ring_shell_from_i = 1; ring_shell_from_i < no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= no_sectors_ring_tr; sector_from_i++)
		{	
	
			Fill_ring(ring_shell_from_i, sector_from_i, W);				
			// G = Zm<
		
			EnergyTr_Compute_nlin();											
			// nlin = Zp.grad Zm<	
	
		
			Ring_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
										*W.V1, *W.V2, *W.V3, 
										*ring_shell_radius, *sector_angle_ring_tr, 
										*temp_ring_tr, kfactor);	
												
			(*W.ring_to_ring_Elsasser)(ring_shell_from_i, sector_from_i, 
										Range::all(), Range::all()) = -*temp_ring_tr;	
			
		}


	IncVF::Elsasser_to_UB_field(W);													
	// Back to U, B vars
	
}



//*********************************************************************************************
// Vector + Scalar
//
void IncVF::Compute_ring_tr(IncVF& W, IncSF& T)
{
	
	// U/W to U/W
	Compute_ring_tr(W);
	
	
	// T to T
	(*ring_to_ring_SF) = 0.0;
	

	
	// skip the last shell -- outer rad = infty				
	for (int ring_shell_from_i = 1; ring_shell_from_i < no_shells; ring_shell_from_i++) 
		for (int sector_from_i = 1; sector_from_i <= no_sectors_ring_tr; sector_from_i++)
		{
						
			Fill_ring(ring_shell_from_i, sector_from_i, T);	
		
			EnergyTr_Compute_nlin(T);												
			// nlin = U.grad Tm	

		
			Ring_mult_all(basis_type, alias_switch, N, *nlin1, *T.F, *ring_shell_radius, 
										*sector_angle_ring_tr, *temp_ring_tr, kfactor);	
											
			(*ring_to_ring_SF)(ring_shell_from_i, sector_from_i, 
										Range::all(), Range::all()) = -*temp_ring_tr;	

		}
}		


//****************************** End of Compute_ring_tr.cc ************************************


