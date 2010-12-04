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

/*! \file  compute_force_feed.cc
 * 
 * @brief Compute energy transferred from the force f to the velocity field Re(f.V*)
 * @sa void IncVF::Compute_force_feed_shell()
 *
 * @note Compute_force_feed() <BR>
 *		-- Computes energy transferred from the force f to the velocity field Re(f.V*) <BR>
 *	
 *	Compute_force_feed(T) <BR>
 *		-- Computes energy transferred from force to the velocity field  Re(f.V*), and from
 *		   scalar-forcing fT to scalar-field Re(fT.T*) <BR>
 *		   
 *	Compute_force_feed(W) <BR>
 *		-- Computes energy transferred from force to the velocity field  Re(f.V*), and 
 *		   vector-forcing fW to vector-field Re(fW.W*)	  <BR>
 *		   
 *	Compute_force_feed(W, T) <BR>
 *		-- Computes energy transferred from force to the velocity field  Re(f.V*), from
 *		   vector-forcing fW to vector-field Re(fW.W*), and from
 *		   scalar-forcing fT to scalar-field Re(fT.T*) <BR>
 *		   
 *	Compute_force_feed(T, Pr_switch) <BR>
 *	    --  If Pr = 0, transfer function is not computed.  <BR>
 *		--  If Pr != 0, same as Compute_force_fee(T)  <BR>
 *		
 *	Compute_force_feed(W, T, Pr_switch)<BR>
 *	    --  If Pr = 0, transfer function is not computed.  <BR>
 *		--  If Pr != 0, same as Compute_force_fee(W, T)	<BR>
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bug
 */
  

#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************


void IncVF::Compute_force_feed_shell()
{

	Shell_mult_all(basis_type, alias_switch, N, *Force1, *Force2, *Force3, *V1, *V2, *V3, 
								*shell_radius, *forceV_shell, kfactor);
												
}


//*********************************************************************************************

void IncVF::Compute_force_feed_shell(IncSF& T)
{

	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Compute_force_feed_shell_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Compute_force_feed_shell_RB(T);	
}



void IncVF::Compute_force_feed_shell_scalar(IncSF& T)
{
	Compute_force_feed_shell();
	
	Shell_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, *shell_radius, 
									*forceSF_shell, kfactor);																
}

//
// RB Convection
//

void IncVF::Compute_force_feed_shell_RB(IncSF& T)
{

	if (globalvar_Pr_switch == "PRZERO")
	{
		Compute_force_feed_ring();
		
		*forceSF_shell = 0.0;
	}
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
		Shell_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, *shell_radius, 
					   *forceSF_shell, kfactor);
		
		*forceV_shell = 0.0;
	}
	
	else
		Compute_force_feed_shell_scalar(T);
		
}

//*********************************************************************************************
void IncVF::Compute_force_feed_shell(IncVF& W)
{
	
	Compute_force_feed_shell();
	
	
	Shell_mult_all(basis_type, alias_switch, N, *W.Force1, *W.Force2, *W.Force3, 
							*W.V1, *W.V2, *W.V3, *shell_radius, *W.forceV_shell, kfactor);
											
}

//*********************************************************************************************

void IncVF::Compute_force_feed_shell(IncVF& W, IncSF& T)
{
	
	Compute_force_feed_shell(W);										
	// Re(f.V*), Re(fw.W*)
	
		
	Shell_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, *shell_radius, 
							*forceSF_shell, kfactor);		

}


/**********************************************************************************************

						FOR THE RINGS
	
***********************************************************************************************/

void IncVF::Compute_force_feed_ring()
{
	
	Ring_mult_all(basis_type, alias_switch, N, *Force1, *Force2, *Force3, *V1, *V2, *V3, 
					*ring_shell_radius, *sector_angle_ring_tr, *forceV_ring, kfactor);		
					
			
}


void IncVF::Compute_force_feed_ring(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Compute_force_feed_ring_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Compute_force_feed_ring_RB(T);
}


void IncVF::Compute_force_feed_ring_scalar(IncSF& T)
{
	
	Compute_force_feed_ring();

	
	Ring_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, *ring_shell_radius, 
								*sector_angle_ring_tr, *forceSF_ring, kfactor);		
}

//
// RB Convection
//

void IncVF::Compute_force_feed_ring_RB(IncSF& T)
{
	if (globalvar_Pr_switch == "PRZERO")
	{
		Compute_force_feed_ring();
		
		*forceSF_ring = 0.0;
	}
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
		Ring_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, *ring_shell_radius, 
					  *sector_angle_ring_tr, *forceSF_ring, kfactor);
		
		*forceV_ring = 0.0;
	}
	
	else
		Compute_force_feed_ring_scalar(T);
}


//
// MHD
//
void IncVF::Compute_force_feed_ring(IncVF& W)
{
	
	Compute_force_feed_ring();

	
	Ring_mult_all(basis_type, alias_switch, N, *W.Force1, *W.Force2, *W.Force3, 
								*W.V1, *W.V2, *W.V3, 
								*ring_shell_radius, *sector_angle_ring_tr, *W.forceV_ring, kfactor);

}

//
// MHD + Scalar
//

void IncVF::Compute_force_feed_ring(IncVF& W, IncSF& T)
{
	
	Compute_force_feed_ring(W);										
	// Re(f.V*), Re(fw.W*)
	

	Ring_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, *ring_shell_radius, 
						*sector_angle_ring_tr, *forceSF_ring, kfactor);		
		
}



/**********************************************************************************************

						FOR THE CYLINDRICAL RINGS
	
***********************************************************************************************/ 

void IncVF::Compute_force_feed_cylinder_ring()
{
	
	Cyl_ring_mult_all(basis_type, alias_switch, N, *Force1, *Force2, *Force3, *V1, *V2, *V3, 
						*cylinder_shell_radius, *cylinder_kpll_array_tr, 
						*forceV_cylinder_ring, kfactor);		
			
}


void IncVF::Compute_force_feed_cylinder_ring(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Compute_force_feed_cylinder_ring_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Compute_force_feed_cylinder_ring_RB(T);
}


void IncVF::Compute_force_feed_cylinder_ring_scalar(IncSF& T)
{
	
	Compute_force_feed_cylinder_ring();
	
	Cyl_ring_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, 
								*cylinder_shell_radius, *cylinder_kpll_array_tr, 
								*forceSF_cylinder_ring, kfactor);		
					
}

//
// RB Convection
//

void IncVF::Compute_force_feed_cylinder_ring_RB(IncSF& T)
{
	if (globalvar_Pr_switch == "PRZERO")
	{
		Compute_force_feed_cylinder_ring();
		
		*forceSF_cylinder_ring = 0.0;
	}
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
		Cyl_ring_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, 
						  *cylinder_shell_radius, *cylinder_kpll_array_tr, 
						  *forceSF_cylinder_ring, kfactor);
		
		*forceV_cylinder_ring = 0.0;
	}
	
	else
		Compute_force_feed_cylinder_ring_scalar(T);
		
}


//
// MHD
//
void IncVF::Compute_force_feed_cylinder_ring(IncVF& W)
{	

	Compute_force_feed_cylinder_ring();
	
	Cyl_ring_mult_all(basis_type, alias_switch, N, *W.Force1, *W.Force2, *W.Force3, 
								*W.V1, *W.V2, *W.V3, 
								*cylinder_shell_radius, *cylinder_kpll_array_tr, 
								*W.forceV_cylinder_ring, kfactor);

}

//
// MHD + Scalar
//

void IncVF::Compute_force_feed_cylinder_ring(IncVF& W, IncSF& T)
{	

	Compute_force_feed_cylinder_ring(W);						
	// Re(f.V*), Re(fw.W*)
	
	Cyl_ring_mult_all(basis_type, alias_switch, N, *T.Force, *T.F, *cylinder_shell_radius, 
								*cylinder_kpll_array_tr, *forceSF_cylinder_ring, kfactor);		
		
}



//******************************  End of compute_force_feed.cc ********************************




