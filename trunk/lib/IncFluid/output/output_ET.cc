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


/*! \file  Output_ET.cc
 * 
 * @brief  Output flux, shell-to-shell, ring-to-ring (spherical & cylinderical), and
 *			Skpq.
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */


#include "../IncFluid.h"   


//*********************************************************************************************

void IncFluid::Output_flux()
{
	
	if (my_id == master_id)
		flux_file << "%% Time = " << Tnow << endl; 	

	Compute_flux();
	Compute_force_feed_shell(); 
	
	if (helicity_flux_switch == 1) 
		if (N[2] > 1) // 3D
			Compute_kinetic_helicity_flux();
		else if (N[2]==1)	// 2D
			Compute_enstrophy_flux();
				
	if (my_id == master_id)
		for (int i=1;  i<= no_spheres; i++) 
		{
			flux_file << i << " " << (*sphere_radius)(i) << " "  
					  << (*flux_self)(i)				<< "  "
					  << (*forceV_shell)(i)				<< "     ";
			
			if (helicity_flux_switch == 1) 
				flux_file << (*flux_hk)(i);
			
			flux_file << endl;
			
		}
						
	if (my_id == master_id)		flux_file << endl;
}  

//*********************************************************************************************
// scalar
  
void IncFluid::Output_flux(IncSF& T)
{
	
	if (my_id == master_id)
		flux_file << "%% Time = " << Tnow << endl; 	

	Compute_flux(T);
	Compute_force_feed_shell(T); 
	
	if (helicity_flux_switch == 1) 
		if (N[2] > 1) // 3D
			Compute_kinetic_helicity_flux();
		else if (N[2]==1)	// 2D
			Compute_enstrophy_flux();
	
	if (my_id == master_id)
		for (int i=1;  i<= no_spheres; i++) 
		{
			flux_file << i << " " << (*sphere_radius)(i) << " "  
					<< (*flux_self)(i)					<< " "   
					<< (*flux_SF)(i)					<< "  " 
					<< (*forceV_shell)(i)				<< " "  
					<< (*forceSF_shell)(i)				<< "     ";
			
			if (helicity_flux_switch == 1) 
				flux_file << (*flux_hk)(i);
			
			flux_file << endl;
		}
		
	if (my_id == master_id)		flux_file << endl;
} 


//*********************************************************************************************
// Vector
  
void IncFluid::Output_flux(IncVF& W)
{
	
	if (my_id == master_id)
		flux_file << "%% Time = " << Tnow << endl; 	

	Compute_flux(W); 
	Compute_force_feed_shell(W);

	if (helicity_flux_switch == 1) 
		if (N[2] > 1) {
			Compute_kinetic_helicity_flux();
			Compute_magnetic_helicity_flux(W);
		}	// 3D
	
		else if (N[2]==1)	{
			Compute_enstrophy_flux();
			Compute_magnetic_enstrophy_flux(W);
		}	// 2D
	
	
	if (my_id == master_id)
		for (int i=1;  i<= no_spheres; i++) 
		{	
			flux_file << i << " " << (*sphere_radius)(i) << " "  
						<< (*flux_self)(i)				<< " "		
						<< (*flux_VF_in_out)(i)			<< " "
						<< (*flux_VF_in_in)(i)			<< "   "	
						<< (*W.flux_self)(i)			<< " " 
						<< (*W.flux_VF_in_out)(i)		<< " "		
						<< (*W.flux_VF_in_in)(i)		<< " " 
						<< (*W.flux_VF_out_out)(i)		<< "    " 
						<< (*flux_Elsasser)(i)			<< " "		
						<< (*W.flux_Elsasser)(i)		<< "    "
						<< (*forceV_shell)(i)			<< " "		
						<< (*W.forceV_shell)(i)			<< "    ";
			
			if (helicity_flux_switch == 1) 
				flux_file << (*flux_hk)(i)	<< " "	
				<< (*W.flux_hk)(i);
			
			flux_file << endl;
		}				
					
	if (my_id == master_id)		flux_file << endl;
}

//*********************************************************************************************
// Magnetoconvection
  
void IncFluid::Output_flux(IncVF& W, IncSF &T)
{
	
	if (my_id == master_id)		flux_file << "%% Time = " << Tnow << endl; 

	Compute_flux(W, T);
	Compute_force_feed_shell(W, T); 
	
	if (helicity_flux_switch == 1) 
		if (N[2] > 1) {
			Compute_kinetic_helicity_flux();
			Compute_magnetic_helicity_flux(W);
		}	// 3D
	
		else if (N[2]==1)	{
			Compute_enstrophy_flux();
			Compute_magnetic_enstrophy_flux(W);
		}	// 2D
	
	if (my_id == master_id)
		for (int i=1;  i<= no_spheres; i++) 
		{
			flux_file << i << " " << (*sphere_radius)(i) << " "  
						<< (*flux_self)(i)				<< " " 
						<< (*flux_VF_in_out)(i)			<< " "
						<< (*flux_VF_in_in)(i)			<< "   " 
						<< (*W.flux_self)(i)			<< " " 
						<< (*W.flux_VF_in_out)(i)		<< " "  
						<< (*W.flux_VF_in_in)(i)		<< " " 
						<< (*W.flux_VF_out_out)(i)		<<  " "  
						<< (*flux_Elsasser)(i)			<< " " 
						<< (*W.flux_Elsasser)(i)		<< "     "
						<< (*flux_SF)(i)				<< "      " 
						<< (*forceV_shell)(i)			<< " " 
						<< (*W.forceV_shell)(i)			<< " " 
						<< (*forceSF_shell)(i)			<< "     ";
			
			if (helicity_flux_switch == 1) 
				flux_file << (*flux_hk)(i)	<< " "	
				<< (*W.flux_hk)(i);
			
			flux_file << endl;
		}
	
	if (my_id == master_id)			flux_file << endl;
}


//*********************************************************************************************
//*********************************************************************************************
// Shell-to-shell

void IncFluid::Output_shell_to_shell()
{
	
	Compute_shell_tr();
	
	if (helicity_shell_ET_switch == 1) 
		if (N[2] > 1) // 3D
			Compute_kinetic_helicity_shell_tr();
		else if (N[2]==1)	// 2D
			Compute_enstrophy_shell_tr();
	
	if (my_id == master_id)
	{
		Range ra(1,no_shells-1);
		Range ra2(1,no_shells);

		shell_to_shell_file << "%% Time = " << Tnow << endl; 	

		
		shell_to_shell_file << "%% U to U: "	<< endl 
			<< (*shelltoshell_self)(ra,ra2)		<< endl;
		
		if (helicity_shell_ET_switch == 1) 
			shell_to_shell_file << "%% Hk to Hk: "	<< endl 
				<< (*shelltoshell_hk)(ra,ra2)		<< endl;
		
		shell_to_shell_file << endl << endl;
		
	}
}  
 
//*********************************************************************************************   
// scalar


void IncFluid::Output_shell_to_shell(IncSF& T)
{
	
	Compute_shell_tr(T);
	
	if (helicity_shell_ET_switch == 1) 
		if (N[2] > 1) // 3D
			Compute_kinetic_helicity_shell_tr();
		else if (N[2]==1)	// 2D
			Compute_enstrophy_shell_tr();
	
	
	if (my_id == master_id)
	{
		Range ra(1,no_shells-1);
		Range ra2(1,no_shells);
		
		shell_to_shell_file << "%% Time = " << Tnow << endl; 	 
		
		shell_to_shell_file << "%% U to U: "			<< endl 
				<< (*shelltoshell_self)(ra,ra2)			<< endl;
				
		shell_to_shell_file << "%% T to T: "			<< endl 
				<< (*shelltoshell_SF)(ra,ra2)			<< endl;
		
		if (helicity_shell_ET_switch == 1) 
			shell_to_shell_file << "%% Hk to Hk: "		<< endl 
				<< (*shelltoshell_hk)(ra,ra2)			<< endl;
		
		shell_to_shell_file << endl << endl;
	
	}	
} 


//*********************************************************************************************
// MHD

void IncFluid::Output_shell_to_shell(IncVF& W)
{

	DP B0mag = W.Get_mag_V0();
	
	Compute_shell_tr(W);
	
	if (B0mag > MYEPS)
		Compute_shell_ET_B0(W); 
	
	if (helicity_shell_ET_switch == 1) 
		if (N[2] > 1) {
			Compute_kinetic_helicity_shell_tr();
			Compute_magnetic_helicity_shell_tr(W);
		}
		else if (N[2]==1)	{
			Compute_enstrophy_shell_tr();
			Compute_magnetic_enstrophy_shell_tr(W);
		}

	
	if (my_id == master_id)
	{
		Range ra(1,no_shells-1);
		Range ra2(1,no_shells);
		
		shell_to_shell_file << "%% Time = "				<< Tnow << endl; 	
		
		shell_to_shell_file << "%% U to U: "			<< endl 
				<< (*shelltoshell_self)(ra,ra2)			<< endl;
				
		shell_to_shell_file << "%% W to W: "			<< endl 
				<< (*W.shelltoshell_self)(ra,ra2)		<< endl;
		
		shell_to_shell_file << "%% U to W: "			<< endl 
				<< (*shelltoshell_VF)(ra,ra2)			<< endl;
		
		shell_to_shell_file << "%% Zp to Zp: "			<< endl 
				<< (*shelltoshell_Elsasser)(ra,ra2)		<< endl;
		
		shell_to_shell_file << "%% Zm to Zm: "			<< endl 
				<< (*W.shelltoshell_Elsasser)(ra,ra2) 
				<< endl << endl << endl;
		 
		if (B0mag > MYEPS)
			shell_to_shell_file << "%% b to u due to B0 " << endl 
					<< (*energy_tr_shell_B0)(ra2)		  << endl;
		
		if (helicity_shell_ET_switch == 1) {
			shell_to_shell_file << "%% Hk to Hk: "		<< endl 
				<< (*shelltoshell_hk)(ra,ra2)			<< endl;
			
			shell_to_shell_file << "%% Hm to Hm: "		<< endl 
				<< (*W.shelltoshell_hk)(ra,ra2)			<< endl;
		}
		
		shell_to_shell_file << endl << endl;
					
	}
				
} 


//*********************************************************************************************
//
//  MHD+Scalar

void IncFluid::Output_shell_to_shell(IncVF& W, IncSF& T)
{

	DP B0mag = W.Get_mag_V0();
	
	Compute_shell_tr(W, T); 
	
	if (B0mag > MYEPS)
		Compute_shell_ET_B0(W);
	
	if (helicity_shell_ET_switch == 1) 
		if (N[2] > 1) {
			Compute_kinetic_helicity_shell_tr();
			Compute_magnetic_helicity_shell_tr(W);
		}
		else if (N[2]==1)	{
			Compute_enstrophy_shell_tr();
			Compute_magnetic_enstrophy_shell_tr(W);
		}
	
	if (my_id == master_id)
	{
		Range ra(1,no_shells-1);
		Range ra2(1,no_shells);
			
		shell_to_shell_file << "%% Time = " << Tnow << endl; 	
		
		shell_to_shell_file << "%% U to U: "			<< endl 
				<< (*shelltoshell_self)(ra,ra2)			<< endl;
				
		shell_to_shell_file << "%% W to W: "			<< endl 
				<< (*W.shelltoshell_self)(ra,ra2)		<< endl;
		
		shell_to_shell_file << "%% U to W: "			<< endl 
				<< (*shelltoshell_VF)(ra,ra2)			<< endl;
		
		shell_to_shell_file << "%% Zp to Zp: "			<< endl 
				<< (*shelltoshell_Elsasser)(ra,ra2)		<< endl;
		
		shell_to_shell_file << "%% Zm to Zm: "			<< endl 
				<< (*W.shelltoshell_Elsasser)(ra,ra2)	<< endl;
			
		shell_to_shell_file << "%% T to T: "			<< endl 
				<< (*shelltoshell_SF)(ra,ra2)			<< endl << endl;
				
		if (B0mag  > MYEPS)
			shell_to_shell_file << "%% b to u due to B0 "<< endl 
					<< (*energy_tr_shell_B0)(ra)		<< endl;
		
		if (helicity_shell_ET_switch == 1) {
			shell_to_shell_file << "%% Hk to Hk: "		<< endl 
				<< (*shelltoshell_hk)(ra,ra2)			<< endl;
			
			shell_to_shell_file << "%% Hm to Hm: "		<< endl 
				<< (*W.shelltoshell_hk)(ra,ra2)			<< endl;
		}
		
		shell_to_shell_file << endl << endl;
			
	}				
} 


//********************************************************************************************* 

// ring-to-ring (spherical)

void IncFluid::Output_ring_to_ring()
{
	if (ET_anisotropic_ring_switch == 1)
	{
		Compute_ring_tr(); 
		
		Compute_force_feed_ring();
		
		
		if (my_id == master_id)
		{
			Range ra1(1, no_ring_shells-1);
			Range ra2(1, no_sectors_ring_tr);
			
			ring_to_ring_file << "%% Time = " << Tnow << endl; 	
			
			ring_to_ring_file << "%% U to U: "						<< endl 
				<< (*ring_to_ring_self)(ra1,ra2,ra1,ra2)			<< endl << endl;
			
			ring_to_ring_file << "%% Force feed: "					<< endl
				<< (*forceV_ring)(ra1,ra2)							<< endl << endl;
				
		}		
	}	
}  

//*********************************************************************************************
// scalar


void IncFluid::Output_ring_to_ring(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_ring_to_ring_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_ring_to_ring_RB(T);
}	

 
void IncFluid::Output_ring_to_ring_scalar(IncSF& T)
{
	if (ET_anisotropic_ring_switch == 1)
	{
		Compute_ring_tr(T);
		
		Compute_force_feed_ring(T);
		
				
		if (my_id == master_id)
		{
			Range ra1(1, no_ring_shells-1);
			Range ra2(1, no_sectors_ring_tr);
					
			ring_to_ring_file << "%% Time = " << Tnow << endl; 	
			
			ring_to_ring_file << "%% U to U: "						<< endl 
				<< (*ring_to_ring_self)(ra1,ra2,ra1,ra2)			<< endl;
			
			ring_to_ring_file << "%% T to T: "						<< endl 
				<< (*ring_to_ring_SF)(ra1,ra2,ra1,ra2)				<< endl << endl;
				
			ring_to_ring_file << "%% Force feed (Fv * V): "			<< endl
				<< (*forceV_ring)(ra1,ra2)							<< endl << endl;
				
			ring_to_ring_file << "%% Force feed (Ftheta * theta): " << endl
				<< (*forceSF_ring)(ra1,ra2)							<< endl << endl;
					
		}	
	}	
}

//	RB Convection	//

void IncFluid::Output_ring_to_ring_RB(IncSF& T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_ring_to_ring();
	
	else
		Output_ring_to_ring_scalar(T);
}

//*********************************************************************************************
void IncFluid::Output_ring_to_ring(IncVF& W)
{

	if (ET_anisotropic_ring_switch == 1)
	{
		Compute_ring_tr(W);
		
		Compute_force_feed_ring(W);
		
		if (W.Get_mag_V0() > MYEPS)
			Compute_ring_ET_B0(W);
			
		
		if (my_id == master_id)
		{
			Range ra1(1, no_ring_shells-1);
			Range ra2(1, no_sectors_ring_tr);
					
			ring_to_ring_file << "%% Time = " << Tnow << endl; 	 
			
			ring_to_ring_file << "%% U to U: "						<< endl 
				<< (*ring_to_ring_self)(ra1,ra2,ra1,ra2)			<< endl;
				
			ring_to_ring_file << "%% W to W: "						<< endl 
				<< (*W.ring_to_ring_self)(ra1,ra2,ra1,ra2)			<< endl;
			
			ring_to_ring_file << "%% U to W: "						<< endl 
				<< (*ring_to_ring_VF)(ra1,ra2,ra1,ra2)				<< endl;
			
			ring_to_ring_file << "%% Zp to Zp: "					<< endl 
				<< (*ring_to_ring_Elsasser)(ra1,ra2,ra1,ra2)		<< endl;
			
			ring_to_ring_file << "%% Zm to Zm: "					<< endl 
				<< (*W.ring_to_ring_Elsasser)(ra1,ra2,ra1,ra2)		<< endl << endl;
			
					
			if (W.Get_mag_V0() > MYEPS)
				ring_to_ring_file << "%% b to u due to B0 "			<< endl 
						<< (*energy_tr_ring_B0)(Range(1,toEnd))		<< endl;

			ring_to_ring_file << "%% Force feed (Fv * V): "				<< endl
				<< (*forceV_ring)(ra1,ra2)								<< endl << endl;
				
			ring_to_ring_file << "%% Force feed (Fw * W): "				<< endl
				<< (*W.forceV_ring)(ra1,ra2)							<< endl << endl;
				
			
		}		
	}		
} 


//*********************************************************************************************
//

void IncFluid::Output_ring_to_ring(IncVF& W, IncSF& T)
{

	if (ET_anisotropic_ring_switch == 1)
	{
		Compute_ring_tr(W, T);
		
		Compute_force_feed_ring(W, T);
		
		if (W.Get_mag_V0() > MYEPS)
			Compute_ring_ET_B0(W);
			
			
		if (my_id == master_id)
		{
			Range ra1(1, no_ring_shells-1);
			Range ra2(1, no_sectors_ring_tr);
					
			ring_to_ring_file << "%% Time = " << Tnow << endl; 	
			
			ring_to_ring_file << "%% U to U: "							<< endl 
				<< (*ring_to_ring_self)(ra1,ra2, ra1,ra2)				<< endl;
			
			ring_to_ring_file << "%% W to W: "							<< endl 
				<< (*W.ring_to_ring_self)(ra1,ra2, ra1,ra2)				<< endl;
			
			ring_to_ring_file << "%% U to W: "							<< endl 
				<< (*ring_to_ring_VF)(ra1,ra2, ra1,ra2)					<< endl;
			
			ring_to_ring_file << "%% Zp to Zp: "						<< endl 
				<< (*ring_to_ring_Elsasser)(ra1,ra2, ra1,ra2)			<< endl;
			
			ring_to_ring_file << "%% Zm to Zm: "						<< endl 
				<< (*W.ring_to_ring_Elsasser)(ra1,ra2, ra1,ra2)			<< endl;
				
			ring_to_ring_file << "%% T to T: "							<< endl 
				<< (*ring_to_ring_SF)(ra1,ra2, ra1,ra2)					<< endl << endl;
				
			if (W.Get_mag_V0() > MYEPS)
				ring_to_ring_file << "%% b to u due to B0 "				<< endl 
						<< (*energy_tr_ring_B0)(Range(1,toEnd))			<< endl;
			
			ring_to_ring_file << "%% Force feed (Fv * V): "				<< endl
				<< (*forceV_ring)(ra1,ra2)								<< endl;
				
			ring_to_ring_file << "%% Force feed (Fw * W): "				<< endl
				<< (*W.forceV_ring)(ra1,ra2)							<< endl;
				
			ring_to_ring_file << "%% Force feed (Ftheta * theta): "		<< endl
				<< (*forceSF_ring)(ra1,ra2)								<< endl;	
			
		}	
	}										
} 


//*********************************************************************************************
//*********************************************************************************************
// ring-to-ring (cylinderical)

void IncFluid::Output_cylinder_ring_to_ring()
{

	if (ET_anisotropic_cylinder_switch == 1)
	{
	
		Compute_cylinder_ring_tr(); 
		
		Compute_force_feed_cylinder_ring();
		
		if (my_id == master_id)
		{
			Range ra1(1, no_cylinder_shells-1);
			Range ra2(1, no_cylinder_kpll_slabs);
			
			cylinder_ring_to_ring_file << "%% Time = " << Tnow << endl; 	
			
			cylinder_ring_to_ring_file << "%% U to U: "							<< endl 
				<< (*cylinder_ring_to_ring_self)(ra1,ra2, ra1,ra2)				<< endl << endl;
				
			cylinder_ring_to_ring_file << "%% Force feed: "						<< endl
				<< (*forceV_cylinder_ring)(ra1,ra2)								<< endl;
			
		}	
	}	
	
}

//*********************************************************************************************
// scalar

void IncFluid::Output_cylinder_ring_to_ring(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_cylinder_ring_to_ring_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_cylinder_ring_to_ring_RB(T);
}


 
void IncFluid::Output_cylinder_ring_to_ring_scalar(IncSF& T)
{
	if (ET_anisotropic_cylinder_switch == 1)
	{
		Compute_cylinder_ring_tr(T); 
		
		Compute_force_feed_cylinder_ring(T);
		
		
		if (my_id == master_id)
		{
			Range ra1(1, no_cylinder_shells-1);
			Range ra2(1, no_cylinder_kpll_slabs);
			
			cylinder_ring_to_ring_file << "%% Time = " << Tnow << endl; 	

			cylinder_ring_to_ring_file << "%% U to U: "							<< endl 
				<< (*cylinder_ring_to_ring_self)(ra1,ra2, ra1,ra2)				<< endl;
			
			cylinder_ring_to_ring_file << "%% T to T: "							<< endl 
				<< (*cylinder_ring_to_ring_SF)(ra1,ra2, ra1,ra2)				<< endl << endl;
			
			cylinder_ring_to_ring_file << "%% Force feed (Fv * V): "			<< endl
				<< (*forceV_cylinder_ring)(ra1,ra2)								<< endl;
				
			cylinder_ring_to_ring_file << "%% Force feed (Ftheta * theta): "	<< endl
				<< (*forceSF_cylinder_ring)(ra1,ra2)							<< endl << endl;
			
		}	
	}	
} 

//	RB Convection	//

void IncFluid::Output_cylinder_ring_to_ring_RB(IncSF& T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_cylinder_ring_to_ring();
	
	else
		Output_cylinder_ring_to_ring_scalar(T);
}

//*********************************************************************************************
void IncFluid::Output_cylinder_ring_to_ring(IncVF& W)
{
	if (ET_anisotropic_cylinder_switch == 1)
	{
		Compute_cylinder_ring_tr(W); 
		
		Compute_force_feed_cylinder_ring(W);
		
		if (W.Get_mag_V0() > MYEPS)
			Compute_cylinder_ring_ET_B0(W);
		
		
		if (my_id == master_id)
		{
			Range ra1(1, no_cylinder_shells-1);
			Range ra2(1, no_cylinder_kpll_slabs);
			

			cylinder_ring_to_ring_file << "%% Time = " << Tnow << endl; 	
			
			cylinder_ring_to_ring_file << "%% U to U: "							<< endl 
				<< (*cylinder_ring_to_ring_self)(ra1,ra2, ra1,ra2)				<< endl;
				
			cylinder_ring_to_ring_file << "%% W to W: "							<< endl 
				<< (*W.cylinder_ring_to_ring_self)(ra1,ra2, ra1,ra2)			<< endl;
			
			cylinder_ring_to_ring_file << "%% U to W: "							<< endl 
				<< (*cylinder_ring_to_ring_VF)(ra1,ra2, ra1,ra2)				<< endl;
			
			cylinder_ring_to_ring_file << "%% Zp to Zp: "						<< endl 
				<< (*cylinder_ring_to_ring_Elsasser)(ra1,ra2, ra1,ra2)			<< endl;
			
			cylinder_ring_to_ring_file << "%% Zm to Zm: "						<< endl 
				<< (*W.cylinder_ring_to_ring_Elsasser)(ra1,ra2, ra1,ra2)		<< endl << endl;
			 
			
			if (W.Get_mag_V0() > MYEPS)
				cylinder_ring_to_ring_file << "%% b to u due to B0 "			<< endl 
						<< (*energy_tr_cylinder_ring_B0)(Range(1,toEnd))			<< endl;
			
			cylinder_ring_to_ring_file << "%% Force feed (Fv * V): "			<< endl
				<< (*forceV_cylinder_ring)(ra1,ra2)								<< endl;
				
			cylinder_ring_to_ring_file << "%% Force feed (Fw * W): "			<< endl
				<< (*W.forceV_cylinder_ring)(ra1,ra2)							<< endl;
				
				
		}	
	}			
} 

//*********************************************************************************************
//

void IncFluid::Output_cylinder_ring_to_ring(IncVF& W, IncSF& T)
{

	if (ET_anisotropic_cylinder_switch == 1)
	{
	
		Compute_cylinder_ring_tr(W, T); 
		
		Compute_force_feed_cylinder_ring(W, T);
		
		if (W.Get_mag_V0() > MYEPS)
			Compute_cylinder_ring_ET_B0(W); 
			
		
		if (my_id == master_id)
		{
			Range ra1(1, no_cylinder_shells-1);
			Range ra2(1, no_cylinder_kpll_slabs);
					
			cylinder_ring_to_ring_file << "%% Time = " << Tnow << endl; 	

			cylinder_ring_to_ring_file << "%% U to U: "							<< endl 
				<< (*cylinder_ring_to_ring_self)(ra1,ra2, ra1,ra2)				<< endl;
			
			cylinder_ring_to_ring_file << "%% W to W: "							<< endl 
				<< (*W.cylinder_ring_to_ring_self)(ra1,ra2, ra1,ra2)			<< endl;
			
			cylinder_ring_to_ring_file << "%% U to W: "							<< endl 
				<< (*cylinder_ring_to_ring_VF)(ra1,ra2, ra1,ra2)				<< endl;
			
			cylinder_ring_to_ring_file << "%% Zp to Zp: "						<< endl 
				<< (*cylinder_ring_to_ring_Elsasser)(ra1,ra2, ra1,ra2)			<< endl;
			
			cylinder_ring_to_ring_file << "%% Zm to Zm: "						<< endl 
				<< (*W.cylinder_ring_to_ring_Elsasser)(ra1,ra2, ra1,ra2)		<< endl;
				
			cylinder_ring_to_ring_file << "%% T to T: "							<< endl 
				<< (*cylinder_ring_to_ring_SF)(ra1,ra2, ra1,ra2)				<< endl << endl;
			
			if (W.Get_mag_V0() > MYEPS)
				cylinder_ring_to_ring_file << "%% b to u due to B0 "			<< endl 
						<< (*energy_tr_cylinder_ring_B0)(Range(1,toEnd))			<< endl;
			
			cylinder_ring_to_ring_file << "%% Force feed (Fv * V): "			<< endl
				<< (*forceV_cylinder_ring)(ra1,ra2)								<< endl;
				
			cylinder_ring_to_ring_file << "%% Force feed (Fw * W): "			<< endl
				<< (*W.forceV_cylinder_ring)(ra1,ra2)							<< endl;
				
			cylinder_ring_to_ring_file << "%% Force feed (Ftheta * theta): "	<< endl
				<< (*forceSF_cylinder_ring)(ra1,ra2)							<< endl;
			
		}	
	}						
} 


//*********************************************************************************************

// Skpq
//
void IncFluid::Output_Skpq()
{
	if (ET_Skpq_switch == 1)
	{
		Compute_Skpq();

		if (my_id == master_id)
			for(int i=1; i<=no_of_Skpq_triads; i++)
				Skpq_file	<< (*triad_array)(i,0) << " " 
							<< (*triad_array)(i,1) << " " 
							<< (*triad_array)(i,2) << "     " 
							<< (*triad_array)(i,3) << " " 
							<< (*triad_array)(i,4) << " "
							<< (*triad_array)(i,5) << "    "
							<< (*Sself_array)(i) 
							<< endl;
			
		if (my_id == master_id)		Skpq_file << endl << endl;
	}
}

//*********************************************************************************************
// Scalar
//

void IncFluid::Output_Skpq(IncSF& T)
{
	if (ET_Skpq_switch == 1)
	{
		Compute_Skpq();

		
		if (my_id == master_id)
			for(int i=1; i<=no_of_Skpq_triads; i++)
				Skpq_file	<< (*triad_array)(i,0)	<< " " 
							<< (*triad_array)(i,1)	<<  " " 
							<< (*triad_array)(i,2)	<< " " 
							<< (*triad_array)(i,3)	<< " " 
							<< (*triad_array)(i,4)	<< " "
							<< (*triad_array)(i,5)	<< " "
							<< (*Sself_array)(i)	<< " " 
							<< (*S_SF_array)(i)
							<< endl; 
			
		if (my_id == master_id)		Skpq_file << endl << endl;	
	}
}


//*********************************************************************************************
// Vector
//

void IncFluid::Output_Skpq(IncVF& W)
{
	if (ET_Skpq_switch == 1)
	{
		Compute_Skpq();

		if (my_id == master_id)
			for(int i=1; i<=no_of_Skpq_triads; i++)
				Skpq_file	<< (*triad_array)(i,0) << " " 
							<< (*triad_array)(i,1) <<  " " 
							<< (*triad_array)(i,2) << " " 
							<< (*triad_array)(i,3) << " " 
							<< (*triad_array)(i,4) << " "
							<< (*triad_array)(i,5) << "   "
							<< (*Sself_array)(i)   << " " 
							<< (*W.Sself_array)(i) << " " 
							<< (*S_to_VF_array)(i)
							<< endl; 
		
		if (my_id == master_id)		Skpq_file << endl << endl;
	}
}


//*********************************************************************************************
// Vector+Scalar
//

void IncFluid::Output_Skpq(IncVF& W, IncSF& T)
{
	if (ET_Skpq_switch == 1)
	{
		Compute_Skpq();

		if (my_id == master_id)
			for(int i=1; i<=no_of_Skpq_triads; i++)
				Skpq_file	<< (*triad_array)(i,0) << " " 
							<< (*triad_array)(i,1) <<  " " 
							<< (*triad_array)(i,2) << " " 
							<< (*triad_array)(i,3) << " " 
							<< (*triad_array)(i,4) << " "
							<< (*triad_array)(i,5) << "    "
							<< (*Sself_array)(i)   << " " 
							<< (*W.Sself_array)(i) << " " 
							<< (*S_to_VF_array)(i) << " " 
							<< (*S_SF_array)(i) 
							<< endl; 
		
		
		if (my_id == master_id)		Skpq_file << endl << endl;
	}

}



//*******************************  End of Output_ET.cc  ***************************************




