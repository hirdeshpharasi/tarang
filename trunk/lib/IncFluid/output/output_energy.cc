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

/*! \file  Output_energy.cc
 * 
 * @brief  Output global, spectrum (shell, rings, cylinderical rings). 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"   


//*********************************************************************************************

void IncFluid::Output_global()
{

	CV_Compute_totalenergy_diss();
	CV_Compute_entropy(); 
	CV_Compute_total_helicity();
	
	if (my_id == master_id) 
	{
		DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
		
		DP total_diss = dissipation_coefficient*CV_total_dissipation;
		
		DP kolm_scale_u = sqrt(sqrt((pow3(dissipation_coefficient) / total_diss)));
		
		DP kmax_eta = kmax*kolm_scale_u;
		
		DP Rlambda = 2*CV_total_energy* sqrt(15/(dissipation_coefficient* total_diss)); 

		
		global_file  << Tnow				<< "	"  
			<<	CV_total_energy				<< " " 
			<<	CV_total_dissipation		<< " " 
			<<	total_diss					<< "    "
			<<	CV_total_helicity1			<< " " 
			<<  CV_total_helicity2			<< "    " 
			<<	CV_total_dissipation_H1		<< " " 
			<<  CV_total_dissipation_H2		<< "    "
			<<	CV_entropy					<< "    "
			<<	kmax_eta					<< "    "
			<<  dissipation_coefficient		<< "    "
			<<	Rlambda						<< "    "
			<<	Tdt							<< "    "
			<<	real((*V1)(0,0,0))			<< " "	
			<<  real((*V2)(0,0,0))			<< " " 
			<<  real((*V3)(0,0,0))
			<<	endl;
	}		
}

//*********************************************************************************************


void IncFluid::Output_global(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_global_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_global_RB(T);
}


void IncFluid::Output_global_scalar(IncSF& T)
{
	CV_Compute_totalenergy_diss(); 
	T.CS_Compute_totalenergy_diss();
	
	CV_Compute_total_helicity();
	
	CV_Compute_entropy(); 
	T.CS_Compute_entropy(); 
	
	if (my_id == master_id) 
	{
		DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
		
		DP total_diss = dissipation_coefficient*CV_total_dissipation;
		
		DP kolm_scale_u = sqrt(sqrt((pow3(dissipation_coefficient) / total_diss)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = dissipation_coefficient/T.diffusion_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
	
		
		DP Rlambda = 2*CV_total_energy* sqrt(15/(dissipation_coefficient* total_diss));  
			
		
		global_file  << Tnow										<< "	"  
			<<	CV_total_energy										<< " " 
			<<	T.CS_total_energy									<< "	" 
			<<	CV_total_dissipation								<< " " 
			<<	T.CS_total_dissipation								<< "	" 
			<<	total_diss											<< " "
			<<	(T.diffusion_coefficient)*(T.CS_total_dissipation)	<< "    "
			<<	CV_total_helicity1									<< " " 
			<<  CV_total_helicity2									<< "    " 
			<<	CV_total_dissipation_H1								<< " " 
			<<  CV_total_dissipation_H2								<< "    "
			<<	CV_entropy											<< " "
			<<	T.CS_entropy										<< "   " 
			<<	kmax_eta1											<< " " 
			<<	kmax_eta2											<< "    "
			<<  dissipation_coefficient								<< " " 
			<<  T.diffusion_coefficient								<< "    "
			<<	Rlambda												<< "    "
			<<	Tdt													<< "    "
			<<	real((*V1)(0,0,0))									<< " "	
			<<  real((*V2)(0,0,0))									<< " " 
			<<  real((*V3)(0,0,0))									<< "	"
			<<	real((*T.F)(0,0,0))
			<<	endl;
	}	
	
}

//*********************************************************************************************

// RB Convection //

void IncFluid::Output_global_RB(IncSF& T)		
{
	
	static DP nusselt_no;
	
	CV_Compute_totalenergy_diss(); 
	T.CS_Compute_totalenergy_diss();
	
	CV_Compute_entropy(); 
	T.CS_Compute_entropy(); 
	
	CV_Compute_total_helicity();
	
	nusselt_no = Get_Nusselt_no(T);
	
	if (my_id == master_id) 
	{
		DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
		
		DP total_diss = dissipation_coefficient*CV_total_dissipation;
		
		DP kolm_scale_u = sqrt(sqrt((pow3(dissipation_coefficient) / total_diss)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = dissipation_coefficient/T.diffusion_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		DP Rlambda = 2*CV_total_energy* sqrt(15/(dissipation_coefficient* total_diss));  
		
		
		global_file  << Tnow										<< "	"  
			<<	CV_total_energy										<< " " 
			<<	T.CS_total_energy									<< "	" 
			<<	CV_total_dissipation								<< " " 
			<<	T.CS_total_dissipation								<< "	" 
			<<	total_diss											<< " "
			<<	(T.diffusion_coefficient)*(T.CS_total_dissipation)	<< "    "
			<<	nusselt_no											<< "    "
			<<	CV_total_helicity1									<< " " 
			<<  CV_total_helicity2									<< "    " 
			<<	CV_total_dissipation_H1								<< " " 
			<<  CV_total_dissipation_H2								<< "    "
			<<	CV_entropy											<< " "
			<<	T.CS_entropy										<< "   " 
			<<	kmax_eta1											<< " " 
			<<	kmax_eta2											<< "    "
			<<  dissipation_coefficient								<< " " 
			<<  T.diffusion_coefficient								<< "    "
			<<	Rlambda												<< "    "
			<<	Tdt													<< "    "
			<<	real((*V1)(0,0,0))									<< " "	
			<<  real((*V2)(0,0,0))									<< " " 
			<<  real((*V3)(0,0,0))									<< "	"
			<<	real((*T.F)(0,0,0))
			<<	endl;	
	}
}

//*********************************************************************************************

void IncFluid::Output_global(IncVF& W)
{
	static DP Hc;
	
		
	CV_Compute_totalenergy_diss(); 
	W.CV_Compute_totalenergy_diss();
	
		
	Hc = Get_cross_helicity(W);
	
	CV_Compute_total_helicity();
	W.CV_Compute_total_helicity();
	
	
	CV_Compute_entropy();
	W.CV_Compute_entropy();
	
		
	if (my_id == master_id) 
	{
		DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
		
		DP total_diss = dissipation_coefficient*CV_total_dissipation;
		
		DP kolm_scale_u = sqrt(sqrt((pow3(dissipation_coefficient) / total_diss)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = dissipation_coefficient/W.dissipation_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		DP Rlambda = 2*CV_total_energy* sqrt(15/(dissipation_coefficient* total_diss)); 
		
				
		global_file  << Tnow											<< "	"  
			<<	CV_total_energy											<< " " 
			<<	W.CV_total_energy										<< "   " 
			<<	CV_total_dissipation									<< " " 
			<<	W.CV_total_dissipation									<< " " 
			<<	total_diss												<< " "
			<<	(W.dissipation_coefficient)*(W.CV_total_dissipation)	<< "    " 
			<<	Hc														<< "    "
			<<	CV_total_helicity1										<< " " 
			<<  CV_total_helicity2										<< " " 
			<<	W.CV_total_helicity1									<< " " 
			<<  W.CV_total_helicity2									<< "    "
			<<	CV_total_dissipation_H1									<< " " 
			<<  CV_total_dissipation_H2									<< "    "
			<<	W.CV_total_dissipation_H1								<< " " 
			<<  W.CV_total_dissipation_H2								<< "    "
			<<	CV_entropy												<< " " 
			<<	W.CV_entropy											<< "   " 
			<<	kmax_eta1												<< " " 
			<<	kmax_eta2												<< "    "
			<<  dissipation_coefficient									<< " " 
			<<  W.dissipation_coefficient								<< "    "
			<<	Rlambda													<< "    "
			<<	Tdt														<< "    "
			<<	real((*V1)(0,0,0))										<< " "	
			<<  real((*V2)(0,0,0))										<< " " 
			<<  real((*V3)(0,0,0))										<< "	"
			<<	real((*W.V1)(0,0,0))									<< " "	
			<<  real((*W.V2)(0,0,0))									<< " " 
			<<  real((*W.V3)(0,0,0)) 
			<<	endl;
	}	
}

//*********************************************************************************************

void IncFluid::Output_global(IncVF& W, IncSF& T)		
{	
	static DP Hc;

	CV_Compute_totalenergy_diss(); 
	W.CV_Compute_totalenergy_diss();
	T.CS_Compute_totalenergy_diss();
	
	Hc = Get_cross_helicity(W);
	
	CV_Compute_entropy();
	W.CV_Compute_entropy();	
	T.CS_Compute_entropy();
	
	CV_Compute_total_helicity();
	W.CV_Compute_total_helicity();
	
	if (my_id == master_id) 
	{
		DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
		
		DP total_diss = dissipation_coefficient*CV_total_dissipation;
		
		DP kolm_scale_u = sqrt(sqrt((pow3(dissipation_coefficient) / total_diss)));
		
		DP kmax_eta1 = kmax * kolm_scale_u;
		
		DP tempvar = dissipation_coefficient/W.dissipation_coefficient;
		DP kmax_eta2 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		tempvar = dissipation_coefficient/T.diffusion_coefficient;
		DP kmax_eta3 = kmax_eta1 * sqrt(sqrt(tempvar))/tempvar;
		
		DP Rlambda = 2*CV_total_energy* sqrt(15/(dissipation_coefficient* total_diss)); 
		 
		global_file  << Tnow											<< "	"  
			<<	CV_total_energy											<< " " 
			<<	W.CV_total_energy										<< "  "
			<<	T.CS_total_energy										<< "    " 
			<<	CV_total_dissipation									<< " " 
			<<	W.CV_total_dissipation									<< " " 
			<<	T.CS_total_dissipation									<< "	"
			<<	total_diss												<< " " 
			<<	(W.dissipation_coefficient)*(W.CV_total_dissipation)	<< " "
			<<	(T.diffusion_coefficient)*(T.CS_total_dissipation)		<< "    " 
			<<	Hc														<< "    "
			<<	CV_total_helicity1										<< " " 
			<<  CV_total_helicity2										<< " " 
			<<	W.CV_total_helicity1									<< " " 
			<<  W.CV_total_helicity2									<< "     "
			<<	CV_total_dissipation_H1									<< " " 
			<<  CV_total_dissipation_H2									<< "    "
			<<	W.CV_total_dissipation_H1								<< " " 
			<<  W.CV_total_dissipation_H2								<< "    "
			<<	CV_entropy												<< " " 
			<<	W.CV_entropy											<<  " " 
			<<	T.CS_entropy											<< "    " 
			<<	kmax_eta1												<< " " 
			<<	kmax_eta2												<< " " 
			<<	kmax_eta3												<< "    "
			<<  dissipation_coefficient									<< " " 
			<<  W.dissipation_coefficient								<< " " 
			<<  T.diffusion_coefficient									<< "    "
			<<	Rlambda													<< "    "
			<<	Tdt														<< "    "
			<<	real((*V1)(0,0,0))										<< " "	
			<<  real((*V2)(0,0,0))										<< " " 
			<<  real((*V3)(0,0,0))										<< "	"
			<<	real((*W.V1)(0,0,0))									<< " "	
			<<  real((*W.V2)(0,0,0))									<< " " 
			<<  real((*W.V3)(0,0,0))									<< "    "
			<<	real((*T.F)(0,0,0))
			<<	endl;
	}
		 
}

//*********************************************************************************************
//*********************************************************************************************
 
void IncFluid::Output_shell_spectrum()
{
	
	if (my_id == master_id) 
		spectrum_file << "%% Time = " << Tnow << endl; 	

	CV_Compute_shell_spectrum(); 
	
	Compute_force_shell_spectrum();

	CV_Compute_shell_spectrum_helicity();
	
	if (my_id == master_id) 
		for (int i=0; i< CV_shell_ek_size; i++)	
			spectrum_file << i										<< "	" 
				<< (*CV_shell_ek1)(i)								<< " " 
				<< (*CV_shell_ek2)(i)								<< " " 
				<< (*CV_shell_ek3)(i)								<< "	" 
				<< dissipation_coefficient*((*CV_shell_dissk)(i))	<< "	" 
				<< (*CV_shell_h1k1)(i)								<< "  "
				<< (*CV_shell_h1k2)(i)								<< "  "
				<< (*CV_shell_h1k3)(i)								<< "	"
				<< (*shell_spectrum_force_Vk)(i)					<< "	" 
				<< endl;	

	if (my_id == master_id)		spectrum_file << endl;
}  

//*********************************************************************************************  
// scalar
void IncFluid::Output_shell_spectrum(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_shell_spectrum_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_shell_spectrum_RB(T);
}


  
void IncFluid::Output_shell_spectrum_scalar(IncSF& T)
{
	
	if (my_id == master_id) 
		spectrum_file << "%% Time = " << Tnow << endl; 	

	CV_Compute_shell_spectrum(); 
	T.CS_Compute_shell_spectrum();
	
	Compute_cross_vT_shell_spectrum(T);
	
	Compute_force_shell_spectrum();
	T.Compute_force_shell_spectrum();

	CV_Compute_shell_spectrum_helicity();	
	
	if (my_id == master_id) 
		for (int i=0; i< CV_shell_ek_size; i++)
			spectrum_file << i										<< "	" 
				<<	(*CV_shell_ek1)(i)								<< " " 
				<<	(*CV_shell_ek2)(i)								<< " " 
				<<	(*CV_shell_ek3)(i)								<< "	"  
				<<	(*T.CS_shell_ek)(i)								<< "	" 
				<<	dissipation_coefficient*((*CV_shell_dissk)(i))	<< " " 
				<<	T.diffusion_coefficient*((*CS_shell_dissk)(i))	<< "	"			
				<<	(*shell_ek_cross_V1T)(i)						<< " " 
				<<	(*shell_ek_cross_V2T)(i)						<< " "
				<<	(*shell_ek_cross_V3T)(i)						<< "	"  
				<<	(*CV_shell_h1k1)(i)								<< "  "
				<<	(*CV_shell_h1k2)(i)								<< "  "
				<<	(*CV_shell_h1k3)(i)								<< "	"
				<<	(*shell_spectrum_force_Vk)(i)					<< " " 
				<<	(*T.shell_spectrum_force_SFk)(i)			
				<<	endl;
	
	if (my_id == master_id)		spectrum_file << endl;
}

//  RB Convection	


void IncFluid::Output_shell_spectrum_RB(IncSF &T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_shell_spectrum();
	
	else
		Output_shell_spectrum_scalar(T);
}
 
//********************************************************************************************* 
// Vector
  
void IncFluid::Output_shell_spectrum(IncVF& W)
{
	
	if (my_id == master_id) 
		spectrum_file << "%% Time = " << Tnow << endl; 	

	CV_Compute_shell_spectrum(); 
	W.CV_Compute_shell_spectrum();
	
	Compute_force_shell_spectrum();
	W.Compute_force_shell_spectrum();
		
	
	CV_Compute_shell_spectrum_helicity();	
	W.CV_Compute_shell_spectrum_helicity();
	
	if (my_id == master_id) 
		for (int i=0; i< CV_shell_ek_size; i++)
			spectrum_file << i											<< "	" 
				<<	(*CV_shell_ek1)(i)									<< " " 
				<<	(*CV_shell_ek2)(i)									<< " " 
				<<	(*CV_shell_ek3)(i)									<< "	"  
				<<	(*W.CV_shell_ek1)(i)								<< " "  
				<<	(*W.CV_shell_ek2)(i)								<< " " 
				<<	(*W.CV_shell_ek3)(i)								<< "	" 
				<<	dissipation_coefficient*((*CV_shell_dissk)(i))		<< " " 
				<<	W.dissipation_coefficient*((*W.CV_shell_dissk)(i))	<< "    "			
				<<	(*shell_ek_cross)(i)								<< "    " 
				<<	(*CV_shell_h1k1)(i)									<< " "
				<<	(*CV_shell_h1k2)(i)									<< " "
				<<	(*CV_shell_h1k3)(i)									<< "	" 
				<<	(*W.CV_shell_h1k1)(i)								<< " "
				<<	(*W.CV_shell_h1k2)(i)								<< " "
				<<	(*W.CV_shell_h1k3)(i)								<< "	"
				<<	(*shell_spectrum_force_Vk)(i)						<< " " 
				<<	(*W.shell_spectrum_force_Vk)(i)						
				<<	endl;
	
	if (my_id == master_id)		spectrum_file << endl;
}


//*********************************************************************************************
// Magneto+scalar
  
void IncFluid::Output_shell_spectrum(IncVF& W, IncSF &T)
{
	
	if (my_id == master_id) 
		spectrum_file << "%% Time = " << Tnow << endl; 	

	CV_Compute_shell_spectrum(); 
	W.CV_Compute_shell_spectrum();
	T.CS_Compute_shell_spectrum();

	Compute_cross_vT_shell_spectrum(T);
	
	Compute_force_shell_spectrum();
	W.Compute_force_shell_spectrum();
	
	T.Compute_force_shell_spectrum();

	CV_Compute_shell_spectrum_helicity();	
	W.CV_Compute_shell_spectrum_helicity();
	
	if (my_id == master_id) 
		for (int i=0; i< CV_shell_ek_size; i++)
			spectrum_file << i											<< "	" 
				<<	(*CV_shell_ek1)(i)									<< " " 
				<<	(*CV_shell_ek2)(i)									<< " " 
				<<	(*CV_shell_ek3)(i)									<< "	"  
				<<	(*W.CV_shell_ek1)(i)								<< " "  
				<<	(*W.CV_shell_ek2)(i)								<< " " 
				<<	(*W.CV_shell_ek3)(i)								<< "	" 
				<<	(*T.CS_shell_ek)(i)									<< "	"
				<<	dissipation_coefficient*((*CV_shell_dissk)(i))		<< " " 
				<<	W.dissipation_coefficient*((*W.CV_shell_dissk)(i))	<< " " 
				<<	T.diffusion_coefficient*((*T.CS_shell_dissk)(i))	<< "    "			  
				<<	(*shell_ek_cross)(i)								<< "    " 
				<<	(*shell_ek_cross_V1T)(i)							<< " " 
				<<	(*shell_ek_cross_V2T)(i)							<< " "
				<<	(*shell_ek_cross_V3T)(i)							<< "   "  
				<<	(*CV_shell_h1k1)(i)									<< " "
				<<	(*CV_shell_h1k2)(i)									<< " "
				<<	(*CV_shell_h1k3)(i)									<< "	" 
				<<	(*W.CV_shell_h1k1)(i)								<< " "
				<<	(*W.CV_shell_h1k2)(i)								<< " "
				<<	(*W.CV_shell_h1k3)(i)								<< "	" 
				<<	(*shell_spectrum_force_Vk)(i)						<< " " 
				<<	(*W.shell_spectrum_force_Vk)(i)						<< " " 
				<<	(*T.shell_spectrum_force_SFk)(i)
				<<	endl;

	if (my_id == master_id)		spectrum_file << endl;
}



 
/**********************************************************************************************
	
						Output_ring_spectrum()

***********************************************************************************************/
 
 
void IncFluid::Output_ring_spectrum()
{
	if (CV_anisotropic_ring_switch == 1)
	{
	
		if (my_id == master_id)
			ring_spectrum_file << "%% Time = " << Tnow << endl; 	

		CV_Compute_ring_spectrum();
		Compute_force_ring_spectrum();

		CV_Compute_ring_spectrum_helicity();
		
		if (my_id == master_id)
		{
			Range ra(1, toEnd);
			
			ring_spectrum_file <<  "%% ek1:"	<<   endl	<<  (*CV_ring_ek1)(ra,ra)	
								<< endl << endl;
								
			ring_spectrum_file <<  "%% ek2: "	<<   endl	<<  (*CV_ring_ek2)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% dissk: " <<   endl	<<  (*CV_ring_dissk)(ra,ra)		
								<< endl << endl;
			
			ring_spectrum_file <<  "%% Hk: "	<<   endl	<<  (*CV_ring_h1k)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Vk*Fk: " <<   endl	<<  (*ring_spectrum_force_Vk)(ra,ra) 
							   << endl << endl;
		}					   
	}
	
}  
 

//*********************************************************************************************  
// scalar
void IncFluid::Output_ring_spectrum(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_ring_spectrum_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_ring_spectrum_RB(T);
}
  
void IncFluid::Output_ring_spectrum_scalar(IncSF& T)
{

	if (CV_anisotropic_ring_switch == 1)
	{
		
		if (my_id == master_id)
			ring_spectrum_file << "%% Time = " << Tnow << endl; 	

		CV_Compute_ring_spectrum(); 
		T.CS_Compute_ring_spectrum();
		
		Compute_cross_vT_ring_spectrum(T);
			
		Compute_force_ring_spectrum();
		T.Compute_force_ring_spectrum();
			
		CV_Compute_ring_spectrum_helicity();
		
		if (my_id == master_id)
		{
			Range ra(1, toEnd);
			
			ring_spectrum_file <<  "%% ek1:"	<<   endl	<<  (*CV_ring_ek1)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% ek2:"	<<   endl	<<  (*CV_ring_ek2)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Tek:"	<<   endl	<<  (*T.CS_ring_ek)(ra,ra)		
								<< endl << endl;
			
			ring_spectrum_file <<  "%% dissk:"	<<   endl	<<	(*CV_ring_dissk)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Tdissk:"	<<   endl	<<  (*T.CS_ring_dissk)(ra,ra)	
								<< endl << endl;
			
			ring_spectrum_file <<  "%% V1k*Tk:"	<<   endl	<<  (*ring_ek_cross_V1T)(ra,ra)	
								<< endl << endl;
								
			ring_spectrum_file <<  "%% V2k*Tk:"	<<   endl	<<  (*ring_ek_cross_V2T)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% V3k*Tk:"	<<   endl	<<  (*ring_ek_cross_V3T)(ra,ra)
								<< endl << endl;
				
			ring_spectrum_file <<  "%% Hk:"		<<   endl	<<	(*CV_ring_h1k)(ra,ra)
								<< endl	<< endl;
			
			ring_spectrum_file <<  "%% Vk*Fk:"	<<   endl	<<  (*ring_spectrum_force_Vk)(ra,ra) 
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Vk*SFk:" <<   endl	
							   <<  (*T.ring_spectrum_force_SFk)(ra,ra)
							   << endl	<< endl;		
		}
	}
}

//  RB Convection	//

void IncFluid::Output_ring_spectrum_RB(IncSF &T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_ring_spectrum();
	
	else
		Output_ring_spectrum_scalar(T);
}
 
//********************************************************************************************* 
// Vector
  
void IncFluid::Output_ring_spectrum(IncVF& W)
{
	if (CV_anisotropic_ring_switch == 1)
	{
		
		if (my_id == master_id)
			ring_spectrum_file << "%% Time = " << Tnow << endl; 	

		CV_Compute_ring_spectrum(); 
		W.CV_Compute_ring_spectrum();
		
		Compute_cross_helicity_ring_spectrum(W);
		
		Compute_force_ring_spectrum();
		W.Compute_force_ring_spectrum();

		CV_Compute_ring_spectrum_helicity();
		W.CV_Compute_ring_spectrum_helicity();
		
		if (my_id == master_id)
		{
			Range ra(1, toEnd);
			
			ring_spectrum_file <<  "%% ek1:"	<<   endl	<<  (*CV_ring_ek1)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% ek2:"	<<   endl	<<  (*CV_ring_ek2)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Wek1:"	<<   endl	<<  (*W.CV_ring_ek1)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Wek2:"	<<   endl	<<  (*W.CV_ring_ek2)(ra,ra)
								<< endl << endl;
			
			ring_spectrum_file <<  "%% dissk:"	<<   endl	<<	(*CV_ring_dissk)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Wdissk:"	<<   endl	<<  (*W.CV_ring_dissk)(ra,ra)	
								<< endl << endl;
			
			ring_spectrum_file <<  "%% Vk*Wk:"	<<   endl	<<  (*ring_ek_cross)(ra,ra)
								<< endl << endl;
				
			ring_spectrum_file <<  "%% Hk:"		<<   endl	<<  (*CV_ring_h1k)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% W.Hk:"	<<   endl	<<  (*W.CV_ring_h1k)(ra,ra)		
								<< endl << endl;
			
			ring_spectrum_file <<  "%% Vk*Fk:"	<<   endl	<<  (*ring_spectrum_force_Vk)(ra,ra)
								<< endl	<< endl;
								
			ring_spectrum_file <<  "%% Wk*Fk:" <<   endl	
								<<  (*W.ring_spectrum_force_Vk)(ra,ra) 
								<< endl	<< endl;	
		}																											
	}

}


//*********************************************************************************************	
// Magneto+scalar
  
void IncFluid::Output_ring_spectrum(IncVF& W, IncSF &T)
{
	if (CV_anisotropic_ring_switch == 1)
	{
		
		if (my_id == master_id)
			ring_spectrum_file << "%% Time = " << Tnow << endl; 	

		CV_Compute_ring_spectrum(); 
		W.CV_Compute_ring_spectrum();
		T.CS_Compute_ring_spectrum();
		
		Compute_cross_helicity_ring_spectrum(W);
		Compute_cross_vT_ring_spectrum(T);
		
		Compute_force_ring_spectrum();
		W.Compute_force_ring_spectrum();
		T.Compute_force_ring_spectrum();
		
		CV_Compute_ring_spectrum_helicity();
		W.CV_Compute_ring_spectrum_helicity();
		
		if (my_id == master_id)
		{
			Range ra(1, toEnd);
			
			ring_spectrum_file <<  "%% ek1:"	<<   endl	<<  (*CV_ring_ek1)(ra,ra)		
								<< endl << endl;
								
			ring_spectrum_file <<  "%% ek2:"	<<   endl	<<  (*CV_ring_ek2)(ra,ra)
								<< endl << endl;
									
			ring_spectrum_file <<  "%% Wek1:"	<<   endl	<<  (*W.CV_ring_ek1)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Wek2:"	<<   endl	<<  (*W.CV_ring_ek2)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Tek:"	<<   endl	<<  (*T.CS_ring_ek)(ra,ra)
								<< endl << endl;
			
			ring_spectrum_file <<  "%% dissk:"	<<   endl	<<	(*CV_ring_dissk)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Wdissk:"	<<   endl	<<  (*W.CV_ring_dissk)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% Tdissk:"	<<   endl	<<  (*T.CS_ring_dissk)(ra,ra)
								<< endl << endl;
			
			ring_spectrum_file <<  "%% Vk*Wk:"	<<   endl	<<  (*ring_ek_cross)(ra,ra)
								<< endl << endl;
			
			ring_spectrum_file <<  "%% V1k*Tk:"	<<   endl	<<  (*ring_ek_cross_V1T)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% V2k*Tk:"	<<   endl	<<  (*ring_ek_cross_V2T)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% V3k*Tk:"	<<   endl	<<  (*ring_ek_cross_V3T)(ra,ra)
								<< endl << endl;
				
			ring_spectrum_file <<  "%% Hk:"		<<   endl	<<  (*CV_ring_h1k)(ra,ra)
								<< endl << endl;
								
			ring_spectrum_file <<  "%% W.Hk:"	<<   endl	<<  (*W.CV_ring_h1k)(ra,ra)
								<< endl << endl;
			
			ring_spectrum_file <<  "%% Vk*Fk:"	<<   endl	<<  (*ring_spectrum_force_Vk)(ra,ra)
								<< endl	<< endl;
								
			ring_spectrum_file <<  "%% Wk*Fk:" <<   endl	
								<<  (*W.ring_spectrum_force_Vk)(ra,ra) 
								<< endl	<< endl;	
								
			ring_spectrum_file <<  "%% Vk*SFk:" <<   endl	
								<<  (*T.ring_spectrum_force_SFk)(ra,ra)
								<< endl	<< endl;
		}																															
	}

}


/**********************************************************************************************
	
						Output_cylinder_ring_spectrum()

***********************************************************************************************/
 
 
void IncFluid::Output_cylinder_ring_spectrum()
{

	if (CV_anisotropic_cylinder_switch == 1)
	{
		CV_Compute_cylinder_ring_spectrum();	
		
		Compute_force_cylinder_ring_spectrum();
		
		CV_Compute_cylinder_ring_spectrum_helicity();
		
		
		if (my_id == master_id)	
		{
			Range ra(1, toEnd);
			
			cylinder_ring_spectrum_file << "%% Time = " << Tnow << endl; 
			
			cylinder_ring_spectrum_file << (*CV_cylinder_ring_ek1)(2,1) << " " 
						<< (*CV_cylinder_ring_ek1)(2,2) << " "	
						<< (*CV_cylinder_ring_ek1)(2,3) << " "	
						<< (*CV_cylinder_ring_ek1)(2,4) << " "
						<< (*CV_cylinder_ring_ek1)(2,5) << endl << endl;
			
			cylinder_ring_spectrum_file <<  "%% ek1:"	<<   endl	
						<<  (*CV_cylinder_ring_ek1)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% ek2: "	<<   endl	
						<<  (*CV_cylinder_ring_ek2)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% dissk: " <<   endl	
						<<  (*CV_cylinder_ring_dissk)(ra,ra)	<< endl << endl;
			
			cylinder_ring_spectrum_file <<  "%% Hk: "	<<   endl	
						<<  (*CV_cylinder_ring_h1k)(ra,ra)		<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Vk*Fk: " <<   endl	
						<<  (*cylinder_ring_spectrum_force_Vk)(ra,ra)  << endl << endl;
		}				
	}								

}  
  
//
//

//
//  scalar
void IncFluid::Output_cylinder_ring_spectrum(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_cylinder_ring_spectrum_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_cylinder_ring_spectrum_RB(T);
}


void IncFluid::Output_cylinder_ring_spectrum_scalar(IncSF& T)
{
	
	if (CV_anisotropic_cylinder_switch == 1)
	{
		CV_Compute_cylinder_ring_spectrum(); 
		T.CS_Compute_cylinder_ring_spectrum();
		
		Compute_cross_vT_cylinder_ring_spectrum(T);
			
		Compute_force_cylinder_ring_spectrum();
		T.Compute_force_cylinder_ring_spectrum();
		
		CV_Compute_cylinder_ring_spectrum_helicity();
		
		if (my_id == master_id)	
		{
			Range ra(1, toEnd);
			
			cylinder_ring_spectrum_file << "%% Time = " << Tnow << endl; 
		
			cylinder_ring_spectrum_file <<  "%% ek1:"		<<   endl	
						<<  (*CV_cylinder_ring_ek1)(ra,ra) 	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% ek2: "		<<   endl	
						<<  (*CV_cylinder_ring_ek2)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Tek:"		<<   endl	
						<<  (*T.CS_cylinder_ring_ek)(ra,ra) << endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% dissk: "	<<   endl	
						<<  (*CV_cylinder_ring_dissk)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Tdissk: "	<<   endl	
						<<  (*T.CS_cylinder_ring_dissk)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<	"%% V1k*Tk: "	<<  endl	
						<<  (*cylinder_ring_ek_cross_V1T)(ra,ra)	<< endl << endl;
																								
			cylinder_ring_spectrum_file <<	"%% V2k*Tk: "	<<  endl	
						<<  (*cylinder_ring_ek_cross_V2T)(ra,ra)	<< endl << endl;
			
			cylinder_ring_spectrum_file <<	"%% V3k*Tk: "	<<  endl			
						<<  (*cylinder_ring_ek_cross_V3T)(ra,ra)	<< endl << endl;							
										
			cylinder_ring_spectrum_file <<  "%% Hk: "	<<   endl	
						<<  (*CV_cylinder_ring_h1k)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Vk*Fk: " <<   endl	
						<<  (*cylinder_ring_spectrum_force_Vk)(ra,ra) << endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Vk*SFk: " <<   endl	
						<<  (*T.cylinder_ring_spectrum_force_SFk)(ra,ra) << endl << endl;
		}				
	}															
								
} 

//  RB Convection	//

void IncFluid::Output_cylinder_ring_spectrum_RB(IncSF &T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_cylinder_ring_spectrum();
	
	else
		Output_cylinder_ring_spectrum_scalar(T);
}

//
//

void IncFluid::Output_cylinder_ring_spectrum(IncVF& W)
{

	if (CV_anisotropic_cylinder_switch == 1)
	{

		CV_Compute_cylinder_ring_spectrum();	
		W.CV_Compute_cylinder_ring_spectrum();
		
		Compute_cross_helicity_cylinder_ring_spectrum(W);

			
		CV_Compute_cylinder_ring_spectrum_helicity();	
		W.CV_Compute_cylinder_ring_spectrum_helicity();
		
		Compute_force_cylinder_ring_spectrum();
		W.Compute_force_cylinder_ring_spectrum();

		if (my_id == master_id)	
		{
			Range ra(1, toEnd);
			
			cylinder_ring_spectrum_file << "%% Time = " << Tnow << endl;  	
		
			cylinder_ring_spectrum_file <<  "%% ek1:"		<<   endl		
						<<  (*CV_cylinder_ring_ek1)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% ek2: "		<<   endl	
						<<  (*CV_cylinder_ring_ek2)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Wek1:"		<<   endl	
						<<  (*W.CV_cylinder_ring_ek1)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Wek2: "		<<   endl	
						<<  (*W.CV_cylinder_ring_ek2)(ra,ra)	<< endl << endl;
																	
																		
			cylinder_ring_spectrum_file <<  "%% dissk: "	<<   endl	
						<<  (*CV_cylinder_ring_dissk)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Wdissk: "	<<   endl	
						<<  (*W.CV_cylinder_ring_dissk)(ra,ra)	<< endl << endl;
																		
										
			cylinder_ring_spectrum_file <<  "%% crosshel:"	<<   endl	
						<<  (*cylinder_ring_ek_cross)(ra,ra)	<< endl << endl;																						
			
			cylinder_ring_spectrum_file <<  "%% Hk: "		<<   endl	
						<<  (*CV_cylinder_ring_h1k)(ra,ra)		<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% W.Hk: "		<<   endl	
						<<  (*W.CV_cylinder_ring_h1k)(ra,ra)	<< endl << endl;							
										
			cylinder_ring_spectrum_file <<  "%% Vk*Fk: " <<   endl	
						<<  (*cylinder_ring_spectrum_force_Vk)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Wk*W.Fk: " <<   endl	
						<<  (*W.cylinder_ring_spectrum_force_Vk)(ra,ra) << endl << endl;
		}					
	}																
									
}  
  
//
//

void IncFluid::Output_cylinder_ring_spectrum(IncVF& W, IncSF& T)
{	
	
	if (CV_anisotropic_cylinder_switch == 1)
	{ 	

		CV_Compute_cylinder_ring_spectrum();	
		W.CV_Compute_cylinder_ring_spectrum();
		T.CS_Compute_cylinder_ring_spectrum();
		
		Compute_cross_helicity_cylinder_ring_spectrum(W);
		Compute_cross_vT_cylinder_ring_spectrum(T);
			
		CV_Compute_cylinder_ring_spectrum_helicity();	
		W.CV_Compute_cylinder_ring_spectrum_helicity();
		
		Compute_force_cylinder_ring_spectrum();
		W.Compute_force_cylinder_ring_spectrum();
		T.Compute_force_cylinder_ring_spectrum();
			
		if (my_id == master_id)	
		{
			Range ra(1, toEnd);
			
			cylinder_ring_spectrum_file << "%% Time = " << Tnow << endl; 
			
			cylinder_ring_spectrum_file <<  "%% ek1:"		<<   endl	
						<<  (*CV_cylinder_ring_ek1)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% ek2: "		<<   endl	
						<<  (*CV_cylinder_ring_ek2)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Wek1:"		<<   endl	
						<<  (*W.CV_cylinder_ring_ek1)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Wek2: "		<<   endl	
						<<  (*W.CV_cylinder_ring_ek2)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Tek:"	<<   endl	
						<<  (*T.CS_cylinder_ring_ek)(ra,ra)		<< endl << endl;								
																		
			cylinder_ring_spectrum_file <<  "%% dissk: "	<<   endl	
						<<  (*CV_cylinder_ring_dissk)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Wdissk: "	<<   endl	
						<<  (*W.CV_cylinder_ring_dissk)(ra,ra)	<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Tdissk: " <<   endl	
						<<  (*T.CS_cylinder_ring_dissk)(ra,ra)	<< endl << endl;								
										
			cylinder_ring_spectrum_file <<  "%% crosshel:"	<<   endl	
						<<  (*cylinder_ring_ek_cross)(ra,ra)	<< endl << endl;	
										
			cylinder_ring_spectrum_file <<	"%% V1k*Tk: "	<<  endl	
						<<  (*cylinder_ring_ek_cross_V1T)(ra,ra)	<< endl << endl;
																								
			cylinder_ring_spectrum_file <<	"%% V2k*Tk: "	<<  endl	
						<<  (*cylinder_ring_ek_cross_V2T)(ra,ra)	<< endl << endl;
			
			cylinder_ring_spectrum_file <<	"%% V3k*Tk: "	<<  endl	
						<<  (*cylinder_ring_ek_cross_V3T)(ra,ra)	<< endl << endl;																																																		
			
			cylinder_ring_spectrum_file <<  "%% Hk: "		<<   endl	
						<<  (*CV_cylinder_ring_h1k)(ra,ra)			<< endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% W.Hk: "		<<   endl	
						<<  (*W.CV_cylinder_ring_h1k)(ra,ra)		<< endl << endl;							
										
			cylinder_ring_spectrum_file <<  "%% Vk*Fk: " <<   endl	
						<<  (*cylinder_ring_spectrum_force_Vk)(ra,ra) << endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Wk*W.Fk: " <<   endl	
						<<  (*W.cylinder_ring_spectrum_force_Vk)(ra,ra) << endl << endl;
										
			cylinder_ring_spectrum_file <<  "%% Vk*SFk: " <<   endl	
						<<  (*T.cylinder_ring_spectrum_force_SFk)(ra,ra) << endl << endl;
		}					
										
	}																
			

} 

 
/**********************************************************************************************
			
						Output_pressure_spectrum
						Output written in *CS_shell_ek of CSF inherited by NLIN

***********************************************************************************************/
 
void IncFluid::Output_pressure_spectrum()
{

	if (my_id == master_id)	
		pressure_file << "%% Time = " << Tnow << endl; 	

	Compute_shell_spectrum_pressure();		
				
	if (my_id == master_id)	
		pressure_file	<<  "%% pressure shell spectrum:"	<<  endl
					<<	*CS_shell_ek << endl << endl;
	
	if (CV_anisotropic_ring_switch == 1)
	{
		Compute_ring_spectrum_pressure();	
		
		if (my_id == master_id)	
			pressure_file	<<  "%% pressure ring spectrum:"	<<  endl
						<<	*CS_ring_ek << endl << endl;
	}

	
	if (CV_anisotropic_cylinder_switch == 1)
	{		
		Compute_cylinder_ring_spectrum_pressure();
		
		if (my_id == master_id)	
			pressure_file	<<  "%% pressure cylinder ring spectrum:"	<<  endl
						<<  *CS_cylinder_ring_ek << endl << endl;
	}
	
											
}


//*******************************   End of output_energy.cc ***********************************


