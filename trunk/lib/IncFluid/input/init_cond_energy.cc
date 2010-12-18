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

/*! \file  init_cond_energy.ccinit_cond_double_para
 * 
 * @brief Construct energy spectrum given parameter values.  Given spectrum construct
 *			random initial condition (phases).
 *
 * @note 	Given  Model energy spectrum for initial condition
 *		Sk(k) = a k^4 exp(-b k^1.1) / (k^4 + q^4)^(2.8/12)
 *		with q = 1.5, b = 0.02
 *
 * @note:   Satisfy reality condition is critical here for kz=0 and N[3]/2 planes.   
 *			Do not remove this function.
 *
 *	Notation:  (Ki =) kki = ki * kfactor[i]
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "../IncFluid.h"

extern Uniform<DP> SPECrand;

//*********************************************************************************************


/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters a = (*init_cond_double_para)(1).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy but the phases
 *			are random. 
 * 
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_spectrum()
{

	DP kkmag, ek;
	DP amp, phase1, phase2, phase3;
	int index, maxN3;

	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_double_para)(1));
	
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if (kkmag > MYEPS)
				{
					index = (int) ceil(kkmag);
					
					ek = (*CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
							
					amp = sqrt(2*ek);
					
					phase1 = 2*M_PI * SPECrand.random();
					phase2 = 2*M_PI * SPECrand.random();
					phase3 = 2*M_PI * SPECrand.random();
	
					Put_vector_amp_phase_comp_conj(l1, l2, l3, N, amp, phase1, phase2, phase3);
				}
				
			}
		
	if (my_id == master_id)
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V3)(0,0,0) = 0.0;

	if (alias_switch == "DEALIAS")		Dealias();
	
	if (N[3] == 2)
	{
		(*V1)(Range::all(), Range::all(), 1) = 0.0;
		(*V2)(Range::all(), Range::all(), 1) = 0.0;
		(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
	}
	
	if (N[2] == 1)
	{
		(*V2)(Range::all(), 0, Range::all()) = 0.0;
	}
	
	Satisfy_reality_condition_field();
	
}


//*********************************************************************************************

/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters: For velocity field a = init_cond_para(1).
 *						For scalar field a = init_cond_para(2).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k), F(k) is constructed so that it has the given energy but the phases
 *			are random.
 *  
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_spectrum(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Init_cond_energy_spectrum_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Init_cond_energy_spectrum_RB(T);
	
}


void  IncFluid::Init_cond_energy_spectrum_scalar(IncSF& T)
{

	DP kkmag, ek, ekT;
	DP amp, phase1, phase2, phase3;
	DP ampT, phaseT;
	int index, maxN3;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_double_para)(1));
	Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_double_para)(2));
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if (kkmag > MYEPS)
				{
					index = (int) ceil(kkmag);
					
					ek = (*CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
						
					ekT= (*T.CS_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
					
					amp = sqrt(2*ek);
					ampT = sqrt(2*ekT);
					
					phase1 = 2*M_PI * SPECrand.random();
					phase2 = 2*M_PI * SPECrand.random();
					phase3 = 2*M_PI * SPECrand.random();
					
					phaseT = 2*M_PI * SPECrand.random();
					
					Put_vector_amp_phase_comp_conj(l1, l2, l3, N, amp, phase1, phase2, phase3);
					T.Put_scalar_amp_phase_comp_conj(l1, l2, l3, N, ampT, phaseT);
				}	
			}	
		
	if (my_id == master_id)
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V3)(0,0,0) =  (*T.F)(0,0,0) = 0.0;
		

	if (alias_switch == "DEALIAS")		Dealias(T);
	
	if (N[3] == 2)
	{
		(*V1)(Range::all(), Range::all(), 1) = 0.0;
		(*V2)(Range::all(), Range::all(), 1) = 0.0;
		(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
		
		(*T.F)(Range::all(), Range::all(), 1) = 0.0;
	}
	
	if (N[2] == 1)
	{
		(*V2)(Range::all(), 0, Range::all()) = 0.0;
	}
	
	Satisfy_reality_condition_field(T);
}


//******************************************************	

void  IncFluid::Init_cond_energy_spectrum_RB(IncSF& T)
{
	
	if (globalvar_Pr_switch == "PRZERO") 
	{
		Init_cond_energy_spectrum();	
		
		*T.F = *V1; 
		Array_divide_ksqr(basis_type, N, *T.F, kfactor);		
	}
	
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
		// first fill T(k) modes and then compute u(k) using the constraint.
		DP kkmag, ekT;
		DP ampT, phaseT;
		int index, maxN3;
		
		if (N[3] > 2)
			maxN3 = N[3]/2;
		
		else	// 2D
			maxN3 = 0;
		
		Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_double_para)(1));
		
		for (int i1=0; i1<N[1]; i1++)
			for (int i2=0; i2<N[2]; i2++) 
				for (int i3=0; i3<=maxN3; i3++) 
				{
					kkmag = Kmagnitude(basis_type, i1, i2, i3, N, kfactor);
					
					if (kkmag > MYEPS)
					{
						index = (int) ceil(kkmag);
						
						ekT= (*T.CS_shell_ek)(index) 
								/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
						
						ampT = sqrt(2*ekT);
						
						phaseT = 2*M_PI * SPECrand.random();
						
						T.Put_scalar_amp_phase_comp_conj(i1, i2, i3, N, ampT, phaseT);
					}	
				}	
		
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V3)(0,0,0) =  (*T.F)(0,0,0) = 0.0;
		
		if (alias_switch == "DEALIAS")		Dealias(T);
		
		if (N[3] == 2)
		{
			(*V1)(Range::all(), Range::all(), 1) = 0.0;
			(*V2)(Range::all(), Range::all(), 1) = 0.0;
			(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
			
			(*T.F)(Range::all(), Range::all(), 1) = 0.0;
		}
		
		if (N[2] == 1)
		{
			(*V2)(Range::all(), 0, Range::all()) = 0.0;
		}
		
		Satisfy_reality_array(basis_type, N, *T.F);
		
		Init_cond_Prinfty(T);
	}
	
	else	// PRLARGE or PRSMALL
		Init_cond_energy_spectrum_scalar(T);
	
	
	Zero_modes_RB_slip(T);
}	




//*********************************************************************************************

/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters: For velocity field a = init_cond_para(1).
 *						For W field a = init_cond_para(2).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k), W(k) is constructed so that it has the given energy but the phases
 *			are random.
 *  
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_spectrum(IncVF& W)
{

	DP kkmag, ek, ekW;
	DP amp, phase1, phase2, phase3;
	DP ampW, phase1W, phase2W, phase3W;
	int index, maxN3;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_double_para)(1));
	Model_initial_shell_spectrum(N, *W.CV_shell_ek, (*init_cond_double_para)(2));
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if (kkmag > MYEPS)
				{
					index = (int) ceil(kkmag);
					
					ek = (*CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
							
					ekW= (*W.CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
					
					amp = sqrt(2*ek);
					ampW = sqrt(2*ekW);
					
					phase1 = 2*M_PI * SPECrand.random();
					phase2 = 2*M_PI * SPECrand.random();
					phase3 = 2*M_PI * SPECrand.random();
					
					phase1W = 2*M_PI * SPECrand.random();
					phase2W = 2*M_PI * SPECrand.random();
					phase3W = 2*M_PI * SPECrand.random();
					
					
					Put_vector_amp_phase_comp_conj(l1, l2, l3, N, amp, phase1, phase2, phase3);
					
					W.Put_vector_amp_phase_comp_conj(l1, l2, l3, N, ampW, 
													  phase1W, phase2W, phase3W);
				}	
			}
	
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
	}	
				

	if (alias_switch == "DEALIAS")		Dealias(W);
	
	if (N[3] == 2)
	{
		(*V1)(Range::all(), Range::all(), 1) = 0.0;
		(*V2)(Range::all(), Range::all(), 1) = 0.0;
		(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
		
		(*W.V1)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V2)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
	}
	
	if (N[2] == 1)
	{
		(*V2)(Range::all(), 0, Range::all()) = 0.0;
		
		(*W.V2)(Range::all(), 0, Range::all()) = 0.0;
	}
	
	Satisfy_reality_condition_field(W);
	
}


//*********************************************************************************************

/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters: For velocity field a = (*init_cond_double_para)(1).
 *						For W field a = (*init_cond_double_para)(2).
 *						For T field a = (*init_cond_double_para)(3).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k), W(k), T.F(k) is constructed so that it has the given energy 
 *			but the phases are random.
 *  
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_spectrum(IncVF& W, IncSF& T)
{

	DP kkmag, ek, ekW, ekT;
	DP amp, phase1, phase2, phase3;
	DP ampW, phase1W, phase2W, phase3W;
	DP ampT, phaseT;
	int index, maxN3;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_double_para)(1));
	Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_double_para)(2));
	Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_double_para)(3));
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if (kkmag > MYEPS)
				{
					index = (int) ceil(kkmag);
					
					ek = (*CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
							
					ekW= (*W.CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
							
					ekT= (*T.CS_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
					
					amp = sqrt(2*ek);
					ampW = sqrt(2*ekW);
					ampT = sqrt(2*ekT);
					
					phase1 = 2*M_PI * SPECrand.random();
					phase2 = 2*M_PI * SPECrand.random();
					phase3 = 2*M_PI * SPECrand.random();
					
					phase1W = 2*M_PI * SPECrand.random();
					phase2W = 2*M_PI * SPECrand.random();
					phase3W = 2*M_PI * SPECrand.random();
					
					phaseT = 2*M_PI * SPECrand.random();
					
					
					Put_vector_amp_phase_comp_conj(l1, l2, l3, N, amp, phase1, phase2, phase3);
					
					W.Put_vector_amp_phase_comp_conj(l1, l2, l3, N, ampW, 
													  phase1W, phase2W, phase3W);
													  
					T.Put_scalar_amp_phase_comp_conj(l1, l2, l3, N, ampT, phaseT);								  
				}	
			}	
	
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
		(*T.F)(0,0,0) = 0.0;	
	}
	

	if (alias_switch == "DEALIAS")		Dealias(W, T);
	
	if (N[3] == 2)
	{
		(*V1)(Range::all(), Range::all(), 1) = 0.0;
		(*V2)(Range::all(), Range::all(), 1) = 0.0;
		(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
		
		(*W.V1)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V2)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
	}
	
	if (N[2] == 1)
	{
		(*V2)(Range::all(), 0, Range::all()) = 0.0;
		
		(*W.V2)(Range::all(), 0, Range::all()) = 0.0;
	}
	
	Satisfy_reality_condition_field(W, T);
}


//*********************************************************************************************

/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters a = IC(1); sk = Hk/(k*ek) = IC(2).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_helicity_spectrum()
{


	DP kkmag, ek, amp, phase1, phase2, phase3;

	complx uperp1, uperp2;
	int index, maxN3;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_double_para)(1));
	DP sk = (*init_cond_double_para)(2);
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if (kkmag > MYEPS)
				{
					index = (int) ceil(kkmag);
				
					ek = (*CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
					amp = sqrt(2*ek);
					
					phase1 = 2*M_PI * SPECrand.random();
					phase2 = phase1 + M_PI/2.0;
					phase3 = asin(sk)/2.0;			// zeta
					
					Put_vector_amp_phase_comp_conj(l1, l2, l3, N, amp, phase1, phase2, phase3);
				} // of if(kkmag > MYEPS)			
			}	
	
	if (my_id == master_id)
		(*V1)(0,0,0) = (*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;

	if (alias_switch == "DEALIAS")		Dealias();
	
	if (N[3] == 2)
	{
		(*V1)(Range::all(), Range::all(), 1) = 0.0;
		(*V2)(Range::all(), Range::all(), 1) = 0.0;
		(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
	}
	
	if (N[2] == 1)
	{
		(*V2)(Range::all(), 0, Range::all()) = 0.0;
	}
	
	Satisfy_reality_condition_field();

}



//*********************************************************************************************
/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters a = IC(1); sk = Hk/(k*ek) = IC(2); a_scalar = IC(3),
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */

void  IncFluid::Init_cond_energy_helicity_spectrum(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Init_cond_energy_helicity_spectrum_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Init_cond_energy_helicity_spectrum_RB(T);
}


void  IncFluid::Init_cond_energy_helicity_spectrum_scalar(IncSF& T)
{

	DP kkmag, ek, amp, phase1, phase2, phase3;
	DP ekT, ampT, phaseT;
	complx uperp1, uperp2;
	int index, maxN3;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_double_para)(1));
	DP sk = (*init_cond_double_para)(2);
	
	Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_double_para)(3));
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if (kkmag > MYEPS)
				{
					index = (int) ceil(kkmag);
				
					ek = (*CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
					amp = sqrt(2*ek);
					
					phase1 = 2*M_PI * SPECrand.random();
					phase2 = phase1 + M_PI/2.0;
					phase3 = asin(sk)/2.0;			// zeta
					
					Put_vector_amp_phase_comp_conj(l1, l2, l3, N, amp, phase1, phase2, phase3);
		
					
					// scalar
					
					ekT = (*T.CS_shell_ek)(index) 
								/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
								
					ampT = sqrt(2*ekT);
					phaseT = 2*M_PI * SPECrand.random();
					
					T.Put_scalar_amp_phase_comp_conj(l1, l2, l3, N, ampT, phaseT);
				}	
			}	
	
	if (my_id == master_id)
	{
		(*V1)(0,0,0)  =  (*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*T.F)(0,0,0) = 0.0;
	}	
	
	if (alias_switch == "DEALIAS")		Dealias(T);
	
	if (N[3] == 2)
	{
		(*V1)(Range::all(), Range::all(), 1) = 0.0;
		(*V2)(Range::all(), Range::all(), 1) = 0.0;
		(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
		
		(*T.F)(Range::all(), Range::all(), 1) = 0.0;
	}
	
	if (N[2] == 1)
	{
		(*V2)(Range::all(), 0, Range::all()) = 0.0;
	}
	
	Satisfy_reality_condition_field(T);
	
}


void  IncFluid::Init_cond_energy_helicity_spectrum_RB(IncSF& T)
{
	
	if (globalvar_Pr_switch == "PRZERO") 
	{
		Init_cond_energy_helicity_spectrum();	
		
		*T.F = *V1; 
		Array_divide_ksqr(basis_type, N, *T.F, kfactor);		
	}
	
	else if (globalvar_Pr_switch == "PRINFTY") 
	{
		cout << "Init_cond_energy_helicity_spectrum_RB  NOT possible for PRINFTY condition" 
		<< endl;
		exit(0);
	}
	
	else
		Init_cond_energy_helicity_spectrum_scalar(T);
	
	Zero_modes_RB_slip(T);		
	
}



//*********************************************************************************************
/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters a = IC(1); sk = IC(2); 
 *			for W: a = IC(3), skW =  HM(k) * k/ EW(k) = IC(4), 
 *					h = 2*Hc / (amp, ampW) = IC(5).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_helicity_spectrum(IncVF& W)
{

	DP kkmag, ek, amp, phase1, phase2, phase3;
	DP ekW, ampW, h, phase1W, phase2W, phase3W;
	DP temp;

	complx uperp1, uperp2, wperp1, wperp2;
	int index, maxN3;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_double_para)(1));
	DP sk = (*init_cond_double_para)(2);
	
	Model_initial_shell_spectrum(N, *W.CV_shell_ek, (*init_cond_double_para)(3));
	
	DP skW = (*init_cond_double_para)(4);			// = HM(k) * k / EW(k)
	h = (*init_cond_double_para)(5);			// = 2*Hc/(amp*ampW)
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
			
				if (kkmag > MYEPS)
				{
					index = (int) ceil(kkmag);
				
					ek = (*CV_shell_ek)(index) 
							/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
					amp = sqrt(2*ek);
					
					phase1 = 2*M_PI * SPECrand.random();
					phase2 = phase1 + M_PI/2.0;
					phase3 = asin(sk)/2.0;			// zeta
					
					Put_vector_amp_phase_comp_conj(l1, l2, l3, N, amp, phase1, phase2, phase3);
					
							
					// W field
					
					ekW = (*W.CV_shell_ek)(index) 
								/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
					ampW = sqrt(2*ekW);

					
					phase3W = asin(skW)/2.0;			// zeta_b				
					temp = h / cos(phase3-phase3W);
					// cos(phase1 - phase1W)
					
					phase1W = phase1 - acos(temp);
					phase2W = phase1W + M_PI/2.0;
					
					W.Put_vector_amp_phase_comp_conj(l1, l2, l3, N, ampW, 
													 phase1W, phase2W, phase3W);
				}  // of if (kkmax > MYEPS)			
			}	
		
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
	}	

	if (alias_switch == "DEALIAS")		Dealias(W);
	
	if (N[3] == 2)
	{
		(*V1)(Range::all(), Range::all(), 1) = 0.0;
		(*V2)(Range::all(), Range::all(), 1) = 0.0;
		(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
		
		(*W.V1)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V2)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
	}
	
	if (N[2] == 1)
	{
		(*V2)(Range::all(), 0, Range::all()) = 0.0;
		(*W.V2)(Range::all(), 0, Range::all()) = 0.0;
	}
	
	Satisfy_reality_condition_field(W);
}


//*********************************************************************************************
/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters a = IC(1); sk = IC(2); 
 *			for W: a = IC(3), skW =  HM(k) * k/ EW(k) = IC(4), 
 *					h = 2*Hc / (amp, ampW) = IC(5).
 *			for T: a = IC(6).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy
 * 
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_helicity_spectrum(IncVF& W, IncSF& T)
{

	DP ekT, ampT, phaseT;
	DP kkmag;
	int index, maxN3;
	
	Init_cond_energy_helicity_spectrum(W);
	
	Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_double_para)(6));
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++) 
			{
				if (!((l1 == 0) && (l2 == 0) && (l3 == 0)))
				{
				
					kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
					index = (int) ceil(kkmag);
					
					ekT = (*T.CS_shell_ek)(index) 
								/ Approx_number_modes_in_shell(basis_type, N, index, kfactor);
								
					ampT = sqrt(2*ekT);
					phaseT = 2*M_PI * SPECrand.random();
					
					T.Put_scalar_amp_phase_comp_conj(l1, l2, l3, N, ampT, phaseT);	
				}		
			}	
			
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
		(*T.F)(0,0,0) = 0.0;
	}		
	
	if (alias_switch == "DEALIAS")		Dealias(W, T);	
	
	if (N[3] == 2)
	{
		(*V1)(Range::all(), Range::all(), 1) = 0.0;
		(*V2)(Range::all(), Range::all(), 1) = 0.0;
		(*V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
		
		(*W.V1)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V2)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
		
		(*T.F)(Range::all(), Range::all(), 1) = 0.0;
	}
	
	if (N[2] == 1)
	{
		(*V2)(Range::all(), 0, Range::all()) = 0.0;
		(*W.V2)(Range::all(), 0, Range::all()) = 0.0;
	}
	
	Satisfy_reality_condition_field(W, T);
					
}

//********************************** init_cond_energy.cc **************************************



  
