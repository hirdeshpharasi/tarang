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

/*! \file  init_cond_energy_2Din3Dgrid.cc
 * 
 * @brief Construct energy spectrum given parameter values.  Given spectrum construct
 *			random initial condition (phases).
 *
 * @note 	Given  Model energy spectrum for initial condition
 *		Sk(k) = a k^4 exp(-b k^1.1) / (k^4 + q^4)^(2.8/12)
 *		with q = 1.5, b = 0.02
 *
 *	Notation:  (Ki =) kki = ki * kfactor[i]
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "../IncFluid.h"

extern Uniform<DP> SPECrand;


//*********************************************************************************************


/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters a = (*init_cond_para)(1).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k) is constructed so that it has the given energy but the phases
 *			are random. 
 * 
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_spectrum_2Din3Dgrid()
{

	DP kkmag, ek;
	DP amp, phase;
	int index;

	
	(*V1) = 0.0;
	(*V2) = 0.0;
	(*V3) = 0.0;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_para)(1));
	
//	cout << "Init energy spectrum U: " << *CV_shell_ek << endl;
	
	for (int lx=0; lx<local_N1; lx++)										
		for (int ly=0; ly<N[2]; ly++)   
		{
			kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
			
			if (kkmag > MYEPS)
			{
				index = (int) ceil(kkmag);
				
				ek = (*CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
						
				amp = sqrt(2*ek);
				
				phase = 2*M_PI * SPECrand.random();
				
				Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, amp, phase);
			}
		}	
		
	if (my_id == master_id)
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V3)(0,0,0) = 0.0;
	
	Satisfy_reality_condition_field();

	if (alias_switch == "DEALIAS")		Dealias();
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
void  IncFluid::Init_cond_energy_spectrum_2Din3Dgrid(IncSF& T)
{

	DP kkmag, ek, ekT;
	DP amp, phase;
	DP ampT, phaseT;
	int index;
	
	
	(*V1) = 0.0;
	(*V2) = 0.0;
	(*V3) = 0.0;
	(*T.F) = 0.0;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_para)(1));
	Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_para)(2));
	
	for (int lx=0; lx<local_N1; lx++)										
		for (int ly=0; ly<N[2]; ly++)  
		{
			kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
			
			if (kkmag > MYEPS)
			{
				index = (int) ceil(kkmag);
				
				ek = (*CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
					
				ekT= (*T.CS_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
				
				amp = sqrt(2*ek);
				ampT = sqrt(2*ekT);
				
				phase  = 2*M_PI * SPECrand.random();
				phaseT = 2*M_PI * SPECrand.random();
				
				Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, amp, phase);
				T.Put_scalar_amp_phase_comp_conj(lx, ly, 0, ampT, phaseT);
					
			}	
		}	
		
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V3)(0,0,0) =  0.0;
		(*T.F)(0,0,0) = 0.0;
	}	
		
	Satisfy_reality_condition_field(T);
		

	if (alias_switch == "DEALIAS")		Dealias(T);
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
void  IncFluid::Init_cond_energy_spectrum_2Din3Dgrid(IncVF& W)
{

	DP kkmag, ek, ekW;
	DP amp, phase;
	DP ampW, phaseW;
	int index;
	
	
	(*V1) = 0.0;
	(*V2) = 0.0;
	(*V3) = 0.0;
	
	(*W.V1) = 0.0;
	(*W.V2) = 0.0;
	(*W.V3) = 0.0;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_para)(1));
	Model_initial_shell_spectrum(N, *W.CV_shell_ek, (*init_cond_para)(2));
	
	for (int lx=0; lx<local_N1; lx++)										
		for (int ly=0; ly<N[2]; ly++)  
		{
			kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
			
			if (kkmag > MYEPS)
			{
				index = (int) ceil(kkmag);
				
				ek = (*CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
					
				ekW= (*W.CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
					
				amp = sqrt(2*ek);
				ampW = sqrt(2*ekW);
				
				phase  = 2*M_PI * SPECrand.random();
				phaseW = 2*M_PI * SPECrand.random();
				
				Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, amp, phase);
				W.Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, ampW, phaseW);
			}	
		}	
		
	
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
	}	
		
	Satisfy_reality_condition_field(W);
				

	if (alias_switch == "DEALIAS")		Dealias(W);
	
}


//*********************************************************************************************

/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters: For velocity field a = (*init_cond_para)(1).
 *						For W field a = (*init_cond_para)(2).
 *						For T field a = (*init_cond_para)(3).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k), W(k), T.F(k) is constructed so that it has the given energy 
 *			but the phases are random.
 *  
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_spectrum_2Din3Dgrid(IncVF& W, IncSF& T)
{

	DP kkmag, ek, ekW, ekT;
	DP amp, phase;
	DP ampW, phaseW;
	DP ampT, phaseT;
	int index;
	
	
	(*V1) = 0.0;
	(*V2) = 0.0;
	(*V3) = 0.0;
	
	(*W.V1) = 0.0;
	(*W.V2) = 0.0;
	(*W.V3) = 0.0;
	
	(*T.F) = 0.0;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_para)(1));
	Model_initial_shell_spectrum(N, *W.CV_shell_ek, (*init_cond_para)(2));
	Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_para)(3));
	
	for (int lx=0; lx<local_N1; lx++)										
		for (int ly=0; ly<N[2]; ly++) 
		{
			kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
			
			if (kkmag > MYEPS)
			{
				index = (int) ceil(kkmag);
				
				ek = (*CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
					
				ekW= (*W.CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
							
				ekT= (*T.CS_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);			
					
				amp = sqrt(2*ek);
				ampW = sqrt(2*ekW);
				ampT = sqrt(2*ekT);
				
				phase  = 2*M_PI * SPECrand.random();
				phaseW = 2*M_PI * SPECrand.random();
				phaseT = 2*M_PI * SPECrand.random();
				
				Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, amp, phase);
				W.Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, ampW, phaseW);
				T.Put_scalar_amp_phase_comp_conj(lx, ly, 0, ampT, phaseT);						
			}	
		}	
		
	
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
		(*T.F)(0,0,0) = 0.0;
	}	
		
	Satisfy_reality_condition_field(W, T);

	if (alias_switch == "DEALIAS")		Dealias(W, T);
	
}


//*********************************************************************************************		

void  IncFluid::Init_cond_energy_spectrum_2Din3Dgrid(string Pr_switch, IncSF& T)
{

	if (Pr_switch == "PRZERO") 
	{
		Init_cond_energy_spectrum_2Din3Dgrid();	
					
		*T.F = *V1; 
		Array_divide_ksqr(basis_type, N, *T.F, kfactor);		
	}
	
	else
		Init_cond_energy_spectrum_2Din3Dgrid(T);
		
	Zero_modes_RB_slip(T);	
	
	if (sincos_horizontal_2D_switch == 1)
		Sincos_horizontal(T);		
}	
	
	
//*********************************************************************************************

void  IncFluid::Init_cond_energy_spectrum_2Din3Dgrid(string Pr_switch, IncVF& W, IncSF& T)
{

	if (Pr_switch == "PRZERO") 
	{
		Init_cond_energy_spectrum_2Din3Dgrid(W);	
					
		*T.F = *V1; 
		Array_divide_ksqr(basis_type, N, *T.F, kfactor);		
	}
	
	else
		Init_cond_energy_spectrum_2Din3Dgrid(W, T);
		
	Zero_modes_RB_slip(W, T);		
	
	if (sincos_horizontal_2D_switch == 1)
		Sincos_horizontal(W, T);

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
void  IncFluid::Init_cond_energy_helicity_spectrum_2Din3Dgrid(IncVF& W)
{
	
	DP kkmag, ek, ekW;
	DP amp, phase;
	DP ampW, phaseW;
	int index;
	
	
	(*V1) = 0.0;
	(*V2) = 0.0;
	(*V3) = 0.0;
	
	(*W.V1) = 0.0;
	(*W.V2) = 0.0;
	(*W.V3) = 0.0;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_para)(1));
	Model_initial_shell_spectrum(N, *W.CV_shell_ek, (*init_cond_para)(2));
	DP h = (*init_cond_para)(3);
	
	for (int lx=0; lx<local_N1; lx++)										
		for (int ly=0; ly<N[2]; ly++)  
		{
			kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
			
			if (kkmag > MYEPS)
			{
				index = (int) ceil(kkmag);
				
				ek = (*CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
				
				ekW= (*W.CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
				
				amp = sqrt(2*ek);
				ampW = sqrt(2*ekW);
				
				phase  = 2*M_PI * SPECrand.random();
				phaseW = phase - acos(h);
				
				Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, amp, phase);
				W.Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, ampW, phaseW);
			}	
		}	
	
	
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
	}	
	
	Satisfy_reality_condition_field(W);
	
	
	if (alias_switch == "DEALIAS")		Dealias(W);
	
}


//*********************************************************************************************

/** @brief Create initial condition based on given spectrum Sk(k). 
 *			Parameters: For velocity field a = (*init_cond_para)(1).
 *						For W field a = (*init_cond_para)(2).
 *						For T field a = (*init_cond_para)(3).
 *
 * @note  The energy Sk(k) is divided equally among all the modes in the shell.
 *			Then V(k), W(k), T.F(k) is constructed so that it has the given energy 
 *			but the phases are random.
 *  
 * @note	The mean mode has zero energy.
 */
void  IncFluid::Init_cond_energy_helicity_spectrum_2Din3Dgrid(IncVF& W, IncSF& T)
{
	
	DP kkmag, ek, ekW, ekT;
	DP amp, phase;
	DP ampW, phaseW;
	DP ampT, phaseT;
	int index;
	
	
	(*V1) = 0.0;
	(*V2) = 0.0;
	(*V3) = 0.0;
	
	(*W.V1) = 0.0;
	(*W.V2) = 0.0;
	(*W.V3) = 0.0;
	
	(*T.F) = 0.0;
	
	Model_initial_shell_spectrum(N, *CV_shell_ek, (*init_cond_para)(1));
	Model_initial_shell_spectrum(N, *W.CV_shell_ek, (*init_cond_para)(2));
	DP h = (*init_cond_para)(3);
	
	Model_initial_shell_spectrum(N, *T.CS_shell_ek, (*init_cond_para)(4));
	
	for (int lx=0; lx<local_N1; lx++)										
		for (int ly=0; ly<N[2]; ly++) 
		{
			kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
			
			if (kkmag > MYEPS)
			{
				index = (int) ceil(kkmag);
				
				ek = (*CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
				
				ekW= (*W.CV_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);
				
				ekT= (*T.CS_shell_ek)(index) 
						/ Approx_number_modes_in_shell(basis_type, index, kfactor);			
				
				amp = sqrt(2*ek);
				ampW = sqrt(2*ekW);
				ampT = sqrt(2*ekT);
				
				phase  = 2*M_PI * SPECrand.random();
				phaseW = phase - acos(h);
				
				phaseT = 2*M_PI * SPECrand.random();
				
				Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, amp, phase);
				W.Put_vector_amp_phase_comp_conj_2Din3Dgrid(lx, ly, ampW, phaseW);
				T.Put_scalar_amp_phase_comp_conj(lx, ly, 0, ampT, phaseT);						
			}	
		}	
	
	
	if (my_id == master_id)
	{
		(*V1)(0,0,0) = 	(*V2)(0,0,0) = (*V2)(0,0,0) = 0.0;
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
		(*T.F)(0,0,0) = 0.0;
	}	
	
	Satisfy_reality_condition_field(W, T);
	
	if (alias_switch == "DEALIAS")		Dealias(W, T);
	
}




//********************************** init_cond_energy.cc **************************************


  
