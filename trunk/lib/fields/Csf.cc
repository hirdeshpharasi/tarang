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



/*! \file  Csf.cc
 * 
 * @brief  Class declaration of Csf, a Complex Scalar Field 
 * @sa Csf.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bugs No known bugs
 */

#include "Csf.h" 


//*********************************************************************************************

CSF::CSF
(
	int *NN, 
	string string_switches[], 
	Array<int,1> switches,
	DP *kfactor, 
	Array<int,1> misc_output_para
) 
{

	CS_basis_type			= string_switches[1];
	CS_alias_switch			= string_switches[2];	
	CS_no_input_field_mode	= string_switches[4];
	CS_no_output_field_mode = string_switches[5];
	
	CS_anisotropic_ring_switch		= switches(4);
	CS_anisotropic_cylinder_switch	= switches(5);
	CS_structure_fn_switch			= switches(6);
	CS_planar_structure_fn_switch	= switches(7);
	
	CS_sector_spectrum_input_scheme = misc_output_para(1);
	CS_no_sectors_spectrum			= misc_output_para(2);
	CS_no_cylinder_slabs_spectrum	= misc_output_para(3);
	CS_structurefn_q_min			= misc_output_para(4);
	CS_structurefn_q_max			= misc_output_para(5);
	
	
	for (int i=1; i<=3; i++) 
	{
		Ncs[i] = NN[i];  
		CS_kfactor[i] = kfactor[i];
		CS_xfactor[i] = 1/kfactor[i];
	}

	CS_total_energy = CS_total_dissipation = 0.0;
 
	// Angular sectors...
	
	if (CS_anisotropic_ring_switch == 1)
	{
		// Angles of the sectors
		
		DP max_theta = Get_max_polar_angle(CS_basis_type);
		
		CS_sector_angle_array_spectrum = new Array<DP,1>(CS_no_sectors_spectrum+1);
		
		if (CS_sector_spectrum_input_scheme == 0)		// using linear spacing
		{ 
			DP dtheta = max_theta/CS_no_sectors_spectrum;
			
			for (int i=0; i<CS_no_sectors_spectrum; i++)
				(*CS_sector_angle_array_spectrum)(i) = i*dtheta;
					
			(*CS_sector_angle_array_spectrum)(CS_no_sectors_spectrum) = max_theta;
		}
		
		else if (CS_sector_spectrum_input_scheme == 1)
		{
				;// work on it
		}
	}
	
	
	F = new Array<complx,3>(local_N1, Ncs[2], Ncs[3]/2+1);  
	*F = 0.0;   // initialize

	CS_shell_ek_size =  Min_radius_outside(CS_basis_type, CS_alias_switch, Ncs, CS_kfactor)+1;
	
	CS_ring_ek_size =  Max_radius_inside(CS_basis_type, CS_alias_switch, Ncs, CS_kfactor)+1;	
	// +1 for k=0.
  
	CS_shell_ek = new Array<DP,1>(CS_shell_ek_size); 
	(*CS_shell_ek)=0.0;
	
	CS_shell_dissk = new Array<double,1>(CS_shell_ek_size);  
	(*CS_shell_dissk)=0.0;

	if (CS_anisotropic_ring_switch == 1)
	{
		CS_ring_ek  = new Array<double,2>(CS_ring_ek_size, CS_no_sectors_spectrum+1);  
		(*CS_ring_ek) = 0.0;
		
		CS_ring_dissk = new Array<double,2>(CS_ring_ek_size, CS_no_sectors_spectrum+1); 
		(*CS_ring_dissk) = 0.0;
	}
	
	//
	//
	
	if (CS_anisotropic_cylinder_switch == 1)
	{
		CS_cylinder_shell_ek_size 
			= Anis_max_Krho_radius_inside(CS_basis_type, CS_alias_switch, Ncs, CS_kfactor)+1;
		
		// slab heights stored in *CV_cylinder_kpll_array_spectrum.
		// using linear spacing
		
		CS_cylinder_kpll_array_spectrum = new Array<DP,1>(CS_no_cylinder_slabs_spectrum+1);
		
		DP	kkpll_min = 0.0;
		DP	kkpll_max = Anis_max_Kpll(CS_basis_type, CS_alias_switch, Ncs, CS_kfactor);
		
		(*CS_cylinder_kpll_array_spectrum)(0) = kkpll_min;
		
		DP dkkpll = (kkpll_max - kkpll_min)/CS_no_cylinder_slabs_spectrum;
			
		for (int i=1; i<CS_no_cylinder_slabs_spectrum; i++)
			(*CS_cylinder_kpll_array_spectrum)(i) = kkpll_min + i*dkkpll;
			
		(*CS_cylinder_kpll_array_spectrum)(CS_no_cylinder_slabs_spectrum) = kkpll_max;
		
		// Allocate memory for cylinder now.		
		CS_cylinder_ring_ek = new Array<double,2>(CS_cylinder_shell_ek_size, 
									CS_no_cylinder_slabs_spectrum+1); 
									
		CS_cylinder_ring_dissk = new Array<double,2>(CS_cylinder_shell_ek_size, 
									CS_no_cylinder_slabs_spectrum+1);	
									
		(*CS_cylinder_ring_ek) = 0.0;
		(*CS_cylinder_ring_dissk) = 0.0;	
	}
	
		
	// Structure function
	
	
	if (CS_structure_fn_switch == 1)
	{
		CS_structurefn_r_max = Max_radius_inside_real_space(CS_basis_type, Ncs, CS_xfactor);
		
		int q_range = CS_structurefn_q_max - CS_structurefn_q_min;
		
		CS_St = new Array<DP,2>(CS_structurefn_r_max+1, q_range+1);
		(*CS_St) = 0.0;
	}
	
	
	if (CS_planar_structure_fn_switch == 1)
	{
#ifdef ANISDIRN1
		CS_structure_fn_rpll_max = Ncs[1];
#endif

#ifdef ANISDIRN2
		CS_structure_fn_rpll_max = Ncs[2];
#endif

#ifdef ANISDIRN3
		CS_structure_fn_rpll_max = Ncs[3];
#endif	
	
		int q_range = CS_structurefn_q_max - CS_structurefn_q_min;
		
		CS_st_planar = new Array<DP,3>(CS_structure_fn_rpll_max,
										CS_structurefn_r_max+1, q_range+1);
		(*CS_st_planar) = 0.0;
	}

}




/**********************************************************************************************

		 Initialize CSF

**********************************************************************************************/

void CSF::CS_Initialize()
{
	
	*F = 0.0;
}



/**********************************************************************************************

		Forward transform(*F) -> *F
			-- Here temp_r is temporary array
		
		Forward_transform_transpose_order(*Ftr) -> *F 
		  -- Here Ftr unchanged.
	
***********************************************************************************************/

void CSF::CS_Forward_transform(Array<complx,3> temp_r)
{
	Forward_transform_array(CS_basis_type, Ncs, *F, temp_r, 1);			// SFT
}

//
//

void CSF::CS_Forward_transform_transpose_order(Array<complx,3> Ftr, Array<complx,3> temp_r)
{
	temp_r = Ftr;
	
	Forward_transform_array_transpose_order(CS_basis_type, Ncs, temp_r, *F, 1);	// SFT
}


/**********************************************************************************************

		Inverse transform(*F) -> *F
			-- Here temp_r is temporary array
		
		Inverse_transform_transpose_order(*F) -> *Ftr    

***********************************************************************************************/

void CSF::CS_Inverse_transform(Array<complx,3> temp_r)
{
	Inverse_transform_array(CS_basis_type, Ncs, *F, temp_r, 1);			// SFT
}


void CSF::CS_Inverse_transform_transpose_order(Array<complx,3> Ftr, Array<complx,3> temp)
{	
	temp = *F;
	
	Inverse_transform_array_transpose_order(CS_basis_type, Ncs, temp, Ftr, 1);	// SFT
}


/**********************************************************************************************

		  F(k) = F(k)/k^2

**********************************************************************************************/

void CSF::CS_divide_ksqr()
{

	Array_divide_ksqr(CS_basis_type, Ncs, *F, CS_kfactor);

}


/**********************************************************************************************

		 Computes the total energy sum[|F(k)|^2/2]
			 dissipation sum[k^2|F(k)|^2]
			 Entropy

**********************************************************************************************/

void CSF::CS_Compute_totalenergy_diss()
{
		
	CS_total_energy = Get_total_energy(CS_basis_type, CS_alias_switch, Ncs, *F);
	
	CS_total_dissipation = 2 * Get_total_Sn(CS_basis_type, CS_alias_switch, 
												Ncs, *F, 2, CS_kfactor);
}


//
//

void CSF::CS_Compute_entropy()
{
	CS_entropy = Get_entropy_scalar(CS_basis_type, CS_alias_switch, Ncs, *F);
}


/**********************************************************************************************

      Computes the shell spectra and ring spectra

**********************************************************************************************/


void CSF::CS_Compute_shell_spectrum()
{
	Compute_shell_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, 0, 
								*CS_shell_ek, CS_kfactor);
								
	Compute_shell_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, 2, 
								*CS_shell_dissk, CS_kfactor);
	
	*CS_shell_dissk = 2.0 * (*CS_shell_dissk);		
}

// Computes the Cross spectrum in  3D

void CSF::CS_Compute_shell_spectrum(Array<complx,3> G, Array<DP,1> result)
{
	Compute_shell_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, G, 0, result, CS_kfactor);		
}

//*********************************************************************************************

void CSF::CS_Compute_ring_spectrum()
{

	if (CS_anisotropic_ring_switch == 1)
	{
		Compute_ring_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, 
							*CS_sector_angle_array_spectrum, 0, *CS_ring_ek, CS_kfactor);
									
		Compute_ring_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, 
							*CS_sector_angle_array_spectrum, 2, *CS_ring_dissk, CS_kfactor);
		
		*CS_ring_dissk = 2.0 * (*CS_ring_dissk);
	}	
}

//
//

void CSF::CS_Compute_ring_spectrum(Array<complx,3> G, Array<DP,2> result)
{

	if (CS_anisotropic_ring_switch == 1)
		Compute_ring_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, G,
							*CS_sector_angle_array_spectrum, 0, result, CS_kfactor);									
}	

//*********************************************************************************************

void CSF::CS_Compute_cylinder_ring_spectrum()
{
	if (CS_anisotropic_cylinder_switch == 1)
	{
		Compute_cylinder_ring_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, 
											*CS_cylinder_kpll_array_spectrum,
											0.0, *CS_cylinder_ring_ek, CS_kfactor);
											
		Compute_cylinder_ring_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, 
											*CS_cylinder_kpll_array_spectrum,
											2.0, *CS_cylinder_ring_dissk, CS_kfactor);
	}																				
}

//
//
void CSF::CS_Compute_cylinder_ring_spectrum(Array<complx,3> G, Array<DP,2> result)
{

	if (CS_anisotropic_cylinder_switch == 1)
		Compute_cylinder_ring_spectrum(CS_basis_type, CS_alias_switch, Ncs, *F, G,
											*CS_cylinder_kpll_array_spectrum,
											0, result, CS_kfactor);
}


/**********************************************************************************************

		 Modal energy

**********************************************************************************************/

DP CSF::CS_Modal_energy(int i1, int i2, int i3)
{
	return Modal_energy(CS_basis_type, *F, i1, i2, i3);	
}


/**********************************************************************************************

		 Structure function

**********************************************************************************************/

/*
void CSF::CS_Compute_structure_function()
{
	if (CS_structure_fn_switch == 1)
		Compute_structure_function(CS_basis_type, Ncs, *F, CS_structurefn_q_min, 
				CS_structurefn_q_max, *CS_St, CS_xfactor);				
}


void CSF::CS_Compute_planar_structure_function()
{
	if (CS_planar_structure_fn_switch == 1)
		Compute_planar_structure_function(CS_basis_type, Ncs, *F, CS_structurefn_q_min, 
				CS_structurefn_q_max, *CS_st_planar, CS_xfactor);		 
}

*/
/**********************************************************************************************

		 Dealias CSF

**********************************************************************************************/

void CSF::CS_Dealias()
{ 
	Dealias_array(CS_basis_type, Ncs, *F);
}


/**********************************************************************************************

		 Input and Output array

**********************************************************************************************/

void CSF::CS_output(ofstream& file_out, Array<complx,3> temp_array)
{
	Write_data_MPI(file_out, Ncs, *F, temp_array);
}

void CSF::CS_output(ofstream& file_out, int Nreduced[], Array<complx,3> temp_array)
{
	Write_data_MPI(file_out, Ncs, Nreduced, *F, temp_array);
}

void CSF::CS_input(ifstream& file_in, Array<complx,3> temp_array)
{
	Read_data_MPI(file_in, Ncs, *F, temp_array);
}

void CSF::CS_input(ifstream& file_in, int Nreduced[], Array<complx,3> temp_array)
{
	Read_data_MPI(file_in, Ncs, Nreduced, *F, temp_array);
}



//******************************** End of CSF.cc **********************************************


