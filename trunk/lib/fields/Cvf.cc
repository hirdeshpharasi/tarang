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

/*! \file  Cvf.cc
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @sa Cvf.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bugs No known bugs
 */

#include "Cvf.h"

					
/**********************************************************************************************

   Class constructor: Allocates for Vi, Ek etc.; Initializes the arrays.
   
**********************************************************************************************/

CVF::CVF
(
	int *NN, 
	string string_switches[], 
	Array<int,1> switches,
	DP *kfactor, 
	Array<int,1> misc_output_para
) 
{

	CV_basis_type			= string_switches[1];
	CV_alias_switch			= string_switches[2];	
	CV_no_input_field_mode	= string_switches[4];
	CV_no_output_field_mode = string_switches[5];
	
	CV_anisotropic_ring_switch		= switches(4);
	CV_anisotropic_cylinder_switch	= switches(5);
	CV_structure_fn_switch			= switches(6);
	CV_planar_structure_fn_switch	= switches(7);
	
	CV_sector_spectrum_input_scheme = misc_output_para(1);
	CV_no_sectors_spectrum			= misc_output_para(2);
	CV_no_cylinder_slabs_spectrum	= misc_output_para(3);
	CV_structurefn_q_min			= misc_output_para(4);
	CV_structurefn_q_max			= misc_output_para(5);
	
	for (int i=1; i<=3; i++) 
	{
		Ncv[i]=NN[i];  
		CV_kfactor[i] = kfactor[i];
		CV_xfactor[i] = 1/kfactor[i];
	}

	CV_total_energy = CV_total_dissipation = CV_entropy = 0.0;


	if (CV_anisotropic_ring_switch == 1)
	{
		// Angles of the sectors (for 2D & 3D)

		DP max_theta = Get_max_polar_angle(CV_basis_type);		// M_PI/2
		
		CV_sector_angle_array_spectrum = new Array<DP,1>(CV_no_sectors_spectrum+1);
		
		if (CV_sector_spectrum_input_scheme == 0)		// using linear spacing
		{ 
			DP dtheta = max_theta/CV_no_sectors_spectrum;
			
			for (int i=0; i<CV_no_sectors_spectrum; i++)
				(*CV_sector_angle_array_spectrum)(i) = i*dtheta;
					
			(*CV_sector_angle_array_spectrum)(CV_no_sectors_spectrum) = max_theta;
		}
		
		else if (CV_sector_spectrum_input_scheme == 1)
		{
				; // work on it
		}
	}
	

	// Memory allocation now..
	
	V1 = new Array<complx,3>(local_N1, Ncv[2], Ncv[3]/2+1);  
	V2 = new Array<complx,3>(local_N1, Ncv[2], Ncv[3]/2+1);  
	V3 = new Array<complx,3>(local_N1, Ncv[2], Ncv[3]/2+1);   

	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0;
	
	
	CV_shell_ek_size =  Min_radius_outside(CV_basis_type, CV_alias_switch, Ncv, CV_kfactor)+1;
	
	CV_ring_ek_size =  Max_radius_inside(CV_basis_type, CV_alias_switch, Ncv, CV_kfactor)+1;
	// +1 for k=0.

	CV_shell_ek1 = new Array<DP,1>(CV_shell_ek_size);  
	(*CV_shell_ek1) = 0.0;
	
	CV_shell_ek2 = new Array<DP,1>(CV_shell_ek_size);  
	(*CV_shell_ek2) = 0.0;
	
	CV_shell_ek3 = new Array<DP,1>(CV_shell_ek_size);  
	(*CV_shell_ek3) = 0.0;
	
	CV_shell_ek  = new Array<DP,1>(CV_shell_ek_size);  
	(*CV_shell_ek) = 0.0;
	
	CV_shell_dissk=new Array<double,1>(CV_shell_ek_size);  
	(*CV_shell_dissk) = 0.0;
	
	CV_shell_h1k1 = new Array<DP,1>(CV_shell_ek_size);  
	(*CV_shell_h1k1) = 0.0;	
	
	CV_shell_h1k2 = new Array<DP,1>(CV_shell_ek_size);  
	(*CV_shell_h1k2) = 0.0;
	
	CV_shell_h1k3 = new Array<DP,1>(CV_shell_ek_size);  
	(*CV_shell_h1k3) = 0.0;
	

	if (CV_anisotropic_ring_switch == 1)
	{
		CV_ring_ek1  = new Array<double,2>(CV_ring_ek_size, CV_no_sectors_spectrum+1);  
		(*CV_ring_ek1) = 0.0;
		
		CV_ring_ek2  = new Array<double,2>(CV_ring_ek_size, CV_no_sectors_spectrum+1);  
		(*CV_ring_ek2) = 0.0;
		
		CV_ring_h1k  = new Array<double,2>(CV_ring_ek_size, CV_no_sectors_spectrum+1);  
		(*CV_ring_h1k) = 0.0;
		
		CV_ring_dissk = new Array<double,2>(CV_ring_ek_size, CV_no_sectors_spectrum+1); 
		(*CV_ring_dissk) = 0.0;
	}

	if (CV_anisotropic_cylinder_switch == 1)
	{
		CV_cylinder_shell_ek_size 
			= Anis_max_Krho_radius_inside(CV_basis_type, CV_alias_switch, Ncv, CV_kfactor)+1;
				// +1 for Krho = 0.		
		
		// slab heights stored in *CV_cylinder_kpll_array_spectrum.
		// using linear spacing
		
		CV_cylinder_kpll_array_spectrum = new Array<DP,1>(CV_no_cylinder_slabs_spectrum+1);
		
		DP	kkpll_min = 0.0;
		DP	kkpll_max = Anis_max_Kpll(CV_basis_type, CV_alias_switch, Ncv, CV_kfactor);
		
		(*CV_cylinder_kpll_array_spectrum)(0) = kkpll_min;
		
		DP dkkpll = (kkpll_max - kkpll_min)/CV_no_cylinder_slabs_spectrum;
			
		for (int i=1; i<CV_no_cylinder_slabs_spectrum; i++)
			(*CV_cylinder_kpll_array_spectrum)(i) = kkpll_min + i*dkkpll;
			
		(*CV_cylinder_kpll_array_spectrum)(CV_no_cylinder_slabs_spectrum) = kkpll_max;
	

		// Allocate space for cylinderical rings etc.	
		CV_cylinder_ring_ek1 = new Array<double,2>(CV_cylinder_shell_ek_size, 
										CV_no_cylinder_slabs_spectrum+1); 												
		
		CV_cylinder_ring_ek2 = new Array<double,2>(CV_cylinder_shell_ek_size, 
										CV_no_cylinder_slabs_spectrum+1); 
										
		CV_cylinder_ring_dissk = new Array<double,2>(CV_cylinder_shell_ek_size, 
										CV_no_cylinder_slabs_spectrum+1); 	
										
		CV_cylinder_ring_h1k = new Array<double,2>(CV_cylinder_shell_ek_size, 
										CV_no_cylinder_slabs_spectrum+1); 																															


		(*CV_cylinder_ring_ek1) = 0.0;
		(*CV_cylinder_ring_ek2) = 0.0;
		(*CV_cylinder_ring_dissk) = 0.0;
		(*CV_cylinder_ring_h1k) = 0.0;
	}
	
	// Structure function
	
	if (CV_structure_fn_switch == 1)
	{
		CV_structurefn_r_max = Max_radius_inside_real_space(CV_basis_type, Ncv, CV_xfactor);
		
		int q_range = CV_structurefn_q_max - CV_structurefn_q_min;
		
		CV_St = new Array<DP,3>(CV_structurefn_r_max+1, q_range+1, 2);
		(*CV_St) = 0.0;
	}
	
	if (CV_planar_structure_fn_switch == 1)
	{
		if (globalvar_anisotropy_switch == 1)
			CV_structure_fn_rpll_max = Ncv[1];

		else if (globalvar_anisotropy_switch == 2)
			CV_structure_fn_rpll_max = Ncv[2];

		else if (globalvar_anisotropy_switch == 3)
			CV_structure_fn_rpll_max = Ncv[3];
		

		int q_range = CV_structurefn_q_max - CV_structurefn_q_min;
		
		CV_st_planar = new Array<DP,4>(CV_structure_fn_rpll_max, 
									CV_structurefn_r_max+1, q_range+1, 2);
		(*CV_st_planar) = 0.0;		
	}
 
}

/**********************************************************************************************

		   Copies Arrays Vi's to to.Vi

**********************************************************************************************/

void CVF::CV_Copy_to(CVF& to)
{
	*to.V1 = *V1;   
	*to.V2 = *V2;   
	*to.V3 = *V3;
}


/**********************************************************************************************

		 Initialize CVF

**********************************************************************************************/


void CVF::CV_Initialize()
{
	*V1 = 0.0;   
	*V2 = 0.0;   
	*V3 = 0.0;
}


/**********************************************************************************************

  	Init_fftw_plan -- creates plan and puts it in
              global vars fftw_plan

***********************************************************************************************/

void CVF::Init_fftw_plan()  
{

	Init_fftw_plan_array(CV_basis_type, Ncv, *V1);
}

		
/**********************************************************************************************

      Inplace Forward transform of all the components   
	 -- temp_r is temporary location 

**********************************************************************************************/


void CVF::CV_Forward_transform(Array<complx,3> temp_r)					
{
	Forward_transform_array(CV_basis_type, Ncv, *V1, temp_r, 1);			// SFT
	Forward_transform_array(CV_basis_type, Ncv, *V2, temp_r, 0);			// CFT
	Forward_transform_array(CV_basis_type, Ncv, *V3, temp_r, 0);			// CFT	
}


/**********************************************************************************************

      Forward_transform_transpose_order(Vitr) = Vi 
		Vir unchanged
		temp_r = N2  x N1 x N3

**********************************************************************************************/


void CVF::CV_Forward_transform_transpose_order
(
		Array<complx,3> V1r, 
		Array<complx,3> V2r,
		Array<complx,3> V3r, 
		Array<complx,3> temp_r
)
{
	temp_r = V1r;
	Forward_transform_array_transpose_order(CV_basis_type, Ncv, temp_r, *V1, 1);		// SFT
	
	temp_r = V2r;
	Forward_transform_array_transpose_order(CV_basis_type, Ncv, temp_r, *V2, 0);		// CFT
	
	temp_r = V3r;
	Forward_transform_array_transpose_order(CV_basis_type, Ncv, temp_r, *V3, 0);		// CFT
}



/**********************************************************************************************

       Inplace Inverse transform of all the components   
	 -- temp_r is temporary space

**********************************************************************************************/


void CVF::CV_Inverse_transform(Array<complx,3> temp_r)
{
	Inverse_transform_array(CV_basis_type, Ncv, *V1, temp_r, 1);			// ISFT
	Inverse_transform_array(CV_basis_type, Ncv, *V2, temp_r, 0);			// ICFT
	Inverse_transform_array(CV_basis_type, Ncv, *V3, temp_r, 0);			// ICFT	
}



/**********************************************************************************************

       Inverse_transform(*Vi) = Vitr 
		*Vi unchanged
		temp = N1 x N2 x N3

**********************************************************************************************/


void CVF::CV_Inverse_transform_transpose_order
(
		Array<complx,3> V1r, 
		Array<complx,3> V2r,
		Array<complx,3> V3r, 
		Array<complx,3> temp
)
{
	temp = *V1;
	Inverse_transform_array_transpose_order(CV_basis_type, Ncv, temp, V1r, 1);		// ISFT
	
	temp = *V2;
	Inverse_transform_array_transpose_order(CV_basis_type, Ncv, temp, V2r, 0);		// ICFT
	
	temp = *V3;
	Inverse_transform_array_transpose_order(CV_basis_type, Ncv, temp, V3r, 0);		// ICFT
}



/**********************************************************************************************

       Computes the total energy sum[|V(k)|^2/2]
			 total dissipation sum[k^2|V(k)|^2]
			 Entropy

**********************************************************************************************/

void CVF::CV_Compute_totalenergy_diss()
{

	CV_total_energy = Get_total_energy(CV_basis_type, CV_alias_switch, Ncv, *V1)
								+ Get_total_energy(CV_basis_type, CV_alias_switch, Ncv, *V2) 
								+ Get_total_energy(CV_basis_type, CV_alias_switch, Ncv, *V3);
									
	CV_total_dissipation  
			= 2 * (Get_total_Sn(CV_basis_type, CV_alias_switch, Ncv, *V1, 2, CV_kfactor)  
				 + Get_total_Sn(CV_basis_type, CV_alias_switch, Ncv, *V2, 2, CV_kfactor) 
				 + Get_total_Sn(CV_basis_type, CV_alias_switch, Ncv, *V3, 2, CV_kfactor));
				 
}


//
//

void CVF::CV_Compute_entropy()
{
	
	 CV_entropy = Get_entropy(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3);
			
}


/**********************************************************************************************

      Computes the shell energy and dissipation spectra in 3D


**********************************************************************************************/
	

void CVF::CVF::CV_Compute_shell_spectrum()
{

	static Array<DP,1> temp(CV_shell_ek_size);  
	temp = 0.0;

	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, 0, 
								*CV_shell_ek1, CV_kfactor);
	
	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V2, 0, 
								*CV_shell_ek2, CV_kfactor);
	
	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V3, 0, 
								*CV_shell_ek3, CV_kfactor);
	

	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, 2, 
								*CV_shell_dissk, CV_kfactor);
	
	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V2, 2, 
								temp, CV_kfactor);
								
	*CV_shell_dissk = *CV_shell_dissk + temp;
	
	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V3, 2, 
								temp, CV_kfactor);
								
	*CV_shell_dissk = *CV_shell_dissk + temp;
	
	*CV_shell_dissk = 2.0*(*CV_shell_dissk);
}

//
// Cross spectrum result(K) = real(A(K). conj(B(K)))/2
//

void CVF::CV_Compute_shell_spectrum
(
		Array<complx,3> W1, Array<complx,3> W2, Array<complx,3> W3, 
		Array<DP,1> result
)
{

	static Array<DP,1> temp(CV_shell_ek_size);	
	temp = 0.0;
	
	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, W1, 0, result, CV_kfactor);
	
	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V2, W2, 0, temp, CV_kfactor);	
	result = result + temp;		
	
	Compute_shell_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V3, W3, 0, temp, CV_kfactor);
	result = result + temp;
}



/**********************************************************************************************

      Computes the ring spectra in 3D


**********************************************************************************************/


void CVF::CV_Compute_ring_spectrum()
{
	if (CV_anisotropic_ring_switch == 1)
	{
		static Array<DP,2> temp(CV_ring_ek_size, CV_no_sectors_spectrum);  
		temp = 0.0;

		Compute_ring_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3, 
									*CV_sector_angle_array_spectrum, 0, 
									*CV_ring_ek1, *CV_ring_ek2, CV_kfactor);

		Compute_ring_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3, 
									*CV_sector_angle_array_spectrum, 2, 
									*CV_ring_dissk, temp, CV_kfactor);
									
		*CV_ring_dissk = *CV_ring_dissk + temp;	
		
		*CV_ring_dissk = 2.0*(*CV_ring_dissk);
	}	

}

//
// Cross spectrum result(K) = A(K). conj(B(K))
//

void CVF::CV_Compute_ring_spectrum
(
	Array<complx,3> W1, Array<complx,3> W2, Array<complx,3> W3, 
	Array<DP,2> result
)
{
	if (CV_anisotropic_ring_switch == 1)
	{
		static Array<DP,2> temp(CV_ring_ek_size, CV_no_sectors_spectrum);  temp = 0.0;
	
		Compute_ring_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3, W1, W2, W3, 
									*CV_sector_angle_array_spectrum, 0, 
									result, temp, CV_kfactor);
									
		result = result + temp;
	}	
	
}


/**********************************************************************************************

						CYLINDRICAL RING SPECTRUM (for 3D only)

**********************************************************************************************/
							

void CVF::CV_Compute_cylinder_ring_spectrum()
{

	if (CV_anisotropic_cylinder_switch == 1)
	{
		static Array<DP,2> temp(CV_cylinder_shell_ek_size, CV_no_cylinder_slabs_spectrum+1);  
		temp = 0.0;

		Compute_cylinder_ring_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3, 
											*CV_cylinder_kpll_array_spectrum, 0, 
											*CV_cylinder_ring_ek1, *CV_cylinder_ring_ek2,
											CV_kfactor);

		Compute_cylinder_ring_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3, 
											*CV_cylinder_kpll_array_spectrum, 2, 
											*CV_cylinder_ring_dissk, temp, CV_kfactor);
									
		*CV_cylinder_ring_dissk = *CV_cylinder_ring_dissk + temp;	
		
		*CV_cylinder_ring_dissk = 2.0*(*CV_cylinder_ring_dissk);
		
	}		
}

//

void CVF::CV_Compute_cylinder_ring_spectrum
(
	Array<complx,3> W1, Array<complx,3> W2, Array<complx,3> W3, 
	Array<DP,2> result
)	
{	
	
	if (CV_anisotropic_cylinder_switch == 1)
	{
		static Array<DP,2> temp(CV_cylinder_shell_ek_size, CV_no_cylinder_slabs_spectrum+1);  
		temp = 0.0;

		Compute_cylinder_ring_spectrum(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3, 
											W1, W2, W3, *CV_cylinder_kpll_array_spectrum, 0, 
											result, temp, CV_kfactor);

		result = result + temp;
	}	

}	



/**********************************************************************************************

       Computes moaal energy, modal helicity, modal vorticity.


**********************************************************************************************/



DP CVF::CV_Modal_energy(int l1, int l2, int l3)
{
	return (Modal_energy(CV_basis_type, *V1, l1, l2, l3) 
				+ Modal_energy(CV_basis_type, *V2, l1, l2, l3)
				+ Modal_energy(CV_basis_type, *V3, l1, l2, l3));
}				

//
//

void CVF::CV_Modal_vorticity(int l1, int l2, int l3, TinyVector<complx,3> &vorticity)
{
	Compute_Modal_vorticity(CV_basis_type, l1, l2, l3, Ncv, *V1, *V2, *V3, 
								CV_kfactor, vorticity);
}

//
//

DP CVF::CV_Modal_helicity(int l1, int l2, int l3)
{
	return Get_Modal_helicity(CV_basis_type, l1, l2, l3, Ncv, *V1, *V2, *V3, CV_kfactor);
}

/**********************************************************************************************

       Computes total helicity, and helicity spectra


**********************************************************************************************/


void CVF::CV_Compute_total_helicity()
{

	Compute_total_helicity(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3,  
										CV_total_helicity1, CV_total_helicity2,
										CV_total_dissipation_H1, CV_total_dissipation_H2, 
										CV_kfactor);
										
}

//
//

void CVF::CV_Compute_shell_spectrum_helicity()
{

	Compute_shell_spectrum_helicity(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3, 
										*CV_shell_h1k1, *CV_shell_h1k2, *CV_shell_h1k3, 
										CV_kfactor);		
}

void CVF::CV_Compute_ring_spectrum_helicity()
{
	if (CV_anisotropic_ring_switch == 1) 
		Compute_ring_spectrum_helicity(CV_basis_type, CV_alias_switch, Ncv, *V1, *V2, *V3, 
							*CV_sector_angle_array_spectrum, 
							*CV_ring_h1k,  CV_kfactor);
}


void CVF::CV_Compute_cylinder_ring_spectrum_helicity()
{

	if (CV_anisotropic_cylinder_switch == 1)
		Compute_cylinder_ring_spectrum_helicity(CV_basis_type, CV_alias_switch, Ncv, 
								*V1, *V2, *V3,
								*CV_cylinder_kpll_array_spectrum, *CV_cylinder_ring_h1k,
								CV_kfactor);									
}


//*********************************************************************************************

void CVF::CV_Compute_structure_function()
{

	if (CV_structure_fn_switch == 1)	
		Compute_structure_function(CV_basis_type, Ncv, *V1, *V2, *V3, CV_structurefn_q_min, 
									CV_structurefn_q_max, *CV_St, CV_xfactor);		
}


//*********************************************************************************************


void CVF::CV_Compute_planar_structure_function()
{
	if (CV_planar_structure_fn_switch == 1)
		Compute_planar_structure_function(CV_basis_type, Ncv, *V1, *V2, *V3, 
									CV_structurefn_q_min, CV_structurefn_q_max, 
									*CV_st_planar, CV_xfactor);	
																			
}

/**********************************************************************************************

       Dealias CVF


**********************************************************************************************/


void CVF::CV_Dealias()
{
	Dealias_array(CV_basis_type, Ncv, *V1);
	Dealias_array(CV_basis_type, Ncv, *V2);
	Dealias_array(CV_basis_type, Ncv, *V3);
}


/**********************************************************************************************

      Input and Output of all the components (Arrays)


**********************************************************************************************/


void CVF::CV_output(ofstream& file_out, Array<complx,3> temp_array, string format)
{
	Write_data_MPI(CV_basis_type, file_out, Ncv, *V1, temp_array, format);
	Write_data_MPI(CV_basis_type, file_out, Ncv, *V2, temp_array, format);
	Write_data_MPI(CV_basis_type, file_out, Ncv, *V3, temp_array, format);
}

/*
void CVF::CV_output_hdf5(DataSet* dataset1R, DataSet* dataset1C, DataSet* dataset2R,
		DataSet* dataset2C, DataSet* dataset3R, DataSet* dataset3C,
	       	Array<complx,3> temp_array)
{
	Write_data_MPI_hdf5(dataset1R, dataset1C, Ncv, *V1, temp_array);
	Write_data_MPI_hdf5(dataset2R, dataset2C, Ncv, *V1, temp_array);
	Write_data_MPI_hdf5(dataset3R, dataset3C, Ncv, *V1, temp_array);
}
*/
void CVF::CV_output(ofstream& file_out, int Nreduced[], string format)
{
	Write_data_MPI(CV_basis_type, file_out, Ncv, Nreduced, *V1, format );
	Write_data_MPI(CV_basis_type, file_out, Ncv, Nreduced, *V2, format);
	Write_data_MPI(CV_basis_type, file_out, Ncv, Nreduced, *V3, format);
}

void CVF::CV_input(ifstream& file_in, Array<complx,3> temp_array, string format)
{
	Read_data_MPI(CV_basis_type, file_in, Ncv, *V1, temp_array, format);
	Read_data_MPI(CV_basis_type, file_in, Ncv, *V2, temp_array, format);
	Read_data_MPI(CV_basis_type, file_in, Ncv, *V3, temp_array, format);	
}


void CVF::CV_input(ifstream& file_in, int Nreduced[], string format)
{
	Read_data_MPI(CV_basis_type, file_in, Ncv, Nreduced, *V1, format);  
	Read_data_MPI(CV_basis_type, file_in, Ncv, Nreduced, *V2, format);	
	Read_data_MPI(CV_basis_type, file_in, Ncv, Nreduced, *V3, format);
	
}




//************************ END of CVF class Definitions ***************************************




