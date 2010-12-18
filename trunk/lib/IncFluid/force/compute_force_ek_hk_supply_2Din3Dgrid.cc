force_int_para/* Tarang-4.0
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

/*! \file  compute_force_ek_hk_supply.cc
 * 
 * @brief Compute force when ek and hk supply rate is given
 *
 * @note 2D:   F(k) = alpha * V(k)
 * @note 3D;   F(k) = alpha * V(k) + beta(k) * Omega(k)
 *				alpha and beta determined from the supply rates
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug sk=1, -1 needs to be handled separately.
 */

#include "../IncFluid.h"


//*********************************************************************************************


void IncFluid::Compute_force_const_energy_helicity_supply_2Din3Dgrid()
{

	DP inner_radius = (*force_int_para)(1);
	DP outer_radius = (*force_int_para)(2);
	DP energy_supply = (*force_int_para)(3);
	
	// cout << *V1 << *V2 << *V3 << endl;
	int nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
	DP energy_supply_per_mode = energy_supply / nf;
	
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
		
	
	int lx, ly;
	DP kkmag, alpha_k;
	
	*Force1 = 0.0; 
	*Force2 = 0.0; 
	*Force3 = 0.0;
	
	int kx_min;
	if (basis_type == "FOUR")
		kx_min = -kx_max;
		
	else if (basis_type == "SCFT")
		kx_min = 0;
		
	else
		kx_min = 0;		// for -Wall
		
				
		
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = -ky_max; ky <= ky_max; ky++)  
		{
			lx = Get_lx(basis_type, kx, N);
			ly = Get_ly3D(basis_type, ky, N);
			
			if ( (lx >= 0) && (lx < local_N1) ) 
			{
				kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
				
				if ((kkmag > inner_radius) && (kkmag <= outer_radius)) 
				{
					
					if (CV_Modal_energy(lx, ly, 0) > MYEPS)
					{
						alpha_k = energy_supply_per_mode/ (2*CV_Modal_energy(lx, ly, 0));
						
						Force_alpha_2Din3Dgrid(kx, ky, alpha_k);
					}		
				}
			}	//  of if ((lx >= 0)...)						
		}		// of for
			
	if ( (alias_switch == "DEALIAS") 
			&& (Is_alias_array(basis_type, N, *V1,	outer_radius, kfactor) == 1) )
		Dealias_force();
						
}

//*********************************************************************************************
//
//	Scalar
//

void IncFluid::Compute_force_const_energy_helicity_supply_2Din3Dgrid(IncSF& T)
{

	DP inner_radius = (*force_int_para)(1);
	DP outer_radius = (*force_int_para)(2);
	DP energy_supply = (*force_int_para)(3);
	DP energy_supply_scalar = (*force_int_para)(4);
	
	
	int nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
	DP energy_supply_per_mode = energy_supply / nf;
	DP energy_supply_scalar_per_mode = energy_supply_scalar / nf;
	
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	
	int lx, ly;	
	DP kkmag, alpha_k, alpha_k_scalar;
	

	*Force1 = 0.0; 
	*Force2 = 0.0; 
	*Force3 = 0.0;
	
	*T.Force = 0.0;
	
	int kx_min;
	if (basis_type == "FOUR")
		kx_min = -kx_max;
		
	else if (basis_type == "SCFT")
		kx_min = 0;
	
	else
		kx_min = 0;		// for -Wall
		
						
		
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = -ky_max; ky <= ky_max; ky++)  
		{
			lx = Get_lx(basis_type, kx, N);
			ly = Get_ly3D(basis_type, ky, N);
			
			if ( (lx >= 0) && (lx < local_N1) ) 
			{
				kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
				
				if ((kkmag > inner_radius) && (kkmag <= outer_radius)) 
				{
					if (energy_supply_per_mode > MYEPS)
						if (CV_Modal_energy(lx, ly, 0) > MYEPS)
						{
							alpha_k = energy_supply_per_mode / (2*CV_Modal_energy(lx, ly, 0)); 
							
							Force_alpha_2Din3Dgrid(kx, ky, alpha_k);
						}	
					
					
					if (energy_supply_scalar_per_mode > MYEPS)
						if (T.CS_Modal_energy(lx, ly, 0) > MYEPS)
						{
							alpha_k_scalar = energy_supply_scalar_per_mode 
											/ T.CS_Modal_energy(lx, ly, 0); 
											
							T.Force_scalar_alpha(kx, ky, 0, alpha_k_scalar);
						}
				}
			}   // of ( (lx >= 0)..)						
		}		// of for
		
	if ( (alias_switch == "DEALIAS") 
			&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
		Dealias_force();
						
}

//*********************************************************************************************
// 
// Vector
//
void IncFluid::Compute_force_const_energy_helicity_supply_2Din3Dgrid(IncVF& W)
{

	DP inner_radius = (*force_int_para)(1);
	DP outer_radius = (*force_int_para)(2);
	DP energy_supply = (*force_int_para)(3);
	DP energy_supply_W = (*force_int_para)(4);
	
	int nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
	DP energy_supply_per_mode = energy_supply / nf;
	DP energy_supply_W_per_mode = energy_supply_W / nf;	
	
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
		
	int lx, ly;
	DP kkmag;
	DP alpha_k;
	DP alpha_k_W;
	
	*Force1 = 0.0; 
	*Force2 = 0.0; 
	*Force3 = 0.0;
	*W.Force1 = 0.0; 
	*W.Force2 = 0.0; 
	*W.Force3 = 0.0;
	
	int kx_min;
	if (basis_type == "FOUR")
		kx_min = -kx_max;
		
	else if (basis_type == "SCFT")
		kx_min = 0;
		
	else
		kx_min = 0;		// for -Wall
		
				
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = -ky_max; ky <= ky_max; ky++)  
		{
			lx = Get_lx(basis_type, kx, N);
			ly = Get_ly3D(basis_type, ky, N);
			
			if ( (lx >= 0) && (lx < local_N1) ) 
			{
				kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
				
				if ((kkmag > inner_radius) && (kkmag <= outer_radius)) 
				{
					
					if (energy_supply_per_mode > MYEPS)
						if (CV_Modal_energy(lx, ly, 0) > MYEPS)
						{
							alpha_k = energy_supply_per_mode / (2*CV_Modal_energy(lx, ly, 0));
							
							Force_alpha_2Din3Dgrid(kx, ky, alpha_k);	
						}	
					
					// for W
					
					if (energy_supply_W_per_mode > MYEPS)
						if (W.CV_Modal_energy(lx, ly, 0) > MYEPS)
						{
							alpha_k_W = energy_supply_W_per_mode 
										/ (2*W.CV_Modal_energy(lx, ly, 0)); 
							
							W.Force_alpha_2Din3Dgrid(kx, ky, alpha_k_W);
						}		
				}  
			}	// of ( (lx >= 0)...)					
		}		// of for
			
	if ( (alias_switch == "DEALIAS") 
			&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
		Dealias_force(W);					
}


//*********************************************************************************************
//
// Vector + Scalar
//

void IncFluid::Compute_force_const_energy_helicity_supply_2Din3Dgrid(IncVF& W, IncSF& T)
{

	DP inner_radius = (*force_int_para)(1);
	DP outer_radius = (*force_int_para)(2);
	DP energy_supply = (*force_int_para)(3);
	DP energy_supply_W = (*force_int_para)(4);
	DP energy_supply_scalar = (*force_int_para)(5);
	
	int nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
	DP energy_supply_per_mode = energy_supply / nf;
	DP energy_supply_W_per_mode = energy_supply_W / nf;	
	DP energy_supply_scalar_per_mode = energy_supply_scalar / nf;	
	
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);

		
	
	int lx, ly;
	DP kkmag;
	DP alpha_k;
	DP alpha_k_W;
	DP alpha_k_scalar;
	
	*Force1 = 0.0; 
	*Force2 = 0.0; 
	*Force3 = 0.0;
	
	*W.Force1 = 0.0; 
	*W.Force2 = 0.0; 
	*W.Force3 = 0.0;
	
	*T.Force = 0.0;
	
	int kx_min;
	if (basis_type == "FOUR")
		kx_min = -kx_max;
		
	else if (basis_type == "SCFT")
		kx_min = 0;
	
	else
		kx_min = 0;		// for -Wall
		
		
			
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = -ky_max; ky <= ky_max; ky++)  
		{
			lx = Get_lx(basis_type, kx, N);
			ly = Get_ly3D(basis_type, ky, N);
			
			if ( (lx >= 0) && (lx < local_N1) ) 
			{
				kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
				
				if ((kkmag > inner_radius) && (kkmag <= outer_radius)) 
				{
					if (energy_supply_per_mode > MYEPS)
						if (CV_Modal_energy(lx, ly, 0) > MYEPS)
						{
							alpha_k = energy_supply_per_mode / (2*CV_Modal_energy(lx, ly, 0));
							
							Force_alpha_2Din3Dgrid(kx, ky, alpha_k);
						}
						
					
					if (energy_supply_scalar_per_mode > MYEPS)
						if (T.CS_Modal_energy(lx, ly, 0) > MYEPS)
						{
							alpha_k_scalar = energy_supply_scalar_per_mode 
											/ T.CS_Modal_energy(lx, ly, 0);
											
							T.Force_scalar_alpha(kx, ky, 0, alpha_k_scalar);
						}
						
					// for W
					
					if (energy_supply_W_per_mode > MYEPS)
						if (W.CV_Modal_energy(lx, ly, 0) > MYEPS)
						{
							alpha_k_W = energy_supply_W_per_mode 
											/ (2*W.CV_Modal_energy(lx, ly, 0)); 
							
							W.Force_alpha_2Din3Dgrid(kx, ky, alpha_k_W);
						}
				}
			}	// of if ( (lx >= 0) ..)	
		}		// of for loop
			
	if ( (alias_switch == "DEALIAS") 
			&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
		Dealias_force(W, T);
					
}


//***********************  End of compute_force_ek_hk_supply.cc *******************************


