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

/*! \file  compute_force_const_ek_hk.cc
 * 
 * @brief Compute force when ek and hk are held constant.
 *
 * @note 2D:   F(k) = alpha * V(k)
 * @note 3D;   F(k) = alpha * V(k) + beta(k) * Omega(k)
 *				alpha and beta determined from the supply rates
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 * 
 * @bug sk=1, -1 needs to be handled separately. Maximum helicity state requires 
 *		V field to be maximum helical, which is zero probability state.
 */ 


#include "../IncFluid.h"


//*********************************************************************************************


void IncFluid::Compute_force_const_energy_helicity()
{

	static DP inner_radius = (*force_field_para)(1);
	static DP outer_radius = (*force_field_para)(2);
	static DP energy_level = (*force_field_para)(3);
	static DP h_by_k_E = (*force_field_para)(4);				// Hk/(k*e)
	
	
	static int nf; 
	static DP energy_per_mode;
	
	static int kx_max, ky_max, kz_max;
	
	if (is_force_field_modes_read == 0)
	{
		nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
		
		energy_per_mode = energy_level / nf;	
		
		if (N[1] > 1)
			kx_max = (int) ceil(outer_radius/kfactor[1]);
		else 
			kx_max = 0;
		
		
		if (N[2] > 1)
			ky_max = (int) ceil(outer_radius/kfactor[2]);
		else
			ky_max = 0;
		
		
		if (N[3] > 2)
			kz_max = (int) ceil(outer_radius/kfactor[3]);
		else	
			kz_max = 0;
		
		is_force_field_modes_read = 1;	
		
	}
	
	
	int lx, ly, lz;
	DP kkmag, alpha_k, beta_k, sk, temp1, temp2, temp3;
	

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
			for (int kz = 0; kz <= kz_max; kz++)
			{
				lx = Get_lx(basis_type, kx, N);
				ly = Get_ly3D(basis_type, ky, N);
				lz = kz;
				
				if ( (lx >= 0) && (lx < local_N1) ) 
				{
					kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
					
					if ((kkmag > inner_radius) && (kkmag <= outer_radius)) 
					{
			
						if (energy_per_mode > MYEPS)
							if (CV_Modal_energy(lx, ly, lz) > MYEPS)
							{
								sk = CV_Modal_helicity(lx, ly, lz) 
									/ (kkmag * CV_Modal_energy(lx, ly, lz));
								
								temp1 = sqrt(energy_per_mode / CV_Modal_energy(lx, ly, lz));
								
								if (abs(sk*sk-1) > MYEPS2)
								{	
									temp2 = sqrt((1+h_by_k_E) / (1+sk));
									temp3 = sqrt((1-h_by_k_E) / (1-sk));
									
									alpha_k = (temp1/2) * (temp2 + temp3);
									beta_k =  (temp1/(2*kkmag)) * (temp2 - temp3);
								}
								
								else
								{
									alpha_k = temp1;
									beta_k = 0.0;
								}	
								
								Force_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
							}		
					}
				}						
			}
			
	if ( (alias_switch == "DEALIAS") 
			&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
		Dealias_force();
	
	Satisfy_reality_condition_force_field();

}


//*********************************************************************************************
//
// scalar
//

void IncFluid::Compute_force_const_energy_helicity(IncSF& T)
{

	static DP inner_radius = (*force_field_para)(1);
	static DP outer_radius = (*force_field_para)(2);
	static DP energy_level = (*force_field_para)(3);
	static DP h_by_k_E	= (*force_field_para)(4);				// Hk/(k*e)
	static DP energy_scalar_level = (*force_field_para)(5);
	
	static int nf; 
	static DP energy_per_mode;
	static DP energy_scalar_per_mode;
	
	static int kx_max, ky_max, kz_max;
	
	if (is_force_field_modes_read == 0)
	{
		nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
		
		energy_per_mode = energy_level / nf;	
		energy_scalar_per_mode = energy_scalar_level / nf;
		
		if (N[1] > 1)
			kx_max = (int) ceil(outer_radius/kfactor[1]);
		else 
			kx_max = 0;
		
		
		if (N[2] > 1)
			ky_max = (int) ceil(outer_radius/kfactor[2]);
		else
			ky_max = 0;
		
		
		if (N[3] > 2)
			kz_max = (int) ceil(outer_radius/kfactor[3]);
		else	
			kz_max = 0;
		
		is_force_field_modes_read = 1;	
		
	}
		
	
	int lx, ly, lz;
	DP kkmag, alpha_k, beta_k, alpha_k_scalar, sk, temp1, temp2, temp3;
	
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
			for (int kz = 0; kz <= kz_max; kz++)
			{
				lx = Get_lx(basis_type, kx, N);
				ly = Get_ly3D(basis_type, ky, N);
				lz = kz;
				
				if ( (lx >= 0) && (lx < local_N1) ) 
				{
					kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
					
					if ((kkmag > inner_radius) && (kkmag <= outer_radius)) 
					{
						if (energy_per_mode > MYEPS)
							if (CV_Modal_energy(lx, ly, lz) > MYEPS)
							{
								sk = CV_Modal_helicity(lx, ly, lz) 
										/ (kkmag * CV_Modal_energy(lx, ly, lz));
								
								temp1 = sqrt(energy_per_mode / CV_Modal_energy(lx, ly, lz));
								
								if (abs(sk*sk-1) > MYEPS2)
								{
									temp2 = sqrt((1+h_by_k_E) / (1+sk));
									temp3 = sqrt((1-h_by_k_E) / (1-sk));
									
									alpha_k = (temp1/2) * (temp2 + temp3);
									beta_k =  (temp1/(2*kkmag)) * (temp2 - temp3);
								}
								
								else
								{
									alpha_k = temp1;
									beta_k = 0.0;
								}
								
								Force_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
							}
							
						if (energy_scalar_per_mode > MYEPS)
							if (T.CS_Modal_energy(lx, ly, lz) > MYEPS)
							{
								alpha_k_scalar = sqrt(energy_scalar_per_mode 
													/ T.CS_Modal_energy(lx, ly, lz));
								
								T.Force_scalar_alpha(kx, ky, kz, alpha_k_scalar);
							}
					}	
				}					
			}
			
	if ( (alias_switch == "DEALIAS") 
			&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
			Dealias_force();
	
	Satisfy_reality_condition_force_field(T);

}


//*********************************************************************************************
//
// Vector
//


void IncFluid::Compute_force_const_energy_helicity(IncVF& W)
{

	static DP inner_radius = (*force_field_para)(1);
	static DP outer_radius = (*force_field_para)(2);
	static DP energy_level = (*force_field_para)(3);
	static DP h_by_k_E = (*force_field_para)(4);				// Hk/(k*e)
	static DP energyW_level = (*force_field_para)(5);
	static DP h_by_k_E_W = (*force_field_para)(6);	
	
	
	static int nf; 
	static DP energy_per_mode;
	static DP energyW_per_mode;
	
	static int kx_max, ky_max, kz_max;
	
	if (is_force_field_modes_read == 0)
	{
		nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
		
		energy_per_mode = energy_level / nf;	
		energyW_per_mode = energyW_level / nf;
		
		if (N[1] > 1)
			kx_max = (int) ceil(outer_radius/kfactor[1]);
		else 
			kx_max = 0;
		
		
		if (N[2] > 1)
			ky_max = (int) ceil(outer_radius/kfactor[2]);
		else
			ky_max = 0;
		
		
		if (N[3] > 2)
			kz_max = (int) ceil(outer_radius/kfactor[3]);
		else	
			kz_max = 0;
		
		is_force_field_modes_read = 1;	
		
	}
		
	
	int lx, ly, lz;
	DP kkmag;
	DP alpha_k, beta_k, sk, temp1, temp2, temp3;
	DP alpha_k_W, beta_k_W, skW, temp1W, temp2W, temp3W;
	

	*Force1 = 0.0; *Force2 = 0.0; *Force3 = 0.0;
	*W.Force1 = 0.0; *W.Force2 = 0.0; *W.Force3 = 0.0;
	
	int kx_min;
	
	if (basis_type == "FOUR")
		kx_min = -kx_max;
		
	else if (basis_type == "SCFT")
		kx_min = 0;
		
	else
		kx_min = 0;		// for -Wall
			
		
	for (int kx = kx_min; kx <= kx_max; kx++)
		for (int ky = -ky_max; ky <= ky_max; ky++)  
			for (int kz = 0; kz <= kz_max; kz++)
			{
				lx = Get_lx(basis_type, kx, N);
				ly = Get_ly3D(basis_type, ky, N);
				lz = kz;
				
				if ( (lx >= 0) && (lx < local_N1) ) 
				{
					kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
					
					if ((kkmag > inner_radius) && (kkmag <= outer_radius)) 
					{
						if (energy_per_mode > MYEPS)
							if (CV_Modal_energy(lx, ly, lz) > MYEPS)
							{
								sk = CV_Modal_helicity(lx, ly, lz) 
										/ (kkmag * CV_Modal_energy(lx, ly, lz));
								
								temp1 = sqrt(energy_per_mode / CV_Modal_energy(lx, ly, lz));
								
								if (abs(sk*sk-1) > MYEPS2)
								{
									temp2 = sqrt((1+h_by_k_E) / (1+sk));
									temp3 = sqrt((1-h_by_k_E) / (1-sk));
									
									alpha_k = (temp1/2) * (temp2 + temp3);
									beta_k =  (temp1/(2*kkmag)) * (temp2 - temp3);
								}
								
								else
								{
									alpha_k = temp1;
									beta_k = 0.0;
								}
								
								Force_alpha_beta(kx, ky, kz, alpha_k, beta_k);
							}
						
						// For W 
						
						if (energyW_per_mode > MYEPS)
							if (W.CV_Modal_energy(lx, ly, lz) > 0)
							{
								skW = W.CV_Modal_helicity(lx, ly, lz) 
										/ (kkmag * W.CV_Modal_energy(lx, ly, lz));
								
								temp1W = sqrt(energyW_per_mode 
											/ W.CV_Modal_energy(lx, ly, lz));
								
								if (abs(skW*skW-1) > MYEPS2)
								{
									temp2W = sqrt((1+h_by_k_E_W) / (1+skW));
									temp3W = sqrt((1-h_by_k_E_W) / (1-skW));
									
									alpha_k_W = (temp1W/2) * (temp2W + temp3W);
									beta_k_W =  (temp1W/(2*kkmag)) * (temp2W - temp3W);
								}
								
								else
								{
									alpha_k_W = temp1W;
									beta_k_W = 0.0;
								}	
								
								W.Force_alpha_beta(kx, ky, kz, alpha_k_W, beta_k_W);
							}			
						
					}		
				}					
			}
			
	if ( (alias_switch == "DEALIAS") 
			&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
			Dealias_force(W);
	
	Satisfy_reality_condition_force_field(W);
					
}


//*********************************************************************************************
//
// Vector + scalar
//

void IncFluid::Compute_force_const_energy_helicity(IncVF& W, IncSF& T)
{

	static DP inner_radius = (*force_field_para)(1);
	static DP outer_radius = (*force_field_para)(2);
	static DP energy_level = (*force_field_para)(3);
	static DP h_by_k_E = (*force_field_para)(4);				// Hk/(k*e)
	static DP energyW_level = (*force_field_para)(5);
	static DP h_by_k_E_W = (*force_field_para)(6);	
	static DP energy_scalar_level = (*force_field_para)(7);
	
	
	static int nf; 
	static DP energy_per_mode;
	static DP energyW_per_mode;
	static DP energy_scalar_per_mode;
	
	static int kx_max, ky_max, kz_max;
	
	if (is_force_field_modes_read == 0)
	{
		nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
		
		energy_per_mode = energy_level / nf;	
		energyW_per_mode = energyW_level / nf;
		energy_scalar_per_mode = energy_scalar_level / nf;
		
		if (N[1] > 1)
			kx_max = (int) ceil(outer_radius/kfactor[1]);
		else 
			kx_max = 0;
		
		
		if (N[2] > 1)
			ky_max = (int) ceil(outer_radius/kfactor[2]);
		else
			ky_max = 0;
		
		
		if (N[3] > 2)
			kz_max = (int) ceil(outer_radius/kfactor[3]);
		else	
			kz_max = 0;
		
		is_force_field_modes_read = 1;	
		
	}
		
	
	int lx, ly, lz;
	DP kkmag;
	DP alpha_k, beta_k, sk, temp1, temp2, temp3;
	DP alpha_k_W, beta_k_W, skW, temp1W, temp2W, temp3W;
	DP alpha_k_scalar;
	

	*Force1 = 0.0; *Force2 = 0.0; *Force3 = 0.0;
	*W.Force1 = 0.0; *W.Force2 = 0.0; *W.Force3 = 0.0;
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
			for (int kz = 0; kz <= kz_max; kz++)
			{
				lx = Get_lx(basis_type, kx, N);
				ly = Get_ly3D(basis_type, ky, N);
				lz = kz;
				
				if ( (lx >= 0) && (lx < local_N1) ) 
				{
					kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
				
					if ((kkmag > inner_radius) && (kkmag <= outer_radius)) 
					{
						
						if (energy_per_mode > MYEPS)
							if (CV_Modal_energy(lx, ly, lz) > MYEPS)
							{
								sk = CV_Modal_helicity(lx, ly, lz) 
										/ (kkmag * CV_Modal_energy(lx, ly, lz));
								
								temp1 = sqrt(energy_per_mode / CV_Modal_energy(lx, ly, lz));
								
								if (abs(sk*sk-1) > MYEPS2)
								{	
									temp2 = sqrt((1+h_by_k_E) / (1+sk));
									temp3 = sqrt((1-h_by_k_E) / (1-sk));
									
									alpha_k = (temp1/2) * (temp2 + temp3);
									beta_k =  (temp1/(2*kkmag)) * (temp2 - temp3);
								}
								else
								{
									alpha_k = temp1;
									beta_k = 0.0;
								}
								
								Force_alpha_beta(kx, ky, kz, alpha_k, beta_k);
							}
					
						
						// scalar
						
						if (energy_scalar_per_mode > MYEPS)
							if (T.CS_Modal_energy(lx, ly, lz) > MYEPS)
							{
								alpha_k_scalar = sqrt(energy_scalar_per_mode 
													/ T.CS_Modal_energy(lx, ly, lz)); 
													
								T.Force_scalar_alpha(kx, ky, kz, alpha_k_scalar);
							}
								
							
						// W 
						
						if (energyW_per_mode > MYEPS)
							if (W.CV_Modal_energy(lx, ly, lz) > MYEPS)
							{
								skW = W.CV_Modal_helicity(lx, ly, lz) 
										/ (kkmag * W.CV_Modal_energy(lx, ly, lz));
								
								temp1W = sqrt(energyW_per_mode 
											/ W.CV_Modal_energy(lx, ly, lz));
								
								if (abs(skW*skW-1) > MYEPS2)
								{
									temp2W = sqrt((1+h_by_k_E_W) / (1+skW));
									temp3W = sqrt((1-h_by_k_E_W) / (1-skW));
									
									alpha_k_W = (temp1W/2) * (temp2W + temp3W);
									beta_k_W =  (temp1W/(2*kkmag)) * (temp2W - temp3W);
								}
								
								else
								{
									alpha_k_W = temp1W;
									beta_k_W = 0.0;
								}
								
								W.Force_alpha_beta(kx, ky, kz, alpha_k_W, beta_k_W);
							}	
						
					} // of if ((kmag > ..) 	
				}	// of if (lx ..)				
			} // of for ()
			
	if ( (alias_switch == "DEALIAS") 
			&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
			Dealias_force(W, T);
					
	Satisfy_reality_condition_force_field(W, T);
}

//*****************************  End of compute_force_const_ek_hk.cc **************************


