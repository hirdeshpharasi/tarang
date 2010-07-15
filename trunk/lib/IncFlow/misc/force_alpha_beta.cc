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

/*! \file  force_alpha_beta.cc
 * 
 * @brief Compute Force given alpha and beta
 *
 * @note 2D:   F(k) = alpha * V(k)
 * @note 3D;   F(k) = alpha * V(k) + beta(k) * Omega(k)
 *
 * @sa void IncVF::Put_random_vector_add_conj(int i1, int ly, DP amp)
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug   No known bugs
 */


#include "../IncVF.h"
#include "../IncSF.h"

extern Uniform<DP> SPECrand;



//*********************************************************************************************
/** @brief Assign a force vector at (kx, ky, kz) given alpha and beta,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha, beta
 *
 *  @return  F(k) = alpha * V(k)  + beta(k) * Omega(k)
 */
void IncVF::Force_alpha_beta(int kx, int ky, int kz, DP alpha, DP beta)
{

	TinyVector<complx,3> vorticity;
	
	int lx = Get_lx(basis_type, kx, N);  
	int ly = Get_ly3D(basis_type, ky, N); 
	int lz = kz;
	
	CV_Modal_vorticity(lx, ly, lz, vorticity);
	
	if (basis_type == "FOUR")	
	{
		if ( (lx >= 0) && (lx < local_N1) ) 
		{
			(*Force1)(lx, ly, lz) = alpha * (*V1)(lx, ly, lz) 
										+ beta * vorticity(0);
										
			(*Force2)(lx, ly, lz) = alpha * (*V2)(lx, ly, lz) 
										+ beta * vorticity(1);
										
			(*Force3)(lx, ly, lz) = alpha * (*V3)(lx, ly, lz) 
										+ beta * vorticity(2);
		}								
	}								
	
	else if (basis_type == "SCFT")
	{
		if ( (lx >= 0) && (lx < local_N1) ) 
		{
			// in FOUR basis
			(*Force1)(lx, ly, lz) = alpha * (-I) * (*V1)(lx, ly, lz) 
										+ beta * vorticity(0);
										
			(*Force2)(lx, ly, lz) = alpha * (*V2)(lx, ly, lz) 
										+ beta * vorticity(1);
										
			(*Force3)(lx, ly, lz) = alpha * (*V3)(lx, ly, lz) 
										+ beta * vorticity(2);
		
			// for SCFT basis conversion
			(*Force1)(lx, ly, lz) = I * (*Force1)(lx, ly, lz);		
			// because of V1(x) = V1(k..) sin(kx*x)*exp(i ky y)
		}		
	}	
							
	
	// Complex conj  is added for kz=0 & Nz/2 modes collectively 
	//  in Compute_force_const_energy_helicity_supply()
	// while calling Satisfy_reality_condition_force_field().

}


//*********************************************************************************************
/** @brief Assign a force vector at (kx, ky, kz) given alpha and beta,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha, beta
 *
 *  @return  F(k) = alpha * V(k)  + beta(k) * Omega(k)
 */
void IncVF::Force_alpha_2Din3Dgrid(int kx, int ky, DP alpha)
{
	
	TinyVector<complx,3> vorticity;
	
	int lx = Get_lx(basis_type, kx, N);  
	int ly = Get_ly3D(basis_type, ky, N); 

	
	if (basis_type == "FOUR")	
	{
		if ( (lx >= 0) && (lx < local_N1) ) 
		{
			(*Force1)(lx, ly, 0) = alpha * (*V1)(lx, ly, 0); 
			(*Force2)(lx, ly, 0) = alpha * (*V2)(lx, ly, 0); 
		}								
	}								
	
	else if (basis_type == "SCFT")
	{
		if ( (lx >= 0) && (lx < local_N1) ) 
		{
			(*Force1)(lx, ly, 0) = alpha * (*V1)(lx, ly, 0); // (-I)* I =1 Compare with 3D
			(*Force2)(lx, ly, 0) = alpha * (*V2)(lx, ly, 0); 
		}		
	}	
	
	
	// Complex conj  is added for kz=0 & Nz/2 modes collectively 
	//  in Compute_force_const_energy_helicity_supply()
	// while calling Satisfy_reality_condition_force_field().
	
}





//*********************************************************************************************

/** @brief Assign a force vector at (kx, ky, kz) given alpha,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha
 *
 *  @return  F(k) = alpha * V(k).
 *  @note in SCFT basis, conversion to FOUR and back to SCFT implies multiplication
 *			of (-I) and (I).  Hence the conversion remains unchanged.
 */
void IncSF::Force_scalar_alpha(int kx, int ky, int kz, DP alpha)
{
	
	int lx = Get_lx(ISF_basis_type, kx, NIs);  
	int ly = Get_ly3D(ISF_basis_type, ky, NIs); 
	int lz = kz;
	
	if ( (lx >= 0) && (lx < local_N1) ) 
		(*Force)(lx, ly, lz) = alpha * (*F)(lx, ly, lz);
						
	
	// Complex conj  is added for kz=0 & Nz/2 modes collectively 
	//  in Compute_force_const_energy_helicity_supply()
	// while calling Satisfy_reality_condition_force_field().					
		
}



//********************************** force_alpha_beta.cc **************************************

	
