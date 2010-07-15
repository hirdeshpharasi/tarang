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

/*! \file  ic_amp_phase.cc
 * 
 * @brief Initial conditions where fields are chosen randomly.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "../IncVF.h"
#include "../IncSF.h"

extern Uniform<DP> SPECrand;


//*********************************************************************************************
/** @brief Put a velocity vector at (lx,ly,lz) given its amp and phases
 *			and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param lx, ly, lz
 *	@param amp		Amplitude of the vector.
 *	@param phase1, phase2, phase3	Phases
 *
 */
void IncVF::Put_vector_amp_phase_comp_conj
(
	int lx, int ly, int lz, 
	DP amp, 
	DP phase1, DP phase2, DP phase3
)
{

	complx vpll, vh1, vh2;
			
	complx uperp1 = amp * exp(I * phase1) * cos(phase3);
	complx uperp2 = amp * exp(I * phase2) * sin(phase3);
	
	DP  theta = AnisKvect_polar_angle(basis_type, lx, ly, lz, N, kfactor);
	DP  phi	  = AnisKvect_azimuthal_angle(basis_type, lx, ly, lz, N, kfactor);
	
	DP kkmag  = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
	DP kkperp = AnisKperp(basis_type, lx, ly, lz, N, kfactor);
	
	if (kkmag > MYEPS)
		if ( kkperp > MYEPS)
		{	
			vpll = -uperp2 * sin(theta);
			vh1 =  uperp2 * cos(theta)*cos(phi) + uperp1 * sin(phi); 
			vh2 =  uperp2 * cos(theta)*sin(phi) - uperp1 * cos(phi); 	
		}
		
		else	// k along x axis.  V on the y-z plane.
		{				
			vpll = 0.0;
			vh1 = uperp1;
			vh2 = uperp2;
		}
	
	
	if ( (lx >= 0) && (lx < local_N1) ) 
	{
#ifdef ANISDIRN1
		(*V1)(lx, ly, lz) = vpll;
		(*V2)(lx, ly, lz) = vh1;
		(*V3)(lx, ly, lz) = vh2;
#endif

#ifdef ANISDIRN2
		(*V2)(lx, ly, lz) = vpll;
		(*V3)(lx, ly, lz) = vh1;
		(*V1)(lx, ly, lz) = vh2;
#endif

#ifdef ANISDIRN3
		(*V3)(lx, ly, lz) = vpll;
		(*V1)(lx, ly, lz) = vh1;
		(*V2)(lx, ly, lz) = vh2;
#endif
	}	
	
				
	if ( (lx >= 0) && (lx < local_N1) ) 									
		if (basis_type == "SCFT")
			(*V1)(lx, ly, lz) = I * (*V1)(lx, ly, lz);		
			// because of V1(x) = V1(k..) sin(kx*x)*exp(i ky y)
		
	// Complex conj  is added for kz=0 & Nz/2 modes collectively 
	//  in Init_cond_energy_spectrum_2Din3Dgrid()
	// while calling Satisfy_reality_condition_field()
	
}



//*********************************************************************************************

/** @brief Put a scalar field at (lx,ly,lz) given its amplitude and phase, 
 *			and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param lx, ly, lz
 *	@param amp  Amplitude of the vector.
 *
 */
void IncSF::Put_scalar_amp_phase_comp_conj(int lx, int ly, int lz, DP amp, DP phase)
{
	
	if ( (lx >= 0) && (lx < local_N1) ) 
		(*F)(lx, ly, lz) = amp * exp(I * phase);
	
	
	if (ISF_basis_type == "SCFT")
		if ( (lx >= 0) && (lx < local_N1) ) 
			(*F)(lx, ly, lz) = I * (*F)(lx, ly, lz);	
			// because of F(x) = F(kx,ky) sin(kx*x)*exp(i ky y)
		
	
	// Complex conj  is added for kz=0 & Nz/2 modes collectively 
	//  in Init_cond_energy_spectrum_2Din3Dgrid()
	// while calling Satisfy_reality_condition_field()		
	
}


//*********************************************************************************************
/** @brief Put a velocity vector at (i1,i2,0) given its amp and phases 
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param i1, i2
 *	@param amp  Amplitude of the vector.
 *
 *  @note  phi is the angle V vector makes with y axis.
 */
void IncVF::Put_vector_amp_phase_comp_conj_2Din3Dgrid(int lx, int ly, DP amp, DP phase)
{
	
	complx vpll, vh;
	
	complx uperp = amp * exp(I * phase);
	
	DP phi = AnisKvect_polar_angle_2Din3Dgrid(basis_type, lx, ly, 0, N, kfactor);
	
	DP kkmag = Kmagnitude(basis_type, lx, ly, 0, N, kfactor);
	
	if (kkmag > MYEPS)
	{
		vpll = -uperp * sin(phi);
		vh = uperp * cos(phi);
	}
	
	if ( (lx >= 0) && (lx < local_N1) ) 
	{
#ifdef ANISDIRN1
		(*V1)(lx, ly, 0) = vpll;
		(*V2)(lx, ly, 0) = vh;
#endif
	
#ifdef ANISDIRN2
		(*V2)(lx, ly, 0) = vpll;
		(*V1)(lx, ly, 0) = vh;
#endif
	}	
	
	if ( (lx >= 0) && (lx < local_N1) ) 
		if (basis_type == "SCFT")
			(*V1)(lx, ly, 0) = I * (*V1)(lx, ly, 0);			
		// because of V1(x) = V1(kx,ky) sin(kx*x)*exp(i ky y)
	
	
	// Complex conj  is added for kz=0 & Nz/2 modes collectively 
	//  in Init_cond_energy_spectrum_2Din3Dgrid()
	// while calling Satisfy_reality_condition_field().
	// For 2Din3D grid, kz=1 is also worked on, but it does not matter since they are all 0s.
	
}



//********************************** ic_random.cc *********************************************

	
