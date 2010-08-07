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


/*! \file four_inline.h 
 * 
 * @brief Inline functions to compute array indices \f$ \vec{i} \f$ given wavenumber  
 * \f$ \vec{k} \f$ and viceversa.
 * Also contains other useful functions like Kmagnitude, Max radius etc.
 *
 * Grid wavenumber \f$ \vec{k} \f$ is integer, that can be computed using the grid index 
 * \f$ \vec{i} \f$.
 * Actual wavenumber \f$ \vec{K} \f$ is computed from the grid wavenumber \f$ \vec{k} \f$ 
 *		using \f$ K_i = k_i * f_i \f$
 * where \f$ f_i \f$ is the kfactor[i]. 
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber 
 *	or gridwavenumber depending on the switch WAVENOACTUAL or WAVENOGRID.   
 *	Actual wavenumber is preferred.
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date	August 2008
 * @bug		No known bugs
 *
 */ 
 

#ifndef _FOUR_INLINE_H
#define _FOUR_INLINE_H

#include <mpi.h>
#include <fftw3-mpi.h>

#include <blitz/array.h>
#include <complex>
#include <cmath>

#include "../basis_basicfn/basis_basicfn_inline.h"
#include "../basis_basicfn/basis_basicfn.h"


using namespace blitz;

//global vars
extern int my_id;							/**<	External: My process id	*/
extern int numprocs;						/**<	External: No of processors	*/
extern const int master_id;					/**<	External: Id of master proc	*/
extern ptrdiff_t local_N1;					/**<	External: local_N1 
															= N1/numprocs in the currentproc */
extern ptrdiff_t local_N1_start;			/**< External: start of i1 in the currentproc */
extern MPI_Status status;


//*********************************************************************************************

/*! @brief	Get grid waveno kx given first local array index l1.
 * 
 *  if i1 <=Ni/2, ki=i1; else kx=i1-Ni.  <BR>
 *	i1= local_N1_start + l1.
 * 
 * \param l1  first local index of an array
 * \return kx corresponding to l1
 */
inline int Get_kx_FOUR(int l1, int N[])  
{ 			
	return  ((local_N1_start + l1) <= N[1]/2) ? 
				(local_N1_start + l1) : (local_N1_start + l1-N[1]); 
} 
 

/*! @brief	Get local array index lx given grid waveno kx.
 * 
 *  If kx>=0, i1=kx; else i1=kx+N1.  <BR>
 *	i1= local_N1_start + l1.
 * 
 * \param	kx  grid wavenumber along x
 * \return	l1  local array index along x
 */
inline int Get_lx_FOUR(int kx, int N[])  
{ 
	return  (kx >= 0) ? (kx-local_N1_start) : (kx + N[1]-local_N1_start);
}	


/*! @brief	Get array index ix given grid waveno kx.
 * 
 *  If kx>=0, i1=kx; else i1=kx+N1.  <BR>
 * 
 * \param	kx  grid wavenumber along x
 * \return	i1  array index along x
 */
inline int Get_ix_FOUR(int kx, int N[])			
{ 
	return  (kx >= 0) ? kx : (kx + N[1]);   
}


/*! @brief	Get grid waveno ky given first local array index l2.
 * 
 *  If l2<=N2/2, k2=l2; else ky=l2-N2.
 * 
 * \param l2  second local index of an array
 * \return kx corresponding to l1
 */
inline int Get_ky3D_FOUR(int l2, int N[])   { return  (l2 <= N[2]/2) ? l2 : (l2-N[2]); } 

/*! @brief	Get local array index ly given grid waveno ky.
 * 
 *  If ky>=0, l2=ky; else l2=ky+N2.
 * 
 * \param	ky  grid wavenumber along y
 * \return	l2  local array index along y
 */
inline int Get_ly3D_FOUR(int ky, int N[]) { return  (ky >= 0) ? ky : (ky + N[2]); } 


/**********************************************************************************************

				 If wavenos computed using actual wavenumber:  
						globalvar_waveno_switch = 0;  
						Ki = Kfactor[i]*grid[i]
				 
				 If wavenos computed using grid wavenumber:  Ki = grid[i]
						globalvar_waveno_switch = 1;  
						Ki = grid[i]

***********************************************************************************************/


///  WAVENOACTUAL: \f$ K = \sqrt{K_x^2 + K_y^2 + K_z^2} \f$
inline DP Kmagnitude_FOUR(int l1, int l2, int l3, int N[], DP kfactor[])
{
	if	(globalvar_waveno_switch == 0)
		return sqrt(  pow2(Get_kx_FOUR(l1,N) * kfactor[1]) 
				+ pow2(Get_ky3D_FOUR(l2,N) * kfactor[2]) + pow2(l3* kfactor[3]) );
	
	else if (globalvar_waveno_switch == 1)
		return sqrt(  pow2(Get_kx_FOUR(l1,N)) + pow2(Get_ky3D_FOUR(l2,N)) + pow2(l3) ); 
	
	else
		return 0;		// for -Wall
}

/// WAVENOACTUAL -- Radius of the smallest sphere that contains the wavenumber K box Ni's. <BR>
/// The range of ki = [-Ni/2+1 : Ni/2].
inline int Min_radius_outside_aliased_FOUR(int N[], DP kfactor[]) 
{
	if	(globalvar_waveno_switch == 0)
		return (int) ceil(sqrt(  pow2((N[1]/2)*kfactor[1]) 
				+ pow2(N[2]/2 * kfactor[2]) + pow2(N[3]/2* kfactor[3]) ));
	
	else if (globalvar_waveno_switch == 1)
		return (int) ceil(sqrt(  pow2(N[1]/2) + pow2(N[2]/2) + pow2(N[3]/2) ));
	
	else
		return 0;		// for -Wall
		
}


/// WAVENOACTUAL -- Radius of the largest sphere that fits inside the wavenumber K box Ni's. <BR>
/// The range of ki = [-Ni/2+1 : Ni/2].
inline int Max_radius_inside_aliased_FOUR(int N[], DP kfactor[]) 
{
	int ans = 1;
	
	if	(globalvar_waveno_switch == 0)
	{
		DP kmag = min( (N[1]/2)*kfactor[1], (N[2]/2)*kfactor[2] );  
		
		if (N[3] > 2)
			kmag = min(kmag, (N[3]/2)*kfactor[3]);
			
		ans = ((int) kmag);	
	}
	
	else if (globalvar_waveno_switch == 1)
	{
		ans = min(N[1]/2, N[2]/2);
		
		if (N[3] > 2)
			ans = min(ans, (N[3]/2));
	}
	
	return ans;
	
}


/// WAVENOACTUAL -- Radius of the smallest sphere that contains the dealiased wavenumber  
///					K box Ni's. <BR>
/// The range of dealiased ki =  [-Ni/3+1 : Ni/3].
inline int Min_radius_outside_dealiased_FOUR(int N[], DP kfactor[]) 
{
	if	(globalvar_waveno_switch == 0)
		return (int) ceil(sqrt(  pow2((N[1]/3)*kfactor[1]) 
				+ pow2(N[2]/3 * kfactor[2]) + pow2(N[3]/3* kfactor[3]) ));
	
	else if (globalvar_waveno_switch == 1)
		return (int) ceil(sqrt(  pow2(N[1]/3) + pow2(N[2]/3) + pow2(N[3]/3) ));
	
	else
		return 0;		// for -Wall
}

/// WAVENOACTUAL -- Radius of the largest sphere that fits inside the dealised wavenumber  
///					K box Ni's. <BR>
/// The range of dealiased ki =  [-Ni/3+1 : Ni/3].
inline int Max_radius_inside_dealiased_FOUR(int N[], DP kfactor[]) 
{
	int ans = 1;
	
	if	(globalvar_waveno_switch == 0)
	{
		DP kmag = min( (N[1]/3)*kfactor[1], (N[2]/3)*kfactor[2] );  
		
		if (N[3] > 2)
			kmag = min(kmag, (N[3]/3)*kfactor[3]);
			
		return ((int) kmag);	
	}
	
	else if (globalvar_waveno_switch == 1)
	{
		ans = min(N[1]/3, N[2]/3);
		
		if (N[3] > 2)
			ans = min(ans, (N[3]/3));
	}
	
	return  ans;
		
}
	


/*! \brief WAVENOACTUAL -- Returns the approximate number of modes in a wavenumber K semishell 
 *							(upper hemisphere) of radius "radius".
 * 
 *  We divide the area in Fourier space by volume of unit lattice \f$ \Pi_i f_i \f$, 
 *	where \f$ f_i \f$ is the factor[i].
 *  In 1D, The area in NOT divided by kfactor.
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber 
 * or gridwavenumber depending on the switch
 * WAVENOACTUAL or WAVENOGRID.   Actual wavenumber is preferred.
 *
 * \param  radius
 * \return The number of modes in a hemispheric shell of radius. 
 *			In 3D, it is hemisphere sphere with (kz>=0).
 */
inline DP Approx_number_modes_in_shell_FOUR(int N[], int radius, DP kfactor[])
{
	DP ans = 0.0;
	
	if ((N[2] > 1) && (N[3] > 2))
	{
		if	(globalvar_waveno_switch == 0)
			ans = (2*M_PI*radius*radius)/(kfactor[1]*kfactor[2]*kfactor[3]);	
	
		else if (globalvar_waveno_switch == 1)
			ans = (2*M_PI*radius*radius);
	}
	else if (N[2] == 1)
	{
		if	(globalvar_waveno_switch == 0)
			ans = (M_PI*radius)/(kfactor[1]*kfactor[3]);
		
		else if	(globalvar_waveno_switch == 1)
			ans = (M_PI*radius);
	}
	else if (N[3] <= 2)
	{
		if	(globalvar_waveno_switch == 0)
			ans = (M_PI*radius)/(kfactor[1]*kfactor[2]);
		
		else if	(globalvar_waveno_switch == 1)
			ans = (M_PI*radius);
	}
	
	return ans;
}

//*********************************************************************************************

inline int Min_radius_outside_FOUR(string alias_switch, int N[], DP kfactor[])
{
	if (alias_switch == "ALIAS")
		return Min_radius_outside_aliased_FOUR(N, kfactor);
		
	else
		return Min_radius_outside_dealiased_FOUR(N, kfactor);	

} 


inline int Max_radius_inside_FOUR(string alias_switch, int N[], DP kfactor[]) 
{
	if (alias_switch == "ALIAS")
		return Max_radius_inside_aliased_FOUR(N, kfactor);
		
	else
		return Max_radius_inside_dealiased_FOUR(N, kfactor);
}


//*********************************************************************************************



/*! \brief Returns multiplication factor for computing enregy spectrum etc. 
 * 
 *  \f$ E(k) = |\vect{u}(\vect{k})|^2 /2 \f$, but in simulation we
 *  double the energy of most of the modes because complex conjugates modes with 
 * -ky are not stored in the simulation.  
 * The modes on the x-axis are not doubled because their c.c. are already counted.
 * 
 * \param  l1, l2
 * \return Multiplication factor to \f$ |\vect{u}(\vect{k})|^2  \f$ for computing 
 *		energy spectrum etc. factor = 1 implies that the modal energy is already doubled.
 */
 
inline DP Multiplicity_factor_FOUR(int l1, int l2, int l3, int N[])
{
	DP factor;

	int kx = Get_kx_FOUR(l1, N);
	int ky = Get_ky3D_FOUR(l2, N);
	
	if (l3 > 0)
		factor = 1.0;
		
	else
		factor = 0.5;
		
	if ( (kx == N[1]/2) || ((ky == N[2]/2)&& (N[2] > 1)) )
		return 2*factor;		// for kx = -N[1]/2 or ky = -N[2]/2 ; 
								// Ignoring corner which would have factor 4
	else
		return factor;
}



/// Modal energy -- \f$ E(k) = |\vect{u}(\vect{k})|^2 /2 \f$
inline DP Modal_energy_FOUR(Array<complx,3> A, int l1, int l2, int l3)
{
	return pow2(abs(A(l1,l2,l3)))/2;	
}



//*********************************************************************************************	

			
/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  \f$ K_{||} = K_1 \f$			
inline DP AnisKpll_FOUR(int l1, int l2, int l3, int N[], DP kfactor[])
{	
	if (globalvar_anisotropy_switch == 1)
		return (Get_kx_FOUR(l1,N) * kfactor[1]); 
	
	else if (globalvar_anisotropy_switch == 2)
		return (Get_ky3D_FOUR(l2,N) * kfactor[2]);
		
	else if (globalvar_anisotropy_switch == 3)
		return (l3 * kfactor[3]);
	
	else
		return 0;		// for -Wall
		
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  \f$ K_\perp =\sqrt{K_2^2 + K_3^2} \f$			
inline DP AnisKperp_FOUR(int l1, int l2, int l3, int N[], DP kfactor[])
{
	if (globalvar_anisotropy_switch == 1)
		return sqrt( pow2(Get_ky3D_FOUR(l2,N) * kfactor[2]) + pow2(l3* kfactor[3]) ); 
	
	else if (globalvar_anisotropy_switch == 2)
		return sqrt( pow2(Get_kx_FOUR(l1,N) * kfactor[1]) + pow2(l3* kfactor[3]) );
		
	else if (globalvar_anisotropy_switch == 3)
		return sqrt( pow2(Get_kx_FOUR(l1,N) * kfactor[1]) 
						+ pow2(Get_ky3D_FOUR(l2,N) * kfactor[2]) );
	else
		return 0;		// for -Wall
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  horizontal direction 1, \f$ K_{h1} = K_2 \f$										
inline DP AnisKh1_FOUR(int l1, int l2, int l3, int N[], DP kfactor[])
{	
	if (globalvar_anisotropy_switch == 1)
		return (Get_ky3D_FOUR(l2,N) * kfactor[2]);  
	
	else if (globalvar_anisotropy_switch == 2)
		return (l3 * kfactor[3]);
		
	else if (globalvar_anisotropy_switch == 3)
		return (Get_kx_FOUR(l1,N) * kfactor[1]);
	
	else
		return 0;		// for -Wall
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  horizontal direction 2, \f$ K_{h2} = K_3 \f$				
inline DP AnisKh2_FOUR(int l1, int l2, int l3, int N[], DP kfactor[])
{	
	if (globalvar_anisotropy_switch == 1)
		return (l3 * kfactor[3]);  
	
	else if (globalvar_anisotropy_switch == 2)
		return (Get_kx_FOUR(l1,N)  * kfactor[1]);
		
	else if (globalvar_anisotropy_switch == 3)
		return (Get_ky3D_FOUR(l2,N) * kfactor[2]); 
	
	else
		return 0;		// for -Wall
}
			
			
/// Cylindrical: Anis_min_Kpll
inline DP Anis_min_Kpll_FOUR(string alias_switch, int N[], DP kfactor[]) 
{ 
	/*
	if (alias_switch == "ALIAS")
	{
		if (globalvar_anisotropy_switch == 1)
			return ((-N[1]/2+1) * kfactor[1]);
		
		else if (globalvar_anisotropy_switch == 2)
			return ((-N[2]/2+1) * kfactor[2]);
			
		else if (globalvar_anisotropy_switch == 3)
			return 0.0;
	}
		 
		
	else
	{
		if (globalvar_anisotropy_switch == 1)
			return ((-N[1]/3+1) * kfactor[1]);
		
		else if (globalvar_anisotropy_switch == 2)
			return ((-N[2]/3+1) * kfactor[2]);
			
		else if (globalvar_anisotropy_switch == 3)
			return 0.0;
	}
	 */
		
	return 0.0;
}
				
/// Cylindrical: Anis_max_Kpll
inline DP Anis_max_Kpll_FOUR(string alias_switch, int N[], DP kfactor[]) 
{ 
	DP maxkpll = 0;
	
	if (alias_switch == "ALIAS")
	{	
		if (globalvar_anisotropy_switch == 1)
			maxkpll = ((N[1]/2) * kfactor[1]); 
		
		else if (globalvar_anisotropy_switch == 2)
			maxkpll = ((N[2]/2) * kfactor[2]); 
			
		else if (globalvar_anisotropy_switch == 3)
			maxkpll = ((N[3]/2) * kfactor[3]);
	}
		
	else
	{
		if (globalvar_anisotropy_switch == 1)
			maxkpll = ((N[1]/3) * kfactor[1]);	
		
		else if (globalvar_anisotropy_switch == 2)
			maxkpll = ((N[2]/3) * kfactor[2]);	 
			
		else if (globalvar_anisotropy_switch == 3)
			maxkpll = ((N[3]/3) * kfactor[3]);
		
	}
	
	return maxkpll;
}

/// 3D Cylindrical: Anis_max_Krho_radius_inside the wavenumber box
inline int Anis_max_Krho_radius_inside_FOUR(string alias_switch, int N[], DP kfactor[]) 	
{
	DP kmag = 0.0;
	
	if (alias_switch == "ALIAS")
	{
		if (globalvar_anisotropy_switch == 1)
			kmag = min( (N[2]/2)*kfactor[2], (N[3]/2)*kfactor[3] );
		
		else if (globalvar_anisotropy_switch == 2)
			kmag = min( (N[1]/2)*kfactor[1], (N[3]/2)*kfactor[3] );  
		
		else if (globalvar_anisotropy_switch == 3)
			kmag = min( (N[1]/2)*kfactor[1], (N[2]/2)*kfactor[2] ); 
	}
	
	else
	{
		if (globalvar_anisotropy_switch == 1)
			kmag = min( (N[2]/3)*kfactor[2], (N[3]/3)*kfactor[3] ); 
		
		else if (globalvar_anisotropy_switch == 2)
			kmag = min( (N[1]/3)*kfactor[1], (N[3]/3)*kfactor[3] ); 
			
		else if (globalvar_anisotropy_switch == 3)
			kmag = min( (N[1]/3)*kfactor[1], (N[2]/3)*kfactor[2] );	
	}
	
		return ((int) kmag);
}

// Max polar angle
inline DP Get_max_polar_angle_FOUR() 
{ 
	
	return M_PI/2;
	
}			

//*********************************************************************************************	


/*! \brief Returns the angle K vector makes with the anisotropic axis 
 * 
 * The range of angle is \f$ [0:\pi] \f$.
 *
 * \param  l1, l2, l3 (3D)
 * \return \f$ \tan^{-1}(K_{\perp}/K_{||}) \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	
inline DP AnisKvect_polar_angle_FOUR(int l1, int l2, int l3, int N[], DP kfactor[])
{
	DP kkpll, kkperp;
	
	kkpll = AnisKpll_FOUR(l1, l2, l3, N, kfactor);
	kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
	
	return Get_polar_angle(kkperp, kkpll);
}


inline DP AnisKvect_azimuthal_angle_2Din3Dgrid_FOUR
(
 int l1, int l2, int l3, 
 int N[], 
 DP kfactor[]
)
{
	// l3 = 0
	DP kkh1 = 1;
	DP kkh2 = 1;
	
	if (N[3]  <= 2)
	{
		if (globalvar_anisotropy_switch == 1)
		{ 
			kkh1 = Get_kx_FOUR(l1,N) * kfactor[1];
			kkh2 = Get_ky3D_FOUR(l2,N) * kfactor[2];
		}
		
		else if (globalvar_anisotropy_switch == 2)	
		{	
			kkh1 = Get_ky3D_FOUR(l2,N) * kfactor[2];
			kkh2 = Get_kx_FOUR(l1,N) * kfactor[1];
		}
	}
	
	return Get_azimuthal_angle(kkh1, kkh2);
}


/*! \brief 3D: Returns the azimutal angle.
 * 
 * The range of angle is \f$ [0:\pi] \f$.
 *
 * \param  l1, l2, l3 (3D)
 * \return \f$ \tan^{-1}(Ky}/Kx \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	
inline DP AnisKvect_azimuthal_angle_FOUR(int l1, int l2, int l3, int N[], DP kfactor[])
{
	
	DP kkh1 = AnisKh1_FOUR(l1, l2, l3, N, kfactor);
	DP kkh2 = AnisKh2_FOUR(l1, l2, l3, N, kfactor);
	
	return Get_azimuthal_angle(kkh1, kkh2);
}			


#endif

//********************************** END of four_inline.h *************************************



