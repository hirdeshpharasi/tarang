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

/*! \file sincosfour_inline.h 
 * 
 * @brief Inline functions to compute array indices \f$ \vec{i} \f$ given wavenumber  
 *		\f$ \vec{k} \f$ and viceversa.
 *		Also contains other useful functions like Kmagnitude, Max radius etc.
 *
 * Grid wavenumber \f$ \vec{k} \f$ is integer, that can be computed using the grid index 
 *		\f$ \vec{i} \f$. Actual wavenumber \f$ \vec{K} \f$ is computed from the grid wavenumber 
 *		\f$ \vec{k} \f$ using \f$ K_i = k_i * f_i \f$  where \f$ f_i \f$ is the kfactor[i]. 
 *
 * l1, l2, l3 = local array indices. <BR>
 * i1, i2, i3 = global or total array indices.
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber 
 *		or gridwavenumber depending on the switch WAVENOACTUAL or WAVENOGRID.
 * 
 * @author  M. K. Verma
 * @version 4.0 Parallel
 * @date	August 2008
 * @bug		No known bugs
 */

#ifndef _SCFT_INLINE_H
#define _SCFT_INLINE_H

#include <blitz/array.h>
#include <complex>
#include <cmath>

#include <mpi.h>
#include <fftw3-mpi.h>

#include "../basis_basicfn/basis_basicfn_inline.h"
#include "../basis_basicfn/basis_basicfn.h"

using namespace blitz;



//global vars
extern int my_id;								/**<	External: My process id	*/
extern int numprocs;							/**<	External: No of processors	*/
extern const int master_id;						/**<	External: Id of master proc	*/
extern ptrdiff_t local_N1;						/**<	External: local_N1 = N1/numprocs in the currentproc */
extern ptrdiff_t local_N1_start;				/**<	External: start of i1 in the currentproc */
extern ptrdiff_t local_N2;						/**<	External: local_N2 = N2/numprocs in the currentproc */
extern ptrdiff_t local_N2_start;				/**<	External: start of i2 in the currentproc */

extern MPI_Status status;

//*********************************************************************************************

/*! @brief	Get grid waveno kx given first local array index l1.
 * 
 *	i1= local_N1_start + l1.
 * 
 * \param l1  first local index of an array
 * \return kx corresponding to l1
 */
inline int Get_kx_SCFT(int l1, int N[]) { return  (local_N1_start + l1); }


/*! @brief	Get local array index lx given grid waveno kx.
 * 
 *	l1= kx - local_N1_start.
 * 
 * \param	kx  grid wavenumber along x
 * \return	l1  local array index along x
 */
inline int Get_lx_SCFT(int kx, int N[])  { return  (kx-local_N1_start); } 



/*! @brief	Get grid waveno ky given first local array index l2.
 * 
 *  If l2<=N2/2, k2=l2; else ky=l2-N2.
 * 
 * \param l2  second local index of an array
 * \return kx corresponding to l1
 */
inline int Get_ky3D_SCFT(int l2, int N[]) { return  (l2 <= N[2]/2) ? l2 : (l2-N[2]); }


/*! @brief	Get local array index ly given grid waveno ky.
 * 
 *  If ky>=0, l2=ky; else l2=ky+N2.
 * 
 * \param	ky  grid wavenumber along y
 * \return	l2  local array index along y
 */
inline int Get_ly3D_SCFT(int ky, int N[]) { return  (ky >= 0) ? ky : (ky + N[2]); }


/**********************************************************************************************

		  If wavenos computed using actual wavenumber:  Ki = Kfactor[i]*grid[i]

***********************************************************************************************/


///  WAVENOACTUAL: \f$ K = \sqrt{K_x^2 + K_y^2 + K_z^2} \f$
inline DP Kmagnitude_SCFT(int l1, int l2, int l3, int N[], DP kfactor[])
{ 
	if	(globalvar_waveno_switch == 0)
		return sqrt( pow2(Get_kx_SCFT(l1,N)*kfactor[1]) 
					+ pow2(Get_ky3D_SCFT(l2,N)*kfactor[2]) + pow2(l3*kfactor[3]) ); 
	
	else if (globalvar_waveno_switch == 1)
		return sqrt( pow2(Get_kx_SCFT(l1,N)) + pow2(Get_ky3D_SCFT(l2,N)) + pow2(l3) );
	
	else
		return 0;		// for -Wall
}


/// WAVENOACTUAL -- Radius of the smallest sphere that contains the wavenumber K box Ni's. <BR>
/// The range of kx=[0:N1-1];  ki = [-Ni/2+1 : Ni/2] along perp directions.
inline int Min_radius_outside_aliased_SCFT(int N[], DP kfactor[]) 
{
	if	(globalvar_waveno_switch == 0)
		return (int) ceil(sqrt(  pow2((N[1]-1)*kfactor[1]) 
					  + pow2(N[2]/2 * kfactor[2]) + pow2(N[3]/2* kfactor[3]) ));
	
	else if (globalvar_waveno_switch == 1)
		return (int) ceil(sqrt(  pow2(N[1]-1) + pow2(N[2]/2) + pow2(N[3]/2) ));
	
	else
		return 0;		// for -Wall
}


/// WAVENOACTUAL -- Radius of the largest sphere that fits inside the wavenumber K box Ni's. <BR>
/// The range of kx=[0:N1-1];  ki = [-Ni/2+1 : Ni/2] along perp directions.
inline int Max_radius_inside_aliased_SCFT(int N[], DP kfactor[]) 
{
	int ans = 1;
	
	if	(globalvar_waveno_switch == 0)
	{
		DP kmag = min( (N[1]-1)*kfactor[1], (N[2]/2)*kfactor[2] );  
		
		if (N[3] > 2)
			kmag = min(kmag, (N[3]/2)*kfactor[3]);
		
		ans = ((int) kmag);
	}
	
	else if (globalvar_waveno_switch == 1)
	{
		ans = min(N[1]-1, N[2]/2);
		
		if (N[3] > 2)
			ans = min(ans, (N[3]/2));
	}
	
	return ans;
}		


/// WAVENOACTUAL -- Radius of the smallest sphere that contains the dealiased wavenumber  
///					K box Ni's. <BR>
/// The range of dealiased  kx=[0:2*N1/3];  ki = [-Ni/3+1 : Ni/3] along perp directions.
inline int Min_radius_outside_dealiased_SCFT(int N[], DP kfactor[]) 
{
	if	(globalvar_waveno_switch == 0)
		return (int) ceil(sqrt(  pow2(2*N[1]/3*kfactor[1]) 
				+ pow2(N[2]/3 * kfactor[2]) + pow2(N[3]/3* kfactor[3]) ));
	
	else if (globalvar_waveno_switch == 1)
		return  (int) ceil(sqrt(  pow2(2*N[1]/3) + pow2(N[2]/3) + pow2(N[3]/3) ));
	
	else
		return 0;		// for -Wall
}


/// WAVENOACTUAL -- Radius of the largest sphere that fits inside the dealised wavenumber  
///					K box Ni's. <BR>
/// The range of dealiased  kx=[0:2*N1/3];  ki = [-Ni/3+1 : Ni/3] along perp directions.
inline int Max_radius_inside_dealiased_SCFT(int N[], DP kfactor[]) 
{
	int ans = 1;
	
	if	(globalvar_waveno_switch == 0)
	{
		DP kmag = min( 2*N[1]/3*kfactor[1], (N[2]/3)*kfactor[2] );  
		
		if (N[3] > 2)
			kmag = min(kmag, (N[3]/3)*kfactor[3]);
		
		ans = ((int) kmag);
	}
	
	else if (globalvar_waveno_switch == 1)
	{
		ans = min(2*N[1]/3, N[2]/3);
		
		if (N[3] > 2)
			ans = min(ans, (N[3]/3));
	}
	
	return  ans;
}	


/*! \brief WAVENOACTUAL -- Returns the approximate number of modes in a wavenumber 
 *			K shell of radius "radius".
 * 
 *  We divide the area in Fourier space by volume of unit lattice \f$ \Pi_i f_i \f$, 
 *			where \f$ f_i \f$ is the factor[i]. In 1D, The area in NOT divided by kfactor.
 *
 * \param  radius
 * \return The number of modes in a shell of radius. In 2D, it is quarter circle (kx, ky>= 0). 
 *			In 3D, it is quarter sphere with (kx,kz>=0).
 */
inline DP Approx_number_modes_in_shell_SCFT(int N[], int radius, DP kfactor[])
{
	DP ans = 0.0;
	
	if ((N[2] > 1) && (N[3] > 2))
	{
		if (globalvar_waveno_switch == 0)
			ans = (M_PI*radius*radius)/(kfactor[1]*kfactor[2]*kfactor[3]);	
		
		else if (globalvar_waveno_switch == 1)
			ans = (M_PI*radius*radius);
	}
	else if (N[2] == 1)
	{
		if	(globalvar_waveno_switch == 0)
			ans = (M_PI*radius/2)/(kfactor[1]*kfactor[3]);
		
		else if	(globalvar_waveno_switch == 1)
			ans = (M_PI*radius/2);
	}
	else if (N[3] <= 2)
	{
		if	(globalvar_waveno_switch == 0)
			ans = (M_PI*radius/2)/(kfactor[1]*kfactor[2]);
		
		else if	(globalvar_waveno_switch == 1)
			ans = (M_PI*radius/2);
	}
	
	return ans;
}


//*********************************************************************************************

inline int Min_radius_outside_SCFT(string alias_switch, int N[], DP kfactor[])
{
	if (alias_switch == "ALIAS")
		return Min_radius_outside_aliased_SCFT(N, kfactor);
		
	else
		return Min_radius_outside_dealiased_SCFT(N, kfactor);	

} 


inline int Max_radius_inside_SCFT(string alias_switch, int N[], DP kfactor[]) 
{
	if (alias_switch == "ALIAS")
		return Max_radius_inside_aliased_SCFT(N, kfactor);
		
	else
		return Max_radius_inside_dealiased_SCFT(N, kfactor);
}


//*********************************************************************************************



/*! \brief Returns multiplication factor for computing enregy spectrum etc. 
 *  
 * Modal energy  = \f$ E(k) =  |A(k_x, \vec{k}_{\perp})|^{2} \f$  if  \f$ (k_x > 0) \f$.  <BR>
 * Modal energy  = \f$ E(k) =  |A(0, \vec{k}_{\perp})|^{2}/2 \f$  if  \f$ (k_x = 0) \f$.  <BR>
 * In simulation we double the energy of most of the modes because complex conjugates modes with 
 * -ky are not stored in the simulation.  
 * The modes on the xy-plane are not doubled because their c.c. are already counted.
 * 
 * \param  l1, l2, l3
 * \return Multiplication factor for computing enregy spectrum etc.
 */
inline DP Multiplicity_factor_SCFT(int l1, int l2, int l3, int N[])
{
	DP factor;

	int m = Get_kx_SCFT(l1, N);
	int ky = Get_ky3D_SCFT(l2, N);
	
	if (l3 > 0)
		if ( m > 0 )
			factor  = 2.0;
		else				// m = 0 plane
			factor = 1.0;
			
	else					// kz = 0 plane
		if (m > 0)
			factor = 1.0;
		else				// kz = 0; m = 0 line
			factor = 0.5;
			
	if ((ky == N[2]/2) && (N[2] > 1))
		return 2*factor;	// for both ky = N[2]/2 and -N[2]/2
	else					
		return factor;
}



/**********************************************************************************************

		Modal energy

***********************************************************************************************/
/// Modal energy  = \f$ E(k) =  |A(k_x, \vec{k}_{\perp})|^{2} \f$  if  \f$ (k_x > 0) \f$.  <BR>
/// Modal energy  = \f$ E(k) =  |A(0, \vec{k}_{\perp})|^{2}/2 \f$  if  \f$ (k_x = 0) \f$.  <BR>
inline DP Modal_energy_SCFT(Array<complx,3> A, int l1, int l2, int l3)
{
	if (l1 == 0)
		return pow2(abs(A(l1,l2,l3)))/2;
	
	else
		return pow2(abs(A(l1,l2,l3)));	
}


//*********************************************************************************************

/// 3D == Anisotropic axis along x1: for anisotropic energy spectrum and 
///			energy transfer calculations,  \f$ K_{||} = K_1 \f$.		
inline DP AnisKpll_SCFT(int l1, int l2, int l3, int N[], DP kfactor[])
{	
	if (globalvar_anisotropy_switch == 1)
		return (Get_kx_SCFT(l1, N)*kfactor[1]); 
	
	else if (globalvar_anisotropy_switch == 2)
		return (Get_ky3D_SCFT(l2,N) * kfactor[2]);
	
	else if (globalvar_anisotropy_switch == 3)
		return (l3 * kfactor[3]);
		
	else
		return 0;		// for -Wall
}

/// 3D == Anisotropic axis along x1: for anisotropic energy spectrum and 
///			energy transfer calculations,  \f$ K_\perp =\sqrt{K_2^2 + K_3^2} \f$.			
inline DP AnisKperp_SCFT(int l1, int l2, int l3, int N[], DP kfactor[])
{
	if (globalvar_anisotropy_switch == 1)
		return sqrt( pow2(Get_ky3D_SCFT(l2,N) * kfactor[2]) + pow2(l3*kfactor[3]) ); 
	
	else if (globalvar_anisotropy_switch == 2)
		return sqrt( pow2(Get_kx_SCFT(l1, N)*kfactor[1]) + pow2(l3*kfactor[3]) );
		
	else if (globalvar_anisotropy_switch == 3)
		return sqrt( pow2(Get_kx_SCFT(l1, N) * kfactor[1]) 
					+ pow2(Get_ky3D_SCFT(l2,N) * kfactor[2]) );
			
	else
		return 0;		// for -Wall
}

/// 3D == Anisotropic axis along x1: for anisotropic energy spectrum and 
///			energy transfer calculations,  horizontal direction 1, \f$ K_{h1} = K_2 \f$.										
inline DP AnisKh1_SCFT(int l1, int l2, int l3, int N[], DP kfactor[])
{	
	if (globalvar_anisotropy_switch == 1)
		return (Get_ky3D_SCFT(l2,N) * kfactor[2]); 
	
	else if (globalvar_anisotropy_switch == 2)
		return (l3 * kfactor[3]);
		
	else if (globalvar_anisotropy_switch == 3)
		return (Get_kx_SCFT(l1, N) * kfactor[1]);
			
	else
		return 0;		// for -Wall
}

/// 3D == Anisotropic axis along x1: for anisotropic energy spectrum and 
///			energy transfer calculations,  horizontal direction 2, \f$ K_{h2} = K_3 \f$.				
inline DP AnisKh2_SCFT(int l1, int l2, int l3, int N[], DP kfactor[])
{	
	if (globalvar_anisotropy_switch == 1)
		return (l3 * kfactor[3]);  
	
	else if (globalvar_anisotropy_switch == 2)
		return (Get_kx_SCFT(l1, N) * kfactor[1]);
		
	else if (globalvar_anisotropy_switch == 3)
		return (Get_ky3D_SCFT(l2,N) * kfactor[2]);
			
	else
		return 0;		// for -Wall
}
			
/// Cylindrical: Anis_min_Kpll
inline DP Anis_min_Kpll_SCFT(string alias_switch, int N[], DP kfactor[]) 
{ 
	/*
	if (globalvar_anisotropy_switch == 1)
		return 0.0; 
	
	else if (globalvar_anisotropy_switch == 3)
		return 0.0;
		
	else if (globalvar_anisotropy_switch == 2)
	{
		if (alias_switch == "ALIAS")
			return ((-N[2]/2+1) * kfactor[2]);
		
		else  
			return ((-N[2]/3+1) * kfactor[2]);
	}		
	
	else
		return 0;		// for -Wall
	 */
	
	return 0.0;
}
				
/// Cylindrical: Anis_max_Kpll
inline DP Anis_max_Kpll_SCFT(string alias_switch, int N[], DP kfactor[]) 
{ 
	
	DP maxkpll = 0.0;
	
	if (alias_switch == "ALIAS")
	{
		if (globalvar_anisotropy_switch == 1)
			maxkpll = ((N[1]-1) * kfactor[1]); 
		
		else if (globalvar_anisotropy_switch == 2)
			maxkpll = ((N[2]/2) * kfactor[2]); 
			
		else if (globalvar_anisotropy_switch == 3)
			maxkpll = ((N[3]/2) * kfactor[3]); 
	}
		
	else
	{
		if (globalvar_anisotropy_switch == 1)
			maxkpll = ((2*N[1]/3) * kfactor[1]);
		
		else if (globalvar_anisotropy_switch == 2)
			maxkpll = ((N[2]/2) * kfactor[2]); 
			
		else if (globalvar_anisotropy_switch == 3)
			maxkpll = ((N[3]/3) * kfactor[3]);	
	}
	
	return maxkpll;
}

	
/// 3D Cylindrical: Anis_max_Krho_radius_inside the wavenumber box.
inline int Anis_max_Krho_radius_inside_SCFT(string alias_switch, int N[], DP kfactor[]) 			
{
	DP kmag = 0.0;
	
	if (alias_switch == "ALIAS")
	{
		if (globalvar_anisotropy_switch == 1)
			kmag = min( (N[2]/2)*kfactor[2], (N[3]/2)*kfactor[3] ); 
		
		else if (globalvar_anisotropy_switch == 2)
			kmag = min( (N[1]-1)*kfactor[1], (N[3]/2)*kfactor[3] );
			
		else if (globalvar_anisotropy_switch == 3)
			kmag = min( (N[1]-1)*kfactor[1], (N[2]/2)*kfactor[2] ); 
	}
		
	else
	{
		if (globalvar_anisotropy_switch == 1)
			kmag = min( (N[2]/3)*kfactor[2], (N[3]/3)*kfactor[3] );
		
		else if (globalvar_anisotropy_switch == 2)
			kmag = min( (2*N[1]/3)*kfactor[1], (N[3]/3)*kfactor[3] );
			
		else if (globalvar_anisotropy_switch == 3)
			kmag = min( (2*N[1]/3)*kfactor[1], (N[2]/3)*kfactor[2] ); 
	}
	
	return ((int) kmag);	
}

// Max polar angle
inline DP Get_max_polar_angle_SCFT() 
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
inline DP AnisKvect_polar_angle_SCFT(int l1, int l2, int l3, int N[], DP kfactor[])
{
	DP kkpll, kkperp;
	
	kkpll = AnisKpll_SCFT(l1, l2, l3, N, kfactor);
	kkperp = AnisKperp_SCFT(l1, l2, l3, N, kfactor);
	
	return Get_polar_angle(kkperp, kkpll);
}


inline DP AnisKvect_azimuthal_angle_2Din3Dgrid_SCFT
(
 int l1, int l2, int l3, 
 int N[], 
 DP kfactor[]
)
{
	
	// i3 = 0
	DP kkh1 = 1;
	DP kkh2 = 1;
	
	if (N[3]  <= 2)
	{
		if (globalvar_anisotropy_switch == 1)
		{ 
			kkh1 = Get_kx_SCFT(l1, N) * kfactor[1];
			kkh2 = Get_ky3D_SCFT(l2,N) * kfactor[2];
		}
		
		else if (globalvar_anisotropy_switch == 2)	
		{	
			kkh1 = Get_ky3D_SCFT(l2,N) * kfactor[2];
			kkh2 = Get_kx_SCFT(l1, N) * kfactor[1];
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
inline DP AnisKvect_azimuthal_angle_SCFT(int l1, int l2, int l3, int N[], DP kfactor[])
{
	
	DP kkh1 = AnisKh1_SCFT(l1, l2, l3, N, kfactor);
	DP kkh2 = AnisKh2_SCFT(l1, l2, l3, N, kfactor);
	
	return Get_azimuthal_angle(kkh1, kkh2);
}			
						

#endif

//=================================== End of inline functions ==========================//



