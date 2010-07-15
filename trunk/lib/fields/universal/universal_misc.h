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



/*! \file  universal_misc.h
 * 
 * @brief  Universal functions based on basis fns
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug satisfy reality condition for SCFT case: shouldn't we set the imag part of some 
 *		parts to zero.
 */


#ifndef _UNIVERAL_MISC_H
#define _UNIVERAL_MISC_H

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif

#include <mpi.h>
#include <fftw3-mpi.h>


#include "universal_inline.h"

#include "../../four/four_basic.h"
#include "../../scft/scft_basic.h"

#include "../array_fn/array_basic.h"


using namespace blitz;

//*********************************************************************************************

																																																																																																																																																																																
/*! @brief 3D:  Ascertain reality condition for array A.
  *
  * Reality condition \f$ A(-\vec{k}) = A^{*}(\vec{k}) \f$ (FOUR). <BR>
  * Reality condition \f$ A(k_{||}, \vec{-k}_{\perp}) = A^*(k_{||}, \vec{k}_{\perp}) \f$. 
  *			(SCFT) <BR> 
  * To satisfy the reality condition, we save only kz>=0 modes. 
  *	In addition we need to do the following <BR>
  * In FOUR basis, \f$ A(-kx,-ky,0) = A^{*}(kx,ky,0) \f$ and  \f$ A(-kx,-ky,N_3/2) 
  *						= A^{*}(kx,ky,N_3/2) \f$ for ky>0. <BR>
  * In SCFT basis, \f$ A(kx,-ky,0) = A^{*}(kx,ky,0) \f$ and  \f$ A(kx,-ky,N_3/2) 
  *						= A^{*}(kx,ky,N_3/2) \f$ for ky>0. <BR>
  *
  * @param  basis_type		FOUR or SCFT
  * @param	 A()			Complex 2D array
  * @param  N[]			size of the array A.	
  *
  * @return The modified A satisfies reality condition.
  */																		
void Satisfy_reality_array(string basis_type,  int N[], Array<complx,3> A);

//*********************************************************************************************


/*! @brief 3D:  Dealiase array A using 2/3 rule (FOUR).
  *
  *  FOUR basis: |kx| >= N1/3 is saved in ix = N1/3:2*N1/3+1. So  A(N1/3:N1/2,:,:) =0. 
  *			Note that kx_max=N1/2.
  *			|ky| >= N2/3 is saved in iy = Ny/3:2*N2/3+1. So  A(:,N2/3:N2/2,:) =0. 
  *			Note that kx_max=N1/2.
  *			Along z, A(:,:,N3/3:N3/2) =0.
  *
  *  SCFT basis: set A(2*N1/3:N1-1,:,:) =0.
  *			|ky| >= N2/3 is saved in iy = Ny/3:2*N2/3+1. So  A(:,N2/3:N2/2,:) =0. 
  *			Note that kx_max=N1/2.
  *			Along z, A(:,:,N3/3:N3/2) =0.
  *
  *  @param  basis_type		FOUR or SCFT
  *  @param	 A()			Complex 2D array
  *  @param  N[]			size of the array A: kx=0:N1/2; ky=-N2/2:N2/2.	
  *
  *  @return Dealiased A.
  */
void Dealias_array(string basis_type, int N[], Array<complx,3> A);


//*********************************************************************************************

/*! @brief 3D:  Test aliasing needed or not in array using 2/3 rule (FOUR). 
 *
 *  @param  basis_type		FOUR
 *  @param	 A()			Complex 2D array
 *  @param  N[]			size of the array A: kx=-N1/2+1:N1/2; ky = 0:N2/2.	
 *  @param  outer_radius	A(kx) is nonzero for 0<=\f$ |\vec{K}| \f$ <= ((int) outer_radius) + 1.
 *
 *  @return FOUR basis: 1 if ((outer_radius/f1 > N[1]/3) || (outer_radius/f2 > N[2]/3) 
 *								|| (outer_radius/f3 > N[3]/3)); 0 otherwise. 
 *  @return SCFT basis: 1 if ((outer_radius/f1 > 2*N[1]/3) || (outer_radius/f2 > N[2]/3) 
 *								|| (outer_radius/f3 > N[3]/3)); 0 otherwise. 
 */
int Is_alias_array
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	DP outer_radius, 
	DP kfactor[]
);



#endif


//******************************** End of universal_misc.h  ***********************************


