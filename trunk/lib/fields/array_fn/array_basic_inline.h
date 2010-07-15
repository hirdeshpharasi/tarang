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

/*! \file array_basic_inline.h
 * 
 * @brief Some basic Array operations: Array_real_mult, Output_asreal, 
 *			Model_initial_shell_spectrum
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug  No know bug
 */ 


#ifndef _FIELD_ARRAY_BASIC_INLINE_H
#define _FIELD_ARRAY_BASIC_INLINE_H

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif


#include "../../basis_basicfn/basis_basicfn_inline.h"
#include "../../basis_basicfn/basis_basicfn.h"


using namespace blitz;

//******************************************************************************

/// return \f$ |\vec{r}_2 - \vec{r}_1| \f$ for 2D.
inline DP Get_abs_dr
(
	TinyVector<DP,2> r1, 
	TinyVector<DP,2> r2
) 
{
	return sqrt(pow2(r2(0)-r1(0)) + pow2(r2(1)-r1(1)));
}


/// return \f$ |\vec{r}_2 - \vec{r}_1| \f$ for 3D.
inline DP Get_abs_dr
(	
	TinyVector<DP,3> r1, 
	TinyVector<DP,3> r2, 
	TinyVector<DP,3> r3
) 
{
	return sqrt(pow2(r2(0)-r1(0)) + pow2(r2(1)-r1(1)) + pow2(r2(2)-r1(2)));
}


//******************************************************************************

/** @brief Maximum radius inside the real sphere.
 *	Along the periodic dirn, max(r)=N[i]/2. \\
 *	Along the nonperiodic dirn, max(r)=N[i]-1. 
 *
 *	@return  min of max(r) along the three axis.
 */
inline int Max_radius_inside_real_space(string basis_type, int N[], DP xfactor[])
{
	if (basis_type == "FOUR")
	{
		DP xmag = min( (N[1]/2)*xfactor[1], (N[2]/2)*xfactor[2] );  
		xmag = min(xmag, (N[3]/2)*xfactor[3]);
		
		return (int) ceil(xmag);	
	}
	
	else if (basis_type == "SCFT")
	{
		DP xmag = min( (N[1]-1)*xfactor[1], (N[2]/2)*xfactor[2] );  
		xmag = min(xmag, (N[3]/2)*xfactor[3]);
		
		return (int) ceil(xmag);	
	}
	
	else
		return 0; // for -Wall
}


//******************************************************************************

inline int Min_radius_outside_real_space(string basis_type, int N[], DP xfactor[])
{
	if (basis_type == "FOUR")
		return (int) ceil( sqrt( pow2((N[1]/2)*xfactor[1]) + pow2((N[2]/2)*xfactor[2]) 
						+ pow2((N[3]/2)*xfactor[3]) ) );
						
	else if (basis_type == "SCFT")
		return (int) ceil( sqrt( pow2((N[1]-1)*xfactor[1]) + pow2((N[2]/2)*xfactor[2]) 
						+ pow2((N[3]/2)*xfactor[3]) ) );	
						
	else
		return 0; // for -Wall																	
}



#endif

//*********************************  End of array_basic_inline.h ******************************



