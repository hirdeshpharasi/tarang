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

/*! \files	planar_struct_fn.h
 * 
 * @brief Structure function calculations
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */ 


#ifndef _PLANAR_STRUCT_FN_H
#define _PLANAR_STRUCT_FN_H

#include "array_basic.h"

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif



using namespace blitz;


//*********************************************************************************************


/*! @brief Compute \f$ (\Delta v)^q \f$ for a pair of points in a plane
 * 
 *	@param	TinyVector<int,3> j1, j2 The two point in the plane; the plane coordinates are
 *				the same.
 *	@param	TinyVector<DP,3> V1, V2	Velocity at the two points.
 *	@param	q_min, q_max: min and max value of q (structure fn index).
 *	@param	kpll  The index of the plane.
 *	@param	xfactor[i]: min_dx[i].
 *			We ignore \f$ 2 \pi \f$ which is the overall factor.
 *
 *	@return \f$ St(kpll,r,q,0) += |\Delta v_{||}|^q \f$  with r=ceil(|r2-r1|).
 *  @return \f$ St(kpll,r,q,0) += |\Delta v_{\perp}|^q \f$.
 *	@return	count(kpll,r) += 1.
 *
 */
void Compute_dSt_planar
(
	string basis_type,
	int N[],
	TinyVector<int,3> j1, TinyVector<int,3> j2,
	TinyVector<DP,3> V1, TinyVector<DP,3> V2,
	int	q_min, int q_max,
	int kpll,
	Array<DP,4> St_planar, 
	Array<int,2> count,
	DP xfactor[]
);

//
//
/*! @brief Compute \f$ (\Delta v)^q \f$ for a pair of points in a plane
 * 
 *	@param	TinyVector<int,3> j1, j2 The two point in the plane; the plane coordinates are
 *				the same.
 *	@param	DP F1, F2	Value of the scalar functions at the two points.
 *	@param	q_min, q_max: min and max value of q (structure fn index).
 *	@param	kpll  The index of the plane.
 *	@param	xfactor[i]: min_dx[i].
 *			We ignore \f$ 2 \pi \f$ which is the overall factor.
 *
 *	@return \f$ St(kpll,r,q) += |\Delta F|^q \f$  with r=ceil(|r2-r1|).
 *	@return	count(kpll,r) += 1.
 *
 */
void Compute_dSt_planar
(
	string basis_type,
	int N[],
	TinyVector<int,3> j1, TinyVector<int,3> j2,
	DP F1, DP F2,
	int	q_min, int q_max,
	int kpll,
	Array<DP,3> St_planar, 
	Array<int,2> count,
	DP xfactor[]
);

//******************************************************************************

/** @brief 3D: Compute the planar structure function of vector field 
 *					\f$ \vec{v}(\vec{x}) \f$.
 *
 * @param	Vx, Vy, Vz  Velocity field; Real values stored in complex arrays.
 * @param	N[]  The size of the arrays A.
 * @param	q_min
 * @param	q_max
 * @param	xfactor[i] = min_dx[i]; 
 *			we ignore \f$ 2 \pi \f$ which is the overall factor.
 *
 * @return	\f$ St(kpll, r, q, 0) = |\Delta v_{||}|^q \f$  with r=ceil(|r2-r1|) in the plane.
 * @return	\f$ St(kpll, r, q, 1) = |\Delta v_{\perp}|^q \f$  
 *
 * @note	Loop over planes, and all points P & Q (Q>=P) in the plane; 
 *			Find r between the two points, 
 *			and compute St(r,q,0-1) by summing over all r's.  
 * @note	When P=Q, the real and imag parts are two different points in R.
 * @note	The final St(r,q,0-1) is divided by total no of points with sep r.
 */

void Compute_planar_structure_function
(
	string basis_type,
	int N[], 
	Array<complx,3> Vx, Array<complx,3> Vy, Array<complx,3> Vz,
	int q_min, int q_max, 
	Array<DP,4> St_planar, 
	DP xfactor[]
);

//
//

/// Same as above but for scalar field F.
void Compute_planar_structure_function
(
	string basis_type,
	int N[], 
	Array<complx,3> F,
	int q_min, int q_max, 
	Array<DP,3> St_planar, 
	DP xfactor[]
);

#endif

//*********************************   End of struct_fn.h  *************************************



