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


/*! \file  Sincos_horizontal.cc
 * 
 * @brief sinccos basis function for horizontal functions for RB convection type problems
 *
 * @note  Vx = sin(m pi x) cos(n k0 y)  [Buoyancy dirn] <BR>
 *		  Vy = cos(m pi x) sin(n k0 y)  <BR>
 *
 * @note  Vx = real with Vx(kx,-ky,kz) = Vx(kx,ky,kz) [Buoyancy dirn] <BR>
 *		  Vy = imag with Vx(kx,-ky,kz) = -Vy(kx,ky,kz)  <BR>
 *		  Vz = imag with Vz(kx,-ky,kz) = Vz(kx,ky,kz)  <BR>
 *			Naturally Vy(kx,0.kz) = Vy(kz,N[2]/2,kz) = 0. <BR>
 *
 *	We set imag(Vx) = 0 and determine Vy using incompressibility conditions. <BR>
 *  For ky (or n) = 0 =>  Vy = 0.  Incompressibility => Vx = 0.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************

void IncVF::free_slip_verticalwall()
{

	imag(*V1) = 0.0;
	for (int i2=1; i2<N[2]/2; i2++)
		(*V1)(Range::all(),N[2]-i2,Range::all()) = (*V1)(Range::all(),i2,Range::all());
	
	real(*V2) = 0.0;
	for (int i2=1; i2<N[2]/2; i2++)
		(*V2)(Range::all(),N[2]-i2,Range::all()) = -(*V2)(Range::all(),i2,Range::all());
	
	imag((*V2)(Range::all(),0,Range::all())) = 0.0;
	imag((*V2)(Range::all(),N[2]/2,Range::all())) = 0.0;
	
	real(*V3) = 0.0;
	for (int i2=1; i2<N[2]/2; i2++)
		(*V3)(Range::all(),N[2]-i2,Range::all()) = (*V3)(Range::all(),i2,Range::all());
	
	Satisfy_reality_condition_field();
	
}	

//*********************************************************************************************

void IncVF::free_slip_verticalwall(IncSF& T)
{
	free_slip_verticalwall();
	
	imag((*T.F)) = 0.0;
}


//*********************************************************************************************

void IncVF::free_slip_verticalwall(IncVF& W)
{
	free_slip_verticalwall();
	
	W.free_slip_verticalwall();
}


//*********************************************************************************************

void IncVF::free_slip_verticalwall(IncVF& W, IncSF& T)
{
	free_slip_verticalwall();
	
	W.free_slip_verticalwall();
	
	imag((*T.F)) = 0.0;
}


//**********************************   End of compute_dt.cc  **********************************





