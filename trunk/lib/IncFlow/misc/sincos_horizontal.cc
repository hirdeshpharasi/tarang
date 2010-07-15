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
 * @note  Vx = sin(m pi x) cos(n k0 y) [Buoyancy dirn] <bR>
 *		  Vy = cos(m pi x) sin(n k0 y)  <BR>
 *
 *	We set imag(Vx) = 0 and determine Vy using incompressibility conditions. <BR>
 *  For ky (or n) = 0 =>  Vy = 0.  Incompressibility => Vx = 0.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************

void IncVF::Sincos_horizontal()
{
	int kx, ky, kz;
	complx Vx, Vy;
	
	
	imag(*V1) = 0.0;
	(*V3) = 0.0;
	
	// kz = 0 plane
	for (int lx=0; lx<local_N1; lx++)										
		for (int ly=0; ly<N[2]; ly++)  
		{													
			kx = Get_kx(basis_type, lx, N);
			ky = Get_ky3D(basis_type, ly, N);
			kz = 0;
			
			if (abs(ky) > 0)
			{
				Vx = (*V1)(lx, ly, 0);
				Vy = DxVx(basis_type, kfactor, kx, Vx) / complex<DP>(0, -ky*kfactor[3]);
				
				(*V2)(lx, ly, 0) = Vy;
			}
			
			else
			{
				(*V1)(lx, ly, 0) = 0.0;  
				(*V2)(lx, ly, 0) = 0.0; 
			}
		}
}	

//*********************************************************************************************

void IncVF::Sincos_horizontal(IncSF& T)
{
	Sincos_horizontal();
	
	imag((*T.F)) = 0.0;
}


//*********************************************************************************************

void IncVF::Sincos_horizontal(IncVF& W)
{
	Sincos_horizontal();
	
	W.Sincos_horizontal();
}


//*********************************************************************************************

void IncVF::Sincos_horizontal(IncVF& W, IncSF& T)
{
	Sincos_horizontal();
	
	W.Sincos_horizontal();
	
	imag((*T.F)) = 0.0;
}


//**********************************   End of compute_dt.cc  **********************************




