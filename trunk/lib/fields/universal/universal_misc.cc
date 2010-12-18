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


/*! \file  universal_misc.cc
 * 
 * @sa	universal_misc.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug No known bugs
 */
 
#include "universal_misc.h"


/**********************************************************************************************

	Work on kz=0 and N[3]/2 plane such that they satisfy reality condition
	
	// For kz =0 and kz = N[3]/2 planes
	// f(-kx, -ky, 0) = conj(f(kx, ky, 0))- do for kz=N[3]/2

***********************************************************************************************/ 

void Satisfy_reality_array(string basis_type,  int N[], Array<complx,3> A)
{
	static Array<complx,2> plane_lower(N[1], N[2]);
	static Array<complx,2> plane_upper(N[1], N[2]);
	
	int kx;
	int conj_nx;		// ix for conj location for nx
	int positive_ix, positive_iy;
	
	// For kz =0 and kz = N[3]/2 planes
	// f(-kx, -ky, 0) = conj(f(kx, ky, 0))- do for kz=N[3]/2
	if (basis_type == "FOUR")
	{
		Get_XY_plane(N, A, plane_lower, 0);
		Get_XY_plane(N, A, plane_upper, N[3]/2);
		
		for (int lx=0; lx<local_N1; lx++)
		{
			kx = Get_kx_FOUR(lx, N);
			conj_nx =  Get_ix_FOUR(-kx, N);				// for ix=0:N[1]-1 range
			
			for (int ly=N[2]/2+1; ly<N[2]; ly++)
			{
				positive_iy = -Get_ky3D_FOUR(ly, N);
				
				A(lx, ly, 0) = conj( plane_lower(conj_nx, positive_iy) );
				A(lx, ly, N[3]/2) = conj( plane_upper(conj_nx, positive_iy) );
			}
		}
		
		
		// For ky =0 and ky = N[2]/2 for the above planes, 
		// A(-kx, 0, ..) = A(kx, 0, ..) for kx >0.
		for (int lx=0; lx<local_N1; lx++)
		{
			kx = Get_kx_FOUR(lx, N);
			
			if (kx <0)			// kx < 0 for lx location.
			{
				positive_ix =  Get_ix_FOUR(-kx, N);	// conj_nx has kx > 0.
				
				A(lx, 0, 0) = conj( plane_lower(positive_ix, 0) );
				A(lx, 0, N[3]/2) = conj( plane_upper(positive_ix, 0) );
				
			}	
		}	
	}	
	
	else if (basis_type == "SCFT")		// f(m, -ky, 0) = conj(f(m, ky, 0))- do for kz=N[3]/2
	{
		for (int ly=N[2]/2+1; ly<N[2]; ly++)
		{
			positive_iy = -Get_ky3D_SCFT(ly, N);			
			A(Range::all(), ly, 0) = conj( A(Range::all(), positive_iy, 0) );
			A(Range::all(), ly, N[3]/2) = conj( A(Range::all(), positive_iy, N[3]/2) );
		}
	}

}


/**********************************************************************************************

		Dealias array using 2/3 rule: fill up outer 1/3 region with zeros for dealising

***********************************************************************************************/ 

void Dealias_array(string basis_type, int N[], Array<complx,3> A)
{
	int kx;
	
	if (basis_type == "FOUR") 
	{
		
		if (N[1] >= 3)
			for (int l1=0;  l1<local_N1; l1++)
			{
				kx = Get_kx(basis_type, l1, N);
				
				if (abs(kx) >= N[1]/3)	
					A(l1, Range::all(), Range::all()) = 0.0;				
			}
		
		if (N[2] >= 3)
			A(Range::all(), Range(N[2]/3, 2*N[2]/3+1), Range::all())= 0.0;
		
		if (N[3] >= 3)	// N[3] = 2 is really a 2D sim
			A(Range::all(), Range::all(), Range(N[3]/3, N[3]/2))= 0.0;
	}
		
		
	else if (basis_type == "SCFT")
	{
		if (N[1] >= 3)
			for (int l1=0;  l1<local_N1; l1++)
			{
				kx = Get_kx(basis_type, l1, N);
				
				if (abs(kx) >= 2*N[1]/3)	
					(A)(l1, Range::all(), Range::all()) = 0.0;				
			}
		
		if (N[2] >= 3)
			A(Range::all(), Range(N[2]/3, 2*N[2]/3+1), Range::all())= 0.0;
		
		if (N[3] >= 3)	// N[3] = 2 is really a 2D sim
			A(Range::all(), Range::all(), Range(N[3]/3, N[3]/2))= 0.0;
	}	
}	


/**********************************************************************************************

		Test aliasing needed or not in array using 2/3 rule.  The array has nonzero values
		only up to r= outer_radius

***********************************************************************************************/ 

int Is_alias_array(string basis_type, int N[], Array<complx,3> A, DP outer_radius, DP kfactor[])
{
	int kx_max = (int) ceil(outer_radius/kfactor[1]);
	int ky_max = (int) ceil(outer_radius/kfactor[2]);
	int kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	if (basis_type == "FOUR") 
	{
		if ((kx_max > N[1]/3) || (ky_max > N[2]/3) || (kz_max > N[3]/3))
			return 1;
		else
			return 0;
	}
		
	else if (basis_type == "SCFT")
	{
		if ((kx_max > 2*N[1]/3) || (ky_max > N[2]/3) || (kz_max > N[3]/3))
			return 1;
		else
			return 0;
	}
	
	else
		return 0;  // for -Wall				
}	


//********************************  End of universal_misc.cc   ********************************


