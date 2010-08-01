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

/*! \file four_basic.cc
 * 
 * @sa four_basic.h
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date	22/08/2008
 * @bug		No known bugs
 */ 

#include "four_basic.h"


/**********************************************************************************************

	-- Number of modes included inside shell (inner_radius, outer_radius]
	-- also count the complex conj
	
***********************************************************************************************/


int Get_local_number_modes_in_shell_FOUR(int N[], DP inner_radius, DP outer_radius, DP kfactor[])
{
	
	int l1, l2, l3;
	
	int kx_max, ky_max, kz_max;
	
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (N[2] > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0.0;
	
	if (N[3] >=2)
		kz_max = (int) ceil(outer_radius/kfactor[3]);
	else
		kz_max = 0.0;
	
	
	int count = 0;
	DP kkmag;
	
	for (int kx = -kx_max; kx <= kx_max; kx++)
		for (int ky = -ky_max; ky <= ky_max; ky++)  
			for (int kz = 0; kz <= kz_max; kz++)
			{
				l1 = Get_lx_FOUR(kx, N);
				l2 = Get_ly3D_FOUR(ky, N);
				l3 = kz;
				
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				if ((kkmag >= inner_radius) && (kkmag < outer_radius))
				{
					if (l3 == 0)
						count++;
					else 
						count = count + 2;
				}		
			}
	
	return count;		

}


int Get_number_modes_in_shell_FOUR(int N[], DP inner_radius, DP outer_radius, DP kfactor[])
{
	int local_count = Get_local_number_modes_in_shell_FOUR(N, inner_radius, outer_radius, 
															kfactor);
	
	int count;
	
	MPI_Allreduce(&local_count, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	return count;
}	

/**********************************************************************************************

	Compute Wavenumber components for grid index (i1,i2).

***********************************************************************************************/

// Real K
void Wavenumber_FOUR(int l1, int l2, int l3, int N[], DP kfactor[], TinyVector<DP,3> &kk)
{
	kk = Get_kx_FOUR(l1, N)*kfactor[1],  Get_ky3D_FOUR(l2, N)*kfactor[2], l3*kfactor[3];
}


// Complex K; The imaginary part is zero.  Written to use cross function of blitz.
// Omega = cross(V,K).
void Wavenumber_FOUR(int l1, int l2, int l3, int N[], DP kfactor[], TinyVector<complx,3> &kk)
{
	kk = complx(Get_kx_FOUR(l1, N)*kfactor[1], 0.0),
		 complx(Get_ky3D_FOUR(l2, N)*kfactor[2], 0.0),
		 complx(l3*kfactor[3], 0.0);
}


/**********************************************************************************************

	Compute Modal helicity for grid index (i1, i2, i3) = K. [Vr x Vi].

***********************************************************************************************/

DP Get_Modal_helicity_FOUR
(
	int l1, int l2, int l3, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP kfactor[]
)
{

	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, kk;
	
	Vreal = real(Ax(l1, l2, l3)), real(Ay(l1, l2, l3)), real(Az(l1, l2, l3));
	Vimag = imag(Ax(l1, l2, l3)), imag(Ay(l1, l2, l3)), imag(Az(l1, l2, l3));
		
	VrcrossVi = cross(Vreal, Vimag);
	Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
	
	return (dot(kk, VrcrossVi));	

}



/**********************************************************************************************=

	Compute Modal Vorticity for grid index (i1, i2, i3) = I K x V(k).

***********************************************************************************************/


void Compute_Modal_vorticity_FOUR
(
	int l1, int l2, int l3, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP kfactor[], 
	TinyVector<complx,3> &vorticity
)
{

	TinyVector<DP,3> kk;
	TinyVector<complx,3> Vi, complx_kk;
	
	Vi = Ax(l1, l2, l3), Ay(l1, l2, l3), Az(l1, l2, l3);

	Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
	
	vorticity(0) = I *  (kk(1)*Vi(2) - kk(2)*Vi(1));
	vorticity(1) = I *  (kk(2)*Vi(0) - kk(0)*Vi(2));
	vorticity(2) = I *  (kk(0)*Vi(1) - kk(1)*Vi(0));
}




void Compute_Modal_vorticity_y_component_FOUR
(
 int l1, int l2, int l3, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
 DP kfactor[], 
 complx &vort_y
 )
{
	
	TinyVector<DP,3> kk;
	TinyVector<complx,3> Vi, complx_kk;
	
	Vi = Ax(l1, l2, l3), Ay(l1, l2, l3), Az(l1, l2, l3);
	
	Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
	
	vort_y = I *  (kk(2)*Vi(0) - kk(0)*Vi(2));
}


/**********************************************************************************************

	Replaces A(k) by A(k)/K^2 with K^2 = sum_i [ki*kfactor(i)]^2 ; A(0)=0

***********************************************************************************************/



void Array_divide_ksqr_FOUR(int N[], Array<complx,3> A, DP kfactor[])
{
	firstIndex l2;
	secondIndex l3;
  
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 
	{
		k1 = Get_kx_FOUR(l1, N);
		
		A(l1, Range(0,N[2]/2), Range::all()) 
			= A(l1, Range(0,N[2]/2), Range::all())
				/( pow2(k1*kfactor[1]) + pow2(l2*kfactor[2]) + pow2(l3*kfactor[3]));
  
		if (N[2] > 1)
			A(l1, Range(N[2]/2+1,N[2]-1), Range::all())  
				= A(l1, Range(N[2]/2+1,N[2]-1), Range::all())
					/( pow2(k1*kfactor[1]) + pow2((l2+1-N[2]/2) * kfactor[2]) 
										   + pow2(l3*kfactor[3]));
	}
	   
	// To avoid division by zero
	if (my_id == master_id)   
		A(0,0,0) = 0.0;   
}

/**********************************************************************************************

       		Replaces A(k) by A(k)*K^2.

***********************************************************************************************/


void Array_mult_ksqr_FOUR(int N[], Array<complx,3> A, DP kfactor[])
{
	firstIndex l2;
	secondIndex l3;
  
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 
	{
		k1 = Get_kx_FOUR(l1, N);	
			
		A(l1, Range(0,N[2]/2), Range::all()) 
			= A(l1, Range(0,N[2]/2), Range::all())
				* ( pow2(k1*kfactor[1]) + pow2(l2*kfactor[2]) + pow2(l3*kfactor[3]));
  
		if (N[2] > 1)
			A(l1, Range(N[2]/2+1,N[2]-1), Range::all())  
				= A(l1, Range(N[2]/2+1,N[2]-1), Range::all())
					* ( pow2(k1*kfactor[1]) + pow2((l2+1-N[2]/2) * kfactor[2]) 
											+ pow2(l3*kfactor[3]));
	}  
}



/**********************************************************************************************

	Replaces A(k) by A(k)*exp(factor*K^2).

***********************************************************************************************/



void Array_exp_ksqr_FOUR(int N[], Array<complx,3> A, DP factor, DP kfactor[])
{
	firstIndex l2;
	secondIndex l3;
  
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 
	{
		k1 = Get_kx_FOUR(l1, N);	
			
		A(l1, Range(0,N[2]/2), Range::all()) 
			= A(l1, Range(0,N[2]/2), Range::all())
				* exp(factor* ( pow2(k1*kfactor[1]) + pow2(l2*kfactor[2]) 
										+ pow2(l3*kfactor[3])) );
  
		if (N[2] > 1)
			A(l1, Range(N[2]/2+1,N[2]-1), Range::all())  
				= A(l1, Range(N[2]/2+1,N[2]-1), Range::all()) 
					 * exp(factor* ( pow2(k1*kfactor[1]) + pow2((l2+1-N[2]/2)*kfactor[2]) 
										+ pow2(l3*kfactor[3])));
	}
}


/**********************************************************************************************

	Replaces A(k) by A(k)*[ exp( factor*K^2+hyper_factor*K^4 )]

***********************************************************************************************/

void Array_exp_ksqr_FOUR(int N[], Array<complx,3> A, DP factor, DP hyper_factor, DP kfactor[])
{
	firstIndex l2;
	secondIndex l3;
  
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 
	{
		k1 = Get_kx_FOUR(l1, N);	
			
		A(l1, Range(0,N[2]/2), Range::all()) 
			= A(l1, Range(0,N[2]/2), Range::all())
				* exp( factor* (pow2(k1*kfactor[1]) + pow2(l2*kfactor[2]) +pow2(l3*kfactor[3])) 
					+ hyper_factor* sqr(pow2(k1*kfactor[1]) + pow2(l2*kfactor[2]) 
										+ pow2(l3*kfactor[3])) );
  
		if (N[2] > 1)
			A(l1, Range(N[2]/2+1,N[2]-1), Range::all())  
				= A(l1, Range(N[2]/2+1,N[2]-1), Range::all()) 
					 * exp(factor* (pow2(k1*kfactor[1]) + pow2((l2+1-N[2]/2)*kfactor[2]) 
										+ pow2(l3*kfactor[3]))
						  + hyper_factor* sqr(pow2(k1*kfactor[1]) + pow2((l2+1-N[2]/2)*kfactor[2]) 
																	+ pow2(l3*kfactor[3])) );
	}
  
}

/**********************************************************************************************

	Replaces A(k) by A(k)*(V0.K)^2 / K^2.

***********************************************************************************************/


void Array_mult_V0_khat_sqr_FOUR(int N[], Array<complx,3> A, TinyVector<DP,3> V0, DP kfactor[])
{
	firstIndex l2;
	secondIndex l3;
  
	DP V0x = V0(0);
	DP V0y = V0(1);
	DP V0z = V0(2);
	
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 
	{
		k1 = Get_kx_FOUR(l1, N);
				
		A(l1, Range(0,N[2]/2), Range::all()) 
			= A(l1, Range(0,N[2]/2), Range::all())
				* sqr( V0x*k1*kfactor[1] + V0y*l2*kfactor[2] + V0z*l3*kfactor[3] );
			
		if (N[2] > 1)
			A(l1, Range(N[2]/2+1,N[2]-1), Range::all())  
				= A(l1, Range(N[2]/2+1,N[2]-1), Range::all())
					*sqr( V0x*k1*kfactor[1] + V0y*(l2+1-N[2]/2)*kfactor[2] + V0z*l3*kfactor[3] );
	}  
	
	Array_divide_ksqr_FOUR(N, A, kfactor);
}


//*****************************  End of four_basic.cc *****************************************	








