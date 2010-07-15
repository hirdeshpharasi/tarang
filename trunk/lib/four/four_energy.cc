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


/*! \file four_energy.cc
 * 
 * @sa four_energy.h
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 22/08/2008
 * @bug Line 1299: Helicity -- do we need factor 2?
 */  

#include "four_energy.h"


/**********************************************************************************************

       		Computes A^2/2 & Re(A.B*)/2 except k=0

***********************************************************************************************/


DP Get_local_energy_FOUR(string alias_switch, int N[], Array<complx,3> A)
{
	DP total = 0.0;
	
		// kz >= 0
	total =   sum(sqr(abs(A)));
	
	if (N[2] > 1)
		total +=  sum(sqr(abs(A(Range::all(), N[2]/2, Range::all()))));			
		// for ky = -N[2]/2 modes that are not stored.
	
	if (numprocs == 1)
		total += sum(sqr(abs(A(N[1]/2, Range::all(), Range::all()))));		
		// for kx = -N[1]/2 modes that are not stored.
	
	else if (my_id == numprocs/2)
		total += sum(sqr(abs(A(0, Range::all(), Range::all()))));			
		// for kx = -N[1]/2 modes that are not stored.
	
	
		// kz = 0: 	 subtract 1/2(...) 
	total += - sum(sqr(abs(A(Range::all(), Range::all(), 0))))/2;
	
	if (N[2] > 1)
		total -= sum(sqr(abs(A(Range::all(), N[2]/2, 0))))/2;				// ky = -N[2]/2
	
	if (numprocs == 1)
		total +=  -sum(sqr(abs(A(N[1]/2, Range::all(), 0))))/2;				// kx = -N[1]/2	
	
	else if (my_id == numprocs/2)
		total +=  -sum(sqr(abs(A(0, Range::all(), 0))))/2;					// kx = -N[1]/2	
	
			 	
	// Origin			
	if (my_id == master_id)  
		total += - pow2(abs(A(0,0,0)))/2;  
	
	total =   sum(sqr(abs(A))) ;
		
	return total;
}

//


DP Get_total_energy_FOUR(string alias_switch, int N[], Array<complx,3> A)
{

	DP local_total = Get_local_energy_FOUR(alias_switch, N, A);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return total;
		
	else
		return 0;
}

//*********************************************************************************************

DP Get_local_energy_FOUR(string alias_switch, int N[],  Array<complx,3> A, Array<complx,3> B)
{
	DP total = 0.0;
	
	// kz >= 0
	total =   real(sum(A*conj(B)))
			+ real( sum( A(Range::all(), N[2]/2, Range::all())
						* conj(B(Range::all(), N[2]/2, Range::all()))) );		
						// for ky = -N[2]/2
	
	if (numprocs == 1)						// for kx =  -N[1]/2
		total += real( sum( A(N[1]/2, Range::all(), Range::all())
							* conj(B(N[1]/2, Range::all(), Range::all()))) );
		
	else if (my_id == numprocs/2)			// for kx = -N[1]/2
		total += real( sum(A(0, Range::all(), Range::all())
							* conj(B(0, Range::all(), Range::all()))) );
		  
	
	// kz = 0: 	 subtract 1/2(...) 
	total += - real( sum(A(Range::all(), Range::all(), 0)
							* conj(B(Range::all(), Range::all(), 0))) )/2
			 - real( sum(A(Range::all(), N[2]/2, 0)* conj(B(Range::all(), N[2]/2, 0))) )/2;
			 
	if (numprocs == 1)						// for kx = -N[1]/2
		total +=  - real( sum(A(N[1]/2, Range::all(), 0)
								* conj(B(N[1]/2, Range::all(), 0))) )/2;
		
	else if (my_id == numprocs/2)			// for kx = -N[1]/2
		total +=  - real( sum(A(0, Range::all(), 0)*conj(B(0, Range::all(), 0))) )/2;
	
	
	// Origin
	if (my_id == master_id) 		 			 
		total += - real( A(0,0,0)*conj(B(0,0,0)) )/2;							
	
	
	return total;	  
}


// 

DP Get_total_energy_FOUR(string alias_switch, int N[],  Array<complx,3> A, Array<complx,3> B)
{

	DP local_total = Get_local_energy_FOUR(alias_switch, N, A, B);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return total;
		
	else
		return 0;	
}


/**********************************************************************************************

       		Computes sum[ K^n |A(k)|^2/2];   k=0 excluded

***********************************************************************************************/


//  Computes \f$ \sum  K^n |A(k)|^2/2 \f$;   k=0 excluded.
DP Get_local_Sn_FOUR(string alias_switch, int N[], Array<complx,3> A, DP n, DP kfactor[])
{
	DP Sn = 0.0;
	DP kkmag;												// kkmag = sqrt(Kx^2+Ky^2+Kz^2)
	
	int	kkmax = Min_radius_outside_FOUR(alias_switch, N, kfactor);
	  
	for (int l1=0; l1<local_N1; l1++)			
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				if (kkmag <= kkmax)
					Sn += Multiplicity_factor_FOUR(l1,l2,l3,N) * pow(kkmag,n) 
														   * pow2(abs(A(l1,l2,l3)));
			}
			
	// Subtract the contribution from the origin
	if (my_id == master_id)
		Sn += -Multiplicity_factor_FOUR(0, 0, 0, N) * pow(0,n) * pow2(abs(A(0,0)));		

	return Sn;
}

//

DP Get_total_Sn_FOUR(string alias_switch, int N[], Array<complx,3> A, DP n, DP kfactor[])
{
	DP local_total;
	local_total = Get_local_Sn_FOUR(alias_switch, N, A, n, kfactor);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return total;
		
	else
		return 0;	
}


/**********************************************************************************************
 
		Computes Sk(K) = sum[ K'^n|A(K')|^2/2 ] with K-1 < K'<= K.
             Sk(0) for n=0 is the mean energy.
			  
			  Sk(A, B, k) = sum[ K'^n Real(A(K')*conj(B(K'))/2 ] 
		
			  For the partial shell: avg* volume of the shell.
			  The volume is 2*pi*K^2.

***********************************************************************************************/


void Compute_local_shell_spectrum_FOUR
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, 
	DP n, 
	Array<DP,1> local_Sk, 
	Array<DP,1> local_Sk_count,  
	DP kfactor[]
)
{
	local_Sk_count = 0.0;
	local_Sk = 0.0;												

	DP kkmag;													
	int index;
	DP factor;
	
	int	kkmax = Min_radius_outside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				index = (int) ceil(kkmag);
				
				if (index <= kkmax)
				{
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
				
					local_Sk(index) = local_Sk(index) + factor * pow(kkmag,n) 
																* pow2(abs(A(l1,l2,l3)));
					local_Sk_count(index) = local_Sk_count(index) + factor;
				}
			}

}


void Compute_shell_spectrum_FOUR
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, 
	DP n, 
	Array<DP,1> Sk,  
	DP kfactor[]
)
{

	static Array<DP,1> local_Sk(Sk.length());
	static Array<DP,1> local_Sk_count(Sk.length());
	
	local_Sk_count = 0.0; 
	local_Sk = 0.0;	
	
	Compute_local_shell_spectrum_FOUR(alias_switch, N, A, n, local_Sk, local_Sk_count, kfactor);
	
	static Array<DP,1> Sk_count(Sk.length());
	
	int data_size = Sk.size();
				
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
					reinterpret_cast<double*>(Sk.data()), 
					data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
					
	MPI_Reduce(reinterpret_cast<double*>(local_Sk_count.data()), 
					reinterpret_cast<double*>(Sk_count.data()), 
					data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 

	// For semi-filled shells
	
	if (my_id == master_id) 
	{
		int kkmax_inside = Max_radius_inside_FOUR(alias_switch, N, kfactor); 	
		int	kkmax = Min_radius_outside_FOUR(alias_switch, N, kfactor);

		if (N[2] > 1)   // N[2]=1 is a 2D case
		{
			for (int index = kkmax_inside+1; index <= kkmax; index++) 
				if (Sk_count(index) >= 1) 
					Sk(index) = Sk(index) * Approx_number_modes_in_shell_FOUR(index, kfactor) 
											/ Sk_count(index); 
					// Approx_modes in shell of radius s is 2*pi*s^2/f1*f2*f3 for WAVENOACTUAL 
					//	(unit volume = f1*f2*f3).
					// For WAVENOGRID, the value is pi*s. 
		}
		
		else if (N[2] == 1)
		{
			DP approx_number_modes_in_shell_2D;
			
			for (int index = kkmax_inside+1; index <= kkmax; index++) 
			{
				approx_number_modes_in_shell_2D = (M_PI*index)/(kfactor[1]*kfactor[3]);
				
				if (Sk_count(index) >= 1) 
					Sk(index) = Sk(index) * approx_number_modes_in_shell_2D	/ Sk_count(index); 
			}	
		}	
	}
	
}
			

//*********************************************************************************************			
//
// Re(A. conj(B))
//


void Compute_local_shell_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> A,  Array<complx,3> B,
	DP n, 
	Array<DP,1> local_Sk, 
	Array<DP,1> local_Sk_count, 
	DP kfactor[]
)
{
	local_Sk_count=0;
	local_Sk = 0.0;															
	
	DP kkmag;																
	int index;
	DP factor;

	int	kkmax = Min_radius_outside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				index = (int) ceil(kkmag);
				
				if (index <= kkmax)
				{
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
				
					local_Sk(index) = local_Sk(index) + factor* pow(kkmag,n) * real(A(l1,l2,l3)
																		 *conj(B(l1,l2,l3)));
					local_Sk_count(index) = local_Sk_count(index) + factor;
				}	
			}

}

//
//
void Compute_shell_spectrum_FOUR
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP n, 
	Array<DP,1> Sk,  
	DP kfactor[]
)
{

	static Array<DP,1> local_Sk(Sk.length());
	static Array<DP,1> local_Sk_count(Sk.length());
	
	local_Sk_count = 0.0; 
	local_Sk = 0.0;	
	Compute_local_shell_spectrum_FOUR(alias_switch, N, A, B, n, local_Sk, local_Sk_count, 
										kfactor);
	
	static Array<DP,1> Sk_count(Sk.length());
	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
						reinterpret_cast<double*>(Sk.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
						
	MPI_Reduce(reinterpret_cast<double*>(local_Sk_count.data()), 
						reinterpret_cast<double*>(Sk_count.data()), 
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 

	// For semi-filled shells
	
	if (my_id == master_id) 
	{
		int kkmax_inside = Max_radius_inside_FOUR(alias_switch, N, kfactor); 	
		int	kkmax = Min_radius_outside_FOUR(alias_switch, N, kfactor);

		if (N[2] > 1)   // N[2]=1 is a 2D case
		{	
			for (int index =  kkmax_inside+1; index <= kkmax; index++) 
				if (Sk_count(index) >= 1) 
					Sk(index) = Sk(index) *Approx_number_modes_in_shell_FOUR(index, kfactor) 
											/ Sk_count(index); 
		}
		
		else if (N[2] == 1)
		{
			DP approx_number_modes_in_shell_2D;
			
			for (int index =  kkmax_inside+1; index <= kkmax; index++) 
			{	
				approx_number_modes_in_shell_2D = (M_PI*index)/(kfactor[1]*kfactor[3]);
				
				if (Sk_count(index) >= 1) 
					Sk(index) = Sk(index) *approx_number_modes_in_shell_2D 	/ Sk_count(index); 
			}	
		}	
			
	}
}



/**********************************************************************************************

	Compute Helicity1 = K . (Vr x Vi)
	Helicity2 = K. (Vr x Vi)/K^2
	Multiplication factor = 1 for kz>0 because of energy spectrum details.

***********************************************************************************************/



void Compute_local_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &local_helicity1, DP &local_helicity2, 
	DP &local_dissipation_H1, DP &local_dissipation_H2,
	DP kfactor[]
)
{
	TinyVector<DP,3> Vreal, Vimag, VrcrossVi;
	TinyVector<DP,3> kk;
	DP modal_helicity, kkmag, kksqr, factor;
	
	local_helicity1 = local_helicity2 = 0.0;
	local_dissipation_H1 =  0.0;
	
	int	kkmax = Min_radius_outside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++)
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				if (kkmag <= kkmax)
				{
					Vreal = real(Ax(l1, l2, l3)), real(Ay(l1, l2, l3)), real(Az(l1, l2, l3));
					Vimag = imag(Ax(l1, l2, l3)), imag(Ay(l1, l2, l3)), imag(Az(l1, l2, l3));
			
					VrcrossVi = cross(Vreal, Vimag);
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);		
					kksqr = pow2(kkmag);
					
					factor = 2*Multiplicity_factor_FOUR(l1, l2, l3, N);
					// factor multiplied by 2 because of the defn  Hk = K . (Vr x Vi).
					// recall the defn of energy spectrum that contains 1/2.
					
					modal_helicity = factor * dot(kk, VrcrossVi);
					local_helicity1 += modal_helicity;
					
					if (kksqr > MYEPS)
						local_helicity2 += modal_helicity / kksqr;
						
					local_dissipation_H1 += kksqr * modal_helicity;
				}	
			}	
			
	local_dissipation_H1 *= 2.0;
	local_dissipation_H2 = 2.0 * local_helicity1;		
}

//

void Compute_total_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &total_helicity1, DP &total_helicity2, 
	DP &total_dissipation_H1, DP &total_dissipation_H2,
	DP kfactor[]
)
{
	DP local_helicity1, local_helicity2;
	DP local_dissipation_H1, local_dissipation_H2;
	
	Compute_local_helicity_FOUR(alias_switch, N, Ax, Ay, Az, 
									local_helicity1, local_helicity2, 
									local_dissipation_H1, local_dissipation_H2, kfactor);
									
	MPI_Reduce(&local_helicity1, &total_helicity1, 1, MPI_DOUBLE, MPI_SUM, master_id, 
								MPI_COMM_WORLD);
								
	MPI_Reduce(&local_helicity2, &total_helicity2, 1, MPI_DOUBLE, MPI_SUM, master_id, 
								MPI_COMM_WORLD);
								
	MPI_Reduce(&local_dissipation_H1, &total_dissipation_H1, 1, MPI_DOUBLE, MPI_SUM, master_id, 
								MPI_COMM_WORLD);
									
	MPI_Reduce(&local_dissipation_H2, &total_dissipation_H2, 1, MPI_DOUBLE, MPI_SUM, master_id, 
								MPI_COMM_WORLD);
}


/**********************************************************************************************

	Compute helicity spectrum
	Helicity1 = K . (Vr x Vi)
	Helicity2 = K. (Vr x Vi)/K^2

***********************************************************************************************/



void Compute_local_shell_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3, 
	Array<DP,1> local_H1k_count, 
	DP kfactor[]
)												
{

	local_H1k_count = 0.0;
	local_H1k1 = 0.0;	
	local_H1k2 = 0.0;
	local_H1k3 = 0.0;
	
	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, kk;
	DP kkmag;															
	int index;
	DP factor;
	
	int	kkmax = Min_radius_outside_FOUR(alias_switch,N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				index = (int) ceil(kkmag);
				
				if (index <= kkmax) 
				{
					factor = 2*Multiplicity_factor_FOUR(l1, l2, l3, N);
					// factor multiplied by 2 because of the defn  Hk = K . (Vr x Vi).
					// recall the defn of energy spectrum that contains 1/2.
					
					Vreal = real(Ax(l1, l2, l3)), real(Ay(l1, l2, l3)), real(Az(l1, l2, l3));
					Vimag = imag(Ax(l1, l2, l3)), imag(Ay(l1, l2, l3)), imag(Az(l1, l2, l3));
			
					VrcrossVi = cross(Vreal, Vimag);
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
					
					// modal_helicity = factor * dot(kk, VrcrossVi);	
					local_H1k1(index) = local_H1k1(index) + factor* (kk(0)*VrcrossVi(0));
					local_H1k2(index) = local_H1k2(index) + factor* (kk(1)*VrcrossVi(1));
					local_H1k3(index) = local_H1k3(index) + factor* (kk(2)*VrcrossVi(2));
					
					local_H1k_count(index) = local_H1k_count(index) + factor;
				}	
			} 

}


void Compute_shell_spectrum_helicity_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3, 
	DP kfactor[]
)	
{

	static Array<DP,1> local_H1k1(H1k1.length());
	static Array<DP,1> local_H1k2(H1k1.length());
	static Array<DP,1> local_H1k3(H1k1.length());
	
	static Array<DP,1> local_H1k_count(H1k1.length());
	
	local_H1k1 = 0.0;
	local_H1k2 = 0.0;
	local_H1k3 = 0.0;
	
	local_H1k_count = 0.0; 
		
	
	Compute_local_shell_spectrum_helicity_FOUR(alias_switch, N, Ax, Ay, Az, 
							local_H1k1, local_H1k2, local_H1k3,
							local_H1k_count, kfactor);
	
	static Array<DP,1> H1k1_count(H1k1.length());	
	int data_size = H1k1.size();
				
	MPI_Reduce(reinterpret_cast<double*>(local_H1k1.data()), 
					reinterpret_cast<double*>(H1k1.data()), 
					data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);
					
	MPI_Reduce(reinterpret_cast<double*>(local_H1k2.data()), 
					reinterpret_cast<double*>(H1k2.data()), 
					data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);
	
	MPI_Reduce(reinterpret_cast<double*>(local_H1k3.data()), 
					reinterpret_cast<double*>(H1k3.data()), 
					data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);								  								  
					
	MPI_Reduce(reinterpret_cast<double*>(local_H1k_count.data()), 
					reinterpret_cast<double*>(H1k1_count.data()), 
					data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 

	// The shells near the edges do not complete half sphere, so normalize the shells.
	
	if (my_id == master_id) 
	{
		int kkmax_inside = Max_radius_inside_FOUR(alias_switch, N, kfactor); 	
		int	kkmax = Min_radius_outside_FOUR(alias_switch,N, kfactor);

		if (N[2] > 1)   // N[2]=1 is a 2D case
		{
			for (int index = kkmax_inside+1; index <= kkmax; index++) 
				if (H1k1_count(index) >= 1) 
				{
					H1k1(index) = H1k1(index) * Approx_number_modes_in_shell_FOUR(index, kfactor) 
											/ H1k1_count(index); 
											
					H1k2(index) = H1k2(index) * Approx_number_modes_in_shell_FOUR(index, kfactor) 
											/ H1k1_count(index); 
											
					H1k3(index) = H1k3(index) * Approx_number_modes_in_shell_FOUR(index, kfactor) 
											/ H1k1_count(index); 
				}
		}
		
		else if (N[2] == 1)
		{
			DP approx_number_modes_in_shell_2D;
			
			for (int index = kkmax_inside+1; index <= kkmax; index++) 
			{	
				approx_number_modes_in_shell_2D = (M_PI*index)/(kfactor[1]*kfactor[3]);
				
				if (H1k1_count(index) >= 1) 
				{
					H1k1(index) = H1k1(index) * approx_number_modes_in_shell_2D 
									/ H1k1_count(index);
					
					H1k2(index) = H1k2(index) * approx_number_modes_in_shell_2D 
									/ H1k1_count(index); 
					
					H1k3(index) = H1k3(index) * approx_number_modes_in_shell_2D 
									/ H1k1_count(index); 
				}
			}	
		}	
			
	}
	
}



/**********************************************************************************************

	Compute entropy = sum(pi log(1/pi) where pi = Ei/Total E
	
	Skip the mean mode k= (0,0,0)

***********************************************************************************************/



DP Get_local_entropy_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
)
{
	
	DP modal_energy, prob, local_entropy;
	
	DP total_energy = Get_total_energy_FOUR(alias_switch, N, Ax) 
						+ Get_total_energy_FOUR(alias_switch, N, Ay) 
						+ Get_total_energy_FOUR(alias_switch, N, Az);
	
	MPI_Bcast( &total_energy, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
	
	local_entropy = 0.0;
			
	for (int l1=0; l1<local_N1; l1++)									
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				modal_energy = Modal_energy_FOUR(Ax, l1, l2, l3) 
								+ Modal_energy_FOUR(Ay, l1, l2, l3) 
								+  Modal_energy_FOUR(Az, l1, l2, l3);
								
				prob = modal_energy / total_energy;		
				
				if (prob > MYEPS) 
					local_entropy += 2 * Multiplicity_factor_FOUR(l1, l2, l3, N) 
										* (-prob * log(prob)/log(2.0));
					// 2*Multiplicity_factor_FOUR for the complex conjugate mode
					// factor 2 is bit convoluted.  
					// This is because of the definition of Multiplicity_factor_FOUR.
					// Also skips origin because modal_energy(k=0) = 0.
			}
			
	return local_entropy;		
}											


DP Get_entropy_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
)
{

	DP local_entropy = Get_local_entropy_FOUR(alias_switch, N, Ax, Ay, Az);
	
	DP entropy;
	MPI_Reduce(&local_entropy, &entropy, 1, MPI_DOUBLE, MPI_SUM, 
				master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return entropy;
		
	else
		return 0;	
}


//
// Scalar
// 


DP Get_local_entropy_scalar_FOUR(string alias_switch, int N[], Array<complx,3> A)
{

	DP modal_energy, prob;
	DP local_entropy = 0.0;
	
	DP total_energy = Get_total_energy_FOUR(alias_switch, N, A);
	
	MPI_Bcast( &total_energy, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD);
			
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=10; l3<=N[3]/2; l3++) 
			{
				modal_energy = Modal_energy_FOUR(A, l1, l2, l3);
				
				prob = modal_energy / total_energy;		
				
				if (prob > MYEPS) 
					local_entropy += 2 * Multiplicity_factor_FOUR(l1, l2, l3, N) 
								 * (-prob * log(prob)/log(2.0));
					// 2*Multiplicity_factor_FOUR for the complex conjugate mode
					// factor 2 is bit convoluted.  
					// This is because of the definition of Multiplicity_factor_FOUR
					// Also skips origin because modal_energy(k=0) = 0.
			}
			
	return local_entropy;		

}				


DP Get_entropy_scalar_FOUR(string alias_switch, int N[], Array<complx,3> A)
{

	DP local_entropy = Get_local_entropy_scalar_FOUR(alias_switch, N, A);
	
	DP entropy;
	MPI_Reduce(&local_entropy, &entropy, 1, MPI_DOUBLE, MPI_SUM,
				master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return entropy;
	
	else
		return 0;		
}



/**********************************************************************************************

	Compute ring spectrum: Craya basis.
	2D: Sk(k, theta)
	3d: Sk1(k, theta) along horzontal direction; S2k(k, theta) 
			in the perpendicular dirn (polar).

***********************************************************************************************/



void Compute_local_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
)											
{

	local_S1k = 0.0;	
	local_S2k = 0.0;
	
	complx anisV1, anisV2;
	TinyVector<complx,3>  kk;
	TinyVector<complx,3> V, VcrossK;
	
	DP kkmag, theta, kkperp;
	DP V2sqr;
	DP factor;
	int shell_index, sector_index;
	
	int	kkmax = Max_radius_inside_FOUR(alias_switch, N, kfactor);
									
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N,  kfactor);
				shell_index = (int) ceil(kkmag);
				
				if ((kkmag > MYEPS) && (shell_index <= kkmax))
				{
					theta = AnisKvect_polar_angle_FOUR(l1, l2, l3, N, kfactor);
					
					sector_index = Get_sector_index(theta, sector_angle_array);
					
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
					
					kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
					
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
					
					if (kkperp < MYEPS)			// k on the anisotropy axis.
					{
					
#ifdef ANISDIRN1
						anisV1 = Ay(l1, l2, l3);
						anisV2 = Az(l1, l2, l3);
#endif

#ifdef ANISDIRN2
						anisV1 = Az(l1, l2, l3);
						anisV2 = Ax(l1, l2, l3);
#endif

#ifdef ANISDIRN3
						anisV1 = Ax(l1, l2, l3);
						anisV2 = Ay(l1, l2, l3);
#endif

						local_S1k(shell_index, sector_index) += factor * pow(kkmag,n)
																* pow2(abs(anisV1));
																
						local_S2k(shell_index, sector_index) += factor * pow(kkmag,n)
																* pow2(abs(anisV2));					
					}
				
				
				
					else
					{

						V = Ax(l1, l2, l3), Ay(l1, l2, l3), Az(l1, l2, l3);
					
						VcrossK = cross(V, kk);
				
#ifdef ANISDIRN1
						anisV1 = VcrossK(0)/kkperp;
#endif

#ifdef ANISDIRN2
						anisV1 = VcrossK(1)/kkperp;
#endif

#ifdef ANISDIRN3
						anisV1 = VcrossK(2)/kkperp;
#endif


						
						local_S1k(shell_index, sector_index) += factor * pow(kkmag,n)
																* pow2(abs(anisV1));
						
						V2sqr = pow2(abs(Ax(l1,l2,l3))) +  pow2(abs(Ay(l1,l2,l3))) 
								+  pow2(abs(Az(l1,l2,l3))) - pow2(abs(anisV1));
						
						local_S2k(shell_index, sector_index) +=  factor * pow(kkmag,n) * V2sqr;
				
					}	// of inner else
				}		// of if
			}			// of for loop

}


//*********************************************************************************************

void Compute_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
)
{

	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_ring_spectrum_FOUR(alias_switch, N, Ax, Ay, Az, sector_angle_array, n, 
										local_S1k, local_S2k, kfactor);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_S1k.data()), 
						reinterpret_cast<double*>(S1k.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);
						
	MPI_Reduce(reinterpret_cast<double*>(local_S2k.data()), 
						reinterpret_cast<double*>(S2k.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);					 

}



//*********************************************************************************************
// 
//  A.B
//


void Compute_local_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
)											
{

	local_S1k = 0.0;	
	local_S2k = 0.0;
	
	
	complx anisV1, anisV2, anisW1, anisW2;
	TinyVector<complx,3>  kk, e1K;
	TinyVector<complx,3> V, VcrossK, W, WcrossK;
	
	DP kkmag, theta, kkh1, kkh2, kkperp;
	DP factor;
	int shell_index, sector_index;
	
	int	kkmax = Max_radius_inside_FOUR(alias_switch,N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				shell_index = (int) ceil(kkmag);
				
				if ((kkmag > MYEPS) && (shell_index <= kkmax))
				{
					theta = AnisKvect_polar_angle_FOUR(l1, l2, l3, N, kfactor);
					
					sector_index = Get_sector_index(theta, sector_angle_array);
					
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
				
					kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
					
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
					
					if (kkperp < MYEPS)			// k on the anisotropy axis.
					{
					
#ifdef ANISDIRN1
						anisV1 = Ay(l1, l2, l3);	anisV2 = Az(l1, l2, l3);
						anisW1 = By(l1, l2, l3);	anisW2 = Bz(l1, l2, l3);
#endif

#ifdef ANISDIRN2
						anisV1 = Az(l1, l2, l3);	anisV2 = Ax(l1, l2, l3);
						anisW1 = Bz(l1, l2, l3);	anisW2 = Bx(l1, l2, l3);
#endif

#ifdef ANISDIRN3
						anisV1 = Ax(l1, l2, l3);	anisV2 = Ay(l1, l2, l3);
						anisW1 = Bx(l1, l2, l3);	anisW2 = By(l1, l2, l3);
#endif

						local_S1k(shell_index, sector_index) += factor * pow(kkmag,n)
																* real( anisV1*conj(anisW1) );
						
						local_S2k(shell_index, sector_index) += factor * pow(kkmag,n)
																* real( anisV2*conj(anisW2) );					
					}
				
				
				
					else
					{

						V = Ax(l1, l2, l3), Ay(l1, l2, l3), Az(l1, l2, l3);
					
						VcrossK = cross(V, kk);
						
						W = Bx(l1, l2, l3), By(l1, l2, l3), Bz(l1, l2, l3);
					
						WcrossK = cross(W, kk);
						
						kkh1 = AnisKh1_FOUR(l1, l2, l3, N, kfactor);
						kkh2 = AnisKh2_FOUR(l1, l2, l3, N, kfactor);
					
#ifdef ANISDIRN1
						anisV1 = VcrossK(0)/kkperp;
						anisW1 = WcrossK(0)/kkperp;
						
						e1K = 0.0, kkh2/kkperp, -kkh1/kkperp;
						
						anisV2 = dot(VcrossK, e1K);
						anisW2 = dot(WcrossK, e1K);
					
#endif

#ifdef ANISDIRN2
						anisV1 = VcrossK(1)/kkperp;
						anisW1 = WcrossK(1)/kkperp;
						
						e1K = -kkh1/kkperp, 0.0, kkh2/kkperp;
						
						anisV2 = dot(VcrossK, e1K);
						anisW2 = dot(WcrossK, e1K);
#endif

#ifdef ANISDIRN3
						anisV1 = VcrossK(2)/kkperp;
						anisW1 = WcrossK(2)/kkperp;
						
						e1K = kkh2/kkperp, -kkh1/kkperp, 0.0;
						
						anisV2 = dot(VcrossK, e1K);
						anisW2 = dot(WcrossK, e1K);
#endif

					
						local_S1k(shell_index, sector_index) += factor * pow(kkmag,n) 
																* real( anisV1*conj(anisW1) );
					
						local_S2k(shell_index, sector_index) += factor * pow(kkmag,n)  
																* real( anisV2*conj(anisW2) );
				
					}	// of inner else
				}		// of if (shell_index <= kkmax)
			}			// of for loop
			
}


//*********************************************************************************************

void Compute_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
)
{

	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_ring_spectrum_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
								sector_angle_array, n, local_S1k, local_S2k, kfactor);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_S1k.data()), 
						reinterpret_cast<double*>(S1k.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);
						
	MPI_Reduce(reinterpret_cast<double*>(local_S2k.data()), 
						reinterpret_cast<double*>(S2k.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);					 

}


//*********************************************************************************************
// scalar
//

void Compute_local_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
)										
{

	local_Sk = 0.0;	
	
	DP kkmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	kkmax = Max_radius_inside_FOUR(alias_switch,N, kfactor);
			
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N,  kfactor);
				shell_index = (int) ceil(kkmag);
				
				if ((kkmag > MYEPS) && (shell_index <= kkmax))
				{
					theta = AnisKvect_polar_angle_FOUR(l1, l2, l3, N, kfactor);
				
					sector_index = Get_sector_index(theta, sector_angle_array);
					
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
				
					local_Sk(shell_index, sector_index) +=  factor * pow(kkmag,n)
															* pow2(abs(F(l1,l2,l3)));
				}						
			}	
	

}


void Compute_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
)
{

	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
			
	Compute_local_ring_spectrum_FOUR(alias_switch, N, F, sector_angle_array, n, 
										local_Sk, kfactor);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
						reinterpret_cast<double*>(Sk.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);

}


//
// F*G
//

void Compute_local_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
)										
{

	local_Sk = 0.0;	
	
	DP kkmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	int	kkmax = Max_radius_inside_FOUR(alias_switch,N, kfactor);
			
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N,  kfactor);
				shell_index = (int) ceil(kkmag);
				
				if ((kkmag > MYEPS) && (shell_index <= kkmax))
				{
					theta = AnisKvect_polar_angle_FOUR(l1, l2, l3, N, kfactor);
				
					sector_index = Get_sector_index(theta, sector_angle_array);
					
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
				
					local_Sk(shell_index, sector_index) +=  factor * pow(kkmag,n)
													* real(F(l1,l2,l3) * conj(G(l1,l2,l3)));
				}						
			}	
	

}


void Compute_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F,
	Array<complx,3> G,  
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
)
{

	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
			
	Compute_local_ring_spectrum_FOUR(alias_switch, N, F, G, sector_angle_array, n, 
										local_Sk, kfactor);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
						reinterpret_cast<double*>(Sk.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);

}



//*********************************************************************************************

void Compute_local_ring_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array,
	Array<DP,2> local_H1k, 
	DP kfactor[]
)												
{

	local_H1k = 0.0;	

	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, kk;
	DP kkmag, theta;														
	DP modal_helicity;
	DP factor;
	int shell_index, sector_index;

	int	kkmax = Max_radius_inside_FOUR(alias_switch,N, kfactor);
		
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N,  kfactor);
				shell_index = (int) ceil(kkmag);
				
				if ((kkmag > MYEPS) && (shell_index <= kkmax))
				{
					theta = AnisKvect_polar_angle_FOUR(l1, l2, l3, N, kfactor); 
					
					sector_index = Get_sector_index(theta, sector_angle_array);
					
					Vreal = real(Ax(l1, l2, l3)), real(Ay(l1, l2, l3)), real(Az(l1, l2, l3));
					Vimag = imag(Ax(l1, l2, l3)), imag(Ay(l1, l2, l3)), imag(Az(l1, l2, l3));
			
					VrcrossVi = cross(Vreal, Vimag);
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
					
					factor = 2*Multiplicity_factor_FOUR(l1, l2, l3, N); 
					// factor multiplied by 2 because of the defn  Hk = K . (Vr x Vi).
					// recall the defn of energy spectrum that contains 1/2.
					
					modal_helicity = factor * dot(kk, VrcrossVi);
					
									
					local_H1k(shell_index, sector_index) +=  modal_helicity;
				}	
			} 

}

//
//

void Compute_ring_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array,
	Array<DP,2> H1k, 
	DP kfactor[]
)	
{

	static Array<DP,2> local_H1k(H1k.length()(0), (H1k.length())(1));
	
	local_H1k = 0.0;	
	
	Compute_local_ring_spectrum_helicity_FOUR(alias_switch, N, Ax, Ay, Az, sector_angle_array, 
												local_H1k, kfactor);
	
	int data_size = (H1k.length())(0) * (H1k.length())(1);
				
	MPI_Reduce(reinterpret_cast<double*>(local_H1k.data()), 
					reinterpret_cast<double*>(H1k.data()), 
					data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
}


/**********************************************************************************************
	
								CYLINDERICAL RING SPECTRUM
				
***********************************************************************************************/


void Compute_local_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
)						
{

	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;

	complx anisV1;
	
	DP kkmag, kkpll, kkperp;										
	DP V2sqr;
	DP factor;
	int shell_index, slab_index;
	
	
	int	kkperp_max = Anis_max_Krho_radius_inside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
				shell_index = (int) ceil(kkperp);
				
				if (shell_index <= kkperp_max)
				{					
					kkpll = AnisKpll_FOUR(l1, l2, l3, N, kfactor);
					
					slab_index = Get_slab_index(kkpll, kkperp, cylinder_kpll_array);
					
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
				
			
#ifdef ANISDIRN1
					anisV1 = Ax(l1, l2, l3);
#endif

#ifdef ANISDIRN2
					anisV1 = Ay(l1, l2, l3);
#endif

#ifdef ANISDIRN3
					anisV1 = Az(l1, l2, l3);
#endif

					
					local_S1k(shell_index, slab_index) += factor * pow(kkmag,n) 
															     * pow2(abs(anisV1));
					
					V2sqr = pow2(abs(Ax(l1,l2,l3))) +  pow2(abs(Ay(l1,l2,l3))) 
							+  pow2(abs(Az(l1,l2,l3))) - pow2(abs(anisV1));
					
					local_S2k(shell_index, slab_index) += factor * pow(kkmag,n) * V2sqr;
				}	
			}
	
}

//
//
void Compute_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
)
{

	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_cylinder_ring_spectrum_FOUR(alias_switch, N, Ax, Ay, Az, 
									cylinder_kpll_array, n, local_S1k, local_S2k, kfactor);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_S1k.data()), 
						reinterpret_cast<double*>(S1k.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);
						
	MPI_Reduce(reinterpret_cast<double*>(local_S2k.data()), 
						reinterpret_cast<double*>(S2k.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);

}

//*********************************************************************************************

void Compute_local_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
)
{
	local_S1k = 0.0;  												// initialize Ek
	local_S2k = 0.0;

	complx V, W, anisV1, anisW1;
	TinyVector<complx,3> Vrho, Wrho_conj;
	TinyVector<DP,3>  e1K;
	
	DP kkmag, kkpll, kkperp;
	DP factor;
	int shell_index, slab_index;

	int	kkperp_max = Anis_max_Krho_radius_inside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
				shell_index = (int) ceil(kkperp);
				
				if (shell_index <= kkperp_max)
				{
					kkpll = AnisKpll_FOUR(l1, l2, l3, N, kfactor);
					
					slab_index = Get_slab_index(kkpll, kkperp, cylinder_kpll_array);
					
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
					
					V = Ax(l1, l2, l3), Ay(l1, l2, l3), Az(l1, l2, l3);
					W = Bx(l1, l2, l3), By(l1, l2, l3), Bz(l1, l2, l3);
			
#ifdef ANISDIRN1
					anisV1 = Ax(l1, l2, l3);
					anisW1 = Bx(l1, l2, l3);
					
					e1K = 1.0, 0.0, 0.0;
					
					Vrho = complx(0.0,0.0), Ay(l1, l2, l3), Az(l1, l2, l3);
					Wrho_conj = complx(0.0,0.0), conj(By(l1, l2, l3)), conj(Bz(l1, l2, l3));					
#endif

#ifdef ANISDIRN2
					anisV1 = Ay(l1, l2, l3);
					anisW1 = By(l1, l2, l3);
					
					e1K = 0.0, 1.0, 0.0;
					
					Vrho = Ax(l1, l2, l3), complx(0.0,0.0), Az(l1, l2, l3);
					Wrho_conj = conj(Bx(l1, l2, l3)), complx(0.0,0.0), conj(Bz(l1, l2, l3));
#endif

#ifdef ANISDIRN3
					anisV1 = Az(l1, l2, l3);
					anisW1 = Bz(l1, l2, l3);
					
					e1K = 0.0, 0.0, 1.0;
					
					Vrho = Ax(l1, l2, l3), Ay(l1, l2, l3), complx(0.0,0.0);
					Wrho_conj = conj(Bx(l1, l2, l3)), conj(By(l1, l2, l3)), complx(0.0,0.0);
#endif
					
					local_S1k(shell_index, slab_index) += factor * pow(kkmag,n) 
														   * real(anisV1*conj(anisW1));
					
					local_S2k(shell_index, slab_index) += factor * pow(kkmag,n) 
														   * real(dot(Vrho, Wrho_conj));
				}											
			}

}	



//
//
void Compute_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
)
{

	static Array<DP,2> local_S1k( (S1k.length())(0), (S1k.length())(1));
	static Array<DP,2> local_S2k( (S1k.length())(0), (S1k.length())(1));
	
	Compute_local_ring_spectrum_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
									cylinder_kpll_array, n, local_S1k, local_S2k, kfactor);
	
	int data_size = (S1k.length())(0) * (S1k.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_S1k.data()), 
						reinterpret_cast<double*>(S1k.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);
						
	MPI_Reduce(reinterpret_cast<double*>(local_S2k.data()), 
						reinterpret_cast<double*>(S2k.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);					 

}


//*********************************************************************************************

void Compute_local_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> cylinder_kpll_array,
	DP n,   
	Array<DP,2> local_Sk, 
	DP kfactor[]
)
{

	local_Sk = 0.0;
	
	DP kkmag, kkperp, kkpll;
	DP factor;
	int shell_index, slab_index;

	int	kkperp_max = Anis_max_Krho_radius_inside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
				shell_index = (int) ceil(kkperp);
				
				if (shell_index <= kkperp_max)
				{
					kkpll = AnisKpll_FOUR(l1, l2, l3, N, kfactor);
					
					slab_index = Get_slab_index(kkpll, kkperp, cylinder_kpll_array);
					
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
					
					local_Sk(shell_index, slab_index) += factor * pow(kkmag,n)  
														  * pow2(abs(F(l1,l2,l3)));
				}	
			}	

}

//
//
void Compute_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> cylinder_kpll_array,
	DP n,   
	Array<DP,2> Sk, 
	DP kfactor[]
)
{

	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_ring_spectrum_FOUR(alias_switch, N, F, cylinder_kpll_array, n, 
										local_Sk, kfactor);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
						reinterpret_cast<double*>(Sk.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);

}


//*********************************************************************************************

void Compute_local_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G,
	Array<DP, 1> cylinder_kpll_array,
	DP n,   
	Array<DP,2> local_Sk, 
	DP kfactor[]
)
{

	local_Sk = 0.0;
	
	DP kkmag, kkperp, kkpll;
	DP factor;
	int shell_index, slab_index;

	int	kkperp_max = Anis_max_Krho_radius_inside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
				shell_index = (int) ceil(kkperp);
				
				if (shell_index <= kkperp_max)
				{
					kkpll = AnisKpll_FOUR(l1, l2, l3, N, kfactor);
					
					slab_index = Get_slab_index(kkpll, kkperp, cylinder_kpll_array);
					
					factor = Multiplicity_factor_FOUR(l1, l2, l3, N);
					
					local_Sk(shell_index, slab_index) += factor * pow(kkmag,n) 
													* real(F(l1,l2,l3) * conj(G(l1,l2,l3)));
				}	
			}	

}

//
//
void Compute_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G,
	Array<DP, 1> cylinder_kpll_array,
	DP n,   
	Array<DP,2> Sk, 
	DP kfactor[]
)
{
	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	Compute_local_ring_spectrum_FOUR(alias_switch, N, F, G, cylinder_kpll_array, n, 
										local_Sk, kfactor);
	
	int data_size = (Sk.length())(0) * (Sk.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
						reinterpret_cast<double*>(Sk.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);

}


//*********************************************************************************************

void Compute_local_cylinder_ring_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_H1k, 
	DP kfactor[]
)
{

	local_H1k = 0.0;  												// initialize Ek

	TinyVector<DP,3> Vreal, Vimag, VrcrossVi, kk;
	DP kkmag, kkpll, kkperp;															
	DP modal_helicity;
	DP factor;
	
	int shell_index, slab_index;
	
	int	kkperp_max = Anis_max_Krho_radius_inside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)										
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{			
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
				
				shell_index = (int) ceil(kkperp);
				
				if (shell_index <= kkperp_max)
				{
					kkpll = AnisKpll_FOUR(l1, l2, l3, N, kfactor);
					
					slab_index = Get_slab_index(kkpll, kkperp, cylinder_kpll_array);
					
					Vreal = real(Ax(l1, l2, l3)), real(Ay(l1, l2, l3)), real(Az(l1, l2, l3));
					Vimag = imag(Ax(l1, l2, l3)), imag(Ay(l1, l2, l3)), imag(Az(l1, l2, l3));
			
					VrcrossVi = cross(Vreal, Vimag);
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
					
					factor = 2*Multiplicity_factor_FOUR(l1, l2, l3, N); 
					// factor multiplied by 2 because of the defn  Hk = K . (Vr x Vi).
					// recall the defn of energy spectrum that contains 1/2.
					
					modal_helicity = factor * dot(kk, VrcrossVi);				
					local_H1k(shell_index, slab_index) +=  modal_helicity;
				}	
			}	

}

//
//
void Compute_cylinder_ring_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> H1k, 
	DP kfactor[]
)
{
	static Array<DP,2> local_H1k(H1k.length()(0), (H1k.length())(1));
	
	local_H1k = 0.0;	
	
	Compute_local_ring_spectrum_helicity_FOUR(alias_switch, N, Ax, Ay, Az, cylinder_kpll_array, 
												local_H1k, kfactor);
	
	int data_size = (H1k.length())(0) * (H1k.length())(1);
				
	MPI_Reduce(reinterpret_cast<double*>(local_H1k.data()), 
					reinterpret_cast<double*>(H1k.data()), 
					data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
	
}

/**********************************************************************************************

			(B0.k) Im[vectA. conj(vectB)]

***********************************************************************************************/

void Compute_local_imag_shell_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> local_Sk, 
	DP kfactor[]
)
{
	DP kkmag;
	int shell_index;
	TinyVector<DP,3> kk;
	
	
	local_Sk = 0.0;
	
	int	kkmax = Min_radius_outside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				
				shell_index = (int) ceil(kkmag);
				
				if (shell_index <= kkmax)
				{
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);	
					
					local_Sk(shell_index) += 2* Multiplicity_factor_FOUR(l1, l2, l3, N) 
											* dot(B0, kk)
											* imag( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
													+ Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
													+ Az(l1,l2,l3)*conj(Bz(l1,l2,l3)) );
													
					// factor multiplied by 2. recall the defn of energy spectrum and factor								
				}									
			}
					
}


//
//
void Compute_imag_shell_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> Sk, 
	DP kfactor[]
)
{

	static Array<DP,1> local_Sk(Sk.length());
	
	local_Sk = 0.0;	
	
	Compute_local_imag_shell_spectrum_B0_FOUR(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz, 
												local_Sk, kfactor);

	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
						reinterpret_cast<double*>(Sk.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
					
}	
				
//******************************** for rings **************************************************

void Compute_local_imag_ring_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
)
{
	DP kkmag, theta;
	int shell_index, sector_index;
	TinyVector<DP,3> kk;
	
	local_Sk = 0.0;
	
	int	kkmax = Max_radius_inside_FOUR(alias_switch,N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
				shell_index = (int) ceil(kkmag);
				
				if ((kkmag > MYEPS) && (shell_index <= kkmax))
				{
					theta = AnisKvect_polar_angle_FOUR(l1, l2, l3, N, kfactor);
				
					sector_index = Get_sector_index(theta, sector_angle_array);	
					
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);	
					
					local_Sk(shell_index, sector_index) += 
								2* Multiplicity_factor_FOUR(l1, l2, l3, N) * dot(B0, kk)
								 * imag( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
										+ Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
										+ Az(l1,l2,l3)*conj(Bz(l1,l2,l3)));										
										
				}									
			}
					
}

//
//
void Compute_imag_ring_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> Sk, 
	DP kfactor[]
)
{

	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	local_Sk = 0.0;	
	
	Compute_local_imag_ring_spectrum_B0_FOUR(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz, 
												sector_angle_array, local_Sk, kfactor);

	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
						reinterpret_cast<double*>(Sk.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
					
}	


//******************************** for cylinder rings *****************************************

void Compute_local_imag_cylinder_ring_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
)
{

	DP kkpll, kkperp;
	int shell_index, slab_index;
	TinyVector<DP,3> kk;
	
	
	local_Sk = 0.0;
	
	int	kkperp_max = Anis_max_Krho_radius_inside_FOUR(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkpll = AnisKpll_FOUR(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_FOUR(l1, l2, l3, N, kfactor);
				
				shell_index = (int) ceil(kkperp);
				
				if (shell_index <= kkperp_max)
				{
					slab_index = Get_slab_index(kkpll, kkperp, cylinder_kpll_array);
					
					Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
					
					local_Sk(shell_index, slab_index) +=  
											   2* Multiplicity_factor_FOUR(l1, l2, l3, N) 
												* dot(B0, kk)
												* imag( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
													+ Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
													+ Az(l1,l2,l3)*conj(Bz(l1,l2,l3)));
														
				}									
			}

}

//
//
void Compute_imag_cylinder_ring_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> Sk, 
	DP kfactor[]
)
{

	static Array<DP,2> local_Sk( (Sk.length())(0), (Sk.length())(1));
	
	local_Sk = 0.0;	
	
	Compute_local_imag_cylinder_ring_spectrum_B0_FOUR(alias_switch, N, B0, Ax, Ay, Az, 
									Bx, By, Bz, cylinder_kpll_array, local_Sk, kfactor);

	
	int data_size = Sk.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_Sk.data()), 
						reinterpret_cast<double*>(Sk.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
					
}	

//*****************************  End of four_energy.cc ****************************************	


