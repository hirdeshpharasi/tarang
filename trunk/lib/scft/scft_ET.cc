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

/*! \file sincosfour_ET.cc 
 * 
 * @brief Contains useful functions for computing products of vector fields over shells and rings (MPI).
 *
 * We compute sum  \f$  \sum[ \Re(A(\vec{K}') * B^{*}(\vec{K}')) ] \f$  for over a region.  
 * This region could be a shell, ring or a cylinderial segment. <BR>
 *
 * @sa sincosfour_ET.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 22/08/2008
 * @bug no known bug
 */


#include "scft_ET.h"


/**********************************************************************************************

					SHELL MULT  Real(A.B*) 

***********************************************************************************************/


DP Local_shell_mult_single_SCFT
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
)
{

	DP kkmag;
	DP result = 0.0;	
	
	for (int l1=0; l1< local_N1; l1++) 
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++)					
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N, kfactor);
				
				if ( (kkmag > inner_radius) && (kkmag <= outer_radius) ) 
					result += 2 * Multiplicity_factor_SCFT(l1, l2, l3, N) 
								* real(A(l1,l2,l3)*conj(B(l1,l2,l3)));
			}
	
	return result;
}

//
//
DP Shell_mult_single_SCFT
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
)
{
	DP local_total;
	
	local_total = Local_shell_mult_single_SCFT(alias_switch, N, A, B, 
									inner_radius, outer_radius, kfactor);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return total;
		
	else
		return 0;	
}


//*********************************************************************************************


void Local_shell_mult_single_real_imag_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP& local_result_real, DP& local_result_imag,
	DP kfactor[] 
)
{
	DP kkmag;
	DP factor;
	
	local_result_real = local_result_imag = 0.0;
	
	for (int l1=0; l1< local_N1; l1++) 
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++)	
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N, kfactor);
				
				if ( (kkmag > inner_radius) && (kkmag <= outer_radius) ) 
				{
					factor = Multiplicity_factor_SCFT(l1, l2, l3, N);
					
					local_result_real += 2*factor * real(A(l1,l2,l3))*real(B(l1,l2,l3)); 
					local_result_imag += 2*factor * imag(A(l1,l2,l3))*imag(B(l1,l2,l3)); 
				}	
			}
	
}

//
//
void Shell_mult_single_real_imag_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP& total_result_real, DP& total_result_imag,
	DP kfactor[]
)
{
	DP local_result_real, local_result_imag;
	
	Local_shell_mult_single_real_imag_SCFT(alias_switch, N, A, B, inner_radius, outer_radius, 
											local_result_real, local_result_imag, kfactor);
	
	MPI_Reduce(&local_result_real, &total_result_real, 1, MPI_DOUBLE, MPI_SUM, 
									master_id, MPI_COMM_WORLD);	
									
	MPI_Reduce(&local_result_imag, &total_result_imag, 1, MPI_DOUBLE, MPI_SUM, 
									master_id, MPI_COMM_WORLD);	

}

/*===========================================================================

	void Shell_termbyterm_mult_SCFT()
	
	-- returns real( A[k]*conj(B[k]) ) for k = [inner_radius, outer_radius)
		
	-- Modes on the surface of inner_radius is included in the shell.

=============================================================================*/ 

void Local_shell_mult_all_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result, 
	DP kfactor[]
)
{
	DP kkmag;
	int shell_index;
	
	local_result = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N,  kfactor);
				
				if (kkmag <= kkmax_inside)
				{
					shell_index = Get_shell_index(kkmag, shell_radius_array);
					
					local_result(shell_index) += 2*Multiplicity_factor_SCFT(l1, l2, l3, N) 
													* real( A(l1,l2,l3)*conj(B(l1,l2,l3)) );
				}									   
			}
		
}


//
//
void Shell_mult_all_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
)
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_all_SCFT(alias_switch, N, A, B, shell_radius_array, local_result,kfactor);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}				


//*********************************************************************************************


void Local_shell_mult_all_real_imag_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result_real, Array<DP,1> local_result_imag, 
	DP kfactor[]
)
{
	DP kkmag;
	DP factor;
	int shell_index;
	
	local_result_real = 0.0;
	local_result_imag = 0.0;

	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N,  kfactor);
				
				if (kkmag <= kkmax_inside)
				{
					shell_index = Get_shell_index(kkmag, shell_radius_array);
					
					factor = Multiplicity_factor_SCFT(l1, l2, l3, N);
		
					local_result_real(shell_index) += 2*factor * real(A(l1,l2,l3))
															   * real(B(l1,l2,l3)); 
					
					local_result_imag(shell_index) += 2*factor * imag(A(l1,l2,l3))
															   * imag(B(l1,l2,l3)); 
				}	
			}
	
}

//
//
void Shell_mult_all_real_imag_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result_real, Array<DP,1> result_imag, 
	DP kfactor[]
)
{
	static Array<DP,1> local_result_real(result_real.length());
	static Array<DP,1> local_result_imag(result_real.length());
	
	Local_shell_mult_all_real_imag_SCFT(alias_switch, N, A, B, shell_radius_array, 
										local_result_real, local_result_imag, kfactor);
	
	int data_size = result_real.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_result_real.data()), 
						reinterpret_cast<double*>(result_real.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
						
	MPI_Reduce(reinterpret_cast<double*>(local_result_imag.data()), 
						reinterpret_cast<double*>(result_imag.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);					

}


//*********************************************************************************************

void Local_shell_mult_all_SCFT
(	
	string alias_switch,
	int N[],  
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result, 
	DP kfactor[]
)
{
	DP kkmag;
	int shell_index;
	
	local_result = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N,  kfactor);
				
				if (kkmag <= kkmax_inside)
				{
					shell_index = Get_shell_index(kkmag, shell_radius_array);
					
					local_result(shell_index) += 2* Multiplicity_factor_SCFT(l1, l2, l3, N) 
												  * real( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
														+ Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
													    + Az(l1,l2,l3)*conj(Bz(l1,l2,l3)) ); 
				}									 
			}
		
}


//
//
void Shell_mult_all_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,   
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
)
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, shell_radius_array, 
								local_result, kfactor);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}				


//*********************************************************************************************


void Local_shell_mult_all_real_imag_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,   
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result_real, Array<DP,1> local_result_imag, 
	DP kfactor[]
)
{
	DP kkmag;
	DP factor;
	int shell_index;
	
	local_result_real = 0.0;
	local_result_imag = 0.0;

	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N,  kfactor);
				
				if (kkmag <= kkmax_inside)
				{
					shell_index = Get_shell_index(kkmag, shell_radius_array);
					factor = Multiplicity_factor_SCFT(l1, l2, l3, N);
		
					local_result_real(shell_index)  
									+= 2*factor *( real(Ax(l1,l2,l3))*real(Bx(l1,l2,l3)) 
												 + real(Ay(l1,l2,l3))*real(By(l1,l2,l3)) 
												 + real(Az(l1,l2,l3))*real(Bz(l1,l2,l3)) ); 
													
					local_result_imag(shell_index) 
									+= 2*factor *( imag(Ax(l1,l2,l3))*imag(Bx(l1,l2,l3)) 
												 + imag(Ay(l1,l2,l3))*imag(By(l1,l2,l3)) 
												 + imag(Az(l1,l2,l3))*imag(Bz(l1,l2,l3)) );
				}
			}
	
}

//
//
void Shell_mult_all_real_imag_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result_real, Array<DP,1> result_imag, 
	DP kfactor[]
)
{
	static Array<DP,1> local_result_real(result_real.length());
	static Array<DP,1> local_result_imag(result_real.length());
	
	Local_shell_mult_all_real_imag_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
						shell_radius_array, local_result_real, local_result_imag, kfactor);
	
	int data_size = result_real.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_result_real.data()), 
						reinterpret_cast<double*>(result_real.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
						
	MPI_Reduce(reinterpret_cast<double*>(local_result_imag.data()), 
						reinterpret_cast<double*>(result_imag.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);					

}


/**********************************************************************************************

			RING MULT

***********************************************************************************************/

void Local_ring_mult_all_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
)
{
	DP kkmag, theta;
	int shell_index, sector_index;
	
	local_result = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N,  kfactor);
				
				if ((kkmag <= kkmax_inside) && (kkmag > MYEPS))
				{
					theta = AnisKvect_polar_angle_SCFT(l1, l2, l3, N, kfactor);
					
					Compute_ring_index(kkmag, theta, shell_radius_array, sector_angle_array, 
														shell_index, sector_index);
					
					local_result(shell_index, sector_index)  
										+=  2 * Multiplicity_factor_SCFT(l1, l2, l3, N)
											  * real( A(l1,l2,l3)*conj(B(l1,l2,l3)) ); 
				}					  
			}
}


//
//
void Ring_mult_all_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_ring_mult_all_SCFT(alias_switch, N, A, B, shell_radius_array, sector_angle_array, 
								local_result, kfactor);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}				



//*********************************************************************************************

void Local_ring_mult_all_real_imag_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result_real, Array<DP,2> local_result_imag, 
	DP kfactor[]
)
{
	DP kkmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	local_result_real = 0.0;
	local_result_imag = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N,  kfactor);
				
				if ((kkmag <= kkmax_inside) && (kkmag > MYEPS))
				{
					theta = AnisKvect_polar_angle_SCFT(l1, l2, l3, N, kfactor);
					
					Compute_ring_index(kkmag, theta, shell_radius_array, sector_angle_array, 
														shell_index, sector_index);
					
					factor = Multiplicity_factor_SCFT(l1, l2, l3, N);
					
					local_result_real(shell_index, sector_index) 
												+= 2 * factor 
													 * ( real(A(l1,l2,l3))*real(B(l1,l2,l3)) );
												
					local_result_imag(shell_index, sector_index) 
												+= 2 * factor 
													 * ( imag(A(l1,l2,l3))*imag(B(l1,l2,l3)) );   
				}							
			}
}


//
//
void Ring_mult_all_real_imag_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result_real( (result_real.length())(0), (result_real.length())(1));
	static Array<DP,2> local_result_imag( (result_real.length())(0), (result_real.length())(1));
	
	Local_ring_mult_all_real_imag_SCFT(alias_switch, N, A, B, shell_radius_array, 
						sector_angle_array, local_result_real, local_result_imag, kfactor);
	
	int data_size = (result_real.length())(0) * (result_real.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result_real.data()), 
						reinterpret_cast<double*>(result_real.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
						
	MPI_Reduce(reinterpret_cast<double*>(local_result_imag.data()), 
						reinterpret_cast<double*>(result_imag.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 											
}		
		

//*********************************************************************************************

void Local_ring_mult_all_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
)
{
	DP kkmag, theta;
	int shell_index, sector_index;
	
	local_result = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N,  kfactor);
				
				if ((kkmag <= kkmax_inside) && (kkmag > MYEPS))
				{
					theta = AnisKvect_polar_angle_SCFT(l1, l2, l3, N, kfactor);
					
					Compute_ring_index(kkmag, theta, shell_radius_array, sector_angle_array, 
														shell_index, sector_index);
					
					local_result(shell_index, sector_index)  
										+=  2 * Multiplicity_factor_SCFT(l1, l2, l3, N)
											  * real( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
													+ Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
													+ Az(l1,l2,l3)*conj(Bz(l1,l2,l3)) ); 
				}									  
			}
		
}


//
//
void Ring_mult_all_SCFT
(
	string alias_switch,
	int N[],  
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_ring_mult_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, shell_radius_array, 
								sector_angle_array, local_result, kfactor);

	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}				


//*********************************************************************************************

void Local_ring_mult_all_real_imag_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result_real, Array<DP,2> local_result_imag, 
	DP kfactor[]
)
{
	DP kkmag, theta;
	DP factor;
	int shell_index, sector_index;
	
	local_result_real = 0.0;
	local_result_imag = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N,  kfactor);
				
				if ((kkmag <= kkmax_inside) && (kkmag > MYEPS))
				{
					theta = AnisKvect_polar_angle_SCFT(l1, l2, l3, N, kfactor);
					
					Compute_ring_index(kkmag, theta, shell_radius_array, sector_angle_array, 
														shell_index, sector_index);
					
					factor = Multiplicity_factor_SCFT(l1, l2, l3, N);
					
					local_result_real(shell_index, sector_index)    
									+= 2 * factor 
										 * ( real(Ax(l1,l2,l3))*real(Bx(l1,l2,l3)) 
										   + real(Ay(l1,l2,l3))*real(By(l1,l2,l3)) 
										   + real(Az(l1,l2,l3))*real(Bz(l1,l2,l3)) );
										
					local_result_imag(shell_index, sector_index)   
									+= 	2 * factor 
										  * ( imag(Ax(l1,l2,l3))*imag(Bx(l1,l2,l3)) 
										    + imag(Ay(l1,l2,l3))*imag(By(l1,l2,l3)) 
										    + imag(Az(l1,l2,l3))*imag(Bz(l1,l2,l3)) );  
				}					
			}
		
}


//
//
void Ring_mult_all_real_imag_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result_real( (result_real.length())(0), (result_real.length())(1));
	static Array<DP,2> local_result_imag( (result_real.length())(0), (result_real.length())(1));
	
	Local_ring_mult_all_real_imag_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										shell_radius_array, sector_angle_array, 
										local_result_real, local_result_imag, kfactor);
	
	int data_size = (result_real.length())(0) * (result_real.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result_real.data()), 
						reinterpret_cast<double*>(result_real.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
						
	MPI_Reduce(reinterpret_cast<double*>(local_result_imag.data()), 
						reinterpret_cast<double*>(result_imag.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 											
}				


/**********************************************************************************************
	
						CYLINDRICAL RING MULTIPLICATION OD ARRAYS
				
***********************************************************************************************/

void Local_cyl_ring_mult_all_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array,  
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
)
{

	DP kkpll, kkperp;
	int shell_index, slab_index;
	
	local_result = 0.0;
	
	int	kkperp_max = Anis_max_Krho_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkpll = AnisKpll_SCFT(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_SCFT(l1, l2, l3, N, kfactor);
				
				if (kkperp <= kkperp_max)
				{
					Compute_cylinder_ring_index(kkpll, kkperp, cylinder_shell_radius_array,
									cylinder_kpll_array, shell_index, slab_index);
					
					local_result(shell_index, slab_index) 
									+=  2 * Multiplicity_factor_SCFT(l1, l2, l3, N) 
										  * real( A(l1,l2,l3)*conj(B(l1,l2,l3)) );
				}										
			}
		
}

//
//

void Cyl_ring_mult_all_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array,  
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_cyl_ring_mult_all_SCFT(alias_switch, N, A, B, cylinder_shell_radius_array, 
								cylinder_kpll_array, local_result, kfactor);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}

//*********************************************************************************************


void Local_cyl_ring_mult_all_real_imag_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array, 
	Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> local_result_real, Array<DP,2> local_result_imag, 
	DP kfactor[]
)
{

	DP kkpll, kkperp;
	DP factor;
	int shell_index, slab_index;
	
	local_result_real = 0.0;
	local_result_imag = 0.0;
	
	int	kkperp_max = Anis_max_Krho_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkpll = AnisKpll_SCFT(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_SCFT(l1, l2, l3, N, kfactor);
				
				if (kkperp <= kkperp_max)
				{
					Compute_cylinder_ring_index(kkpll, kkperp, cylinder_shell_radius_array,
									cylinder_kpll_array, shell_index, slab_index);
					
					factor = Multiplicity_factor_SCFT(l1, l2, l3, N);
					
					local_result_real(shell_index, slab_index)  
								+=  2 * factor * ( real(A(l1,l2,l3))*real(B(l1,l2,l3)) );
													
					local_result_imag(shell_index, slab_index)  
								+=  2 * factor * ( imag(A(l1,l2,l3))*imag(B(l1,l2,l3)) );
				}				
			}
					
}

//
//
void Cyl_ring_mult_all_real_imag_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array, 
	Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result_real((result_real.length())(0), (result_real.length())(1));
	static Array<DP,2> local_result_imag((result_real.length())(0), (result_real.length())(1));
	
	Local_cyl_ring_mult_all_real_imag_SCFT(alias_switch, N, A, B, cylinder_shell_radius_array, 
						cylinder_kpll_array, local_result_real, local_result_imag, kfactor);
	
	int data_size = (result_real.length())(0) * (result_real.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result_real.data()), 
						reinterpret_cast<double*>(result_real.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
						
	MPI_Reduce(reinterpret_cast<double*>(local_result_imag.data()), 
						reinterpret_cast<double*>(result_imag.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 											
}	

//*********************************************************************************************

void Local_cyl_ring_mult_all_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
)
{

	DP kkpll, kkperp;
	int shell_index, slab_index;
	
	local_result = 0.0;
	
	int	kkperp_max = Anis_max_Krho_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkpll = AnisKpll_SCFT(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_SCFT(l1, l2, l3, N, kfactor);
				
				if (kkperp <= kkperp_max)
				{
					Compute_cylinder_ring_index(kkpll, kkperp, cylinder_shell_radius_array,
									cylinder_kpll_array, shell_index, slab_index);
					
					local_result(shell_index, slab_index) 
									+=  2 * Multiplicity_factor_SCFT(l1, l2, l3, N) 
										  * real( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
												+ Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
												+ Az(l1,l2,l3)*conj(Bz(l1,l2,l3)) );
				}									
			}

}


void Cyl_ring_mult_all_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_cyl_ring_mult_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
					cylinder_shell_radius_array, cylinder_kpll_array, local_result, kfactor);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}


//*********************************************************************************************

void Local_cyl_ring_mult_all_real_imag_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array, 
	Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> local_result_real, Array<DP,2> local_result_imag, 
	DP kfactor[]
)
{
	
	DP kkpll, kkperp, factor;
	int shell_index, slab_index;
	
	local_result_real = 0.0;
	local_result_imag = 0.0;
	
	int	kkperp_max = Anis_max_Krho_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkpll = AnisKpll_SCFT(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_SCFT(l1, l2, l3, N, kfactor);
				
				if (kkperp <= kkperp_max)
				{
					Compute_cylinder_ring_index(kkpll, kkperp, cylinder_shell_radius_array,
									cylinder_kpll_array, shell_index, slab_index);
					
					factor = Multiplicity_factor_SCFT(l1, l2, l3, N);
					
					local_result_real(shell_index, slab_index) 
											+=  2 * factor 
												  * ( real(Ax(l1,l2,l3))*real(Bx(l1,l2,l3)) 
													+ real(Ay(l1,l2,l3))*real(By(l1,l2,l3)) 
													+ real(Az(l1,l2,l3))*real(Bz(l1,l2,l3)) );
					
					local_result_imag(shell_index, slab_index) 
											+=  2 * factor 
												  * ( imag(Ax(l1,l2,l3))*imag(Bx(l1,l2,l3)) 
													+ imag(Ay(l1,l2,l3))*imag(By(l1,l2,l3)) 
													+ imag(Az(l1,l2,l3))*imag(Bz(l1,l2,l3)) );
				}									
			}
		
			
}


//
//

void Cyl_ring_mult_all_real_imag_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array, 
	Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result_real((result_real.length())(0), (result_real.length())(1));
	static Array<DP,2> local_result_imag((result_real.length())(0), (result_real.length())(1));
	
	Local_cyl_ring_mult_all_real_imag_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										cylinder_shell_radius_array,  cylinder_kpll_array, 
										local_result_real, local_result_imag, kfactor);
	
	int data_size = (result_real.length())(0) * (result_real.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result_real.data()), 
						reinterpret_cast<double*>(result_real.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 
						
	MPI_Reduce(reinterpret_cast<double*>(local_result_imag.data()), 
						reinterpret_cast<double*>(result_imag.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 											
}
/**********************************************************************************************

			(B0.k) Im[vectA. conj(vectB)]

***********************************************************************************************/


void Local_shell_mult_all_imagVW_B0_SCFT
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result, 
	DP kfactor[]
)
{
	DP kkmag;
	int shell_index;
	TinyVector<DP,3> kk;
	
	local_result = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N, kfactor);
				
				if (kkmag <= kkmax_inside)
				{
					shell_index = Get_shell_index(kkmag, shell_radius_array);
					
					Wavenumber_SCFT(l1, l2, l3, N, kfactor, kk);	
					
					local_result(shell_index) += 2* Multiplicity_factor_SCFT(l1, l2, l3, N) 
											* dot(B0, kk)
											* imag( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
												  + Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
												  + Az(l1,l2,l3)*conj(Bz(l1,l2,l3)) );
				}									
			}
					
}

//
//

void Shell_mult_all_imagVW_B0_SCFT
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
)
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_all_imagVW_B0_SCFT(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz, 
								shell_radius_array, local_result, kfactor);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}

//
//
void Local_ring_mult_all_imagVW_B0_SCFT
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> ring_shell_radius_array,  Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
)
{
	DP kkmag, theta;
	int shell_index, sector_index;
	TinyVector<DP,3> kk;
	
	local_result = 0.0;
		
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude_SCFT(l1, l2, l3, N, kfactor);
				
				if ((kkmag <= kkmax_inside) && (kkmag > MYEPS))
				{
					theta = AnisKvect_polar_angle_SCFT(l1, l2, l3, N, kfactor);
					
					Compute_ring_index(kkmag, theta, ring_shell_radius_array, 
								sector_angle_array, shell_index, sector_index);
					
					Wavenumber_SCFT(l1, l2, l3, N, kfactor, kk);
					
					local_result(shell_index, sector_index)   
							+=	 2* Multiplicity_factor_SCFT(l1, l2, l3, N) 
								  * dot(B0, kk) 
								  * imag( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
										+ Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
										+ Az(l1,l2,l3)*conj(Bz(l1,l2,l3)) );
				}									
			}
					
}

//
//
//

void Ring_mult_all_imagVW_B0_SCFT
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> ring_shell_radius_array, Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_ring_mult_all_imagVW_B0_SCFT(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz, 
						 ring_shell_radius_array, sector_angle_array, local_result, kfactor);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}


//
//
void Local_cyl_ring_mult_all_imagVW_B0_SCFT
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
)
{

	DP kkpll, kkperp;
	int shell_index, slab_index;
	TinyVector<DP,3> kk;
	
	local_result = 0.0;
	
	int	kkperp_max = Anis_max_Krho_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int l1=0; l1<local_N1; l1++)				
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=N[3]/2; l3++) 
			{
				kkpll = AnisKpll_SCFT(l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp_SCFT(l1, l2, l3, N, kfactor);
				
				if (kkperp <= kkperp_max)
				{
					Compute_cylinder_ring_index(kkpll, kkperp, cylinder_shell_radius_array,
									cylinder_kpll_array, shell_index, slab_index);

					Wavenumber_SCFT(l1, l2, l3, N, kfactor, kk);

					local_result(shell_index, slab_index) += 
												2* Multiplicity_factor_SCFT(l1, l2, l3, N) 
												 * dot(B0, kk)
												 * real( Ax(l1,l2,l3)*conj(Bx(l1,l2,l3)) 
													   + Ay(l1,l2,l3)*conj(By(l1,l2,l3)) 
													   + Az(l1,l2,l3)*conj(Bz(l1,l2,l3)) );
				}									
			}

}

//
//

void Cyl_ring_mult_all_imagVW_B0_SCFT
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{
	static Array<DP,2> local_result( (result.length())(0), (result.length())(1));
	
	Local_cyl_ring_mult_all_imagVW_B0_SCFT(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz, 
					cylinder_shell_radius_array, cylinder_kpll_array, local_result, kfactor);
	
	int data_size = (result.length())(0) * (result.length())(1);
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
						reinterpret_cast<double*>(result.data()),
						data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}


	//*********************************************************************************************


DP Local_shell_mult_vorticity_SCFT
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 )
{
	
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP kkmag;
	DP result = 0.0;
	
	for (int i1=0; i1<local_N1; i1++)				
		for (int i2=0; i2<N[2]; i2++) 
			for (int i3=0; i3<=N[3]/2; i3++) 
			{
				kkmag = Kmagnitude_SCFT(i1, i2, i3, N, kfactor);
				
				if ( (kkmag > inner_radius) && (kkmag <= outer_radius) ) 
				{
					if (N[2] > 1)
					{	
						Compute_Modal_vorticity_SCFT(i1, i2, i3, N, 
													 Bx, By, Bz, kfactor, vorticity);
						
							// factor 2 below to account for the absent complex conj modes. 
						result += 2 * Multiplicity_factor_SCFT(i1, i2, i3, N)
						* real( Ax(i1,i2,i3)*conj(vorticity(0)) 
							   + Ay(i1,i2,i3)*conj(vorticity(1)) 
							   + Az(i1,i2,i3)*conj(vorticity(2)) );
						
					}	
					
					else if (N[2] ==1)		// 2D
					{
						Compute_Modal_vorticity_y_component_SCFT(i1, 0, i3, N, 
																 Bx, By, Bz, kfactor, vort_y);
						
							// nlin2 contains U.grad omega
							// factor 2 to account for the absent complex conj modes.
						result += 2 * Multiplicity_factor_SCFT(i1, 0, i3, N)
						* real( Ay(i1,i2,i3)*conj(vort_y) );
						
					}	
					
				}	
			}
	
	return result;
}

	//
	//


DP Shell_mult_vorticity_SCFT
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 )
{
	DP local_total;
	local_total = Local_shell_mult_vorticity_SCFT(alias_switch, N, Ax, Ay, Az,
												  Bx, By, Bz, inner_radius, outer_radius, kfactor);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);	
	
	if (my_id == master_id) 
		return total;
	
	else
		return 0;	
}	

	//
	//

DP Local_shell_mult_vector_potential_SCFT
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 )
{
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP kkmag;
	DP result = 0.0;
	
	
	for (int i1=0; i1<local_N1; i1++)				
		for (int i2=0; i2<N[2]; i2++) 
			for (int i3=0; i3<=N[3]/2; i3++) 
			{
				kkmag = Kmagnitude_SCFT(i1, i2, i3, N, kfactor);
				
				if ( (kkmag > inner_radius) && (kkmag <= outer_radius) ) 
				{
					if (N[2] > 1)
					{	
						Compute_Modal_vorticity_SCFT(i1, i2, i3, N, 
													 Bx, By, Bz, kfactor, vorticity);
						
						result += 2 * Multiplicity_factor_SCFT(i1, i2, i3, N)/pow(kkmag,2)
						* real( Ax(i1,i2,i3)*conj(vorticity(0)) 
							   + Ay(i1,i2,i3)*conj(vorticity(1)) 
							   + Az(i1,i2,i3)*conj(vorticity(2)) );
							// factor 2 to account for the absent complex conj modes. 
					}	
					
					else if (N[2] ==1)		// 2D
					{
						Compute_Modal_vorticity_y_component_SCFT(i1, 0, i3, N, 
																 Bx, By, Bz, kfactor, vort_y);
						
						result += 2 * Multiplicity_factor_SCFT(i1, 0, i3, N)
						* real( Ay(i1,i2,i3)*conj(vort_y) )/pow(kkmag,2);
							// nlin2 contains U.grad a; a(k) = curl(b)/k^2
							// factor 2 to account for the absent complex conj modes. 
						
					}	
					
				}	
			}
	
	return result;
}


	//
	//
DP Shell_mult_vector_potential_SCFT
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 )
{
	DP local_total;
	local_total = Local_shell_mult_vector_potential_SCFT(alias_switch, N, Ax, Ay, Az,
														 Bx, By, Bz, inner_radius, outer_radius, kfactor);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);	
	
	if (my_id == master_id) 
		return total;
	
	else
		return 0;	
}	


	//
	//
void Local_shell_mult_vorticity_all_SCFT
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 )
{
	
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP kkmag;
	int shell_index;
	
	result = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int i1=0; i1<local_N1; i1++)				
		for (int i2=0; i2<N[2]; i2++) 
			for (int i3=0; i3<=N[3]/2; i3++) 
			{
				kkmag = Kmagnitude_SCFT(i1, i2, i3, N, kfactor);
				
				if (kkmag <= kkmax_inside)
				{
					shell_index = Get_shell_index(kkmag, shell_radius_array);
					
					if (N[2] > 1)
					{	
						Compute_Modal_vorticity_SCFT(i1, i2, i3, N, 
													 Bx, By, Bz, kfactor, vorticity);
						
							// factor 2 to account for the absent complex conj modes.
						result(shell_index) += 2 * Multiplicity_factor_SCFT(i1, i2, i3, N)
						* real( Ax(i1,i2,i3)*conj(vorticity(0)) 
							   + Ay(i1,i2,i3)*conj(vorticity(1)) 
							   + Az(i1,i2,i3)*conj(vorticity(2)) );
						
					}	
					
					else if (N[2] ==1)		// 2D
					{
						Compute_Modal_vorticity_y_component_SCFT(i1, 0, i3, N, 
																 Bx, By, Bz, kfactor, vort_y);
						
							// nlin2 contains U.grad omega
							// factor 2 to account for the absent complex conj modes. 
						result(shell_index) += 2 * Multiplicity_factor_SCFT(i1, 0, i3, N)
						* real( Ay(i1,i2,i3)*conj(vort_y) );
						
					}
				}							
			}
	
}

	//
	//
void Shell_mult_vorticity_all_SCFT
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 )
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_vorticity_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										shell_radius_array, local_result,kfactor);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
			   reinterpret_cast<double*>(result.data()),
			   data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}	


	//
	//
void Local_shell_mult_vector_potential_all_SCFT
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 )
{
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP kkmag;
	int shell_index;
	
	result = 0.0;
	
	int	kkmax_inside = Max_radius_inside_SCFT(alias_switch, N, kfactor);
	
	for (int i1=0; i1<local_N1; i1++)				
		for (int i2=0; i2<N[2]; i2++) 
			for (int i3=0; i3<=N[3]/2; i3++) 
			{
				kkmag = Kmagnitude_SCFT(i1, i2, i3, N, kfactor);
				
				if (kkmag <= kkmax_inside)
				{
					shell_index = Get_shell_index(kkmag, shell_radius_array);
					
					if (N[2] > 1)
					{	
						Compute_Modal_vorticity_SCFT(i1, i2, i3, N, 
													 Bx, By, Bz, kfactor, vorticity);
						
							// factor 2 to account for the absent complex conj modes.
						result(shell_index) += 2 * Multiplicity_factor_SCFT(i1, i2, i3, N)
						* real( Ax(i1,i2,i3)*conj(vorticity(0)) 
							   + Ay(i1,i2,i3)*conj(vorticity(1)) 
							   + Az(i1,i2,i3)*conj(vorticity(2)) )
						/pow(kkmag,2);
						
					}	
					
					else if (N[2] ==1)		// 2D
					{
						Compute_Modal_vorticity_y_component_SCFT(i1, 0, i3, N, 
																 Bx, By, Bz, kfactor, vort_y);
						
							// nlin2 contains U.grad omega
							// factor 2 to account for the absent complex conj modes. 
						result(shell_index) += 2 * Multiplicity_factor_SCFT(i1, 0, i3, N)
						* real( Ay(i1,i2,i3)*conj(vort_y) )
						/ pow(kkmag,2);
						
					}
				}							
			}
	
}


	//
	//
void Shell_mult_vector_potential_all_SCFT
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 )
{
	static Array<DP,1> local_result(result.length());
	
	Local_shell_mult_vector_potential_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
											   shell_radius_array, local_result,kfactor);
	
	int data_size = result.size();
	
	MPI_Reduce(reinterpret_cast<double*>(local_result.data()), 
			   reinterpret_cast<double*>(result.data()),
			   data_size, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD); 						
}	



//****************************  End of four_ET.cc   *******************************************




