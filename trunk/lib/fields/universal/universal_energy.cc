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

/*! \file  universal_energy.cc
 * 
 * @sa	universal_energy.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug  No known bugs
 */
 
#include "universal_energy.h"



//*********************************************************************************************

// Get total fluctuating energy
DP Get_total_energy(string basis_type, string alias_switch, int N[], Array<complx,3> A)
{
	if (basis_type == "FOUR")
		return Get_total_energy_FOUR(alias_switch, N, A);
		
	else if (basis_type == "SCFT") 
		return Get_total_energy_SCFT(alias_switch, N, A);
		
	else
		return 0;		// Safe exit	
}


//

DP Get_total_energy
(
	string basis_type, 
	string alias_switch, 
	int N[],  
	Array<complx,3> A, 
	Array<complx,3> B
)
{

	if (basis_type == "FOUR")
		return  Get_total_energy_FOUR(alias_switch, N, A, B);
		
	else if (basis_type == "SCFT") 
		return  Get_total_energy_SCFT(alias_switch, N, A, B);
		
	else
		return 0;		// Safe exit	
}




//*********************************************************************************************

// Returns sum[ k^n |A(k)|^2/2]


DP Get_total_Sn
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> A, 
	DP n, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR")
		return  Get_total_Sn_FOUR(alias_switch, N, A, n, kfactor);
		
	else if (basis_type == "SCFT") 
		return  Get_total_Sn_SCFT(alias_switch, N, A, n, kfactor);
		
	else
		return 0;		// Safe exit	
}


//*********************************************************************************************


void Compute_shell_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> A, 
	DP n, 
	Array<DP,1> Sk, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR")  
		Compute_shell_spectrum_FOUR(alias_switch, N, A, n, Sk, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_shell_spectrum_SCFT(alias_switch, N, A, n, Sk, kfactor);
}

//
// A*B
//

void Compute_shell_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> B, 
	DP n, 
	Array<DP,1> Sk, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR")  
		Compute_shell_spectrum_FOUR(alias_switch, N, A, B, n, Sk, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_shell_spectrum_SCFT(alias_switch, N, A, B, n, Sk, kfactor);
}



// Computes total of helicity1, helicity2, dissipation_helicity1, dissipation_helicity2.
void Compute_total_helicity
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &total_helicity1, DP &total_helicity2, 
	DP &total_dissipation_H1, DP &total_dissipation_H2,
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Compute_total_helicity_FOUR(alias_switch, N, Ax, Ay, Az, 
										total_helicity1,  total_helicity2,
										total_dissipation_H1, total_dissipation_H2, kfactor);
		
	else if (basis_type == "SCFT")
		Compute_total_helicity_SCFT(alias_switch, N, Ax, Ay, Az,  
										total_helicity1,  total_helicity2,
										total_dissipation_H1, total_dissipation_H2, kfactor);
}	

// Computes shell spectrum \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$.
void Compute_shell_spectrum_helicity
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3,
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Compute_shell_spectrum_helicity_FOUR(alias_switch, N, Ax, Ay, Az, 
												H1k1, H1k2, H1k3, kfactor);
		
	else if (basis_type == "SCFT")
		Compute_shell_spectrum_helicity_SCFT(alias_switch, N, Ax, Ay, Az, 
												H1k1, H1k2, H1k3, kfactor);
}	

//*********************************************************************************************
//
// Entropy
//

DP Get_entropy
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
)
{
	if (basis_type == "FOUR") 
		return Get_entropy_FOUR(alias_switch, N, Ax, Ay, Az);
		
	else if (basis_type == "SCFT")
		return Get_entropy_SCFT(alias_switch, N, Ax, Ay, Az);
		
	else
		return 0;		// Safe exit	
}


DP Get_entropy_scalar(string basis_type, string alias_switch, int N[], Array<complx,3> A)
{
	if (basis_type == "FOUR") 
		return Get_entropy_scalar_FOUR(alias_switch, N, A);
		
	else if (basis_type == "SCFT")
		return Get_entropy_scalar_SCFT(alias_switch, N, A);
		
	else
		return 0;		// Safe exit	
}	

//*********************************************************************************************
//
//		Ring spectrum
// 

void Compute_ring_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR")  
		Compute_ring_spectrum_FOUR(alias_switch, N, Ax, Ay, Az, sector_angle_array, 
										n, S1k, S2k, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_ring_spectrum_SCFT(alias_switch, N, Ax, Ay, Az, sector_angle_array, 
										n, S1k, S2k, kfactor);
}

//
// A*B
//

void Compute_ring_spectrum
(
	string basis_type, 
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
	if (basis_type == "FOUR")  
		Compute_ring_spectrum_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										sector_angle_array, n, S1k, S2k, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_ring_spectrum_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										sector_angle_array, n, S1k, S2k, kfactor);
}


//
// scalar
//

void Compute_ring_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR")  
		Compute_ring_spectrum_FOUR(alias_switch, N, F, sector_angle_array, n, Sk, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_ring_spectrum_SCFT(alias_switch, N, F, sector_angle_array, n, Sk, kfactor);
}

//
//


void Compute_ring_spectrum
(
	string basis_type, 
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
	if (basis_type == "FOUR")  
		Compute_ring_spectrum_FOUR(alias_switch, N, F, G, sector_angle_array, n, Sk, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_ring_spectrum_SCFT(alias_switch, N, F, G, sector_angle_array, n, Sk, kfactor);
}

//
//

void Compute_ring_spectrum_helicity
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> H1k, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR")  
		Compute_ring_spectrum_helicity_FOUR(alias_switch, N, Ax, Ay, Az, 
										sector_angle_array, H1k, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_ring_spectrum_helicity_SCFT(alias_switch, N, Ax, Ay, Az, 
										sector_angle_array, H1k, kfactor);
}		



//*********************************************************************************************

void Compute_cylinder_ring_spectrum
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR")  
		Compute_cylinder_ring_spectrum_FOUR(alias_switch, N, Ax, Ay, Az, cylinder_kpll_array, 
										n, S1k, S2k, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_cylinder_ring_spectrum_SCFT(alias_switch, N, Ax, Ay, Az, cylinder_kpll_array, 
										n, S1k, S2k, kfactor);

}


//*********************************************************************************************

void Compute_cylinder_ring_spectrum
(
	string basis_type,
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

	if (basis_type == "FOUR")  
		Compute_cylinder_ring_spectrum_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										cylinder_kpll_array, n, S1k, S2k, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_cylinder_ring_spectrum_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										cylinder_kpll_array, n, S1k, S2k, kfactor);

}
	
//*********************************************************************************************

void Compute_cylinder_ring_spectrum
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> cylinder_kpll_array,  
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR")  
		Compute_cylinder_ring_spectrum_FOUR(alias_switch, N, F, cylinder_kpll_array, 
							n, Sk, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_cylinder_ring_spectrum_SCFT(alias_switch, N, F, cylinder_kpll_array, 
							n, Sk, kfactor);

}



//*********************************************************************************************

void Compute_cylinder_ring_spectrum
(
	string basis_type,
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

	if (basis_type == "FOUR")  
		Compute_cylinder_ring_spectrum_FOUR(alias_switch, N, F, G, cylinder_kpll_array, 
							n, Sk, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_cylinder_ring_spectrum_SCFT(alias_switch, N, F, G, cylinder_kpll_array, 
							n, Sk, kfactor);

}


//*********************************************************************************************

void Compute_cylinder_ring_spectrum_helicity
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> H1k, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR")  
		Compute_cylinder_ring_spectrum_helicity_FOUR(alias_switch, N, Ax, Ay, Az, 
										cylinder_kpll_array, H1k, kfactor);
		
	else if (basis_type == "SCFT")  
		Compute_cylinder_ring_spectrum_helicity_SCFT(alias_switch, N, Ax, Ay, Az, 
										cylinder_kpll_array, H1k, kfactor);

}

//*********************************************************************************************	


void Compute_imag_shell_spectrum_B0
(
	string basis_type,
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> Sk, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Compute_imag_shell_spectrum_B0_FOUR(alias_switch, N, B0, Ax, Ay, Ax, 
										Bx, By, Bz, Sk, kfactor);
		
	else if (basis_type == "SCFT")
		Compute_imag_shell_spectrum_B0_SCFT(alias_switch, N, B0, Ax, Ay, Ax, 
										Bx, By, Bz, Sk, kfactor);
}



//*********************************************************************************************	

void Compute_imag_ring_spectrum_B0
(
	string basis_type,
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> Sk, DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Compute_imag_ring_spectrum_B0_FOUR(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz,
										sector_angle_array, Sk, kfactor);

		
	else if (basis_type == "SCFT")
		Compute_imag_ring_spectrum_B0_SCFT(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz,
										sector_angle_array, Sk, kfactor);
}


//*********************************************************************************************	

void Compute_imag_cylinder_ring_spectrum_B0
(
	string basis_type,
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> Sk, DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Compute_imag_cylinder_ring_spectrum_B0_FOUR(alias_switch, N, B0, Ax, Ay, Az, 
										Bx, By, Bz, cylinder_kpll_array, Sk, kfactor);
		
	else if (basis_type == "SCFT")
		Compute_imag_cylinder_ring_spectrum_B0_SCFT(alias_switch, N, B0, Ax, Ay, Az, 
										Bx, By, Bz, cylinder_kpll_array, Sk, kfactor);
}


//********************************  End of universal_energy.cc   ******************************



