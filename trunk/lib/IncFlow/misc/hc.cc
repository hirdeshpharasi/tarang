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

/*! \file  hc.cc
 * 
 * @brief Computes total cross helicity, cross helicity spectrum, force spectrum
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

/** @brief Computes cross helicity V.W/2.
 * 
 *  @param IncVF W
 *
 * @return \f$ Hc = \sum V.W/2 \f$ <BR>
 */
DP IncVF::Get_cross_helicity(IncVF& W)
{
	return (Get_total_energy(basis_type, alias_switch, N, *V1, *W.V1)  
			  + Get_total_energy(basis_type, alias_switch, N, *V2, *W.V2) 
			  + Get_total_energy(basis_type, alias_switch, N, *V3, *W.V3));
}



//*********************************************************************************************

/** @brief Computes cross helicity shell spectrum
 * 
 *  @param IncVF W
 *
 * @return \f$ *shell_ek_cross = real(V(k).W^*(k))/2 \f$ 
 */
void IncVF::Compute_cross_helicity_shell_spectrum(IncVF& W)
{

	CV_Compute_shell_spectrum(*W.V1, *W.V2, *W.V3, *shell_ek_cross);
}



//*********************************************************************************************
/** @brief Computes cross helicity ring spectrum
 * 
 *  @param IncVF W
 *
 * @return \f$ *shell_ek_cross = real(V(k).W^*(k))/2 \f$ 
 */
void IncVF::Compute_cross_helicity_ring_spectrum(IncVF& W)
{
	if (CV_anisotropic_ring_switch == 1)
		CV_Compute_ring_spectrum(*W.V1, *W.V2, *W.V3, *ring_ek_cross);	
}


//*********************************************************************************************
/** @brief Computes cross helicity ring spectrum for cylindrical rings
 * 
 *  @param IncVF W
 *
 * @return \f$ *shell_ek_cross = real(V(k).W^*(k))/2 \f$ 
 */
void IncVF::Compute_cross_helicity_cylinder_ring_spectrum(IncVF& W)
{
	if 	(CV_anisotropic_cylinder_switch == 1)
		CV_Compute_cylinder_ring_spectrum(*W.V1, *W.V2, *W.V3, *cylinder_ring_ek_cross);
}


//*********************************************************************************************
//*********************************************************************************************

/** @brief Computes  shell spectrum \f$ \Real(v_i(k) * T^*(k)/2) \f$ 
 * 
 *  @param IncSF T
 *
 * @return \f$ *shell_ek_cross = real(V_i(k) T^*(k))/2 \f$ 
 */
void IncVF::Compute_cross_vT_shell_spectrum(IncSF& T)
{
	
	Compute_shell_spectrum(basis_type, alias_switch, N, *V1, *T.F, 0, 
										*shell_ek_cross_V1T, kfactor);
										
	Compute_shell_spectrum(basis_type, alias_switch, N, *V2, *T.F, 0, 
										*shell_ek_cross_V2T, kfactor);
										
	Compute_shell_spectrum(basis_type, alias_switch, N, *V3, *T.F, 0, 
										*shell_ek_cross_V3T, kfactor);									
}


//*********************************************************************************************
/** @brief Computes ring spectrum \f$ \Real(v_i(k) * T^*(k)/2) \f$ 
 * 
 *  @param IncSF T
 *
 * @return \f$ *shell_ek_cross = real(V_i(k) T^*(k))/2 \f$ 
 */
void IncVF::Compute_cross_vT_ring_spectrum(IncSF& T)
{
	if (CV_anisotropic_ring_switch == 1)
	{
		Compute_ring_spectrum(basis_type, alias_switch, N, *V1, *T.F,
											*sector_angle_array_spectrum, 0,
											*ring_ek_cross_V1T, kfactor);
											
		Compute_ring_spectrum(basis_type, alias_switch, N, *V2, *T.F,
											*sector_angle_array_spectrum, 0,
											*ring_ek_cross_V2T, kfactor);
											
		Compute_ring_spectrum(basis_type, alias_switch, N, *V3, *T.F,	
											*sector_angle_array_spectrum, 0,
											*ring_ek_cross_V3T, kfactor);
	}																			
}


//*********************************************************************************************
/** @brief Computes ring spectrum \f$ \Real(v_i(k) * T^*(k)/2) \f$ for cylindrical geometry.
 * 
 *  @param IncSF T
 *
 * @return \f$ *shell_ek_cross = real(V_i(k) T^*(k))/2 \f$ 
 */
void IncVF::Compute_cross_vT_cylinder_ring_spectrum(IncSF& T)
{

	if 	(CV_anisotropic_cylinder_switch == 1)
	{
		Compute_cylinder_ring_spectrum(basis_type, alias_switch, N, *V1, *T.F,		
											*cylinder_kpll_array_spectrum, 0,
											*cylinder_ring_ek_cross_V1T, kfactor);
											
		Compute_cylinder_ring_spectrum(basis_type, alias_switch, N, *V2, *T.F, 
											*cylinder_kpll_array_spectrum, 0,
											*cylinder_ring_ek_cross_V2T, kfactor);
											
		Compute_cylinder_ring_spectrum(basis_type, alias_switch, N, *V3, *T.F, 
											*cylinder_kpll_array_spectrum, 0,
											*cylinder_ring_ek_cross_V3T, kfactor);	
	}																		
}


//***********************************  End of hc.cc *******************************************







