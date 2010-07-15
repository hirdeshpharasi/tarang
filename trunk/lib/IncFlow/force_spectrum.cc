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

/*! \file  force_spectrum.cc
 * 
 * @brief Computes force spectrum
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "IncVF.h"
#include "IncSF.h"

//*********************************************************************************************

/** @brief Computes Force spectrum (Force(K).V(K)^*), no factor of 1/2
 *
 * @return \f$ (*shell_spectrum_force_Vk) = \Re(F(K) \cdot V(K)^*)\f$ <BR>
 */
void IncVF::Compute_force_shell_spectrum()
{
	CVF::CV_Compute_shell_spectrum(*Force1, *Force2, *Force3, *shell_spectrum_force_Vk);	
	
	*shell_spectrum_force_Vk = 2*(*shell_spectrum_force_Vk);
}


//*********************************************************************************************

/** @brief Computes Force spectrum \f$ \Re(F^{\psi}(K).\psi(K)^*)\f$, no factor of 1/2
*			for IncSF
 *
 * @return \f$ (*shell_spectrum_force_Vk) = \Re(F^{\psi}(K).\psi(K)^*)\f$ <BR>
 */
void IncSF::Compute_force_shell_spectrum()							
{
	CSF::CS_Compute_shell_spectrum(*Force, *shell_spectrum_force_SFk);
	
	*shell_spectrum_force_SFk = 2*(*shell_spectrum_force_SFk);	
}



//*********************************************************************************************

void IncVF::Compute_force_ring_spectrum()
{

	CVF::CV_Compute_ring_spectrum(*Force1, *Force2, *Force3, *ring_spectrum_force_Vk);	
	
	*ring_spectrum_force_Vk = 2*(*ring_spectrum_force_Vk);
	
}

//*********************************************************************************************

void IncSF::Compute_force_ring_spectrum()					// Note IncSF function
{

	CSF::CS_Compute_ring_spectrum(*Force, *ring_spectrum_force_SFk);
	
	*ring_spectrum_force_SFk = 2*(*ring_spectrum_force_SFk);
	
}

//*********************************************************************************************

void IncVF::Compute_force_cylinder_ring_spectrum()
{

	CVF::CV_Compute_cylinder_ring_spectrum(*Force1, *Force2, *Force3, 
												*cylinder_ring_spectrum_force_Vk);
											
	*cylinder_ring_spectrum_force_Vk = 2*(*cylinder_ring_spectrum_force_Vk);
}

//*********************************************************************************************

void IncSF::Compute_force_cylinder_ring_spectrum()					// Note IncSF function
{

	CSF::CS_Compute_cylinder_ring_spectrum(*Force, *cylinder_ring_spectrum_force_SFk);
	
	*cylinder_ring_spectrum_force_SFk = 2*(*cylinder_ring_spectrum_force_SFk);
		
}



//***********************************  End of force_spectrum.cc *******************************





