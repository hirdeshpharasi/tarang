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


/*! \file  zero_mode_RB_slip.cc
 * 
 * @brief For RB convection only <BR>
 *			In 2D: set V(0,j) = theta(0,j) = 0; In 3D V(0,j, k) = theta(0,j, k) = 0;
 *
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 */

#include "../IncVF.h"
#include "../IncSF.h"

//**********************************************************************************************


/** @brief Set the modes to zerp for (0,ky,kz) in 3D. 
 *
 * @return In 3D: \f$ V_1(0,ky,kz) = 0 \f$  because sin(0) =0 [kx =0].
 * @return In 2D; \f$ V_2(0,ky,kz) = 0, V_3(0,ky,kz) = 0 \f$  by choice.
 */
void  IncVF::Zero_modes_RB_slip(IncSF& T)
{

	if (my_id == 0)			// i = 0 reside in master node
		for (int j=0; j<N[2]; j++)
			for (int k=0; k<=N[3]/2; k++)
				(*V1)(0,j,k) = (*T.F)(0,j,k) = 0.0;
		
}	

//
// Magnetoconvection RB slip
//

void  IncVF::Zero_modes_RB_slip(IncVF& W, IncSF& T)
{
	if (my_id == 0)
	{
		cout << "Implement zero mode " << endl;
		exit(1);
	}	
		
}				




//******************************** End of zero_mode_RB_slip.cc ********************************


