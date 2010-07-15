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


/*! \file  Rsf.cc
 * 
 * @brief  Class declaration of Rsf, a Real Scalar Field 
 *
 * @sa	Rsf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug No known bugs
 */
 

#include "Rsf.h"


//*********************************************************************************************

//  Class constructor; Allocation for F, Ek etc. Initialize the arrays.

RSF::RSF(int NN[], string string_switches[]) 
{

	RS_basis_type = string_switches[1];
	
    for (int i=1; i<=3; i++) 
		Nrs[i] = NN[i];

#ifdef TRANSPOSE
	Fr = new Array<complx,3>(local_N2, Nrs[1], Nrs[3]/2+1); 
#else    
    Fr = new Array<complx,3>(local_N1, Nrs[2], Nrs[3]/2+1);  
#endif	

    *Fr = 0.0;   // initialize
	
}

/**********************************************************************************************

					Inplace Forward Fourier transform
   
**********************************************************************************************/

void RSF::RS_Forward_transform(Array<complx,3> temp_r)
{
	Forward_transform_array(RS_basis_type, Nrs, *Fr, temp_r, 1);			// SFT
}



/**********************************************************************************************

					Forward_transform_transpose_order(*Fr) = F 
						*Fr unchanged
   
**********************************************************************************************/

void RSF::RS_Forward_transform_transpose_order(Array<complx,3> F, Array<complx,3> temp_r)
{
	temp_r = *Fr;
	
	Forward_transform_array_transpose_order(RS_basis_type, Nrs, temp_r, F, 1);		// SFT
}



/**********************************************************************************************

					Inplace Inverse Fourier transform
   
**********************************************************************************************/

void RSF::RS_Inverse_transform(Array<complx,3> temp_r)
{
	Inverse_transform_array(RS_basis_type, Nrs, *Fr, temp_r, 1);			// ISFT
}

/**********************************************************************************************

					Inverse_transform_transpose_order(F) = *Fr 
						F unchanged
   
**********************************************************************************************/

void RSF::RS_Inverse_transform_transpose_order(Array<complx,3> F, Array<complx,3> temp)
{
	temp = F;
	
	Inverse_transform_array_transpose_order(RS_basis_type, Nrs, temp, *Fr, 1);		// ISFT
}


/**********************************************************************************************

				Dealias RSF
   
**********************************************************************************************/

void RSF::RS_Dealias()
{
	Dealias_array(RS_basis_type, Nrs, *Fr);
}




/**********************************************************************************************

    Outputs F in real space
   
**********************************************************************************************/

void RSF::RS_Output(ofstream& fileout, Array<complx,3> temp_array)
{
  Output_asreal(fileout, Nrs, *Fr, temp_array);
}


//************************ END of RSF class Definitions ***************************************



