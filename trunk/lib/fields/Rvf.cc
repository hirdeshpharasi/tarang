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



/*! \file  Rvf.cc
 * 
 * @brief  Class declaration of Rvf, a Real Vector Field 
 *
 * @sa	Rvf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug  No known bugs
 */
 
 

#include "Rvf.h"

/*
#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif
 */

//*********************************************************************************************

RVF::RVF(int NN[], string string_switches[]) 
{ 

	RV_basis_type = string_switches[1];
	
    for (int i=1; i<=3; i++) 
      Nrv[i]=NN[i];
	  
	// If Transpose switch on, FT(*Vtr) = *V
	
#ifdef TRANSPOSE
	V1r = new Array<complx,3>(local_N2, Nrv[1], Nrv[3]/2+1);  
	V2r = new Array<complx,3>(local_N2, Nrv[1], Nrv[3]/2+1);  
	V3r = new Array<complx,3>(local_N2, Nrv[1], Nrv[3]/2+1);  
	
#else 
	V1r = new Array<complx,3>(local_N1, Nrv[2], Nrv[3]/2+1);			
	V2r = new Array<complx,3>(local_N1, Nrv[2], Nrv[3]/2+1);  
	V3r = new Array<complx,3>(local_N1, Nrv[2], Nrv[3]/2+1); 
#endif	

	(*V1r) = 0.0; 
	(*V2r) = 0.0; 
	(*V3r) = 0.0;
}

/**********************************************************************************************

				Forward transform: Inplace 
   
**********************************************************************************************/

void RVF::RV_Forward_transform(Array<complx,3> temp_r)
{
	Forward_transform_array(RV_basis_type, Nrv, *V1r, temp_r, 1);			// SFT
	
	Forward_transform_array(RV_basis_type, Nrv, *V2r, temp_r, 0);			// CFT
	
	Forward_transform_array(RV_basis_type, Nrv, *V3r, temp_r, 0);			// CFT	
}

/**********************************************************************************************

			Forward_transform(*Vir) = Vi 
				Vir unchanged.
   
**********************************************************************************************/


void RVF::RV_Forward_transform_transpose_order
(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp_r
)
{
	temp_r = *V1r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V1, 1);			// SFT
	
	temp_r = *V2r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V2, 0);			// CFT
	
	temp_r = *V3r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V3, 0);			// CFT	
}


/**********************************************************************************************

		Inplace Inverse Fourier transform.
   
**********************************************************************************************/

void RVF::RV_Inverse_transform(Array<complx,3> temp_r)
{
	Inverse_transform_array(RV_basis_type, Nrv, *V1r, temp_r, 1);			// ISFT
	
	Inverse_transform_array(RV_basis_type, Nrv, *V2r, temp_r, 0);			// ICFT
	
	Inverse_transform_array(RV_basis_type, Nrv, *V3r, temp_r, 0);			// ICFT
}


/**********************************************************************************************

			Inverse_transform(Vi) = *Vir 
				Keeping Vi unchanged....
				temp = N1 x N2 x N3
   
**********************************************************************************************/

void RVF::RV_Inverse_transform_transpose_order
(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp
)
{
	temp = V1;
	Inverse_transform_array_transpose_order(RV_basis_type, Nrv, temp, *V1r, 1);	// ISFT
	
	temp = V2;
	Inverse_transform_array_transpose_order(RV_basis_type, Nrv, temp, *V2r, 0);	// ICFT
	
	temp = V3;
	Inverse_transform_array_transpose_order(RV_basis_type, Nrv, temp, *V3r, 0);	// ICFT

}


/**********************************************************************************************

			 Forward  transform RSprod 
			  3D: *Vr[i] With SCFT -- 1,3 SFT;  2 CFT
   
**********************************************************************************************/


void RVF::RV_Forward_transform_RSprod(Array<complx,3> temp_r)
{	
	Forward_transform_array(RV_basis_type, Nrv, *V1r, temp_r, 1);			// SFT
	
	Forward_transform_array(RV_basis_type, Nrv, *V2r, temp_r, 0);			// CFT
	
	Forward_transform_array(RV_basis_type, Nrv, *V3r, temp_r, 1);			// SFT	
}
	

// *Vir unchanged..
void RVF::RV_Forward_transform_RSprod_transpose_order
(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp_r
)
{
	temp_r = *V1r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V1, 1);			// SFT
	
	temp_r = *V2r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V2, 0);			// CFT
	
	temp_r = *V3r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V3, 1);			// SFT	
}


/**********************************************************************************************

			Dealiase RVF
   
**********************************************************************************************/

void RVF::RV_Dealias()
{
	Dealias_array(RV_basis_type, Nrv, *V1r);
	
	Dealias_array(RV_basis_type, Nrv, *V2r);
	
	Dealias_array(RV_basis_type, Nrv, *V3r);
}


/**********************************************************************************************

			 Outputs Vi in real space
   
**********************************************************************************************/

void RVF::RV_Output(ofstream& fileout, Array<complx,3> temp_array)
{	

	Output_asreal(fileout, Nrv, *V1r, temp_array);
	
	Output_asreal(fileout, Nrv, *V2r, temp_array); 
	
	Output_asreal(fileout, Nrv, *V3r, temp_array);
}

/*
void RVF::RV_Output_hdf5(DataSet* dataset1, DataSet* dataset2, DataSet* dataset3, 
		 Array<complx,3> temp_array)
{	

	Output_asreal_hdf5(dataset1, Nrv, *V1r, temp_array);
	Output_asreal_hdf5(dataset2, Nrv, *V2r, temp_array);
	Output_asreal_hdf5(dataset3, Nrv, *V3r, temp_array);
	
	
}*/
//**************************  End of RVF class definitions ************************************






