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

/*! \file  Output_field.cc
 * 
 * @brief  Output_field   
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"   


//*********************************************************************************************

void IncFluid::Output_field()
{
	if (my_id == master_id)	
		field_out_file << "%% Time = " << Tnow << endl; 
	
	if (output_vx_vy_switch == 0)
		CV_output(field_out_file, *VF_temp, output_field_format);
	
	else {
		Write_data_MPI(basis_type, field_out_file, N, *V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_out_file, N, *V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_file, N, *V3, output_field_format);
	}
	// *VF_temp is the temporary array useful in this operation
}
  
/*
 void IncFluid::Output_field_hdf5()
 {
 CV_output_hdf5(field_dataset1R, field_dataset1C, field_dataset2R,
 field_dataset2C, field_dataset3R, field_dataset3C,
 *VF_temp);		
 // *VF_temp is the temporary array useful in this operation
 }
 */

void IncFluid::Output_field(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_field_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_field_RB(T);
}

void IncFluid::Output_field_scalar(IncSF& T)
{
	if (my_id == master_id)	
		field_out_file << "%% Time = " << Tnow << endl; 
	
	if (output_vx_vy_switch == 0)
		CV_output(field_out_file, *VF_temp, output_field_format);
	
	else {
		Write_data_MPI(basis_type, field_out_file, N, *V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_out_file, N, *V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_file, N, *V3, output_field_format);
	}
	
	T.CS_output(field_out_file, *VF_temp, output_field_format);		
	// *VF_temp is the temporary array useful in this operation
}


void IncFluid::Output_field_RB(IncSF& T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_field();
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
		if (my_id == master_id)	
			field_out_file << "%% Time = " << Tnow << endl; 
		
		T.CS_output(field_out_file, *VF_temp, output_field_format);
	}
	
	else
		Output_field_scalar(T);
}

void IncFluid::Output_field(IncVF& W)
{
	if (my_id == master_id)	
		field_out_file << "%% Time = " << Tnow << endl; 
	
	if (output_vx_vy_switch == 0) {
		CV_output(field_out_file, *VF_temp, output_field_format);
		W.CV_output(field_out_file, *VF_temp, output_field_format);
	}
	
	else {
		Write_data_MPI(basis_type, field_out_file, N, *V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_out_file, N, *V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_file, N, *V3, output_field_format);
		
		Write_data_MPI(basis_type, field_out_file, N, *W.V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_out_file, N, *W.V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_file, N, *W.V3, output_field_format);
	}
	// *VF_temp is the temporary array useful in this operation
}


void IncFluid::Output_field(IncVF& W, IncSF& T)
{
	if (my_id == master_id)	
		field_out_file << "%% Time = " << Tnow << endl; 
	
	if (output_vx_vy_switch == 0) {
		CV_output(field_out_file, *VF_temp, output_field_format);
		W.CV_output(field_out_file, *VF_temp, output_field_format);
	}
	
	else {
		Write_data_MPI(basis_type, field_out_file, N, *V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_out_file, N, *V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_file, N, *V3, output_field_format);
		
		Write_data_MPI(basis_type, field_out_file, N, *W.V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_out_file, N, *W.V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_file, N, *W.V3, output_field_format);
	}
	
	T.CS_output(field_out_file, *VF_temp, output_field_format);	
	// *VF_temp is the temporary array useful in this operation
}


//*********************************************************************************************

//			IncFluid::Output_field_frequent()

//*********************************************************************************************

void IncFluid::Output_field_frequent()
{
	if (my_id == master_id)	
	{
		string filename = "/out/field_frequent_out.d";
		filename = data_dir_name+ filename;   
		field_frequent_out_file.open(filename.c_str(), ios::trunc);	
		
	//	field_frequent_out_file << "%% Time = " << Tnow << endl; 
	}
	
	if (output_vx_vy_switch == 0)
		CV_output(field_frequent_out_file, *VF_temp, output_field_format);
	
	else {
		Write_data_MPI(basis_type, field_frequent_out_file, N, *V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_frequent_out_file, N, *V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_frequent_out_file, N, *V3, output_field_format);
	}

	if (my_id == master_id)	
		field_frequent_out_file.close();
}
  

void IncFluid::Output_field_frequent(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_field_frequent_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_field_frequent_RB(T);
} 


void IncFluid::Output_field_frequent_scalar(IncSF& T)
{
	if (my_id == master_id)	
	{
		string filename = "/out/field_frequent_out.d";
		filename = data_dir_name+ filename;   
		field_frequent_out_file.open(filename.c_str(), ios::trunc);	
		
//		field_frequent_out_file << "%% Time = " << Tnow << endl; 
	}
	
	
	if (output_vx_vy_switch == 0)
		CV_output(field_frequent_out_file, *VF_temp, output_field_format);
	
	else {
		Write_data_MPI(basis_type, field_frequent_out_file, N, *V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_frequent_out_file, N, *V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_frequent_out_file, N, *V3, output_field_format);
	}
	
	T.CS_output(field_frequent_out_file, *VF_temp, output_field_format);	
	

	if (my_id == master_id)	
		field_frequent_out_file.close();
}
  
void IncFluid::Output_field_frequent_RB(IncSF& T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_field_frequent();
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
		if (my_id == master_id)	
		{
			string filename = "/out/field_frequent_out.d";
			filename = data_dir_name+ filename;   
			field_frequent_out_file.open(filename.c_str(), ios::trunc);	
		}
		
//		field_frequent_out_file << "%% Time = " << Tnow << endl; 
		
		T.CS_output(field_frequent_out_file, *VF_temp, output_field_format);
		
		if (my_id == master_id)	
			field_frequent_out_file.close();
	}
	
	else
		Output_field_frequent_scalar(T);
}



void IncFluid::Output_field_frequent(IncVF& W)
{
	if (my_id == master_id)	
	{
		string filename = "/out/field_frequent_out.d";
		filename = data_dir_name+ filename;   
		field_frequent_out_file.open(filename.c_str(), ios::trunc);	
		
//		field_frequent_out_file << "%% Time = " << Tnow << endl; 
	}
	
	if (output_vx_vy_switch == 0) {
		CV_output(field_out_file, *VF_temp, output_field_format);
		W.CV_output(field_frequent_out_file, *VF_temp, output_field_format);
	}
	
	else {
		Write_data_MPI(basis_type, field_frequent_out_file, N, *V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_frequent_out_file, N, *V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_frequent_out_file, N, *V3, output_field_format);
		
		Write_data_MPI(basis_type, field_frequent_out_file, N, *W.V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_frequent_out_file, N, *W.V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_frequent_out_file, N, *W.V3, output_field_format);
	}

	if (my_id == master_id)	
		field_frequent_out_file.close();
}

void IncFluid::Output_field_frequent(IncVF& W, IncSF& T)
{
	if (my_id == master_id)	
	{
		string filename = "/out/field_frequent_out.d";
		filename = data_dir_name+ filename;   
		field_frequent_out_file.open(filename.c_str(), ios::trunc);	
		
//		field_frequent_out_file << "%% Time = " << Tnow << endl; 
	} 
	
	if (output_vx_vy_switch == 0) {
		CV_output(field_frequent_out_file, *VF_temp, output_field_format);
		W.CV_output(field_frequent_out_file, *VF_temp, output_field_format);
	}
	
	else {
		Write_data_MPI(basis_type, field_frequent_out_file, N, *V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_frequent_out_file, N, *V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_frequent_out_file, N, *V3, output_field_format);
		
		Write_data_MPI(basis_type, field_frequent_out_file, N, *W.V1, *VF_temp, output_field_format);
		Write_data_MPI(basis_type, field_frequent_out_file, N, *W.V2, *VF_temp, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_frequent_out_file, N, *W.V3, output_field_format);
	}
	
	T.CS_output(field_frequent_out_file, *VF_temp, output_field_format);

	if (my_id == master_id)	
		field_frequent_out_file.close();
}


//*********************************************************************************************

//			Output_realfield()

//*********************************************************************************************
 

void IncFluid::Output_realfield()
{
	if (my_id == master_id)
		realfield_out_file << "%% Time = " << Tnow << endl; 	
	

#ifdef TRANSPOSE
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp); 
	RV_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format); 

#else	
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	RV_Inverse_transform(*VF_temp_r);

	RV_Output(realfield_out_file, *VF_temp, output_field_format);  
#endif	
}

//*********************************************************************************************
// Scalar

void IncFluid::Output_realfield(IncSF& T)
{	
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_realfield_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_realfield_RB(T);
}


void IncFluid::Output_realfield_scalar(IncSF& T)
{	
	if (my_id == master_id)		
		realfield_out_file << "%% Time = " << Tnow << endl; 	
	
#ifdef TRANSPOSE
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp); 
	T.RS_Inverse_transform_transpose_order(*T.F, *VF_temp); 
	
	RV_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format); 
	T.RS_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format);
	
#else	

	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	RV_Inverse_transform(*VF_temp_r);
	
	*T.Fr = *T.F;   
	T.RS_Inverse_transform(*VF_temp_r);
  
	RV_Output(realfield_out_file, *VF_temp, output_field_format); 
	T.RS_Output(realfield_out_file, *VF_temp, output_field_format);
#endif	
}

//	RB Convection	//

void IncFluid::Output_realfield_RB(IncSF& T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_realfield();
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
#ifdef TRANSPOSE		
		T.RS_Inverse_transform_transpose_order(*T.F, *VF_temp);
		T.RS_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format);
		
#else
		*T.Fr = *T.F;
		T.RS_Inverse_transform(*VF_temp_r);	
		T.RS_Output(realfield_out_file, *VF_temp, output_field_format);
#endif
		
	}
	
	else
		Output_realfield_scalar(T);
}

//*********************************************************************************************
// Vector
void IncFluid::Output_realfield(IncVF& W)
{
	if (my_id == master_id)		
		realfield_out_file << "%% Time = " << Tnow << endl; 	
	
#ifdef TRANSPOSE
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp); 
	W.RV_Inverse_transform_transpose_order(*W.V1, *W.V2, *W.V3, *W.VF_temp);
	
	RV_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format);
	W.RV_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format);
	
#else
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	RV_Inverse_transform(*VF_temp_r);
	
	*W.V1r = *W.V1;  
	*W.V2r = *W.V2;  
	*W.V3r = *W.V3; 
	W.RV_Inverse_transform(*VF_temp_r);
	
	RV_Output(realfield_out_file, *VF_temp, output_field_format); 
	W.RV_Output(realfield_out_file, *VF_temp, output_field_format);
#endif
}

//*********************************************************************************************
// Vectpr+scalar
void IncFluid::Output_realfield(IncVF& W, IncSF& T)
{
	if (my_id == master_id)		
		realfield_out_file << "%% Time = " << Tnow << endl; 	
	
#ifdef TRANSPOSE
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp); 
	W.RV_Inverse_transform_transpose_order(*W.V1, *W.V2, *W.V3, *W.VF_temp); 
	T.RS_Inverse_transform_transpose_order(*T.F, *VF_temp); 
	
	RV_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format);
	W.RV_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format);
	T.RS_Output_transpose_order(realfield_out_file, *VF_temp_r, output_field_format);
	
#else	
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	RV_Inverse_transform(*VF_temp_r);
	
	*W.V1r = *W.V1;  
	*W.V2r = *W.V2;  
	*W.V3r = *W.V3;
	W.RV_Inverse_transform(*VF_temp_r);
	
	*T.Fr = *T.F;   
	T.RS_Inverse_transform(*VF_temp_r);
	
	RV_Output(realfield_out_file, *VF_temp, output_field_format); 
	W.RV_Output(realfield_out_file, *VF_temp, output_field_format);
	T.RS_Output(realfield_out_file, *VF_temp, output_field_format);
#endif
}  



//*********************************************************************************************

//			IncFluid::Output_field_reduced()
//			IN ASCII FORMAT

//********************************************************************************************* 

void IncFluid::Output_field_reduced()
{

	if (my_id == master_id)	
		field_out_reduced_file << "%% Time = " << Tnow << endl;

	if (output_vx_vy_switch == 0)
		CV_output(field_out_reduced_file, N_out_reduced, output_field_format);
	
	else {
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *V1, output_field_format);
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *V2, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, 
								*V3, output_field_format);
	}
} 

//*********************************************************************************************
// scalar


void IncFluid::Output_field_reduced(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Output_field_reduced_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Output_field_reduced_RB(T);
}


void IncFluid::Output_field_reduced_scalar(IncSF& T)
{

	if (my_id == master_id)	
		field_out_reduced_file << "%% Time = " << Tnow << endl;

	if (output_vx_vy_switch == 0)
		CV_output(field_out_reduced_file, N_out_reduced, output_field_format);
	
	else {
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *V1, output_field_format);
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *V2, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, 
								*V3, output_field_format);
	}	
	
	T.CS_output(field_out_reduced_file, N_out_reduced, output_field_format);
}

//		RB-Convection		//

void IncFluid::Output_field_reduced_RB(IncSF& T)
{
	if (globalvar_Pr_switch == "PRZERO")
		Output_field_reduced();
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
		if (my_id == master_id)
			field_out_reduced_file << "%% Time = " << Tnow << endl; 
		
		T.CS_output(field_out_reduced_file, N_out_reduced, output_field_format);
	}
	
	else
		Output_field_reduced_scalar(T);
}

//*********************************************************************************************
//		MHD		
void IncFluid::Output_field_reduced(IncVF& W)
{
	if (my_id == master_id)	
		field_out_reduced_file << "%% Time = " << Tnow << endl;

	if (output_vx_vy_switch == 0) {
		CV_output(field_out_reduced_file, N_out_reduced, output_field_format);
		W.CV_output(field_out_reduced_file, N_out_reduced, output_field_format);
	}
	
	else {
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *V1, output_field_format);
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *V2, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, 
								*V3, output_field_format);
		
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *W.V1, output_field_format);
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *W.V2, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, 
								*W.V3, output_field_format);
	}

} 


//*********************************************************************************************
//		Convective MHD		

void IncFluid::Output_field_reduced(IncVF& W, IncSF& T)
{
	if (my_id == master_id)	
		field_out_reduced_file << "%% Time = " << Tnow << endl;

	if (output_vx_vy_switch == 0) {
		CV_output(field_out_reduced_file, N_out_reduced, output_field_format);
		W.CV_output(field_out_reduced_file, N_out_reduced, output_field_format);
	}
	
	else {
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *V1, output_field_format);
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *V2, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, 
								*V3, output_field_format);
		
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *W.V1, output_field_format);
		Write_data_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, *W.V2, output_field_format);
		
		Write_data_kz0plane_MPI(basis_type, field_out_reduced_file, N, N_out_reduced, 
								*W.V3, output_field_format);
	}

	T.CS_output(field_out_reduced_file, N_out_reduced, output_field_format);

}

/**********************************************************************************************

							IncFluid::Output_field_k()
			
***********************************************************************************************/


void IncFluid::Output_field_k()
{	
	
	static Array<DP,1> field_k_buf(MAXSIZE_K_R_BUFFER);
	field_k_buf  = 0;
	int	 tag = 123;
	
	int lx, ly, lz, ix;			// storage loc
	int kx, ky, kz;			// wavenumber
	DP Tuk;
	
	for (int i=1; i<= nos_output_waveno; i++)
	{
		kx = (*output_k_array)(i,1);
		ky = (*output_k_array)(i,2);			
		kz = (*output_k_array)(i,3);
		
		lx = Get_lx(basis_type, kx, N);	 
		ly = Get_ly3D(basis_type, ky, N); 
		lz = kz;	
		
		ix = Get_ix(basis_type, kx, N);
		
		if ((lx >= 0) && (lx < local_N1)) 
		{
			Tuk = Get_Tk(kx, ky, kz);
		
			field_k_buf(1) = real((*V1)(lx,ly,lz));		
			field_k_buf(2) = imag((*V1)(lx,ly,lz));
			field_k_buf(3) = real((*V2)(lx,ly,lz));		
			field_k_buf(4) = imag((*V2)(lx,ly,lz));
			field_k_buf(5) = real((*V3)(lx,ly,lz));		
			field_k_buf(6) = imag((*V3)(lx,ly,lz));
			
			field_k_buf(7) = Tuk;
		}
		
		// the point lies in the master node.
		if (ix < (N[1]/numprocs)) 
		{	
			if (my_id == master_id)
				field_k_out_file << Tnow << "  " << kx << " " << ky << " " << kz << "   " 
					<< field_k_buf(1) << " " << field_k_buf(2) << " " << field_k_buf(3) << " "
					<< field_k_buf(4) << " " << field_k_buf(5) << " " << field_k_buf(6) << " " 
					<< field_k_buf(7) << endl;
		}
		else
		{				
			if ((my_id != master_id)  && ((lx >= 0) && (lx < local_N1)))  // send to master
			{
				MPI_Send( reinterpret_cast<double*>(field_k_buf.data()), MAXSIZE_K_R_BUFFER, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );	
			}	
			
			else  if (my_id == master_id)   // Receive data from the source
			{
				MPI_Recv( reinterpret_cast<double*>(field_k_buf.data()), MAXSIZE_K_R_BUFFER, 
							MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);	
				
				field_k_out_file << Tnow << "  " << kx << " " << ky << " " << kz << "   " 
						<< field_k_buf(1) << " " << field_k_buf(2) << " " << field_k_buf(3) << " "
						<< field_k_buf(4) << " " << field_k_buf(5) << " " << field_k_buf(6) << " " 
						<< field_k_buf(7) << endl;
			}									
																				
			MPI_Barrier(MPI_COMM_WORLD);		// Sync procs
		}	
	}		// of for loop

}
				
					
//*********************************************************************************************
// SCALAR
//

void IncFluid::Output_field_k(IncSF& T)
{
	
	static Array<DP,1> field_k_buf(MAXSIZE_K_R_BUFFER);
	field_k_buf  = 0;
	int	 tag = 123;
	
	int lx, ly, lz, ix;			// storage loc
	int kx, ky, kz;			// wavenumber
	DP Tuk, TFk;
	
	for (int i=1; i<= nos_output_waveno; i++) 
	{
		kx = (*output_k_array)(i,1);
		ky = (*output_k_array)(i,2);			
		kz = (*output_k_array)(i,3);
		
		lx = Get_lx(basis_type, kx, N);	 
		ly = Get_ly3D(basis_type, ky, N); 
		lz = kz;
		
		ix = Get_ix(basis_type, kx, N);
		
		if ((lx >= 0) && (lx < local_N1)) 
		{
			Tuk = Get_Tk(kx, ky, kz);  
			TFk = T.Get_Tk(kx, ky, kz);
			
			field_k_buf(1) = real((*V1)(lx,ly,lz));		
			field_k_buf(2) = imag((*V1)(lx,ly,lz));
			field_k_buf(3) = real((*V2)(lx,ly,lz));		
			field_k_buf(4) = imag((*V2)(lx,ly,lz));
			field_k_buf(5) = real((*V3)(lx,ly,lz));		
			field_k_buf(6) = imag((*V3)(lx,ly,lz));
			
			field_k_buf(7) = real((*T.F)(lx,ly,lz));	
			field_k_buf(8) = imag((*T.F)(lx,ly,lz));
			
			field_k_buf(9) = Tuk;	
			field_k_buf(10) = TFk;
		}	
			
		// the point lies in the master node.
		if (ix < (N[1]/numprocs)) 
		{	
			if (my_id == master_id)
				field_k_out_file << Tnow << "  " << kx << " " << ky << " " << kz << "   " 
					<< field_k_buf(1) << " " << field_k_buf(2) << " " << field_k_buf(3) << " "
					<< field_k_buf(4) << " " << field_k_buf(5) << " " << field_k_buf(6) << " " 
					<< field_k_buf(7) << " " << field_k_buf(8) << " " << field_k_buf(9) << " "
					<< field_k_buf(10) << " "<< endl;
		}
		else
		{	
			if ((my_id != master_id)  && ((lx >= 0) && (lx < local_N1)))  // send to master 
			{
				MPI_Send( reinterpret_cast<double*>(field_k_buf.data()), MAXSIZE_K_R_BUFFER,
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );
			}	
		
			else  if (my_id == master_id)   // Receive data from the source
			{
				MPI_Recv( reinterpret_cast<double*>(field_k_buf.data()), MAXSIZE_K_R_BUFFER,
							MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status );
							
				field_k_out_file << Tnow << "  " << kx << " " << ky << " " << kz << "   " 
						<< field_k_buf(1) << " " << field_k_buf(2) << " " << field_k_buf(3) << " "
						<< field_k_buf(4) << " " << field_k_buf(5) << " " << field_k_buf(6) << " " 
						<< field_k_buf(7) << " " << field_k_buf(8) << " " << field_k_buf(9) << " "
						<< field_k_buf(10) << endl;
			}
			
			MPI_Barrier(MPI_COMM_WORLD);		// Sync procs
		}	
	}		// of for loop
}		

		
		
		
//*********************************************************************************************														
//		MHD			
//

void IncFluid::Output_field_k(IncVF& W)
{

	static Array<DP,1> field_k_buf(MAXSIZE_K_R_BUFFER);
	field_k_buf  = 0;
	int	 tag = 123;
	
	int lx, ly, lz, ix;			// storage loc
	int kx, ky, kz;			// wavenumber
	DP Tuk, TWk;
	
	for (int i=1; i<= nos_output_waveno; i++) 
	{
		kx = (*output_k_array)(i,1);
		ky = (*output_k_array)(i,2);			
		kz = (*output_k_array)(i,3);
		
		ix = Get_ix(basis_type, kx, N);
		
		lx = Get_lx(basis_type, kx, N);	 
		ly = Get_ly3D(basis_type, ky, N); 
		lz = kz;	 
		
		if ((lx >= 0) && (lx < local_N1)) 
		{
			Tuk = Get_Tk(kx, ky, kz);  
			TWk = W.Get_Tk(kx, ky, kz); 	  
			
			field_k_buf(1) = real((*V1)(lx,ly,lz));		
			field_k_buf(2) = imag((*V1)(lx,ly,lz));
			field_k_buf(3) = real((*V2)(lx,ly,lz));		
			field_k_buf(4) = imag((*V2)(lx,ly,lz));
			field_k_buf(5) = real((*V3)(lx,ly,lz));		
			field_k_buf(6) = imag((*V3)(lx,ly,lz));
			
			field_k_buf(7) = real((*W.V1)(lx,ly,lz));	
			field_k_buf(8) = imag((*W.V1)(lx,ly,lz));
			field_k_buf(9) = real((*W.V2)(lx,ly,lz));	
			field_k_buf(10) = imag((*W.V2)(lx,ly,lz));
			field_k_buf(11) = real((*W.V3)(lx,ly,lz));	
			field_k_buf(12) = imag((*W.V3)(lx,ly,lz));
			
			field_k_buf(13) = Tuk;	
			field_k_buf(14) = TWk;
		}	
			
		// the point lies in the master node.
		if (ix < (N[1]/numprocs)) 
		{	
			if (my_id == master_id)
				field_k_out_file << Tnow << "  " << kx << " " << ky << " " << kz << "   " 
					<< field_k_buf(1) << " " << field_k_buf(2) << " " << field_k_buf(3) << " "
					<< field_k_buf(4) << " " << field_k_buf(5) << " " << field_k_buf(6) << " " 
					<< field_k_buf(7) << " " << field_k_buf(8) << " " << field_k_buf(9) << " "
					<< field_k_buf(10) << " " << field_k_buf(11) << " " << field_k_buf(12) << " "
							<< field_k_buf(13) << " " << field_k_buf(14) << endl;
		}
		else
		{	
			if ((my_id != master_id)  && ((lx >= 0) && (lx < local_N1)))  // send to master 
			{
				MPI_Send( reinterpret_cast<double*>(field_k_buf.data()), MAXSIZE_K_R_BUFFER, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );
			}	
		
			else  if (my_id == master_id)   // Receive data from the source
			{
				MPI_Recv( reinterpret_cast<double*>(field_k_buf.data()), MAXSIZE_K_R_BUFFER, 
							MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status );
							
				field_k_out_file << Tnow << "  " << kx << " " << ky << " " << kz << "   " 
						<< field_k_buf(1) << " " << field_k_buf(2) << " " << field_k_buf(3) << " "
						<< field_k_buf(4) << " " << field_k_buf(5) << " " << field_k_buf(6) << " " 
						<< field_k_buf(7) << " " << field_k_buf(8) << " " << field_k_buf(9) << " "
						<< field_k_buf(10) << " " << field_k_buf(11) << " " << field_k_buf(12) << " "
						<< field_k_buf(13) << " " << field_k_buf(14) << endl;
			}
			
			MPI_Barrier(MPI_COMM_WORLD);		// Sync procs
		}	
	}		// of for loop
		
}



//*********************************************************************************************
// VF + Scalar
//

void IncFluid::Output_field_k(IncVF& W, IncSF& T)
{
	static Array<DP,1> field_k_buf(MAXSIZE_K_R_BUFFER);
	field_k_buf  = 0;
	int	 tag = 123;
	
	int lx, ly, lz, ix;			// storage loc
	int kx, ky, kz;			// wavenumber
	DP Tuk, TWk, TFk;
	
	
	for (int i=1; i<= nos_output_waveno; i++) 
	{
		kx = (*output_k_array)(i,1);
		ky = (*output_k_array)(i,2);			
		kz = (*output_k_array)(i,3);
		
		lx = Get_lx(basis_type, kx, N);	 
		ly = Get_ly3D(basis_type, ky, N); 
		lz = kz;	
		
		ix = Get_ix(basis_type, kx, N);
		
		if ((lx >= 0) && (lx < local_N1)) 
		{
		
			Tuk = Get_Tk(kx, ky, kz);
			TWk = W.Get_Tk(kx, ky, kz);
			TFk = T.Get_Tk(kx, ky, kz);
			
			field_k_buf(1) = real((*V1)(lx,ly,lz));		
			field_k_buf(2) = imag((*V1)(lx,ly,lz));
			field_k_buf(3) = real((*V2)(lx,ly,lz));		
			field_k_buf(4) = imag((*V2)(lx,ly,lz));
			field_k_buf(5) = real((*V3)(lx,ly,lz));		
			field_k_buf(6) = imag((*V3)(lx,ly,lz));
			
			field_k_buf(7) = real((*W.V1)(lx,ly,lz));	
			field_k_buf(8) = imag((*W.V1)(lx,ly,lz));
			field_k_buf(9) = real((*W.V2)(lx,ly,lz));	
			field_k_buf(10) = imag((*W.V2)(lx,ly,lz));
			field_k_buf(11) = real((*W.V3)(lx,ly,lz));	
			field_k_buf(12) = imag((*W.V3)(lx,ly,lz));
			
			field_k_buf(13) = real((*T.F)(lx,ly,lz));	
			field_k_buf(14) = imag((*T.F)(lx,ly,lz));
			
			field_k_buf(15) = Tuk;	
			field_k_buf(16) = TWk; 	
			field_k_buf(16) = TFk;
		}	
		
		// the point lies in the master node.
		if (ix < (N[1]/numprocs)) 
		{
			if (my_id == master_id)
				field_k_out_file << Tnow << "  " << kx << " " << ky << " " << kz << "   " 
					<< field_k_buf(1) << " " << field_k_buf(2) << " " << field_k_buf(3) << " "
					<< field_k_buf(4) << " " << field_k_buf(5) << " " << field_k_buf(6) << " " 
					<< field_k_buf(7) << " " << field_k_buf(8) << " " << field_k_buf(9) << " "
					<< field_k_buf(10) << " " << field_k_buf(11) << " " << field_k_buf(12) << " "
					<< field_k_buf(13) << " " << field_k_buf(14) << field_k_buf(15) << " " 
					<< field_k_buf(16) << " " << field_k_buf(17)  << endl;
		}	
		else
		{	
			if ((my_id != master_id)  && ((lx >= 0) && (lx < local_N1)))  // send to master 
			{
				MPI_Send( reinterpret_cast<double*>(field_k_buf.data()), MAXSIZE_K_R_BUFFER, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );	
			}
			else  if (my_id == master_id)   // Receive data from the source
			{
				MPI_Recv( reinterpret_cast<double*>(field_k_buf.data()), MAXSIZE_K_R_BUFFER, 
							MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status );
				field_k_out_file << Tnow << "  " << kx << " " << ky << " " << kz << "   " 
						<< field_k_buf(1) << " " << field_k_buf(2) << " " << field_k_buf(3) << " "
						<< field_k_buf(4) << " " << field_k_buf(5) << " " << field_k_buf(6) << " " 
						<< field_k_buf(7) << " " << field_k_buf(8) << " " << field_k_buf(9) << " "
						<< field_k_buf(10) << " " << field_k_buf(11) << " " << field_k_buf(12) << " "
						<< field_k_buf(13) << " " << field_k_buf(14) << field_k_buf(15) << " " 
						<< field_k_buf(16) << " " << field_k_buf(17)  << endl;
			}
			
			MPI_Barrier(MPI_COMM_WORLD);		// Sync procs
		}	
	}		// of for loop
}


/**********************************************************************************************

							IncFluid::Output_field_r()
			
***********************************************************************************************/


void IncFluid::Output_field_r()
{	

	static Array<DP,1> field_r_buf(MAXSIZE_K_R_BUFFER);
	field_r_buf  = 0;
	int	 tag = 123;
	
	int ix, iy, iz;			// array index
	int lx, ly, lz;			// storage loc in local processor
	int kx, ky, kz;
	int is_real;
			
	
#ifdef TRANSPOSE
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp); 

#else	
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	RV_Inverse_transform(*VF_temp_r);
#endif

		
	for (int i=1; i<= nos_output_position; i++)
	{
		ix = (*output_position_array)(i,1);
		iy = (*output_position_array)(i,2);	
		iz = (*output_position_array)(i,3);
		
		kx = ix;
		ky = iy;
		Real_to_fourier_index(iz, kz, is_real);
		
		lx = Get_lx(basis_type, kx, N);	 
		ly = Get_ly3D(basis_type, ky, N); 
		lz = kz;	
		
		if ((lx >= 0) && (lx < local_N1)) 
		{
		
			if (is_real == 1)
			{
				field_r_buf(1) =  real((*V1r)(lx,ly,lz)); 
				field_r_buf(2) =  real((*V2r)(lx,ly,lz));	
				field_r_buf(3) =  real((*V3r)(lx,ly,lz));	
			}
							 
			else
			{
				field_r_buf(1) =  imag((*V1r)(lx,ly,lz)); 
				field_r_buf(2) =  imag((*V2r)(lx,ly,lz));	
				field_r_buf(3) =  imag((*V3r)(lx,ly,lz));	
			}
		}	
		
		// Printing now...
		if (ix < (N[1]/numprocs))	// data in the master node itself
		{
			if (my_id == master_id)
				field_r_out_file << Tnow			<< "  " 
								<<	ix				<< " "
								<<	iy				<< " "
								<<  iz				<< " "
								<<  field_r_buf(1)	<< " " 
								<<  field_r_buf(2)	<< " " 
								<<  field_r_buf(3)	<< endl;
		}	
		else	
		{
			if ((my_id != master_id)  && ((lx >= 0) && (lx < local_N1)))  // send to master
			{	
				MPI_Send( reinterpret_cast<double*>(field_r_buf.data()), MAXSIZE_K_R_BUFFER, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );
			}
			else  if (my_id == master_id)   // master receives data and prints
			{
				MPI_Recv( reinterpret_cast<double*>(field_r_buf.data()), MAXSIZE_K_R_BUFFER, 
							MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
							
				field_r_out_file << Tnow			<< "  " 
									<<	ix				<< " "
									<<	iy				<< " "
									<<  iz				<< " "
									<<  field_r_buf(1)	<< " " 
									<<  field_r_buf(2)	<< " " 
									<<  field_r_buf(3)	<< endl;
			}	
			
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
}

//*********************************************************************************************

void IncFluid::Output_field_r(IncSF& T)
{	
	static Array<DP,1> field_r_buf(MAXSIZE_K_R_BUFFER);
	field_r_buf  = 0;
	int	 tag = 123;
	
	int ix, iy, iz;			// array index
	int lx, ly, lz;			// storage loc in local processor
	int kx, ky, kz;
	int is_real;
	
#ifdef TRANSPOSE
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp); 
	T.RS_Inverse_transform_transpose_order(*T.F, *VF_temp); 
	
#else	

	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	RV_Inverse_transform(*VF_temp_r);
	
	*T.Fr = *T.F;   
	T.RS_Inverse_transform(*VF_temp_r);
#endif
		
	for (int i=1; i<= nos_output_position; i++)
	{
		
		ix = (*output_position_array)(i,1);
		iy = (*output_position_array)(i,2);	
		iz = (*output_position_array)(i,3);
		
		kx = ix;
		ky = iy;	
		Real_to_fourier_index(iz, kz, is_real);
	
		lx = Get_lx(basis_type, kx, N);	 
		ly = Get_ly3D(basis_type, ky, N); 
		lz = kz;	
				
		if ((lx >= 0) && (lx < local_N1)) 
		{
		
			if (is_real == 1)
			{
				field_r_buf(1) =  real((*V1r)(lx,ly,lz)); 
				field_r_buf(2) =  real((*V2r)(lx,ly,lz));	
				field_r_buf(3) =  real((*V3r)(lx,ly,lz));
				
				field_r_buf(4) =  real((*T.Fr)(lx,ly,lz));			
			}
							 
			else
			{
				field_r_buf(1) =  imag((*V1r)(lx,ly,lz)); 
				field_r_buf(2) =  imag((*V2r)(lx,ly,lz));	
				field_r_buf(3) =  imag((*V3r)(lx,ly,lz));	
				
				field_r_buf(4) =  imag((*T.Fr)(lx,ly,lz));	
			}
		}
	
		// the point lies in the master node.
		if (ix < (N[1]/numprocs)) 
		{
			if (my_id == master_id)  
				field_r_out_file << Tnow	<< "  " 
						<<	ix				<< " "
						<<	iy				<< " "
						<<  iz				<< " "
						<<  field_r_buf(1)	<< " " 
						<<  field_r_buf(2)	<< " " 
						<<  field_r_buf(3)	<< " "
						<<  field_r_buf(4)	<< endl;
		}
		else 
		{	
			if ((my_id != master_id)  && ((lx >= 0) && (lx < local_N1)))  // send to master 
			{	
				MPI_Send( reinterpret_cast<double*>(field_r_buf.data()), MAXSIZE_K_R_BUFFER, 
					MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );
			}	
			
			else if (my_id == master_id)	// master receives the data and prints
			{
				MPI_Recv( reinterpret_cast<double*>(field_r_buf.data()), MAXSIZE_K_R_BUFFER, 
						 MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
				
				field_r_out_file << Tnow			<< "  " 
				<<	ix				<< " "
				<<	iy				<< " "
				<<  iz				<< " "
				<<  field_r_buf(1)	<< " " 
				<<  field_r_buf(2)	<< " " 
				<<  field_r_buf(3)	<< " "
				<<  field_r_buf(4)	<< endl;
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
		}
			
	}
			

}

//*********************************************************************************************


void IncFluid::Output_field_r(IncVF& W)
{	
	static Array<DP,1> field_r_buf(MAXSIZE_K_R_BUFFER);
	field_r_buf  = 0;
	int	 tag = 123;
	
	int ix, iy, iz;			// array index
	int lx, ly, lz;			// storage loc in local processor
	int kx, ky, kz;
	int is_real;
			
	
#ifdef TRANSPOSE
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp); 
	W.RV_Inverse_transform_transpose_order(*W.V1, *W.V2, *W.V3, *W.VF_temp); 
	
#else
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	RV_Inverse_transform(*VF_temp_r);
	
	*W.V1r = *W.V1;  
	*W.V2r = *W.V2;  
	*W.V3r = *W.V3; 
	W.RV_Inverse_transform(*VF_temp_r);
#endif
	
		
	for (int i=1; i<= nos_output_position; i++)
	{
		ix = (*output_position_array)(i,1);
		iy = (*output_position_array)(i,2);	
		iz = (*output_position_array)(i,3);
		
		kx = ix;
		ky = iy;
		Real_to_fourier_index(iz, kz, is_real);
		
		lx = Get_lx(basis_type, kx, N);	 
		ly = Get_ly3D(basis_type, ky, N); 
		lz = kz;	
		
		if ((lx >= 0) && (lx < local_N1)) 
		{		
			if (is_real == 1)
			{
				field_r_buf(1) =  real((*V1r)(lx,ly,lz)); 
				field_r_buf(2) =  real((*V2r)(lx,ly,lz));	
				field_r_buf(3) =  real((*V3r)(lx,ly,lz));
				
				field_r_buf(4) =  real((*W.V1r)(lx,ly,lz)); 
				field_r_buf(5) =  real((*W.V2r)(lx,ly,lz));	
				field_r_buf(6) =  real((*W.V3r)(lx,ly,lz));		
			}
							 
			else
			{
				field_r_buf(1) =  imag((*V1r)(lx,ly,lz)); 
				field_r_buf(2) =  imag((*V2r)(lx,ly,lz));	
				field_r_buf(3) =  imag((*V3r)(lx,ly,lz));	
				
				field_r_buf(4) =  imag((*W.V1r)(lx,ly,lz)); 
				field_r_buf(5) =  imag((*W.V2r)(lx,ly,lz));	
				field_r_buf(6) =  imag((*W.V3r)(lx,ly,lz));
			}
		}
		
		
		// the point lies in the master node.
		if (ix < (N[1]/numprocs))
		{	
			if (my_id == master_id)  
				field_r_out_file << Tnow			<< "  " 
								<<	ix				<< " "
								<<	iy				<< " "
								<<  iz				<< " "
								<<  field_r_buf(1)	<< " " 
								<<  field_r_buf(2)	<< " " 
								<<  field_r_buf(3)	<< " "
								<<  field_r_buf(4)	<< " " 
								<<  field_r_buf(5)	<< " " 
								<<  field_r_buf(6)	<< endl;
		}	
		else
		{	
			if ((my_id != master_id)  && ((lx >= 0) && (lx < local_N1)))  // send to master 
			{	
				MPI_Send( reinterpret_cast<double*>(field_r_buf.data()), MAXSIZE_K_R_BUFFER, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );
			}	
		
		
			else  if (my_id == master_id)   // master receives the data and prints
			{
				MPI_Recv( reinterpret_cast<double*>(field_r_buf.data()), MAXSIZE_K_R_BUFFER, 
							MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
							
				field_r_out_file << Tnow			<< "  " 
									<<	ix				<< " "
									<<	iy				<< " "
									<<  iz				<< " "
									<<  field_r_buf(1)	<< " " 
									<<  field_r_buf(2)	<< " " 
									<<  field_r_buf(3)	<< " "
									<<  field_r_buf(4)	<< " " 
									<<  field_r_buf(5)	<< " " 
									<<  field_r_buf(6)	<< endl;
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
	
		}

	}
}

//*********************************************************************************************


void IncFluid::Output_field_r(IncVF& W, IncSF& T)
{	
	static Array<DP,1> field_r_buf(MAXSIZE_K_R_BUFFER);
	field_r_buf  = 0;
	int	 tag = 123;
	
	int ix, iy, iz;			// array index
	int lx, ly, lz;			// storage loc in local processor
	int kx, ky, kz;
	int is_real;
	
#ifdef TRANSPOSE
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp); 
	W.RV_Inverse_transform_transpose_order(*W.V1, *W.V2, *W.V3, *W.VF_temp); 
	T.RS_Inverse_transform_transpose_order(*T.F, *VF_temp); 
	
#else	
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	RV_Inverse_transform(*VF_temp_r);
	
	*W.V1r = *W.V1;  
	*W.V2r = *W.V2;  
	*W.V3r = *W.V3;
	W.RV_Inverse_transform(*VF_temp_r);
	
	*T.Fr = *T.F;   
	T.RS_Inverse_transform(*VF_temp_r);
#endif

		
	for (int i=1; i<= nos_output_position; i++)
	{
		ix = (*output_position_array)(i,1);
		iy = (*output_position_array)(i,2);	
		iz = (*output_position_array)(i,3);
		
		kx = ix;
		ky = iy;
		Real_to_fourier_index(iz, kz, is_real);
		
		lx = Get_lx(basis_type, kx, N);	 
		ly = Get_ly3D(basis_type, ky, N); 
		lz = kz;	
		
		if ((lx >= 0) && (lx < local_N1)) 
		{
		
			if (is_real == 1)
			{
				field_r_buf(1) =  real((*V1r)(lx,ly,lz)); 
				field_r_buf(2) =  real((*V2r)(lx,ly,lz));	
				field_r_buf(3) =  real((*V3r)(lx,ly,lz));
				
				field_r_buf(4) =  real((*W.V1r)(lx,ly,lz)); 
				field_r_buf(5) =  real((*W.V2r)(lx,ly,lz));	
				field_r_buf(6) =  real((*W.V3r)(lx,ly,lz));
				
				field_r_buf(7) =  real((*T.Fr)(lx,ly,lz));			
			}
							 
			else
			{
				field_r_buf(1) =  imag((*V1r)(lx,ly,lz)); 
				field_r_buf(2) =  imag((*V2r)(lx,ly,lz));	
				field_r_buf(3) =  imag((*V3r)(lx,ly,lz));	
				
				field_r_buf(4) =  imag((*W.V1r)(lx,ly,lz)); 
				field_r_buf(5) =  imag((*W.V2r)(lx,ly,lz));	
				field_r_buf(6) =  imag((*W.V3r)(lx,ly,lz));
				
				field_r_buf(7) =  imag((*T.Fr)(lx,ly,lz));
			}
		}
		
		// the point lies in the master node.
		if (ix < (N[1]/numprocs)) 
		{
			if (my_id == master_id)
				field_r_out_file << Tnow			<< "  " 
								<<	ix				<< " "
								<<	iy				<< " "
								<<  iz				<< " "
								<<  field_r_buf(1)	<< " " 
								<<  field_r_buf(2)	<< " " 
								<<  field_r_buf(3)	<< " "
								<<  field_r_buf(4)	<< " " 
								<<  field_r_buf(5)	<< " " 
								<<  field_r_buf(6)	<< " "
								<<  field_r_buf(7)	<< endl;
		}
		else
		{
			if ((my_id != master_id)  && ((lx >= 0) && (lx < local_N1)))  // send to master 
			{	
				MPI_Send( reinterpret_cast<double*>(field_r_buf.data()), MAXSIZE_K_R_BUFFER, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );
			}
			else  if (my_id == master_id)  // master receives the data and prints
			{
				MPI_Recv( reinterpret_cast<double*>(field_r_buf.data()), MAXSIZE_K_R_BUFFER, 
							MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
							
				field_r_out_file << Tnow			<< "  " 
									<<	ix				<< " "
									<<	iy				<< " "
									<<  iz				<< " "
									<<  field_r_buf(1)	<< " " 
									<<  field_r_buf(2)	<< " " 
									<<  field_r_buf(3)	<< " "
									<<  field_r_buf(4)	<< " " 
									<<  field_r_buf(5)	<< " " 
									<<  field_r_buf(6)	<< " "
									<<  field_r_buf(7)	<< endl;
			}																									
			
			MPI_Barrier(MPI_COMM_WORLD);
		}	
	}

}


//******************************  End of output_field.cc  **************************************

