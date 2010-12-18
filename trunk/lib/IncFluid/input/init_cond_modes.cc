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

/*! \file  init_cond_modes.cc
 * 
 * @brief Input several modes as initial condition
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 */


#include "../IncFluid.h"

extern Uniform<DP> SPECrand;

//*********************************************************************************************


/** @brief Read wavenos k and V(k) from file field_in_file. 
 * 
 *  @param k_dim	Dimension of k vector, e.g., 3 in 3D.
 *	@param V_n		No of components of vector to be read.
 *	@param num_field_waveno		Number of k vector read from the file.
 */ 	
void IncFluid::Read_waveno_field(int k_dim, int V_n, int& num_field_waveno)
{
	if (my_id == master_id) 
	{
		int i, ki;
		complx Vi;
		
		cout << "******** Reading field wavenumbers and Vx  *********" << endl;		
		
		int iter = 0;
		while (! field_in_file.eof()) 
		{	
			field_in_file >> ki;
			  
			if (ki < N[1])  
			{								// Stopper.. 
				iter++;	
				(*field_modes_k_array)(iter,1) = ki;
				
				for (i=2; i<=k_dim; i++) 
				{
					field_in_file >> ki;	
					(*field_modes_k_array)(iter,i) = ki;  
				}
		
				for (i=1; i<= V_n; i++) 
				{
					field_in_file >> Vi;	
					(*field_modes_V_array)(iter,i) = Vi; 
				}
			}
			else
				break;	
		}
			
		num_field_waveno = iter;
		cout << "******* End of Reading field wavenumbers and Vx, Vy **********" << endl;
		cout  << "Reading of field configurations ended successfully" << endl;
	}	
}	
	
//*********************************************************************************************

/** @brief Read V(k) for some modes as initial condition.  
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ V_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Vy computed using \f$ V_y = -D_x V_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * Complex conjugate modes are automatically added.
 */	
void  IncFluid::Init_cond_modes_SIMPLE()
{
	
	Input_prefix(field_in_file);
	
	(*field_modes_k_array) = 0;  
	(*field_modes_V_array) = 0.0;

	int  space_dim = 3;
	int  no_of_vectors_read = 2;
	
	if (my_id == master_id)	
		Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
	
	MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
	int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
	MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
	data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
	MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
									
	
	int kx, ky, kz;	
	complex<DP> Vx, Vy, Vz;
	
	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0;
	
	for (int i = 1; i <= num_field_waveno; i++) 
	{
		kx = (*field_modes_k_array)(i,1); 
		ky = (*field_modes_k_array)(i,2); 
		kz = (*field_modes_k_array)(i,3);
			
		Vx = (*field_modes_V_array)(i,1); 
		Vy = (*field_modes_V_array)(i,2); 
		
		Last_component(kx, ky, kz, Vx, Vy, Vz);
		
		Assign_field_comp_conj(kx, ky, kz, Vx, Vy, Vz);
		// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
		
		if (my_id == master_id)
			cout	<< "k, V : (" << kx << " " << ky  << " " <<kz << ") " 
					<<	Vx	<<	" "		<< Vy	<<	" "		<<	Vz	<<  endl;  
	}
	
	if (my_id == master_id)
		cout << "**********************************************" << endl;

}


//*********************************************************************************************


/** @brief Read V(k) and T.F(k) for some modes as initial condition.
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ V_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Vy computed using \f$ V_y = -D_x V_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * @note Complex conjugate modes are automatically added.
 *
 * @note Scalar field is added simply.
 */	 

void  IncFluid::Init_cond_modes_SIMPLE(IncSF& T)
{
	
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Init_cond_modes_SIMPLE_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Init_cond_modes_SIMPLE_RB(T);
}

void  IncFluid::Init_cond_modes_SIMPLE_scalar(IncSF& T)
{
		
	Input_prefix(field_in_file);
	
	(*field_modes_k_array) = 0;  
	(*field_modes_V_array) = 0.0;
	
	int  space_dim = 3;
	int  no_of_vectors_read = 3;
	
	if (my_id == master_id)	
		Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
	
	
	MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
	int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
	MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
	data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
	MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
	
	int kx, ky, kz;
	complex<DP> Vx, Vy, Vz, G;
	
	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0; 
	(*T.F) = 0.0;
	
	for (int i = 1; i <= num_field_waveno; i++) 
	{
		kx = (*field_modes_k_array)(i,1); 
		ky = (*field_modes_k_array)(i,2); 
		kz = (*field_modes_k_array)(i,3);
			
		Vx = (*field_modes_V_array)(i,1); 
		Vy = (*field_modes_V_array)(i,2); 
		G = (*field_modes_V_array)(i,3); 
		
		Last_component(kx, ky, kz, Vx, Vy, Vz);
		
		Assign_field_comp_conj(kx, ky, kz, Vx, Vy, Vz);
		T.Assign_field_comp_conj(kx, ky, kz, G);
		
		// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
		if (my_id == master_id)
			cout	<< "k, V, G : (" << kx << " " << ky  << " " <<kz << ") " 
					<< Vx  << " "  << Vy << " " << Vz << " "  
					<< G << endl;  
	}
	
	if (my_id == master_id)
		cout << "*******************************************" << endl;
	
}

// RBC
void  IncFluid::Init_cond_modes_SIMPLE_RB(IncSF& T)
{
	
	if (globalvar_Pr_switch == "PRZERO") 
	{
		Init_cond_modes_SIMPLE();		
		
		*T.F = *V1; 
		Array_divide_ksqr(basis_type, N, *T.F, kfactor);		
	}
	
	
	
	else if (globalvar_Pr_switch == "PRINFTY") 
	{
		
		Input_prefix(field_in_file);
		
		(*field_modes_k_array) = 0;  
		(*field_modes_V_array) = 0.0;
		
		int  space_dim = 3;
		int  no_of_vectors_read = 1;
		
		if (my_id == master_id)	
			Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
		
		
		MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
				  MPI_INT, master_id, MPI_COMM_WORLD); 
		
		data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
				  MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
		
		int kx, ky, kz;
		complex<DP> G;
		
		(*V1) = 0.0; 
		(*V2) = 0.0; 
		(*V3) = 0.0; 
		(*T.F) = 0.0;
		
		for (int i = 1; i <= num_field_waveno; i++) 
		{
			kx = (*field_modes_k_array)(i,1); 
			ky = (*field_modes_k_array)(i,2); 
			kz = (*field_modes_k_array)(i,3);
			
			G = (*field_modes_V_array)(i,1); 
			
			T.Assign_field_comp_conj(kx, ky, kz, G);
			
			// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
			if (my_id == master_id)
				cout	<< "k, G : (" << kx << " " << ky  << " " <<kz << ") " 
						<< G << endl;  
		}
		
		if (my_id == master_id)
			cout << "*******************************************" << endl;
		
		Init_cond_Prinfty(T);
	}
	
	else
		Init_cond_modes_SIMPLE_scalar(T);
	
	Zero_modes_RB_slip(T);
}	
//*********************************************************************************************

/** @brief Read V(k) and W(k) for some modes as initial condition.
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ V_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Vy computed using \f$ V_y = -D_x V_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * @note Complex conjugate modes are automatically added.
 *
 * @note Vector field W is added in a similar manner.
 */	
void  IncFluid::Init_cond_modes_SIMPLE(IncVF& W)
{

	Input_prefix(field_in_file);
	
	(*field_modes_k_array) = 0;  
	(*field_modes_V_array) = 0.0;

	int  space_dim = 3;
	int  no_of_vectors_read = 4;
	
	if (my_id == master_id)	
		Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
	
	
	MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
	int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
	MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
	data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
	MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
	
	int kx, ky, kz;
	complex<DP> Vx, Vy, Vz, Wx, Wy, Wz;
	
	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0;
	
	(*W.V1) = 0.0; 
	(*W.V2) = 0.0;   
	(*W.V3) = 0.0; 
	
	for (int i = 1; i <= num_field_waveno; i++) 
	{
		kx = (*field_modes_k_array)(i,1); 
		ky = (*field_modes_k_array)(i,2); 
		kz = (*field_modes_k_array)(i,3);
			
		Vx = (*field_modes_V_array)(i,1); 
		Vy = (*field_modes_V_array)(i,2); 
		 
		Wx = (*field_modes_V_array)(i,3); 
		Wy = (*field_modes_V_array)(i,4);
		
		Last_component(kx, ky, kz, Vx, Vy, Vz);
		Last_component(kx, ky, kz, Wx, Wy, Wz);
		
		Assign_field_comp_conj(kx, ky, kz, Vx, Vy, Vz);
		W.Assign_field_comp_conj(kx, ky, kz, Wx, Wy, Wz);
		// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
		
		if (my_id == master_id)
			cout	<< "k, V, W : (" << kx << " " << ky  << " " <<kz << ") " 
					<< Vx << " "	<< Vy << " "  << Vz << " " 
					<< Wx << " "	<< Wy << " "  << Wz << endl;
	}
	
	if (my_id == master_id)
		cout << "******************************************"  << endl;
	
}


//*********************************************************************************************

/** @brief Read V(k), W(k), and T.F(k) for some modes as initial condition.
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ V_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Vy computed using \f$ V_y = -D_x V_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * @note Complex conjugate modes are automatically added.
 *
 * @note Vector field W and scalar field T are added in a similar manner.
 */
void  IncFluid::Init_cond_modes_SIMPLE(IncVF& W, IncSF& T)
{

	Input_prefix(field_in_file);
	
	(*field_modes_k_array) = 0;  
	(*field_modes_V_array) = 0.0;
	
	int  space_dim = 3;
	int  no_of_vectors_read = 5;
	
	if (my_id == master_id)	
		Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
	
	MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
	int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
	MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
	data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
	MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
									
									
	int kx, ky, kz;	
	complex<DP> Vx, Vy, Vz, Wx, Wy, Wz, G;
	
	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0;
	
	(*W.V1) = 0.0; 
	(*W.V2) = 0.0;   
	(*W.V3) = 0.0; 
	
	(*T.F) = 0.0;
	
	for (int i = 1; i <= num_field_waveno; i++) 
	{
		kx = (*field_modes_k_array)(i,1); 
		ky = (*field_modes_k_array)(i,2); 
		kz = (*field_modes_k_array)(i,3);
			
		Vx = (*field_modes_V_array)(i,1); 
		Vy = (*field_modes_V_array)(i,2); 
		 
		Wx = (*field_modes_V_array)(i,3); 
		Wy = (*field_modes_V_array)(i,4); 
		
		G = (*field_modes_V_array)(i,5);
		
		Last_component(kx, ky, kz, Vx, Vy, Vz);
		Last_component(kx, ky, kz, Wx, Wy, Wz);
		
		Assign_field_comp_conj(kx, ky, kz, Vx, Vy, Vz);
		W.Assign_field_comp_conj(kx, ky, kz, Wx, Wy, Wz);
		T.Assign_field_comp_conj(kx, ky, kz, G);
		// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
		
		if (my_id == master_id)
			cout	<< "k, V, G : (" << kx << " " << ky  << " " <<kz << ") " 
					<< Vx << " " << Vy << " " << Vz << " " 
					<< Wx << " " << Wy << " " << Wz << " " 
					<< G << endl;
	}
	
	if (my_id == master_id)
		cout << "**********************************************" << endl;
	
}


//*********************************************************************************************
//*********************************************************************************************

/** @brief Read V(k) for some modes as initial condition.
 *				We provide V[1] and Omega[1] (vorticity along 1),
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), \Omega_x(\vec{k}), \f$ 
 *		In typical case Vz is computed using incompressible cond. <BR>
 *		kz = 0 cases are handled as a special case. See fn IncVF::Compute_VyVz(..).
 *
 * @sa void IncVF::Compute_VyVz(int kx, int ky, int kz, complx& Vx, complx Omega, 
 *		complx &Vy, complx &Vz)
 *
 * @note Complex conjugate modes are automatically added.
 */
void  IncFluid::Init_cond_modes_VORTICITY()
{
	
	Input_prefix(field_in_file);
	
	(*field_modes_k_array) = 0;  
	(*field_modes_V_array) = 0.0;

	int  space_dim = 3;
	int  no_of_vectors_read = 2;
	
	if (my_id == master_id)	
		Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
	
	MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
	int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
	MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
	data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
	MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
									
									
	
	int kx, ky, kz;	
	complex<DP> Vx, Vy, Vz, Omega;
	
	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0;
	
	for (int i = 1; i <= num_field_waveno; i++) 
	{
		kx = (*field_modes_k_array)(i,1); 
		ky = (*field_modes_k_array)(i,2); 
		kz = (*field_modes_k_array)(i,3);	
		
		Vx = (*field_modes_V_array)(i,1); 
		Omega = (*field_modes_V_array)(i,2); 
		
		Compute_VyVz(kx, ky, kz, Vx, Omega, Vy, Vz); 
		Assign_field_comp_conj(kx, ky, kz, Vx, Vy, Vz);
		// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
		
		if (my_id == master_id)
			cout	<< "k, V : (" << kx << " " << ky  << " " <<kz << ") " 
				<<	Vx << " " << Vy << " " << Vz <<  endl;  
	}

	if (my_id == master_id)
		cout << "**********************************************" << endl;

}

//*********************************************************************************************

/** @brief Read V(k) and T.F(k) for some modes as initial condition.
 *				We provide V[1], F, and Omega[1] (vorticity along 1),
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), \Omega_x(\vec{k}), \f$ 
 *		In typical case Vz is computed using incompressible cond. <BR>
 *		kz = 0 cases are handled as a special case. See fn IncVF::Compute_VyVz(..).
 *
 * @sa void IncVF::Compute_VyVz(int kx, int ky, int kz, complx& Vx, complx Omega, 
 *	complx &Vy, complx &Vz)
 *
 * @note Complex conjugate modes are automatically added.
 */

void  IncFluid::Init_cond_modes_VORTICITY(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Init_cond_modes_VORTICITY_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Init_cond_modes_VORTICITY_RB(T);
}


void  IncFluid::Init_cond_modes_VORTICITY_scalar(IncSF& T)
{
	
	Input_prefix(field_in_file);
	
	(*field_modes_k_array) = 0;  
	(*field_modes_V_array) = 0.0;
	
	int  space_dim = 3;
	int  no_of_vectors_read = 3;
	
	if (my_id == master_id)	
		Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
	
	MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
	int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
	MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
	data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
	MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
									
	
	int kx, ky, kz;	
	complex<DP> Vx, Vy, Vz, Omega, G;
	
	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0; 
	(*T.F) = 0.0;
	
	for (int i = 1; i <= num_field_waveno; i++) 
	{
		kx = (*field_modes_k_array)(i,1); 
		ky = (*field_modes_k_array)(i,2); 
		kz = (*field_modes_k_array)(i,3);	
		
		Vx = (*field_modes_V_array)(i,1); 
		Omega = (*field_modes_V_array)(i,2); 
		
		G = (*field_modes_V_array)(i,3); 
		
		Compute_VyVz(kx, ky, kz, Vx, Omega, Vy, Vz);
		 
		Assign_field_comp_conj(kx, ky, kz, Vx, Vy, Vz);
		T.Assign_field_comp_conj(kx, ky, kz, G);		
		// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
		
		if (my_id == master_id)
			cout	<< "k, V, G : (" << kx << " " << ky  << " " <<kz << ") " 
				<< Vx  << " "  << Vy << " " << Vz << " "  
				<< G << endl;  
	}
	
	if (my_id == master_id)
		cout << "**********************************************" << endl;


}



void  IncFluid::Init_cond_modes_VORTICITY_RB(IncSF& T)
{
	
	if (globalvar_Pr_switch == "PRZERO") 
	{
		Init_cond_modes_VORTICITY();	
		
		*T.F = *V1; 
		Array_divide_ksqr(basis_type, N, *T.F, kfactor);		
	}
	
	
	else if (globalvar_Pr_switch == "PRINFTY") 
	{
		Input_prefix(field_in_file);
		
		(*field_modes_k_array) = 0;  
		(*field_modes_V_array) = 0.0;
		
		int  space_dim = 3;
		int  no_of_vectors_read = 1;
		
		if (my_id == master_id)	
			Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
		
		MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
				  MPI_INT, master_id, MPI_COMM_WORLD); 
		
		data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
				  MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
		
		int kx, ky, kz;	
		complex<DP> G;
		
		(*V1) = 0.0; 
		(*V2) = 0.0; 
		(*V3) = 0.0; 
		(*T.F) = 0.0;
		
		for (int i = 1; i <= num_field_waveno; i++) 
		{
			kx = (*field_modes_k_array)(i,1); 
			ky = (*field_modes_k_array)(i,2); 
			kz = (*field_modes_k_array)(i,3);	
			
			G = (*field_modes_V_array)(i,1); 
			
			T.Assign_field_comp_conj(kx, ky, kz, G);		
			// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, G : (" << kx << " " << ky  << " " <<kz << ") " 
						<< G << endl;  
		}
		
		if (my_id == master_id)
			cout << "============================================" << endl;
		
		Init_cond_Prinfty(T);
	}
	
	else
		Init_cond_modes_VORTICITY_scalar(T);
	
	Zero_modes_RB_slip(T);	
}	
		

//*********************************************************************************************

/** @brief Read V(k) and W(k) for some modes as initial condition.
 *				We provide V[1], W[1], and Omega[1], Omega_W[1] (vorticity along 1),
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), \Omega_x(\vec{k}), \f$ 
 *		In typical case Vz is computed using incompressible cond. <BR>
 *		kz = 0 cases are handled as a special case. See fn IncVF::Compute_VyVz(..).
 *
 * @sa void IncVF::Compute_VyVz(int kx, int ky, int kz, complx& Vx, complx Omega, 
 *	complx &Vy, complx &Vz)
 *
 * @note Complex conjugate modes are automatically added.
 * @note W(k) are handled similarly.
 */
void  IncFluid::Init_cond_modes_VORTICITY(IncVF& W)
{

	Input_prefix(field_in_file);
	
	(*field_modes_k_array) = 0;  
	(*field_modes_V_array) = 0.0;

	int  space_dim = 3;
	int  no_of_vectors_read = 4;
	
	if (my_id == master_id)	
		Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
	
	MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
	int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
	MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
	data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
	MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
									
									
	int kx, ky, kz;	
	complex<DP> Vx, Vy, Vz, Wx, Wy, Wz, Omega, OmegaW;
	
	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0;
	
	(*W.V1) = 0.0; 
	(*W.V2) = 0.0;   
	(*W.V3) = 0.0; 
	
	for (int i = 1; i <= num_field_waveno; i++) 
	{
		kx = (*field_modes_k_array)(i,1); 
		ky = (*field_modes_k_array)(i,2); 
		kz = (*field_modes_k_array)(i,3);	
		
		Vx = (*field_modes_V_array)(i,1); 
		Omega = (*field_modes_V_array)(i,2);  
		
		Wx = (*field_modes_V_array)(i,3); 
		OmegaW = (*field_modes_V_array)(i,4);
		
		Compute_VyVz(kx, ky, kz, Vx, Omega, Vy, Vz); 
		Compute_VyVz(kx, ky, kz, Wx, OmegaW, Wy, Wz);
		
		Assign_field_comp_conj(kx, ky, kz, Vx, Vy, Vz);
		W.Assign_field_comp_conj(kx, ky, kz, Wx, Wy, Wz);
		// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
		
		if (my_id == master_id)
			cout	<< "k, V, W : (" << kx << " " << ky  << " " <<kz << ") " 
					<< Vx << " " << Vy << " "  << Vz << "   " 
					<< Wx << " " << Wy << " "  << Wz << endl;
	}
	
	if (my_id == master_id)
		cout << "**********************************************" << endl;

}


//*********************************************************************************************
// Vector + Scalar
/** @brief Read V(k) and W(k) for some modes as initial condition.
 *				We provide V[1], W[1], F[1], and Omega[1], Omega_W[1] (vorticity along 1),
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ V_x(\vec{k}), \Omega_x(\vec{k}), \f$ 
 *		In typical case Vz is computed using incompressible cond. <BR>
 *		kz = 0 cases are handled as a special case. See fn IncVF::Compute_VyVz(..).
 *
 * @sa void IncVF::Compute_VyVz(int kx, int ky, int kz, complx& Vx, complx Omega, 
 *	complx &Vy, complx &Vz)
 *
 * @note Complex conjugate modes are automatically added.
 * @note W(k) are handled similarly.
 */
void  IncFluid::Init_cond_modes_VORTICITY(IncVF& W, IncSF& T)
{

	Input_prefix(field_in_file);
	
	(*field_modes_k_array) = 0;  
	(*field_modes_V_array) = 0.0;
	
	int  space_dim = 3;
	int  no_of_vectors_read = 5;
	
	if (my_id == master_id)	
		Read_waveno_field(space_dim, no_of_vectors_read, num_field_waveno);
	
	
	MPI_Bcast(&num_field_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
	int data_size = INIT_COND_MAX_NO_FIELD_MODES * 4;
	MPI_Bcast( reinterpret_cast<int*>((*field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
	data_size = INIT_COND_MAX_NO_FIELD_MODES * 16;
	MPI_Bcast( reinterpret_cast<double*>((*field_modes_V_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
									
	
	int kx, ky, kz;	
	complex<DP> Vx, Vy, Vz, Wx, Wy, Wz, Omega, OmegaW, G;
	
	(*V1) = 0.0; 
	(*V2) = 0.0; 
	(*V3) = 0.0;
	
	(*W.V1) = 0.0; 
	(*W.V2) = 0.0;   
	(*W.V3) = 0.0; 
	
	(*T.F) = 0.0;
	
	for (int i = 1; i <= num_field_waveno; i++) 
	{
		kx = (*field_modes_k_array)(i,1); 
		ky = (*field_modes_k_array)(i,2); 
		kz = (*field_modes_k_array)(i,3);	
		
		Vx = (*field_modes_V_array)(i,1); 
		Omega = (*field_modes_V_array)(i,2);  
		
		Wx = (*field_modes_V_array)(i,3); 
		OmegaW = (*field_modes_V_array)(i,4); 
		
		G = (*field_modes_V_array)(i,5);
		
		Compute_VyVz(kx, ky, kz, Vx, Omega, Vy, Vz); 
		Compute_VyVz(kx, ky, kz, Wx, OmegaW, Wy, Wz);
		
		Assign_field_comp_conj(kx, ky, kz, Vx, Vy, Vz);
		W.Assign_field_comp_conj(kx, ky, kz, Wx, Wy, Wz);
		T.Assign_field_comp_conj(kx, ky, kz, G);
		// *V(k)=(Vx,Vy,Vz); Add complex conj if kz=0.
		
		if (my_id == master_id)
			cout	<< "k, V, G : (" << kx << " " << ky  << " " <<kz << ") " 
					<< Vx << " " << Vy << " " << Vz << " " 
					<< Wx << " " << Wy << " " << Wz << " " 
					<< G << endl;
	}
	
	if (my_id == master_id)
		cout << "**********************************************" << endl;
	
}


//******************************  End of Init_cond_modes.cc  **********************************



