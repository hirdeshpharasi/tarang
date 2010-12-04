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
/*! \file  compute_force_modes.cc
 * 
 * @brief Force given set of modes only
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"


//*********************************************************************************************


/** @brief Read wavenos k and F(k) from file force_field_in_file. 
 * 
 *  @param k_dim	Dimension of k vector, e.g., 3 in 3D.
 *	@param F_n		No of components of vector to be read.
 *	@param num_force_waveno		Number of k vector read from the file.
 */
void IncFluid::Read_waveno_force_field(int k_dim, int F_n, int& num_force_waveno)
{
	
	if (my_id == master_id) 
	{
		int i, ki;
		complx Fi;
		
		cout << "****** Reading force field wavenumbers and F  ******" << endl;		
		
		Open_force_waveno_input_files();
		
		int iter = 0;
		while ( ! force_field_in_file.eof()) 
		{	
			force_field_in_file >> ki;  
			
			if (ki < N[1])  
			{								// Stopper.. 
				iter++;	
				(*force_field_modes_k_array)(iter,1) = ki;
				
				for (i=2; i<=k_dim; i++) 
				{
					force_field_in_file >> ki;	
					(*force_field_modes_k_array)(iter,i) = ki;  
				}
		
				for (i=1; i<= F_n; i++) 
				{
					force_field_in_file >> Fi;	
					(*force_field_modes_F_array)(iter,i) = Fi; 
				}
			}
			else break;	
		}
			
		num_force_waveno = iter;
		
		cout << "****** End of Reading force field wavenumbers and components ******" << endl;
		cout  << "Reading of force field configurations ended successfully" << endl;
		
		Close_force_waveno_input_files();
	}	
}	
	
//*********************************************************************************************

/** @brief Read F(k) for some modes as initial condition.  
 *				We provide V[1],..,V[D-1] modes.
 * 
 * @note 2D:   we read (kx, ky) and \f$ F_x(\vec{k}) \f$ except when Ky =0 (x axis).  
 *		Here the read component is Vy, not Vx (and Vx = 0). <BR>
 *		Fy computed using \f$ F_y = -D_x F_x / K_y \f$ <BR>
 *
 * @note 3D:   we read (kx, ky, kz) and \f$ F_x(\vec{k}), V_y(\vec{k}), \f$ 
 *		except when Kz =0 (x-y plane), for which case the read component is Vx and Vz. 
 *		Here we obtain Vy using incompressibility criteria. <BR>
 *		In typical case Vz is computed using incompressible cond. <BR>
 *
 * Complex conjugate modes are automatically added.
 */	
void  IncFluid::Compute_force_given_modes_SIMPLE()
{

	if (is_force_field_para_read == 0)
	{
		Input_force_prefix(force_field_in_file);
		
		(*force_field_modes_k_array) = 0;  
		(*force_field_modes_F_array) = 0.0;

		int  space_dim = 3;
		int  no_of_vectors_read = 2;
		
		if (my_id == master_id)	
			Read_waveno_force_field(space_dim, no_of_vectors_read, num_force_waveno);
			
		MPI_Bcast(&num_force_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = MAX_NO_FORCE_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*force_field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
		data_size = MAX_NO_FORCE_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*force_field_modes_F_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
									
									
		
		int kx, ky, kz;			
		complex<DP> Fx, Fy, Fz;
		
		(*Force1) = 0.0; 
		(*Force2) = 0.0; 
		(*Force3) = 0.0;
		
		for (int i = 1; i <= num_force_waveno; i++) 
		{
			kx = (*force_field_modes_k_array)(i,1); 
			ky = (*force_field_modes_k_array)(i,2); 
			kz = (*force_field_modes_k_array)(i,3);	
			
			Fx = (*force_field_modes_F_array)(i,1); 
			Fy = (*force_field_modes_F_array)(i,2); 
			
			Last_component(kx, ky, kz, Fx, Fy, Fz);
			
			Assign_force_and_comp_conj(kx, ky, kz, Fx, Fy, Fz);
			// *F(k)=(Fx,Fy,Fz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, Force : (" << kx << " " << ky  << " " <<kz << ") " 
						<<  Fx << " " << Fy << " " << Fz <<  endl;  
		}
		
		if (my_id == master_id) 
			cout << "******************************************" << endl;	
			
		
		is_force_field_para_read = 1;		// force field modes read now
	}
	
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


void IncFluid::Compute_force_given_modes_SIMPLE(IncSF& T)
{

	if (is_force_field_para_read == 0)
	{
		Input_force_prefix(force_field_in_file);
		
		(*force_field_modes_k_array) = 0;  
		(*force_field_modes_F_array) = 0.0;

		int  space_dim = 3;
		int  no_of_vectors_read = 3;
		
		if (my_id == master_id)	
			Read_waveno_force_field(space_dim, no_of_vectors_read, num_force_waveno);
			
		MPI_Bcast(&num_force_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = MAX_NO_FORCE_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*force_field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
		data_size = MAX_NO_FORCE_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*force_field_modes_F_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
									
									
		
		int kx, ky, kz;			
		complex<DP> Fx, Fy, Fz, G;
		
		(*Force1) = 0.0; 
		(*Force2) = 0.0; 
		(*Force3) = 0.0;
		(*T.Force) = 0.0;
		
		for (int i = 1; i <= num_force_waveno; i++) 
		{
			kx = (*force_field_modes_k_array)(i,1); 
			ky = (*force_field_modes_k_array)(i,2); 
			kz = (*force_field_modes_k_array)(i,3);	
			
			Fx = (*force_field_modes_F_array)(i,1); 
			Fy = (*force_field_modes_F_array)(i,2); 
			G = (*force_field_modes_F_array)(i,3); 
			
			Last_component(kx, ky, kz, Fx, Fy, Fz);
			
			Assign_force_and_comp_conj(kx, ky, kz, Fx, Fy, Fz);
			
			T.Assign_force_and_comp_conj(kx, ky, kz, G);
			// *F(k)=(Fx,Fy,Fz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, Force : (" << kx << " " << ky  << " " <<kz << ") " 
						<< Fx << " " << Fy << " " << Fz << "    " 
						<< G << endl;  
		}
		
		if (my_id == master_id) 
			cout << "******************************************" << endl;	
		
		
		is_force_field_para_read = 1;		// force field modes read now
	}
	
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
void  IncFluid::Compute_force_given_modes_SIMPLE(IncVF& W)
{

	if (is_force_field_para_read == 0)
	{
		Input_force_prefix(force_field_in_file);
		
		(*force_field_modes_k_array) = 0;  
		(*force_field_modes_F_array) = 0.0;

		int  space_dim = 3;
		int  no_of_vectors_read = 4;
		
		if (my_id == master_id)	
			Read_waveno_force_field(space_dim, no_of_vectors_read, num_force_waveno);
			
		MPI_Bcast(&num_force_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = MAX_NO_FORCE_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*force_field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
		data_size = MAX_NO_FORCE_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*force_field_modes_F_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
									
									
		
		int kx, ky, kz;			
		complex<DP> Fx, Fy, Fz;
		complex<DP> FxW, FyW, FzW;
		
		(*Force1) = 0.0; 
		(*Force2) = 0.0; 
		(*Force3) = 0.0;
		
		(*W.Force1) = 0.0; 
		(*W.Force2) = 0.0; 
		(*W.Force3) = 0.0;
		
		for (int i = 1; i <= num_force_waveno; i++) 
		{
			kx = (*force_field_modes_k_array)(i,1); 
			ky = (*force_field_modes_k_array)(i,2); 
			kz = (*force_field_modes_k_array)(i,3);	
			
			Fx = (*force_field_modes_F_array)(i,1); 
			Fy = (*force_field_modes_F_array)(i,2);
			FxW = (*force_field_modes_F_array)(i,3); 
			FyW = (*force_field_modes_F_array)(i,4); 
			
			Last_component(kx, ky, kz, Fx, Fy, Fz);
			Last_component(kx, ky, kz, FxW, FyW, FzW);
			
			Assign_force_and_comp_conj(kx, ky, kz, Fx, Fy, Fz);
			W.Assign_force_and_comp_conj(kx, ky, kz, FxW, FyW, FzW);
			// *F(k)=(Fx,Fy,Fz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, Force : (" << kx << " " << ky  << " " <<kz << ") " 
						<< Fx	<< " " << Fy  << " " << Fz << " " 
						<< FxW	<< " " << FyW << " " << FzW  <<  endl;  
		}
		
		if (my_id == master_id) 
			cout << "******************************************" << endl;	
		
			
		is_force_field_para_read = 1;		// force field modes read now
	}
	
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
void  IncFluid::Compute_force_given_modes_SIMPLE(IncVF& W, IncSF& T)
{

	if (is_force_field_para_read == 0)
	{
		Input_force_prefix(force_field_in_file);
		
		(*force_field_modes_k_array) = 0;  
		(*force_field_modes_F_array) = 0.0;

		int  space_dim = 3;
		int  no_of_vectors_read = 5;
		
		if (my_id == master_id)	
			Read_waveno_force_field(space_dim, no_of_vectors_read, num_force_waveno);
			
		MPI_Bcast(&num_force_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = MAX_NO_FORCE_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*force_field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
		data_size = MAX_NO_FORCE_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*force_field_modes_F_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
									
									
		
		int kx, ky, kz;			
		complex<DP> Fx, Fy, Fz;
		complex<DP> FxW, FyW, FzW, G;
		
		(*Force1) = 0.0; 
		(*Force2) = 0.0; 
		(*Force3) = 0.0;
		
		(*W.Force1) = 0.0; 
		(*W.Force2) = 0.0; 
		(*W.Force3) = 0.0;
		
		(*T.Force) = 0.0;
		
		for (int i = 1; i <= num_force_waveno; i++) 
		{
			kx = (*force_field_modes_k_array)(i,1); 
			ky = (*force_field_modes_k_array)(i,2); 
			kz = (*force_field_modes_k_array)(i,3);	
			
			Fx = (*force_field_modes_F_array)(i,1); 
			Fy = (*force_field_modes_F_array)(i,2);
			
			FxW = (*force_field_modes_F_array)(i,3); 
			FyW = (*force_field_modes_F_array)(i,4); 
			
			G   = (*force_field_modes_F_array)(i,3);
			
			Last_component(kx, ky, kz, Fx, Fy, Fz);
			Last_component(kx, ky, kz, FxW, FyW, FzW);
			
			Assign_force_and_comp_conj(kx, ky, kz, Fx, Fy, Fz);
			W.Assign_force_and_comp_conj(kx, ky, kz, FxW, FyW, FzW);
			T.Assign_force_and_comp_conj(kx, ky, kz, G);
			// *F(k)=(Fx,Fy,Fz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, Force : (" << kx << " " << ky  << " " <<kz << ") " 
						<< Fx << " " << Fy << " " << Fz << " " 
						<< FxW << " " << FyW << " " << FzW  <<  " " 
						<< G << endl;  
		}
		
		if (my_id == master_id) 
			cout << "******************************************" << endl;
			
			
		is_force_field_para_read = 1;		// force field modes read now
	}
	
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
void  IncFluid::Compute_force_given_modes_VORTICITY()
{

	if (is_force_field_para_read == 0)
	{
		Input_force_prefix(force_field_in_file);
		
		(*force_field_modes_k_array) = 0;  
		(*force_field_modes_F_array) = 0.0;

		int  space_dim = 3;
		int  no_of_vectors_read = 2;
		
		if (my_id == master_id)	
			Read_waveno_force_field(space_dim, no_of_vectors_read, num_force_waveno);
			
		MPI_Bcast(&num_force_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = MAX_NO_FORCE_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*force_field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
		data_size = MAX_NO_FORCE_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*force_field_modes_F_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
									
									
		
		int kx, ky, kz;			
		complex<DP> Fx, Fy, Fz, F_omega;
		
		(*Force1) = 0.0; 
		(*Force2) = 0.0; 
		(*Force3) = 0.0;
		
		for (int i = 1; i <= num_force_waveno; i++) 
		{
			kx = (*force_field_modes_k_array)(i,1); 
			ky = (*force_field_modes_k_array)(i,2); 
			kz = (*force_field_modes_k_array)(i,3);	
			
			Fx = (*force_field_modes_F_array)(i,1); 
			F_omega = (*force_field_modes_F_array)(i,2); 
			
			Compute_VyVz(kx, ky, kz, Fx, F_omega, Fy, Fz);
			Assign_force_and_comp_conj(kx, ky, kz, Fx, Fy, Fz);
			// *Force(k)=(Fx,Fy,Fz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, Force : (" << kx << " " << ky  << " " <<kz << ") " 
						<< Fx << " " << Fy << " " << Fz <<  endl;  
		}
		
		if (my_id == master_id) 
			cout << "******************************************" << endl;
			
		is_force_field_para_read = 1;		// force field modes read now
	}

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
void  IncFluid::Compute_force_given_modes_VORTICITY(IncSF& T)
{

	if (is_force_field_para_read == 0)
	{
		Input_force_prefix(force_field_in_file);
		
		(*force_field_modes_k_array) = 0;  
		(*force_field_modes_F_array) = 0.0;

		int  space_dim = 3;
		int  no_of_vectors_read = 3;
		
		if (my_id == master_id)	
			Read_waveno_force_field(space_dim, no_of_vectors_read, num_force_waveno);
			
		MPI_Bcast(&num_force_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = MAX_NO_FORCE_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*force_field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
		data_size = MAX_NO_FORCE_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*force_field_modes_F_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
		
		int kx, ky, kz;			
		complex<DP> Fx, Fy, Fz, F_omega, G;
		
		(*Force1) = 0.0; 
		(*Force2) = 0.0; 
		(*Force3) = 0.0;
		(*T.Force) = 0.0;
		
		for (int i = 1; i <= num_force_waveno; i++) 
		{
			kx = (*force_field_modes_k_array)(i,1); 
			ky = (*force_field_modes_k_array)(i,2); 
			kz = (*force_field_modes_k_array)(i,3);	
			
			Fx = (*force_field_modes_F_array)(i,1); 
			F_omega = (*force_field_modes_F_array)(i,2);
			
			G = (*force_field_modes_F_array)(i,3); 
			
			Compute_VyVz(kx, ky, kz, Fx, F_omega, Fy, Fz);
			
			Assign_force_and_comp_conj(kx, ky, kz, Fx, Fy, Fz);
			T.Assign_force_and_comp_conj(kx, ky, kz, G);
			// *Force(k)=(Fx,Fy,Fz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, Force, G : (" << kx << " " << ky  << " " <<kz << ") " 
						<< Fx << " " << Fy << " " << Fz << " " 
						<< G  <<  endl;  
		}
		
		if (my_id == master_id) 
			cout << "******************************************" << endl;
						
		is_force_field_para_read = 1;		// force field modes read now
	}

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
void  IncFluid::Compute_force_given_modes_VORTICITY(IncVF& W)
{
	
	if (is_force_field_para_read == 0)
	{
		Input_force_prefix(force_field_in_file);
		
		(*force_field_modes_k_array) = 0;  
		(*force_field_modes_F_array) = 0.0;

		int  space_dim = 3;
		int  no_of_vectors_read = 4;
		
		if (my_id == master_id)	
			Read_waveno_force_field(space_dim, no_of_vectors_read, num_force_waveno);
			
		MPI_Bcast(&num_force_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = MAX_NO_FORCE_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*force_field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
		data_size = MAX_NO_FORCE_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*force_field_modes_F_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
									
									
		
		int kx, ky, kz;			
		complex<DP> Fx, Fy, Fz, F_omega;
		complex<DP> FxW, FyW, FzW, F_omegaW;
		
		(*Force1) = 0.0; 
		(*Force2) = 0.0; 
		(*Force3) = 0.0;
		
		(*W.Force1) = 0.0; 
		(*W.Force2) = 0.0; 
		(*W.Force3) = 0.0;
		
		for (int i = 1; i <= num_force_waveno; i++) 
		{
			kx = (*force_field_modes_k_array)(i,1); 
			ky = (*force_field_modes_k_array)(i,2); 
			kz = (*force_field_modes_k_array)(i,3);	
			
			Fx = (*force_field_modes_F_array)(i,1); 
			F_omega = (*force_field_modes_F_array)(i,2); 
			
			FxW = (*force_field_modes_F_array)(i,3); 
			F_omegaW = (*force_field_modes_F_array)(i,4); 
			
			Compute_VyVz(kx, ky, kz, Fx, F_omega, Fy, Fz);
			Compute_VyVz(kx, ky, kz, FxW, F_omegaW, FyW, FzW);
			
			Assign_force_and_comp_conj(kx, ky, kz, Fx, Fy, Fz);
			W.Assign_force_and_comp_conj(kx, ky, kz, FxW, FyW, FzW);
			// *Force(k)=(Fx,Fy,Fz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, Force : (" << kx << " " << ky  << " " <<kz << ") " 
						<< Fx << " "  << Fy  << " " << Fz <<  " " 
						<< FxW << " " << FyW << " " << FzW << endl;  
		}
		
		if (my_id == master_id) 
			cout << "******************************************" << endl;	
		
		is_force_field_para_read = 1;		// force field modes read now
	}

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
void  IncFluid::Compute_force_given_modes_VORTICITY(IncVF& W, IncSF& T)
{

	if (is_force_field_para_read == 0)
	{
		Input_force_prefix(force_field_in_file);
		
		(*force_field_modes_k_array) = 0;  
		(*force_field_modes_F_array) = 0.0;

		int  space_dim = 3;
		int  no_of_vectors_read = 5;
		
		if (my_id == master_id)	
			Read_waveno_force_field(space_dim, no_of_vectors_read, num_force_waveno);
			
		MPI_Bcast(&num_force_waveno, 1, MPI_INT, master_id, MPI_COMM_WORLD); 
		
		int data_size = MAX_NO_FORCE_FIELD_MODES * 4;
		MPI_Bcast( reinterpret_cast<int*>((*force_field_modes_k_array).data()), data_size, 
									MPI_INT, master_id, MPI_COMM_WORLD); 
	
		data_size = MAX_NO_FORCE_FIELD_MODES * 16;
		MPI_Bcast( reinterpret_cast<double*>((*force_field_modes_F_array).data()), data_size, 
									MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
									
									
		
		int kx, ky, kz;			
		complex<DP> Fx, Fy, Fz, F_omega;
		complex<DP> FxW, FyW, FzW, F_omegaW, G;
		
		(*Force1) = 0.0; 
		(*Force2) = 0.0; 
		(*Force3) = 0.0;
		
		(*W.Force1) = 0.0; 
		(*W.Force2) = 0.0; 
		(*W.Force3) = 0.0;
		
		(*T.Force) = 0.0;
		
		for (int i = 1; i <= num_force_waveno; i++) 
		{
			kx = (*force_field_modes_k_array)(i,1); 
			ky = (*force_field_modes_k_array)(i,2); 
			kz = (*force_field_modes_k_array)(i,3);	
			
			Fx = (*force_field_modes_F_array)(i,1); 
			F_omega = (*force_field_modes_F_array)(i,2); 
			
			FxW = (*force_field_modes_F_array)(i,3); 
			F_omegaW = (*force_field_modes_F_array)(i,4); 
			
			G = (*force_field_modes_F_array)(i,5); 
			
			Compute_VyVz(kx, ky, kz, Fx, F_omega, Fy, Fz);
			Compute_VyVz(kx, ky, kz, FxW, F_omegaW, FyW, FzW);
			
			Assign_force_and_comp_conj(kx, ky, kz, Fx, Fy, Fz);
			W.Assign_force_and_comp_conj(kx, ky, kz, FxW, FyW, FzW);
			T.Assign_force_and_comp_conj(kx, ky, kz, G);
			// *Force(k)=(Fx,Fy,Fz); Add complex conj if kz=0.
			
			if (my_id == master_id)
				cout	<< "k, Force : (" << kx << " " << ky  << " " <<kz << ") " 
						<< Fx  << " " << Fy << " "  << Fz <<  " " 
						<< FxW << " " << FyW << " " << FzW << " " 
						<< G << endl;  
		}
		
		if (my_id == master_id) 
			cout << "******************************************" << endl;	
		
		is_force_field_para_read = 1;		// force field modes read now
	}

}



//*******************************  End of compute_force_modes.cc  *****************************


