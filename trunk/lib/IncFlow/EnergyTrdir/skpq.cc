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


/*! \file  skpq.cc
 * 
 * @brief Compute Skpq for a given set of tirads.
 * @sa DP IncVF::Get_Suu(TinyVector<DP,2> k, TinyVector<DP,2> p)
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bugs  Lots of problems...
 */


#include "../IncVF.h"
#include "../IncSF.h"

//***************************************  3D  ************************************************

// u(p) to u(k)
DP IncVF::Get_Suu(TinyVector<int,3> k, TinyVector<int,3> p)
{
	int lx, ly;
	
	TinyVector<int,3> q;
	TinyVector<complx,3>  uk_conj, uk, up, uq, kk;
	
	// get u(k)^*
	lx = Get_lx(basis_type, k(0), N);	
	ly = Get_ly3D(basis_type, k(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		uk =  (*V1)(lx, ly, k(2)), (*V2)(lx, ly, k(2)), (*V3)(lx, ly, k(2));
	
	int data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(uk.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
	
	// get u(p)
	lx = Get_lx(basis_type, p(0), N);
	ly = Get_ly3D(basis_type, p(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		up = (*V1)(lx, ly, p(2)), (*V2)(lx, ly, p(2)), (*V3)(lx, ly, p(2));
	
	data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(up.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
					
	// get u(q)						
	q = k(0)-p(0), k(1)-p(1), k(2)-p(2);
	
	lx = Get_lx(basis_type, q(0), N);	
	ly = Get_ly3D(basis_type, q(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		uq = (*V1)(lx, ly, q(2)), (*V2)(lx, ly, q(2)), (*V3)(lx, ly, q(2));
		
	data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(uq.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
					
					
	if (my_id == master_id)
	{
		if (basis_type == "SCFT")
		{
			uk_conj = conj(uk(0)*(-I)), conj(uk(1)), conj(uk(2));
			up = up(0)*(-I), up(1), up(2);
			uq = uq(0)*(-I), uq(1), uq(2);
		}
		
		kk = complx(k(0)* kfactor[1], 0.0),  
			 complx(k(1)* kfactor[2], 0.0), 
			 complx(k(2)* kfactor[3], 0.0); 
			 
			
		return imag( dot(kk,uq) * dot(uk_conj, up) );
	}	
	
	else
		return 0;		// for -Wall

}

//*********************************************************************************************
// w(p) to w(k)
DP IncVF::Get_Sww(IncVF& W, TinyVector<int,3> k, TinyVector<int,3> p)
{

	int lx, ly;

	TinyVector<int,3> q;
	TinyVector<complx,3>  wk_conj, wk, wp, uq, kk;
	
	
	// get u(k)^*
	lx = Get_lx(basis_type, k(0), N);	
	ly = Get_ly3D(basis_type, k(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		wk =  (*W.V1)(lx, ly, k(2)), (*W.V2)(lx, ly, k(2)), (*W.V3)(lx, ly, k(2));
	
	int data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(wk.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
	
	// get u(p)
	lx = Get_lx(basis_type, p(0), N);
	ly = Get_ly3D(basis_type, p(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		wp = (*W.V1)(lx, ly, p(2)), (*W.V2)(lx, ly, p(2)), (*W.V3)(lx, ly, p(2));
	
	data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(wp.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
					
	// get u(q)						
	q = k(0)-p(0), k(1)-p(1), k(2)-p(2);
	
	lx = Get_lx(basis_type, q(0), N);	
	ly = Get_ly3D(basis_type, q(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		uq = (*V1)(lx, ly, q(2)), (*V2)(lx, ly, q(2)), (*V3)(lx, ly, q(2));
		
	data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(uq.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
					
	
	
	if (my_id == master_id)
	{
		if (basis_type == "SCFT")
		{
			wk_conj = conj(wk(0)*(-I)), conj(wk(1)), conj(wk(2));
			wp = wp(0)*(-I), wp(1), wp(2);
			uq = uq(0)*(-I), uq(1), uq(2);
		}
		
		kk = complx(k(0)* kfactor[1], 0.0),  
			 complx(k(1)* kfactor[2], 0.0), 
			 complx(k(2)* kfactor[3], 0.0); 
			 
			
		return imag( dot(kk,uq) * dot(wk_conj, wp) );
	}	
	
	else
		return 0;   // for -Wall
	
	
}

//*********************************************************************************************
// u(p) to w(k)
DP IncVF::Get_Swu(IncVF& W, TinyVector<int,3> k, TinyVector<int,3> p)
{

	int lx, ly;

	TinyVector<int,3> q;
	TinyVector<complx,3>  wk_conj, wk, up, wq, kk;
	
	
	// get u(k)^*
	lx = Get_lx(basis_type, k(0), N);	
	ly = Get_ly3D(basis_type, k(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		wk =  (*W.V1)(lx, ly, k(2)), (*W.V2)(lx, ly, k(2)), (*W.V3)(lx, ly, k(2));
	
	int data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(wk.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
	
	// get u(p)
	lx = Get_lx(basis_type, p(0), N);
	ly = Get_ly3D(basis_type, p(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		up = (*V1)(lx, ly, p(2)), (*V2)(lx, ly, p(2)), (*V3)(lx, ly, p(2));
	
	data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(up.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
					
	// get u(q)						
	q = k(0)-p(0), k(1)-p(1), k(2)-p(2);
	
	lx = Get_lx(basis_type, q(0), N);
	ly = Get_ly3D(basis_type, q(1), N);	
	
	if ((lx >= 0) && (lx < local_N1))
		wq = (*W.V1)(lx, ly, q(2)), (*W.V2)(lx, ly, q(2)), (*W.V3)(lx, ly, q(2));
		
	data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(wq.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
					
	
	
	if (my_id == master_id)
	{
		if (basis_type == "SCFT")
		{
			wk_conj = conj(wk(0)*(-I)), conj(wk(1)), conj(wk(2));
			up = up(0)*(-I), up(1), up(2);
			wq = wq(0)*(-I), wq(1), wq(2);
		}
		
		kk = complx(k(0)* kfactor[1], 0.0),  
			 complx(k(1)* kfactor[2], 0.0), 
			 complx(k(2)* kfactor[3], 0.0); 
			 
			
		return imag( dot(kk,wq) * dot(wk_conj, up) );
	}	
	
	else
		return 0;   // for -Wall

}

//*********************************************************************************************
// T(p) to T(k)
DP IncVF::Get_Stt(IncSF& T, TinyVector<int,3> k, TinyVector<int,3> p)
{

	int lx, ly;
	
	TinyVector<int,3>	  q;
	
	TinyVector<complx,1>	tk_conj, tk, tp;
	TinyVector<complx,3>  uq, kk;
	
	
	// get u(k)^*
	lx = Get_lx(basis_type, k(0), N);
	ly = Get_ly3D(basis_type, k(1), N);	
	
	if ((lx >= 0) && (lx < local_N1))
		tk =  (*T.F)(lx, ly, k(2));
	
	int data_size = 2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(tk.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
	
	// get u(p)
	lx = Get_lx(basis_type, p(0), N);
	ly = Get_ly3D(basis_type, p(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		tp = (*T.F)(lx, ly, p(2));
	
	data_size = 2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(tp.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
					
	// get u(q)						
	q = k(0)-p(0), k(1)-p(1), k(2)-p(2);
	
	lx = Get_lx(basis_type, q(0), N);	
	ly = Get_ly3D(basis_type, q(1), N);
	
	if ((lx >= 0) && (lx < local_N1))
		uq = (*V1)(lx, ly, q(2)), (*V2)(lx, ly, q(2)), (*V3)(lx, ly, q(2));
		
	data_size = 3*2;  // factor 2 due to complex to double conversion
	MPI_Bcast(reinterpret_cast<double*>(uq.data()), data_size, MPI_DOUBLE, 
					master_id, MPI_COMM_WORLD);
					
	
	
	if (my_id == master_id)
	{
		if (basis_type == "SCFT")
		{
			tk_conj = conj(tk(0)*(-I));
			tp = tp(0)*(-I);
			uq = uq(0)*(-I), uq(1), uq(2);
		}
		
		kk = complx(k(0)* kfactor[1], 0.0),  
			 complx(k(1)* kfactor[2], 0.0), 
			 complx(k(2)* kfactor[3], 0.0); 
			 
			
		return imag( dot(kk,uq) * dot(tk_conj,tp) );
	}	
	
	else
		return 0;   // for -Wall
		
}


//*********************************************************************************************


/** @brief Compute Suu for set of triads specified in triad_array. 
 * 
 *  @note	Triads are fed in triad_array.  triad_array(i,:) = kx, ky, px, py, 0 , 0 in 2D.
 *			while triad_array(i,:) = kx, ky, kz, px, py, pz in 3D.
 *	
 * @return (*Sself_array)(i) = Suu(i) for all i triads 
 */
void IncVF::Compute_Skpq()
{
	no_of_Skpq_triads=1;
			
	TinyVector<int,3> k,p;
	
	// Place k & p of triads here.
	(*triad_array)(1,Range::all()) = 3,3,2, 1,1,1;
	
	for (int i=1; i<=no_of_Skpq_triads; i++)
	{
		k = (*triad_array)(i,0), (*triad_array)(i,1), (*triad_array)(i,2);
		p = (*triad_array)(i,3), (*triad_array)(i,4), (*triad_array)(i,5);
		
		(*Sself_array)(i) = Get_Suu(k, p);
	}
		
}	


//*********************************************************************************************

/** @brief Compute Suu and Stt for set of triads specified in triad_array. 
 * 
 *  @note	Triads are fed in triad_array.  triad_array(i,:) = kx, ky, px, py, 0 , 0 in 2D.
 *			while triad_array(i,:) = kx, ky, kz, px, py, pz in 3D.
 *	
 * @return (*Sself_array)(i) = Suu(i) for all i triads 
 * @return (*S_SF_array)(i)  = Stt(i) for all i triads 
 */
void IncVF::Compute_Skpq(IncSF& T)
{
	no_of_Skpq_triads=1;
		
	TinyVector<int,3> k,p;
	
	// Place k & p of triads here.
	(*triad_array)(1,Range::all()) = 3,3,2, 2,2, 1;
	
	for (int i=1; i<=no_of_Skpq_triads; i++)
	{
		k = (*triad_array)(i,0), (*triad_array)(i,1), (*triad_array)(i,2);
		p = (*triad_array)(i,3), (*triad_array)(i,4), (*triad_array)(i,5);
		
		(*Sself_array)(i) = Get_Suu(k, p);
		(*S_SF_array)(i) = Get_Stt(T, k, p);
	}
		
}



//*********************************************************************************************

/** @brief Compute Suu, Sww, Swu [u to w] for set of triads specified in triad_array. 
 * 
 *  @note	Triads are fed in triad_array.  triad_array(i,:) = kx, ky, px, py, 0 , 0 in 2D.
 *			while triad_array(i,:) = kx, ky, kz, px, py, pz in 3D.
 *	
 * @return (*Sself_array)(i) = Suu(i) for all i triads 
 * @return (*W.Sself_array)(i) = Sww(i) for all i triads 
 * @return (*S_to_VF_array)(i)  = Swu(i) for all i triads [ u to w]
 */
void IncVF::Compute_Skpq(IncVF& W)
{
	no_of_Skpq_triads=1;
	
	TinyVector<int,3> k,p;
	
	// Place k & p of triads here.
	(*triad_array)(1,Range::all()) = 3,2,1, 1,1,1;
	
	for (int i=1; i<=no_of_Skpq_triads; i++)
	{
		k = (*triad_array)(i,0), (*triad_array)(i,1), (*triad_array)(i,2);
		p = (*triad_array)(i,3), (*triad_array)(i,4), (*triad_array)(i,5);
		
		(*Sself_array)(i) = Get_Suu(k, p);
		(*W.Sself_array)(i) = Get_Sww(W, k, p);
		(*S_to_VF_array)(i) = Get_Swu(W, k, p);
	}	
	
}

//*********************************************************************************************

/** @brief Compute Suu, Sww, Swu [u to w], Stt for set of triads specified in triad_array. 
 * 
 *  @note	Triads are fed in triad_array.  triad_array(i,:) = kx, ky, px, py, 0 , 0 in 2D.
 *			while triad_array(i,:) = kx, ky, kz, px, py, pz in 3D.
 *	
 * @return (*Sself_array)(i)	= Suu(i) for all i triads 
 * @return (*W.Sself_array)(i)	= Sww(i) for all i triads 
 * @return (*S_to_VF_array)(i)  = Swu(i) for all i triads [ u to w]
 * @return (*S_SF_array)(i)		= Stt(i) for all i triads 
 */
void IncVF::Compute_Skpq(IncVF& W, IncSF& T)
{
	no_of_Skpq_triads=1;

	TinyVector<int,3> k,p;
	
	// Place k & p of triads here.
	(*triad_array)(1,Range::all()) = 3,2,1, 1,1,1;
	
	for (int i=1; i<=no_of_Skpq_triads; i++)
	{
		k = (*triad_array)(i,0), (*triad_array)(i,1), (*triad_array)(i,2);
		p = (*triad_array)(i,3), (*triad_array)(i,4), (*triad_array)(i,5);
		
		(*Sself_array)(i) = Get_Suu(k, p);
		(*W.Sself_array)(i) = Get_Sww(W, k, p);
		(*S_to_VF_array)(i) = Get_Swu(W, k, p);
		(*S_SF_array)(i) = Get_Stt(T, k, p);
	}
	
}



//*******************************  End of Skpq.cc  ********************************************



