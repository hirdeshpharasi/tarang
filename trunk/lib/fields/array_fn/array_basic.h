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

/*! \file field_basic.h
 * 
 * @brief Some basic Array operations: Array_real_mult, Output_asreal, Model_initial_shell_spectrum
 *
 * @author  M. K. Verma
 * @date Sept 2008
 *
 * @bugs Out of date functions Read_data_split_arrays_MPI, Write_data_split_arrays_MPI contain bugs.
 */ 
 

#ifndef _FIELD_BASIC
#define _FIELD_BASIC

#include <mpi.h>
#include <fftw3-mpi.h>

#include <blitz/array.h>
#include <complex>
#include <cmath>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif


#include "array_basic_inline.h"

using namespace blitz;

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

/** @brief Multiply real arrays A and B term by term, and put the result in C.  
 *
 *  The multiplication is done in each process independently 
 *	without any requirement of communications.
 *
 * @param  A  Real array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$. 
 * @param  B  Real array \f$ B [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$. 
 * @param  N[]  The size of the arrays A, B.
 * 
 * @return   C(i,j,k) = A(i,j,k)*B(i,j,k).  
 */
void Array_real_mult
(
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<complx,3> C
);

void Array_real_divide
(
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<complx,3> C
);


//*********************************************************************************************

/** @brief Output real array complete A to file_out
 *
 * The master_id first outputs its A.  Then it receives the data from other processors 
 *		in sequence and outputs them.
 *
 * @param  A  Real array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$. 
 * @param  N[]  The size of the array A.
 * @param  file_out	Output file to which A is to be written.
 * @param  temp_array Temporary array \f$  [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$ 
 *			for receiving array data from other processors.
 * 
 * @return   Complete Array A is written in file_out. 
 */
void Output_asreal
(
	ofstream& file_out, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_array,
	string format
);


void Output_asreal_transpose_order
(
	 ofstream& file_out, 
	 int N[], 
	 Array<complx,3> A, 
	 Array<complx,3> temp_r,
	string format
 );

/*
void Output_asreal_hdf5
(
	DataSet* dataset,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_array
);
*/
//*********************************************************************************************

/** @brief Read complex array complete A from file_in
 *
 * Complete A is stored in file_in in the unstructured format 
 *		(series of blitz++ complex data format) <BR>
 * The master_id first reads its A.  Then it reads the data for other processors in sequence, 
 *	and then sends the data to the processor.
 *
 * @param  A  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$. 
 * @param  N[]  The size of the array A.
 * @param  file_in  The file containing complete A data in complex data format.
 * @param  temp_array Temporary array \f$  [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$ 
 *			for receiving array data from other processors.
 * 
 * @return   Array A in each processor read in.
 */
void Read_data_MPI
(
	string basis_type,
	ifstream& file_in, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_array,
	string format
);

void Read_data_kz0plane_MPI
(
	string basis_type,
	ifstream& file_in, 
	int N[], 
	Array<complx,3> A,
	string format
);


//*********************************************************************************************

/** @brief Read complex array reducedA from file_in
 *
 * reducedA is stored in file_in in the unstructured format 
 *		(series of blitz++ complex data format) <BR>
 * 
 * Steps: <BR>
 * 
 * (1) Compute last processor that would contain reduced A:    
 *     If Nreduced[1] <= local_N1, all of reducedA is stored in master_id.  
 *     Otherwise it is distribute till processor last_dist. <BR>
 *		last_dest = Nreduced[1]/local_N1-1 (if local_N1 divisible by Nreduced[1]). <BR>
 *		last_dest = Nreduced[1]/local_N1 (if local_N1 indivisible by Nreduced[1]). <BR>
 *
 * (2) Read data into master_id.  The last_l1 = min(local_N1, Nreduced[1])-1. <BR>
 * (3) If last_dist > 0, distribute data to other processors in sequence.  
 *		We calculate for each processor <BR>
 *			remaining_i = Nreduced[1] - dest*local_N1. <BR>
 *			last_i = min(local_N1, remaining_i) - 1. <BR>
 *     and keep filling A in each processor.
 *     Note that data is first read in temp_array of master_id, and then it is transferred 
 *			to the processor.
 *
 * @param  A  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$. 
 * @param  Nreduced[]  The size of the reduced array A.
 * @param  file_in	The file containing reduced A data in complex data format.
 * @param  temp_array Temporary array \f$  [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$ 
 *			for receiving array data from other processors.
 * 
 * @return   Array A in each processor. 
 */
void Read_data_MPI
(
	string basis_type,
	ifstream& file_in, 
	int N[], 
	int Nreduced[], 
	Array<complx,3> A,
	string format
);

void Read_data_kz0plane_MPI
(
 string basis_type,
 ifstream& file_in, 
 int N[], 
 int Nreduced[], 
 Array<complx,3> A,
 string format
);

//*********************************************************************************************

/** @brief Write complex array complete A to file_out
 *
 * The master_id first write its A to file_out.  Then it recieves A from other processors 
 *			in sequence, and then writes to file_out.
 *
 * @param  A  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$. 
 * @param  N[]  The size of the array A. 
 * @param  Nreduced[]  The size of the array reducedA.
 * @param  file_out  Output file to which A is to be written.
 * @param  temp_array Temporary array \f$  [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$ 
 *			for receiving array data from other processors.
 * 
 * @return   file_out contains complte A in unstructured format 
 *				(series of blitz++ complex data format).
 */
void Write_data_MPI
(
	string basis_type,
	ofstream& file_out, 
	int N[], 
	Array<complx,3> A,  
	Array<complx,3> temp_array,
	string format
);

void Write_data_kz0plane_MPI
(
 string basis_type,
 ofstream& file_out, 
 int N[], 
 Array<complx,3> A,
 string format
);

/*
void Write_data_MPI_hdf5
(
	DataSet* datasetR,
        DataSet* datasetC,	
	int N[], 
	Array<complx,3> A,  
	Array<complx,3> temp_array
);
*/
//*********************************************************************************************

/** @brief Write complex array reducedA to file_out
 *
 * reducedA is sent to file_out in the unstructured format 
 *		(series of blitz++ complex data format) <BR>
 * 
 * Steps: <BR>
 * 
 * (1) Compute last processor that would contain reduced A:  <BR>
 *     If Nreduced[1] <= local_N1, all of reducedA is stored in master_id.  
 *     Otherwise it is distribute till processor last_dist. <BR>
 *		last_dest = Nreduced[1]/local_N1-1 (if local_N1 divisible by Nreduced[1]). <BR>
 *		last_dest = Nreduced[1]/local_N1 (if local_N1 indivisible by Nreduced[1]). <BR>
 *
 * (2) Master_id writes its data to file_out.  The last_l1 = min(local_N1, Nreduced[1])-1.<BR>
 * (3) If last_dist > 0, distribute data to other processors in sequence.  
 *		We calculate for each processor <BR>
 *			remaining_i = Nreduced[1] - dest*local_N1. <BR>
 *			last_i = min(local_N1, remaining_i) - 1. <BR>
 *     and keep filling A in each processor.
 *     Note that data is first read in temp_array of master_id, 
 *			and then it is transferred to the processor.
 *
 * @param  A  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$. 
 * @param  Nreduced[]  The size of the reduced array A.
 * @param  file_in	The file containing reduced A data in complex data format.
 * @param  temp_array Temporary array \f$  [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$ 
 *			for receiving array data from other processors.
 * 
 * @return   Array A in each processor. 
 */
void Write_data_MPI
(
	string basis_type,
	ofstream& file_out, 
	int N[], 
	int Nreduced[], 
	Array<complx,3> A,
	string format
);

void Write_data_kz0plane_MPI
(
 string basis_type,
 ofstream& file_out, 
 int N[], 
 int Nreduced[], 
 Array<complx,3> A,
 string format
);

//*********************************************************************************************

/** @brief Collect pieced of kz planes from all the procs and put them in plane
 *
 * @param  A  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$. 
 * @param  N[]  The size of the array A.
 * @param  Nreduced[]  The size of the array reducedA.
 * @param  plane  2D array for storing completeA(:,:,kz).
 * @param  kz	grid wavenumber along z direction whose planar data is to be put in plane.
 * 
 * @return  plane() = 2D array for storing completeA(:,:,kz).
 */
void Get_XY_plane(int N[], Array<complx,3> A, Array<complx,2> plane, int kz);



/** @brief Model initial shell spectrum.
 *
 * @param  N[]  The size of the arrays A.
 * @param  a  A parameter.
 * 
 * @return  \f$ Sk(k) = a k^4 exp(-0.02 k^1.1) / (k^4 + q^4)^(1+2.8/12) \f$ with q=1.5.
 */
void Model_initial_shell_spectrum(int N[],  Array<DP,1> Sk, DP a);


#endif

//*********************************   End of array_basic.h  ***********************************




