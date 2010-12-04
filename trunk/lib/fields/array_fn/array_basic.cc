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


/*! \file array_basic.cc
 * 
 * @sa field_basic.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */ 
 
 
#include "array_basic.h"

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
/**********************************************************************************************

     		C = A.B (assuming A and B are real)
			Term by term multiplication.

***********************************************************************************************/

void Array_real_mult
(
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<complx,3> C
)
{
	real(C) = real(A) * real(B);
	imag(C) = imag(A) * imag(B);
}


/**********************************************************************************************

     		Output A assuiming it to be real.

***********************************************************************************************/ 

void Output_asreal
(
	ofstream& file_out, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_array
)
{

	int source_local_N1, source_local_N1_start;			// that of worker
	int data_size;
	int tag = 123;
	
	if (my_id == master_id) 
	{
		file_out << "%% Array " <<  N[1] << " x "  << N[2] << " x " << N[3] << endl;
	
		for (int i=0; i<local_N1; i++) 
		{
			for (int j=0; j<N[2]; j++) 
				for (int k=0; k<N[3]/2; k++) 
					file_out << real(A(i,j,k)) << " " << imag(A(i,j,k)) << " " ;
			file_out << endl; 
		}
			
		for (int source = 1; source <= numprocs-1; source++) 
		{
			MPI_Recv( &source_local_N1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			MPI_Recv( &source_local_N1_start, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			MPI_Recv( &data_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			
			MPI_Recv( reinterpret_cast<double*>((temp_array).data()), data_size, 
						MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );	
			
			for (int i=0; i<local_N1; i++) 
			{
				for (int j=0; j<N[2]; j++) 
					for (int k=0; k<N[3]/2; k++) 
						file_out << real(temp_array(i,j,k)) << " " << imag(temp_array(i,j,k)) 
								 << " " ;
				file_out << endl; 
			}
		}
		
		file_out << endl << endl;
	}
	
	else		// process is not master
	{
		MPI_Send( &local_N1, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		MPI_Send( &local_N1_start, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		data_size = 2* local_N1 * N[2] * (N[3]/2 + 1);
		MPI_Send( &data_size, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		
		MPI_Send( reinterpret_cast<double*>(A.data()), data_size, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );								
	}						
	
	MPI_Barrier(MPI_COMM_WORLD);


}


void Output_asreal_transpose_order
(
	 ofstream& file_out, 
	 int N[], 
	 Array<complx,3> A, 
	 Array<complx,3> temp_r
)
{
	
	int source_local_N2, source_local_N2_start;			// that of worker
	int data_size;
	int tag = 123;
	
	if (my_id == master_id) 
	{
		file_out << "%% Array " <<  N[2] << " x "  << N[1] << " x " << N[3] << endl;
		
		for (int i=0; i<local_N2; i++) 
		{
			for (int j=0; j<N[1]; j++) 
				for (int k=0; k<N[3]/2; k++) 
					file_out << real(A(i,j,k)) << " " << imag(A(i,j,k)) << " " ;
			file_out << endl; 
		}
		
		for (int source = 1; source <= numprocs-1; source++) 
		{
			MPI_Recv( &source_local_N2, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			MPI_Recv( &source_local_N2_start, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			MPI_Recv( &data_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			
			MPI_Recv( reinterpret_cast<double*>((temp_r).data()), data_size, 
					 MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );	
			
			for (int i=0; i<local_N2; i++) 
			{
				for (int j=0; j<N[1]; j++) 
					for (int k=0; k<N[3]/2; k++) 
						file_out << real(temp_r(i,j,k)) << " " << imag(temp_r(i,j,k)) 
						<< " " ;
				file_out << endl; 
			}
		}
		
		file_out << endl << endl;
	}
	
	else		// process is not master
	{
		MPI_Send( &local_N2, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		MPI_Send( &local_N2_start, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		data_size = 2* local_N2 * N[1] * (N[3]/2 + 1);
		MPI_Send( &data_size, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		
		MPI_Send( reinterpret_cast<double*>(A.data()), data_size, 
				 MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );								
	}						
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
}

/*
void Output_asreal_hdf5
(
 	DataSet* dataset,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_array
)
{

	int source_local_N1, source_local_N1_start;			// that of worker
	int data_size;
	int tag = 123;
	

	if (my_id == master_id) 
	{
		hsize_t start[3]; //Start of hyperslab
		hsize_t stride[3]; //Stride of hyperslab
		hsize_t count[3]; //Block count
		hsize_t block[3]; //Block sizes
		start[1] = 0; start[2] = 0;
		stride[0] = 1; stride[1] = 1; stride[2] = 1;
		count[0] = 1; count[1] = 1; count[2] = 1;
		block[0] = local_N1; block[1] = N[2]; block[2] = N[3];
		
		start[0] = 0;

        hsize_t realfield_dim[] = {N[1], N[2], N[3]}; // dim sizes of ds (on disk)
        DataSpace fspace( 3, realfield_dim );

		hsize_t dim[] = {local_N1, N[2], N[3]};
		DataSpace mspace(3, dim);

		fspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		mspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

		dataset->write(real(A).data(), PredType::NATIVE_FLOAT, mspace, fspace);

//		file_out << "%% Array " <<  N[1] << " x "  << N[2] << " x " << N[3] << endl;
//	
//		for (int i=0; i<local_N1; i++) 
//		{
//			for (int j=0; j<N[2]; j++) 
//				for (int k=0; k<N[3]/2; k++) 
//					file_out << real(A(i,j,k)) << " " << imag(A(i,j,k)) << " " ;
//			file_out << endl; 
//		}
			
		for (int source = 1; source <= numprocs-1; source++) 
		{
			MPI_Recv( &source_local_N1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			MPI_Recv( &source_local_N1_start, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			MPI_Recv( &data_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			
			MPI_Recv( reinterpret_cast<double*>((temp_array).data()), data_size, 
						MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );	
			
			start[0] = source * local_N1;
			fspace.selectNone();
			fspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
			start[0] = 0; // Is mspace selection needed again?
			mspace.selectNone();
			mspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

			dataset->write(real(temp_array).data(), PredType::NATIVE_FLOAT, mspace, fspace);
//
//			for (int i=0; i<local_N1; i++) 
//    		{
//				for (int j=0; j<N[2]; j++) 
//					for (int k=0; k<N[3]/2; k++) 
//						file_out << real(temp_array(i,j,k)) << " " << imag(temp_array(i,j,k)) 
//								 << " " ;
//				file_out << endl; 
//			}

		}
		
//		file_out << endl << endl;

	}
	
	else		// process is not master
	{
		MPI_Send( &local_N1, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		MPI_Send( &local_N1_start, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		data_size = 2* local_N1 * N[2] * (N[3]/2 + 1);
		MPI_Send( &data_size, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		
		MPI_Send( reinterpret_cast<double*>(A.data()), data_size, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );								
	}						
	
	MPI_Barrier(MPI_COMM_WORLD);

}
*/

/**********************************************************************************************

		Master reads data from file_in and distributes to A in all procs 

***********************************************************************************************/ 

void Read_data_MPI
(
	ifstream& file_in, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_array
)
{

	int dest_local_N1[numprocs];					// that of worker
	int dest_local_N1_start[numprocs];
	int dest_data_size[numprocs];
	
	int temp1, temp2, temp3;
	
	int data_size;
	int tag1 = 111;
	int tag2 = 222;
	int tag3 = 333;
	
	//
	// First collect the local_N1, local_N1_start, and data_size from the destinations
	//
	if (my_id == master_id) 
	{
		for (int dest = 1; dest <= numprocs-1; dest++) 
		{
			MPI_Recv( &temp1, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD, &status );
			MPI_Recv( &temp2, 1, MPI_INT, dest, tag2, MPI_COMM_WORLD, &status );
			MPI_Recv( &temp3, 1, MPI_INT, dest, tag3, MPI_COMM_WORLD, &status );	
			
			dest_local_N1[dest] = temp1;
			dest_local_N1_start[dest] = temp2;
			dest_data_size[dest] =  temp3;														
		}
	}
	
	else		// process is not master
	{
		MPI_Send( &local_N1, 1, MPI_INT, master_id, tag1, MPI_COMM_WORLD );
		MPI_Send( &local_N1_start, 1, MPI_INT, master_id, tag2, MPI_COMM_WORLD );
		data_size =  2* local_N1 * N[2] * (N[3]/2 + 1);
		MPI_Send( &data_size, 1, MPI_INT, master_id, tag3, MPI_COMM_WORLD );
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//
	// Now master reads the data and sends to all procs.
	//
	
	if (my_id == master_id) 
	{
		for (int i=0; i<local_N1; i++)
			for (int j=0; j<N[2]; j++)
				for (int k=0; k<=N[3]/2; k++)
					file_in >> A(i,j,k);
				
		
		for (int dest = 1; dest <= numprocs-1; dest++) 
		{
			for (int i=0; i<local_N1; i++)
				for (int j=0; j<N[2]; j++)
					for (int k=0; k<=N[3]/2; k++)
						file_in >> temp_array(i,j,k);
		
			data_size = dest_data_size[dest];			
			MPI_Send( reinterpret_cast<double*>((temp_array).data()), data_size, 
						MPI_DOUBLE, dest, tag1, MPI_COMM_WORLD );																							
		}
	}
	
	else		// process is not master
	{
		MPI_Recv( reinterpret_cast<double*>(A.data()), data_size, 
						MPI_DOUBLE, master_id, tag1, MPI_COMM_WORLD, &status);	
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
}



/**********************************************************************************************

		Master reads data of REDUCED SIZE from file_in and distributes to A in all procs 

***********************************************************************************************/ 

void Read_data_MPI
(
	ifstream& file_in, 
	int N[], 
	int Nreduced[], 
	Array<complx,3> A,  
	Array<complx,3> temp_array
)
{
	
	int dest_local_N1[numprocs];					// that of worker
	int dest_local_N1_start[numprocs];
	int dest_data_size[numprocs];
	
	int temp1, temp2, temp3;
	
	int data_size;
	int tag1 = 111;
	int tag2 = 222;
	int tag3 = 333;
	
	//
	// First collect the local_N1, local_N1_start, and data_size from the destinations
	//
	if (my_id == master_id) 
	{
		for (int dest = 1; dest <= numprocs-1; dest++) 
		{
			MPI_Recv( &temp1, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD, &status );
			MPI_Recv( &temp2, 1, MPI_INT, dest, tag2, MPI_COMM_WORLD, &status );
			MPI_Recv( &temp3, 1, MPI_INT, dest, tag3, MPI_COMM_WORLD, &status );	
			
			dest_local_N1[dest] = temp1;
			dest_local_N1_start[dest] = temp2;
			dest_data_size[dest] =  temp3;														
		}
	}
	
	else		// process is not master
	{
		MPI_Send( &local_N1, 1, MPI_INT, master_id, tag1, MPI_COMM_WORLD );
		MPI_Send( &local_N1_start, 1, MPI_INT, master_id, tag2, MPI_COMM_WORLD );
		data_size =  2* local_N1 * N[2] * (N[3]/2 + 1);
		MPI_Send( &data_size, 1, MPI_INT, master_id, tag3, MPI_COMM_WORLD );
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	//
	// Now send the data
	//
	
	int last_dest;					// The last processor that will contain reducedA.
	int	remaining_i;
	
		
	if (Nreduced[1] > local_N1)
	
		if ((Nreduced[1] % local_N1) == 0)
			last_dest = Nreduced[1]/local_N1 -1;
		else
			last_dest = Nreduced[1]/local_N1;
			
	else
		last_dest = 0;		
		
						
	A = 0.0;		// Initialize for every proc
	
	if (my_id == master_id) 
	{
		int last_i = minimum(local_N1, Nreduced[1]) - 1;			
		// lower value of local_N1-1 or Nreducec[1]-1.
		
		for (int i=0; i<=last_i; i++)
			for (int j=0; j<Nreduced[2]; j++) 
				for (int k=0; k<=Nreduced[3]/2; k++) 
					file_in >> A(i,j,k);
	
		for (int dest = 1; dest <= last_dest; dest++) 
		{
			remaining_i = Nreduced[1] - dest*local_N1;
			last_i = minimum(local_N1, remaining_i) - 1;
		
			temp_array = 0.0;
			
			for (int i=0; i<=last_i; i++)
				for (int j=0; j<Nreduced[2]; j++) 
					for (int k=0; k<=Nreduced[3]/2; k++) 
						file_in >> temp_array(i,j,k);
			
			data_size = dest_data_size[dest];			
			MPI_Send( reinterpret_cast<double*>((temp_array).data()), data_size, 
							MPI_DOUBLE, dest, tag1, MPI_COMM_WORLD );	
		}
	}	
				
	else if ((my_id > 0) && (my_id <= last_dest))			// process is not master
	{
		MPI_Recv( reinterpret_cast<double*>(A.data()), data_size, 
						MPI_DOUBLE, master_id, tag1, MPI_COMM_WORLD, &status);	
	}
	
	else
		A = 0.0;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
}


/**********************************************************************************************

		Master colllects  A  from all procs and writes in file_out 

***********************************************************************************************/ 
 

void Write_data_MPI
(
	ofstream& file_out, 
	int N[], 
	Array<complx,3> A,  
	Array<complx,3> temp_array
)
{

	int source_local_N1, source_local_N1_start;			// that of worker
	int data_size;
	int tag = 123;
	
	if (my_id == master_id) 
	{
		for (int i=0; i<local_N1; i++) 
			for (int j=0; j<N[2]; j++) 
				for (int k=0; k<=N[3]/2; k++) 
					file_out << A(i,j,k) << " " ;
		
		for (int source = 1; source <= numprocs-1; source++) 
		{
			MPI_Recv( &source_local_N1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			MPI_Recv( &source_local_N1_start, 1, MPI_INT, source, tag, MPI_COMM_WORLD,&status);
			MPI_Recv( &data_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			
			MPI_Recv( reinterpret_cast<double*>((temp_array).data()), data_size, 
						MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );	
		
			for (int i=0; i<local_N1; i++) 
				for (int j=0; j<N[2]; j++) 
					for (int k=0; k<=N[3]/2; k++) 
						file_out << temp_array(i,j,k) << " " ;				
		}
		
		file_out << endl << endl;
	}
	
	else		// process is not master
	{
		MPI_Send( &local_N1, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		MPI_Send( &local_N1_start, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		data_size = 2* local_N1 * N[2] * (N[3]/2 + 1);
		MPI_Send( &data_size, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		
		MPI_Send( reinterpret_cast<double*>(A.data()), data_size, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );								
	}						
	
	MPI_Barrier(MPI_COMM_WORLD);
	
}

/*
void Write_data_MPI_hdf5
(
	DataSet* datasetR,
    DataSet* datasetC,	
	int N[], 
	Array<complx,3> A,  
	Array<complx,3> temp_array
)
{

	int source_local_N1, source_local_N1_start;			// that of worker
	int data_size;
	int tag = 123;
	
	if (my_id == master_id) 
	{
		hsize_t start[3]; //Start of hyperslab
		hsize_t stride[3]; //Stride of hyperslab
		hsize_t count[3]; //Block count
		hsize_t block[3]; //Block sizes
		start[1] = 0; start[2] = 0;
		stride[0] = 1; stride[1] = 1; stride[2] = 1;
		count[0] = 1; count[1] = 1; count[2] = 1;
		block[0] = local_N1; block[1] = N[2]; block[2] = (N[3]/2 + 1);
		
		start[0] = 0;

         	hsize_t field_dim[] = {N[1], N[2], (N[3]/2 + 1)}; // dim sizes of ds (on disk)
         	DataSpace fspace( 3, field_dim );

		hsize_t dim[] = {local_N1, N[2], (N[3]/2 + 1)};
		DataSpace mspace(3, dim);

		fspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
		mspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

		datasetR->write(real(A).data(), PredType::NATIVE_FLOAT, mspace, fspace);
		datasetC->write(imag(A).data(), PredType::NATIVE_FLOAT, mspace, fspace);

//		for (int i=0; i<local_N1; i++) 
//			for (int j=0; j<N[2]; j++) 
//				for (int k=0; k<=N[3]/2; k++) 
//					file_out << A(i,j,k) << " " ;
//		
		for (int source = 1; source <= numprocs-1; source++) 
		{
			MPI_Recv( &source_local_N1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			MPI_Recv( &source_local_N1_start, 1, MPI_INT, source, tag, MPI_COMM_WORLD,&status);
			MPI_Recv( &data_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			
			MPI_Recv( reinterpret_cast<double*>((temp_array).data()), data_size, 
						MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );	

			start[0] = source * local_N1;
			fspace.selectNone();
			fspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
			start[0] = 0; // Is mspace selection needed again?
			mspace.selectNone();
			mspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

			datasetR->write(real(temp_array).data(), PredType::NATIVE_FLOAT, mspace, fspace);
			datasetC->write(imag(temp_array).data(), PredType::NATIVE_FLOAT, mspace, fspace);
		
//			for (int i=0; i<local_N1; i++) 
//				for (int j=0; j<N[2]; j++) 
//					for (int k=0; k<=N[3]/2; k++) 
//						file_out << temp_array(i,j,k) << " " ;				
		}
		
//		file_out << endl << endl;
	}
	
	else		// process is not master
	{
		MPI_Send( &local_N1, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		MPI_Send( &local_N1_start, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		data_size = 2* local_N1 * N[2] * (N[3]/2 + 1);
		MPI_Send( &data_size, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		
		MPI_Send( reinterpret_cast<double*>(A.data()), data_size, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );								
	}						
	
	MPI_Barrier(MPI_COMM_WORLD);
	
}
*/
/**********************************************************************************************

		Master outputs data of REDUCED SIZE into file_out
		-- We assume that the Nreduced[1] > local_N1, so the Reduced data is to be
				fit in the master node.

***********************************************************************************************/ 


void Write_data_MPI
(
	ofstream& file_out, 
	int N[], 
	int Nreduced[], 
	Array<complx,3> A, 
	Array<complx,3> temp_array
)
{

	int source_local_N1, source_local_N1_start;			// that of worker
	int data_size;
	int last_source;
	int remaining_i;
	int tag = 123;

	
	if (Nreduced[1] > local_N1)
	
		if ((Nreduced[1] % local_N1) == 0)
			last_source = Nreduced[1]/local_N1 -1;
		else
			last_source = Nreduced[1]/local_N1;
			
	else
		last_source = 0;
	
		
				
	if (my_id == master_id) 
	{
		
		int last_i = minimum(local_N1, Nreduced[1]) - 1;
		
		for (int i=0; i<=last_i; i++)
			for (int j=0; j<Nreduced[2]; j++) 
				for (int k=0; k<=Nreduced[3]/2; k++) 
					file_out << A(i,j,k);
				
		for (int source = 1; source <= last_source; source++) 
		{
			MPI_Recv( &source_local_N1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			MPI_Recv( &source_local_N1_start, 1, MPI_INT, source, tag, MPI_COMM_WORLD,&status);
			MPI_Recv( &data_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			
			MPI_Recv( reinterpret_cast<double*>((temp_array).data()), data_size, 
						MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );	
			
			
			remaining_i = Nreduced[1] - source*local_N1;
			last_i = minimum(local_N1, remaining_i) - 1;
			
			for (int i=0; i<=last_i; i++)
				for (int j=0; j<Nreduced[2]; j++) 
					for (int k=0; k<=Nreduced[3]/2; k++) 
						file_out << temp_array(i,j,k);
		}
		
		file_out << endl << endl;
	}	
	
	else		// process is not master
	{
		if ((my_id > 0) && (my_id <= last_source))
		{
			MPI_Send( &local_N1, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
			MPI_Send( &local_N1_start, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
			data_size = 2* local_N1 * N[2] * (N[3]/2 + 1);
			MPI_Send( &data_size, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
			
			MPI_Send( reinterpret_cast<double*>(A.data()), data_size, 
							MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );
		}												
	}						
	
	MPI_Barrier(MPI_COMM_WORLD);
	
}


/**********************************************************************************************

		Collect pieced of kz planes from all the procs and put them in plane

***********************************************************************************************/ 


void Get_XY_plane(int N[], Array<complx,3> A, Array<complx,2> plane, int kz)
{
	
	int source_local_N1, source_local_N1_start;			// that of worker
	int data_size;
	int tag = 999;
	
	static Array<complx,2> plane_part(local_N1, N[2]);
	
	if (my_id == master_id) 
	{
		plane( Range(0,local_N1-1), Range::all() ) = A(Range::all(), Range::all(), kz);	

		for (int source = 1; source <= numprocs-1; source++) 
		{
			MPI_Recv( &source_local_N1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			MPI_Recv( &source_local_N1_start, 1, MPI_INT, source, tag, MPI_COMM_WORLD,&status);
			MPI_Recv( &data_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
			
			MPI_Recv( reinterpret_cast<double*>(plane_part.data()), data_size, 
						MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );	
									
			plane( Range(source_local_N1_start, source_local_N1_start+source_local_N1-1), 
							Range::all(), kz )  = plane_part;	
		}
	}
	else		// process is not master
	{
		MPI_Send( &local_N1, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		MPI_Send( &local_N1_start, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
		data_size = 2* local_N1 * N[2];
		MPI_Send( &data_size, 1, MPI_INT, master_id, tag, MPI_COMM_WORLD );
				
		plane_part = A(Range::all(), Range::all(), kz);		
		MPI_Send( reinterpret_cast<double*>(plane_part.data()), data_size, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD );								
	}						
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Now master sends the plane to all the nodes
	
	tag = 899;
	data_size = 2 * N[1] * N[2];
	
	if (my_id == master_id)		
		for (int dest = 1; dest <= numprocs-1; dest++) 
			MPI_Send( reinterpret_cast<double*>(plane.data()), data_size, 
						MPI_DOUBLE, dest, tag, MPI_COMM_WORLD );
							
	else	
		MPI_Recv( reinterpret_cast<double*>(plane.data()), data_size, 
						MPI_DOUBLE, master_id, tag, MPI_COMM_WORLD, &status );
						
	
	MPI_Barrier(MPI_COMM_WORLD);
													
}


/**********************************************************************************************

		Model energy spectrum for initial condition
		Sk(k) = a k^4 exp(-b k^1.1) / (k^4 + q^4)^(1+2.8/12)
		with q = 1.5, b = 0.02 

***********************************************************************************************/ 


void Model_initial_shell_spectrum(int N[],  Array<DP,1> Sk, DP a)
{
	
	DP q = 1.5;
	DP b = 0.02;
	DP c = 1+2.8/12;
	DP k;	// radius
	
//	TinyVector<int, 1> len = Sk.length();
	
	Sk = 0.0;						
									
	for (int i=0; i < (Sk.length())(0); i++)
	{	
		k = 1.0*i;
		Sk(i) = a * pow(k,4) * exp(- b * pow(k,1.1)) /pow( (pow(k,4)+pow(q,4)), c);
    }

}



//*********************************   End of array_basic.cc ***********************************





