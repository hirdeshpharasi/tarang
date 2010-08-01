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

/*! \file  Output_main.cc
 * 
 * @brief  File operations
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 */	

#include "IncFluid.h"

/*
#include <sstream>

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

const H5std_string DATASET_REALFIELD_NAME1( "RealField1" );
const H5std_string DATASET_REALFIELD_NAME2( "RealField2" );
const H5std_string DATASET_REALFIELD_NAME3( "RealField3" );
const H5std_string DATASET_FIELD_NAME1R( "ComplexField1R" );
const H5std_string DATASET_FIELD_NAME1C( "ComplexField1C" );
const H5std_string DATASET_FIELD_NAME2R( "ComplexField2R" );
const H5std_string DATASET_FIELD_NAME2C( "ComplexField2C" );
const H5std_string DATASET_FIELD_NAME3R( "ComplexField3R" );
const H5std_string DATASET_FIELD_NAME3C( "ComplexField3C" );
*/
//*********************************************************************************************

/// If the input file field_in_file is not open, then error message.
void IncFluid::Input_prefix(ifstream& field_in_file)
{
	if (my_id == master_id) 
	{
		if (! field_in_file.is_open())  
		{
			cout << "Unable to open field_in_file: Exiting Program: Exiting Program " << endl;
			exit(1);
		}   
		
		else 
			cout << "Reading field configurations from field_in_file " << endl ;
	}
}

//*********************************************************************************************


/// If the input file force_field_in_file is not open, then error message.
void IncFluid::Input_force_prefix(ifstream& force_field_in_file)
{
	if (my_id == master_id) 
	{
		if (! force_field_in_file.is_open())  
		{
			cout << "MYERROR: Unable to open force_field_in_file: Exiting Program: Exiting Program" 
				 << endl;
			exit(1);
		}   
		
		else 
			cout << "Reading field configurations from force_field_in_file " << endl ;
	}		

}


//*********************************************************************************************

/// In file X, output prefix_str containing diffusive parameters, and
/// Grid size, basis & Integration schemes, Tinit, Tfinal, and Tdt.
void IncFluid::Output_prefix(ofstream& X, string prefix_str)
{
	if (my_id == master_id) 
	{
		X << prefix_str << endl;
		X << "%% Grid size " <<  N[1] << " "  << N[2] << " " << N[3] << endl;
		X << "%% Basis, Integration scheme:  " << basis_type << " " << integ_scheme << endl;
		X << "%% Tinit, Tfinal, Tdt_fixed:  " << Tinit << "    " << Tfinal 
			<< "    " << Tdt_fixed << endl << endl;
		
		X << "%% globalvar_anisotropy_switch: " << globalvar_anisotropy_switch << endl;
		X << "%% globalvar_waveno_switch: " << globalvar_waveno_switch << endl;
	}		
}

//*********************************************************************************************

std::string FloatToString(float number)
{
	std::ostringstream buff;
	buff<<number;
	return buff.str();
}

/// Open input file "field_in.d".
void IncFluid::Open_field_input_files()
{
	if (my_id == master_id) 
	{
		string filename;
		
		filename = "/in/field_in.d";
		filename = data_dir_name+ filename;   
		field_in_file.open(filename.c_str());
	}	
}	


//*********************************************************************************************
/// Open force input file (force wavenos).
void IncFluid::Open_force_waveno_input_files()
{
	if (my_id == master_id) 
	{
		string filename;
		
		filename = "/in/force_field_in.d";
		filename = data_dir_name+ filename;   
		force_field_in_file.open(filename.c_str());
	}	
}	

//*********************************************************************************************
/// Open output files
void IncFluid::Open_output_files(string prefix_str)
{
	if (my_id == master_id) 
	{
		string filename;

		filename = "/out/glob.d"; 
		filename = data_dir_name+ filename;   global_file.open(filename.c_str());
		
		filename = "/out/field_out.d";
		filename = data_dir_name+ filename;   field_out_file.open(filename.c_str());
		
		filename = "/out/field_out_reduced.d";
		filename = data_dir_name+ filename;   field_out_reduced_file.open(filename.c_str());
		
		filename = "/out/realfield_out.d";
		filename = data_dir_name+ filename;   realfield_out_file.open(filename.c_str());
		
		filename = "/out/field_k_out.d";
		filename = data_dir_name+ filename;	  field_k_out_file.open(filename.c_str());
		
		filename = "/out/field_r_out.d";
		filename = data_dir_name+ filename;	  field_r_out_file.open(filename.c_str());
		
		filename = "/out/spectrum.d";
		filename = data_dir_name+ filename;   spectrum_file.open(filename.c_str());
		
		filename = "/out/pressure.d";
		filename = data_dir_name+ filename;	  pressure_file.open(filename.c_str());
		
		filename = "/out/flux.d";
		filename = data_dir_name+ filename;   flux_file.open(filename.c_str());
		
		filename = "/out/shell.d";
		filename = data_dir_name+ filename;   shell_to_shell_file.open(filename.c_str());
		
		filename = "/out/ring_specctrum.d";
		filename = data_dir_name+ filename;	  ring_spectrum_file.open(filename.c_str());
		
		filename = "/out/ring_to_ring.d";
		filename = data_dir_name+ filename;	  ring_to_ring_file.open(filename.c_str());
		
		filename = "/out/cyl_ring_spectrum.d";
		filename = data_dir_name+ filename;	  cylinder_ring_spectrum_file.open(filename.c_str());
		
		filename = "/out/cyl_ring_to_ring.d";
		filename = data_dir_name+ filename;	  cylinder_ring_to_ring_file.open(filename.c_str());
		
		filename = "/out/struct_fn.d";
		filename = data_dir_name+ filename;	  structure_fn_file.open(filename.c_str());
	
		filename = "/out/planar_structure_fn.d";
		filename = data_dir_name+ filename;	  planar_structure_fn_file.open(filename.c_str());
		
		filename = "/out/moment.d";
		filename = data_dir_name+ filename;	  moment_file.open(filename.c_str());
		
		filename = "/out/Skpq.d";
		filename = data_dir_name+ filename;	  Skpq_file.open(filename.c_str());
		
		// Set precision digit
		global_file.setf(ios::fixed, ios::floatfield);
		global_file.precision(MY_PRECISION); // sets number of decimal places
		
		field_k_out_file.setf(ios::fixed, ios::floatfield);
		field_k_out_file.precision(MY_PRECISION); 
		
		field_r_out_file.setf(ios::fixed, ios::floatfield);
		field_r_out_file.precision(MY_PRECISION); 
		
		cout.setf(ios::fixed, ios::floatfield);
		cout.precision(MY_PRECISION); 
		
		
		filename = "/out/readme"; 
		filename = data_dir_name+ filename;   readme_file.open(filename.c_str());
		Output_prefix(readme_file, prefix_str);
	
		// output size and other information to the file
		
					
		//	Output_prefix(field_out_file, prefix_str);	
		// prefix_str removed to make reading of the file easier.
		
		/*
		Output_prefix(global_file, prefix_str);	
		Output_prefix(field_out_reduced_file, prefix_str);
		Output_prefix(realfield_out_file, prefix_str);	
		Output_prefix(field_k_out_file, prefix_str);	
		Output_prefix(field_r_out_file, prefix_str);
		Output_prefix(spectrum_file, prefix_str);	
		Output_prefix(pressure_file, prefix_str);
		Output_prefix(flux_file, prefix_str);	
		*/

		// shell-to-shell 
		Output_prefix(shell_to_shell_file, prefix_str);		
		
		shell_to_shell_file << "%% Shellradius = " 
				<< (*shell_radius)(Range(1,toEnd)) << endl;
		
		// ring-spectrum
		
		if (CV_anisotropic_ring_switch == 1)
		{
		// Output_prefix(ring_spectrum_file, prefix_str);
			ring_spectrum_file << "%% Sectors = " << (*sector_angle_array_spectrum) << endl;
		}	
		
		// ring-to-ring transfer
		
		if (ET_anisotropic_ring_switch == 1)
		{
		//	Output_prefix(ring_to_ring_file, prefix_str);	
			ring_to_ring_file << "%% Ring radius = " 
					<< (*ring_shell_radius)(Range(1, no_ring_shells-1)) << endl;
					
			ring_to_ring_file << "%% Sectors = " << (*sector_angle_ring_tr) << endl;
		}
		// cylinderical-ring-spectrum
		
		if (CV_anisotropic_cylinder_switch == 1)
		{
		//  Output_prefix(cylinder_ring_spectrum_file, prefix_str);
			cylinder_ring_spectrum_file << "%% Slabs = " 
					<< (*cylinder_kpll_array_spectrum) << endl;
		}	
		
		// cylinderical-ring-to-ring transfer
		
		if (ET_anisotropic_cylinder_switch == 1)
		{
		//  Output_prefix(cylinder_ring_to_ring_file, prefix_str);	
			cylinder_ring_to_ring_file << "%% Ring radius = " 
					<< (*cylinder_shell_radius)(Range(1, no_cylinder_shells-1)) << endl;
					
			cylinder_ring_to_ring_file << "%% Slabs = " << (*cylinder_kpll_array_tr) << endl;
		}
		//
		
		/*
		Output_prefix(structure_fn_file, prefix_str);	
		Output_prefix(planar_structure_fn_file, prefix_str);
		Output_prefix(moment_file, prefix_str);	
		Output_prefix(Skpq_file, prefix_str);
		 */
	}
}

/*
/// Open output files in hdf5 format
void IncFluid::Open_output_files_hdf5(string prefix_str, float timeNow, int N[], int rank)
{
	if (my_id == master_id)
	{
		filename_hdf5 = "/out/hdf/Out"+ FloatToString(timeNow)+ ".h5";
		filename_hdf5 = data_dir_name+ filename_hdf5;
         	file = new H5File( filename_hdf5, H5F_ACC_TRUNC );

		if (Tnow >= Trealfield_save_next) 
		{
         		//
         		//Create property list for a dataset and set up fill values.
         		//
         		float fillvalue = 0.0;   // Fill value for the dataset
         		DSetCreatPropList plist;
         		plist.setFillValue(PredType::NATIVE_FLOAT, &fillvalue);
	        	
         		
			    //
          		// Create dataspace for the dataset in the file.
          		//
         		hsize_t realfield_dim[] = {N[1], N[2], N[3]}; // dim sizes of ds (on disk)
         		DataSpace realfield_space( rank, realfield_dim );
 		
         		//
          		// Create dataset and write it into the file.
          		//
         		realfield_dataset1 = new DataSet(file->createDataSet( DATASET_REALFIELD_NAME1,
						PredType::NATIVE_FLOAT, realfield_space, plist));
         		realfield_dataset2 = new DataSet(file->createDataSet( DATASET_REALFIELD_NAME2,
						PredType::NATIVE_FLOAT, realfield_space, plist));
         		realfield_dataset3 = new DataSet(file->createDataSet( DATASET_REALFIELD_NAME3,
						PredType::NATIVE_FLOAT, realfield_space, plist));
		}
		if (Tnow >= Tfield_save_next) 
		{
         		//
         		// Create property list for a dataset and set up fill values.
         		//
         		float fillvalue = 0.0;   // Fill value for the dataset
         		DSetCreatPropList plist;
         		plist.setFillValue(PredType::NATIVE_FLOAT, &fillvalue);
	        	
         		
			    //
          		// Create dataspace for the dataset in the file.
          		//
         		hsize_t field_dim[] = {N[1], N[2], (N[3]/2 + 1)}; // dim sizes of ds (on disk)
         		DataSpace field_space( rank, field_dim );
 		
         		//
          		// Create dataset and write it into the file.
          		//
         		field_dataset1R = new DataSet(file->createDataSet( DATASET_FIELD_NAME1R,
						PredType::NATIVE_FLOAT, field_space, plist));
         		field_dataset1C = new DataSet(file->createDataSet( DATASET_FIELD_NAME1C,
						PredType::NATIVE_FLOAT, field_space, plist));
         		field_dataset2R = new DataSet(file->createDataSet( DATASET_FIELD_NAME2R,
						PredType::NATIVE_FLOAT, field_space, plist));
         		field_dataset2C = new DataSet(file->createDataSet( DATASET_FIELD_NAME2C,
						PredType::NATIVE_FLOAT, field_space, plist));
         		field_dataset3R = new DataSet(file->createDataSet( DATASET_FIELD_NAME3R,
						PredType::NATIVE_FLOAT, field_space, plist));
         		field_dataset3C = new DataSet(file->createDataSet( DATASET_FIELD_NAME3C,
						PredType::NATIVE_FLOAT, field_space, plist));
		}
	}
}

/// Close output files in hdf5 format
void IncFluid::Close_output_files_hdf5(string key)
{
	if (my_id == master_id) 
	{
		if(key == "realField")
		{
			delete realfield_dataset1;
			delete realfield_dataset2;
			delete realfield_dataset3;
		}
		if(key == "complexField")
		{
			delete field_dataset1R;
			delete field_dataset1C;
			delete field_dataset2R;
			delete field_dataset2C;
			delete field_dataset3R;
			delete field_dataset3C;
		}
		if(key == "final")
		{
			delete file;
		}
	}
}
*/

//*********************************************************************************************
/// Close input file.
void IncFluid::Close_field_input_files()
{
	if (my_id == master_id) 
		field_in_file.close();	
}

//*********************************************************************************************
/// Close force_field input file.
void IncFluid::Close_force_waveno_input_files()
{
	if (my_id == master_id) 
		force_field_in_file.close();	
}

//*********************************************************************************************
/// Close output files.
void IncFluid::Close_output_files()
{	
	if (my_id == master_id) 
	{
		global_file.close(); 
		field_out_file.close(); 
		field_out_reduced_file.close();
		realfield_out_file.close(); 
		field_k_out_file.close();
		field_r_out_file.close();
		spectrum_file.close(); 
		pressure_file.close();
		flux_file.close(); 
		shell_to_shell_file.close(); 
		ring_spectrum_file.close();
		ring_to_ring_file.close();
		cylinder_ring_spectrum_file.close();
		cylinder_ring_to_ring_file.close();
		structure_fn_file.close();
		planar_structure_fn_file.close();
		moment_file.close();
		Skpq_file.close();
	}	
}


//********************************** End of file_operations.cc ********************************

