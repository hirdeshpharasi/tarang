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


/*! \file  read_prog_para.cc
 * 
 *	@brief  Reads prog prameters: prog_para_file, prog_kind, data_dir_name, D.
 *
 *	@author  M. K. Verma
 *	@version 4.0  MPI
 *	@date Sept 2008
 */


#include "main.h"


//*********************************************************************************************					

void Read_prog_para(ifstream& prog_para_file, string& prog_kind, string& data_dir_name)
{

	string s; 

	if (my_id == master_id) 
		cout << endl << "======= Reading program parameters =========" << endl << endl << endl;

	if (! prog_para_file.is_open())  
	{
		cout << "Unable to open prog_para_file: MY_ID  " << my_id <<  endl;
		exit(1);
	}

	//
	//

	getline(prog_para_file, s); getline(prog_para_file, s); getline(prog_para_file, s);
	getline(prog_para_file, s); getline(prog_para_file, s); getline(prog_para_file, s);

	prog_para_file >> prog_kind;
	if (my_id == master_id)		
		cout << "	Program: " << prog_kind << endl;
	getline(prog_para_file, s); getline(prog_para_file, s); getline(prog_para_file, s);

	prog_para_file >> data_dir_name;
	if (my_id == master_id)		
		cout << "	data_dir_name: " << data_dir_name << endl;

	cout  << "===========================================" << endl << endl;
	
}	


//********************************** End of read_prog_para..h *********************************




