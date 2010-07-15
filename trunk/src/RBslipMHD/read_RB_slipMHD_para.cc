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

/*! \file RB_slip_main.h
 * 
 * @brief Reads parameters related to RB convection (Pr, r, q/k0. eta, Pr_switch, and Uscaling).
 *
 * @author  M. K. Verma
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 

//********************************************************************************************* 


#include "../main.h"

extern int my_id;								// My process id
extern int numprocs;							// No of processors
extern const int master_id;						// Id of master proc
extern ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
extern ptrdiff_t local_N2;
extern MPI_Status status;


//********************************************************************************************* 

void Read_RB_slipMHD_para
(	
	ifstream& RB_slipMHD_para_file,  
	double& Pr, 
	double& r, 
	double& eta,		// magnetic diffusivity
	double& qbyk0, 
	string& Pr_switch, 
	string& Uscaling
)
{

	string s; 
	
	if (my_id == master_id) 
		cout  << "************** Reading RB parameters ******************" << endl << endl;	
		
	if (! RB_slip_para_file.is_open())  
	{
		cout << "Unable to open RB_slip_para_file;  Exiting Program: MY_ID " << my_id << endl;
		exit(1);
	}
	
	
	getline(RB_slip_para_file, s); getline(RB_slip_para_file, s); getline(RB_slip_para_file, s);
		
	RB_slip_para_file >> Pr >> r >> eta >> qbyk0;	
	
	if (my_id == master_id)		
		cout << "	Pr, r, eta. q/k0: " << Pr << " "  << r << " " 
			 << eta  << " " << qbyk0  << endl << endl;
		
	getline(RB_slip_para_file, s); getline(RB_slip_para_file, s); getline(RB_slip_para_file, s);
	
	RB_slip_para_file >>  Pr_switch >> Uscaling ;
	
	if (my_id == master_id)		
		cout << "	 Pr_switch & Uscaling : " << Pr_switch   << " "  << Uscaling << endl << endl;
	
	if (my_id == master_id) 
		cout  << "******************************************************" << endl << endl;
}	



//******************************** End of RB_slipMHD_para.cc ********************************** 



