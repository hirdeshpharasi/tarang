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

/*! \file RB_slipMHD_main.h
 * 
 * @brief Main program executing RB magnetoconvection under no-slip boundary condition.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 

//********************************************************************************************* 

 #include "../main.h"
// main.h contains IncFluid.h (all declarations)

void Read_RB_slipMHD_para
(	
	ifstream& RB_para_file,  
	double& Pr, 
	double& r, 
	double& eta, 
	double& qbyk0, 
	string& Pr_switch, 
	string& Uscaling
);
																					
		
//******************************** End of RB_slip_main.h **************************************


																								

