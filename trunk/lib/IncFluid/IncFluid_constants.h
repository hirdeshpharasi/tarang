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


/*! \file IncFluid_constants.h
 * 
 * @brief Constants for Incompressible flows
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 

//*********************************************************************************************

#ifndef _H_IF_constants
#define _H_IF_constants

const int MAXSIZE_STRING_SWITCHES_ARRAY = 20;
const int MAXSIZE_SWITCHES_ARRAY = 20;

const int MAXSIZE_SOLVER_META_PARA = 10;
const int MAXSIZE_SOLVER_INT_PARA = 10;
const int MAXSIZE_SOLVER_DOUBLE_PARA = 10;
const int MAXSIZE_SOLVER_STRING_PARA = 10;

const int MAXSIZE_DIAGNOSTIC_PARA = 20;

const int MAXSIZE_TIME_PARA_ARRAY = 10;
const int MAXSIZE_TIME_SAVE_ARRAY = 30;

const int MAXSIZE_MISC_OUTPUT_PARA = 20;

const int MAXSIZE_OUT_K_R_ARRAY = 200;	
const int MAXSIZE_K_R_BUFFER = 20;                // for MPI data transfer of u(k) etc.


const int MAXSIZE_ET_PARAMETERS = 20;				// For energy transfer 
const int MAXSIZE_ET_RADII_SECTOR_ARRAY = 200;							

const int MAXSIZE_INIT_COND_META_PARA = 10;
const int MAXSIZE_INIT_COND_INT_PARA = 10;
const int MAXSIZE_INIT_COND_DOUBLE_PARA = 10;
const int MAXSIZE_INIT_COND_STRING_PARA = 10;

const int MAXSIZE_FORCE_META_PARA = 10;
const int MAXSIZE_FORCE_INT_PARA = 10;
const int MAXSIZE_FORCE_DOUBLE_PARA = 10;
const int MAXSIZE_FORCE_STRING_PARA = 10;



// When we input init conditions using modes.
const int INIT_COND_MAX_NO_FIELD_MODES = 51;
const int MAX_NO_VARIABLES = 10;				// at each k: e.g., Vx, Vy, Vz Wx, Wy, Wz

// When we provide force using modes.  The number of vars per mode is MAX_NO_VARIABLES.
const int MAX_NO_FORCE_FIELD_MODES = 51;

#endif

 
//********************************** End of IncFluid_constants.h ******************************


