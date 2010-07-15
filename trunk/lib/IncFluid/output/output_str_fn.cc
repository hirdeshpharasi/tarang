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

/*! \file  Output_str_fn.cc
 * 
 * @brief  Output Structure function, Planar structure function.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"   

//*********************************************************************************************


void IncFluid::Output_structure_fn()
{
	if (my_id == master_id)
		structure_fn_file << "%% Time = " << Tnow << endl; 
		
	CV_Compute_structure_function();
	
	if (my_id == master_id)
		for (int r=1; r<=CV_structurefn_r_max; r++)
		{
			for (int q=CV_structurefn_q_min; q<=CV_structurefn_q_max; q++)
				structure_fn_file	<< r				<< "  " 
									<< (*CV_St)(r,q,0)	<< " " 
									<< (*CV_St)(r,q,1);
								
			structure_fn_file << endl;
		}	
}
  
  
//
//

void IncFluid::Output_structure_fn(IncSF& T)
{
/*
	if (my_id == master_id)
		structure_fn_file << "%% Time = " << Tnow << endl; 
		
	CV_Compute_structure_function();
	T.CS_Compute_structure_function();
	
	if (my_id == master_id)
		for (int r=1; r<=CV_structurefn_r_max; r++)
		{
			for (int q=CV_structurefn_q_min; q<=CV_structurefn_q_max; q++)
				structure_fn_file	<< r					<< "  " 
									<< (*CV_St)(r,q,0)		<< " " 
									<< (*CV_St)(r,q,1)		<< "   "
									<< (*T.CS_St)(r,q);
								
			structure_fn_file << endl;
		}	
	*/	
}
    
//
//

void IncFluid::Output_structure_fn(IncVF& W)
{
	if (my_id == master_id)
		structure_fn_file << "%% Time = " << Tnow << endl; 
		
	CV_Compute_structure_function();
	W.CV_Compute_structure_function();
	
	if (my_id == master_id)
		for (int r=1; r<=CV_structurefn_r_max; r++)
		{
			for (int q=CV_structurefn_q_min; q<=CV_structurefn_q_max; q++)
				structure_fn_file	<< r					<< "  " 
									<< (*CV_St)(r,q,0)		<< " " 
									<< (*CV_St)(r,q,1)		<< "    "
									<< (*W.CV_St)(r,q,0)	<< " " 
									<< (*W.CV_St)(r,q,1);
								
			structure_fn_file << endl;
		}		
} 

//
//

void IncFluid::Output_structure_fn(IncVF& W, IncSF& T)
{
/*
	if (my_id == master_id)
		structure_fn_file << "%% Time = " << Tnow << endl; 
		
	CV_Compute_structure_function();
	W.CV_Compute_structure_function();
	T.CS_Compute_structure_function();
	
	
	if (my_id == master_id)
		for (int r=1; r<=CV_structurefn_r_max; r++)
		{
			for (int q=CV_structurefn_q_min; q<=CV_structurefn_q_max; q++)
				structure_fn_file	<< r					<< "  " 
									<< (*CV_St)(r,q,0)		<< " " 
									<< (*CV_St)(r,q,1)		<< "    "
									<< (*W.CV_St)(r,q,0)	<< " " 
									<< (*W.CV_St)(r,q,1)	<< "   "
									<< (*T.CS_St)(r,q);
								
			structure_fn_file << endl;
		}	
 */		
} 

//
//

void IncFluid::Output_structure_fn(IncSF& T, string Pr_switch)
{
	if (Pr_switch == "PRZERO")
		Output_structure_fn();
	
	else
		Output_structure_fn(T);
}

void IncFluid::Output_structure_fn(IncVF& W, IncSF& T, string Pr_switch)
{
	if (Pr_switch == "PRZERO")
		Output_structure_fn(W);
	
	else
		Output_structure_fn(W, T);
}


//*********************************************************************************************

// Planar structure function
//
void IncFluid::Output_planar_structure_fn()
{
	if (my_id == master_id)
		planar_structure_fn_file << "%% Time = " << Tnow << endl; 
		
	CV_Compute_planar_structure_function();
	
	if (my_id == master_id)
		planar_structure_fn_file	<< "%% Planar struncture function:  V " << endl 
								<< (*CV_st_planar)		<< endl << endl;
}
  
  
//
//

void IncFluid::Output_planar_structure_fn(IncSF& T)
{
/*
	if (my_id == master_id)
		planar_structure_fn_file << "%% Time = " << Tnow << endl; 
		
	CV_Compute_planar_structure_function();
	T.CS_Compute_planar_structure_function();
	
	if (my_id == master_id)
	{
		planar_structure_fn_file	<< "%% Planar struncture function:  V " << endl 
									<< (*CV_st_planar)		<< endl << endl;
						
		planar_structure_fn_file	<< "%% Planar struncture function: Scalar field T " << endl 
									<< (*T.CS_st_planar)	<< endl << endl;	
	}	
*/											
}

//
//

void IncFluid::Output_planar_structure_fn(IncVF& W)
{
	if (my_id == master_id)
		planar_structure_fn_file << "%% Time = " << Tnow << endl; 
	
		
	CV_Compute_planar_structure_function();
	W.CV_Compute_planar_structure_function();
	
	if (my_id == master_id)
	{
		planar_structure_fn_file	<< "%% Planar struncture function:  V " << endl 
									<< (*CV_st_planar)		<< endl << endl;
						
		planar_structure_fn_file	<< "%% Planar struncture function: Vector field W " 
									<< endl 
									<< (*W.CV_st_planar)		<< endl << endl;	
	}														
}
  
//
//

void IncFluid::Output_planar_structure_fn(IncVF& W, IncSF& T)
{
/*
	if (my_id == master_id)
		planar_structure_fn_file << "%% Time = " << Tnow << endl; 
		

	CV_Compute_planar_structure_function();
	W.CV_Compute_planar_structure_function();
	T.CS_Compute_planar_structure_function();
		
	if (my_id == master_id)
	{			
		planar_structure_fn_file	<< "%% Planar struncture function:  V " << endl 
									<< (*CV_st_planar)		<< endl << endl;
						
		planar_structure_fn_file	<< "%% Planar struncture function: Vector field W " << endl 
									<< (*W.CV_st_planar)		<< endl << endl;	
						
		planar_structure_fn_file	<< "%% Planar struncture function: Scalar field T " << endl 
									<< (*T.CS_st_planar)		<< endl << endl;					
	}	
	*/													
}
//
//

void IncFluid::Output_planar_structure_fn(IncSF& T, string Pr_switch)
{
	if (Pr_switch == "PRZERO")
		Output_planar_structure_fn();
	
	else
		Output_planar_structure_fn(T);
}

void IncFluid::Output_planar_structure_fn(IncVF& W, IncSF& T, string Pr_switch)
{
	if (Pr_switch == "PRZERO")
		Output_planar_structure_fn(W);
	
	else
		Output_planar_structure_fn(W, T);
}


//********************************* End of output_str_fn **************************************


