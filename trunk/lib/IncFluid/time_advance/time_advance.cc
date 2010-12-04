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


/*! \file  time_advance.cc
 * 
 * @brief Time advances velocity field in Incompressible NS by one unit
 *
 * Status before entering the fn   <BR>
 *		nlin = FT[U.grad U]; U.F = p(k);  <BR>
 *		(computed by Compute_nlin and compute_pressure)  <BR>
 *	
 *	Procedure (EULER)	 <BR>											
 *		1. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()     <BR>          
 *		2. V(t+dt) = V(t) + dt*nlin in fn single_time_step() <BR>
 *		
 *	Procedure (RK2)	 <BR>
 *		1. Save V[i] in Vcopy[i] <BR>
 *		2. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()  <BR>
 *		3. V(t+dt) = V(t) + (dt/2)*nlin using fn single_time_step(): MID POINT <BR>
 *		4. compute_nlin(V(t+dt/2)) <BR>
 *		5. compute p(t+dt/2) using compute_pressure() <BR>
 *		6. V[i] = Vcopy[i] (copy back) <BR>
 *		7. V(t+dt) = V(t) + dt*nlin(t+dt/2) using fn single_time_step()  <BR>
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug No known bugs
 */                
   
#include "../IncFluid.h"

//*********************************************************************************************
 
void IncFluid::Time_advance()
{

//*********************************************************************************************
//	EULER	
//

	if (integ_scheme == "EULER") 
	{
		Compute_rhs();  
		Single_time_step(Tdt, 1, 1, 1);   
	}
	
	
//*********************************************************************************************		
//	RK2		
//
	
	else if (integ_scheme == "RK2") 
	{
		// Allocated once and for all- Vcopy location (static) 
		static PlainCVF Vcopy(N);
							
		Copy_field_to(Vcopy);
  
		Compute_rhs();  
		Single_time_step(Tdt, 0.5, 0.5, 0.5);					// Goto the mid point
  
		Compute_force_TO_rhs();
  
		Copy_field_from(Vcopy);
		
		Single_time_step(Tdt, 1, 0.5, 1);						
		// Time-step by Tdt now using the mid-point slopes 
	}
	
//*********************************************************************************************
//		RK4
//
	
	else if (integ_scheme == "RK4") 
	{
		// Allocated once and for all- Vcopy location (static)
		static PlainCVF Vcopy(N);									
		static PlainCVF tot_Vrhs(N);
											
		Copy_field_to(Vcopy);
		tot_Vrhs.PlainCV_Initialize();
		
		
		// Step 1
		
		Compute_rhs();    
		Single_time_step(Tdt, 0.5, 0.5, 0.5);					// u1: Goto the mid point

		// tot_Vrhs += (Tdt/6)*RHS(t)*exp(-nu k^2 b*dt) with b=1.0
		//	RHS is contained in *nlin
		Compute_RK4_parts(tot_Vrhs, Tdt, 1.0, Tdt/6);
		
		// Step 2
		
		Compute_force_TO_rhs();
		
		Copy_field_from(Vcopy);
		
		Single_time_step(Tdt, 0.5, 0, 0.5);						// u2: Goto the mid point
		
		Compute_RK4_parts(tot_Vrhs, Tdt, 0.5, Tdt/3);
		
		// Step 3
		
		Compute_force_TO_rhs();
		
		Copy_field_from(Vcopy);
		
		Single_time_step(Tdt, 1, 0.5, 1);						// u3: Go to the end point
		
		Compute_RK4_parts(tot_Vrhs, Tdt, 0.5, Tdt/3);
		
		// Step 4
		
		Compute_force_TO_rhs();
										
		Compute_RK4_parts(tot_Vrhs, Tdt, 0, Tdt/6);				// tot_rhs += rhs(t+dt, u3)
		
		// Final result
		
		
		Copy_field_from(Vcopy);
		Mult_field_exp_ksqr_dt(Tdt, 1.0);
			
		*V1 = *V1 + (*tot_Vrhs.V1);
		*V2 = *V2 + (*tot_Vrhs.V2);
		*V3 = *V3 + (*tot_Vrhs.V3);
		 
	}		// Of RK4
	
	if (alias_switch == "DEALIAS")		Dealias();					// Keep V(k) dealiased	
}

//**********************************   End of time_advance.cc  ********************************
   

