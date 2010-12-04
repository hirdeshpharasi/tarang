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
/*! \file  time_advance_SF.cc
 * 
 * @brief Time advances velocity and passive scalar fields in Incompressible flow by one unit
 *
 * Status before entering the fn   <BR>
 *		nlin = FT[U.grad U]; U.F = p(k);  <BR>
 *		T.nlin = FT[U.grad T]
 *		(computed by Compute_nlin and compute_pressure)  <BR>
 *	
 *	Procedure (EULER)	 <BR>											
 *		1. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()     <BR>          
 *		2. V,T(t+dt) = V,T(t) + dt*nlin(V,T) in fn single_time_step() <BR>
 *		
 *	Procedure (RK2)	 <BR>
 *		1. Save V[i] in Vcopy[i] <BR>
 *		2. Compute_rhs(T): nlin[i] = -nlin[i] - grad(pressure) and T.nlin = -T.nlin <BR>
 *		3. V,T(t+dt) = V,T(t) + (dt/2)*nlin(V,T) using fn single_time_step(): MID POINT <BR>
 *		4. Compute_nlin(V(t+dt/2),T(t+dt/2)) <BR>
 *		5. compute p(t+dt/2) using compute_pressure() <BR>
 *		6. V,T[i] = V,Tcopy[i] (copy back) <BR>
 *		7. V,T(t+dt) = V,T(t) + dt*nlinU,T(t+dt/2) using fn single_time_step()  <BR>
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
                
   
#include "../IncFluid.h"

//********************************************************************************************* 


void IncFluid::Time_advance(IncSF& T)
{

//*********************************************************************************************
//		EULER

	if (integ_scheme == "EULER") 
	{ 
		Compute_rhs(T);              
		Single_time_step_EULER(T, Tdt);   
	}
		
		
//*********************************************************************************************		
//		RK2			
	
	else if (integ_scheme == "RK2") 
	{
		
		// Allocated once and for all- Vcopy,Scopy  (static)
		static PlainCVF Vcopy(N);							
		static PlainCSF Scopy(N); 

		
		Copy_field_to(Vcopy); 
		T.Copy_field_to(Scopy);								// Vcopy[i] <- V[i] ; Scopy <- T.F     
  
		Compute_rhs(T);  
		  
		Single_time_step_EULER(T,Tdt/2);					// Goto the mid point
		
		Compute_force(T);
  
		Compute_nlin(T);									// Compute nlin using V(t+dt/2)
		
		Satisfy_reality_condition_nlin(T);
		
		Add_force(T);										// nlin = nlin - f
		
		Compute_pressure();									// Compute pressure using V(t+dt/2)
		
		Compute_rhs(T);										// compute rhs using V(t+dt/2) 
  

		Copy_field_from(Vcopy); 
		T.Copy_field_from(Scopy);							// Copy back into V[i],T

		IncFluid::Single_time_step_RK2(T,Tdt);							
		// Time-step by Tdt now using the mid-point slopes 
	}
	
//*********************************************************************************************		
//		RK4		
	
	else if (integ_scheme == "RK4") 
	{
		
		static PlainCVF Vcopy(N);													
		static PlainCSF Scopy(N);
							  
		static PlainCVF tot_Vrhs(N);								
		static PlainCSF tot_Srhs(N);
								
		
		Copy_field_to(Vcopy);	
		T.Copy_field_to(Scopy);
		
		tot_Vrhs.PlainCV_Initialize();  
		tot_Srhs.PlainCS_Initialize();
		
		// Step 1
		
		Compute_rhs(T);    
		Single_time_step_EULER(T, Tdt/2);						// u1: Goto the mid point

		Mult_nlin_exp_ksqr_dt(Tdt);								// *nlin = *nlin x exp(-nu k^2 dt/d)
		Add_nlin_to_field(tot_Vrhs, 1.0);						// rhs(t0, u0) in nlin
		T.Mult_nlin_exp_ksqr_dt(Tdt);	
		T.Add_nlin_to_field(tot_Srhs, 1.0);
		
		// Step 2
		
		Compute_force(T);
		Compute_nlin(T);										// Compute nlin using V(t+dt/2)
		Satisfy_reality_condition_nlin(T);
		Add_force(T);											// nlin = nlin - f
		Compute_pressure();										// Compute pressure using V(t+dt/2)
		Compute_rhs(T);
		Copy_field_from(Vcopy);	T.Copy_field_from(Scopy);
		Single_time_step_Semi_implicit(T, Tdt/2);				// u2: Goto the mid point
		
		Mult_nlin_exp_ksqr_dt(Tdt/2);							// rhs(t_mid, u1) in nlin	
		Add_nlin_to_field(tot_Vrhs, 2.0);
		T.Mult_nlin_exp_ksqr_dt(Tdt/2);	
		T.Add_nlin_to_field(tot_Srhs, 2.0);
		
		// Step 3
		
		Compute_force(T);
		Compute_nlin(T);	
		Satisfy_reality_condition_nlin(T);												
		Add_force(T);													
		Compute_pressure();												
		Compute_rhs(T);
		Copy_field_from(Vcopy);  T.Copy_field_from(Scopy);
		Single_time_step_RK2(T, Tdt);							// u3: Go to the end point
		
		Mult_nlin_exp_ksqr_dt(Tdt/2);									
		Add_nlin_to_field(tot_Vrhs, 2.0);						// rhs(t_mid, u2) in nlin
		T.Mult_nlin_exp_ksqr_dt(Tdt/2);	
		T.Add_nlin_to_field(tot_Srhs, 2.0);
		
		// Step 4
		
		Compute_force(T);
		Compute_nlin(T);
		Satisfy_reality_condition_nlin(T);													
		Add_force(T);													
		Compute_pressure();												
		Compute_rhs(T);
										
		Add_nlin_to_field(tot_Vrhs, 1.0);						// rhs(t_mid, u3) in nlin
		T.Add_nlin_to_field(tot_Srhs, 1.0);
		
		// Final result
		
		Copy_field_from(Vcopy);  T.Copy_field_from(Scopy);
		Mult_field_exp_ksqr_dt(Tdt);  T.Mult_field_exp_ksqr_dt(Tdt);
			
		*V1 = *V1 + (Tdt/6) * (*tot_Vrhs.V1);
		*V2 = *V2 + (Tdt/6) * (*tot_Vrhs.V2);
		*V3 = *V3 + (Tdt/6) * (*tot_Vrhs.V3);
		
		*T.F = *T.F + (Tdt/6) * (*tot_Srhs.F);
		 
	}		// of RK4
	
	//
	// For all schemes
	if (alias_switch == "DEALIAS")		Dealias(T);					// Keep V(k) dealiased	

}

//**********************************   End of time_advance_SF.cc  *****************************


