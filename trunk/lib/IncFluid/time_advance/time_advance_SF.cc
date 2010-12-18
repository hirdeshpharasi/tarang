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
		Single_time_step(T, Tdt, 1, 1, 1);     
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
		  
		Single_time_step(T, Tdt, 0.5, 0.5, 0.5);			// Goto the mid point
		
		Compute_force_TO_rhs(T);
  

		Copy_field_from(Vcopy); 
		T.Copy_field_from(Scopy);							// Copy back into V[i],T

		Single_time_step(T, Tdt, 1, 0.5, 1);							
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
		Single_time_step(T, Tdt, 0.5, 0.5, 0.5);				// u1: Goto the mid point

		// tot_Vrhs += (Tdt/6)*RHS(t)*exp(-nu k^2 b*dt) with b=1.0
		//	RHS is contained in *nlin
		Compute_RK4_parts(T, tot_Vrhs, tot_Srhs, Tdt, 1.0, Tdt/6); 
		
		// Step 2
		
		Compute_force_TO_rhs(T);
		
		Copy_field_from(Vcopy);	
		T.Copy_field_from(Scopy);
		
		Single_time_step(T, Tdt, 0.5, 0, 0.5);					// u2: Goto the mid point
		
		Compute_RK4_parts(T, tot_Vrhs, tot_Srhs, Tdt, 0.5, Tdt/3);
		
		// Step 3
		
		Compute_force_TO_rhs(T);
		
		Copy_field_from(Vcopy);  
		T.Copy_field_from(Scopy);
		
		Single_time_step(T, Tdt, 1, 0.5, 1);					// u3: Go to the end point
		
		Compute_RK4_parts(T, tot_Vrhs, tot_Srhs, Tdt, 0.5, Tdt/3);
		
		// Step 4
		
		Compute_force_TO_rhs(T);
										
		Compute_RK4_parts(T, tot_Vrhs, tot_Srhs, Tdt, 0, Tdt/6);
		
		// Final result
		
		Copy_field_from(Vcopy);  
		T.Copy_field_from(Scopy);
		
		Mult_field_exp_ksqr_dt(Tdt, 1.0);  
		T.Mult_field_exp_ksqr_dt(Tdt, 1.0);
		
		if  ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		{
			*V1 = *V1 + (*tot_Vrhs.V1);
			*V2 = *V2 + (*tot_Vrhs.V2);
			*V3 = *V3 + (*tot_Vrhs.V3);
			
			*T.F = *T.F + (*tot_Srhs.F);
		}
		
		else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		{
			
			if (globalvar_Pr_switch == "PRZERO") 
			{
					
				*V1 = *V1 + (*tot_Vrhs.V1);
				*V2 = *V2 + (*tot_Vrhs.V2);
				*V3 = *V3 + (*tot_Vrhs.V3);
				
				*T.F = *V1;
				Array_divide_ksqr_SCFT(N, *T.F, kfactor);
			}
			
			else if (globalvar_Pr_switch == "PRINFTY") 
			{
				
				*T.F = *T.F + (*tot_Srhs.F);
				
				*nlin1 = 0.0;
				*nlin2 = 0.0;
				*nlin3 = 0.0;
				
				Compute_force(T);
				Add_force(T);
				Compute_pressure();									// Compute pressure using V(t+dt/2)
				Compute_rhs();
				
				*V1 = *nlin1;
				*V2 = *nlin2;
				*V3 = *nlin3;
				
				Array_divide_ksqr_SCFT(N, *V1, kfactor);				// u(k) = nlin(k)/k^2	
				Array_divide_ksqr_SCFT(N, *V2, kfactor);
				Array_divide_ksqr_SCFT(N, *V3, kfactor);
			}
			
			else
			{	
				*V1 = *V1 + (*tot_Vrhs.V1);
				*V2 = *V2 + (*tot_Vrhs.V2);
				*V3 = *V3 + (*tot_Vrhs.V3);
				
				*T.F = *T.F + (*tot_Srhs.F);
			}
		}
		 
	}		// of RK4
	
	// if free_slip_verticalwall condition is initialized as IC, it should
	// be satisfied at all times, yet we set it again just to make sure everytime step.
	if ((free_slip_verticalwall_switch == 1) && (basis_type == "SCFT"))
		free_slip_verticalwall(T);
	
	//
	// For all schemes
	if (alias_switch == "DEALIAS")		Dealias(T);					// Keep V(k) dealiased	

}

//**********************************   End of time_advance_SF.cc  *****************************


