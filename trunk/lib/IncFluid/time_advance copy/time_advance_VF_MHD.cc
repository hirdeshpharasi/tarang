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

/*! \file  time_advance_VF_MHD.cc
 * 
 * @brief Time advances velocity and magnetic fields in MHD flow by one unit
 *
 * Status before entering the fn   <BR>
 *		nlin = FT[U.grad U-B.grad B]; U.F = p(k);  <BR>
 *		W.nlin = FT[U.grad B - B.grad U]
 *		(computed by Compute_nlin and compute_pressure)  <BR>
 *	
 *	Procedure (EULER)	 <BR>											
 *		1. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()     <BR>          
 *		2. V,W(t+dt) = V,W(t) + dt*nlin(V,W) in fn single_time_step() <BR>
 *		
 *	Procedure (RK2)	 <BR>
 *		1. Save V[i] in Vcopy[i] <BR>
 *		2. Compute_rhs(T): nlin[i] = -nlin[i] - grad(pressure) and T.nlin = -T.nlin <BR>
 *		3. V,T(t+dt) = V,T(t) + (dt/2)*nlin(V,T) using fn single_time_step(): MID POINT <BR>
 *		4. Compute_nlin(V(t+dt/2),T(t+dt/2)) <BR>
 *		5. compute p(t+dt/2) using compute_pressure() <BR>
 *		6. V,T[i] = V,Tcopy[i] (copy back) <BR>
 *		7. V,T(t+dt) = V,T(t) + dt*nlinU,T(t+dt/2) using fn single_time_step()  <BR> *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 */
          
   
#include "../IncFluid.h"

//*********************************************************************************************


void IncFluid::Time_advance(IncVF& W)
{

//*********************************************************************************************
//		EULER


	if (integ_scheme == "EULER") 
	{
		Compute_rhs(W);              
		Single_time_step_EULER(W, Tdt);						// Single-time-step both U and W 
	}


//*********************************************************************************************
//		RK2

	
	else if (integ_scheme == "RK2") 
	{
		
		// Allocated once and for all- Vcopy, Wcopy (static)
		static PlainCVF Vcopy(N);		
		static PlainCVF Wcopy(N);     
		
		Copy_field_to(Vcopy); 
		W.Copy_field_to(Wcopy);
     
		Compute_rhs(W);
		Single_time_step_EULER(W,Tdt/2);					// Goto the mid point
  
		Compute_force(W);
  
		Compute_nlin(W);									// Compute nlin(Zp,Zm) at mid point
		
		Satisfy_reality_condition_nlin(W);
			
		Add_force(W);										// nlin = nlin - f (for V and B)
    
		Compute_pressure();									// Compute pressure using V(t+dt/2)
		Compute_rhs(W);	
		
		Copy_field_from(Vcopy); 
		W.Copy_field_from(Wcopy);							 							 
  
		Single_time_step_RK2(W, Tdt);							
		// Time-step by Tdt now using the mid-point slopes
	}
	
	
//*********************************************************************************************
//		RK4


	
	else if (integ_scheme == "RK4") 
	{
		
		static PlainCVF Vcopy(N);										
		static PlainCVF Wcopy(N);
							
		static PlainCVF tot_Vrhs(N);								
		static PlainCVF tot_Wrhs(N);
							
	
		Copy_field_to(Vcopy); 
		W.Copy_field_to(Wcopy);
		
		tot_Vrhs.PlainCV_Initialize(); 
		tot_Wrhs.PlainCV_Initialize();
		
		Compute_rhs(W);    
		Single_time_step_EULER(W, Tdt/2);						// u1: Goto the mid point

		Mult_nlin_exp_ksqr_dt(Tdt);								// *nlin = *nlin x exp(-nu k^2 dt/d)
		Add_nlin_to_field(tot_Vrhs, 1.0);						// rhs(t0, u0) in nlin
		W.Mult_nlin_exp_ksqr_dt(Tdt);	
		W.Add_nlin_to_field(tot_Wrhs, 1.0);
		
		// Step 2
		
		Compute_force(W);
		Compute_nlin(W);										// Compute nlin using V(t+dt/2)
		Satisfy_reality_condition_nlin(W);
		Add_force(W);											// nlin = nlin - f
		Compute_pressure();										// Compute pressure using V(t+dt/2)
		Compute_rhs(W);
		Copy_field_from(Vcopy);	W.Copy_field_from(Wcopy);
		Single_time_step_Semi_implicit(W, Tdt/2);				// u2: Goto the mid pt 
		
		Mult_nlin_exp_ksqr_dt(Tdt/2);	
		W.Mult_nlin_exp_ksqr_dt(Tdt/2);							
		Add_nlin_to_field(tot_Vrhs, 2.0);						// rhs(t_mid, u1) in nlin		
		W.Add_nlin_to_field(tot_Wrhs, 2.0);	
		
		// Step 3
		
		Compute_force(W);
		Compute_nlin(W);										// Compute nlin using u2
		Satisfy_reality_condition_nlin(W);
		Add_force(W);											
		Compute_pressure();										
		Compute_rhs(W);
		Copy_field_from(Vcopy);	W.Copy_field_from(Wcopy);
		Single_time_step_RK2(W, Tdt);							// u3 : Goto the mid pt 
		
		Mult_nlin_exp_ksqr_dt(Tdt/2);	
		W.Mult_nlin_exp_ksqr_dt(Tdt/2);							
		Add_nlin_to_field(tot_Vrhs, 2.0);						// rhs(t_mid, u2) in nlin
		W.Add_nlin_to_field(tot_Wrhs, 2.0);	
		
		// Step 4
		
		Compute_force(W);
		Compute_nlin(W);										// Compute nlin using u3
		Satisfy_reality_condition_nlin(W);
		Add_force(W);											
		Compute_pressure();										
		Compute_rhs(W);
						
		Add_nlin_to_field(tot_Vrhs, 1.0);						// rhs(t_mid, u3) in nlin
		W.Add_nlin_to_field(tot_Wrhs, 1.0);
		
		// Final result
		
		Copy_field_from(Vcopy);  W.Copy_field_from(Wcopy);
		Mult_field_exp_ksqr_dt(Tdt); W.Mult_field_exp_ksqr_dt(Tdt);
				
		 *V1 = *V1 + (Tdt/6) * (*tot_Vrhs.V1);
		 *V2 = *V2 + (Tdt/6) * (*tot_Vrhs.V2);
		 *W.V1 = *W.V1 + (Tdt/6) * (*tot_Wrhs.V1);
		 *W.V2 = *W.V2 + (Tdt/6) * (*tot_Wrhs.V2);
		 
		 *V3 = *V3 + (Tdt/6) * (*tot_Vrhs.V3);
		 *W.V3 = *W.V3 + (Tdt/6) * (*tot_Wrhs.V3);	

	}   // Of RK4
	
	//
	//  For all schemes
	if (alias_switch == "DEALIAS")		Dealias(W);					// Keep V(k) dealiased	
	
}

//**********************************   End of time_advance_VF_MHD.cc  *************************



