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


/*! \file  single_time_step.cc
 * 
 * @brief single time step fields using Euler, RK2, and semi-implicit scheme.
 * @sa void IncFluid::Single_time_step_EULER(DP dt)
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "../IncFluid.h"


//*********************************************************************************************

/** @brief Single time step using Euler's scheme.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) + dt \vec{N}(t)] 
 *					*\exp(-K^2*\nu*dt)
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$
 */
void IncFluid::Single_time_step_EULER(DP dt)
{	

	Add_nlin_dt(dt);										// V = V + dt*nlin
	Mult_field_exp_ksqr_dt(dt);
}

	
//*********************************************************************************************	
//	Passive scalar 
/** @brief Passive scalar Single time step using Euler's scheme.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) + dt \vec{N}(t)] 
 *					*\exp(-K^2*\nu*dt).
 * @return \f$ F(t+dt) = [F(t) + dt T.N(t)] 
 *					*\exp(-K^2*\kappa*dt).
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$. Similarly for scalar
 */
void IncFluid::Single_time_step_EULER(IncSF& T, DP dt)
{						
	Single_time_step_EULER(dt);	
	
	T.Add_nlin_dt(dt);
	T.Mult_field_exp_ksqr_dt(dt);	
}  


//*********************************************************************************************
//	MHD		
/** @brief Passive scalar Single time step using Euler's scheme.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) + dt \vec{N}(t)] 
 *					*\exp(-K^2*\nu*dt).
 * @return \f$ \vec{W}(t+dt) = [\vec{W}(t) + dt W.\vec{N}(t)] 
 *					*\exp(-K^2*W.\nu*dt).
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$. Similarly for vector W.
 */
void IncFluid::Single_time_step_EULER(IncVF& W, DP dt)
{  
	Single_time_step_EULER(dt);									

	W.Add_nlin_dt(dt);
	W.Mult_field_exp_ksqr_dt(dt);
}  


//*********************************************************************************************
//	MHD + passive scalar  
/** @brief Passive scalar Single time step using Euler's scheme.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) + dt \vec{N}(t)] 
 *					*\exp(-K^2*\nu*dt).
 * @return \f$ \vec{W}(t+dt) = [\vec{W}(t) + dt W.\vec{N}(t)] 
 *					*\exp(-K^2*W.\nu*dt).
 * @return \f$ F(t+dt) = [F(t) + dt T.N(t)] 
 *					*\exp(-K^2*\kappa*dt).
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$. Similarly for vector W &scalar
 */
void IncFluid::Single_time_step_EULER(IncVF& W, IncSF& T,  DP dt)
{
	Single_time_step_EULER(W, dt);
																										
	T.Add_nlin_dt(dt);
	T.Mult_field_exp_ksqr_dt(dt);
}



//*********************************************************************************************
//	RB Convection	
/** @brief RB convection: Single time step using Euler's scheme.
 *
 * @return If Pr != 0, time-advance passive scalar.
 *			Else time-advance velocity field, and \f$ F = V_x /K^2 \f$.
 */
void IncFluid::Single_time_step_EULER(IncSF& T, string Pr_switch, DP dt)
{

	if (Pr_switch == "PRZERO") 
	{
		Single_time_step_EULER(dt);								// V advance
		*T.F = *V1;
		Array_divide_ksqr_SCFT(N, *T.F, kfactor);				// theta(k) = v1(k)/k^2
	}  

	else
		Single_time_step_EULER(T, dt);

}


//*********************************************************************************************
//	Magnetoconvection	
/** @brief RB magnetoconvection: Single time step using Euler's scheme.
 *
 * @return If Pr != 0, time-advance passive scalar & vector W.
 *			Else time-advance velocity & W fields, and \f$ F = V_x /K^2 \f$.
 */
void IncFluid::Single_time_step_EULER(IncVF& W, IncSF& T, string Pr_switch, DP dt)
{

	if (Pr_switch == "PRZERO") 
	{
		Single_time_step_EULER(W, dt);
		*T.F = *V1;
		Array_divide_ksqr_SCFT(N, *T.F, kfactor);							// theta(k) = v1(k)/k^2
	}  

	else 
		Single_time_step_EULER(W, T, dt);
}



//*********************************************************************************************
//*********************************************************************************************

/** @brief Single time step using RK2 scheme.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) *\exp(-K^2*\nu*dt) + dt \vec{N}(t)] 
 *					*\exp(-K^2*\nu*dt)
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$
 */
void IncFluid::Single_time_step_RK2(DP dt)
{	
	Mult_field_exp_ksqr_dt(dt / 2);
	Add_nlin_dt(dt);	
	Mult_field_exp_ksqr_dt(dt / 2);
}

	
//*********************************************************************************************	
//	Passive scalar  
/** @brief Single time step using RK2 scheme.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) *\exp(-K^2*\nu*dt) + dt \vec{N}(t)] 
 *					*\exp(-K^2*\nu*dt)
 * @return \f$ F(t+dt) = [F(t) *\exp(-K^2*\kappa*dt)  + dt T.N(t)] 
 *					*\exp(-K^2*\kappa*dt).
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$. Similarly for scalar.
 */
void IncFluid::Single_time_step_RK2(IncSF& T, DP dt)
{
	Single_time_step_RK2(dt);								// V time-advance
 
	T.Mult_field_exp_ksqr_dt(dt / 2);
	T.Add_nlin_dt(dt);
	T.Mult_field_exp_ksqr_dt(dt / 2);	
}  


//*********************************************************************************************
//	MHD		
/** @brief Single time step using RK2 scheme for V and W.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) *\exp(-K^2*\nu*dt) + dt \vec{N}(t)] 
 *					*\exp(-K^2*\nu*dt)
 * @return \f$ \vec{W}(t+dt) = [\vec{W}(t) *\exp(-K^2*W.\nu*dt) + dt W.\vec{N}(t)] 
 *					*\exp(-K^2*W.\nu*dt)
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$.  Similar for W.
 */
void IncFluid::Single_time_step_RK2(IncVF& W, DP dt)
{  
	Single_time_step_RK2(dt);									// V advance
    
	W.Mult_field_exp_ksqr_dt(dt / 2);
	W.Add_nlin_dt(dt);
	W.Mult_field_exp_ksqr_dt(dt / 2);
}  


//*********************************************************************************************
//	MHD + passive scalar  
/** @brief Single time step using RK2 scheme for V, W, and T.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) *\exp(-K^2*\nu*dt) + dt \vec{N}(t)] 
 *					*\exp(-K^2*\nu*dt)
 * @return \f$ \vec{W}(t+dt) = [\vec{W}(t) *\exp(-K^2*W.\nu*dt) + dt W.\vec{N}(t)] 
 *					*\exp(-K^2*W.\nu*dt)
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$.  Similar for W.
 */
void IncFluid::Single_time_step_RK2(IncVF& W, IncSF& T,  DP dt)
{
	Single_time_step_RK2(W, dt);
	
	T.Mult_field_exp_ksqr_dt(dt / 2);																										
	T.Add_nlin_dt(dt);
	T.Mult_field_exp_ksqr_dt(dt / 2);
}



//*********************************************************************************************
//	RB Convection	
/** @brief RB convection: Single time step using RK2 scheme.
 *
 * @return If Pr != 0, time-advance passive scalar.
 *			Else time-advance velocity field, and \f$ F = V_x /K^2 \f$.
 */
void IncFluid::Single_time_step_RK2(IncSF& T, string Pr_switch, DP dt)
{
	if (Pr_switch == "PRZERO") 
	{
		Single_time_step_RK2(dt);								// V advance
		*T.F = *V1;
		Array_divide_ksqr_SCFT(N, *T.F, kfactor);				// theta(k) = v1(k)/k^2
	}  

	else
		Single_time_step_RK2(T, dt);
}


//*********************************************************************************************
//	Magnetoconvection	
/** @brief RB magnetoconvection: Single time step using Euler's scheme.
 *
 * @return If Pr != 0, time-advance passive scalar & vector W.
 *			Else time-advance velocity & W fields, and \f$ F = V_x /K^2 \f$.
 */
void IncFluid::Single_time_step_RK2(IncVF& W, IncSF& T, string Pr_switch, DP dt)
{

	if (Pr_switch == "PRZERO") 
	{
		Single_time_step_RK2(W, dt);
		*T.F = *V1;
		Array_divide_ksqr_SCFT(N, *T.F, kfactor);				// theta(k) = v1(k)/k^2
	}  

	else 
		Single_time_step_RK2(W, T, dt);
}


//*********************************************************************************************
//*********************************************************************************************

/** @brief Single time step using Semi_implicit scheme.
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) *\exp(-K^2*\nu*dt) + dt \vec{N}(t+dt)] 
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$
 */
void IncFluid::Single_time_step_Semi_implicit(DP dt)
{	
	Mult_field_exp_ksqr_dt(dt);	
	Add_nlin_dt(dt);										// V = V + dt*nlin
}

	
//*********************************************************************************************	
//	Passive scalar 
void IncFluid::Single_time_step_Semi_implicit(IncSF& T, DP dt)
{						
	Single_time_step_Semi_implicit(dt);	
	
	T.Mult_field_exp_ksqr_dt(dt);	
	T.Add_nlin_dt(dt);
}  


//*********************************************************************************************
//	MHD	
void IncFluid::Single_time_step_Semi_implicit(IncVF& W, DP dt)
{  
	Single_time_step_Semi_implicit(dt);									

	W.Mult_field_exp_ksqr_dt(dt);
	W.Add_nlin_dt(dt);
}  


//*********************************************************************************************
//	MHD + passive scalar 
void IncFluid::Single_time_step_Semi_implicit(IncVF& W, IncSF& T,  DP dt)
{
	Single_time_step_Semi_implicit(W, dt);
																										
	T.Mult_field_exp_ksqr_dt(dt);
	T.Add_nlin_dt(dt);
}



//*********************************************************************************************
//	RB Convection	

void IncFluid::Single_time_step_Semi_implicit(IncSF& T, string Pr_switch, DP dt)
{

	if (Pr_switch == "PRZERO")
	{
		Single_time_step_Semi_implicit(dt);						// V advance
		*T.F = *V1;
		Array_divide_ksqr_SCFT(N, *T.F, kfactor);				// theta(k) = v1(k)/k^2
	}  

	else
		Single_time_step_Semi_implicit(T, dt);
}


//*********************************************************************************************
//	Magnetoconvection

void IncFluid::Single_time_step_Semi_implicit(IncVF& W, IncSF& T, string Pr_switch, DP dt)
{

	if (Pr_switch == "PRZERO") 
	{
		Single_time_step_Semi_implicit(W, dt);
		*T.F = *V1;
		Array_divide_ksqr_SCFT(N, *T.F, kfactor);				// theta(k) = v1(k)/k^2
	}  

	else 
		Single_time_step_Semi_implicit(W, T, dt);
}


//**********************************   End of Single_time_step.cc  ****************************



