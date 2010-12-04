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

/*! \file  mult_exp_visc.cc
 * 
 * @brief Multiply fields with \f$ \exp(-\nu K^2 dt) \f$ or like.
 *   When hyper_dissipation_switch == 1 (on), then multiply by 
 *	\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 *
 * @sa void IncVF::Mult_field_exp_ksqr_dt(DP dt)
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************

/** @brief Multiply vector field by \f$ \exp(-\nu K^2 dt) \f$ or 
 *			\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 * 
 *  @param dt
 *
 * @return \f$ V(k) = V(k) * \exp(-\nu K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ V(k) = V(k) * \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 */
void IncVF::Mult_field_exp_ksqr_dt(DP dt, DP a)
{
	
	if (hyper_dissipation_switch == 0) 
	{
		Array_exp_ksqr(basis_type, N,  *V1, -dissipation_coefficient*a*dt, kfactor);
		Array_exp_ksqr(basis_type, N,  *V2, -dissipation_coefficient*a*dt, kfactor);
		Array_exp_ksqr(basis_type, N,  *V3, -dissipation_coefficient*a*dt, kfactor);
	}
	
	else
	{
		Array_exp_ksqr(basis_type, N,  *V1, -dissipation_coefficient*a*dt, 
						-hyper_dissipation_coefficient*a*dt, kfactor);
						
		Array_exp_ksqr(basis_type, N,  *V2, -dissipation_coefficient*a*dt, 
						-hyper_dissipation_coefficient*a*dt, kfactor);
						
		Array_exp_ksqr(basis_type, N,  *V3, -dissipation_coefficient*a*dt, 
						-hyper_dissipation_coefficient*a*dt, kfactor);
	}

}



//*********************************************************************************************
/** @brief Multiply nonlinear field by \f$ \exp(-\nu K^2 dt) \f$ or 
 *			\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 * 
 *  @param dt
 *
 * @return \f$ N(k) = N(k) * \exp(-\nu K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ V(k) = V(k) * \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 */

void IncVF::Mult_nlin_exp_ksqr_dt(DP dt, DP a)
{

	if (hyper_dissipation_switch == 0) 
	{
		Array_exp_ksqr(basis_type, N,  *nlin1, -dissipation_coefficient*a*dt, kfactor);
		Array_exp_ksqr(basis_type, N,  *nlin2, -dissipation_coefficient*a*dt, kfactor);
		Array_exp_ksqr(basis_type, N,  *nlin3, -dissipation_coefficient*a*dt, kfactor);
	}
	
	else
	{
		Array_exp_ksqr(basis_type, N,  *nlin1, -dissipation_coefficient*a*dt, 
						-hyper_dissipation_coefficient*a*dt, kfactor);
						
		Array_exp_ksqr(basis_type, N,  *nlin2, -dissipation_coefficient*a*dt, 
						-hyper_dissipation_coefficient*a*dt, kfactor);
						
		Array_exp_ksqr(basis_type, N,  *nlin3, -dissipation_coefficient*a*dt, 
						-hyper_dissipation_coefficient*a*dt, kfactor);
	}

}



//*********************************************************************************************
/** @brief Multiply scalar field by \f$ \exp(-\kappa K^2 dt) \f$ or 
 *			\f$ \exp(-\kappa K^2 dt) + \exp(-\kappa_h)  K^4 dt) \f$
 * 
 *  @param dt
 *
 * @return \f$ F(k) = F(k) * \exp(-\kappa K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ F(k) = F(k) * \exp(-\kappa K^2 dt) + \exp(-\kappa_h)  K^4 dt) \f$
 */
void IncSF::Mult_field_exp_ksqr_dt(DP dt, DP a)
{

	 if (hyper_diffusion_switch == 0) 
		Array_exp_ksqr(ISF_basis_type, NIs,  *F, -diffusion_coefficient*a*dt, ISF_kfactor); 
		
	 else
		Array_exp_ksqr(ISF_basis_type, NIs,  *F, -diffusion_coefficient*a*dt, 
						-hyper_diffusion_coefficient*a*dt, ISF_kfactor); 	
	
} 
 


//*********************************************************************************************
/** @brief Multiply nlin of scalar field by \f$ \exp(-\kappa K^2 dt) \f$ 
 *			\f$ \exp(-\kappa K^2 dt) + \exp(-\kappa_h)  K^4 dt) \f$
 * 
 *  @param dt
 *
 * @return \f$ N(k) = N(k) * \exp(-\kappa K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ F(k) = F(k) * \exp(-\kappa K^2 dt) + \exp(-\kappa_h)  K^4 dt) \f$
 */
void IncSF::Mult_nlin_exp_ksqr_dt(DP dt, DP a)
{ 

	if (hyper_diffusion_switch == 0) 						  
		Array_exp_ksqr(ISF_basis_type, NIs,  *nlin, -diffusion_coefficient*a*dt, ISF_kfactor); 
		
	else
		Array_exp_ksqr(ISF_basis_type, NIs,  *nlin, -diffusion_coefficient*a*dt, 
						-hyper_diffusion_coefficient*a*dt, ISF_kfactor);
} 

//******************************* End of mult_exp_visc.cc  ************************************


	
