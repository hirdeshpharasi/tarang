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

/*! \file  ic_model.cc
 * 
 * @brief Initial condition where low k modes are forced: TG, ABC flows etc.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug	 No known bugs
 */


#include "../IncVF.h"
#include "../IncSF.h"

extern Uniform<DP> SPECrand;



//*********************************************************************************************

/** @brief Init condition for Taylor Green flow 
 * 
 * @note Vx = amp*sin(k0 x) cos(k0 y) cos(k0 z)
 *		Vy = -amp*cos(k0 x) sin(k0 y) cos(k0 z)
 *		Vz = 0
 *
 *  @param k0
 *	@param amp  
 *
 */
void IncVF::Setup_Taylor_Green_field(int k0, DP amp)
{

	*V1 = 0.0; 
	*V2 = 0.0;  
	*V3 = 0.0;				// initialize

	if (basis_type == "FOUR")
	{
		int lx_k0 = Get_lx("FOUR", k0, N);
		int lx_minus_k0 = Get_lx("FOUR", -k0, N);
		int ly_k0 = Get_ly3D("FOUR", k0, N);
		int ly_minus_k0 = Get_ly3D("FOUR", -k0, N);
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V1)(lx_k0, ly_k0, k0)  = -I * (amp/8);
			(*V1)(lx_k0, ly_minus_k0, k0) = -I * (amp/8);
			
			(*V2)(lx_k0, ly_k0, k0)  = I * (amp/8);
			(*V2)(lx_k0, ly_minus_k0, k0) = -I * (amp/8);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{
			(*V1)(lx_minus_k0, ly_k0, k0) = I * (amp/8);
			(*V1)(lx_minus_k0, ly_minus_k0, k0)= I * (amp/8);
			
			(*V2)(lx_minus_k0, ly_k0, k0) = I * (amp/8);
			(*V2)(lx_minus_k0, ly_minus_k0, k0)= -I * (amp/8);
		}
	}
	
	else if (basis_type == "SCFT")
	{
		int lx_k0 = Get_lx("SCFT", k0, N);
		int ly_k0 = Get_ly3D("SCFT", k0, N);
		int ly_minus_k0 = Get_ly3D("SCFT", -k0, N);
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V1)(lx_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V1)(lx_k0, ly_minus_k0, k0) = complx(amp/8, 0.0);
			
			(*V2)(lx_k0, ly_k0, k0) = I *(amp/8);
			(*V2)(lx_k0, ly_minus_k0, k0) = -I *(amp/8);
		}
	}

}


//*********************************************************************************************

/** @brief Init condition for Taylor Green flow 
 * 
 * @note  Vx = amp*sin(k0 x) sin(k0 y) cos(k0 z)
 *		Vy = amp*cos(k0 x) cos(k0 y) cos(k0 z)
 *		Vz = 0
 *
 *  @param k0
 *	@param amp  
 *
 */
void IncVF::Setup_Taylor_Green_field_mag_field(int k0, DP amp)
{

	*V1 = 0.0; 
	*V2 = 0.0;  
	*V3 = 0.0;				// initialize

	if (basis_type == "FOUR")
	{
		int lx_k0 = Get_lx("FOUR", k0, N);
		int lx_minus_k0 = Get_lx("FOUR", -k0, N);
		int ly_k0 = Get_ly3D("FOUR", k0, N);
		int ly_minus_k0 = Get_ly3D("FOUR", -k0, N);
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V1)(lx_k0, ly_k0, k0)  = complx(-amp/8, 0.0);
			(*V1)(lx_k0, ly_minus_k0, k0) = complx(amp/8, 0.0);
			
			(*V2)(lx_k0, ly_k0, k0)  = complx(amp/8, 0.0);
			(*V2)(lx_k0, ly_minus_k0, k0) = complx(amp/8, 0.0);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{
			(*V1)(lx_minus_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V1)(lx_minus_k0, ly_minus_k0, k0)= complx(-amp/8, 0.0);
			
			(*V2)(lx_minus_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V2)(lx_minus_k0, ly_minus_k0, k0)= complx(amp/8, 0.0);
		}
	}
	
	else if (basis_type == "SCFT")
	{
		int lx_k0 = Get_lx("SCFT", k0, N);
		int ly_k0 = Get_ly3D("SCFT", k0, N);
		int ly_minus_k0 = Get_ly3D("SCFT", -k0, N);
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V1)(lx_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V1)(lx_k0, ly_minus_k0, k0) = complx(-amp/8, 0.0);
			
			(*V2)(lx_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V2)(lx_k0, ly_minus_k0, k0) = complx(amp/8, 0.0);
		}
	}

}



//*********************************************************************************************

/** @brief Init condition for Taylor Green flow 
 * 
 * @note \f$ V_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ V_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ V_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 *  @param k0
 *	@param amp  
 *
 */
void IncVF::Setup_ABC_field(int k0, DP amp, DP A, DP B, DP C)
{

	*V1 = 0.0; 
	*V2 = 0.0;  
	*V3 = 0.0;			// initialize

	if (basis_type == "FOUR")
	{
		int lx_k0 = Get_lx("FOUR", k0, N);
		int lx_minus_k0 = Get_lx("FOUR", -k0, N);
		int ly_k0 = Get_ly3D("FOUR", k0, N);
		int ly_minus_k0 = Get_ly3D("FOUR", -k0, N);
		
		if (my_id == master_id)  
		{
			(*V1)(0, ly_k0, 0)  =	complx(B * amp/2, 0.0);
			(*V1)(0, ly_minus_k0, 0) =  complx(B * amp/2, 0.0);
			(*V1)(0, 0, k0)  =  -I * (C * amp/2);
			
			(*V2)(0, 0, k0)  =   complx(C * amp/2, 0.0);
			
			(*V3)(0, ly_k0, 0)  =   -I * (B * amp/2);
			(*V3)(0, ly_minus_k0, 0) =    I * (B * amp/2);
		}
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V2)(lx_k0, 0, 0)  =  -I * (A * amp/2);
			(*V3)(lx_k0, 0, 0)  =   complx(A * amp/2, 0.0);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{	
			(*V2)(lx_minus_k0, 0, 0) =   I * (A * amp/2);
			(*V3)(lx_minus_k0, 0, 0) =   complx(A * amp/2, 0.0);
		}	
	}
	
	else if ((basis_type == "SCFT") && (my_id == master_id))
	{
		cout << "ABC field not allowed for SCFT " << endl;
		exit(0);
	}
}



//*********************************************************************************************


void IncVF::Setup_SIX_MODE_field(int k0, DP amp101, DP amp011, DP amp112, DP h)
{

	*V1 = 0.0; 
	*V2 = 0.0;  
	*V3 = 0.0; 
	// initialize

	if (basis_type == "FOUR")
	{
		int lx_k0 = Get_lx("FOUR", k0, N);
		int lx_minus_k0 = Get_lx("FOUR", -k0, N);
		int ly_k0 = Get_ly3D("FOUR", k0, N);
		int ly_minus_k0 = Get_ly3D("FOUR", -k0, N);
		
		DP factor1 = 2.0/sqrt(2+h*h) * amp101;
		DP factor2 = 2.0/sqrt(2+h*h) * amp011;
		DP factor3 = 4.0/sqrt(6+10*h*h) * amp112;
		
		if (my_id == master_id)  
		{
			(*V1)(0, ly_k0, k0) = complx(-h*(factor2/4), 0.0);
			(*V1)(0, ly_minus_k0, k0) = complx(h*(factor2/4), 0.0);
			
			(*V2)(0, ly_k0, k0)  = (I) * (factor2/4);
			(*V2)(0, ly_minus_k0, k0) = (-I) * (factor2/4);
			
			(*V3)(0, ly_k0, k0)  = (-I) * (factor2/4);
			(*V3)(0, ly_minus_k0, k0) = (-I) * (factor2/4); 
		}
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) )
		{
			(*V1)(lx_k0, 0, k0)  = I * (factor1/4);
			(*V2)(lx_k0, 0, k0)  = complx(-h*(factor1/4), 0.0);
			(*V3)(lx_k0, 0, k0)  =  (-I) * (factor1/4);
			
			(*V1)(lx_k0, ly_k0, 2*k0)  = I * (factor3/8);
			(*V1)(lx_k0, ly_minus_k0, 2*k0) = I * (factor3/8);
			
			(*V2)(lx_k0, ly_k0, 2*k0)  = I * (factor3/8) + complx(h*factor3/4, 0.0);
			(*V2)(lx_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) + complx(h*factor3/4, 0.0);
			
			(*V3)(lx_k0, ly_k0, 2*k0)  = (-I) * (factor3/8) - complx(h*factor3/8, 0.0);
			(*V3)(lx_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) + complx(h*factor3/8, 0.0);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{
			(*V1)(lx_minus_k0, 0, k0)  = -I * (factor1/4);
			(*V2)(lx_minus_k0, 0, k0)  = complx(h*(factor1/4), 0.0);
			(*V3)(lx_minus_k0, 0, k0)  =  (-I) * (factor1/4);
			
			(*V1)(lx_minus_k0, ly_k0, 2*k0)  = (-I) * (factor3/8);
			(*V1)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8);	
					
			(*V2)(lx_minus_k0, ly_k0, 2*k0)  = I * (factor3/8) - complx(h*factor3/4, 0.0);
			(*V2)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) 
													- complx(h*factor3/4, 0.0);
			
			(*V3)(lx_minus_k0, ly_k0, 2*k0)  = (-I) * (factor3/8) + complx(h*factor3/8, 0.0);
			(*V3)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) 
													- complx(h*factor3/8, 0.0);
		}	
		
	}
	
	else if ((basis_type == "SCFT") && (my_id == master_id))
	{
		cout << "Setup_SIX_MODE_field Not allowed in SCFT basis_type " << endl;
		exit(0);
	}

}


//*******************************  End of ic_model.cc *****************************************




	
