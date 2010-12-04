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

/*! \file  init_cond_RB_misc.cc
 * 
 * @brief Initial conditions as RB convection: Lorenz equations.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FOUR basis; not for SCFT basis.
 *
 * @note \f$ V_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ V_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ V_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 */

#include "../IncFluid.h"


//**********************************************************************************************

/** @brief Set the modes according to Lorenz equations for SCFT mode only.
 * 
 * @return In 2D: set \f$ V_1(1,1), \theta(1,1), \theta(2,0) \f$.  Usually all real value <BR>
 *		Note that (1,-1) is not stored in our method.  
 *		\f$ V_2(1,1) \f$ determined using incompressibility cond.
 *
 * @return In 3D: \f$ V_1(1,0,1), \theta(1,0,1), \theta(2,0,0) \f$. Usually all real value <BR>
 *		Note that (1,0,-1) is not stored in our method.  
 *		\f$ V_2(1,0,1) \f$ determined using incompressibility cond.
 */
void  IncFluid::Init_cond_RB_Lorenz(IncSF& T)
{

	 DP		W101,T101,T200;
	W101 = (*init_cond_double_para)(1);
	T101 = (*init_cond_double_para)(2);
	T200 = (*init_cond_double_para)(3);
	  

	if (my_id == master_id) 
	{
		complx Vz;	
		
		(*V1)=0.0; 
		(*V2)=0.0;  
		(*V3)=0.0;
		(*T.F)=0.0;
	 
		(*V1)(1,0,1) = complex<DP>(W101,0); 
		(*V2)(1,0,1) = complex<DP>(0, 0);
		
		Last_component(1,0,1,(*V1)(1,0,1),(*V2)(1,0,1), Vz);	
		(*V3)(1,0,1) = Vz; 
		
		(*T.F)(2,0,0) = complex<DP>(T200,0); 
		(*T.F)(1,0,1) = complex<DP>(T101,0); 
		
		cout << "init cond V: Lorenz vars V1(101), V2(101), V3(101), T(101), T(200): " << endl 
			<< (*V1)(1,0,1)		<< " " 
			<< (*V2)(1,0,1)		<< " " 
			<< (*V3)(1,0,1)		<< " " 
			<< (*T.F)(1,0,1)	<< " " 
			<< (*T.F)(2,0,0) << endl;
	}		

}



//******************************** End of Init_cond_RB_misc.cc ********************************


