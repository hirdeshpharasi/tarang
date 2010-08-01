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

/*! \file  field_k_fns.cc
 * 
 * @brief Functions related to waven numbers.
 * @sa void IncVF::Last_component(int kx, int ky, complx &Vx, complx &Vy)
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  Add_complex_conj does not work across the processors.
 */

#include "../IncVF.h"
#include "../IncSF.h"

extern Uniform<DP> SPECrand;

//*********************************************************************************************

/** @brief Compute the last component Vy 
 * 
 *  @param kx, ky
 *	@param Vx
 *
 * @return \f$ V_y = -D_x V_x / K_y \f$ <BR>
 *  except when Ky =0 (x axis).  Here the read component is Vy, not Vx.  For this case
 *   \f$ V = (0, V_x) \f$ (input).
 */
void IncVF::Last_component(int kx, int ky, complx &Vx, complx &Vy)
{

	if (ky == 0) 
	{
		Vy = Vx; Vx = complex<DP>(0,0);				// for SCFT: Vx must be real
	}
	
	else 
		Vy = DxVx(basis_type, kfactor, kx, Vx) / complex<DP>(0, -ky*kfactor[2]);
}


//*********************************************************************************************

/** @brief Compute the last component Vz 
 * 
 *  @param kx, ky, kz
 *	@param Vx, Vy
 *
 * @return \f$ V_z = -(D_x V_x + D_y V_y) / K_z \f$ <BR>
 *		except when Kz =0 (xy plane).  Here the read components are Vx and Vz (in place of Vy).
 *		Obtain Vy using incompressibility criteria.
 */
void IncVF::Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz)
{
	complx temp;	
	
	if (kz != 0) 
	{
		temp = DxVx_plus_DyVy(basis_type, kfactor, kx, ky, Vx, Vy);  // Dx Vx + Dy Vy
		Vz = temp / complex<DP>(0, -kz*kfactor[3]);
	}
	
	// kz= 0; The read components are Vx and Vz	
	else if (N[3] > 2)			
	{
		Vz =  Vy;
		Last_component(kx, ky, Vx, Vy);
	}	
	
	// 2D emulation:  ignore read Vy.  Compute Vy from Vx using incompressible cond
	// Set Vz = 0. 
	else
	{
		if (ky != 0)
		{
			Vy = DxVx(basis_type, kfactor, kx, Vx) / complex<DP>(0, -ky*kfactor[3]);
			Vz = 0.0;
		}
		else
		{
			Vx = 0.0; 
			Vy = 0.0; 
			Vz = 0.0;
		}	
	}
}		


//*********************************************************************************************

/** @brief Compute Vy, Vz given \f$ V_x \f$ and \f$ \Omega_x \f$
 * 
 *  @param kx, ky, kz
 *	@param \f$ V_x , \Omega \f$
 *
 *  @return \f$ V_y = (D_y D_x V_x + D_x \Omega_x)/ K_{\perp}^2 \f$
 *  @return \f$ V_z = (D_z D_x V_x - D_x \Omega_x)/ K_{\perp}^2 \f$
 *  @return  Exception: Points on the kx axis (\f$ Ky=kz=0 \f$).  Here \f$ V =0 \f$.
 */
 void IncVF::Compute_VyVz
 (
	int kx, int ky, int kz, 
	complx& Vx, 
	complx Omega, 
	complx &Vy, complx &Vz
)
{
	complx delxvx, numr; 
	DP kperp_sqr, Kx, Ky, Kz;
	
	Kx = kx* kfactor[1]; 
	Ky = ky* kfactor[2]; 
	Kz = kz* kfactor[3];
	kperp_sqr = pow2(Ky) + pow2(Kz);
	
	if	( (ky == 0) && (kz == 0) )   		// kperp == 0
	{
		cout << "Problematic case-- ky = kz =0: So setting Vx = Vy = Vz = 0" << endl;
		Vx = Vy = Vz = complex<DP>(0, 0); 
	}
	
	else 
	{	
		delxvx = DxVx(basis_type, kfactor, kx, Vx);
		numr = complex<DP>(0, Ky)* delxvx + complex<DP>(0, Kz)* Omega;
		Vy = numr / kperp_sqr;
		
		numr = complex<DP>(0, Kz)* delxvx + complex<DP>(0, -Ky)* Omega;
		Vz = numr / kperp_sqr;
	}	
}

//*********************************************************************************************

/** @brief D=3: Add complex conjugate at \f$ \vec{k} \f$ 
*				once \f$ \vec{V}(kx, ky, kz) \f$ is given.
 *  
 *	@bug  Add_complex_conj does not work across the processors.
 * @note  For ky = 0 and ky = N(2)/2 planes only.
 * 
 * @return  FOUR: \f$ \vec{V}(-kx, -ky, kz) = \vec{V}^*(kx, ky, kz)\f$
 * @return  SCFT: For (ky != 0) -> \f$ \vec{V}(kx, -ky, kz) = \vec{V}^*(kx, ky, kz)\f$.
 * @return  SCFT: For (ky = 0) -> \f$ imag(V_y(kx, 0, kz)) = imag(V_z(kx, 0, kz)) = 0\f$,
 *					and \f$ V_x(kx, 0, kz) = 0 \f$.
 *				
 */
void IncVF::Add_complex_conj(int kx, int ky, int kz, complx Vx, complx Vy, complx Vz)			
{	
	
	int lx, ly, lz;

	// On kz=0 or kz=N[3]/2 plane
	if ((basis_type == "FOUR") && ((kz == 0) || (kz == N[3]/2)) )				
	{
		lx = Get_lx(basis_type, -kx, N);  
		ly = Get_ly3D(basis_type, -ky, N); 
		lz = kz;	
		
		if ( (lx >= 0) && (lx < local_N1) ) 
		{
			(*V1)(lx, ly, lz) = conj(Vx);		
			(*V2)(lx, ly, lz) = conj(Vy);
			(*V3)(lx, ly, lz) = conj(Vz);
			
			cout << "Complex-conj(V) added for k = (" << -kx << "," << -ky << "," << kz << ")" 
					<< endl;	
		}

	}	
		
	else if ((basis_type == "SCFT") && ((kz == 0) || (kz == N[3]/2)) ) 		
	{
		lx = Get_lx(basis_type, kx, N);   
		ly = Get_ly3D(basis_type, -ky, N);   
		lz = kz;
		
		if ( (lx >= 0) && (lx < local_N1) ) 
		{
			if (ky != 0)
			{
				(*V1)(lx, ly, lz) = conj(Vx);		
				(*V2)(lx, ly, lz) = conj(Vy);
				(*V3)(lx, ly, lz) = conj(Vz);
				
				cout << "Complex-conj(V) added for k = (" << kx << "," << -ky << "," 
						<< kz << ")" << endl;	
			}
			
			else		// ky =  0
			{
				(*V1)(lx, ly, lz) = 0.0;		
				imag((*V2)(lx, ly, lz)) = 0.0;
				imag((*V3)(lx, ly, lz)) = 0.0;
				
					cout << "imag(V2, V3) set to zero for k = (" << kx <<  ",0,0 ); & V1 = 0.0" 
						<< endl;
			}		
		}
	}
	
}	

//*********************************************************************************************

/** @brief D=3: Add complex conjugate at \f$ \vec{k} \f$ 
 *				once \f$ F(kx, ky, kz) \f$ is given.
 *  
 * @bug  Add_complex_conj does not work across the processors.
 * @note  For ky = 0 and ky = N(2)/2 planes only.
 * 
 * @return  FOUR: \f$ F(-kx, -ky, kz) = F^*(kx, ky, kz)\f$
 * @return  SCFT: For (ky != 0) -> \f$ F(kx, -ky, kz) = F^*(kx, ky, kz)\f$.
 * @return  SCFT: For (ky = 0) -> \f$ imag(F(kx, 0, kz))  = 0\f$,
 *				
 */
 void IncSF::Add_complex_conj(int kx, int ky, int kz, complx G)					
{	

	int lx, ly, lz;
	
	// On kz=0 or kz=N[3]/2 plane
	if ((ISF_basis_type == "FOUR") && ((kz == 0) || (kz == NIs[3]/2)) ) 			
	{
		lx = Get_lx(ISF_basis_type, -kx, NIs);  
		ly = Get_ly3D(ISF_basis_type, -ky, NIs); 
		lz = kz;	
		
		if ( (lx >= 0) && (lx < local_N1) )		
			(*F)(lx, ly, lz) = conj(G);	
	//	cout << "Complex-conj(G) added for k = " << -kx << "," << -ky << "," << kz << ")" 
	//			<< endl;	
	}	
	
	else if ((ISF_basis_type == "SCFT") && ((kz == 0) || (kz == NIs[3]/2)) ) 	
	{
		lx = Get_lx(ISF_basis_type, kx, NIs); 
		ly = Get_ly3D(ISF_basis_type, -ky, NIs);   
		lz = kz;
		
		if ( (lx >= 0) && (lx < local_N1) )		
		{
			if (ky != 0)
			{
				(*F)(lx, ly, lz) = conj(G);	
				
				cout << "Complex-conj(G) added for k = (" << -kx << "," << -ky << "," 
						<< kz << ")" << endl;
			}	
			
			else // ky =  0	
			{
				imag((*F)(lx, ly, lz)) = 0.0;
				
				cout << "imag(G) set to zero for k = (" << kx << ", 0, 0)" << endl;
			}
		}	
	}
}



//*********************************************************************************************

/** @brief D=3:  \f$ \vec{V}(kx, ky) = (V_x, V_y, V_z) \f$ and add complex conjugate.
 * 
 * @bug  Assign_field_comp_conj does not work across the processors.
 * @sa Add_complex_conj(int kx, int ky, int kz, complx Vx, complx Vy, complx Vz)
 */
void IncVF::Assign_field_comp_conj(int kx, int ky, int kz, complx Vx, complx Vy, complx Vz)
{
	
	int lx = Get_lx(basis_type, kx, N);	
	int ly = Get_ly3D(basis_type, ky, N); 
	int lz = kz;	
	
	if ( (lx >= 0) && (lx < local_N1) ) 
	{																																			
		(*V1)(lx, ly, lz) = Vx;		
		(*V2)(lx, ly, lz) = Vy; 
		(*V3)(lx, ly, lz) = Vz;
	}
	
	Add_complex_conj(kx, ky, kz, Vx, Vy, Vz);	
}


//*********************************************************************************************

/** @brief D=3:  \f$ F(kx, ky, kz) = G \f$ and add complex conjugate.
 *  
 * @bug  Assign_field_comp_conj does not work across the processors.
 * @sa Add_complex_conj(int kx, int ky,int kz,  complx G)
 */
void IncSF::Assign_field_comp_conj(int kx, int ky, int kz, complx G)
{
	int lx = Get_lx(ISF_basis_type, kx, NIs);	
	int ly = Get_ly3D(ISF_basis_type, ky, NIs); 
	int lz = kz;	
	
	if ( (lx >= 0) && (lx < local_N1) )   
		(*F)(lx, ly, lz) = G;	
							
	Add_complex_conj(kx, ky, kz, G);		
}



//*********************************************************************************************
/** @brief D=3: Add complex conjugate at \f$ \vec{k} \f$ 
 *				once force field \f$ \vec{F}(kx, ky, kz) \f$ is given.
 *  
 * @bug Does not work across the processors.
 *
 * @note  For ky = 0 and ky = N(2)/2 planes only.
 * 
 * @return  FOUR: \f$ \vec{F}(-kx, -ky, kz) = \vec{F}^*(kx, ky, kz)\f$
 * @return  SCFT: For (ky != 0) -> \f$ \vec{F}(kx, -ky, kz) = \vec{F}^*(kx, ky, kz)\f$.
 * @return  SCFT: For (ky = 0) -> \f$ imag(F_y(kx, 0, kz)) = imag(F_z(kx, 0, kz)) = 0\f$,
 *					and \f$ F_x(kx, 0, kz) = 0 \f$.
 *				
 */
void IncVF::Add_complex_conj_force(int kx, int ky, int kz, complx Fx, complx Fy, complx Fz)	
{	
	
	int lx, ly, lz;

	if ((basis_type == "FOUR") && ((kz == 0) || (kz == N[3]/2)) )				
	{
		lx = Get_lx(basis_type, -kx, N);  
		ly = Get_ly3D(basis_type, -ky, N); 
		lz = kz;	
		
		if ( (lx >= 0) && (lx < local_N1) ) 
		{
			(*Force1)(lx, ly, lz) = conj(Fx);		
			(*Force2)(lx, ly, lz) = conj(Fy);
			(*Force3)(lx, ly, lz) = conj(Fz);
			
			// cout << "Complex-conj(Force) added for k = " << -kx << "," << -ky << "," << kz 
			//		<< ")" << endl;	
		}
	}
			
	else if ((basis_type == "SCFT") && ((kz == 0) || (kz == N[3]/2)) ) 		
	{
		lx = Get_lx(basis_type, kx, N);   
		ly = Get_ly3D(basis_type, -ky, N);   
		lz = kz;
		
		
		if ( (lx >= 0) && (lx < local_N1) ) 
		{
			if (ky != 0)
			{
				(*Force1)(lx, ly, lz) = conj(Fx);		
				(*Force2)(lx, ly, lz) = conj(Fy);
				(*Force3)(lx, ly, lz) = conj(Fz);
				
				// cout << "Complex-conj(Force) added for k = " << kx << "," << -ky << ","  
				//		<< kz << ")" << endl;		
			}
			
			else		// ky =  0
			{
				(*Force1)(lx, ly, lz) = 0.0;		
				imag((*Force2)(lx, ly, lz)) = 0.0;
				imag((*Force3)(lx, ly, lz)) = 0.0;
				
				//	cout << "imag(Force2, Force3) set to zero for k = (" << kx <<  "0,0 );   
				//				& V1 = 0.0" << endl;
			}		
		}
	}
	
}	

//
// Scalar 
//

void IncSF::Add_complex_conj_force(int kx, int ky, int kz, complx ForceG)						
{	

	int lx, ly, lz;
	
	if ((ISF_basis_type == "FOUR") && ((kz == 0) || (kz == NIs[3]/2)) ) 			
	{
		lx = Get_lx(ISF_basis_type, -kx, NIs);  
		ly = Get_ly3D(ISF_basis_type, -ky, NIs); 
		lz = kz;	
		
		if ( (lx >= 0) && (lx < local_N1) )		
			(*Force)(lx, ly, lz) = conj(ForceG);	
			
		//  cout << "Complex-conj(G) added for k = " << -kx << "," << -ky << "," << kz 
		//			<< ")" << endl;	
	}	
	
	else if ((ISF_basis_type == "SCFT") && ((kz == 0) || (kz == NIs[3]/2)) ) 	
	{
		lx = Get_lx(ISF_basis_type, kx, NIs); 
		ly = Get_ly3D(ISF_basis_type, -ky, NIs);   
		lz = kz;
		
		if ( (lx >= 0) && (lx < local_N1) )	
		{
			if (ky != 0)
			{	
				(*Force)(lx, ly, lz) = conj(ForceG);	
			
				// cout << "Complex-conj(ForceG) added for k = " << -kx << "," << -ky << "," 
				//		<< kz << ")" << endl;
			}
			
			else  // ky = 0
			{
				imag((*Force)(lx, ly, lz)) = 0.0;
				
				cout << "imag(Force of scalar) set to zero for k = (" 
					 << kx << "," << ky  << "," << kz << ")" << endl;
			}
		}	
	}
}


//*********************************************************************************************

/** @brief D=3:  force field \f$ \vec{F}(kx, ky, kz) = (F_x, F_y, F_z) \f$ 
 *				and add complex conjugate.
 *  
 * @bug Does not work across the processors.
 *
 * @sa Add_complex_conj_force(int kx, int ky, kz, complx Fx, complx Fy, complx Fz)
 */
void IncVF::Assign_force_and_comp_conj(int kx, int ky, int kz, complx Fx, complx Fy, complx Fz)
{
	
	int lx = Get_lx(basis_type, kx, N);	
	int ly = Get_ly3D(basis_type, ky, N); 
	int lz = kz;	
	
	if ( (lx >= 0) && (lx < local_N1) ) 
	{																																			
		(*Force1)(lx, ly, lz) = Fx;		
		(*Force2)(lx, ly, lz) = Fy; 
		(*Force3)(lx, ly, lz) = Fz;
	}
	
	Add_complex_conj(kx, ky, kz, Fx, Fy, Fz);	
}

//
//
//*********************************************************************************************

/** @brief D=2:  force field for scalar \f$ F(kx, ky, kz) = FG \f$ and add complex conjugate.
 *  
 * @bug Does not work across the processors.
 *
 * @sa Add_complex_conj_force(int kx, int ky, int kz, complx ForceG)
 */
void IncSF::Assign_force_and_comp_conj(int kx, int ky, int kz, complx FG)
{
	int lx = Get_lx(ISF_basis_type, kx, NIs);	
	int ly = Get_ly3D(ISF_basis_type, ky, NIs); 
	int lz = kz;	
	
	if ( (lx >= 0) && (lx < local_N1) )   
		(*Force)(lx, ly, lz) = FG;	
							
	Add_complex_conj(kx, ky, kz, FG);		
}



//*********************************************************************************************
//*********************************************************************************************

/// 3D: Return Tk = Real(-nlin(k). conj(V(k))  for vector V
DP IncVF::Get_Tk(int kx, int ky, int kz)
{
	int lx = Get_lx(basis_type, kx, N);	
	int ly = Get_ly3D(basis_type, ky, N); 
	int lz = kz;	

	if ( (lx >= 0) && (lx < local_N1) ) 
	{
		return (	-real((*nlin1)(lx, ly, lz) * conj((*V1)(lx, ly, lz)))
				-real((*nlin2)(lx, ly, lz) * conj((*V2)(lx, ly, lz)))
				-real((*nlin3)(lx, ly, lz) * conj((*V3)(lx, ly, lz))) );	
	}	
	
	else 
		return 0; // for -Wall
}



//
//
/// 3D:  Return Tk = Real(-nlin(k). conj(F(k)) for scalar T
DP IncSF::Get_Tk(int kx, int ky, int kz)
{
	int lx = Get_lx(ISF_basis_type, kx, NIs);	
	int ly = Get_ly3D(ISF_basis_type, ky, NIs); 
	int lz = kz;	
		
	if ( (lx >= 0) && (lx < local_N1) )	
			return ( -real((*nlin)(lx, ly, lz) * conj((*F)(lx, ly, lz))) );
			
	else
		return 0;  // for -Wall
}


//*********************************************************************************************
//*********************************************************************************************

/// 3D: Generate random vector at (kx,ky,kz) with range = rand_range
void IncVF::Get_random_vector
(
	int kx, int ky, int kz, 
	DP rand_range, 
	complx& Vx, complx& Vy, complx& Vz
)
{

	if (kz != 0) 
	{
		Vx = complx( 2*rand_range*(SPECrand.random()-0.5), 
					 2*rand_range*(SPECrand.random()-0.5) );
					 
		Vy = complx( 2*rand_range*(SPECrand.random()-0.5), 
					 2*rand_range*(SPECrand.random()-0.5) );
		
		Last_component(kx, ky, kz, Vx, Vy, Vz);
	}
	
	else 				// kz= 0;
	{
		Vx = complx( 2*rand_range*(SPECrand.random()-0.5), 
					 2*rand_range*(SPECrand.random()-0.5) );
					 
		Last_component(kx, ky, Vx, Vy);	// div = on in kz=0 plane
		
		Vz = complx( 2*rand_range*(SPECrand.random()-0.5), 
					 2*rand_range*(SPECrand.random()-0.5) );
	}	
}

//
// SCALAR
// 
/// 3D: Generate random vector at (kx,ky,kz) with range = rand_range 
void IncVF::Get_random_scalar(int kx, int ky, int kz, DP rand_range, complx& F)
{
	F = complx( 2*rand_range*(SPECrand.random()-0.5), 2*rand_range*(SPECrand.random()-0.5) );
}



//**********************************   End of field_k_fns.cc  *********************************



