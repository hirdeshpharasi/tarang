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

/*! \file  EnergyTr.h
 * 
 * @brief  Class declaration of energy transfer vars and functions. 
 * 
 *	We omputes energy transfers from region A1 of Giver field G, to the receiver field R
 *	with the help of helper field H.
 *
 *	The Giver field is filled in k-space appropriately (region A1).  <BR>
 *	We do the same for the receiver field (region A2). <BR>
 * 
 *	Definitions:
 *		\f$ T^{G,A1}_{R,A2} = Sign \sum_{A2} R_i^*(\vec{K})  N_i(\vec{K}) \f$ <BR>
 *		where the nonlinear term 
 *		\f$ \mathcal{N} = \mathcal{F}( (\vec{H}.\nabla) \vec{G}^{A1} ) \f$. <BR>
 *
 *	For u-to-u:	Sign = -; H, G, R = \f$ \vec{V} \f$. <BR>
 *	For b-to-b:	Sign = -; H = \f$ \vec{V} \f$; G, R = \f$ \vec{B} \f$. <BR> 
 *	For u-to-b: Sign = +; G = \f$ \vec{V} \f$; H, R = \f$ \vec{B} \f$. <BR>
 *  For b-to-u: Sign = +; R = \f$ \vec{V} \f$; G, H = \f$ \vec{B} \f$. <BR>
 *	
 *	For \f$\psi\f$ to \f$ \psi \f$: Sign = -; H = \f$ \vec{V} \f$, and G,R= \f$ \psi \f$.
 *
 *	
 *	Shell radii 0,r(1),r(2), .., r(n-1), r(noshells)=Infty. <BR>
 *	Shell(n) = (r(n-1),r(n)]. <BR> 
 *	There are n+1 shells with Shell(0) being the origin. <BR>
 *	The last shell contains modes beyond the largest sphere that fits inside the box.
 *
 *	For rings, we divide angles into sectors, with angles being 
 *		\f$ 0, \theta_{max}/n, 2*\theta_{max}/n, .., \theta_{max} \f$. <BR>
 *	thetamax can be either \f$ \pi \f$ or \f$ pi/2 \f$.
 *	Sector(n) = \f$ (\theta(n-1), \theta(n)] \f$.  <BR>
 *	Sector(1) contains \f$ \theta=0 \f$ points as well.
 *	Sector(0) does not contain any point.
 *
 *	For cylinderical rings, we divide the anisotropic directions into noslabs slabs.
 *	Slab(n) = (H(n-1), H(n)] with H(n)=0, .., kmax. 
 *
 *	ET_shell_radii_sector_array() -- entries <BR>
 *		Shell: 0, r(1), r(2), .. r(n-1); Inf skipped.  <BR>
 *		Ring-shell: 0, r(1),  r(2), .. r(n-1); Inf skipped.  <BR>
 *		Ring-sector: 0, s(1), s(2), ... s(n).  <BR>
 *		cylinder-ring-shell: 0, r(1), r(2), ..., r(n-1); Inf skipped.<BR>
 *
 *	Shell_to_shell_transfer is computed for shells 1:noshells. Note that the last shell
 *	contains modes outside the largest sphere.
 *
 *	For ring-to-ring transfer we ignore the last shell since the ring-to-ring transfer
 *	takes too much time.  <BR>
 *	For ring-to-ring transfer it is better to use-defined r(i).
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec 2008
 *
 * @bug   No known bugs
 */


//*********************************************************************************************

#ifndef _EnergyTr
#define _EnergyTr

#include "../fields/fields.h"

const int MAX_NO_SKPQ_TRIADS = 40;			// for S(k,p,q) calculataions

//*********************************************************************************************

class EnergyTr  
{ 
 public:
	//! G contains energy giver in triad interaction
	Array<complx,3> *G1;   
	Array<complx,3> *G2;  
	Array<complx,3> *G3;   
	

public:

	int		NET[4];
	
	
	string	ET_basis_type;
	string	ET_alias_switch;
	
	// wavenumber K[i] = ET_kfactor*k[i]; k[i] = index
	DP		ET_kfactor[4];						
	
	//! Number of spheres
	int		no_spheres;							

	//!  Number of shells
	int		no_shells;							
	
	//! No of shells in ring transfer
	int		no_ring_shells;						
	
	//!  Number of sectors for energy transfers
	int		no_sectors_ring_tr;					
	
	//! No of cylinderical shells
	int		no_cylinder_shells;					
	
	//! No of kpll slabs
	int		no_cylinder_kpll_slabs;	
	
	
	//!  output separately real(V)*real(nlin*) and imag(V)*imag(nlin*) 
	//!  0: no, 1: yes
	int		ET_real_imag_switch;
	
	int		ET_anisotropic_ring_switch;
	int		ET_anisotropic_cylinder_switch;
	int		ET_Skpq_switch;
	
	
	//! Array containing the radius of spheres 
	Array<DP,1>  *sphere_radius;				

	//! shell_n = ( shell_radius(n-1),shellradis(n) ]
	Array<DP,1>  *shell_radius;					
	
	//! ring_shell_n = ( ring_shell_radius(n-1), ring_shell_radis(n) ]
	Array<DP,1>  *ring_shell_radius;			
	
	//! ring_n = ( sector_angle(n-1),sector_angle(n) ];
	//! the last ring includes both angles
	Array<DP,1>  *sector_angle_ring_tr;			
  
	Array<DP,1>  *cylinder_shell_radius;
	Array<DP,1>	 *cylinder_kpll_array_tr;
	
	//! for energy transfer to be passed Shell_mult_all
	Array<DP,1>  *temp_shell_tr;				
	
	//! for energy transfer to be passed Shell_mult_all
	Array<DP,2>  *temp_ring_tr;					

	//! for energy transfer to be passed Shell_mult_all
	Array<DP,2>  *temp_cylinder_ring_tr;		

	
	// S(k,p,q) calculations
	
	int				no_of_Skpq_triads;
	
	Array<int,2>	*triad_array;					// grid index e.g., 1,2,...
	Array<DP,1>		*Sself_array;
	Array<DP,1>		*S_to_VF_array;
	Array<DP,1>		*S_SF_array;

	// fluxes
	
	Array<DP,1>  *flux_self;	
	Array<DP,1>  *flux_VF_in_in;				// from Ufrom to extern field in_to_in
	Array<DP,1>  *flux_VF_in_out;				// from Ufrom to extern field in_to_out
	Array<DP,1>  *flux_VF_out_out;				// from Ufrom to extern field out_to_out
	Array<DP,1>  *flux_Elsasser;				//  From Z to Z
	Array<DP,1>  *flux_SF;						//  T< to T>
	
	Array<DP,1>  *flux_hk;
	
	// force*V
	
	Array<DP,1>  *forceV_shell;					// fu(k)*u(k)
	Array<DP,1>  *forceSF_shell;				// fT(k)*T(k)

	
	// Isotropic shell-to-shell
	
	Array<DP,2>  *shelltoshell_self;
	Array<DP,2>  *shelltoshell_VF;				// from Ufrom to extern VF B
	Array<DP,2>  *shelltoshell_Elsasser;		// from Z to Z
	Array<DP,2>  *shelltoshell_SF;				// theta->theta
	
	Array<DP,2>  *shelltoshell_hk;
  
	// Ring to ring
	
	Array<DP,4>  *ring_to_ring_self;
	Array<DP,4>  *ring_to_ring_VF;				// from Ufrom to extern VF B
	Array<DP,4>  *ring_to_ring_Elsasser;			// from Z to Z
	Array<DP,4>  *ring_to_ring_SF;				// theta->theta
	
	
	Array<DP,2>  *forceV_ring;					// fu(k)*u(k)
	Array<DP,2>  *forceSF_ring;					// fT(k)*T(k)
	
	
	// Cylinderical ring to ring
	
	Array<DP,4>  *cylinder_ring_to_ring_self;
	Array<DP,4>  *cylinder_ring_to_ring_VF;				// from Ufrom to extern VF B
	Array<DP,4>  *cylinder_ring_to_ring_Elsasser;		// from Z to Z
	Array<DP,4>  *cylinder_ring_to_ring_SF;				// theta->theta
	
	Array<DP,2>  *forceV_cylinder_ring;					// fu(k)*u(k)
	Array<DP,2>  *forceSF_cylinder_ring;					// fT(k)*T(k)
	
	
	// B0 effect
	//! Energy transfer from b to u due to B0.
	//! \f$ \sum (\vec{B0} \cdot \vec{K}) \Im(\vec{u}(\vec{K'}) \vec{b}^*(\vec{K})) \f$
	Array<DP,1>	*energy_tr_shell_B0;
	Array<DP,2>	*energy_tr_ring_B0;
	Array<DP,2>	*energy_tr_cylinder_ring_B0;	
	
public:
	EnergyTr
	(
		int *NN, 
		string string_switches[], 
		Array<int,1> switches, 
		DP *prog_kfactor, 
		Array<int,1> ET_parameters, 
		Array<DP,1> ET_shell_radii_sector_array
	);

};

#endif

//***************************   End Class decl of EnergyTr   **********************************


 

