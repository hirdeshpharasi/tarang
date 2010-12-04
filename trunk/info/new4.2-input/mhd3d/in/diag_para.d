#diagnostic_para.d
#Only change the numerical fields -- Don't delete or add any line
# N[i]	
8	8	8


(1)solver-meta-para(int,float,string)		1	3	1
(2)Solver-para-int		2
(3)Solver-para-float		1	1	1
(4)Solver-para-string		TEST


Diagnostic-procedure	0	10000


(1)basis_type		FOUR
(2)alias_switch		DEALIAS
(3)Integration_scheme	RK4
(4)Input_number_mode	ASCII
(5)Output_no_mode	ASCII


(1)hyperdissipation 	1
(2)Low-dim		0
(3)LES			0
(4)Anisotropic-ring	0
(5)Anisotropic-cylinder-ring 	0
(6)structure-fn		0
(7)planar-structure-fn	0
(8)skpq-switch		0 
(9)Output-real-space-r 	1
(10)sincos-horizintal2d   0
(11)free-slip-3d	0
(12)helicity-flux-shells 0
(13)anisotropy-switch(1,2,3)		3
(14)waveno(Actual(0)_or_grid(1))	0
 

#diss_coeff + diffusion_coefficient (If RB -- ignore it..)
0.1	0.2

#Hyper_dissipation_switch + hyper dissipation coefficient
0	0


U.Tinit 	0
U.Tfinal 	0.01
U.Tdt		0.001	


(1)Ring-spectrum-input-scheme 	0
(2)Number-sectors-ring-spectrum	5
(3)Number-slabs-cyl 		5
(4)qmin 			1
(5)qmax 			10
(6-to-5+D)N_out_reduced		4  4  2


ET:(1)real-imag-switch	0
(2)no-spheres 		10
(3)no-shells 		10	
(4)no-ring-shells 	9
(5)no-sectors 		5
(6)cylinder-no-shells	6
(7)cylinder-no-slabs	5
(8)Shell-input-scheme		0
(9)Ring-shell-input-scheme	0
(10)Ring-ector-input-scheme	0
(11)Cylinder-shell-input-scheme	0
(12)Cylinder-slab-input-scheme	0


(1)no-output-k-modes	1
(2)no-output-r-points	1

1	0	1

1	0	1
	

(1)input-meta-para(int,float,string)		4	3	1
(2)input-para-int:field_input_proc,N_in_reduceced	3	64 	64	2
(3)input-para-float		1.0	1.0	1.0
(4)input-para-string		TEST


(1)force-meta-para(int,float,string)		1	2	1
(2)force-para-int:field_input_para,N_in_reduceced	0
(3)force-para-float		5.0	5.0
(4)force-para-string		TEST



