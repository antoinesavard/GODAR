!=======================================================================
!             Common block: bond properties
!=======================================================================

    double precision :: 	 	& ! global physical parameters
				eb			 , 	& ! elastic stiffness
				lambda_rb	 , 	& ! radius parameter
				lambda_lb	 , 	& ! lenght parameter
				sigmatb_crit , 	&
				sigmacb_crit ,	&
				tau_crit	 ,	&
				gamma_d			  ! prevents spurious oscillations

	double precision ::			& ! individual physical properties
				knb 	(n,n), 	& ! bond normal stiffness
				ktb 	(n,n), 	& ! bond shear stiffness
				rb		(n,n),	& ! bond half width
				hb		(n,n),	& ! bond thicknesses
				lb		(n,n),	& ! bond lenght
				sb		(n,n),	& ! bond cross sectionnal area
				ib		(n,n),  & ! bond moment of inertia
                thetarelb(n,n), & ! bond relative angle
                deltanb (n,n),  & ! elongation
                deltatb (n,n)     ! deflection

	double precision ::			& ! forces in bonds
				fbn		(n,n),	& ! force due to floes sliding
				fbt		(n,n),	& ! force due to compresion/elongation
				mbb		(n,n),  & ! moment du to bending
				mb      (n)  , 	& ! total moment due to bonds on i
				fbx     (n)  ,  & ! forces in x
				fby     (n)       ! forces in y

	double precision ::			& ! stresses in bonds
				taub	(n,n),	&
				sigmacb	(n,n),	&
				sigmatb	(n,n)

    integer :: 					& ! bond presence or not
				bond 	(n,n)

    common/bond_bool/           &
                bond              ! bond between disks i and j  [0 or 1]

    common/bond_param/			& ! global physical parameters
                eb			 ,	& ! elastic modulus of bonds	 [N/m^2]
				lambda_rb	 , 	& ! radius parameter                 [-]
				lambda_lb	 , 	& ! lenght parameter                 [-]
				sigmatb_crit , 	&
				sigmacb_crit ,	&
				tau_crit	 ,	&
				gamma_d

	common/bond_var/			& ! individual physical properties
				knb 		 , 	& ! bond normal stiffness
				ktb			 , 	& ! bond shear stiffness
				rb			 ,	& ! bond half width
				hb			 ,	& ! bond thicknesses
				lb			 ,	& ! bond lenght
				sb			 ,	& ! bond cross sectionnal area
				ib		     ,  & ! bond moment of inertia
                thetarelb    ,  & ! bond relative angle
                deltanb      ,  & ! elongation
                deltatb           ! deflection

	common/bond_var/			& ! forces in bonds
				fbn			 ,	& ! force due to floes sliding
				fbt			 ,	& ! force due to compresion/elongation
				mbb		     ,  & ! moment du to bending
				mb			 ,	& ! total moment due to bonds on i
				fbx			 ,	& ! forces in x due to bonds
				fby			      ! forces in y due to bonds

	common/bond_var/			& ! stresses in bonds
				taub		 ,	& ! shear stress
				sigmacb		 ,	& ! compressive stress
				sigmatb           ! tensile stress
