!=======================================================================
!             COMMAND BLOCK: bond properties
!=======================================================================

    double precision :: 		& ! global physical parameters
				eb			, 	&
				lambda_rb	, 	&
				lambda_lb	, 	&
				sigmat_max	, 	&
				sigman_max	,	&
				tau_max		,	&
				gamma_d

	double precision ::			& ! individual physical properties
				knb 	(n,n), 	&
				ktb 	(n,n), 	&
				rb		(n,n),	&
				hb		(n,n),	&
				lb		(n,n),	&
				sb		(n,n),	&
				ib		(n,n)

	double precision ::			& ! forces in bonds
				fbn		(n,n),	&
				fbt		(n,n),	&
				mbb		(n,n)

	double precision ::			& ! stresses in bonds
				taub	(n,n),	&
				sigmanb	(n,n),	&
				sigmatb	(n,n)

    integer :: 					& ! bond presence or not
				bond 	(n,n)

    common/bond_param/			& ! global physical parameters
                bond        , 	& ! bond between disks i and j  [0 or 1]
                eb			,	& ! elastic modulus of bonds	 [N/m^2]
				lambda_rb	, 	&
				lambda_lb	, 	&
				sigmat_max	, 	&
				sigman_max	,	&
				tau_max		,	&
				gamma_d

	common/bond_var/			& ! individual physical properties
				knb 		, 	&
				ktb			, 	&
				rb			,	&
				hb			,	&
				lb			,	&
				sb			,	&
				ib		

	common/bond_var/			& ! forces in bonds
				fbn			,	&
				fbt			,	&
				mbb		

	common/bond_var/			& ! stresses in bonds
				taub		,	&
				sigmanb		,	&
				sigmatb
