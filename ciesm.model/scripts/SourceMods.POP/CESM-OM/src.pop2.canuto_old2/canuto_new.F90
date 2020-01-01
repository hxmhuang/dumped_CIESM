!-----------------------------------------------------------------------
!     start of subroutine nasa_giss_mixing_column 
!-----------------------------------------------------------------------
      subroutine nasa_giss_mixing_column(Coriol, na, z, &
                        t, s, an2, s2, ri, rid, &
                        ustar_, buoytur, buoysol, &
                        amld, akm, akh, aks, mytid)


      use const_for_nasagissmixing 
      !use param_mod, only: mytid
      implicit none
      save


            integer, intent(in) :: mytid
            real(r8),   intent(in)  :: Coriol !cry_modify 
            integer,    intent(in)  :: na !cry_modify 
            real(r8),   dimension(na+1),  intent(in)  ::    z !cry_modify
            real(r8),   dimension(na),    intent(in)  ::    t,s,an2,rid !cry_modify
            real(r8),   dimension(na),    intent(inout)     ::    s2,ri
            real(r8),   intent(in)  ::    ustar_,buoytur,buoysol !cry_modify
            real(r8),   intent(out) ::    amld !cry_modify
            real(r8),   dimension(nmax),  intent(out) :: akm,akh,aks !cry_modify

!-----------------------------------------------------------------------
!     added by hxs and cry 201709141200
!-----------------------------------------------------------------------
            integer :: i,ibg
            integer :: idifs,idif
            integer :: iri,irisign,iristep
            integer :: irid,iridsign,iridstep
            integer :: ind2on2,ind2on2sign,ind2on2step
            integer :: ipoint,icall,iproblem,inegproblem
            integer :: ira_r,ipenra_r,itheta_r,itheta_r0,itheta_r1
            integer :: ifunreal,ifbelow,ifnofsmall,ifpureshear,ifrafglt

            integer :: jtheta_r,jtheta_r0,jtheta_r1

            integer :: kbot

            integer :: mt0s,mtm1s
            integer :: m,imax,ilmldpoint,ilmldpointneg
            
            real (r8) :: al0deep          
            real (r8) :: al,al0,al2
            real (r8) :: amtaun2,akz
            real (r8) :: an,and2,and2on2
            real (r8) :: anlq2,anlq2_back
            real (r8) :: al_back,al0_back,al2_back


            real (r8) :: b1,buoytot
            real (r8) :: back_ra_r1,back_rit1,back_ric1

            real (r8) :: c_y001

            real (r8) :: deltanum,deltaden,delta
            real (r8) :: deltheta_r1,delback_ra_r
            real (r8) :: delra_r,deltheta_r,delsisa_r
            real (r8) :: delz,delrh,del2rh,dzrh,d2zrh
            real (r8) :: delsm_back,delsh_back,delss_back,delslq2_back
            real (r8) :: dsm_back_o_dtheta,dsh_back_o_dtheta,dss_back_o_dtheta
            real (r8) :: dsisa_o_dtheta,dback_ra_r_o_dtheta,dslq2_back_o_dtheta


            real (r8) :: eplatidepend,eplatidepend_

            real (r8) :: rimax
            real (r8) :: ri1,back_ri1
            real (r8) :: rit,rit0,rit1
            real (r8) :: ric,ric0,ric1
            real (r8) :: rid1,back_rid1
            real (r8) :: ri_r,ri_r0,ri_r1
            real (r8) :: ra_r,ra_r0,ra_r1
            real (r8) :: rid_r,rid_r0,rid_r1
            real (r8) :: rrcrn,rrcrp,rdzlndzrh
            real (r8) :: rl,rlomega,rlblackadar

            real (r8) :: sm,sh,ss
            real (r8) :: sisa,sisa1
            real (r8) :: sm_r,sh_r,ss_r
            real (r8) :: sm_r0,sh_r0,ss_r0
            real (r8) :: smosloq_r,smosloq_0r
            real (r8) :: smb_00,shb_00,ssb_00
            real (r8) :: sm_last,sh_last,ss_last
            real (r8) :: sm_back,sh_back,ss_back,s2_back
            real (r8) :: slq2,slq2_r,slq2_r0,slq2_back,slq2b_00


            real (r8) :: tmp,tmp_back
            real (r8) :: theta_rcrn,theta_rcrn_deg
            real (r8) :: theta_rcrp,theta_rcrp_deg
            real (r8) :: theta_r,theta_r_deg,theta_r0,theta_r1

            external eplatidepend_
!-----------------------------------------------------------------------
!     added by hxs and cry 201709141200
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     level 2 turbulence models to be used in ogcm
!     sm and sh depends only on the richardson number
!     In salinity model case level 2 means S_M,S_H,S_C depend only on Ri,Ri_d
!-----------------------------------------------------------------------
!     inputs:
!     
!     nmax       = total number of vertical grid points
!     n          = number of vertical grid points above sea bottom
!                  minus one
!     Need the level number of the bottom ocean level for enhanced bottom mixing schemes.
!     kbot       = number of vertical grid points above sea bottom (= n+1)
!
!	  EXTEND TO CASES OF ONE AND ZERO OCEAN LEVEL. (n = 0,-1)
!     na         = Max(1,number of vertical grid points above sea bottom -1)
!     nb         = Max(0,number of vertical grid points above sea bottom -1)
!
!	  Increase number of vertical levels passed in by 1 because need to know 
!	  the bottom depth for purposes of enhanced bottom mixing option.
!     z(n+1)     = depth of vertical grid (cm, positive)
!
!     t(n)       = mean temperature (c)
!	  s(n)       = mean salinity (parts per thousand minus 35)
!	  rh(n)      = mean density (g cm**(-3))
!     ri(n)      = richardson number
!	  rid(n)     = temperature component minus salinity component of Ri
!     s2(n)      = (dU/dz)**2+(dV/dz)**2, shear squared, (1/s)**2
!	  Shear**2 units were mislabelled in turb_2, they should be frequency^2.
!     v_back(n)  = background turbulent viscosity (cm**2/sec)
!     t_back(n)  = background turbulent heat diffusivity (cm**2/sec)
!     s_back(n)  = background turbulent salt diffusivity (cm**2/sec)
!	  Add square of Brunt Vaisala frequency to the inputs.
!	  an2(n)     = N^2 (1/s)**2
!	  s_back(n) becomes a quantity recalculated in the subroutine
!	  when ifsalback>0.

!	  v_back,t_back,s_back become quantities calculated in the subroutine 
!     when ifback or ifsalback > 2

!	  SURFACE FLUXES:
!     ustar_        = sfc friction velocity                    (cm/s)
!     buoytur       = sfc turbulent buoyancy forcing         (cm2/s3)
!     buoysol       = sfc radiative buoyancy forcing         (cm2/s3)

!	  Switch for use of surface fluxes:
!    	isurfuse  = 0 for don't use surface fluxes for dimensionalization of turbulence model. 
!		    		  Rely on l_MLD^2 S, where l_MLD is the Blackadar lengthscale calculated using MLD.
!              	  = 1 for use total surface buoyancy, buoytot=buoytur+buoysol, 
!                  	  to calculate dissipation, epsilon. epsilon=-buoytot/(1-(1/2)yS_M). y=(\tau S)^2.


!	  CORIOLIS PARAMETER
!     Coriol	 = 2 \Omega sin(latitude) 		     	(1/s)	


!	  Choice of rotation influence on lengthscale.
!     ilomega    = 0 for l=l_Blackadar traditional Blackadar lengthscale without rotation
!
!                = 1 for l=(l_Blackadar^{-1} + l_\Omega^{-1})^{-1} ; when MLD>=500meters
!	                 l_\Omega \equiv \sqrt{-B*/f^3} for B*<0, and \infinity for B*>0 .
!     
!
!     fricmx     = max viscosity (cm**2/sec)
!     wndmix     = min viscosity at bottom of 1st level to simulate
!                  missing high frequency windstress components
!                  (cm**2/sec)
!     visc_cbu_limit = largest "visc_cbu" in regions of gravitational
!                      instability (cm**2/sec)
!     diff_cbt_limit = largest "diff_cbt" in regions of gravitational
!                      instability (cm**2/sec)
!
!	  For ifextermld = {0,1} calculate MLD {here, externally in calling unit).
!		  ifextermld = MLD switch [0 or 1]

!	  For ifoutput={0,1} don't do output fields to a file.

!	  For icondear = 0 K_X = (1/2) B_1 {l_Deardorff}^2 S/(Sl/q)_Deardorff S_X_P=eps (traditional)
!	      icondear =-1 K_X = (1/2) B_1 {l_Blackadar}^2 S/{(Sl/q)_P=eps} S_X_P=eps (No Deardorff)
!	      icondear =+1 K_X = (1/2) B_1 {l_Deardorff}^2 S/{(Sl/q)_P=eps} S_X_P=eps (Deardorff length only) 


! 	  For ifepson2=0  K_X = (1/2) B_1 l^2 S/(Sl/q) S_X
!  	  Background  ifepson2=1  K_X = (1/2) B_1^2 Ri (Sl/q)^2 (\epsilon/N^2) S_X
! 	  Bg.&deep Fg. ifepson2=2  K_X = (1/2) B_1^2 Ri (Sl/q)^2 (\epsilon/N^2) S_X
!	  Do NOT use Deardorff lengthscale limitation with (\epsilon/N^2) case.


!	  For ifdeeplat=0	pelagic background (\epsilon/N^2) constant(for ifepson2>0)
!	      ifdeeplat=1      pelagic bg.(\epsilon/N^2)~L(latitude,N) fr. Gregg et al.03

!
!	  FOR ifrafgmax = {0,1} {Don't,Do} limit background ra_r to be at maximum
!		  foreground ra_r at \theta_r where turbulence exists as Ri=> +infinity.

!		For ifsalback=3 case, v_back(n),t_back(n),s_back(n) become outputs.



!     outputs:
!     amld       = mixed layer depth (cm)
!     akm(nmax)  = turbulent viscosity (cm**2/sec)
!	  akh(nmax)  = turbulent heat diffusivity (cm**2/sec)
!	  aks(nmax)  = turbulent salt diffusivity (cm**2/sec)
!  
!     internal quantities:
!
!     nmodel     = 1: improved second order closure model
!                  2: parameter-free stochastic model
!     ntbl       = number of points in the look-up table
!     ria(ntbl)  = richardson number in the look-up table
!     slq2(ntbl) = shear number squared (s*l/q)**2 in the look-up table
!     sma(ntbl)  = sm in the look-up table
!     sha(ntbl)  = sh in the look-up table

!	  Quantities for the temperature-salinity model. Enlarge the table.
!     rib(-mt:mt)      = Ri for the 2D look-up table with enlarged dimensions
!     ridb(-mt:mt)     = difference of Tem and Sal Ri's for the 2D look-up table
!     slq2b(-mt:mt,-mt:mt) = shear number squared (s*l/q)**2 in the 2D table
!     smb(-mt:mt,-mt:mt)   = sm in the look-up table
!     shb(-mt:mt,-mt:mt)   = sh in the look-up table
!     ssb(-mt:mt,-mt:mt)   = ss in the look-up table



!     ri1        = Ri at the given level to be used in table interpolation
!     rid1       = Ri_d at the given level to be used in table interpolation
!     rimax      = maximum richardson number
!     visc_cbu_limit = largest "visc_cbu" in regions of gravitational
!                      instability (cm**2/sec)
!     diff_cbt_limit = largest "diff_cbt" in regions of gravitational
!                       instability (cm**2/sec)
!
!	  Quantities for dimensionalization of foreground diffusivities using surface forcing.
!	  buoytot	 = total buoyancy forcing into ocean (includes penetrative radiation) (cm**2/sec***3)
!	  epsy(na)       = dissipation, epsilon, times (\tau S)^2  (cm**2/sec***3)			
!     lifupper    = LOGICAL variable true in upper ocean where can use l_MLD S^2 dimensionalization.
!	  lifepsy(na) = LOGICAL array for epsy dimensionaliation points (i.e. lifupper .TRUE. and epsy>=0).


!
!	   Mixed Layer Depth definition
!	   idefmld    = 0: MLD = 0.1 degrees centigrade off surface pot. temp. .
!	         	  = 1: MLD = deltemld deg. centigrade off surface pot. temp. .
!   	          = 2: MLD = delrhmld g cm^{-3} off surface pot. density .
!      deltemld   = pot. temperature difference criterion for idefmld = 1
!	   delrhmld   = pot. density difference criterion     for idefmld = 2

!	   Switch for the 1 pt. closure model to allow use of the same set of model
!	   constants which have been used in the atmosphere.
!	   \gamma_{1 to 8} = [.96,.06,.16,.10,7.66,0.4,0,0.29] .
!	   ifchengcon	= 0: old ocean constants,near-surface profile assumption
!	   ifchengcon 	= 1: constants used in atmosphere,matched to experiments


!	   Switch for changing background diffusivities in old model with heat and
!		salt diffusivities equal
!		ifback=0   : Use input background values
!		ifback=4   : Background diffusivities from our model with S = N/sqrt(Ri)
!	       		     and fixed Ri = ri_internal
!		ifback=5   : Background diffusivities from our model with S = N/sqrt(Ri)
!			    	 and fixed Ri = backfrac*Ri_Critical


!		SWITCH TO TURN ON SALINITY MODEL
!		ifsali	= 0: old model with heat and salt diffusivities equal
!		ifsali	= 1: Canuto's salinity model as worked out by the summer of '98


!		Minimum on shear to avoid singularity in ifzeroshear = .FALSE case.
!		ifshearmin = .FALSE.:leave out minimum (to allow my older runs to be reproduced).
!		ifshearmin = .TRUE. :minimum foreground shear of s2min to avoid shear=>0 problems.



!	SWITCH TO ENABLE ZERO SHEAR FORMULA (for Ri more negative than table)
!	ifzeroshear = .FALSE. : Old method: always use table of Ri and Ri_d.
!	ifzeroshear = .TRUE. : New method: use table of N_d^2/N^2 for Ri < Ri_tableminimum


!	Switch for making salinity and heat backgrounds in ratio of S_S/S_H
!	using S_H and S_S from lowest level where they are nonzero for
!	levels where they are zero.(Leave them equal if S's zero at 1st level.)
!	ifsalback = 0: Use input salinity background values.
!	ifsalback = 1: Alter salinity background values to fit S_S/S_H ratio.
!	ifsalback = 2: As above but ratio taken at subcrit Ri at point's R_r.
!	ifsalback = 3: ALL backgrounds from our model using internal wave est.S.
!	ifsalback = 4: ALL bgs our model int. wave S=N/(Ri_i^(1/2)),Ri_i const..
!	ifsalback = 5: Like 4 but ra_r_i = constant factor * ra_r_crit.(theta_r)
!	ifsalback = 6: Like 4 but (S_M/(Sl/q))(ra_r_i)=constant*(S_M/(Sl/q))(0) 



!	Constants used in ifsalback=3,4 cases.
!	back_l_0 = estimated lengthscale, l_0, for int. wave-generated turb.(cm)
!	back_k_0 = `minimum turbulent' wavenumber which yields back_l_0  (cm^-1)


!	Constants used in ifsalback=3 case.
!	back_s2  = estimated shear squared due to internal waves (sec^{-2})
!	back_sm2 = 1/back_s2 (sec^2)
!	v_back0 = residual background constant viscosity(may need for stability)
!	t_back0 = residual background constant heat diffusivity 
!	s_back0 = residual background constant salt diffusivity 



!	Richardson numbers for the background turbulence in the ifsalback=3 case
!	back_ri1  = N^2/{S_{internal}}^2 = ({N_T}^2+{N_C}^2)/{S_{internal}^2
!	back_rid1 = {N_d}^2/(S_{internal})^2 = ({N_T}^2-{N_C}^2)/{S_{internal}^2
!	{S_internal}^2 \equiv back_s2



!	Ri and S^2 for the background turbulence in the ifsalback>=4 cases
!	ri_internal = Dubovikov model constant internal wave Richardson number 
!  	back_rid1 = {N_d}^2/(S_{internal})^2 = ({N_T}^2-{N_C}^2)/{S_{internal}^2
!	N^2 / {S_internal}^2  = ri_internal ; {S_internal}^2 = N^2 / ri_internal
!	back_rid1 = ri_internal {{N_T}^2 - {N_C}^2  \over N^2}, dividing by S^2,
!	back_rid1 = ri_internal {Ri_T - Ri_C \over Ri}
!	back_rid1 = ri_internal {Ri_d \over Ri}

!	(N^2 / S_internal^2) = ri_internal ;  S_internal^2 = (N^2 / ri_internal) . 
!	s2_back \equiv S_internal^2 (sec^(-2)).
!	s2_back   = N^2 / ri_internal
!	ifsalback > 4 case :
!	back_ra_r = array for background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!	Add array for c_y at ra_r_max. [See NBp.000319-3,4 (Vol.IX)].
!	c_y_r0 = array for c_y at maximum ({Ri_T^2 + {Ri_C}^2)^{1/2} at \theta_r
!	Add arrays for dimensionless turbulence functions at background ra_r
!	sm_r1 = array for S_M at background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!	sh_r1 = array for S_H at background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!	ss_r1 = array for S_S at background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r


!	Add an array for (Sl/q)^2 at background ra_r also.
!	slq2_r1 = array for (Sl/q)^2 at background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!	Introduce  a switch for interpolation of 1D arrays versus theta_r 
!	for background dimensionless functions instead of 2D arrays vs. Ri,Ri_d.
!	ifbg_theta_interp=0: Interpolate 2D array with (Ri,Ri_d) indices.
!	ifbg_theta_interp=1: Interpolate 1D array with theta_r index.




!	ifback = 5 :
!	backfrac = constant fraction of Ri_critical used for background
!	ifsalback =5 :
!	backfrac = constant fraction of ra_r_crit.(\theta_r) used for background
!	ifsalback =6 :
!	backfact = constant fraction of K_M(ra_r=0) which K_M background equals.


!	For case with ifdeeplat=1 interior value of (epson/N^2) is multiplied by
!	a latitude and stratification dependent factor equal to unity when
!	N =N_0= 5.24E-3 sec^{-1} and latitude = 30^o as per Gregg et al. Nature Vol.423.
!	epson2__ is the constant =pelagic val. for ifdeeplat=0;30^o,N_0 val. for ifdeeplat=1. 
!	epson2_ is the pelagic value = (\epsilon/N^2) for ifbotenhance=0.
!	epson2 is the value of \epsilon/N^2 used. 

!	For case with both ifepson2>0 and ifbotenhance>0 background epsilon/(N^2) increased 
!	from interior constant to higher values near the bottom.
!   epson2_    = non-enhanced backgound value for epson2 [epsilon/N^2)] (cm**2/sec).
!   eps_bot    = enhanced bottom dissipation of turbulent kinetic energy [epsilon] (cm**2 sec**(-3)).
!   epson2_bot = enhanced bottom epson2 [epsilon/N^2] (cm**2/sec). 



!	For case with "ifepson2=2", background AND foreground epsilon/(N^2)
!	taken as constant at levels below the highest in which foreground dies.
!	For case with "ifepson2=1", background epsilon/(N^2) taken as constant.
!	epson2     = (dissipation of turbulent kinetic energy)/(N^2), (cm**2/sec)


!     internal subroutines and functions:
!     formld,formld_te,formld_rh
!     for nmodel=1: ccoeff, ourl2 
!     for nmodel=2: mcoeff, mikel2, mksmsh, fct, rtwi

!     oursal2,mikesal2a,interp2d
!-----------------------------------------------------------------------
! For variables declares
!-----------------------------------------------------------------------
      	real(r8), parameter 	:: e=2.71828182845904509D0  	!yzf 201709111100
      	real(r8), parameter 	:: pi=3.14159265358979312D0  	!yzf 201709111100
      	integer, 	parameter 	:: nmodel=1        				!yzf 201709111100
      	integer, 	parameter 	:: ntbl= 501         			!yzf 201709111100


!	   	Temperature=Salt diffusivity model background model swith. 
      	integer, parameter :: ifback=5   !yzf 201709111100


!		Salinity model switch
      	integer, parameter ::  ifsali=1  			!yzf 201709111100

      	integer, parameter :: icondear=0  			!yzf 201709111100


!		Background (epsilon/N^2) dimensionalization of diffusivities switch.
      	integer, parameter :: ifepson2=2 			!yzf 201709111100



! 	      Value of (epsilon/N^2)/(1 cm/sec^2) used. See NBp.000203-2, Vol.VIII .
!           Introduce parameters for bottom enhancement.
!   	      ifbotenhance = 0 : no bottom enhancemant
!              		 = 1 : epsilon exponentially decreasing with height above bottom to min.


      	integer, parameter :: ifbotenhance=0   		!yzf 201709111100



!		Parameters for a Coriolis-based latitude dependence of \epsilon/N^2.
!		ifdeeplat = 0 : no Coriolis parameter dependence        
!				  = 1 : depends on Coriolis parameter and stratification as in Gregg et al.
      	integer, parameter :: ifdeeplat=1   		!yzf 201709111100
      	integer, parameter :: ifcheckbottomeps=0   	!yzf 201709111100



!		Switch for limiting BackGround ra_r to at most Foreground ra_r when Ri>0
! 		for R_r in the [R_r_crit_DoubleDiffusion,R_r_crit_SaltFingers] regime.
      	integer, parameter :: ifrafgmax=1    		!yzf 201709111100


!		Salinity background modification switch. 
      	integer, parameter :: ifsalback=5  			!yzf 201709111100


!		Extend table to avoid use of analytic limiting behavior for sal model.
!		Introduce option to have the nonlinear part of the table increase
! 		exponentially with the absolute value of the table index.
      	integer, parameter :: nextrtbl0= 62   		!yzf 201709111100


! 		ifexpabstable = {0,1} for {don't,do} use ~e^|i| nonlinear table spacing.
      	integer, parameter :: ifexpabstable=1   	!yzf 201709111100


!		Introduce option ifast=1, yielding ifastexpabs=1, to allow use of alternate interpolation scheme
!   	      "interp2d_expabs" tailored to exponential absolute value case, should be faster than "interp2d".
      	integer, parameter  :: ifast=1       		!yzf 201709111100
      	integer, parameter  :: ifastexpabs=ifast*ifexpabstable   	!yzf 201709111100


!     	Increase table size for e^|i| spacing to achieve sufficiently high values.
      	integer, parameter :: nextrtbl1=1000
      	integer, parameter :: nextrtbl=nextrtbl0+ifexpabstable*nextrtbl1    	!yzf 201709111100

!		mt0 sets the bounds for the part of the table with constant stepsize.
      	integer, parameter :: nposapprox=101       !yzf 201709111100
      	integer, parameter :: mt0=ntbl-nposapprox  !yzf 201709111100
      	integer, parameter :: mt=mt0+nextrtbl      !yzf 201709111100


      	integer, parameter :: ifchengcon=0          !yzf 201709111100  
      	integer, parameter :: idefmld=0             !yzf 201709111100


!		Introduce an integer parameter for the effect of rotation on the lengthscale.**
!		Introduce the complex variable zlomega for \sqrt{-B*/f^3} for diagnosis.
      	integer, parameter :: ilomega=0    			!yzf 201709111100



!		Introduce array for a deep lengthscale used in ifepson2=2 case.
! 		Do not want to make aldeep a parameter, but want it to be 
!           a big enough array to accomodate the model's number of levels.
      	integer, parameter :: nbig=100          	!yzf 201709111400



!		K_S/K_H (= "sisamax") as a function of angle in Ri_C,Ri_T space 
! 		at a radius just before the realizability limit.
      	integer, parameter :: mt_ra_r=nposapprox-1                    !yzf 201709111400


!		sisamax ranges from theta_r of -pi/4 to pi/4 to more than 
! 		cover the unrealizable region.
      	integer, parameter :: n_theta_r_oct=INT(((pi/4.D0)*mt_ra_r)/15.D0)*15   	!yzf 201709111400


!		Switch to write out polar 2D turbulence table. 
      	integer, parameter :: ifpolartablewrite=0   	!yzf 201709111400

!           Introduce flag for use of \theta_r arrays to interpolate background.
      	integer, parameter :: ifbg_theta_interp=1  		!yzf 201709111400

!        integer, parameter :: ifoutput = 0 


!		Minimum of linear range in lookup table.
      	real(r8), parameter :: ri0=- 4.D0   			!yzf 201709111100


      	real(r8), parameter :: s2min=1.D-14   			!yzf 201709111100
      	real(r8), parameter :: epson2__=.288D0    		!yzf 201709111100


!		The value of epsilon at the bottom. St. Laurent et al. JPO2001 give 
!       for slopes epsilon=(3--9)E-9 W/kg and decay scale = (150 +or- 50) meters
!  		and for crests and canyons (2--5)E-9 W/kg and decay scale = (500 +or- 100) meters .
!       1W/kg = 10^7 erg /(kg s) = 10^7 erg /(10^3 g s) = 10^{7-3} erg/(g s)=10^{4} cm^2/(s^3) .
      	real(r8), parameter :: eps_bot0=2.D-5    		!yzf 201709111100


!		The scale of decrease of bottom mixing with heigh in centimeters.
      	real(r8), parameter :: scale_bot=5.D4   		!yzf 201709111100


!		Gregg et al. admit their formula (equation(2)) breaks down right at the equator where
!		it would predict zero epsilon. Figure 1 of Gregg et al. suggests to me that the 
!		value at the equator of (\epsilon/\epsilon_reference) is between 0.02 and 0.05 .
!		In the v_3x3 NCOM the nearest tracer point to the equator is at 0.91^o and at 
!		N_0=5.24e-3 sec^{-1} at this latitude their L(\theta,f) equals \approx 0.0538.
!		Introduce a minimum on L(\theta,f_Coriolis), called eplatidepend in the program.
!		Current guess is best value for eplatidependmin between 0.1 and 0.05 .
      	real(r8), parameter :: eplatidependmin=7.D-2  	!yzf 201709111100


!	  	For the mixed layer depth criterion
      	real(r8), parameter :: deltemld=0.5D0       		!yzf 201709111100
      	real(r8), parameter :: delrhmld=0.125D-3    		!yzf 201709111100



!		amldminlom is the minimum MLD for use of the lengthscale \sqrt{-B*/f^3}.
      	real(r8), parameter :: amldminlom=5.D4   			!yzf 201709111400


!		Parameters for ifsalback=3 case.
!     	Gargett et. al. JPO Vol.11 p.1258-71 gives for "the deep record",
! 		\phi_0 = 6 \times 10^{-5} s^{-2} cpm^{-1} . "cpm" is 'cycles per meter'.
! 		\phi_0 = 6 \times 10^{-5} s^{-2} (2 pi/100)^{-1} cm 
      	real(r8), parameter :: back_ph_0=(6.D-5)*(1.D2/(2.D0*pi))     	!yzf 201709111400



!		Gargett et. al. favor the value, k_0 = 0.1 cpm. 
! 		But k_0=0.05-0.2 cpm might be viable, see section 5 of their paper.
! 		Take k_0 = 0.1 cpm * adjust_gargett, where adjust_gargett is adjustable.
!   	Convert to radians per cm: k_0 = 0.1 (2pi/100cm) * adjust_gargett.
!		used for ifsalback=4 case also, but set adjust_gargett=1 for ifsalback=4 
      	real(r8), parameter :: adjust_gargett=1.D0                        			!yzf 201709111400
      	real(r8), parameter :: back_k_0=(0.1D0)*(2.D0)*pi*(1.D-2)*adjust_gargett  	!yzf 201709111400


!       Introduce the lengthscale \Delta_0 \equiv pi/k_0 .
! 		The units of \Delta_0 are centimeters, with k_0 in radians per cm.
      	real(r8), parameter :: back_del_0=pi/back_k_0   			!yzf 201709111400
      	real(r8), parameter :: back_s2=back_ph_0*back_k_0   		!yzf 201709111400
      	real(r8), parameter :: back_sm2=1.D0/back_s2      			!yzf 201709111400


! 		Residual constant background diffusivities to be added to model ones. 
      	real(r8), parameter :: v_back0=0.01D0         !yzf 201709111400
      	real(r8), parameter :: t_back0=0.01D0         !yzf 201709111400
      	real(r8), parameter :: s_back0=0.01D0         !yzf 201709111400


!		Parameter for ifsalback=4 case.
      	real(r8), parameter :: ri_internal=1.D0     	!yzf 201709111400

      	
!       Parameter for ifback or ifsalback=5 case.
        real(r8), parameter :: backfrac = 85.D-2    !yzf 201709111400

!       Parameter for ifsalback=6 case.
        real(r8), parameter :: backfact = e**(-1)  !yzf 201709111400

!       Minimum shear^2 and switch for its use in ifzeroshear=.FALSE. case.
!       LOGICAL ifshearmin
        logical, parameter :: ifshearmin=.TRUE.     !yzf 201709111100
      
!       Zero shear parameterization for stongly unstable case switch.
!       LOGICAL ifzeroshear
        logical, parameter :: ifzeroshear=.TRUE.   !yzf 201709111100

!       Need timescale ratios to calculate R_r_crit .
!       Common block with ratios of timescales [See NBp.030403-8to10.]
        real(r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot     !yzf 201709111100
        COMMON /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot

        complex(kind=8) zlomega
        integer :: ifirst
        data ifirst/0/

!       Table realizability limit array.
!       DIMENSION irimax(-mt:mt)
        integer, dimension(-mt:mt) :: irimax    !yzf 201709111400

!       Change the parameter's name from "epson2" to "epson2_" 
!       because epson2 can vary with height above the bottom when ifbotenhance>0.
!       Further change the parameter's name from "epson2_" to "epson2__"
!       because epson2_ can vary with latitude and stratification when ifdeeplat>0.
        real(r8)  ::  epson2,epson2_,epson2_bot    !yzf 201709111100
!       The bottom enhancement epsilon value at the given height.[See NBp020214-1(Vol.XIII).] 
        real(r8) :: eps_bot             !yzf 201709111100
!       eps_bot at the level beneath (which is the bottom level for k=n).
        real(r8) :: eps_bot__under  !yzf 201709111100

!       Take variables to be used in unstable case with zero shear approximation 
!       from my code for HYCOM in common_blocks_giss_fixed2 [See NBp.030401-2to3.].
!       Add variables for use in Ri => -infinity case.
!       Add shearless table as a function of the one variable (N_d^2/N^2) .
!       Ratios needed for the fast interpolation of the exponential absolute value tables.
        real(r8), dimension(-mt:mt) :: and2on2a1,amtaun2a1,sma1,sha1,ssa1  !yzf 201709111100
        real(r8) :: dand2on2, rri, rnd2on2                                 !yzf 201709111100

!       For safety's sake make internal arrays' dimension at least 1 . 
!       Introduce height array for Ri_d.
!       Add z value for lowest ocean level to accomodate bottom enhancement.
!       Introduce array for an2.
        real(r8), dimension(na)   ::  rh                         !yzf 201709111100
        real(r8), dimension(nmax) ::  v_back,t_back,s_back       !yzf 201709111100

!       (First Day of Autumn)X Introduce an array for [\epsilon y], that is, dissipation times {\tau Shear}^2 .
        real(r8), dimension(na)   :: epsy
        real(r8), dimension(nbig) :: aldeep      !yzf 201709111400

!       INTEGER ifextermld   ! MLD switch for outside calculation
        real(r8), dimension(ntbl) :: ria,slq2a,sma,sha

!       Arrays for 2D table used in temperature-salinity model.
        real(r8), dimension(-mt:mt) :: rib,ridb                     !yzf 201709111400
        real(r8), dimension(-mt:mt,-mt:mt) :: slq2b,smb,shb,ssb     !yzf 201709111400

!       Critical ra_r[({Ri_T}^2 + {Ri_C}^2)^(1/2)] array for ifsalback>4 .
!       c_y(ra_rmax) array for ifsalback>4 to use as input guess for background.
!       Background Ri constant for ifback >= 4
!       Background ra_r as a function of \theta_r array for ifsalback>4
!       Background S_M,S_H,S_S at ra_r as functions of \theta_r arrays for ifsalback>4
!       Add (Sl/q)^2 as function of \theta_r array.
        real(r8), dimension(-n_theta_r_oct:3*n_theta_r_oct) :: sisamax, ra_rmax, c_y_r0, back_ra_r, sm_r1, sh_r1, ss_r1, slq2_r1         !yzf 201709111400
!       Introduce logical variable lifupper for region where l_MLD^2 S dimensionalization makes sense
!       and logical array lifepsy for points in column where use epsilon y dimensionalization.
        logical :: lifupper                  !yzf 201709111100
        logical, dimension(na) :: lifepsy    !yzf 201709111100

        integer :: n, nb, k
        real(r8) :: visc_cbu_limit, diff_cbt_limit, ako, back_l_0, dri

        real (r8) :: c_y0, c_y00
        data c_y0   /0.0d0/
        data c_y00  /0.0d0/

        character (500) :: bid_string
        character (500) :: run_file
        character (500) :: error_file
        integer :: nytid

!         integer, parameter :: write_option = 0
        common /mytid_table/ run_file, error_file
        common /mytid_table2/ nytid


        nytid     =     mytid
        write(bid_string,"(I4.4)") mytid
        run_file        =   "giss_runtime_report_"//trim(adjustl(bid_string))//".txt"
        error_file      =   "giss_error_report_"//trim(adjustl(bid_string))//".txt"

!       REFERENCE NOTE: START OF EXECUTABLE PROGRAM.
!       Raise flag and stop if number of levels exceeds nbig.
    n   =   na    

    v_back  =   0.0D0
    t_back  =   0.0D0
    s_back  =   0.0D0

    do k = 1, na
        call density_mcdougall2003(rh(k), t(k), s(k)*1000+35, z(k)*0.01)
    end do
    rh  =   rh*0.001

    if(nmax.gt.nbig) then
        open(unit=10,position='Append',file=error_file)
			write(10,*) 	" "
			write(10,*) 	" "
			write(10,*) 	"****************************"
			write(10,*) 	"*problem in turb_2 routine.*"
			write(10,*) 	"number of model levels exceeds nbig."
			write(10,*) 	"nbig=",nbig,"	nmax=",nmax
			write(10,*) 	"if want to use this many levels"// &
		             	" increase nbig to be bigger than nmax."
			write(10,*) 	"program is stopping now."
        close(10)
        stop
    end if

    nb=max(0,n)
    if(ifirst.eq.0) then
!	Add nmodel to writeout. nmodel=1,2 for 1,2 pt. closure.
       	open(unit=11,position='Append',file=run_file)
			write(11,*) " "
			write(11,*) " "
			write(11,*) "************************************"
			write(11,*) "************************************"
			write(11,*) "turbulence calculated by turb_2gi1b 040217 version."
			write(11,*) "************************************"
			write(11,*) "nmodel=",nmodel
			write(11,*) "************************************"
			write(11,*) "************************************"
			write(11,*) " "
			write(11,*) " "
        close(11)
!    Take zero shear unstable case *foreground* calculation from my code for HYCOM
!	 in inigiss_fixed2.fs0 . [See NBp.030401-2to3.] 
!	 When the Richardson number is more negative than the most negative table value
!	 it is more accurate to use the zero shear approximation derived from Canuto's
!	 analytic pure convection formula [See NBp.030326-3to9.].
!    Amend to make 1D table vs. (N_d^2/N^2) for zero shear unstable case.
!    Choose (N_d^2/N^2) table values to be the same as Ri_d table values.

!   Maximum diffusivity, fricmx, and surface minimum, wndmix, set externally
        visc_cbu_limit=fricmx
        diff_cbt_limit=fricmx
!   *SET THE VALUE OF THE KOLMOGOROV CONSTANT, K_0, HERE NAMED "ako".*
		ako = 1.6D0
!   START OF SALINITY MODEL BACKGROUND LENGTHSCALE CALCULATION SECTION.
		if(ifsali.eq.1) then
!   Calculate constant lengthscale for the background for ifsalback=3,4,5
!  	\Delta_0 = {B_1 pi \over (3 Ko)^{3/2}} l_0
!	l_0 = {(3 Ko)^{3/2} \over B_1 pi} \Delta_0
!	"back_l_0" is the constant background l_0 in centimeters.
          	if(nmodel.eq.1) then  
                call oursal2_2a(b1,0.d0,0.d0,slq2b_00, &
                		smb_00,shb_00,ssb_00, &
                		c_y0,c_y00,0,0)
	  		else if(nmodel.eq.2) then
	    		call mikesal2a(b1,0.d0,0.d0,slq2b_00, &
                     			smb_00,shb_00,ssb_00, &
                     			c_y0,c_y00,0,0)
	  		end if
			back_l_0 = (((3.D0*ako)**(3.D0/2.D0))/(b1*pi))*back_del_0
			if(ifsalback.eq.3) then
                open(unit=11,position='Append',file=run_file)
					write(11,*) " "
					write(11,*) "************************************"
					write(11,*) "internal wave constants for background."
					write(11,*) "residual constant background diffusivities:"
					write(11,*) "k_m /(cm^2 sec^{-1})",v_back0
					write(11,*) "k_h /(cm^2 sec^{-1})",t_back0
					write(11,*) "k_s /(cm^2 sec^{-1})",s_back0
					write(11,*) "."
					write(11,*) " "
					write(11,*) "shear^2/(sec^{-2}) =",back_s2
					write(11,*) "lengthscale, del_0/(cm) =",back_del_0
					write(11,*) "lengthscale, l_0/(cm) =",back_l_0
					write(11,*) '"adjust_gargett="',adjust_gargett
					write(11,*) "************************************"
					write(11,*) " "
                close(11)
			else if(ifsalback.ge.4) then
                open(unit=11,position='Append',file=run_file)
					write(11,*) " "
					write(11,*) "************************************"
					write(11,*) "dubovikov internal wave constants for background."
                close(11)
          		if(ifsalback.eq.4) then	
                    open(unit=11,position='Append',file=run_file)
	    				write(11,*) "internal wave richardson number=",ri_internal
                    close(11)
	  			else if(ifsalback.eq.5) then
                    open(unit=11,position='Append',file=run_file)
                        write(11,*) "ratio of background to critical ra_r"// &
                                    " [\\equiv ({ri_t}^2 + {ri_c}^2)^(1/2)]",backfrac
                    close(11)
	  			else if(ifsalback.eq.6) then
                    open(unit=11,position='Append',file=run_file)
                        write(11,*) "ratio of background dimensionless k_m ="// &
                                    " s_m/(s l/q) to its value at ri_t=ri_c=0", &
                                    backfact
                    close(11)
	  			end if
                open(unit=11,position='Append',file=run_file)
                    write(11,*) "lengthscale, del_0/(cm) =",back_del_0
                    write(11,*) "lengthscale, l_0/(cm) =",back_l_0
                    write(11,*) "************************************"
                    write(11,*) " "
                close(11)
			end if

			go to 10
      	end if     ! if(ifsali.eq.1) 
!   end of salinity model background lengthscale calculation section.
      	

        if(nmodel.eq.1) then
!   choose which set of model constants to use for cheng model.
        	if(ifchengcon.eq.0) then
          		call ccoeff(b1,rimax)
        	else if(ifchengcon.eq.1) then
          		call ccoeff1(b1,rimax)
        	end if
      	else
          	call mcoeff(b1,rimax)
      	endif


!   Calculate constant lengthscale for the background for ifback > 2
!  	\Delta_0 = {B_1 pi \over (3 Ko)^{3/2}} l_0
!	l_0 = {(3 Ko)^{3/2} \over B_1 pi} \Delta_0
!	"back_l_0" is the constant background l_0 in centimeters.
      	back_l_0 = (((3.D0*ako)**(3.D0/2.D0))/(b1*pi))*back_del_0

!   Temperature=Salinity Diffusivity Models Writeouts for background
      	if(ifback.ge.4) then
            open(unit=11,position='Append',file=run_file)
			write(11,*) " "
			write(11,*) "************************************"
			write(11,*) "dubovikov internal wave constants for background."
        	if(ifback.eq.4) then	
      			write(11,*) "internal wave richardson number=",ri_internal
        	else if(ifback.eq.5) then
      			write(11,*) "ratio of background to critical ri = ",backfrac
        	end if
			write(11,*) "************************************"
			write(11,*) " "
            close(11)
      	end if

        open(unit=11,position='Append',file=run_file)
        	write(11,*) "************************************"
			write(11,*) "isurfuse=",isurfuse
        	write(11,*) "************************************"
        close(11)

!   building the look-up tables of slq2,sm and sh vs. ri
        dri=(rimax-ri0)/float(ntbl-1)
        do k=1,ntbl
            ria(k)=ri0+DFLOAT(k-1)*dri
            if(k.eq.ntbl) ria(k)=rimax

            if(nmodel.eq.1) then
                call ourl2(b1,ria(k),slq2a(k),sma(k),sha(k))
            else
                call mikel2(b1,ria(k),slq2a(k),sma(k),sha(k))
            end if
        end do
!       end of building look-up tables


        open(unit=11,position='Append',file=run_file)
            write(11,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(11,*) "nmodel=",nmodel," the table's ntbl=",ntbl
			write(11,*) "temperature=salinity diffusivity model"
			write(11,*) "rimax  =",rimax
			write(11,*) "ifback =",ifback
			write(11,*) " "
			write(11,*) "icondear=",icondear
			
            if(icondear.eq.-1) then
	  			write(11,*) "do *not* use deardorff lengthscale modification."
			else if(icondear.eq.0) then
	  			write(11,*) "ye cheng's old deardorff:"// &
     						"modify l and \tau n leaving s_x unmodified."
			else if(icondear.eq.1) then
	  			write(11,*) "ye cheng's new deardorff:"// &
     						"modify l but leave *both* \tau n and s_x unmodified."
			end if

			write(11,*) " "
			write(11,*) "ifepson2=",ifepson2
        	
            if(ifepson2.eq.2) then
        		write(11,*) "epsilon/(n^2) used even for strong mixing beneath weak mixing" 
        	end if
			
            write(11,*) "ifdeeplat=",ifdeeplat
			
            if(ifdeeplat.gt.0) then
          		write(11,*) &
         			"use latitude dependence of interior ocean value of \\epsilon/n^2"
	  			write(11,*) "eplatidependmin=",eplatidependmin
			end if

			if(ifepson2.gt.0) write(11,*) "epson2__=",epson2__
			write(11,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        close(11)
		go to 17
   10   continue
!	yzf_Align_declaration_201709111750

        if(ifsali.eq.1) dri = -ri0/dfloat(mt0)

!   BUILD SALINITY MODEL TABLES VS. "Ri = Ri_T + Ri_C" AND "Ri_d = Ri_T - Ri_C".
!   Use separate loops for calculation of independent table variables.**
        do iridsign=0,1
            iridstep=(-1)**iridsign 
            do irid= 0,mt*iridstep,iridstep
!   set ri_d table values. (see nbp59,63=p#a27,30.)
                if(abs(irid).le.mt0) then
                    ridb(irid) = dfloat(irid)*dri
                else
                    mt0s = mt0*iridstep
                    mtm1s = (mt0-1)*iridstep
!   1st day of spring introduction of exponential absolute val table option.
                    if(ifexpabstable.eq.0) then
                        idifs = (abs(irid)-mt0)*iridstep
                        ridb(irid) = ridb(mt0s)*((ridb(mt0s)/ &
                                        ridb(mtm1s))**(idifs**2))
                    else if(ifexpabstable.eq.1) then
                        idif = abs(irid)-mt0
                        ridb(irid) = ridb(mt0s)*((ridb(mt0s)/ &
                                        ridb(mtm1s))**(idif))
                    end if
                end if
            end do 
       end do



       do irisign=0,1
           iristep=(-1)**irisign 
           do iri= 0,mt*iristep,iristep
!   set ri table values. (see nbp59,63=p#a27,30.)
                if(abs(iri).le.mt0) then
                    rib(iri) = dfloat(iri)*dri
                else
                    mt0s = mt0*iristep
                    mtm1s = (mt0-1)*iristep
!   1st day of spring introduction of exponential absolute val table option.
                    if(ifexpabstable.eq.0) then
                        idifs = (abs(iri)-mt0)*iristep
                        rib(iri) = rib(mt0s)*((rib(mt0s)/ &
                                    rib(mtm1s))**(idifs**2))
                    else if(ifexpabstable.eq.1) then
                        idif = abs(iri)-mt0
                        rib(iri) = rib(mt0s)*((rib(mt0s)/ &
                                    rib(mtm1s))**(idif))
                    end if
                end if
            end do 
        end do


!   if using interp2d_expabs introduce ratio between adjacent richardson numbers in nonlinear part of table.***
!   rnd2on2 is the ratio between adjacent n_d^2/n^2 in the nonlinear part of the zero shear 1d table.
        if(ifastexpabs.eq.1) then
            rri = rib(mt0)/rib(mt0-1)
!   only calculate rnd2on2 when zero shear parameterization is enabled.
            if(ifzeroshear) rnd2on2 = rri   
        end if

        do iridsign=0,1
           iridstep=(-1)**iridsign 
            do irid= 0,mt*iridstep,iridstep
                do irisign=0,1
                    iristep=(-1)**irisign
                    do iri= 0,mt*iristep,iristep
                        if(nmodel.eq.1) then
    !   need to pass back the value of b_1 from oursal2 for use here. 
    !   call version of submodule oursal2 which has an option to use b1={\tau s}^{3/2} (0,0).
                            call oursal2_2a(b1,rib(iri),ridb(irid),slq2b(iri,irid), &
                                   smb(iri,irid),shb(iri,irid),ssb(iri,irid), &
                                   c_y0,c_y00,iri,irid)
                        else if(nmodel.eq.2) then
    !   need to pass back the value of b_1 from mikesal2 for use here. 
    !   replace mikesal2 by mikesal2a
                            call mikesal2a(b1,rib(iri),ridb(irid),slq2b(iri,irid), &
                                       smb(iri,irid),shb(iri,irid),ssb(iri,irid), &
                                       c_y0,c_y00,iri,irid)
                        end if
                        if(slq2b(iri,irid).lt.0) then
                            irimax(irid) = iri - 1 
                            go to 15
                        end if
                    end do
   15      continue
                end do
            end do  
        end do
!   yzf_Align_201709131020


!   only calculate the zero shear table when the zero shear option is enabled.
        if(ifzeroshear) then
!   calculation of table for zero shear approximation from my hycom inigiss_fixed2.fs0 .
!   make 1 dimensional table of turbulence functions of n_d^2/n^2
!   to be used for the unstable case with negligible shear.
!   n_d^2 \equiv n_heat^2 - n_salt^2 . n_d^2/n^2 ranges from - to + infinity.
!   n_d^2 is analogous to ri_d^2, so i try making its table values exactly the same.
!   oursal2_zeroshear assumes shear^2=0 and n^2 < 0.
            dand2on2 = dri
!   set n_d^2/n^2 table values.
!   calculate moving out from zero index first postive indices and then negative ones.[see nbp.030407-08.]
            do ind2on2sign=0,1
                ind2on2step=(-1)**ind2on2sign
                do ind2on2 = 0,mt*ind2on2step,ind2on2step
                    if(abs(ind2on2).le.mt0) then
                        and2on2a1(ind2on2) = dfloat(ind2on2)*dand2on2
                    else
                        mt0s  = sign(mt0,ind2on2)
                        mtm1s = sign(mt0-1,ind2on2)
    !   introduction of exponential absolute val table option.
                        if(ifexpabstable.eq.0) then
                            idifs = sign((abs(ind2on2)-mt0),ind2on2)
                            and2on2a1(ind2on2) = and2on2a1(mt0s)*((and2on2a1(mt0s)/ &
                                                    and2on2a1(mtm1s))**(idifs**2))
                        else if(ifexpabstable.eq.1) then
                            idif = abs(ind2on2)-mt0
                            and2on2a1(ind2on2) = and2on2a1(mt0s)*((and2on2a1(mt0s)/ &
                                                    and2on2a1(mtm1s))**(idif))
                        end if
                    end if
                end do
            end do


            do ind2on2 = -mt,mt
                call oursal2_zeroshear &
                        (and2on2a1(ind2on2),amtaun2a1(ind2on2) &
                        ,sma1(ind2on2),sha1(ind2on2),ssa1(ind2on2))
            end do
        end if

!   Add writes in salinity model case.
       
        open(unit=11,position='append',file=run_file)
            write(11,*) " "
            write(11,*) " "
            write(11,*) " "
            write(11,*) "************************************************"
            write(11,*) " "
            write(11,*) "new temperature-salinity model"
            write(11,*) " "
            write(11,*) "************************************"
            write(11,*) "isurfuse=",isurfuse
            write(11,*) "************************************"
            write(11,*) "ifsali=",ifsali
            write(11,*) "ifsalback=",ifsalback

            if(.not.ifzeroshear) then
                write(11,*) "ifshearmin=",ifshearmin
                if(ifshearmin) write(11,*) "s2min=",s2min
            end if
	
            write(11,*) "ifzeroshear=",ifzeroshear
            write(11,*) "ilomega=",ilomega
            write(11,*) "amldminlom=",amldminlom
            write(11,*) "icondear=",icondear
        
            if(icondear.eq.-1) then
                write(11,*) "do *not* use deardorff lengthscale modification."
            else if(icondear.eq.0) then
                write(11,*) "ye cheng's old deardorff:"// &
                            "modify l and \tau n leaving s_x unmodified."
            else if(icondear.eq.1) then
                write(11,*) "ye cheng's new deardorff:"// &
                            "modify l but leave *both* \tau n and s_x unmodified."
            end if

            write(11,*) "ifepson2=",ifepson2

            if(ifepson2.gt.0) then 
                if(ifepson2.eq.2) write(11,*) &
                    "epsilon/(n^2) used even for strong mixing beneath weak mixing" 
                write(11,*) "ifdeeplat=",ifdeeplat
                if(ifdeeplat.gt.0) then
                    write(11,*) &
                        "use latitude dependence of interior ocean value of \\epsilon/n^2"
                    write(11,*) "eplatidependmin=",eplatidependmin
                end if
                write(11,*) "epson2__=",epson2__
                write(11,*) "ifbotenhance=",ifbotenhance
                if(ifbotenhance.eq.1) then
                    write(11,*) "eps_bot0=",eps_bot0
                    write(11,*) "scale_bot=",scale_bot
                end if
            end if

            write(11,*)"ifrafgmax=",ifrafgmax
            write(11,*) " "
            write(11,*)"ifbg_theta_interp=",ifbg_theta_interp
            write(11,*) " "
            write(11,*) " "
            write(11,*) "***********************************************"
            write(11,*) "ifexpabstable=",ifexpabstable
            write(11,*) "ifastexpabs=",ifastexpabs
            write(11,*) "***********************************************"
            write(11,*) " "
            write(11,*) &
                        "    i      ", &
                        "    rib(i)      ","    ridb(i)     ", &
                        "irimax(i)  "
            do i= -mt,mt
                write(11,"(i8,2x,2e16.4,i8,2x)") i,rib(i),ridb(i),irimax(i)
            end do

            write(11,*) " "
            write(11,*) "irid       ri_d        ri(irimax)  " &
                    // "s_m        s_h        s_s        " &
                  // "s_m/s_h    s_s/s_h    "
            do irid= -mt,mt
               write(11,"(i8,2x,2e12.4,3f11.6,2f11.4)") irid,ridb(irid),rib(irimax(irid)), &
                   smb(irimax(irid),irid), &
                   shb(irimax(irid),irid), &
                   ssb(irimax(irid),irid), &
                   smb(irimax(irid),irid)/shb(irimax(irid),irid), &
                   ssb(irimax(irid),irid)/shb(irimax(irid),irid)
            end do
        close(11)

!   yzf_Align_201709131120



!   CALCULATE "R_r_Critical" USING CANUTO'S 000228 ANALYTIC FORMULA
!	FOR "R_rho_Critical". See NBp.000229-3 and 000316-4.
!	R_rho_Canuto \equiv -Ri_C/Ri_T \equiv -R_r .
!	In a sheet dated 000228 Canuto gave me:
!	"R_\rho^{cr} = {1 \over \Deta} [1 {+\over-} \sqrt{1 - \Delta^2}] 
!	\Delta \equiv {{\pi_2(1 + {15 \over 7} \pi_3)} \over
!	{\pi_3 - \pi_2 + (15 \over 14} \pi_3^2}} ".
!	Note that the + and - choices are reciprocals so this covers
!	both the Salt Fingering and Double Diffusive Critical R_\rho's.
!	From Ocean Turbulence III paper have: 
!	\pi_{1,2,3,4,5} = 
!	(\tau_pc,\tau_c\theta,\tau_c,\tau_p\theta,\tau_\theta)/\tau 
!	R_r_Crit = [-1 -/+ \sqrt{1 - \Delta^2}]/Delta
!	\Delta = {{{\tau_c\theta \over \tau} ( 1 + (15/7)*{\tau_c \over \tau})}
!	        \over {{\tau_c \over \tau} - {\tau_c\theta \over \tau} + 
!	                (15/14) {\tau_c \over \tau}^2}}
        deltanum = tctot*(1.d0 + ((15.d0/7.d0)*tcot))
        deltaden = tcot - tctot + ((15.d0/14.d0)*(tcot**2))
        delta = deltanum/deltaden
        rrcrn = (-1.d0 - sqrt(1.d0 - (delta**2)))/delta
        rrcrp = (-1.d0 + sqrt(1.d0 - (delta**2)))/delta
        theta_rcrn = atan(rrcrn)
        theta_rcrp = atan(rrcrp)
!   Make sure the right choice of arctan(R_r)=[\theta_r] is made.
!	Arctan covers the range (-pi/2,pi/2) while 
!	\theta_r_Crit must be in the range (-pi/4,3pi/4) (The range of Ri>0.)
        if(theta_rcrn.lt.-pi/4.d0) theta_rcrn = theta_rcrn + pi
        if(theta_rcrp.lt.-pi/4.d0) theta_rcrp = theta_rcrp + pi
        theta_rcrn_deg = theta_rcrn*(180.d0/pi)
        theta_rcrp_deg = theta_rcrp*(180.d0/pi)

        open(unit=11,position='append',file=run_file)
            write(11,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(11,*) "r_r_crit+ =",rrcrp
            write(11,*) "r_r_crit- =",rrcrn
            write(11,*) "\\theta_r_crit+ =",theta_rcrp
            write(11,*) "\\theta_r_crit- =",theta_rcrn
            write(11,*) "\\theta_r_crit+ in degrees =",theta_rcrp_deg
            write(11,*) "\\theta_r_crit- in degrees =",theta_rcrn_deg
            write(11,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        close(11)
!	Increments in radial and angular coordinates in (Ri_T,Ri_C) plane.
    	delra_r = 1.D0/DFLOAT(mt_ra_r)
    	deltheta_r = (pi/4.D0)/DFLOAT(n_theta_r_oct)
!	Calculate the ratio \sigma_sa_max \equiv S_S/S_H as a function 
!	of the angle \theta_r in Ri_T,Ri_C space,
!       \theta_r \equiv arctan(Ri_C/Ri_T). 
!   The range of angles where unrealizability occurs is 
!	a subset of theta_r = -pi/4 to 3pi/4.


!   For Ri_T and Ri_C positive find the realizability limits  
!	in polar coordinates in the (Ri_T,Ri_C) plane : (ra_r,theta_r).
        
!         if(ifpolartablewrite.eq.1) then
!             if (mytid.eq.0) then
!                 open(unit=68,file="turb_ra_th",status="new")
!             end if
!         end if
	

        do itheta_r = -n_theta_r_oct,3*n_theta_r_oct
            theta_r = dfloat(itheta_r)*deltheta_r
            theta_r_deg = theta_r*(180.d0/pi)
!   introduce jtheta_r, an angle index that begins at zero   
!   for the purposes of letting oursal2 know it starts at the origin.
            jtheta_r = itheta_r + n_theta_r_oct
!   initialize sisamax to the impossible negative value of -99.999 to 
!   let places where the realizability limit is not reached stand out.
            sisamax(itheta_r) = -99.999
!   initialize sm_r0,sh_r0,ss_r0 to the inconsistent absurd value -9.999999.
            sm_r0 = -9.999999
            sh_r0 = -9.999999
            ss_r0 = -9.999999
!   flag ibg determines if the background value of ra_r has been calculated.
!   option had never worked. george halliwell points out i erroneously had "isalback" here before.
            if(ifsalback.eq.6) ibg=0
!   flag ifunreal determines if realizability limit has been found.
            ifunreal=0
!   write to file "back_ra_r_values", the allowed values of background ra_r.
            if(itheta_r.eq.-n_theta_r_oct) then
                if (mytid.eq.0) then
                    open(unit=69,file="back_ra_r_values")
                endif
            endif

!   Make the ra_r max value not too large to try to avoid numerical trouble.
            do ira_r = 0,(mt_ra_r**2)/4
                if(ira_r.le.mt_ra_r) then
                    ra_r = dfloat(ira_r)*delra_r
                else
                    ra_r = ((1.d0+delra_r)**(ira_r - mt_ra_r)) &
                            *(dfloat(mt_ra_r)*delra_r)
                end if

                rit = ra_r*cos(theta_r)
                ric = ra_r*sin(theta_r)
                ri_r  = rit + ric
                rid_r = rit - ric
                if(nmodel.eq.1) then
                    call oursal2_2a(b1,ri_r,rid_r,slq2_r, &
                           sm_r,sh_r,ss_r, &
                           c_y0,c_y00,ira_r,jtheta_r)
                else if(nmodel.eq.2) then
                    call mikesal2a(b1,ri_r,rid_r,slq2_r, &
                               sm_r,sh_r,ss_r, &
                               c_y0,c_y00,ira_r,jtheta_r)
                end if



                if(ifsalback.eq.6) then
                    smosloq_r = sm_r/sqrt(slq2_r)
                    if(ira_r.eq.0) smosloq_0r = smosloq_r
                    if((smosloq_r.le.backfact*smosloq_0r).and.(ibg.eq.0))then
                        ra_r1                   = ra_r
                        rit1                    = rit
                        ric1                    = ric
                        ri_r1                   = ri_r
                        rid_r1                  = rid_r
                        slq2_r1(itheta_r)       = slq2_r
                        sm_r1(itheta_r)         = sm_r
                        sh_r1(itheta_r)         = sh_r
                        ss_r1(itheta_r)         = ss_r
                        ibg                     = 1
                    end if
                end if
                if(slq2_r.le.0.d0) then 
                   sisamax(itheta_r) = ss_r0/sh_r0 
                   ra_rmax(itheta_r) = ra_r0
                   if(ifsalback.eq.5) then
                       back_ra_r(itheta_r) = backfrac*ra_rmax(itheta_r)
                   else if(ifsalback.eq.6) then
                       back_ra_r(itheta_r) = ra_r1
                   end if
                   ifunreal = 1 
                   go to 16
                end if
                ra_r0   = ra_r
                rit0    = rit
                ric0    = ric
                ri_r0   = ri_r
                rid_r0  = rid_r
                slq2_r0 = slq2_r
                sm_r0   = sm_r
                sh_r0   = sh_r
                ss_r0   = ss_r

!   Store c_y as c_y_0 for possible use as a  guess in background calc. .
                c_y_r0(itheta_r) = c_y0
            end do
!   close file with background ra_r values.
            if (mytid.eq.0) then
               if(itheta_r.eq.-n_theta_r_oct) close(69)
            endif
!   write out stability functions, the s's and sisamax.
   16  continue

!   set background ra_r large at angles where unrealizability doesn't occur.
!   make the ra_r max value not too large to try to avoid numerical trouble.
            if(ifunreal.eq.0) then
               ipenra_r = (mt_ra_r**2)/4-1
               back_ra_r(itheta_r) = ((1.d0+delra_r)**(ipenra_r - mt_ra_r)) &
                                        *(dfloat(mt_ra_r)*delra_r) 
            end if
!   for ifsalback=5 case get value for initialization of c_y calculation. 
            if(ifsalback.eq.5) then
                if(jtheta_r.eq.0) then
                    c_y001 = c_y0
                end if
            end if
        end do

!   yzf_Align_201709131450




!         if (mytid.eq.0) then
!             if(ifpolartablewrite.eq.1) close(68)
!         endif

        if(ifsalback.gt.4) then
            do itheta_r = -n_theta_r_oct,3*n_theta_r_oct
                theta_r = dfloat(itheta_r)*deltheta_r
                theta_r_deg = theta_r*(180.d0/pi)
                rit1 = back_ra_r(itheta_r)*cos(theta_r)
                ric1 = back_ra_r(itheta_r)*sin(theta_r)
                ri_r1  = rit1 + ric1
                rid_r1 = rit1 - ric1
                if(ifsalback.eq.5) then
                    jtheta_r = itheta_r + n_theta_r_oct
                    if(nmodel.eq.1) then
                        call oursal2_2a(b1,ri_r1,rid_r1,slq2_r1(itheta_r), &
                                sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
                                c_y_r0(itheta_r),c_y001,jtheta_r,1)
                    else if(nmodel.eq.2) then
                        call mikesal2a(b1,ri_r1,rid_r1,slq2_r1(itheta_r), &
                         sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
                                c_y_r0(itheta_r),c_y001,jtheta_r,1)
                    end if
                end if

                sisa1 = ss_r1(itheta_r)/sh_r1(itheta_r)

                if(slq2_r1(itheta_r).lt.0.d0) then
                    open(unit=10,position='Append',file=error_file)
                        write(10,*) &
                            "negative (sl/q)^2 in table of background vs. \\theta_r."
                        write(10,*) "itheta_r=",itheta_r, &
                             "   slq2_r1(itheta_r)=",slq2_r1(itheta_r)
                        write(10,*) "program is stopping in turb_2."
                    close(10)
                    stop
                end if
            end do
        end if

        if(ifzeroshear) then 
!                 close(68)
            open(unit=11,position='append',file=run_file)
            write(11,*) "************************************************"
            write(11,*)          "index (n_d^2)/(n^2)  -(\\tau n)^2    "// &
                 "s_m            s_h            s_s            "
            write(11,*) " "
            do ind2on2 = -mt,mt
                write(11,"(i6,5f15.9)") ind2on2,and2on2a1(ind2on2), &
                            amtaun2a1(ind2on2),sma1(ind2on2),sha1(ind2on2),ssa1(ind2on2)
            end do
            close(11)
            write(*,*) "************************************************"
        end if
   17   continue
        ifirst=1
    end if

!   yzf_Align_201709131520



!   REFERENCE NOTE: END OF INITIALIZATION.
!   Bottom level for use with bottom enhancement.
!	                               (Must be set outside initialization if block.)
    kbot=n+1
!   Surface Buoyancy Flux (can be used for dimensionalization of turbulence model) ** 
!   Total Buoyancy Flux = Sum of "turbulent" and "solar" contributions**
	buoytot = buoytur + buoysol
!   Choose the definition of MLD. 
!   If ifextermld=1, keep the one that was defined externally.
    if(ifextermld.eq.0) then
!   use mixed-layer routine only when there are at least two levels of sea.
        if(n.gt.0) then
            if(idefmld.eq.0) then
                call formld(z,t,amld,n)
            else if(idefmld.eq.1) then
                call formld_te(z,t,deltemld,amld,n)
            else if(idefmld.eq.2) then
                call formld_rh(z,rh,delrhmld,amld,n)
            end if
        else if(n.eq.0) then
            amld = z(1)
        else 
            amld = 0.d0
        end if
    end if

    al0=0.17*amld
!   Write internal turbulence quantities to fort.91 when writing enabled.
!   Headers for each outputstep.
!   Add S_M,H,S to outputs.
!   Add (epsilon (\tau S}^2) to outputs in the case where it is calculated.
!   START OF FIRST LOOP THROUGH LEVELS
    if(ifepson2.eq.2) then
!   initialize switch for sub(background-only) depth. 
        ifbelow=0      
    end if




    do 22 k=1,n
        if((ifepson2.eq.2).and.(ifdeeplat.gt.0)) then
            ifnofsmall=0
            if(an2(k).ge.0.d0) an = sqrt(an2(k))
            if((an/abs(coriol)).lt.1.d0) then   !arccosh(n/f) is undefined can't use gregg et al.
                ifnofsmall=1
            end if
        end if
        ri1=ri(k)
        if(ifsali.eq.0) then
            if(ri1.ge.rimax) then
                ri(k)=rimax    ! ad hoc
                if(ifback.eq.0) then
                    akm(k)=v_back(k)
                    akh(k)=t_back(k)
                    aks(k)=akh(k)
                    goto 22
                else
                    sm   = sma(ntbl)
                    sh   = sha(ntbl)
                    slq2 = slq2a(ntbl)  
                    ss   = sh
                end if
            elseif(ri1.le.ri0) then !asymptotically
                sm=sma(1)
                sh=sha(1)
                ss=sh
                slq2=slq2a(1)*ria(1)/ri1
            else ! linearly interpolate the look-up tables
                m=int((ri1-ri0)/dri)+1
                tmp=(ri1-ria(m))/(ria(m+1)-ria(m))
                sm=sma(m)+(sma(m+1)-sma(m))*tmp
                sh=sha(m)+(sha(m+1)-sha(m))*tmp
                slq2=slq2a(m)+(slq2a(m+1)-slq2a(m))*tmp
                ss=sh
            endif
        else if(ifsali.eq.1) then
            rid1=rid(k)
            and2 = (rid(k)/ri(k))*an2(k)
            and2on2 = and2/an2(k)
            if((ifzeroshear).and.(an2(k).lt.rib(-mt)*s2(k)))then
              ifpureshear=1
            else
              ifpureshear=0
            end if

            if(ifpureshear.eq.1) then
                imax = mt
                if(ifastexpabs.eq.0) then
                    call interp1d(and2on2, &
                          and2on2a1,amtaun2a1,sma1,sha1,ssa1, &
                          amtaun2,sm,sh,ss, &
                          imax,mt,mt0,dand2on2)
                else if(ifastexpabs.eq.1) then
                    call interp1d_expabs(and2on2, &
                          and2on2a1,amtaun2a1,sma1,sha1,ssa1, &
                          amtaun2,sm,sh,ss, &
                          imax,mt,mt0,dand2on2,rnd2on2)
                end if
                slq2 = (-amtaun2)/((b1**2)*ri(k))
                go to 5               !skip 2d interpolation.
            end if


            if(ifastexpabs.eq.0) then
                call interp2d(ri1,rid1, &
                      rib,ridb,slq2b,smb,shb,ssb, &
                      slq2,sm,sh,ss, &
                      irimax,mt,mt0,dri)
            else if(ifastexpabs.eq.1) then
                call interp2d_expabs(ri1,rid1, &
                      rib,ridb,slq2b,smb,shb,ssb, &
                      slq2,sm,sh,ss, &
                      irimax,mt,mt0,dri,rri)
            end if
        end if

        if(slq2.lt.0.d0) then
            open(unit=10,position='Append',file=error_file)
                write(10,*) "************************************************"
                write(10,*) "error detected in turbulence module." 
                write(10,*) "'slq2' negative in turb_2 subroutine" &
                            //" after interpolation."
                write(10,*) "k=",k,"     slq2=",slq2
                write(10,*) "sm=",sm,"   sh=",sh,"   ss=",ss
                write(10,*) "ri1=",ri1,"    rid1=",rid1
                write(10,*) "dri=",dri
                write(10,*) "program will stop."
                write(10,*) "************************************************"
            close(10)
            stop
        end if

!   Assume region contiguous with surface where foreground model is
!	realizable has ended when get 0 "slq2".
        IF(slq2.EQ.0.D0) ifbelow = 1
!*****C
!030401 Skipped from 1D table interpolation to here for unstable zero shear approximation.
    5   CONTINUE

!020912-24,25X ******Calculate epsy(k) \equiv epsilon * (\tau S)^2 for surface forcing dimensionalization.******
!X	       Take epsilon = -buoytot/(1 - ((1/2)(\tau S)^2)S_M))
!X	       epsy = -buoytot/((b1^2 (Sl/q)^2)^{-1} - (1/2)S_M) [See NBp020912-8,12 
!X	       AND NBp020925-2{extension on cover}.(error of (1/4) for (1/2) had to be corrected.) ]
!020919X Keep track of total number of points to compare with number of problem points.
!020923X Define problem points as points where epsy<0 *AND* l^2 S dimensionalization is used 
!X       in isurfuse=0 case AND l is based on MLD. 
!020923X Define negative problem points: problem points with Ri<0 ,
!X	 Keep track of number of points where l^2 S dimensionalization is used AND l is based on MLD.
!020924,26X,030504Z1 Introduce lifupper, logical which is 1 only where l_MLD^2 S would be used for isurfuse=0 .
!26x-030803zi1a      set lifupper false when unrealizable
        lifupper= (((ifepson2.lt.2).or.(ifbelow.eq.0)   &
                  .or.((ifdeeplat.gt.0).and.(ifnofsmall.eq.1)).or.& !030504 revert to l^2 s for n/f<1 in gregg et. al. case.  
                ((ri1.lt.0.d0).and.(k.le.2)))   &
                .and.(slq2.gt.0.d0)) 
        if(lifupper) then
            ilmldpoint=ilmldpoint+1
            if(ri1.lt.0.d0) ilmldpointneg=ilmldpointneg+1
        end if
        ipoint=ipoint+1
        if(k.eq.1) icall=icall+1
        epsy(k)=0.d0
        lifepsy(k)=.false.
        if((isurfuse.eq.1).and.(n.gt.0)) then       
            if(slq2.eq.0.d0) then
                epsy(k)=0.d0
            else
                epsy(k) = -buoytot/((1.d0/((b1**2)*(slq2))) - 0.5d0*sm)
            end if
            lifepsy(k)= ((epsy(k).ge.0.d0).and.lifupper)
            if((epsy(k).lt.0.d0).and.lifupper) then
                iproblem=iproblem+1
                if(ri1.lt.0.d0) inegproblem=inegproblem+1
            end if
        end if

!   yzf_align_201709131615

        akz=0.4*z(k)
        al=akz*al0/(al0+akz)
        if(ilomega.eq.1) then
            if((buoytot.lt.0.d0).and.(amld.ge.amldminlom)) then
                rlomega = sqrt((coriol**3)/(-buoytot))
            else
                rlomega = 0.d0
            end if
            rlblackadar = 1.d0/al  
            rl = rlblackadar + rlomega
            al = 1.d0/rl
        end if

        al2=al*al
        if(.not.(((ifepson2.eq.2).and.(ifbelow.eq.1)).or.lifepsy(k)))then
            if(ri1.gt.0.d0) then
                anlq2=slq2*ri1
                if((anlq2.gt.0.281d0).and.(icondear.ge.0)) then  !0.281=0.53**2
                    al2=0.281d0/anlq2*al2
                    if(icondear.eq.0) slq2=0.281d0/(ri1+1.d-20)
                endif
            endif
        end if





        if(an2(k).lt.0.d0) then
            epson2_ = epson2__      !just to "hold the place". value irrelevent here.
        else 
            if(ifdeeplat.eq.0) then
              epson2_ = epson2__
            else if((ifdeeplat.eq.1).or.(ifdeeplat.eq.2)) then !030803zi1a (see nbp.030803-4,5.)
                if(ifnofsmall.eq.1) then    
                    eplatidepend = 0.d0
                else
                    eplatidepend = eplatidepend_(abs(coriol),an)
                end if
                eplatidepend = max(eplatidepend,eplatidependmin)
                epson2_ = epson2__*eplatidepend
            end if
        end if

        if(ifepson2.ge.1) then
            if(ifbotenhance.eq.0) then 
                epson2 = epson2_ 
            else if(ifbotenhance.eq.1) then
                eps_bot = eps_bot0 * exp((z(k) - z(kbot))/scale_bot)
                epson2_bot = eps_bot/(ri(k)*s2(k))
                epson2 = max(epson2_,epson2_bot)
            end if
        end if


        if((ifback.ge.4).and.(ifsali.eq.0)) then
            if(ifback.eq.4) then
                back_ri1  = ri_internal
            else
                back_ri1 = backfrac*rimax    
            end if

            m=int((back_ri1-ri0)/dri)+1
            tmp=(back_ri1-ria(m))/(ria(m+1)-ria(m))
            sm_back=sma(m)+(sma(m+1)-sma(m))*tmp
            sh_back=sha(m)+(sha(m+1)-sha(m))*tmp
            slq2_back=slq2a(m)+(slq2a(m+1)-slq2a(m))*tmp
            ss_back=sh_back


            if(slq2_back.lt.0) then
                open(unit=10,position='append',file=error_file)
                    write(10,*) "************************************************"
                    write(10,*) "error detected in turbulence module."
                    write(10,*) "'slq2_back' negative in turb_2 subroutine" &
                                //" after 1d interpolation of background."
                    write(10,*) "k=",k,"     slq2_back=",slq2_back
                    write(10,*) "sm_back=",sm_back,"   sh_back=",sh_back
                    write(10,*) "back_ri1=",back_ri1
                    write(10,*) "dri=",dri
                    write(10,*) "program will stop."
                    write(10,*) "************************************************"
                close(10)
                stop
            end if
!yzf_align_201709131650



            s2_back = (ri1/back_ri1)*s2(k)
            if(ri1.le.0.d0) s2_back = 0.d0
            if(ri1.lt.0.d0) then
            sm_back = 0.d0
            sh_back = 0.d0
            ss_back = 0.d0
            end if
            if(ifepson2.eq.0) then
                al0_back = back_l_0
                akz=0.4d0*z(k)
                al_back=akz*al0_back/(al0_back+akz)
                al2_back=al_back*al_back
                if(back_ri1.gt.0.d0) then
                    anlq2_back=slq2_back*back_ri1
                    if((anlq2_back.gt.0.281d0).and.(icondear.ge.0)) then  !0.281=0.53**2
                        al2_back=0.281d0/anlq2_back*al2_back
                        if(icondear.eq.0) slq2_back=0.281d0/(back_ri1+1.d-20)
                    endif
                endif

                tmp_back=0.5d0*b1*al2_back*sqrt(s2_back/(slq2_back+1.d-40))
            else if(ifepson2.ge.1) then
                al_back=0.d0
                tmp_back=0.5d0*b1**2*back_ri1*slq2_back*epson2
            end if
            v_back(k)=tmp_back*sm_back
            t_back(k)=tmp_back*sh_back
            s_back(k)=tmp_back*ss_back
        end if



        if(ifsali.gt.0) then
            if(ifsalback.eq.1) then
                if(slq2.ne.0.d0) then
                    sm_last = sm
                    sh_last = sh
                    ss_last = ss
                else if(k.eq.1) then
                    s_back(k) = t_back(k)
                    go to 20 
                end if
                s_back(k) = (ss_last/sh_last)*t_back(k)
            else if(ifsalback.eq.2) then
                if(slq2.ne.0.d0) then
                    sisa = ss/sh
                else 
                    rit = (ri(k) + rid(k))/2.d0
                    ric = (ri(k) - rid(k))/2.d0
                    if(rit.eq.0.d0) then
                        if(ric.eq.0.d0) then
                            theta_r = atan(1.d0)
                        else
                            theta_r = pi/2.d0   ! arctangent of infinity.
                        end if
                    else
                        theta_r = atan(ric/rit)
                    end if


                    if(abs(theta_r).gt.(pi/2.d0)) stop
                    if(theta_r.lt.(-pi)/4.d0) theta_r = theta_r + pi
                    jtheta_r0 = int((theta_r + (pi/4.d0))/deltheta_r)
                    itheta_r0 = jtheta_r0 - n_theta_r_oct
                    itheta_r1 = itheta_r0+1
                    theta_r0 = itheta_r0*deltheta_r
                    theta_r1 = itheta_r1*deltheta_r
                    theta_r_deg = theta_r*180.d0/pi
	    
                    if((itheta_r1.gt.3*n_theta_r_oct).or.  &
                        (itheta_r0.lt.-n_theta_r_oct)) then
                        open(unit=10,position='Append',file=error_file)
                            write(10,*) "************************************************"
                            write(10,*) "problem in turbulence module!"
                            write(10,*) "unrealizability outside ri>0 region. "
                            write(10,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
                            write(10,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
                            write(10,*) "rit=",rit,"ric=",ric,"    theta_r=",theta_r
                            write(10,*) "theta_r_deg=",theta_r_deg
                            write(10,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
                            write(10,*) "n_theta_r_oct=",n_theta_r_oct 
                            write(10,*) "program will stop."
                            write(10,*) "************************************************"
                        close(10)
                        stop
                    end if

                    deltheta_r1 = theta_r - theta_r0
                    delsisa_r = sisamax(itheta_r1) - sisamax(itheta_r0)
                    dsisa_o_dtheta = delsisa_r/deltheta_r
                    sisa = sisamax(itheta_r0)+deltheta_r1*dsisa_o_dtheta
                end if
                s_back(k) = sisa*t_back(k)


            else if(ifsalback.eq.3) then
                back_ri1  = ri1*s2(k)*back_sm2
                back_rid1 = rid1*s2(k)*back_sm2
                if(ifastexpabs.eq.0) then
                    call interp2d(back_ri1,back_rid1, &
                          rib,ridb,slq2b,smb,shb,ssb, &
                          slq2_back,sm_back,sh_back,ss_back, &
                          irimax,mt,mt0,dri)
                else if(ifastexpabs.eq.1) then
                    call interp2d_expabs(back_ri1,back_rid1, &
                          rib,ridb,slq2b,smb,shb,ssb, &
                          slq2_back,sm_back,sh_back,ss_back, &
                          irimax,mt,mt0,dri,rri)
                end if


                if(slq2_back.lt.0) then
                    open(unit=10,position='append',file=error_file)
                        write(10,*) "************************************************"
                        write(10,*) "error detected in turbulence module."
                        write(10,*) "'slq2_back' negative in turb_2 subroutine" &
                                    //" after interpolation of background."
                        write(10,*) "k=",k,"     slq2_back=",slq2_back
                        write(10,*) & 
                            "sm_back=",sm_back,"   sh_back=",sh_back,"   ss_back=",ss_back
                        write(10,*) "back_ri1=",back_ri1,"   back_rid1=",back_rid1
                        write(10,*) "dri=",dri
                        write(10,*) "program will stop."
                        write(10,*) "************************************************"
                    close(10)
                    stop
                end if


                if(ifepson2.eq.0) then
                    al0_back = back_l_0
                    akz=0.4d0*z(k)
                    al_back=akz*al0_back/(al0_back+akz)
                    al2_back=al_back*al_back
                    if(back_ri1.gt.0.d0) then
                        anlq2_back=slq2_back*back_ri1
                        if((anlq2_back.gt.0.281d0).and.(icondear.ge.0)) then  !0.281=0.53**2
                            al2_back=0.281d0/anlq2_back*al2_back
                            if(icondear.eq.0) slq2_back=0.281d0/(back_ri1+1.d-20)
                        endif
                    endif

                    tmp_back=0.5d0*b1*al2_back*sqrt(back_s2/(slq2_back+1.d-40))

                else if(ifepson2.gt.0) then
                    tmp_back=0.5d0*b1**2*back_ri1*slq2_back*epson2
                end if
                v_back(k)=tmp_back*sm_back+v_back0
                t_back(k)=tmp_back*sh_back+t_back0
                s_back(k)=tmp_back*ss_back+s_back0

            else if(ifsalback.ge.4) then

                if(ifsalback.eq.4) then
                    back_ri1  = ri_internal
                    back_rid1 = (rid1/ri1)*ri_internal
                else
                    if(ri(k).le.0.d0) then
                        back_ra_r1 = 0.d0
                        back_rit1  = 0.d0
                        back_ric1  = 0.d0
                        back_ri1   = 0.d0
                        back_rid1  = 0.d0
                        go to 19
                    end if

                    rit = (ri(k) + rid(k))/2.D0
                    ric = (ri(k) - rid(k))/2.D0
                    ra_r = SQRT((rit**2) + (ric**2))

                    if(rit.eq.0.d0) then
                        if(ric.eq.0.d0) then
                            theta_r = atan(1.d0)
                        else
                            theta_r = pi/2.d0       ! arctangent of infinity.
                        end if
                    else
                        theta_r = atan(ric/rit)
                    end if

                    theta_r = ATAN(ric/rit)

                    if(abs(theta_r).gt.(pi/2.d0)) stop
                    if(theta_r.lt.(-pi)/4.d0) theta_r = theta_r + pi
                    jtheta_r0 = int((theta_r + (pi/4.d0))/deltheta_r)
                    jtheta_r1 = jtheta_r0+1
                    itheta_r0 = jtheta_r0 - n_theta_r_oct
                    itheta_r1 = itheta_r0+1
                    theta_r0 = itheta_r0*deltheta_r
                    theta_r1 = itheta_r1*deltheta_r
                    if((theta_r0.le.theta_rcrp).and.(theta_r.gt.theta_rcrp)) then
                        theta_r = theta_r1
                        theta_r0 = theta_r1
                        itheta_r0 = itheta_r1 
                        itheta_r1 = itheta_r1+1
                        theta_r1 = theta_r1 + deltheta_r
                    else if((theta_r1.ge.theta_rcrn).and.  &
                                    (theta_r.lt.theta_rcrn)) then
                        theta_r = theta_r0
                        theta_r1 = theta_r0
                        itheta_r1 = itheta_r0 
                        itheta_r0 = itheta_r0-1
                        theta_r0 = theta_r0 - deltheta_r
                    end if

                    theta_r_deg = theta_r*180.D0/pi

                    IF((itheta_r1.GT.3*n_theta_r_oct).OR.  &
                                (itheta_r0.LT.-n_theta_r_oct)) THEN
       
                        open(unit=10,position='append',file=error_file)
                            write(10,*) "************************************************"
                            write(10,*) "problem in turbulence module!"
                            write(10,*) "unrealizability outside ri>0 region. "
                            write(10,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
                            write(10,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
                            write(10,*) "rit=",rit,"ric=",ric,"    theta_r=",theta_r
                            write(10,*) "theta_r_deg =",theta_r_deg
                            write(10,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
                            write(10,*) "n_theta_r_oct=",n_theta_r_oct 
                            write(10,*) "program will stop."
                            write(10,*) "************************************************"
                        close(10)
                        stop
                    end if
                    deltheta_r1 = theta_r - theta_r0
                    delback_ra_r = back_ra_r(itheta_r1) - back_ra_r(itheta_r0)
                    dback_ra_r_o_dtheta = delback_ra_r/deltheta_r
                    back_ra_r1 = back_ra_r(itheta_r0) + & 
                              deltheta_r1*dback_ra_r_o_dtheta




                    ifrafglt=0
                    if(ifrafgmax.eq.1) then
                        if((theta_r.le.theta_rcrp).or.(theta_r.ge.theta_rcrn)) then
                            if(back_ra_r1.gt.ra_r) then
                                ifrafglt=1
                                back_ra_r1=ra_r
                            end if
                        end if
                    end if
   18    continue 
     
                    if(back_ra_r1.lt.0.d0) then
                        open(unit=10,position='append',file=error_file)
                            write(10,*) "************************************************"
                            write(10,*) "problem in turbulence module!"
                            write(10,*) "negative bg ra_r \\equiv (ri_t^2+ri_c^2)^(1/2)"
                            write(10,*) "back_ra_r1 =", back_ra_r1
                            write(10,*) "theta_r =", theta_r
                            write(10,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
                            write(10,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
                            write(10,*) "rit=",rit,"ric=",ric
                            write(10,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
                            write(10,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
                            write(10,*) "theta_r_deg =",theta_r_deg
                            write(10,*) "n_theta_r_oct=",n_theta_r_oct 
                            write(10,*) "program will stop."
                            write(10,*) "************************************************"

                        close(10)
                        stop
                    end if 

                    back_rit1 = cos(theta_r)*back_ra_r1
                    back_ric1 = sin(theta_r)*back_ra_r1
                    back_ri1  = back_rit1 + back_ric1
                    back_rid1 = back_rit1 - back_ric1
   
                end if

                if((ifbg_theta_interp.eq.0).or.(ifrafglt.eq.1)) then
                    if(ifastexpabs.eq.0) then
                        call interp2d(back_ri1,back_rid1, &
                              rib,ridb,slq2b,smb,shb,ssb, &
                              slq2_back,sm_back,sh_back,ss_back, &
                              irimax,mt,mt0,dri)
                    else if(ifastexpabs.eq.1) then
                        call interp2d_expabs(back_ri1,back_rid1, &
                              rib,ridb,slq2b,smb,shb,ssb, &
                              slq2_back,sm_back,sh_back,ss_back, &
                              irimax,mt,mt0,dri,rri)
                    end if
                else if(ifbg_theta_interp.eq.1) then
                    deltheta_r1 = theta_r - itheta_r0*deltheta_r
                    delsm_back = sm_r1(itheta_r1) - sm_r1(itheta_r0)
                    dsm_back_o_dtheta = delsm_back/deltheta_r
                    sm_back = sm_r1(itheta_r0) + & 
                              deltheta_r1*dsm_back_o_dtheta
                    delsh_back = sh_r1(itheta_r1) - sh_r1(itheta_r0)
                    dsh_back_o_dtheta = delsh_back/deltheta_r
                    sh_back = sh_r1(itheta_r0) + & 
                              deltheta_r1*dsh_back_o_dtheta
                    delss_back = ss_r1(itheta_r1) - ss_r1(itheta_r0)
                    dss_back_o_dtheta = delss_back/deltheta_r
                    ss_back = ss_r1(itheta_r0) + & 
                              deltheta_r1*dss_back_o_dtheta
                    delslq2_back = slq2_r1(itheta_r1) - slq2_r1(itheta_r0)
                    dslq2_back_o_dtheta = delslq2_back/deltheta_r
                    slq2_back = slq2_r1(itheta_r0) + &
                               deltheta_r1*dslq2_back_o_dtheta


                else
                    open(unit=10,position='append',file=error_file)
                        write(10,*) "problem with choice of background interpolation."
                        write(10,*) "ifbg_theta_interp=",ifbg_theta_interp
                        write(10,*) "ifrafglt=",ifrafglt
                        write(10,*) "program is stopping."
                    close(10)
                    stop
                end if



                if(slq2_back.lt.0) then
                    open(unit=10,position='append',file=error_file)
                        write(10,*) "************************************************"
                        write(10,*) "error detected in turbulence module."
                        write(10,*) "'slq2_back' negative in turb_2 subroutine" &
                                    //" after interpolation of background."
                        write(10,*) "k=",k,"     slq2_back=",slq2_back
                        write(10,*) &
                            "sm_back=",sm_back,"   sh_back=",sh_back,"   ss_back=",ss_back
                        write(10,*) "back_ri1=",back_ri1,"   back_rid1=",back_rid1
                        write(10,*) "dri=",dri
                        write(10,*) "program will stop."
                        write(10,*) "************************************************"
                    close(10)
                    stop
                end if


                s2_back = (ri1/back_ri1)*s2(k)
   19           if(ri1.le.0.d0) s2_back = 0.d0


                if(ri1.lt.0.d0) then
                    sm_back = 0.d0
                    sh_back = 0.d0
                    ss_back = 0.d0
                end if
!*****C
	  IF((sm_back.LT.0.D0).OR.  &
       (sh_back.LT.0.D0).OR.  &
       (ss_back.LT.0.D0)) THEN
       if (mytid.eq.0) then
	         WRITE(*,*) & 
        "************************************************"
	      WRITE(*,*) "Problem in turbulence module!"
	      WRITE(*,*) "Negative Structure Function in Background."
!        write(*,*) 'i=',ii,'j=',jj
        write(*,*) 'temperature'
        write(*,'(5d15.5)') t
        write(*,*)
        write(*,*) 'salinity'
        write(*,'(5d15.5)') s
        write(*,*)
        write(*,*) 'density'
        write(*,'(5d15.5)') rh
        write(*,*)
        write(*,*) 'richardson'
        write(*,'(5d15.5)') ri
        write(*,*)
        write(*,*) 'rid'
        write(*,'(5d15.5)') rid
        write(*,*)
        write(*,*) 'shear'
        write(*,'(5d15.5)') s2
        write(*,*)
        write(*,*) 'BV frequency'
        write(*,'(5d15.5)') an2
        write(*,*)
        write(*,*) 'ustart=',ustar_
        write(*,*)
        write(*,*) 'buoytur=',buoytur
        write(*,*)
        write(*,*) 'buoysol=',buoysol
	      WRITE(*,*) "slq2_back=",slq2_back
	      WRITE(*,*) "sm_back=",sm_back, &
                   " sh_back=",sh_back, &
                   " ss_back=",ss_back
	      WRITE(*,*) " "
	      WRITE(*,*) "back_ra_r1 =", back_ra_r1
	      WRITE(*,*) "theta_r =", theta_r
	      WRITE(*,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
	      WRITE(*,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
	      WRITE(*,*) " "
	      WRITE(*,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
	      WRITE(*,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
	      WRITE(*,*) "theta_r_deg=",theta_r_deg
	      WRITE(*,*) "n_theta_r_oct=",n_theta_r_oct 
	      WRITE(*,*) " "
	      WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	      WRITE(*,*) "rit=",rit,"ric=",ric
	      WRITE(*,*) " "
	      WRITE(*,*) " "
	      WRITE(*,*) "Program will stop."
	      WRITE(*,*) & 
        "************************************************"
         endif
	      STOP
	    END IF

               if(ifepson2.eq.0) then
                    al0_back = back_l_0
                    akz=0.4d0*z(k)
                    al_back=akz*al0_back/(al0_back+akz)
                    al2_back=al_back*al_back

                    if(back_ri1.gt.0.d0) then
                        anlq2_back=slq2_back*back_ri1
                        if((anlq2_back.gt.0.281d0).and.(icondear.ge.0)) then  !0.281=0.53**2
                            al2_back=0.281d0/anlq2_back*al2_back
                            if(icondear.eq.0) slq2_back=0.281d0/(back_ri1+1.d-20)
                        endif
                    endif
                    tmp_back=0.5d0*b1*al2_back*sqrt(s2_back/(slq2_back+1.d-40))
                else if(ifepson2.gt.0) then
                    tmp_back=0.5d0*b1**2*back_ri1*slq2_back*epson2
                end if
                v_back(k)=tmp_back*sm_back
                t_back(k)=tmp_back*sh_back
                s_back(k)=tmp_back*ss_back



	  IF((v_back(k).LT.0.D0).OR.  &
       (t_back(k).LT.0.D0).OR.  &
       (s_back(k).LT.0.D0)) THEN
       if (mytid.eq.0) then
	         WRITE(*,*) &
         "************************************************"
	      WRITE(*,*) "Problem in turbulence module!"
	      WRITE(*,*) "Negative Background Diffusivity."
!        write(*,*) 'i=',ii,'j=',jj
        write(*,*) 'temperature'
        write(*,'(5d15.5)') t
        write(*,*)
        write(*,*) 'salinity'
        write(*,'(5d15.5)') s
        write(*,*)
        write(*,*) 'density'
        write(*,'(5d15.5)') rh
        write(*,*)
        write(*,*) 'richardson'
        write(*,'(5d15.5)') ri
        write(*,*)
        write(*,*) 'rid'
        write(*,'(5d15.5)') rid
        write(*,*)
        write(*,*) 'shear'
        write(*,'(5d15.5)') s2
        write(*,*)
        write(*,*) 'BV frequency'
        write(*,'(5d15.5)') an2
        write(*,*)
        write(*,*) 'ustart=',ustar_
        write(*,*)
        write(*,*) 'buoytur=',buoytur
        write(*,*)
        write(*,*) 'buoysol=',buoysol
	      WRITE(*,*) "v_back=",v_back, &
                   " t_back=",t_back, &
                   " s_back=",s_back
	      WRITE(*,*) " "
	      WRITE(*,*) "slq2_back=",slq2_back
	      WRITE(*,*) "sm_back=",sm_back, &
                   " sh_back=",sh_back, &
                   " ss_back=",ss_back
	      WRITE(*,*) " "
	      WRITE(*,*) "back_ra_r1 =", back_ra_r1
	      WRITE(*,*) "theta_r =", theta_r, &
                   "   theta_r_deg=",theta_r_deg
	      WRITE(*,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
	      WRITE(*,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
	      WRITE(*,*) " "
	      WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	      WRITE(*,*) "rit=",rit,"ric=",ric
	      WRITE(*,*) " "
	      WRITE(*,*) " "
	      WRITE(*,*) "Program will stop."
	      WRITE(*,*) &
         "************************************************"
        endif
	      STOP
	    END IF
!000314 Stop if background diffusivities are zero at positive Ri.
!                   v_back(k)=0.01D0
!                   t_back(k)=0.01D0
!                   s_back(k)=0.01D0
	  IF((ri(k).GT.0.D0).AND.((v_back(k).EQ.0.D0).OR.  &
       (t_back(k).EQ.0.D0).OR.  &
       (s_back(k).EQ.0.D0))) THEN
       if (mytid.eq.0) then
	         WRITE(*,*) & 
        "************************************************"
	      WRITE(*,*) "Problem in turbulence module!"
	      WRITE(*,*) "Zero Background Diffusivity in stable case."
!        write(*,*) 'i=',ii,'j=',jj
        write(*,*) 'temperature'
        write(*,'(5d15.5)') t
        write(*,*)
        write(*,*) 'salinity'
        write(*,'(5d15.5)') s
        write(*,*)
        write(*,*) 'density'
        write(*,'(5d15.5)') rh
        write(*,*)
        write(*,*) 'richardson'
        write(*,'(5d15.5)') ri
        write(*,*)
        write(*,*) 'rid'
        write(*,'(5d15.5)') rid
        write(*,*)
        write(*,*) 'shear'
        write(*,'(5d15.5)') s2
        write(*,*)
        write(*,*) 'BV frequency'
        write(*,'(5d15.5)') an2
        write(*,*)
        write(*,*) 'ustart=',ustar_
        write(*,*)
        write(*,*) 'buoytur=',buoytur
        write(*,*)
        write(*,*) 'buoysol=',buoysol
	      WRITE(*,*) "v_back=",v_back(k), &
                   " t_back=",t_back(k), &
                   " s_back=",s_back(k)
	      WRITE(*,*) " "
	      WRITE(*,*) "slq2_back=",slq2_back
	      WRITE(*,*) "sm_back=",sm_back, &
                   " sh_back=",sh_back, &
                   " ss_back=",ss_back
	      WRITE(*,*) " "
	      WRITE(*,*) "slq2_r1(itheta_r0)=",slq2_r1(itheta_r0), &
                   " slq2_r1(itheta_r1)=",slq2_r1(itheta_r1)
	      WRITE(*,*) "sm_r1(itheta_r0)=",sm_r1(itheta_r0), &
                   " sm_r1(itheta_r1)=",sm_r1(itheta_r1)
	      WRITE(*,*) "sh_r1(itheta_r0)=",sh_r1(itheta_r0), &
                   " sh_r1(itheta_r1)=",sh_r1(itheta_r1)
	      WRITE(*,*) "ss_r1(itheta_r0)=",ss_r1(itheta_r0), &
                   " ss_r1(itheta_r1)=",ss_r1(itheta_r1)
	      WRITE(*,*) " "
	      WRITE(*,*) "back_ra_r1 =", back_ra_r1
	      WRITE(*,*) "theta_r =", theta_r
	      WRITE(*,*) "theta_r_deg =", theta_r_deg
	      WRITE(*,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
	      WRITE(*,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
	      WRITE(*,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
	      WRITE(*,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
	      WRITE(*,*) "n_theta_r_oct=",n_theta_r_oct 
	      WRITE(*,*) "deltheta_r=",deltheta_r
	      WRITE(*,*) " "
	      WRITE(*,*) " "
	      WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	      WRITE(*,*) "rit=",rit,"ric=",ric
	      WRITE(*,*) " "
	      WRITE(*,*) " "
	      WRITE(*,*) "Program will stop."
	      WRITE(*,*) &
         "************************************************"
          endif
	      STOP
	    END IF
!************************************************************************
	END IF
 20     CONTINUE
!*****C
        END IF



        aldeep(k)=0.d0

        if((ifepson2.eq.2).and.(ifbelow.eq.1)) then
            if(ri1.ge.0.d0.or.((ifdeeplat.eq.2).and.(ifnofsmall.eq.1))) then
                tmp=0.5d0*b1**2*ri1*slq2*epson2
            else if(k.gt.2) then
                if(k.eq.n) then
                    delz = z(k) - z(k-1)
                    delrh = rh(k) - rh(k-1)
                    del2rh = rh(k) - 2.d0*rh(k-1) + rh(k-2)
                else
                    delz = z(k+1) - z(k-1)
                    delrh = rh(k+1) - rh(k-1)
                    del2rh = rh(k+1) - 2.d0*rh(k) + rh(k-1)
                end if
                dzrh = delrh/delz
                d2zrh = 4.d0*del2rh/(delz**2)
                rdzlndzrh = dzrh/d2zrh
                al0deep=0.17d0*abs(rdzlndzrh)
                akz=0.4d0*z(k)
                aldeep(k)=akz*al0deep/(al0deep+akz)
                al2=aldeep(k)*aldeep(k)
                if(ifpureshear.eq.1) then 
                    go to 21
                else if(ifshearmin) then
                    s2(k) = max(s2(k),s2min)
                end if
                tmp=0.5d0*b1*al2*sqrt(s2(k)/(slq2+1.d-40))
            else
                if(ifpureshear.eq.1) then 
                    go to 21
                else if(ifshearmin) then
                    s2(k) = max(s2(k),s2min)
                end if
                if(lifepsy(k)) then
                    tmp=0.5d0*epsy(k)/(s2(k)+1.d-40)
                else 
                    tmp=0.5d0*b1*al2*sqrt(s2(k)/(slq2+1.d-40))
                end if
            end if
        else
            if(ifpureshear.eq.1) then 
                go to 21
            else if(ifshearmin) then
                s2(k) = max(s2(k),s2min)
            end if
            if(lifepsy(k)) then
                tmp=0.5d0*epsy(k)/(s2(k)+1.d-40)
            else 
                tmp=0.5d0*b1*al2*sqrt(s2(k)/(slq2+1.d-40))
            end if
        end if
   21   if(ifpureshear.eq.1) tmp=0.5d0*(b1**2)*al2*sqrt(-an2(k)/amtaun2)
            akm(k)=min(tmp*sm+v_back(k),visc_cbu_limit)
            akh(k)=min(tmp*sh+t_back(k),diff_cbt_limit)
            aks(k)=min(tmp*ss+s_back(k),diff_cbt_limit)
 22   continue
        do k=nb+1,nmax 
            akm(k)=0.d0
            akh(k)=0.d0
            aks(k)=0.d0
        end do
        if(n.gt.0) then
            if(akm(1).lt.wndmix) akm(1)=wndmix
            if(akh(1).lt.wndmix) akh(1)=wndmix
            if(aks(1).lt.wndmix) aks(1)=wndmix
        end if
        zlomega = sqrt(cmplx(-buoytot/(coriol**3)))


        DO k =1,n
           IF((akm(k).LT.0.D0).OR.(akh(k).LT.0.D0).OR.(akm(k).LT.0.D0)) THEN 
                open(unit=10,position='Append',file=error_file)
                    WRITE(10,*) "Diffusivity is negative."
                    WRITE(10,*) "k=",k
                    WRITE(10,*) "z[cm]      tem[C]     sal[ppt]   rho[g/cm3] "// &
                                 "Ri         Ri_d	   S^2[/s2]   "// &
                                 "K_M[cm2/s] K_H[cm2/s] K_S[cm2/s] "
                    WRITE(10,"(12(1pe11.3))") z(k),t(k),s(k),rh(k), &
                                    ri(k),rid(k),s2(k), &
                                    akm(k),akh(k),aks(k)
                    WRITE(10,*) "Program will stop."
                close(10)
                STOP
            END IF
        END DO


      return
      end subroutine nasa_giss_mixing_column
!-----------------------------------------------------------------------
!     end of subroutine nasa_giss_mixing_column 
!-----------------------------------------------------------------------



      function eplatidepend_(f,an) !hxs_modified
!030429-30 function for use in turb_2 to calculate a latitude, and brunt vaisala frequency,
!     dependent factor by which to multiply the constant we had been using for 
!     (\epsilon/ n^2) . use formula arccosh = ln(x+\sqrt{x^2 - 1}) in place of "acosh".
!     adapted from program gregglatifunc.f:
!030428 to calculate the latitude dependence of the paper of gregg,sanford & winkel
!     in nature vol.422, "reduced mixing from the breaking of internal waves 
!     in equatorial waters", equation (2).
!     f     = coriolis parameter, 2 \omega sin(latitude),   [sec^{-1}]
!     an    = brunt vaisala frequency, n,                   [sec^{-1}]
!     eplatidepend_ = l(\theta,n) from gregg et. al nature vol.422,3 april 2003, eqn.(2).
!     an0   = garerett and munk reference stratification,   [sec^{-1}]
!     f_30  = coriolis parameter at 30^o reference latitude,[sec^{-1}]
      
            use const_for_nasagissmixing, only: r8
            implicit none
            save

            real (r8),intent(in)    :: an,f

            real (r8), parameter    :: an0=5.24d-3 !see gregg et. al in methods section below eqn.(4).
            real (r8)                     :: f_30,eplatidepend_
            real (r8)                     :: anum,den
            real (r8)                     :: pi,omega
            integer                       :: icalled

            real (r8)                     :: acosh1,x
            real (r8)                     :: wavelat,xf,yn

            data icalled /0/
            save icalled,den
      
            acosh1(x) = log(x + sqrt((x**2) - 1.d0))
            wavelat(xf,yn) = xf*acosh1(yn/xf)

            pi     = 4.d0*atan(1.d0)
            omega  = pi/43082.0d0
            f_30=omega

      if(icalled.eq.0) then
            den=wavelat(f_30,an0)
            write(*,*) "f_30deg arccosh(n_0/f_0)  =",den
            icalled=icalled+1
      end if

      anum = wavelat(f,an)
      eplatidepend_ = anum/den

      end function eplatidepend_


	
!-----------------------------------------------------------------------
!     finds mixed layer depth
!-----------------------------------------------------------------------
!980501 Choice of definitions.
	subroutine formld(z,t,amld,n)

		use const_for_nasagissmixing, only: r8
		implicit none
        save
		integer,intent(in) :: n
		real(r8),dimension(n),intent(in) :: z,t
		real(r8),intent(out) ::amld

		integer :: k
		real(r8) :: tm

		!     z and amld are both positive

		do k=1,n
			if (abs(t(k) - t(1)).gt.0.1D0) then
				tm = t(1) - sign(0.1D0,t(1) - t(k))
				amld = z(k) + (z(k-1) - z(k))* &
					(tm - t(k))/(t(k-1) - t(k) + 1.e-20)
				return
			end if
		end do
		amld=z(n)
		return
 	end subroutine formld

!0501 Temperature difference criterion used, but chose outside \Delta T = delte.
	subroutine formld_te(z,t,delte,amld,n)
!980501 Make double precision to conform to calling cctmix routine.
		use const_for_nasagissmixing, only: r8
		implicit none
        save
		integer,intent(in) :: n
		real(r8),intent(in) :: delte
		real(r8),dimension(n),intent(in) :: z,t
		real(r8),intent(out) :: amld

		integer :: k
		real(r8) :: tm

!z and amld are both positive
		do k=1,n
			if (abs(t(k) - t(1)).gt.delte) then
		    	tm = t(1) - sign(delte,t(1) - t(k))
		    	amld = z(k) + (z(k-1) - z(k))* &
					(tm - t(k))/(t(k-1) - t(k) + 1.e-20)
		    	return
		  end if
		end do
		amld=z(n)
		return
	end subroutine formld_te

!0501 Density difference criterion used, but chose outside \Delta \rho = delrh.
	subroutine formld_rh(z,t,delrh,amld,n)
!980501 Make double precision to conform to calling cctmix routine.
		use const_for_nasagissmixing, only: r8
		implicit none
        save
		integer,intent(in) :: n
		real(r8),intent(in) :: delrh
		real(r8),dimension(n),intent(in) :: z,t
		real(r8),intent(out) :: amld

		integer :: k
		real(r8) :: tm

!     z and amld are both positive
		do k=1,n
			if (abs(t(k) - t(1)).gt.delrh) then
				tm = t(1) - sign(delrh,t(1) - t(k))
				amld = z(k) + (z(k-1) - z(k))* &
					(tm - t(k))/(t(k-1) - t(k) + 1.e-20)
				return
			end if
		end do
		amld=z(n)
		return
	end subroutine formld_rh
!******


!-----------------------------------------------------------------------
!     beginning of improved turbulence model subroutines (nmodel=1)
!-----------------------------------------------------------------------
	subroutine ccoeff(b1,rimax)
!980501 Make double precision to conform to calling cctmix routine.
		use const_for_nasagissmixing, only: r8
		implicit none	
        save
!slq=s*l/q, slq2=slq**2, 1/2*q**2=e
		real (r8), intent(out) :: b1,rimax 	!out parameters

		real (r8) :: g1,g2,g3,g4,g5,g6,g7,g8 	!(g1-g8)
		real (r8) :: d0,d1,d2,d3,d4,d5 			!(d0-d5)
		real (r8) :: s0,s1,s2,s3,s4,s5,s6,s7,s8,s9 	!(s0-s9)

		!temp variables
		real (r8) :: aa,bb,cc
		real (r8) :: temp1,temp2
		real (r8) :: prt0,qus,beta5,qty,cthem1,c7

		common/const0/g1,g2,g3,g4,g5,g6,g7,g8 	!(g1-g8)
		common/const1/d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9 	!(d0-d5) (s0-s9)

!----------------------------------------------------------------------
!     defines constants and model coefficients
!     b1 == (q/ustar)**3,  qus == q/ustar
!     g2=(1.-3.*v2dq2)/2.
!     g3=(3.*u2dq2-1.-g2)/3.
!     w2dq2=(1+g2-3.*g3)/3.
!     shih & shabbir: g2=0.00336, g3=0.0906
!     apbl: g2=0.06, g3=0.16
!     beta5=0.6*(1-c5), if c5=.3 then beta5=.42
!     qty == us**2*theta2/(w*theta)**2
!     cthem1 == 1/(c sub theta)
!-----------------------------------------------------------------------
		prt0=1.0D0
		b1=16.6d0
		qus=b1**(1./3.)
		g2 = 0.00336D0
		g3 = 0.0906D0
		g6 = 0.4D0
		g7 = 0.0D0
		beta5=0.42D0
		qty=3.1D0
		g1=4./3.*(3.*g3**2-g2**2)+4.*qus**(-4)
		g4=15./8.*beta5*g1
		temp1=1./6.*prt0*(1.+g2-3.*g3)*qus**4
		temp2=9.*(g6-g7)*(g6+g7+2.*prt0)/ &
			( (1.+g2-3.*g3)*prt0 )**2 *qus**(-4)
		g5=temp1*(1.+sqrt(1.+temp2))
		cthem1=1./(prt0*qus**2)* qty  ! 0.952766
		c7=1./3.
		g8 = (1.-c7)*cthem1
		!
		s0=3./2.*g1*g5*g5
		s1=-(g6+g7)*g4 + 2.*(g1-g3-1./3.*g2)*g4*g5 + 3./2.*g1*g5*g8
		s2=-3./8.*(g6*g6-g7*g7)*g1
		s3=3./2.*(g6 + g7)*g4 + (3.*g3 + g2)*g4*g5
		s4=2.*g5
		s5=2.*g4
		s6=-2./3.*g5*(g2*g2-3.*g3*g3)+1./2.*g1*g5*(g2-3.*g3)+3./4.*g1*(g6-g7)
		s7=3.*g5
		s8=3.*g4
		s9=g5*(3.*g3*g3-g2*g2)
		!
		d0=3.*g5*g5
		d1=g5*(7.*g4 + 3.*g8)
		d2=g5*g5*(3.*g3*g3 - g2*g2) - 3./4.*(g6*g6 - g7*g7)
		d3=g4*(4.*g4 + 3.*g8)
		d4=g4*(g2*g6 - 3.*g3*g7 - g5*(g2*g2 - g3*g3)) + g5*g8*(3.*g3*g3 - g2*g2)
		d5=1./4.*(g2*g2 - 3.*g3*g3)*(g6*g6 - g7*g7)
		!     find rimax:
		aa=(d3+g4)
		bb=(d4-s1/2.+s6/2.)
		cc=d5-s2/2.
		rimax=(-bb+sqrt(bb**2-4.*aa*cc))/(2*aa)
		rimax=rimax*0.999  !ad hoc
		return
	end subroutine ccoeff


!---------------------------------------------------------------------
!     start of subroutine ccoeff1
!---------------------------------------------------------------------

subroutine ccoeff1(b1,rimax)
      !13/9/17 
      !980912 Created to implement the option of using the same constants Ye Cheng
      !use precision_mod !cry_delete
      !use param_mod, only: mytid !cry_delete
      use const_for_nasagissmixing, only: r8 !cry_add 
      !implicit real(r8) (a-h,o-z)
      real (r8), intent(out) :: b1,rimax !cry check b1

      real(r8) :: g1,g2,g3,g4,g5,g6,g7,g8 !cry_add 
      common/const0/g1,g2,g3,g4,g5,g6,g7,g8

      real(r8) :: d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9 !cry_add 
      common/const1/d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9

      real(r8) :: prt0,qus !cry_add 
      real(r8) :: aa,bb,cc !cry_add 

!-----------------------------------------------------------------------
      prt0=1.0
      b1=16.6d0
      qus=b1**(1./3.)
!980912 Values used in atmosphere given me by Ye Cheng on 980911.
      g1=0.192D0/2.D0
      g2=0.06D0
      g3=0.16D0
      g4=0.1D0
      g5=7.66D0
      g6=0.4D0
      g7=0.D0
      g8=0.29D0
!
      s0=3./2.*g1*g5*g5
      s1=-(g6+g7)*g4 + 2.*(g1-g3-1./3.*g2)*g4*g5 + 3./2.*g1*g5*g8
      s2=-3./8.*(g6*g6-g7*g7)*g1
      s3=3./2.*(g6 + g7)*g4 + (3.*g3 + g2)*g4*g5
      s4=2.*g5
      s5=2.*g4
      s6=-2./3.*g5*(g2*g2-3.*g3*g3)+1./2.*g1*g5*(g2-3.*g3) +3./4.*g1*(g6-g7)
      s7=3.*g5
      s8=3.*g4
      s9=g5*(3.*g3*g3-g2*g2)

      d0=3.*g5*g5
      d1=g5*(7.*g4 + 3.*g8)
      d2=g5*g5*(3.*g3*g3 - g2*g2) - 3./4.*(g6*g6 - g7*g7)
      d3=g4*(4.*g4 + 3.*g8)
      d4=g4*(g2*g6 - 3.*g3*g7 - g5*(g2*g2 - g3*g3)) + g5*g8*(3.*g3*g3 - g2*g2)
      d5=1./4.*(g2*g2 - 3.*g3*g3)*(g6*g6 - g7*g7)
!     find rimax:
      aa=(d3+g4)
      bb=(d4-s1/2.+s6/2.)
      cc=d5-s2/2.
      rimax=(-bb+sqrt(bb**2-4.*aa*cc))/(2*aa)
      rimax=rimax*0.999  !ad hoc
      
      return

      end subroutine ccoeff1 
!---------------------------------------------------------------------
!     end of subroutine ccoeff1
!---------------------------------------------------------------------




	subroutine ourl2(b1,ri,slq2,sm,sh)
!980501 Make double precision to conform to calling cctmix routine.
		use const_for_nasagissmixing, only: r8
! 		use param_mod, only: mytid   !will be removed
		implicit none	
        save
!slq=s*l/q, slq2=slq**2, 1/2*q**2=e
		real (r8), intent(in) :: b1,ri 	!in parameters
		real (r8), intent(out) :: slq2,sm,sh 	!out parameters

		real (r8) :: g1,g2,g3,g4,g5,g6,g7,g8 	!(g1-g8)
		real (r8) :: d0,d1,d2,d3,d4,d5 			!(d0-d5)
		real (r8) :: s0,s1,s2,s3,s4,s5,s6,s7,s8,s9 	!(s0-s9)

		!temp variables
		real (r8) :: aa,bb,cc,dd
		real (r8) :: temp,q,y
		common/const0/g1,g2,g3,g4,g5,g6,g7,g8
		common/const1/d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9

        character (500) :: run_file
        character (500) :: error_file
        common /mytid_table/ run_file, error_file


		aa = (d3+g4)*ri*ri+(d4-s1/2.+s6/2.)*ri+d5-s2/2.
		bb = (d1+g5)*ri+d2-s0/2.
		cc = d0
		temp = bb**2-4.*aa*cc
		if(temp.le.0.) then
            open(unit=10,position='Append',file=error_file)
    			write(10,*) 'in turb_2, temp <= 0, stop, ri =',ri
            close(10)
			stop
		end if
		q = -1./2.*(bb+sign(1.D0,bb)*sqrt(temp))
		if(bb.lt.0.) then
			y=cc/q
		else
			y=q/aa
		endif
		if(y.lt.0.) then
            open(unit=10,position='Append',file=error_file)
				write(10,*) 'in turb_2, y < 0, stop; ri =', ri
            close(10)
			stop
		end if
		dd=d0+(d1*ri+d2)*y+(d3*ri*ri+d4*ri+d5)*y*y
		sm=(s0+(s1*ri+s2)*y)/dd
		sh=(s4+(s5*ri+s6)*y)/dd
		slq2=y/(b1*b1)
		return
	end subroutine ourl2

!---------------------------------------------------------------------
!     end of improved turbulence model subroutines (nmodel=1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     beginning of parameter-free turbulence model subroutines(nmodel=2)
!-----------------------------------------------------------------------


!---------------------------------------------------------------------
!     start of subroutine mcoeff
!---------------------------------------------------------------------

subroutine mcoeff(b1,rimax)
      !13/9/17
      !980501 Make double precision to conform to calling cctmix routine.
      !use param_mod, only: mytid !cry_delete
      use const_for_nasagissmixing, only: r8 !cry_add 
      !implicit real(r8) (a-h,o-z) !cry_delete
      real (r8), intent(out) :: b1,rimax !cry check b1

      b1=24.7D0
      rimax=1.68D0

      return
      
      end subroutine mcoeff
 
!---------------------------------------------------------------------
!     end of subroutine mcoeff
!---------------------------------------------------------------------

      



!---------------------------------------------------------------------
!     start of subroutine mikel2
!---------------------------------------------------------------------

    subroutine mikel2(b1,ri,slq2,sm,sh)
        use const_for_nasagissmixing, only: r8
!         use param_mod, only: mytid
        implicit none
        save

        real (r8),intent(in)  :: b1,ri
        real (r8),intent(out) :: slq2,sm,sh

        real (r8) :: x,ri1,y,val,yst,eps
        integer :: ier,iend,ifirst

        character (500) :: run_file
        character (500) :: error_file
        common /mytid_table/ run_file, error_file

        common /mikebb/ri1
        data ifirst/0/
        external fct

        ri1=ri
        eps=1.e-6                                              
        iend=10000
        y = 0.0d0                                             
                    
        if(ifirst.eq.0) then
            yst=0.7
            ifirst=1
        else
            yst=y                                
        endif

        call rtwi(y,val,fct,yst,eps,iend,ier)  
        if (ier.ne.0) then
            open(unit=10,position='Append',file=error_file)
                write(10,*) "rtwi call problem, ier=",ier
            close(10)
            stop
        endif

        x=ri*y
        call mksmsh(x,y,sm,sh)
        slq2=y/(b1*b1)
        return
    end subroutine mikel2

!---------------------------------------------------------------------
!     end of subroutine mikel2
!---------------------------------------------------------------------




!---------------------------------------------------------------------
!     start of function fct
!---------------------------------------------------------------------

    function fct(y)
        use const_for_nasagissmixing, only: r8
        implicit none
        save
        real (r8),intent(in)    :: y
        real (r8)               :: fct
        real (r8)               :: x,ri,sm,sh

        common /mikebb/ri
        x=ri*y
        call mksmsh(x,y,sm,sh)
        fct=2./(sm-ri*sh)
        return
    end function fct

!---------------------------------------------------------------------
!     end of function fct
!---------------------------------------------------------------------



                                          

!---------------------------------------------------------------------
!     start of subroutine mksmsh
!---------------------------------------------------------------------

      subroutine mksmsh(x,y,sm,sh)
      !980717 Make double precision to conform to calling cctmix routine.
            use const_for_nasagissmixing, only: r8
            implicit none
            save

            real (r8),intent(in)  :: x
            real (r8),intent(inout)  :: y
            real (r8),intent(out) :: sm,sh

            real (r8), parameter :: x1min=15.d0/8.d0
            real (r8) :: x1,fc,ri,a,eta,xi

            x1=0.5d0*sqrt(y)

            if(x1.lt.x1min) then
                  fc=5.d0/9.d0
            else
                  fc=1.d0-5.d0/(3.d0*x1)+25.d0/(16.d0*x1**2)
            end if

            y=(2.d0*x1)**2
            ri=x/y
            a=1.d0+1.d0/(1.+0.16d0*x1**2*ri)
            eta=0.05d0*x1**2*ri*a
            sm=1.d0/25.d0*fc/(1.d0+5.d0/(9.d0*fc)*eta)
            xi=x/60.d0
            sh=0.056d0/(1.d0+2.4d0*xi)
            return
      end subroutine mksmsh

!---------------------------------------------------------------------
!     end of subroutine mksmsh
!---------------------------------------------------------------------


    subroutine oursal2_2a(b1_arg,ri,rid,slq2,sm,sh,sc,c_y0,c_y00,iri,irid)
      
      !11/9/17
    !use precision_mod  !cry_delete
!     use param_mod, only: mytid 
    use const_for_nasagissmixing, only: r8 !cry_add
    
    implicit none !cry_add
    save

    integer, intent(in) :: iri,irid !cry_add
    real (r8), intent(in) :: ri,rid  !cry_add
    real (r8), intent(out) :: b1_arg,slq2,sm,sh,sc   !cry_add
    real (r8), intent(inout) :: c_y0,c_y00 !cry_add need check    
      real (r8), parameter :: b1_0=16.6d0,c_yst0=8.527882d0 !cry_modify
      logical, parameter :: ifchengb1=.false. !cry_modify
      
      !temp variables
      integer ib1set 
      data ib1set / 0 /
      
      integer ib1set0

      real (r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot !cry_add
    common /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot
    
    real (r8) :: rit,ric !cry_add
    common /bb/rit,ric

    character (500) :: run_file
    character (500) :: error_file
    common /mytid_table/ run_file, error_file


    external fct_sal

    real (r8) :: b1,taus2,taus3  !cry_add
    real (r8) :: val,c_yst  !cry_add
    real (r8) :: c_y,c_n,c_c  !cry_add
    real (r8) :: eps=1.d-6  !cry_modify                                          
    integer :: iend=300   !cry_modify 
    integer :: ier !cry_add                                       
  

      ib1set0=ib1set
      if(.not.ifchengb1) then
        if(ib1set0.eq.0) then
        b1 = b1_0
          b1_arg=b1
          ib1set=ib1set+1
        end if
      end if

    call smshsc_a3(0.d0,0.d0,0.d0,sm,sh,sc)
      
    if(iri.eq.0.and.irid.eq.0) then
            c_yst = c_yst0
    else if(iri.eq.0) then
            c_yst = c_y00
    else 
            c_yst = c_y0
    end if

      rit = (ri + rid)/2.d0
      ric = (ri - rid)/2.d0
    
    call rtwi(c_y,val,fct_sal,c_yst,eps,iend,ier) !need check
    
    if(ier.ne.0) then
        open(unit=10,position='Append',file=error_file)
            write(10,*) "in oursal2 subroutine"
            write(10,*) "c_y00=",c_y00,"   c_y0=",c_y0
            write(10,*) "ri=",ri,"   rid=",rid
            write(10,*) "rit=",rit," ric=",ric
            write(10,*) "initial guess for rtwi c_yst=",c_yst
            write(10,*) "rtwi call problem, ier=",ier
        close(10)
        stop
    endif

    if(c_y.ge.0) then
        c_y0=c_y
    else 
        if(ri.lt.0) then
            open(unit=10,position='Append',file=error_file)
                write(10,*) "c_y negative at negative ri"
                write(10,*) "ri=",ri,"   c_y=",c_y
                write(10,*) "unstable realizability limit unexpected:" 
                write(10,*) "stopping in oursal2."
            close(10)
            stop
        end if
    end if
       
      if(iri.eq.0) then 
            c_y00=c_y
      end if 

      if((iri.eq.0).and.(irid.eq.0).and.  &
                  (abs(c_y - c_yst0).gt.1.d-6)) then
            open(unit=10,position='Append',file=error_file)
                  write(10,*) "inconsistency in neutral value of c_y"
                  write(10,*) "value used =",c_yst0
                  write(10,*) "value calculated =",c_y
                  write(10,*) "program stopping in oursal2"
            close(10)
            stop
      end if


      if(ifchengb1) then
            if(ib1set.eq.0) then
                  taus2 = c_y/(tpvot**2)
                  taus3 = taus2*sqrt(taus2)
                  b1 = sqrt(taus3)
                  b1_arg = b1
                  ib1set=ib1set+1
            end if
      end if
      
      if(b1_arg.le.0.d0) then

            open(unit=10,position='Append',file=error_file)
                  write(10,*) "b1 <= zero."
                  write(10,*) "b1=",b1_arg
                  write(10,*) "something must be wrong."
                  write(10,*) "program is stopping in oursal2."
            close(10)
            stop
      end if
      
      if(ib1set.ne.1) then
            open(unit=10,position='Append',file=error_file)
                  write(10,*) "problem in oursal2; b1 not properly set."
                  write(10,*) "number of times b1 set=",ib1set
                  write(10,*) "b1_arg=",b1_arg
                  write(10,*) "program is stopping."
            close(10)
            stop
      end if
      
      if(ib1set0.eq.0) then
            open(unit=11,position='Append',file=run_file)

                  write(11,*) " "
                  write(11,*) "ifchengb1=",ifchengb1
                  write(11,*) "b_1=",b1_arg
                  write(11,*) " "
            close(11)
      end if
      
      ib1set0=ib1set

      slq2 = c_y/((b1*tpvot)**2)

    c_n = -(tcot*tctot/(tpvot**2))*c_y*rit
    c_c = -((tcot**2)/(tpvot**2))*c_y*ric
         
    call smshsc_a3(c_y,c_n,c_c,sm,sh,sc)


1003 format(12(i8))
1004 format(12(1pe14.5))
      
    end subroutine oursal2_2a



!---------------------------------------------------------------------
!     start of function fct_sal
!---------------------------------------------------------------------

    function fct_sal(c_y)                              
    
    !use precision_mod !cry_delete
    !use param_mod, only: mytid !cry_delete
    use const_for_nasagissmixing, only: r8 !cry_add

    implicit none !cry_add
    save
    
    real (r8),intent(in) :: c_y !cry_add
    real (r8) :: fct_sal !cry_add

    real (r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot !cry_add
    common /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot
    
    real (r8) :: rit,ric !cry_add
    common /bb/ rit,ric

    real (r8) :: c_n,c_c,sm,sh,sc !cry_add

    c_n = -((tcot*tctot)/(tpvot**2))*c_y*rit
    c_c = -((tcot**2)/(tpvot**2))*c_y*ric
    
    call smshsc_a3(c_y,c_n,c_c,sm,sh,sc) 
    
    fct_sal=(2.d0*(tpvot**2))/(sm-rit*sh-ric*sc)
    
    return                                          
    
    end function fct_sal

!---------------------------------------------------------------------
!     end of function fct_sal
!---------------------------------------------------------------------





                                        
                                                    

subroutine smshsc_a3(y,n,c,sm,sh,sc)
      !12/9/17  in: y,n,c   out: sm,sh,sc
      !use precision_mod  !cry_delelte
!       use param_mod, only: mytid
      use const_for_nasagissmixing, only: r8 !cry_add 
      !implicit real(r8) (a-h,o-z) !cry_delelte
      implicit none !cry_add
      save

      real(r8),intent(in) :: y,n,c !cry_add
      real(r8),intent(out) :: sm,sh,sc !cry_add

      real(r8) :: d,nm,nh,nc !cry_add
      real(r8) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p1m,p2m !cry_add
      real(r8) :: a0,a1,a2,a3,a4,a5 !cry_add
      real(r8) :: d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15 !cry_add

      integer :: ifrecall !cry_add

      real(r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot !cry_add
      common /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot

      integer :: mytid 
      common /mytid_table2/ mytid

      integer,parameter :: ifmodelconstout=0 !cry_add
      real(r8),parameter :: tpvot0 = 0.4d0, sgmt=0.72d0, tptot0=(1.d0/5.d0)*(1.d0/(1.d0+(1.d0/sgmt))),&
                                      tpcot0 = tptot0,ttot0=sgmt,tcot0=ttot0,tctot0=1.d0/3.d0 !cry_add
!****************************************
!000125 set common block passable timescale ratios to parameter statement values
!000127 include \tau_pv \over \tau .
      tpvot = tpvot0
      tptot = tptot0
      tpcot = tpcot0
      ttot  = ttot0
      tcot  = tcot0
      tctot = tctot0
!*****c
!000125 calculate the p's.
      p1  = 0.832d0
      p2  = 0.545d0
      p3  = (5.d0/2.d0)*tpcot
      p4  = (1.d0/5.d0)*tpcot*(tcot**(-2))
      p5  = tpcot*tctot*(tcot**(-2))
      p6  = (1.d0/5.d0)*(tcot**(-1))*(tctot**(-1))*tptot
      p7  = 5.d0*tctot
      p8  = (5.d0/2.d0)*tptot
      p9  = ttot*tptot*((tcot*tctot)**(-1))
      p10 = tctot*tptot*(tcot**(-2))
      p11 = tpcot*(tcot**(-1))
      p1m = 1.d0 - p1
      p2m = 1.d0 - p2
!*****c
!results.2_1
!990513 values of a's and d's calculated from p's using cheng's fortran code
!     to do so, from today's email from him, cheng990513.results.2_1 .
!results.2_1
!##########################
!##  fortran code:
!##########################
      a0 = 12
      a1 = p11*(12*p9+8*p6-30*p6*p8-5*p6*(p1m+3*p2m))
      a2 = 5*(2*p4*p6*p7-p4*p9-p6*p11)*(p1m+3*p2m)+8*p6*p11+8*p4*p9-16*p4&
                  *p6*p7+12*p11*p9+12*p11*p10-12*p4*p7**2*p6-30*p6*p11*p8+30*p4*p6*&
                  p7*p8+30*p6*p4*p7*p3-30*p4*p9*p3
      a3 = p10*(12*p11+8*p4-30*p3*p4-5*p4*(p1m+3*p2m))
      a4 = -p6*(8-30*p8-5*p1m-15*p2m)-12*p9-12*p11
      a5 = -p4*(8-30*p3-5*p1m-15*p2m)-12*p10-12*p11
      
      d0 = 24
      d1 = p11*((-p6-2*p9)*p1m**2+(p6+6*p9)*p2m**2+2*p6*p8*(p1m-3*p2m))
      d2 = (2*p4*p6*p7-p4*p9-p6*p11)*(p1m**2-p2m**2)+2*(-p11*p10-p11*p9+&
                  p4*p7**2*p6)*(p1m**2-3*p2m**2)+2*(-p6*p4*p7*p3-p4*p6*p7*p8+p4*p9*p3&
                  +p6*p11*p8)*(p1m-3*p2m)
      d3 = p10*((-p4-2*p11)*p1m**2+2*p4*p3*(p1m-3*p2m)+(6*p11+p4)*p2m**2)
      d4 = -4*p6*p11*(3*p9+2*p6)
      d5 = 4*p4*p6**2*p7*(4+3*p7)-4*p4*p9*(3*p11+2*p6)-4*p6*p11*(3*p9+3*&
                  p10+2*p4+2*p6)
      d6 = 4*p4**2*p6*p7*(4+3*p7)-4*p4*p9*(2*p4+3*p11)-8*p4*p6*(p11+p10)&
                  -12*p10*p11*(p4+p6)
      d7 = -4*p4*p10*(2*p4+3*p11)
      d8 = (2*p9+2*p11+p6)*p1m**2-2*p6*p8*(p1m-3*p2m)-(p6+6*p9+6*p11)*p2m**2
      d9 = (2*p10+p4+2*p11)*p1m**2-2*p4*p3*(p1m-3*p2m)-(p4+6*p10+6*p11)*p2m**2
      d10 = 8*p6**2+4*(7*p11+3*p9)*p6+24*p11*p9
      d11 = -8*(4+3*p7)*p4*p6*p7+4*p4*(4*p6+7*p9+3*p11)+4*p6*(3*p10+7*p11)&
                  +24*p11*(p10+p9)
      d12 = 4*p10*(7*p4+6*p11)+4*p4*(2*p4+3*p11)
      d13 = 6*p2m**2-2*p1m**2
      d14 = -28*p6-24*p9-24*p11
      d15 = -24*p10-28*p4-24*p11
!results.2_1
!****************************************

!980728     write out the p's.
!000125 writeout the timescale ratios as well.
      ifrecall=1

      if(ifrecall.eq.0) then
            if (mytid.eq.0) then
                  write(*,*) "tau_pv/tau     =",tpvot 
                  write(*,*) "tau_ptheta/tau =",tptot
                  write(*,*) "tau_pc/tau =",tpcot
                  write(*,*) "tau_theta/tau  =",ttot
                  write(*,*) "tau_c/tau  =",tcot
                  write(*,*) "tau_ctheta/tau  =",tctot
                  write(*,*) " "
                  write(*,*) "p1 =",p1
                  write(*,*) "p2 =",p2
                  write(*,*) "p3 =",p3
                  write(*,*) "p4 =",p4
                  write(*,*) "p5 =",p5
                  write(*,*) "p6 =",p6
                  write(*,*) "p7 =",p7
                  write(*,*) "p8 =",p8
                  write(*,*) "p9 =",p9
                  write(*,*) "p10=",p10
                  write(*,*) "p11=",p11
!990513 write out the a's and d's as well.
                  write(*,*) "a0=",a0
                  write(*,*) "a1=",a1
                  write(*,*) "a2=",a2
                  write(*,*) "a3=",a3
                  write(*,*) "a4=",a4
                  write(*,*) "a5=",a5
                  write(*,*) "d0=",d0
                  write(*,*) "d1=",d1
                  write(*,*) "d2=",d2
                  write(*,*) "d3=",d3
                  write(*,*) "d4=",d4
                  write(*,*) "d5=",d5
                  write(*,*) "d6=",d6
                  write(*,*) "d7=",d7
                  write(*,*) "d8=",d8
                  write(*,*) "d9=",d9
                  write(*,*) "d10=",d10
                  write(*,*) "d11=",d11
                  write(*,*) "d12=",d12
                  write(*,*) "d13=",d13
                  write(*,*) "d14=",d14
                  write(*,*) "d15=",d15
!990513 output p#, a# and d# to the file model_constants if the switch is set.
!000125 writeout the timescale ratios as well.
                  if(ifmodelconstout.eq.1) then
                        open(unit=66,file='model_constants',status='unknown')
                              write(*,*) "tau_pv/tau     =",tpvot 
                              write(*,*) "tau_ptheta/tau =",tptot
                              write(*,*) "tau_pc/tau =",tpcot
                              write(*,*) "tau_theta/tau  =",ttot
                              write(*,*) "tau_c/tau  =",tcot
                              write(*,*) "tau_ctheta/tau  =",tctot
                              write(*,*) " "
                              write(66,*) "p1 =",p1
                              write(66,*) "p2 =",p2
                              write(66,*) "p3 =",p3
                              write(66,*) "p4 =",p4
                              write(66,*) "p5 =",p5
                              write(66,*) "p6 =",p6
                              write(66,*) "p7 =",p7
                              write(66,*) "p8 =",p8
                              write(66,*) "p9 =",p9
                              write(66,*) "p10=",p10
                              write(66,*) "p11=",p11
                              write(66,*) "a0 =",a0
                              write(66,*) "a1 =",a1
                              write(66,*) "a2 =",a2
                              write(66,*) "a3 =",a3
                              write(66,*) "a4 =",a4
                              write(66,*) "a5 =",a5
                              write(66,*) "d0 =",d0
                              write(66,*) "d1 =",d1
                              write(66,*) "d2 =",d2
                              write(66,*) "d3 =",d3
                              write(66,*) "d4 =",d4
                              write(66,*) "d5 =",d5
                              write(66,*) "d6 =",d6
                              write(66,*) "d7 =",d7
                              write(66,*) "d8 =",d8
                              write(66,*) "d9 =",d9
                              write(66,*) "d10=",d10
                              write(66,*) "d11=",d11
                              write(66,*) "d12=",d12
                              write(66,*) "d13=",d13
                              write(66,*) "d14=",d14
                              write(66,*) "d15=",d15
                  close(66)
                  end if
!*****c
            endif
      end if
      
      ifrecall = 1
!******

!980610 modification of section of "sx" containing the den and nums of the "s"'s

!###############################################

      d = d0 + d1*y*n**2 + d2*y*n*c + d3*y*c**2 + d4*n**3 + d5*n**2*c &
            + d6*n*c**2 + d7*c**3 &
            + d8*y*n + d9*y*c + d10*n**2 + d11*n*c + d12*c**2 + d13*y &
            + d14*n + d15*c

!########################################################################
      nm = a0 + a1*n**2 + a2*n*c + a3*c**2 + a4*n + a5*c

!###########################################################################


      nh = - (30.d0*n*p6 + 30.d0*c*p4 - 60.d0 &
            - ( 2.d0*p1m  + 15.d0*p2m**2 - 6.d0*p2m  - 5.d0*p1m**2 ) &
            * y ) &
            * (c*p4*p7 - c*p11 - n*p11 + 1.d0)


      nc = (30.d0*n*p6 + 30.d0*c*p4 - 60.d0 &
            - ( 2.d0*p1m  + 15.d0*p2m**2 - 6.d0*p2m  - 5.d0*p1m**2 ) &
            * y ) &
            * (c*p10 - 1.d0 - n*p6*p7 + n*p9)


!980610-15 modification of section of "sx" containing sm, sh, sc

!*******************************************************************************
      sm = (4.d0/15.d0) * tpvot * nm/d

      sh = (4.d0/15.d0) * tptot * nh/d

      sc = (4.d0/15.d0) * tpcot * nc/d
!*******************************************************************************
      return

1004 format(12(1pe14.5))

      end subroutine smshsc_a3







subroutine rtwi(x,val,fct,xst,eps,iend,ier)
!     to solve general nonlinear equations of the form x=fct(x)
!     by means of wegsteins iteration method
!     prepare iteration
      use const_for_nasagissmixing, only: r8
      implicit none
      save

      real (r8),intent(in) :: xst,eps  !in parameters
      integer,intent(in) :: iend
      integer,intent(out) :: ier    !out parameters
      real (r8), intent(out) :: x   !out parameters
      real (r8), intent(inout) :: val

      !temp variables
      real (r8) :: tol,a,b,d
      real (r8) :: fct
      integer :: i

      ier=0
      tol=xst
      x=fct(tol)
      a=x-xst
      b=-a
      tol=x
      val=x-fct(tol)
!     start iteration loop
      do i=1,iend
!981103 crude fix to avoid mysterious problem which 
!occurred with a close but not too close guess.
            if(dabs(val).lt.1.d-12) then
                  val =0.d0
                  return
            else
                  b=b/val-1.d0                                      
                  if (b .eq. 0.0d0) then      !2,8,2
                        ier=2
                        return
                  else
                        a=a/b
                        x=x+a
                        b=val
                        tol=x
                        val=x-fct(tol)
                        tol=eps                                 
                        d=abs(x)
                        if(d .gt. 1.d0) tol=tol*d
                        if(abs(a) .le. tol) then !5,5,6
                              if(abs(val) .le. 10.d0*tol) then !7,7,6
                                    return
                              else
                                    continue
                              end if
                        else
                              continue
                        end if
                  end if
            end if
      end do
!     no convergence after iend iteration steps. error return. 
      ier=1
      return
end subroutine rtwi
                                                                        


!030401-04Y Adaptation to the NCAR CSM Ocean Model of my calculation of the turbulence
!Y	functions of the one variable (N_d^2/N^2) for the zero shear unstable case
!Y      written for HYCOM [See NBp.030401-2to3]. I must make variables double precision
!Y	and remove the inclusion of the HYCOM file common_blocks_giss.h and instate the
!Y	common block bb0 which carries the timescale ratios [See NBp.030402-2].
!Y	[See NBp.030403-1.] Introduce a check on positivity of S_X in response to an 
!Y	error that I found in the calculation[See NBp.030404-5,6,8 .]
!030403Y Include "ttot'=`{\tau_\theta \over \tau}' in the common block with timescale 
!Y       ratios, /bb0/. See NBp030403-8to12.
!
!_______________________________________________________________________
!-----------------------------------------------------------------------
      subroutine oursal2_zeroshear(and2on2,amtaun2,sm,sh,sc)
      !12/9/17
      !use precision_mod !cry_delelte
!       use param_mod, only: mytid
      use const_for_nasagissmixing, only: r8 !cry_add

      !implicit real(r8) (a-h,o-z) cry_delelte
      implicit none !cry_add
      save

      real (r8), intent(in) :: and2on2 !cry_add
      real (r8), intent(out) :: amtaun2,sm,sh,sc !cry_add

      real (r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot !cry_add
      common /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot

      real (r8) :: eps=1.e-6  !cry_modify                                          
      integer :: iend=300   !cry_modify 

      real (r8) :: anh2on2,ans2on2,taunh2,tauns2 !cry_add 
      real (r8) :: c_y,c_n,c_c  !cry_add

      integer :: ier !cry_add 

      integer :: mytid
      common /mytid_table2/ mytid

      call smshsc_a3(0.d0,0.d0,0.d0,sm,sh,sc)     

      anh2on2 = (1. + and2on2)/2.
      ans2on2 = (1. - and2on2)/2.

      call quad_sal_pureconv(anh2on2,ans2on2,amtaun2,ier)

      if(ier.ne.0) then
            if (mytid.eq.0) then
                  write(*,*) "in oursal2_zeroshear subroutine"
                  write(*,*) "anh2on2",anh2on2,"      ans2on2=",ans2on2
                  write(*,*) "error returned by quad_sal_pureconv ier=",ier
                  write(*,*) "ier=1 means choice of root must be reconsidered" &
                                    //"for these model constants."
            end if
            stop
      end if

      taunh2 = -amtaun2*anh2on2
      tauns2 = -amtaun2*ans2on2
      
      c_y = 0.
      c_n = -(tcot*tctot)*taunh2
      c_c = -(tcot**2)*tauns2
         
      call smshsc_a3(c_y,c_n,c_c,sm,sh,sc)

!030304y stop in case any of the s_x is negative.
      if((sm.lt.0.d0).or.(sh.lt.0.d0).or.(sc.lt.0.d0)) then !cry_170913_change_'ss'_to_'sc'
            if (mytid.eq.0) then
                  write(*,*) " "
                  write(*,*) "in subroutine oursal2_zeroshear after call to smshsc_a3"
                  write(*,*) "error: negative structure function."
                  write(*,*) "sm=",sm
                  write(*,*) "sh=",sh
                  write(*,*) "sc=",sc
                  write(*,*) " "
                  write(*,*) "c_y=",c_y
                  write(*,*) "c_n=",c_n
                  write(*,*) "c_c=",c_c
                  write(*,*) " "
                  write(*,*) "taunh2=",taunh2
                  write(*,*) "tauns2=",tauns2
                  write(*,*) "anh2on2=",anh2on2
                  write(*,*) "ans2on2=",ans2on2
                  write(*,*) "and2on2=",and2on2
                  write(*,*) " "
                  write(*,*) "program is stopping."
            endif
            stop
      end if


1003 format(12(i8))
1004 format(12(1pe14.5))
      
      end subroutine oursal2_zeroshear

!_______________________________________________________________________
!-----------------------------------------------------------------------
!030403Y    I *MUST* add "ttot"=`{\tau_\theta \over \tau}' to the "/bb0/" common block with 
!Y	    timescales, because I use it in this routine. Trouble with this routine made me
!Y	    realize that the common block was missing it. [See NBp.030407-12.]
!030401-02Y Adaptation of the hycom routine to the NCOM must make explicity double precision.
!Y	    Fix error of hyphen for minus sign in amtaun2_negri_0rid test and
!Y	    remove common_blocks_giss.h inclusion since I don't have that file in current NCOM
!Y	    and instate the common block bb0 with the timescale ratios [See NBp.030402-2.]


  subroutine quad_sal_pureconv(anh2on2,ans2on2,amtaun2,ier) 
  !12/9/17                              
  !use precision_mod !cry_delete
!   use param_mod, only: mytid
  use const_for_nasagissmixing, only: r8 !cry_add
  !implicit real(r8) (a-h,o-z) !cry_delete
  real (r8), intent(in) :: anh2on2,ans2on2 !cry_add
  real (r8), intent(out) :: amtaun2 !cry_modify
  integer, intent(out) :: ier !cry_modify

  real (r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot !cry_add
  common /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot

  real (r8), parameter :: errorbound=1.d-12,amtaun2_negri_0rid = (16.6**2)* &
    (-(6.5489662174209907d-04)*(-56.94668746578522)) !cry_modify

  real (r8) :: eta3p,amu3p
  real (r8) :: braasal,braatem,braa,para,a3p
  real (r8) :: brabsal,brabtem,brab,parb,b3p
  real (r8) :: anum,bnum,cnum,radical,qnum,rootplus,rootminus
  real (r8) :: errorplus,errorminus,errorproplus,errorprominus

integer :: mytid
common /mytid_table2/ mytid

  amtaun2 = 0.
  ier = 0


  eta3p = tpcot*(tctot*anh2on2 + tcot*ans2on2)
  amu3p = tptot*(tctot*ans2on2 + ttot*anh2on2)


  braasal = tpcot*(amu3p  - tptot*tctot*anh2on2)*ans2on2
  braatem = tptot*(eta3p - tpcot*tctot*ans2on2)*anh2on2
  braa = braasal+braatem
  para = (eta3p*amu3p - tptot*tpcot*(tctot**2)*anh2on2*ans2on2)
  a3p   = -(7./15.)*braa - para

  brabsal = tpcot*ans2on2
  brabtem = tptot*anh2on2
  brab = brabsal+brabtem
  parb = (eta3p + amu3p)
  b3p   = -(7./15.)*brab - parb


  anum = a3p
  bnum = -b3p
  cnum = -1.
  radical = sqrt(bnum**2 - 4*anum*cnum)
  qnum = (-1./2.)*(bnum + sign(radical,bnum))
  rootplus  = qnum/anum
  rootminus = cnum/qnum
  !030325-27 check that calculated solutions actually satisfy the quadratic equation.
  !030327 make acceptable error a small fraction of the root size.
  errorplus  = a3p*(rootplus**2)-b3p*(rootplus)-1.
  errorminus = a3p*(rootminus**2)-b3p*(rootminus)-1.
  errorproplus  = errorplus/rootplus
  errorprominus = errorminus/rootminus
  
  if(max(abs(errorproplus),abs(errorprominus)).gt.errorbound) then
    if (mytid.eq.0) then
      write(*,*) "problem: error too great in calculating amtaun2 ."
      write(*,*) "in quad_sal_pureconv:"
      write(*,*) "-(\\tau n)^2=q/a \\equiv first root:",rootplus
      write(*,*) "a''' (- \\tau n)^2)^2 - b'''((-\\tau n)^2) - 1 =", &
      errorplus
      write(*,*) "proportional error= ",errorproplus
      write(*,*) "-(\\tau n)^2=c/q \\equiv second root:",rootminus
      write(*,*) "a''' (- \\tau n)^2)^2 - b'''((-\\tau n)^2) - 1 =", &
      errorminus
      write(*,*) "proportional error= ",errorprominus
      write(*,*) "theoretically both roots should give zero."
      write(*,*) "maximum proportional error we permit =",errorbound
      write(*,*) "a=",anum," b=",bnum," c=",cnum
      write(*,*) "q=",qnum
      write(*,*) "program is stopping."
    endif
    stop
  end if
  !******
  !030326-0402y     i choose "rootminus" because it approximates the solution with shear
  !     at very negative ri and because it bends downward as (n_h^2/n_d^2)
  !     departs from zero which is consistent with an additional driving of
  !     mixing due to double-diffusivity even in the unstable case.
  !     if change model constants *must* revisit the choice of root.
  amtaun2 = rootminus
  
  if(abs(anh2on2-ans2on2).lt.errorbound) then
    if(abs(rootminus-amtaun2_negri_0rid).gt.1.) then  
      if (mytid.eq.0) then
        write(*,*) "check if have the right root."
        write(*,*) "rootminus=",rootminus
        write(*,*) "amtaun2_negri_0rid=",amtaun2_negri_0rid
      endif
      ier=1
    end if
  end if
  
  return 

  end subroutine quad_sal_pureconv                                            
                                                    



!_______________________________________________________________________







!---------------------------------------------------------------------
!     start of subroutine mikesal2a
!---------------------------------------------------------------------

!990701-02 submodule intended to calculate mikhail dubovikov's model with salt
!   for the ogcm turbulence module in analogy with oursal2 for cheng's. 
! introduce arguments for last guess and last row's zero guess for c_y.  
!************************************************************c

    subroutine mikesal2a(b1_arg,ri,rid,slq2,sm,sh,ss,c_y0,c_y00,index0,index1)
        use const_for_nasagissmixing, only: r8
!         use param_mod, only: mytid
        implicit none
        save

        real (r8),intent(in)  :: ri,rid
        real (r8),intent(out) :: b1_arg,slq2,sm,sh,ss
        real (r8),intent(inout) :: c_y0,c_y00
        integer,intent(in)  :: index0,index1

        real (r8), parameter :: tpvot = 0.4d0
        real (r8), parameter :: b1=24.7d0

        real (r8) :: ri1
        common /mikebb/ri1

        real (r8) :: rid1,ric1,rit1,rr1
        common /mikebbs/rid1,ric1,rit1,rr1

        integer :: ifirst
        data ifirst/0/

      integer :: mytid
      common /mytid_table2/ mytid

        external fctysal

        real (r8) :: eps,rty_00,c_yst0,yst,c_yst,c_y
        real (r8) :: y,val,fctysal
        real (r8) :: x,z
        integer   :: iend,ier

        b1_arg=b1
        ri1=ri
        rid1=rid
        rit1=0.5d0*(ri+rid)
        ric1=0.5d0*(ri-rid)
        rr1 = ric1/rit1
        
        eps=(1.e-6)*c_y0
        iend=50000

        if((ifirst.eq.0).or.(index0.eq.0.and.index1.eq.0)) then
            rty_00 = (5.d0/6.d0)*(2.d0 + sqrt(67.d0))
            yst=rty_00**2
            c_yst0 = yst*(tpvot**2)
            c_yst = c_yst0
            ifirst=1
        else if(index0.eq.0) then
            c_yst = c_y00
            yst=c_yst/(tpvot**2)                                
        else
            c_yst = c_y0
            yst=c_yst/(tpvot**2)                                
        end if

        call rtwi(y,val,fctysal,yst,eps,iend,ier)  

        if(ier.ne.0) then
            if (mytid.eq.0) then
                write(*,*) "rtwi call problem, ier=",ier
                write(*,*) "could not solve turbulence model at:"
                write(*,*) "first index=",index0,"second index=",index1
                write(*,*) "in mikesal2a subroutine"
                write(*,*) "c_y00=",c_y00," c_y0=",c_y0
                write(*,*) "ri=",ri,"       rid=",rid
                write(*,*) "rit=",rit1,"     ric=",ric1
                write(*,*) "initial guess for rtwi c_yst=",c_yst
                write(*,*) " "
            end if
        end if

        slq2=y/(b1*b1)
        c_y = y*(tpvot**2)

        if(c_y.ge.0) then
            c_y0=c_y
        else
            if(ri.lt.0) then
                if (mytid.eq.0) then
                    write(*,*) "c_y negative at negative ri"
                    write(*,*) "ri=",ri,"      c_y=",c_y
                    write(*,*) "first index=",index0,"second index=",index1
                    write(*,*) "unstable realizability limit unexpected:"
                    write(*,*) "stopping in mikesal2a."
                end if
                stop
            end if
        end if
        if(index0.eq.0) c_y00=c_y
        if((index0.eq.0).and.(index1.eq.0).and.  &
        (abs(c_y - c_yst0).gt.1.d-6)) then
            if (mytid.eq.0) then
                write(*,*) "inconsistency in neutral value of c_y"
                write(*,*) "value used =",c_yst0
                write(*,*) "value calculated =",c_y
                write(*,*) "program stopping in mikesal2a"
            end if
            stop
        end if

        x=rit1*y
        z=ric1*y
        call mksmshss1(x,y,z,sm,sh,ss)
        return
    end subroutine mikesal2a

!---------------------------------------------------------------------
!     end of subroutine mikesal2a
!---------------------------------------------------------------------






!---------------------------------------------------------------------
!     start of function fctysal
!---------------------------------------------------------------------

      function fctysal(y)                              

            use const_for_nasagissmixing, only: r8
            implicit none
            save
            real (r8),intent(in)    :: y
            real (r8)               :: fctysal

            real (r8)               :: x,z,sm,sh,ss
            real (r8)               :: ri,rid,ric,rit,rr

            common /mikebb/ri
            common /mikebbs/rid,ric,rit,rr

            x=rit*y
            z=ric*y
            call mksmshss1(x,y,z,sm,sh,ss)
           
            fctysal=2.D0/(sm-rit*sh-ric*ss)
            return                                          
      end function fctysal

!---------------------------------------------------------------------
!     start of function fctysal
!---------------------------------------------------------------------



         
!---------------------------------------------------------------------
!     start of subroutine mksmshss1
!---------------------------------------------------------------------

    subroutine mksmshss1(x,y,z,sm,sh,ss)
        use const_for_nasagissmixing, only: r8
!         use param_mod, only: mytid
        implicit none
        save

        real (r8),intent(inout)  :: x,y,z
        real (r8),intent(out) :: sm,sh,ss

        real (r8), parameter :: x1min=15.d0/8.d0
        real (r8), parameter :: calculation_epsilon=1.d-60
        integer, parameter :: ifnew=0

        real (r8) :: phi,x1,fc_twid,fc,psi
        real (r8) :: ri,rit,ric
        real (r8) :: eta_t,eta_s,eta_sum
        real (r8) :: sigma_m,st
        integer :: ifrecall

      integer :: mytid
      common /mytid_table2/ mytid

        ifrecall = 0

        if(ifrecall.eq.0) then
            if (mytid.eq.0) then
                write(*,*) " "
                write(*,*)"************************************************"
                write(*,*) " "
                write(*,*)"regularization used for dubovikov salinity model"
                write(*,*) "ifnew=",ifnew
                if (ifnew.eq.0) then
                    write(*,*) "1/1+x regularization"
                else if(ifnew.eq.1) then
                    write(*,*) "e^-x regularization"
                end if  
                write(*,*) " "
                write(*,*)"************************************************"
                write(*,*) " "
            end if
        end if
!*****************************************************************************c
        rit=x/y
        ric=z/y
        y = max(y,calculation_epsilon)
        x = rit*y
        z = ric*y

        x1=0.5d0*sqrt(y)
        if(x1.lt.x1min) then
            fc_twid=5.d0/9.d0
        else
            fc_twid=1.d0-5.d0/(3.d0*x1)+25.d0/(16.d0*x1**2)
        end if

        ri = rit + ric
        fc = 1.8d0*fc_twid
  
        if(ifnew.eq.0) then
            phi = 1.d0/(1.d0 + 0.08*x1**2*ri)
            psi = 1.d0/ &
            (  1.d0 & 
            + 0.16d0*x1**2*ri &
             + 0.013d0*x1**4*rit*ric*phi )
        else if(ifnew.eq.1) then    !990923 version
            phi =  exp(-0.08*x1**2*ri)
            psi = exp(-0.16d0*x1**2*ri  - 0.013d0*x1**4*rit*ric*phi)
        end if

        eta_t = 0.05d0*x1**2*rit*(1.d0 + (1.d0 +0.16d0*ric*x1**2)*psi) 
        eta_s = 0.05d0*x1**2*ric*(1.d0 + (1.d0 +0.16d0*rit*x1**2)*psi) 
        eta_sum = eta_t + eta_s
        if(ifnew.eq.0) then
            sigma_m = 1.d0/(1.d0 + eta_sum/fc)
        else if(ifnew.eq.1) then        !990923 version
            sigma_m = exp(-eta_sum/fc)
        end if

        sm=1.d0/25.d0*fc_twid*sigma_m
        sh=0.056d0*(1.d0 + 0.08d0*x1**2*ric*phi)*psi
        ss=0.056d0*(1.d0 + 0.08d0*x1**2*rit*phi)*psi

        if(ifnew.eq.1) then
            st=0.056d0*(1.d0 + 0.08d0*x1**2*ri*phi)*psi
        end if
    
        return
    end subroutine mksmshss1
!---------------------------------------------------------------------
!     end of subroutine mksmshss1
!---------------------------------------------------------------------                                                              
 


!-----------------------------------------------------------------------
!     end of parameter-free turbulence model subroutines (nmodel=2)
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
!     end of parameter-free turbulence model subroutines (nmodel=2)
!-----------------------------------------------------------------------



!030401-02Y Adaptation to the NCAR CSM Ocean Model of my 1D interpolation routine 
!     written for HYCOM [See NBp.030401-2to3]. I must make variables double precision.
      subroutine interp1d(x, &
                    x_1,slq2_1,sm_1,sh_1,ss_1, &
                    slq2,sm,sh,ss, &
                    ixmax,m,m0,delta)
            use const_for_nasagissmixing, only: r8
!             use param_mod, only: mytid   !will be removed
            implicit none
            save

            real (r8), intent(inout) :: x       !in parameters
            real (r8), intent(in) :: delta

            integer,intent(in) :: m,m0,ixmax
            real(r8),dimension(-m:m),intent(in) :: x_1,slq2_1,sm_1,sh_1,ss_1

            real (r8), intent(out) :: slq2,sm,sh,ss   !out parameters
            
            !temp variables
            integer :: lx1,l,lx0
            real (r8) :: deltaxta,deltax,dslq2_x,dsm_x,dsh_x,dss_x

            integer :: mytid
            common /mytid_table2/ mytid

! --- 6-990112-C030327  Subroutine for a ONE VARIABLE modular interpolation calculation.
!     x is the independent variable in the turbulence model calculation.
!
! --- 1D input array with table spacing:  x_1
! --- table limit value:                  ixmax
! --- 1D input arrays with table values:  slq2_1,sm_1,sh_1,ss_1
! --- Output interpolated values:         slq2,sm,sh,ss

!     print *,'start interp1d'
! --- 030326 Take values off the edge of the table back to the table limits.
            if(x.gt.x_1(ixmax)) then
            x = x_1(ixmax)
            else if(x.lt.x_1(-m)) then
            x = x_1(-m)
            end if 
!
! --- Interpolate points within the table range.
! --- Table index ranges from -m to m with equal spacing for -m0 to m0.
!
            if(abs(x).lt.x_1(m0)) then
! --- Find Interpolation points in the equally spaced x part of the table.
            lx1 = int(x/delta)+nint(sign(dfloat(1),x))
            else
! --- Find Interpolation points in unequally spaced x part of the table.
                  do l=m0,m
! --- Keep stepping until pass input value.
                  if(abs(x).lt.abs(x_1(nint(sign(dfloat(l),x))))) then
! --- Cover both positive and negative index cases.
                              lx1 = nint(sign(dfloat(l),x))
                              go to 250
                        end if
                  end do
!
! --- Special case where have a value which falls at the limit of the table.
                  lx0 = nint(sign(dfloat(m),x))
                  lx1 = lx0
                  go to 252
            end if 
            250   continue 
!
! --- Make lx0 one less or greater than lx1 according to sgn(x).
        lx0 = lx1 - nint(sign(dfloat(1),x)) 
!
        if(x.eq.0.) lx1 = 1
            252   continue
!
! --- Check that the x value falls within the interpolation interval.
            if((x.gt.0..and.  &
            (x.lt.x_1(lx0).or.x.gt.x_1(lx1))).or.  &
            (x.lt.0..and.  &
            (x.gt.x_1(lx0).or.x.lt.x_1(lx1)))) then
                  if (mytid.eq.0) then
                        write(*,*) "x is outside interpolation range in interp1d."
                        write(*,*) "x=  ",x,"lx0= ",lx0,"lx1= ",lx1
                        write(*,*) "x_1(lx0)=  ",x_1(lx0),"   x_1(lx1)= ",x_1(lx1)
                        write(*,*) "Program is stopping."
                  endif
                  stop
            end if

            270   continue
! --- Interpolate turbulence fields.
! --- Introduce table spacing variables.
            deltaxta = x_1(lx1) - x_1(lx0)
            deltax = x - x_1(lx0)

! --- slq2
! --- Set delta field to zero in special cases falling at limit of the table. 
            if(lx1.eq.lx0) then
                  dslq2_x = 0.
            else
                  dslq2_x = (slq2_1(lx1) - slq2_1(lx0))/ &
            deltaxta
            end if
            slq2   = slq2_1(lx0) + &
        dslq2_x*deltax

! --- sm
            if(lx1.eq.lx0) then
                  dsm_x   = 0.
            else
                  dsm_x = (sm_1(lx1) - sm_1(lx0))/ &
            deltaxta
            end if
            sm     = sm_1(lx0) + &
        dsm_x*deltax
!
! --- sh
            if(lx1.EQ.lx0) then
                  dsh_x   = 0.
            else
                  dsh_x = (sh_1(lx1) - sh_1(lx0))/ &
            deltaxta
            end if
            sh     = sh_1(lx0) + &
        dsh_x*deltax
!
! --- ss
            if(lx1.EQ.lx0) then
                  dss_x   = 0.
            else
                  dss_x = (ss_1(lx1) - ss_1(lx0))/ &
            deltaxta
            end if
            ss     = ss_1(lx0) + &
        dss_x*deltax
!     print *,'finish interp1d'
      return
      end subroutine interp1d




!-----------------------------------------------------------
!030401 Adaptation to the NCAR CSM Ocean Model of my 1D interpolation routine
!       written for HYCOM [See NBp.030401-2to3]. 
      subroutine interp1d_expabs(x, &
                    x_1,slq2_1,sm_1,sh_1,ss_1, &
                    slq2,sm,sh,ss, &
                    ixmax,m,m0,delta,rat)

            use const_for_nasagissmixing, only: r8
!             use param_mod, only: mytid   !will be removed
            implicit none
            save

            real (r8), intent(inout) :: x
            real (r8), intent(in) :: delta,rat

            integer,intent(in) :: m,m0,ixmax          
            real(r8),dimension(-m:m),intent(in) :: x_1,slq2_1,sm_1,sh_1,ss_1

            real (r8), intent(out) :: slq2,sm,sh,ss
            
            !temp variables
            integer :: lx1,lx0 
            real (r8) :: deltaxta,deltax
            real (r8) :: dslq2_x,dsm_x,dsh_x,dss_x
            real (r8) :: tabindx

            integer :: mytid
            common /mytid_table2/ mytid
!030129-030328 Subroutine for a faster interpolation calculation in the ifexpabstable=1 case.
! --- 6-990112-C030326  Subroutine for a ONE VARIABLE modular interpolation calculation.
!     x is the independent variable in the turbulence model calculation.
!
! --- 1D input array with table spacing:  x_1
! --- table limit value:                  ixmax
! --- 1D input arrays with table values:  slq2_1,sm_1,sh_1,ss_1
! --- Output interpolated values:         slq2,sm,sh,ss
!01107yXI Input value of ratio between entries  rat


!030326 Take values off the edge of the table back to the table.
!     print *,'start interp1d_expabs'
            if(x.gt.x_1(ixmax)) then
                  x = x_1(m)
            else if(x.lt.x_1(-m)) then
                  x = x_1(-m)
            end if
!*****C
!981019 Interpolate points within the table range.
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
            if(abs(x).lt.x_1(m0)) then
!981019-27 Find Interpolation points in the equally spaced x part of the table.
                  lx1 = int(x/delta)+nint(sign(dfloat(1),x))
!011107-08yXI Find Interpolation points in exponential absolute value spaced x part of the table.
            else if((abs(x)).ge.(x_1(m))) then
!981103 Special case where have a value which falls at the limit of the table.
                  lx0 = nint(sign(dfloat(m),x))
                  lx1 = lx0
                  go to 252
!*****C
            else
                  tabindx = sign( dfloat(m0) + ((log(abs(x)) - log(x_1(m0)))/log(rat)), x)
                  lx1 = int(tabindx)+nint(sign(dfloat(1),x))
!yXI
            end if
!011108yXI It is conceivable that rounding errors may in borderline cases 
!        throw the calculated table indices for x off by one.
!        Check and allow moving one to either side to take care of this.
            if((abs(x_1(lx1))).lt.(abs(x))) then
                  lx1 = lx1 + sign(1,lx1)
            else if((abs(x_1(lx1-sign(1,lx1)))).gt.(abs(x))) then
                  lx1 = lx1 - sign(1,lx1)
            end if
!yXI
            250   continue
!981019-27 Make lx0 one less or greater than lx1 according to sgn(x).
        lx0 = lx1 - nint(sign(dfloat(1),x)) 
!*****C
        if(x.eq.0.D0) lx1 = 1
            252   continue
!981019-28 Check that the x value falls within the interpolation interval.
            if((x.gt.0.D0.and.  &
            (x.lt.x_1(lx0).or.x.gt.x_1(lx1))).or.  &
            (x.lt.0.D0.and.  &
            (x.gt.x_1(lx0).or.x.lt.x_1(lx1)))) then
                  if (mytid.eq.0) then
                        write(*,*) "x is outside interpolation range in interp1d_expabs."
                        write(*,*) "delta=",delta
                        write(*,*) "m0=",m0," m=",m," rat=",rat
                        write(*,*) "x=  ",x," lx0= ",lx0," lx1= ",lx1
                        write(*,*) "x_1(lx0)=  ",x_1(lx0), "   x_1(lx1)= ",x_1(lx1)
                        write(*,*) "Program is stopping."
                  endif
                  stop
            end if
!*****C
!981019-27 Interpolate turbulence fields.
!981027-990112 Introduce table spacing variables.
            deltaxta = x_1(lx1) - x_1(lx0)
            deltax = x - x_1(lx0)
!     slq2
!981103 Set delta field to zero in special cases falling at limit of the table. 
            if(lx1.eq.lx0) then
                  dslq2_x = 0.D0
            else
                  dslq2_x = (slq2_1(lx1) - slq2_1(lx0))/ &
            deltaxta
            end if
            slq2   = slq2_1(lx0) + &
            dslq2_x*deltax
!     sm
            if(lx1.eq.lx0) then
                  dsm_x   = 0.D0
            else
                  dsm_x = (sm_1(lx1) - sm_1(lx0))/ &
                  deltaxta
            end if
            sm     = sm_1(lx0) + &
            dsm_x*deltax
!     sh
            if(lx1.eq.lx0) then
                  dsh_x   = 0.D0
            else
                  dsh_x = (sh_1(lx1) - sh_1(lx0))/ &
                  deltaxta
            end if
            sh     = sh_1(lx0) + &
            dsh_x*deltax
!     ss
            if(lx1.eq.lx0) then
                  dss_x   = 0.D0
            else
                  dss_x = (ss_1(lx1) - ss_1(lx0))/ &
                  deltaxta
            end if
            ss     = ss_1(lx0) + &
            dss_x*deltax


!     print *,'finish interp1d_expabs'
      return
      end subroutine interp1d_expabs



!-----------------------------------------------------------------------------
      subroutine interp2d(ri,rid, &
                    ri_1,rid_1,slq2_2,sm_2,sh_2,ss_2, &
                    slq2,sm,sh,ss, &
                    irimax,m,m0,delta)
            use const_for_nasagissmixing, only: r8
!             use param_mod, only: mytid   !will be removed
            implicit none
            save

            real (r8), intent(inout) :: ri,rid
            real (r8), intent(in) :: delta      

            integer,intent(in) :: m,m0
            integer,dimension(-m:m),intent(in) :: irimax
            real(r8),dimension(-m:m),intent(in) :: ri_1,rid_1
            real(r8),dimension(-m:m,-m:m),intent(in) :: slq2_2,sm_2,sh_2,ss_2

            real (r8), intent(out) :: slq2,sm,sh,ss
            
            !temp variables
            integer :: lrid1,l,lrid0,lri1,lri0
            real (r8) :: deltaridta,deltarita,deltarid,deltari
            real (r8) :: dslq2_rid,dslq2_ri,dsm_rid,dsm_ri,dsh_rid,dsh_ri,dss_rid,dss_ri

            integer :: mytid
            common /mytid_table2/ mytid

!981016-990112    Subroutine for a modular interpolation calculation.
!*
!     1D input arrays with table spacing:       ri_1,rid_1
!     1D input array with table limits:   irimax
!     2D input arrays with table values:  slq2_2,sm_2,sh_2,ss_2
!     Output interpolated values:         slq2,sm,sh,ss

!     Take values off the edge of the table back to the table on radial lines.
!981103 Must use ratio of Ri's taken *before* the cut-off has taken place.
            if(ri.gt.ri_1(m)) then
                  if(abs(rid).le.ri) then
                        rid = ri_1(m)*(rid/ri)
                        ri  = ri_1(m)
                  else if(rid.gt.ri) then
                        ri  = rid_1(m)*(ri/rid)
                        rid = rid_1(m)
                  else if(rid.lt.-ri) then
                        ri  = rid_1(-m)*(ri/rid)
                        rid = rid_1(-m)
                  end if
            else if(ri.lt.ri_1(-m)) then
                  if(abs(rid).le.-ri) then
                        rid = ri_1(-m)*(rid/ri)
                        ri  = ri_1(-m)
                  else if(rid.gt.-ri) then
                        ri  = rid_1(m)*(ri/rid)
                        rid = rid_1(m)
                  else if(rid.lt.ri) then
                        ri  = rid_1(-m)*(ri/rid)
                        rid = rid_1(-m)
                  end if
            else if(rid.gt.rid_1(m)) then
                  ri  = rid_1(m)*(ri/rid)
                  rid = rid_1(m)
            else if(rid.lt.rid_1(-m)) then
                  ri  = rid_1(-m)*(ri/rid)
                  rid = rid_1(-m)
            end if
!*****C
!981019 Interpolate points within the table range.
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
            if(abs(rid).lt.rid_1(m0)) then
!981019-27 Find Interpolation points in the equally spaced Ri_d part of the table.
                  lrid1 = int(rid/delta)+nint(sign(dfloat(1),rid))
            else
!981019 Find Interpolation points in unequally spaced Ri_d part of the table.
                  do l=m0,m
!981027 Keep stepping until pass input value.
                        if(abs(rid).lt.abs(rid_1(nint(sign(dfloat(l),rid))))) then
!981027 Cover both positive and negative index cases.
                              lrid1 = nint(sign(dfloat(l),rid))
!*****C
                        go to 250
                        end if
                  end do
!981103 Special case where have a value which falls at the limit of the table.
                  lrid0 = nint(sign(dfloat(m),rid))
                  lrid1 = lrid0
            go to 252
!*****C
            end if
            250   continue
!981019-27 Make lrid0 one less or greater than lrid1 according to sgn(rid).
        lrid0 = lrid1 - nint(sign(dfloat(1),rid)) 
!*****C
        if(rid.eq.0.D0) lrid1 = 1
            252   continue
!981019-27 Check that the Ri_d value falls within the interpolation interval.
            if((rid.gt.0.D0.and.  &
            (rid.lt.rid_1(lrid0).or.rid.gt.rid_1(lrid1))).or.  &
            (rid.lt.0.D0.and.  &
            (rid.gt.rid_1(lrid0).or.rid.lt.rid_1(lrid1)))) then
                  if (mytid.eq.0) then
                  write(*,*) "Ri_d is outside interpolation range in interp2d."
                  write(*,*) "rid=  ",rid,"lrid0= ",lrid0,"lrid1= ",lrid1
                  write(*,*) "rid_1(lrid0)=  ",rid_1(lrid0),"   rid_1(lrid1)= ",rid_1(lrid1)
                  write(*,*) "Program is stopping."
                  endif
                  stop
            end if
!*****C
!C981022    Artificially reduce Ri if it threatens to surpass Ri_max(Ri_d).
!C    This is to conform to the 1D table's realizability limit treatment. 
!     IF(ri.GT.MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) THEN
!       ri = MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))
!     END IF
!C*****C
!981110 Set turbulence to zero if Ri threatens to surpass the realizability limit.
        if(ri.gt.min(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) then
          slq2=0.D0
          sm = 0.D0
          sh = 0.D0
          ss = 0.D0
          return
        end if
!*****C
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
            if(abs(ri).lt.ri_1(m0)) then
!981019-27 Find Interpolation points in the equally spaced Ri part of the table.
                  lri1 = int(ri/delta)+nint(sign(dfloat(1),ri)) 
            else
!981019 Find Interpolation points in unequally spaced Ri part of the table.
                  do l=m0,m
!981027 Keep stepping until pass input value.
                        if(abs(ri).lt.abs(ri_1(nint(sign(dfloat(l),ri))))) then
!981027 Cover both positive and negative index cases.
                              lri1 = nint(sign(dfloat(l),ri))
                        go to 270
                        end if
                  end do
!981103 Special case where have a value which falls at the limit of the table.
                  lri0 = nint(sign(dfloat(m),ri))
                  lri1 = lri0
                  go to 272
!*****C
            270   continue
            end if
!981019-27 Make lri0 one less or greater than lri1 according to sgn(ri).
        lri0 = lri1 - nint(sign(dfloat(1),ri)) 
!*****C
        if(ri.eq.0.D0) lri1 = 1
            272   continue
!981019-27 Check that the Ri_d value falls within the interpolation interval.
            if((ri.gt.0.D0.and.(ri.lt.ri_1(lri0).or.ri.gt.ri_1(lri1))).or.  &
            (ri.lt.0.D0.and.(ri.gt.ri_1(lri0).or.ri.lt.ri_1(lri1)))) then
                  if (mytid.eq.0) then
                        write(*,*) "Ri is outside interpolation range in interp2d."
                        write(*,*) "ri=  ",ri,"lri0= ",lri0,"lri1= ",lri1
                        write(*,*) "ri_1(lri0)=  ",ri_1(lri0), "   ri_1(lri1)= ",ri_1(lri1)
                        write(*,*) "Program is stopping."
                  endif
                  stop
            end if
!*****C
!981019-27 Interpolate turbulence fields.
!981027-990112 Introduce table spacing variables.
            deltaridta = rid_1(lrid1) - rid_1(lrid0)
            deltarita  = ri_1(lri1)  - ri_1(lri0)
            deltarid = rid - rid_1(lrid0)
            deltari  = ri - ri_1(lri0)
!     slq2
!981103 Set delta field to zero in special cases falling at limit of the table. 
            if(lrid1.eq.lrid0) then
                  dslq2_rid = 0.D0
            else
                  dslq2_rid = (slq2_2(lri0,lrid1) - slq2_2(lri0,lrid0))/ &
            deltaridta
            end if            
            if(lri1.eq.lri0) then
                  dslq2_ri  = 0.D0
            else
                  dslq2_ri = (slq2_2(lri1,lrid0) - slq2_2(lri0,lrid0))/ &
                  deltarita
            end if            
            slq2   = slq2_2(lri0,lrid0) + &
        dslq2_ri*deltari + dslq2_rid*deltarid
            
!     sm
            if(lrid1.eq.lrid0) then
                  dsm_rid   = 0.D0
            else
                  dsm_rid = (sm_2(lri0,lrid1) - sm_2(lri0,lrid0))/ &
                  deltaridta
            end if
            if(lri1.eq.lri0) then
                  dsm_ri    = 0.D0
            else
                  dsm_ri = (sm_2(lri1,lrid0) - sm_2(lri0,lrid0))/ &
                  deltarita
            end if
            sm     = sm_2(lri0,lrid0) + &
        dsm_ri*deltari + dsm_rid*deltarid
            
!     sh
            if(lrid1.eq.lrid0) then
                  dsh_rid   = 0.D0
            else
                  dsh_rid = (sh_2(lri0,lrid1) - sh_2(lri0,lrid0))/ &
                  deltaridta
            end if
            if(lri1.eq.lri0) then
                  dsh_ri    = 0.D0
            else
                  dsh_ri = (sh_2(lri1,lrid0) - sh_2(lri0,lrid0))/ &
                  deltarita
            end if
            sh     = sh_2(lri0,lrid0) + &
        dsh_ri*deltari + dsh_rid*deltarid
            
!     ss
            if(lrid1.eq.lrid0) then
                  dss_rid   = 0.D0
            else
                  dss_rid = (ss_2(lri0,lrid1) - ss_2(lri0,lrid0))/ &
                  deltaridta
            end if
            if(lri1.eq.lri0) then
                  dss_ri    = 0.D0
            else
                  dss_ri = (ss_2(lri1,lrid0) - ss_2(lri0,lrid0))/ &
                  deltarita
            end if
            ss     = ss_2(lri0,lrid0) + &
        dss_ri*deltari + dss_rid*deltarid

      return
      end subroutine interp2d





      subroutine interp2d_expabs(ri,rid, &
                    ri_1,rid_1,slq2_2,sm_2,sh_2,ss_2, &
                    slq2,sm,sh,ss, &
                    irimax,m,m0,delta,rat)

            use const_for_nasagissmixing, only: r8
!             use param_mod, only: mytid   !will be removed
            implicit none
            save

            real (r8), intent(inout) :: ri,rid
            real (r8), intent(in) :: delta,rat

            integer,intent(in) :: m,m0
            integer,dimension(-m:m),intent(in) :: irimax
            real(r8),dimension(-m:m),intent(in) :: ri_1,rid_1
            real(r8),dimension(-m:m,-m:m),intent(in) :: slq2_2,sm_2,sh_2,ss_2

            real (r8), intent(out) :: slq2,sm,sh,ss
            
            !temp variables
            integer :: lrid1,lrid0,lri1,lri0
            real (r8) :: tabindrid,tabindri,deltaridta,deltarita,deltarid,deltari
            real (r8) :: dslq2_rid,dslq2_ri,dsm_rid,dsm_ri,dsh_rid,dsh_ri,dss_rid,dss_ri

            integer :: mytid
            common /mytid_table2/ mytid

!030123z Remove the "sys_flush" from the model E interpolation routine that avoids
!      in the exponential absolute value case the super-time-consuming stepping 
!      through the nonlinear part of the table by estimating the table index near
!      the value to be interpolated in order to use this more efficient routine
!      in my traditional NCOM runs.
!011107-08yXI Subroutine for an interpolation calculation for a table whose nonlinear part 
!yXI          is exponential in the absolute value of the index. Based on:
!981016-990112    Subroutine for a modular interpolation calculation.
!*
!01107yXI Input value of ratio between entries  rat
!     1D input arrays with table spacing:       ri_1,rid_1
!     1D input array with table limits:   irimax
!     2D input arrays with table values:  slq2_2,sm_2,sh_2,ss_2
!     Output interpolated values:         slq2,sm,sh,ss

!     Take values off the edge of the table back to the table on radial lines.
!981103 Must use ratio of Ri's taken *before* the cut-off has taken place.
            if(ri.gt.ri_1(m)) then
                  if(abs(rid).le.ri) then
                        rid = ri_1(m)*(rid/ri)
                        ri  = ri_1(m)
                  else if(rid.gt.ri) then
                        ri  = rid_1(m)*(ri/rid)
                        rid = rid_1(m)
                  else if(rid.lt.-ri) then
                        ri  = rid_1(-m)*(ri/rid)
                        rid = rid_1(-m)
                  end if
            else if(ri.lt.ri_1(-m)) then
                  if(abs(rid).le.-ri) then
                        rid = ri_1(-m)*(rid/ri)
                        ri  = ri_1(-m)
                  else if(rid.gt.-ri) then
                        ri  = rid_1(m)*(ri/rid)
                        rid = rid_1(m)
                  else if(rid.lt.ri) then
                        ri  = rid_1(-m)*(ri/rid)
                        rid = rid_1(-m)
                  end if
            else if(rid.gt.rid_1(m)) then
                  ri  = rid_1(m)*(ri/rid)
                  rid = rid_1(m)
            else if(rid.lt.rid_1(-m)) then
                  ri  = rid_1(-m)*(ri/rid)
                  rid = rid_1(-m)
            end if
!*****C
!981019 Interpolate points within the table range.
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
            if(abs(rid).lt.rid_1(m0)) then
!981019-27 Find Interpolation points in the equally spaced Ri_d part of the table.
                  lrid1 = int(rid/delta)+nint(sign(dfloat(1),rid))
!011107-08yXI Find Interpolation points in exponential absolute value spaced Ri_d part of the table.
            else if((abs(rid)).ge.(rid_1(m))) then
!981103 Special case where have a value which falls at the limit of the table.
                  lrid0 = nint(sign(dfloat(m),rid))
                  lrid1 = lrid0
                  go to 252
!*****C
            else
                  tabindrid = sign(dfloat(m0) + ((log(abs(rid)) - log(rid_1(m0)))/log(rat)), rid)
                  lrid1 = int(tabindrid)+nint(sign(dfloat(1),rid))
!yXI
            end if
!011108yXI It is conceivable that rounding errors may in borderline cases 
!        throw the calculated table indices for Ri_d off by one.
!        Check and allow moving one to either side to take care of this.
            if((abs(rid_1(lrid1))).lt.(abs(rid))) then
                  lrid1 = lrid1 + sign(1,lrid1)
            else if((abs(rid_1(lrid1-sign(1,lrid1)))).gt.(abs(rid))) then
                  lrid1 = lrid1 - sign(1,lrid1)
            end if
!yXI
            250   continue
!981019-27 Make lrid0 one less or greater than lrid1 according to sgn(rid).
        lrid0 = lrid1 - nint(sign(dfloat(1),rid)) 
!*****C
        if(rid.eq.0.D0) lrid1 = 1
            252   continue
!981019-27 Check that the Ri_d value falls within the interpolation interval.
            if((rid.gt.0.D0.and.  &
            (rid.lt.rid_1(lrid0).or.rid.gt.rid_1(lrid1))).or.  &
            (rid.lt.0.D0.and.  &
            (rid.gt.rid_1(lrid0).or.rid.lt.rid_1(lrid1)))) then
                  if (mytid.eq.0) then
                        write(*,*) "Ri_d is outside interpolation range in interp2d."
                        write(*,*) "rid=  ",rid,"lrid0= ",lrid0,"lrid1= ",lrid1
                        write(*,*) "rid_1(lrid0)=  ",rid_1(lrid0),"   rid_1(lrid1)= ",rid_1(lrid1)
                        write(*,*) "Program is stopping."
                  endif
                  stop
            end if
!*****C
!C981022    Artificially reduce Ri if it threatens to surpass Ri_max(Ri_d).
!C    This is to conform to the 1D table's realizability limit treatment. 
!     IF(ri.GT.MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) THEN
!       ri = MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))
!     END IF
!C*****C
!981110 Set turbulence to zero if Ri threatens to surpass the realizability limit.
        if(ri.gt.min(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) then
                  slq2=0.D0
                  sm = 0.D0
                  sh = 0.D0
                  ss = 0.D0
                  return
        end if
!*****C
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
            if(abs(ri).lt.ri_1(m0)) then
!981019-27 Find Interpolation points in the equally spaced Ri part of the table.
                  lri1 = int(ri/delta)+nint(sign(dfloat(1),ri)) 
!011107-08yXI Find Interpolation points in exponential absolute value spaced Ri part of the table.
            else if((abs(ri)).ge.(ri_1(m))) then
!981103 Special case where have a value which falls at the limit of the table.
                  lri0 = nint(sign(dfloat(m),ri))
                  lri1 = lri0
                  go to 272
!*****C
            else
                  tabindri = sign( dfloat(m0) + ((log(abs(ri)) - log(ri_1(m0)))/log(rat)), ri)
                  lri1 = int(tabindri)+nint(sign(dfloat(1),ri))
!yXI
            270   continue
            end if
!011108yXI It is conceivable that rounding errors will in borderline cases 
!        throw the calculated table indices for Ri off by one.
!        Check and allow moving one to either side to take care of this.
            if((abs(ri_1(lri1))).lt.(abs(ri))) then
                  lri1 = lri1 + sign(1,lri1)
            else if((abs(ri_1(lri1-sign(1,lri1)))).gt.(abs(ri))) then
                  lri1 = lri1 - sign(1,lri1)
            end if
!yXI
!981019-27 Make lri0 one less or greater than lri1 according to sgn(ri).
        lri0 = lri1 - nint(sign(dfloat(1),ri)) 
!*****C
        if(ri.eq.0.D0) lri1 = 1
            272   continue
!981019-27 Check that the Ri_d value falls within the interpolation interval.
            if((ri.gt.0.D0.and.(ri.lt.ri_1(lri0).or.ri.gt.ri_1(lri1))).or.  &
            (ri.lt.0.D0.and.(ri.gt.ri_1(lri0).or.ri.lt.ri_1(lri1)))) then
                  if (mytid.eq.0) then
                        write(*,*) "Ri is outside interpolation range in interp2d."
                        write(*,*) "ri=  ",ri,"lri0= ",lri0,"lri1= ",lri1
                        write(*,*) "ri_1(lri0)=  ",ri_1(lri0),"   ri_1(lri1)= ",ri_1(lri1)
                        write(*,*) "Program is stopping."
                  endif
                  stop
            end if
!*****C
!981019-27 Interpolate turbulence fields.
!981027-990112 Introduce table spacing variables.
            deltaridta = rid_1(lrid1) - rid_1(lrid0)
            deltarita  = ri_1(lri1)  - ri_1(lri0)
            deltarid = rid - rid_1(lrid0)
            deltari  = ri - ri_1(lri0)
!     slq2
!981103 Set delta field to zero in special cases falling at limit of the table. 
            if(lrid1.eq.lrid0) then
                  dslq2_rid = 0.D0
            else
                  dslq2_rid = (slq2_2(lri0,lrid1) - slq2_2(lri0,lrid0))/ &
            deltaridta
            end if
            if(lri1.eq.lri0) then
                  dslq2_ri  = 0.D0
            else
                  dslq2_ri = (slq2_2(lri1,lrid0) - slq2_2(lri0,lrid0))/ &
            deltarita
            end if
            slq2   = slq2_2(lri0,lrid0) + &
        dslq2_ri*deltari + dslq2_rid*deltarid
            
!     sm
            if(lrid1.eq.lrid0) then
                  dsm_rid   = 0.D0
            else
                  dsm_rid = (sm_2(lri0,lrid1) - sm_2(lri0,lrid0))/ &
            deltaridta
            end if
            if(lri1.eq.lri0) then
                  dsm_ri    = 0.D0
            else
                  dsm_ri = (sm_2(lri1,lrid0) - sm_2(lri0,lrid0))/ &
            deltarita
            end if
            sm     = sm_2(lri0,lrid0) + &
        dsm_ri*deltari + dsm_rid*deltarid
            
!     sh
            if(lrid1.eq.lrid0) then
                  dsh_rid   = 0.D0
            else
                  dsh_rid = (sh_2(lri0,lrid1) - sh_2(lri0,lrid0))/ &
            deltaridta
            end if
            if(lri1.eq.lri0) then
                  dsh_ri    = 0.D0
            else
                  dsh_ri = (sh_2(lri1,lrid0) - sh_2(lri0,lrid0))/ &
            deltarita
            end if
            sh     = sh_2(lri0,lrid0) + &
        dsh_ri*deltari + dsh_rid*deltarid
            
!     ss
            if(lrid1.EQ.lrid0) then
                  dss_rid   = 0.D0
            else
                  dss_rid = (ss_2(lri0,lrid1) - ss_2(lri0,lrid0))/ &
            deltaridta
            end if
            if(lri1.EQ.lri0) then
                  dss_ri    = 0.D0
            else
                  dss_ri = (ss_2(lri1,lrid0) - ss_2(lri0,lrid0))/ &
            deltarita
            end if
            ss     = ss_2(lri0,lrid0) + &
        dss_ri*deltari + dss_rid*deltarid

      return
      end subroutine interp2d_expabs
