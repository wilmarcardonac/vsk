Program vsk 

!####################
! LOAD NEEDED MODULES 
!####################

    use fiducial
    use arrays
    use functions 

!#################################
! DECLARE VARIABLES AND PARAMETERS
!#################################

    Implicit none

!!$    Integer*4 :: m,n,i                                      ! INTERGER FOR SHORT LOOPS 
!!$    Integer*4 :: seed1,seed2                                      ! SEEDS FOR RANDOM NUMBER GENERATOR 
!!$    Integer*4 :: number_accepted_points,number_rejected_points    ! MCMC PARAMETERS
!!$    Integer*4 :: weight                                           ! IT COUNTS THE NUMBER OF STEPS TAKEN BEFORE MOVING TO A NEW POINT IN MCMC 
!!$
!!$    Real*4 :: average_acceptance_probability           ! SAVES ACCEPTANCE PROBABILITY 
!!$    Real*4 :: genunf                            ! RANDOM UNIFOR DEVIATES 
!!$
!!$    Real*8 :: random_uniform                           ! SAVES RANDOM UNIFORM DEVIATE BETWEEN 0 AND 1 
!!$    Real*8 :: old_loglikelihood,current_loglikelihood  ! TEMPORALY SAVES LIKELIHOOD VALUES 
!!$    Real*8 :: chi2SNIabestfit                          ! SAVES CHI2 OF SNIA AT THE BEST FIT
!!$    Real*4,dimension(number_of_parameters*(number_of_parameters+3)/2 + 1) :: parm     ! ARRAY NEEDED BY RANDON NUMBER GENERATOR 
!!$    Real*4,dimension(number_of_parameters) :: work,x_old,x_new                        ! ARRAYS NEEDED BY RANDOM NUMBER GENERATOR 
!!$    Real*8,dimension(number_of_parameters) :: bestfit,means                           ! SAVES BESTFIT AND MEANS OF PARAMETERS 
!!$
!!$    Logical :: not_good_aap,non_plausible_parameters   ! IT CONTROLS PLAUSIBLE VALUES OF COSMOLOGICAL PARAMETERS
!!$    Logical,dimension(number_of_parameters) :: plausibility  
!!$
!!$    Character(len=10) :: string ! STORES STRINGS FOR INTEGERS
!!$    Character(len=12),dimension(number_hyperparameters) :: alpha_string
!!$    Character(len=12),dimension(number_of_parameters) :: paramnames,latexname
!!$    Character(len=5) :: galaxy

!##########################################################
! Assignmentsc AND INITIALIZATION OF RANDOM NUMBER GENERATOR
!##########################################################

!!$    galaxy = host(1)
!!$
!!$    weight = 1
!!$
!!$    number_rejected_points = 0
!!$
!!$    number_accepted_points = 0
!!$
!!$    call initialize()                               ! INITIALIZE RANDOM NUMBER GENERATOR
!!$
!!$    call phrtsd(phrase,seed1,seed2)                 ! GENERATE SEEDS FOR RANDOM NUMBERS FROM PHRASE
!!$
!!$    call set_initial_seed(seed1,seed2)              ! SET INITIAL SEEDS FOR RANDOM NUMBER GENERATOR
!!$
!!$    ! ALLOCATING MEMORY FOR POINTS IN PARAMETER SPACE AND ACCEPTANCE PROBABILITY
!!$    allocate (old_point(1:number_of_parameters),current_point(1:number_of_parameters),&
!!$    acceptance_probability(number_iterations),Covgauss(number_of_parameters,number_of_parameters),&
!!$    Covguess(number_of_parameters,number_of_parameters),stat = status1)
!!$
!!$    call set_covariance_matrix()
!!$
!!$    call read_data()
    
!##################################
! MARKOV CHAIN MONTE CARLO ANALYSIS
!##################################

    call generate_gaussian_cmb_map()

!!$    write(UNIT_EXE_FILE,*) 'STARTING MCMC ANALYSIS'
!!$
!!$    If (testing_Gaussian_likelihood) then
!!$
!!$       write(UNIT_EXE_FILE,*) 'TESTING CODE WITH GAUSSIAN LIKELIHOOD'
!!$
!!$       open(UNIT_MCMC_FINAL_FILE,file='./output/chains/mcmc_final_output.txt')
!!$
!!$       open(UNIT_RANGES_FILE,file='./output/chains/mcmc_final_output.ranges')    !    OPEN FILE WITH HARD BOUNDS NEEDED BY GETDIST
!!$
!!$       write(UNIT_RANGES_FILE,*) 'A    N    N '
!!$
!!$       write(UNIT_RANGES_FILE,*) 'bw    N    N '
!!$
!!$       !    write(17,*) 'sigma_int    N    N '
!!$
!!$       close(UNIT_RANGES_FILE)
!!$
!!$       Do i=1,number_of_parameters
!!$
!!$          x_old(i) = genunf(-1.,1.)
!!$            
!!$       End Do
!!$
!!$       Do m=1,number_of_parameters
!!$    
!!$          If (m .eq. 4) then
!!$
!!$             old_point(m) = dble(x_old(m)) !exp(dble(x_old(m)))/(1.d1**1.d1)
!!$
!!$          else
!!$
!!$             old_point(m) = dble(x_old(m))
!!$
!!$          End If
!!$
!!$       End Do
!!$
!!$       old_loglikelihood = log_Gaussian_likelihood(old_point)
!!$
!!$    Else 
!!$       !###########################################
!!$       ! GENERATE A RANDOM POINT IN PARAMETER SPACE
!!$       ! RANDOM NUMBER GENERATOR WORKS WITH SINGLE PRECISION WHEREAS THIS CODES USES DOUBLE PRECISION; CHANGES ARE CORRESPONDINGLY MADE.
!!$       !################################################################################################################################
!!$
!!$       If (start_from_fiducial) then
!!$
!!$          write(UNIT_EXE_FILE,*) 'STARTING FROM FIDUCIAL POINT'
!!$
!!$          If (doing_R11_analysis) then
!!$
!!$             If (include_only_cepheids) then
!!$
!!$                If (all_R11_hosts) then
!!$
!!$                   old_point(1) = prior_mu1
!!$
!!$                   old_point(2) = prior_mu2
!!$
!!$                   old_point(3) = prior_mu3 
!!$
!!$                   old_point(4) = prior_mu4
!!$
!!$                   old_point(5) = prior_mu5
!!$
!!$                   old_point(6) = prior_mu6
!!$
!!$                   old_point(7) = prior_mu7
!!$
!!$                   old_point(8) = prior_mu8
!!$
!!$                   old_point(9) = prior_mu9
!!$
!!$                   old_point(10) = prior_zpw
!!$
!!$                   old_point(11) = prior_bw
!!$
!!$                   old_point(12) = prior_Zw
!!$
!!$                Else
!!$
!!$                   old_point(1) = prior_zpw
!!$
!!$                   old_point(2) = prior_bw
!!$
!!$                   old_point(3) = prior_Zw
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         old_point(1) = prior_mu1
!!$
!!$                         old_point(2) = prior_mu2
!!$
!!$                         old_point(3) = prior_mu3 
!!$
!!$                         old_point(4) = prior_mu4
!!$
!!$                         old_point(5) = prior_mu5
!!$
!!$                         old_point(6) = prior_mu6
!!$
!!$                         old_point(7) = prior_mu7
!!$
!!$                         old_point(8) = prior_mu8
!!$
!!$                         old_point(9) = prior_mu9
!!$
!!$                         old_point(10) = prior_mu10
!!$
!!$                         old_point(11) = prior_Mw
!!$
!!$                         old_point(12) = prior_bw
!!$
!!$                         old_point(13) = prior_H0
!!$
!!$                         old_point(14) = prior_Zw
!!$
!!$                         old_point(15) = a_v
!!$
!!$                         old_point(16) = a_cal
!!$                        
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         old_point(1) = prior_mu1
!!$
!!$                         old_point(2) = prior_mu2
!!$
!!$                         old_point(3) = prior_mu3 
!!$
!!$                         old_point(4) = prior_mu4
!!$
!!$                         old_point(5) = prior_mu5
!!$
!!$                         old_point(6) = prior_mu6
!!$
!!$                         old_point(7) = prior_mu7
!!$
!!$                         old_point(8) = prior_mu8
!!$
!!$                         old_point(9) = prior_mu9
!!$
!!$                         old_point(10) = prior_mu10
!!$
!!$                         old_point(11) = prior_Mw
!!$
!!$                         old_point(12) = prior_bw
!!$
!!$                         old_point(13) = prior_H0
!!$
!!$                         old_point(14) = prior_Zw
!!$
!!$                         old_point(15) = a_v
!!$
!!$                         old_point(16) = a_cal
!!$                        
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         old_point(1) = prior_mu1
!!$
!!$                         old_point(2) = prior_mu2
!!$
!!$                         old_point(3) = prior_mu3 
!!$
!!$                         old_point(4) = prior_mu4
!!$
!!$                         old_point(5) = prior_mu5
!!$
!!$                         old_point(6) = prior_mu6
!!$
!!$                         old_point(7) = prior_mu7
!!$
!!$                         old_point(8) = prior_mu8
!!$
!!$                         old_point(9) = prior_mu9
!!$
!!$                         old_point(10) = prior_Mw
!!$
!!$                         old_point(11) = prior_bw
!!$
!!$                         old_point(12) = prior_H0
!!$
!!$                         old_point(13) = prior_Zw
!!$
!!$                         old_point(14) = a_v
!!$
!!$                         old_point(15) = a_cal
!!$                        
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         old_point(1) = prior_mu1
!!$
!!$                         old_point(2) = prior_mu2
!!$
!!$                         old_point(3) = prior_mu3 
!!$
!!$                         old_point(4) = prior_mu4
!!$
!!$                         old_point(5) = prior_mu5
!!$
!!$                         old_point(6) = prior_mu6
!!$
!!$                         old_point(7) = prior_mu7
!!$
!!$                         old_point(8) = prior_mu8
!!$
!!$                         old_point(9) = prior_mu9
!!$
!!$                         old_point(10) = prior_mu10
!!$
!!$                         old_point(11) = prior_Mw
!!$
!!$                         old_point(12) = prior_bw
!!$
!!$                         old_point(13) = prior_H0
!!$
!!$                         old_point(14) = prior_Zw
!!$
!!$                         old_point(15) = a_v
!!$
!!$                         old_point(16) = a_cal
!!$                        
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         old_point(1) = prior_mu1
!!$
!!$                         old_point(2) = prior_mu2
!!$
!!$                         old_point(3) = prior_mu3 
!!$
!!$                         old_point(4) = prior_mu4
!!$
!!$                         old_point(5) = prior_mu5
!!$
!!$                         old_point(6) = prior_mu6
!!$
!!$                         old_point(7) = prior_mu7
!!$
!!$                         old_point(8) = prior_mu8
!!$
!!$                         old_point(9) = prior_mu9
!!$
!!$                         old_point(10) = prior_Mw
!!$
!!$                         old_point(11) = prior_bw
!!$
!!$                         old_point(12) = prior_H0
!!$
!!$                         old_point(13) = prior_Zw
!!$
!!$                         old_point(14) = a_v
!!$
!!$                         old_point(15) = a_cal
!!$                        
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         old_point(1) = prior_mu1
!!$
!!$                         old_point(2) = prior_mu2
!!$
!!$                         old_point(3) = prior_mu3 
!!$
!!$                         old_point(4) = prior_mu4
!!$
!!$                         old_point(5) = prior_mu5
!!$
!!$                         old_point(6) = prior_mu6
!!$
!!$                         old_point(7) = prior_mu7
!!$
!!$                         old_point(8) = prior_mu8
!!$
!!$                         old_point(9) = prior_mu9
!!$
!!$                         old_point(10) = prior_mu10
!!$
!!$                         old_point(11) = prior_zpwLMC
!!$
!!$                         old_point(12) = prior_bw
!!$
!!$                         old_point(13) = prior_H0
!!$
!!$                         old_point(14) = prior_Zw
!!$
!!$                         old_point(15) = a_v
!!$
!!$                         old_point(16) = a_cal
!!$                        
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         old_point(1) = prior_mu1
!!$
!!$                         old_point(2) = prior_mu2
!!$
!!$                         old_point(3) = prior_mu3 
!!$
!!$                         old_point(4) = prior_mu4
!!$
!!$                         old_point(5) = prior_mu5
!!$
!!$                         old_point(6) = prior_mu6
!!$
!!$                         old_point(7) = prior_mu7
!!$
!!$                         old_point(8) = prior_mu8
!!$
!!$                         old_point(9) = prior_mu9
!!$
!!$                         old_point(10) = prior_zpw
!!$
!!$                         old_point(11) = prior_bw
!!$
!!$                         old_point(12) = prior_H0
!!$
!!$                         old_point(13) = prior_Zw
!!$
!!$                         old_point(14) = a_v
!!$
!!$                      Else
!!$
!!$                         old_point(1) = prior_mu1
!!$
!!$                         old_point(2) = prior_mu2
!!$
!!$                         old_point(3) = prior_mu3 
!!$
!!$                         old_point(4) = prior_mu4
!!$
!!$                         old_point(5) = prior_mu5
!!$
!!$                         old_point(6) = prior_mu6
!!$
!!$                         old_point(7) = prior_mu7
!!$
!!$                         old_point(8) = prior_mu8
!!$
!!$                         old_point(9) = prior_mu9
!!$
!!$                         old_point(10) = prior_zpw
!!$
!!$                         old_point(11) = prior_bw
!!$
!!$                         old_point(12) = prior_H0
!!$
!!$                         old_point(13) = prior_Zw
!!$
!!$                         old_point(14) = a_v
!!$                        
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         old_point(1) = prior_zpH
!!$
!!$                         old_point(2) = prior_zpH
!!$
!!$                         old_point(3) = prior_zpH 
!!$
!!$                         old_point(4) = prior_zpH
!!$
!!$                         old_point(5) = prior_zpH
!!$
!!$                         old_point(6) = prior_zpH
!!$
!!$                         old_point(7) = prior_zpH
!!$
!!$                         old_point(8) = prior_zpH
!!$
!!$                         old_point(9) = prior_zpH
!!$
!!$                         old_point(10) = prior_bH
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'
!!$
!!$                   stop
!!$
!!$                End If ! OF ANCHORS
!!$    
!!$             End If ! OF CEPHEIDS
!!$
!!$          Else
!!$
!!$             old_point(1) = prior_A         ! A 
!!$             old_point(2) = prior_bw        ! bw 
!!$             !old_point(3) = prior_sigma_int ! sigma_int 
!!$ 
!!$          End If
!!$
!!$          If (hyperparameters_as_mcmc) then
!!$             ! SETTING INITIAL POINT FOR HYPER-PARAMETERS
!!$             Do m=number_model_parameters+1,number_of_parameters
!!$
!!$                old_point(m) = prior_alpha_j
!!$
!!$             End Do
!!$            
!!$          End If
!!$    
!!$          Do m=1,number_of_parameters
!!$
!!$             If (m .gt. number_model_parameters) then
!!$                    
!!$                If (using_jeffreys_prior) then
!!$
!!$                   x_old(m) = real(log10(old_point(m)))    ! CONVERT TO log10(\alpha_j)
!!$
!!$                Else
!!$
!!$                   x_old(m) = real(old_point(m))
!!$
!!$                End If
!!$
!!$             else
!!$
!!$                x_old(m) = real(old_point(m))
!!$
!!$             End If
!!$
!!$          End Do
!!$
!!$       Else
!!$
!!$          If (doing_R11_analysis) then
!!$
!!$             If (include_only_cepheids) then
!!$
!!$                If (all_R11_hosts) then
!!$
!!$                   x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                   x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                   x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                   x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                   x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                   x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                   x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                   x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                   x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                   x_old(10) = genunf(real(prior_zpw - sigma_zpw),real(prior_zpw + sigma_zpw))
!!$
!!$                   x_old(11) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$                   
!!$                   x_old(12) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$                   
!!$                Else
!!$
!!$                   x_old(1) = genunf(real(prior_zpw - sigma_zpw),real(prior_zpw + sigma_zpw))
!!$
!!$                   x_old(2) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                   x_old(3) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$                     
!!$                         x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                         x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                         x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                         x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                         x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                         x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                         x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                         x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                         x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                         x_old(10) = genunf(real(prior_mu10 - sigma_mu10),real(prior_mu10 + sigma_mu10))
!!$
!!$                         x_old(11) = genunf(real(prior_Mw - sigma_Mw),real(prior_Mw + sigma_Mw))
!!$
!!$                         x_old(12) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                         x_old(13) = genunf(real(prior_H0 - sigma_H0),real(prior_H0 + sigma_H0))
!!$
!!$                         x_old(14) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                         x_old(15) = genunf(real(a_v - sigma_a_v),real(a_v + sigma_a_v))
!!$
!!$                         x_old(16) = genunf(real(a_cal - sigma_a_cal),real(a_cal + sigma_a_cal))
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$                   
!!$                Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$                     
!!$                         x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                         x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                         x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                         x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                         x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                         x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                         x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                         x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                         x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                         x_old(10) = genunf(real(prior_mu10 - sigma_mu10),real(prior_mu10 + sigma_mu10))
!!$
!!$                         x_old(11) = genunf(real(prior_Mw - sigma_Mw),real(prior_Mw + sigma_Mw))
!!$
!!$                         x_old(12) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                         x_old(13) = genunf(real(prior_H0 - sigma_H0),real(prior_H0 + sigma_H0))
!!$
!!$                         x_old(14) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                         x_old(15) = genunf(real(a_v - sigma_a_v),real(a_v + sigma_a_v))
!!$
!!$                         x_old(16) = genunf(real(a_cal - sigma_a_cal),real(a_cal + sigma_a_cal))
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$                     
!!$                         x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                         x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                         x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                         x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                         x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                         x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                         x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                         x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                         x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                         x_old(10) = genunf(real(prior_Mw - sigma_Mw),real(prior_Mw + sigma_Mw))
!!$
!!$                         x_old(11) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                         x_old(12) = genunf(real(prior_H0 - sigma_H0),real(prior_H0 + sigma_H0))
!!$
!!$                         x_old(13) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                         x_old(14) = genunf(real(a_v - sigma_a_v),real(a_v + sigma_a_v))
!!$
!!$                         x_old(15) = genunf(real(a_cal - sigma_a_cal),real(a_cal + sigma_a_cal))
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$                     
!!$                         x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                         x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                         x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                         x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                         x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                         x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                         x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                         x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                         x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                         x_old(10) = genunf(real(prior_mu10 - sigma_mu10),real(prior_mu10 + sigma_mu10))
!!$
!!$                         x_old(11) = genunf(real(prior_Mw - sigma_Mw),real(prior_Mw + sigma_Mw))
!!$
!!$                         x_old(12) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                         x_old(13) = genunf(real(prior_H0 - sigma_H0),real(prior_H0 + sigma_H0))
!!$
!!$                         x_old(14) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                         x_old(15) = genunf(real(a_v - sigma_a_v),real(a_v + sigma_a_v))
!!$
!!$                         x_old(16) = genunf(real(a_cal - sigma_a_cal),real(a_cal + sigma_a_cal))
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$                     
!!$                         x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                         x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                         x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                         x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                         x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                         x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                         x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                         x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                         x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                         x_old(10) = genunf(real(prior_Mw - sigma_Mw),real(prior_Mw + sigma_Mw))
!!$
!!$                         x_old(11) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                         x_old(12) = genunf(real(prior_H0 - sigma_H0),real(prior_H0 + sigma_H0))
!!$
!!$                         x_old(13) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                         x_old(14) = genunf(real(a_v - sigma_a_v),real(a_v + sigma_a_v))
!!$
!!$                         x_old(15) = genunf(real(a_cal - sigma_a_cal),real(a_cal + sigma_a_cal))
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$                     
!!$                         x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                         x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                         x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                         x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                         x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                         x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                         x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                         x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                         x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                         x_old(10) = genunf(real(prior_mu10 - sigma_mu10),real(prior_mu10 + sigma_mu10))
!!$
!!$                         x_old(11) = genunf(real(prior_zpwLMC - sigma_zpwLMC),real(prior_zpwLMC + sigma_zpwLMC))
!!$
!!$                         x_old(12) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                         x_old(13) = genunf(real(prior_H0 - sigma_H0),real(prior_H0 + sigma_H0))
!!$
!!$                         x_old(14) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                         x_old(15) = genunf(real(a_v - sigma_a_v),real(a_v + sigma_a_v))
!!$
!!$                         x_old(16) = genunf(real(a_cal - sigma_a_cal),real(a_cal + sigma_a_cal))
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$
!!$                         x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                         x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                         x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                         x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                         x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                         x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                         x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                         x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                         x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                         x_old(10) = genunf(real(prior_zpw - sigma_zpw),real(prior_zpw + sigma_zpw))
!!$
!!$                         x_old(11) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                         x_old(12) = genunf(real(prior_H0 - sigma_H0),real(prior_H0 + sigma_H0))
!!$
!!$                         x_old(13) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                         x_old(14) = genunf(real(a_v - sigma_a_v),real(a_v + sigma_a_v))
!!$
!!$                      Else
!!$                     
!!$                         x_old(1) = genunf(real(prior_mu1 - sigma_mu1),real(prior_mu1 + sigma_mu1))
!!$
!!$                         x_old(2) = genunf(real(prior_mu2 - sigma_mu2),real(prior_mu2 + sigma_mu2))
!!$
!!$                         x_old(3) = genunf(real(prior_mu3 - sigma_mu3),real(prior_mu3 + sigma_mu3))
!!$
!!$                         x_old(4) = genunf(real(prior_mu4 - sigma_mu4),real(prior_mu4 + sigma_mu4))
!!$
!!$                         x_old(5) = genunf(real(prior_mu5 - sigma_mu5),real(prior_mu5 + sigma_mu5))
!!$
!!$                         x_old(6) = genunf(real(prior_mu6 - sigma_mu6),real(prior_mu6 + sigma_mu6))
!!$
!!$                         x_old(7) = genunf(real(prior_mu7 - sigma_mu7),real(prior_mu7 + sigma_mu7))
!!$
!!$                         x_old(8) = genunf(real(prior_mu8 - sigma_mu8),real(prior_mu8 + sigma_mu8))
!!$
!!$                         x_old(9) = genunf(real(prior_mu9 - sigma_mu9),real(prior_mu9 + sigma_mu9))
!!$
!!$                         x_old(10) = genunf(real(prior_zpw - sigma_zpw),real(prior_zpw + sigma_zpw))
!!$
!!$                         x_old(11) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw))
!!$
!!$                         x_old(12) = genunf(real(prior_H0 - sigma_H0),real(prior_H0 + sigma_H0))
!!$
!!$                         x_old(13) = genunf(real(prior_Zw - sigma_Zw),real(prior_Zw + sigma_Zw))
!!$
!!$                         x_old(14) = genunf(real(a_v - sigma_a_v),real(a_v + sigma_a_v))
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         x_old(1) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(2) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(3) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(4) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(5) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(6) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(7) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(8) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(9) = genunf(real(prior_zpH - sigma_zpH),real(prior_zpH + sigma_zpH))
!!$
!!$                         x_old(10) = genunf(real(prior_bH - sigma_bH),real(prior_bH + sigma_bH))
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$                        
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'
!!$
!!$                   stop
!!$
!!$                End If ! OF ANCHORS
!!$
!!$             End If ! OF ONLY CEPHEIDS
!!$    
!!$          Else
!!$
!!$             x_old(1) = genunf(real(prior_A - sigma_A),real(prior_A + sigma_A))         ! A
!!$             x_old(2) = genunf(real(prior_bw - sigma_bw),real(prior_bw + sigma_bw)) ! bw
!!$             !x_old(3) = genunf(real(-10.d0),real(0.d0)) ! log10(sigma_int)
!!$
!!$          End If
!!$
!!$          If (hyperparameters_as_mcmc) then
!!$             ! SETTING INITIAL POINT FOR HYPER-PARAMETERS
!!$             Do m=number_model_parameters+1,number_of_parameters
!!$                    
!!$                If (using_jeffreys_prior) then
!!$
!!$                   x_old(m) = genunf(real(-1.d0),real(1.d0)) ! log10(\alpha_j)
!!$
!!$                Else
!!$
!!$                   x_old(m) = genunf(real(0.d0),real(1.d0))
!!$
!!$                End If
!!$
!!$             End Do
!!$            
!!$          End If
!!$
!!$          Do m=1,number_of_parameters
!!$
!!$             If (m .gt. number_model_parameters) then
!!$                  
!!$                If (using_jeffreys_prior) then
!!$
!!$                   old_point(m) = 10**(dble(x_old(m)))    !    CONVERT TO \alpha_j
!!$
!!$                Else
!!$
!!$                   old_point(m) = dble(x_old(m)) 
!!$
!!$                End If
!!$
!!$             else
!!$
!!$                old_point(m) = dble(x_old(m))
!!$
!!$             End If
!!$
!!$          End Do
!!$
!!$       End If
!!$
!!$       If (using_hyperparameters) then
!!$          ! OPEN FILE TO STORE MCMC COMPUTATION
!!$          open(UNIT_MCMC_FINAL_FILE,file='./output/chains/mcmc_final_output_HP.txt')
!!$
!!$          write(UNIT_EXE_FILE,*) 'WORKING WITH HYPER-PARAMETERS'
!!$
!!$          write(UNIT_EXE_FILE,*) 'COMPUTING LOG_LIKELIHOOD FOR INITIAL POINT'
!!$            
!!$          ! COMPUTE INITIAL LIKELIHOOD
!!$          If (doing_R11_analysis) then
!!$
!!$             If (include_only_cepheids) then
!!$
!!$                If (all_R11_hosts) then
!!$                    
!!$                   old_loglikelihood = log_R11_likelihood_W_cepheids(old_point(1:number_model_parameters-3),&
!!$                        old_point(number_model_parameters-2),old_point(number_model_parameters-1),&
!!$                        old_point(number_model_parameters),prior_sigma_int)
!!$
!!$                Else
!!$
!!$                   old_loglikelihood = log_likelihood_only_cepheids(galaxy,old_point(1),old_point(2),&
!!$                        old_point(3),prior_sigma_int)
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then !!!HERE
!!$    
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED USING LMC+MW+NGC4258 AS ANCHORS AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$                     
!!$                      Else
!!$                          
!!$                         old_loglikelihood = log_R11_likelihood_W_LMC_MW_NGC4258(old_point(1:number_model_parameters-6),&
!!$                              old_point(number_model_parameters-5),old_point(number_model_parameters-4),&
!!$                              old_point(number_model_parameters-3),old_point(number_model_parameters-2),&
!!$                              old_point(number_model_parameters-1),old_point(number_model_parameters),&
!!$                              prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED USING LMC AND NGC4258 AS ANCHORS AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$                     
!!$                      Else
!!$                          
!!$                         old_loglikelihood = log_R11_likelihood_W_LMC_NGC4258(old_point(1:number_model_parameters-6),&
!!$                              old_point(number_model_parameters-5),old_point(number_model_parameters-4),&
!!$                              old_point(number_model_parameters-3),old_point(number_model_parameters-2),&
!!$                              old_point(number_model_parameters-1),old_point(number_model_parameters),&
!!$                              prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED USING MW+NGC4258 AS ANCHORS AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$                     
!!$                      Else
!!$                          
!!$                         old_loglikelihood = log_R11_likelihood_W_MW_NGC4258(old_point(1:number_model_parameters-6),&
!!$                              old_point(number_model_parameters-5),old_point(number_model_parameters-4),&
!!$                              old_point(number_model_parameters-3),old_point(number_model_parameters-2),&
!!$                              old_point(number_model_parameters-1),old_point(number_model_parameters),&
!!$                              prior_sigma_int,prior_sigma_int_MW)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED USING LMC AND MW AS ANCHORS AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$                     
!!$                      Else
!!$                          
!!$                         old_loglikelihood = log_R11_likelihood_W_LMC_MW(old_point(1:number_model_parameters-6),&
!!$                              old_point(number_model_parameters-5),old_point(number_model_parameters-4),&
!!$                              old_point(number_model_parameters-3),old_point(number_model_parameters-2),&
!!$                              old_point(number_model_parameters-1),old_point(number_model_parameters),&
!!$                              prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                         stop
!!$
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED USING MW AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$                     
!!$                      Else
!!$                          
!!$                         old_loglikelihood = log_R11_likelihood_W_MW(old_point(1:number_model_parameters-6),&
!!$                              old_point(number_model_parameters-5),old_point(number_model_parameters-4),&
!!$                              old_point(number_model_parameters-3),old_point(number_model_parameters-2),&
!!$                              old_point(number_model_parameters-1),old_point(number_model_parameters),&
!!$                              prior_sigma_int,prior_sigma_int_MW)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      End If
!!$
!!$                   End If
!!$                     
!!$                Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         print *,'H BAND NOT IMPLEMENTED USING LMC AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$                     
!!$                      Else
!!$                          
!!$                         old_loglikelihood = log_R11_likelihood_W_LMC(old_point(1:number_model_parameters-6),&
!!$                              old_point(number_model_parameters-5),old_point(number_model_parameters-4),&
!!$                              old_point(number_model_parameters-3),old_point(number_model_parameters-2),&
!!$                              old_point(number_model_parameters-1),old_point(number_model_parameters),&
!!$                              prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$
!!$                         print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                         stop
!!$
!!$                      End If
!!$
!!$                   End If
!!$
!!$                Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                   If (use_metallicity) then 
!!$
!!$                      If (use_H_band) then
!!$
!!$                         old_loglikelihood = log_R11_likelihood_H_NGC4258(old_point(1:number_model_parameters-5),&
!!$                              old_point(number_model_parameters-4),old_point(number_model_parameters-3),&
!!$                              old_point(number_model_parameters-2),old_point(number_model_parameters-1),&
!!$                              old_point(number_model_parameters),prior_sigma_int)
!!$                     
!!$                      Else
!!$
!!$                         old_loglikelihood = log_R11_likelihood_W(old_point(1:number_model_parameters-5),&
!!$                              old_point(number_model_parameters-4),old_point(number_model_parameters-3),&
!!$                              old_point(number_model_parameters-2),old_point(number_model_parameters-1),&
!!$                              old_point(number_model_parameters),prior_sigma_int)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If (use_H_band) then
!!$                     
!!$                         old_loglikelihood = log_R11_likelihood_H(old_point(1:number_model_parameters-1),&
!!$                              old_point(number_model_parameters),prior_sigma_int)
!!$
!!$                      Else
!!$
!!$                         print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                         stop
!!$
!!$                      End If
!!$
!!$                   End If
!!$
!!$                End If ! OF ANCHORS 
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$!             old_loglikelihood = log_Efstathiou_likelihood_hyperparameters(old_point(1),old_point(2),prior_sigma_int)
!!$             old_loglikelihood = log_Efstathiou_likelihood(old_point(1),old_point(2),prior_sigma_int_LMC)
!!$
!!$          End If
!!$
!!$          open(UNIT_PARAMNAMES_FILE,file='./output/chains/mcmc_final_output_HP.paramnames')    !    OPEN FILE WITH PARAMETER NAMES NEEDED BY GETDIST
!!$
!!$          open(UNIT_RANGES_FILE,file='./output/chains/mcmc_final_output_HP.ranges')    !    OPEN FILE WITH HARD BOUNDS NEEDED BY GETDIST
!!$
!!$       Else
!!$          ! OPEN FILE TO STORE MCMC COMPUTATION
!!$          open(UNIT_MCMC_FINAL_FILE,file='./output/chains/mcmc_final_output.txt')
!!$
!!$          write(UNIT_EXE_FILE,*) 'NOT WORKING WITH HYPER-PARAMETERS' 
!!$
!!$          write(UNIT_EXE_FILE,*) 'COMPUTING LOG_LIKELIHOOD FOR INITIAL POINT'
!!$
!!$          ! COMPUTE INITIAL LIKELIHOOD
!!$          old_loglikelihood = log_Efstathiou_likelihood(old_point(1),old_point(2),prior_sigma_int_LMC) 
!!$
!!$          open(UNIT_PARAMNAMES_FILE,file='./output/chains/mcmc_final_output.paramnames')    !    OPEN FILE WITH PARAMETER NAMES NEEDED BY GETDIST
!!$
!!$          open(UNIT_RANGES_FILE,file='./output/chains/mcmc_final_output.ranges')    !    OPEN FILE WITH HARD BOUNDS NEEDED BY GETDIST 
!!$
!!$       End If
!!$
!!$       If (doing_R11_analysis) then
!!$
!!$          If (include_only_cepheids) then
!!$
!!$             If (all_R11_hosts) then
!!$
!!$                paramnames(1) = 'mu01'
!!$                latexname(1) = '\mu_{0,1}'
!!$
!!$                paramnames(2) = 'mu02'
!!$                latexname(2) = '\mu_{0,2}'
!!$
!!$                paramnames(3) = 'mu03'
!!$                latexname(3) = '\mu_{0,3}'
!!$                 
!!$                paramnames(4) = 'mu04'
!!$                latexname(4) = '\mu_{0,4}'
!!$
!!$                paramnames(5) = 'mu05'
!!$                latexname(5) = '\mu_{0,5}'
!!$
!!$                paramnames(6) = 'mu06'
!!$                latexname(6) = '\mu_{0,6}'
!!$
!!$                paramnames(7) = 'mu07'
!!$                latexname(7) = '\mu_{0,7}'
!!$
!!$                paramnames(8) = 'mu08'
!!$                latexname(8) = '\mu_{0,8}'
!!$
!!$                paramnames(9) = 'mu04258'
!!$                latexname(9) = '\mu_{0,4258}'
!!$
!!$                paramnames(10) = 'zpw4258'
!!$                latexname(10) = 'zp_{w,4258}'
!!$
!!$                paramnames(11) = 'bw'
!!$                latexname(11) = 'b_w'
!!$
!!$                paramnames(12) = 'Zw'
!!$                latexname(12) = 'Z_w'
!!$
!!$                Do m=1,number_model_parameters
!!$
!!$                   write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                End Do
!!$
!!$             Else
!!$
!!$                paramnames(1) = 'zpw'
!!$                latexname(1) = 'zp_{w}'
!!$
!!$                paramnames(2) = 'bw'
!!$                latexname(2) = 'b_w'
!!$              
!!$                paramnames(3) = 'Zw'
!!$                latexname(3) = 'Z_w'
!!$
!!$                Do m=1,number_model_parameters
!!$
!!$                   write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                End Do
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$             If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      paramnames(1) = 'mu01'
!!$                      latexname(1) = '\mu_{0,1}'
!!$
!!$                      paramnames(2) = 'mu02'
!!$                      latexname(2) = '\mu_{0,2}'
!!$
!!$                      paramnames(3) = 'mu03'
!!$                      latexname(3) = '\mu_{0,3}'
!!$
!!$                      paramnames(4) = 'mu04'
!!$                      latexname(4) = '\mu_{0,4}'
!!$
!!$                      paramnames(5) = 'mu05'
!!$                      latexname(5) = '\mu_{0,5}'
!!$
!!$                      paramnames(6) = 'mu06'
!!$                      latexname(6) = '\mu_{0,6}'
!!$
!!$                      paramnames(7) = 'mu07'
!!$                      latexname(7) = '\mu_{0,7}'
!!$
!!$                      paramnames(8) = 'mu08'
!!$                      latexname(8) = '\mu_{0,8}'
!!$
!!$                      paramnames(9) = 'mu04258'
!!$                      latexname(9) = '\mu_{0,4258}'
!!$
!!$                      paramnames(10) = 'muLMC'
!!$                      latexname(10) = '\mu_{0,LMC}'
!!$
!!$                      paramnames(11) = 'Mw'
!!$                      latexname(11) = 'M_w'
!!$
!!$                      paramnames(12) = 'bw'
!!$                      latexname(12) = 'b_w'
!!$
!!$                      paramnames(13) = 'H0'
!!$                      latexname(13) = 'H_0'
!!$
!!$                      paramnames(14) = 'Zw'
!!$                      latexname(14) = 'Z_w'
!!$
!!$                      paramnames(15) = 'av'
!!$                      latexname(15) = 'a_v'
!!$
!!$                      paramnames(16) = 'acal'
!!$                      latexname(16) = 'a_{cal}'
!!$
!!$                      Do m=1,number_model_parameters
!!$
!!$                         write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                      End Do
!!$
!!$                   End If
!!$
!!$                Else
!!$                   
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      paramnames(1) = 'mu01'
!!$                      latexname(1) = '\mu_{0,1}'
!!$
!!$                      paramnames(2) = 'mu02'
!!$                      latexname(2) = '\mu_{0,2}'
!!$
!!$                      paramnames(3) = 'mu03'
!!$                      latexname(3) = '\mu_{0,3}'
!!$
!!$                      paramnames(4) = 'mu04'
!!$                      latexname(4) = '\mu_{0,4}'
!!$
!!$                      paramnames(5) = 'mu05'
!!$                      latexname(5) = '\mu_{0,5}'
!!$
!!$                      paramnames(6) = 'mu06'
!!$                      latexname(6) = '\mu_{0,6}'
!!$
!!$                      paramnames(7) = 'mu07'
!!$                      latexname(7) = '\mu_{0,7}'
!!$
!!$                      paramnames(8) = 'mu08'
!!$                      latexname(8) = '\mu_{0,8}'
!!$
!!$                      paramnames(9) = 'mu04258'
!!$                      latexname(9) = '\mu_{0,4258}'
!!$
!!$                      paramnames(10) = 'muLMC'
!!$                      latexname(10) = '\mu_{0,LMC}'
!!$
!!$                      paramnames(11) = 'Mw'
!!$                      latexname(11) = 'M_w'
!!$
!!$                      paramnames(12) = 'bw'
!!$                      latexname(12) = 'b_w'
!!$
!!$                      paramnames(13) = 'H0'
!!$                      latexname(13) = 'H_0'
!!$
!!$                      paramnames(14) = 'Zw'
!!$                      latexname(14) = 'Z_w'
!!$
!!$                      paramnames(15) = 'av'
!!$                      latexname(15) = 'a_v'
!!$
!!$                      paramnames(16) = 'acal'
!!$                      latexname(16) = 'a_{cal}'
!!$
!!$                      Do m=1,number_model_parameters
!!$
!!$                         write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                      End Do
!!$
!!$                   End If
!!$
!!$                Else
!!$                   
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      paramnames(1) = 'mu01'
!!$                      latexname(1) = '\mu_{0,1}'
!!$
!!$                      paramnames(2) = 'mu02'
!!$                      latexname(2) = '\mu_{0,2}'
!!$
!!$                      paramnames(3) = 'mu03'
!!$                      latexname(3) = '\mu_{0,3}'
!!$
!!$                      paramnames(4) = 'mu04'
!!$                      latexname(4) = '\mu_{0,4}'
!!$
!!$                      paramnames(5) = 'mu05'
!!$                      latexname(5) = '\mu_{0,5}'
!!$
!!$                      paramnames(6) = 'mu06'
!!$                      latexname(6) = '\mu_{0,6}'
!!$
!!$                      paramnames(7) = 'mu07'
!!$                      latexname(7) = '\mu_{0,7}'
!!$
!!$                      paramnames(8) = 'mu08'
!!$                      latexname(8) = '\mu_{0,8}'
!!$
!!$                      paramnames(9) = 'mu04258'
!!$                      latexname(9) = '\mu_{0,4258}'
!!$
!!$                      paramnames(10) = 'Mw'
!!$                      latexname(10) = 'M_w'
!!$
!!$                      paramnames(11) = 'bw'
!!$                      latexname(11) = 'b_w'
!!$
!!$                      paramnames(12) = 'H0'
!!$                      latexname(12) = 'H_0'
!!$
!!$                      paramnames(13) = 'Zw'
!!$                      latexname(13) = 'Z_w'
!!$
!!$                      paramnames(14) = 'av'
!!$                      latexname(14) = 'a_v'
!!$
!!$                      paramnames(15) = 'acal'
!!$                      latexname(15) = 'a_{cal}'
!!$
!!$                      Do m=1,number_model_parameters
!!$
!!$                         write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                      End Do
!!$
!!$                   End If
!!$
!!$                Else
!!$                   
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      paramnames(1) = 'mu01'
!!$                      latexname(1) = '\mu_{0,1}'
!!$
!!$                      paramnames(2) = 'mu02'
!!$                      latexname(2) = '\mu_{0,2}'
!!$
!!$                      paramnames(3) = 'mu03'
!!$                      latexname(3) = '\mu_{0,3}'
!!$
!!$                      paramnames(4) = 'mu04'
!!$                      latexname(4) = '\mu_{0,4}'
!!$
!!$                      paramnames(5) = 'mu05'
!!$                      latexname(5) = '\mu_{0,5}'
!!$
!!$                      paramnames(6) = 'mu06'
!!$                      latexname(6) = '\mu_{0,6}'
!!$
!!$                      paramnames(7) = 'mu07'
!!$                      latexname(7) = '\mu_{0,7}'
!!$
!!$                      paramnames(8) = 'mu08'
!!$                      latexname(8) = '\mu_{0,8}'
!!$
!!$                      paramnames(9) = 'mu04258'
!!$                      latexname(9) = '\mu_{0,4258}'
!!$
!!$                      paramnames(10) = 'muLMC'
!!$                      latexname(10) = '\mu_{0,LMC}'
!!$
!!$                      paramnames(11) = 'Mw'
!!$                      latexname(11) = 'M_w'
!!$
!!$                      paramnames(12) = 'bw'
!!$                      latexname(12) = 'b_w'
!!$
!!$                      paramnames(13) = 'H0'
!!$                      latexname(13) = 'H_0'
!!$
!!$                      paramnames(14) = 'Zw'
!!$                      latexname(14) = 'Z_w'
!!$
!!$                      paramnames(15) = 'av'
!!$                      latexname(15) = 'a_v'
!!$
!!$                      paramnames(16) = 'acal'
!!$                      latexname(16) = 'a_{cal}'
!!$
!!$                      Do m=1,number_model_parameters
!!$
!!$                         write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                      End Do
!!$
!!$                   End If
!!$
!!$                Else
!!$                   
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      paramnames(1) = 'mu01'
!!$                      latexname(1) = '\mu_{0,1}'
!!$
!!$                      paramnames(2) = 'mu02'
!!$                      latexname(2) = '\mu_{0,2}'
!!$
!!$                      paramnames(3) = 'mu03'
!!$                      latexname(3) = '\mu_{0,3}'
!!$
!!$                      paramnames(4) = 'mu04'
!!$                      latexname(4) = '\mu_{0,4}'
!!$
!!$                      paramnames(5) = 'mu05'
!!$                      latexname(5) = '\mu_{0,5}'
!!$
!!$                      paramnames(6) = 'mu06'
!!$                      latexname(6) = '\mu_{0,6}'
!!$
!!$                      paramnames(7) = 'mu07'
!!$                      latexname(7) = '\mu_{0,7}'
!!$
!!$                      paramnames(8) = 'mu08'
!!$                      latexname(8) = '\mu_{0,8}'
!!$
!!$                      paramnames(9) = 'mu04258'
!!$                      latexname(9) = '\mu_{0,4258}'
!!$
!!$                      paramnames(10) = 'Mw'
!!$                      latexname(10) = 'M_w'
!!$
!!$                      paramnames(11) = 'bw'
!!$                      latexname(11) = 'b_w'
!!$
!!$                      paramnames(12) = 'H0'
!!$                      latexname(12) = 'H_0'
!!$
!!$                      paramnames(13) = 'Zw'
!!$                      latexname(13) = 'Z_w'
!!$
!!$                      paramnames(14) = 'av'
!!$                      latexname(14) = 'a_v'
!!$
!!$                      paramnames(15) = 'acal'
!!$                      latexname(15) = 'a_{cal}'
!!$
!!$                      Do m=1,number_model_parameters
!!$
!!$                         write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                      End Do
!!$
!!$                   End If
!!$
!!$                Else
!!$                   
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      paramnames(1) = 'mu01'
!!$                      latexname(1) = '\mu_{0,1}'
!!$
!!$                      paramnames(2) = 'mu02'
!!$                      latexname(2) = '\mu_{0,2}'
!!$
!!$                      paramnames(3) = 'mu03'
!!$                      latexname(3) = '\mu_{0,3}'
!!$
!!$                      paramnames(4) = 'mu04'
!!$                      latexname(4) = '\mu_{0,4}'
!!$
!!$                      paramnames(5) = 'mu05'
!!$                      latexname(5) = '\mu_{0,5}'
!!$
!!$                      paramnames(6) = 'mu06'
!!$                      latexname(6) = '\mu_{0,6}'
!!$
!!$                      paramnames(7) = 'mu07'
!!$                      latexname(7) = '\mu_{0,7}'
!!$
!!$                      paramnames(8) = 'mu08'
!!$                      latexname(8) = '\mu_{0,8}'
!!$
!!$                      paramnames(9) = 'mu04258'
!!$                      latexname(9) = '\mu_{0,4258}'
!!$
!!$                      paramnames(10) = 'muLMC'
!!$                      latexname(10) = '\mu_{0,LMC}'
!!$
!!$                      paramnames(11) = 'zpwLMC'
!!$                      latexname(11) = 'zp_{w,LMC}'
!!$
!!$                      paramnames(12) = 'bw'
!!$                      latexname(12) = 'b_w'
!!$
!!$                      paramnames(13) = 'H0'
!!$                      latexname(13) = 'H_0'
!!$
!!$                      paramnames(14) = 'Zw'
!!$                      latexname(14) = 'Z_w'
!!$
!!$                      paramnames(15) = 'av'
!!$                      latexname(15) = 'a_v'
!!$
!!$                      paramnames(16) = 'acal'
!!$                      latexname(16) = 'a_{cal}'
!!$
!!$                      Do m=1,number_model_parameters
!!$
!!$                         write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                      End Do
!!$
!!$                   End If
!!$
!!$                Else
!!$                   
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      paramnames(1) = 'mu01'
!!$                      latexname(1) = '\mu_{0,1}'
!!$
!!$                      paramnames(2) = 'mu02'
!!$                      latexname(2) = '\mu_{0,2}'
!!$
!!$                      paramnames(3) = 'mu03'
!!$                      latexname(3) = '\mu_{0,3}'
!!$
!!$                      paramnames(4) = 'mu04'
!!$                      latexname(4) = '\mu_{0,4}'
!!$
!!$                      paramnames(5) = 'mu05'
!!$                      latexname(5) = '\mu_{0,5}'
!!$
!!$                      paramnames(6) = 'mu06'
!!$                      latexname(6) = '\mu_{0,6}'
!!$
!!$                      paramnames(7) = 'mu07'
!!$                      latexname(7) = '\mu_{0,7}'
!!$
!!$                      paramnames(8) = 'mu08'
!!$                      latexname(8) = '\mu_{0,8}'
!!$
!!$                      paramnames(9) = 'mu04258'
!!$                      latexname(9) = '\mu_{0,4258}'
!!$
!!$                      paramnames(10) = 'zpH4258'
!!$                      latexname(10) = 'zp_{H,4258}'
!!$
!!$                      paramnames(11) = 'bH'
!!$                      latexname(11) = 'b_H'
!!$
!!$                      paramnames(12) = 'H0'
!!$                      latexname(12) = 'H_0'
!!$
!!$                      paramnames(13) = 'ZH'
!!$                      latexname(13) = 'Z_H'
!!$
!!$                      paramnames(14) = 'aH'
!!$                      latexname(14) = 'a_H'
!!$
!!$                      Do m=1,number_model_parameters
!!$
!!$                         write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                      End Do
!!$
!!$                   Else
!!$
!!$                      paramnames(1) = 'mu01'
!!$                      latexname(1) = '\mu_{0,1}'
!!$
!!$                      paramnames(2) = 'mu02'
!!$                      latexname(2) = '\mu_{0,2}'
!!$
!!$                      paramnames(3) = 'mu03'
!!$                      latexname(3) = '\mu_{0,3}'
!!$
!!$                      paramnames(4) = 'mu04'
!!$                      latexname(4) = '\mu_{0,4}'
!!$
!!$                      paramnames(5) = 'mu05'
!!$                      latexname(5) = '\mu_{0,5}'
!!$
!!$                      paramnames(6) = 'mu06'
!!$                      latexname(6) = '\mu_{0,6}'
!!$
!!$                      paramnames(7) = 'mu07'
!!$                      latexname(7) = '\mu_{0,7}'
!!$
!!$                      paramnames(8) = 'mu08'
!!$                      latexname(8) = '\mu_{0,8}'
!!$
!!$                      paramnames(9) = 'mu04258'
!!$                      latexname(9) = '\mu_{0,4258}'
!!$
!!$                      paramnames(10) = 'zpw4258'
!!$                      latexname(10) = 'zp_{w,4258}'
!!$
!!$                      paramnames(11) = 'bw'
!!$                      latexname(11) = 'b_w'
!!$
!!$                      paramnames(12) = 'H0'
!!$                      latexname(12) = 'H_0'
!!$
!!$                      paramnames(13) = 'Zw'
!!$                      latexname(13) = 'Z_w'
!!$
!!$                      paramnames(14) = 'av'
!!$                      latexname(14) = 'a_v'
!!$
!!$                      Do m=1,number_model_parameters
!!$
!!$                         write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$                      End Do
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH1    zp_{H1}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH2    zp_{H2}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH3    zp_{H3}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH4    zp_{H4}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH5    zp_{H5}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH6    zp_{H6}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH7    zp_{H7}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH8    zp_{H8}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'zpH4258    zp_{H4258}'
!!$
!!$                      write(UNIT_PARAMNAMES_FILE,*) 'bH    b_H'
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             End If ! OF ANCHOR
!!$
!!$          End If ! OF ONLY CEPHEIDS
!!$    
!!$       Else
!!$
!!$          paramnames(1) = 'A'
!!$          latexname(1) = 'A'
!!$
!!$          paramnames(2) = 'bw'
!!$          latexname(2) = 'b_w'
!!$
!!$          Do m=1,number_model_parameters
!!$
!!$             write(UNIT_PARAMNAMES_FILE,*) ''//trim(paramnames(m))//'    '//trim(latexname(m))//''
!!$
!!$          End Do
!!$
!!$          !write(16,*) 'sigma_int    \sigma_{int}'
!!$
!!$       End If
!!$
!!$       If (hyperparameters_as_mcmc) then
!!$          ! WRITING PARAMNAMES FOR HYPER-PARAMETERS
!!$          Do m=number_model_parameters+1,number_of_parameters
!!$
!!$             write(string,'(i2.2)') m-number_model_parameters
!!$
!!$             write(UNIT_PARAMNAMES_FILE,*) 'alpha_'//trim(string)//'    \alpha_{'//trim(string)//'}'
!!$
!!$             alpha_string(m-number_model_parameters) = 'alpha_'//trim(string)//'    '
!!$
!!$          End Do
!!$            
!!$       End If

!!$       close(UNIT_PARAMNAMES_FILE)
!!$
!!$       If (doing_R11_analysis) then
!!$
!!$          If (include_only_cepheids) then
!!$              
!!$             If (all_R11_hosts) then
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    25.    34.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.2    -2.5'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -1.    1.'
!!$
!!$             Else
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    25.    34.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    -3.2    -2.5'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    -1.    1.'
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$             If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    0.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(16))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    0.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(16))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'   -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    0.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(16))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'   -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    0.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    15.    25.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(16))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    25.    34.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.2    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    -1.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    0.    1.'
!!$
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    25.    34.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    -1.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    0.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH1    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH2    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH3    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH4    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH5    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH6    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH7    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH8    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH4258    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'bH    -20.    0.'
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             End If ! OF ANCHORS
!!$
!!$          End If ! OF ONLY CEPHEIDS
!!$    
!!$       Else
!!$
!!$          write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    0.    50.'
!!$
!!$          write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    -20.    0.'
!!$
!!$          !    write(17,*) 'sigma_int    1.e-10    1 '
!!$
!!$       End If
!!$
!!$       If (hyperparameters_as_mcmc) then
!!$          ! WRITING HARD BOUNDS FOR HYPER-PARAMETERS
!!$          Do m=number_model_parameters+1,number_of_parameters
!!$
!!$             write(string,'(i2.2)') m-number_model_parameters
!!$                    
!!$             If (using_jeffreys_prior) then
!!$
!!$                print *,'MODIFY THIS PART ACCORDING TO PRIOR'
!!$
!!$                stop
!!$
!!$             Else
!!$
!!$                write(UNIT_RANGES_FILE,*) 'alpha_'//trim(string)//'    0.    1.'
!!$
!!$             End If
!!$
!!$          End Do
!!$            
!!$       End If
!!$
!!$       close(UNIT_RANGES_FILE)
!!$
!!$    End If
!!$
!!$    ! OPEN TEMPORARY FILE TO SAVE CHAIN
!!$    open(UNIT_MCMC_FILE,file='./output/mcmc_output.txt')     
!!$
!!$    write(UNIT_EXE_FILE,*) '# NUMBER OF ITERATIONS IN MCMC : ', number_iterations - steps_taken_before_definite_run
!!$
!!$    If (start_from_fiducial .and. .not.testing_Gaussian_likelihood) then
!!$
!!$        write(UNIT_EXE_FILE,*) '# FIDUCIAL MODEL IS (PARAMETERS ARE ORDERED AS IN CHAINS FILES) :', old_point
!!$
!!$        write(UNIT_EXE_FILE,'(a37,es18.10)') '# ln(L/L_max) AT THE FIDUCIAL MODEL :', old_loglikelihood
!!$
!!$    End If
!!$
!!$    If (hyperparameters_as_mcmc) then
!!$
!!$        write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    A    bw   ', alpha_string(1:number_hyperparameters)
!!$ 
!!$    Else
!!$
!!$       If (doing_R11_analysis) then
!!$
!!$          If (include_only_cepheids) then
!!$
!!$             write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$          Else
!!$
!!$             If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+NGC4258+MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+NGC4258+MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    zpH1    zpH2    zpH3'//trim(&
!!$                           '    zpH4    zpH5    zpH6    zpH7    zpH8    zpH4258    bH')//'' 
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             End If ! OF ANCHORS
!!$
!!$          End If ! ONLY CEPHEIDS
!!$    
!!$       Else
!!$
!!$          !write(13,*) '# Weight   -ln(L/L_{max})    A    bw    sigma_int ' 
!!$          write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})   ', paramnames(1:number_model_parameters) 
!!$
!!$       End If
!!$
!!$    End If
!!$
!!$    write(UNIT_EXE_FILE,*)'STARTING SAMPLING OF PARAMETER SPACE'
!!$
    ! LOOP TO EXPLORE PARAMETER SPACE STARTS HERE
!!$    Do m=1,number_iterations
!!$
!!$       ! GENERATE NEW POINT IN PARAMETER SPACE FROM MULTIVARIATE DISTRIBUTION; CODE USES RANLIB LIBRARY (BE CAREFUL WITH X_OLD AND OLD_POINT DEFINITIONS)
!!$       If (testing_Gaussian_likelihood) then
!!$
!!$          call setgmn(x_old,real(Covgauss),number_of_parameters,parm) 
!!$ 
!!$          call genmn(parm,x_new,work)
!!$
!!$       Else
!!$          
!!$          call setgmn(x_old,real(Covguess),number_of_parameters,parm) 
!!$
!!$          call genmn(parm,x_new,work)
!!$
!!$       End If
!!$
!!$       If (doing_R11_analysis) then
!!$
!!$          If (include_only_cepheids) then
!!$
!!$             If (all_R11_hosts) then
!!$
!!$                plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$
!!$                plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                plausibility(10) =  (x_new(10) .le. real(25.d0)) .or. (x_new(10) .ge. real(34.d0)) 
!!$
!!$                plausibility(11) =  (x_new(11) .le. real(-3.2d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                plausibility(12) =  (x_new(12) .le. real(-1.d0)) .or. (x_new(12) .ge. real(1.d0)) 
!!$
!!$             Else
!!$
!!$                plausibility(1) =  (x_new(1) .le. real(25.d0)) .or. (x_new(1) .ge. real(34.d0)) 
!!$
!!$                plausibility(2) =  (x_new(2) .le. real(-3.2d0)) .or. (x_new(2) .ge. real(-2.5d0)) 
!!$
!!$                plausibility(3) =  (x_new(3) .le. real(-1.d0)) .or. (x_new(3) .ge. real(1.d0)) 
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$             If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(40.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-7.d0)) .or. (x_new(11) .ge. real(-4.d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(-3.5d0)) .or. (x_new(12) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(55.d0)) .or. (x_new(13) .ge. real(95.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(-2.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(0.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                      plausibility(16) =  (x_new(16) .le. real(-1.d0)) .or. (x_new(16) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(40.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-7.d0)) .or. (x_new(11) .ge. real(-4.d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(-3.5d0)) .or. (x_new(12) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(55.d0)) .or. (x_new(13) .ge. real(95.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(-2.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(0.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                      plausibility(16) =  (x_new(16) .le. real(-1.d0)) .or. (x_new(16) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(-7.d0)) .or. (x_new(10) .ge. real(-4.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-3.5d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(55.d0)) .or. (x_new(12) .ge. real(95.d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(-2.d0)) .or. (x_new(13) .ge. real(1.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(0.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(-1.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(40.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-7.d0)) .or. (x_new(11) .ge. real(-4.d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(-3.5d0)) .or. (x_new(12) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(55.d0)) .or. (x_new(13) .ge. real(95.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(-2.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(0.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                      plausibility(16) =  (x_new(16) .le. real(-1.d0)) .or. (x_new(16) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(-7.d0)) .or. (x_new(10) .ge. real(-4.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-3.5d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(55.d0)) .or. (x_new(12) .ge. real(95.d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(-2.d0)) .or. (x_new(13) .ge. real(1.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(0.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(-1.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(40.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(15.d0)) .or. (x_new(11) .ge. real(25.d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(-3.5d0)) .or. (x_new(12) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(55.d0)) .or. (x_new(13) .ge. real(95.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(-2.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(0.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                      plausibility(16) =  (x_new(16) .le. real(-1.d0)) .or. (x_new(16) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$                           
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) =  (x_new(10) .le. real(25.d0)) .or. (x_new(10) .ge. real(34.d0)) 
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-3.2d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(55.d0)) .or. (x_new(12) .ge. real(95.d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(-1.d0)) .or. (x_new(13) .ge. real(1.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(0.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$                           
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) =  (x_new(10) .le. real(25.d0)) .or. (x_new(10) .ge. real(34.d0)) 
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-3.5d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(55.d0)) .or. (x_new(12) .ge. real(95.d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(-1.d0)) .or. (x_new(13) .ge. real(1.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(0.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))
!!$
!!$                      plausibility(2) = (x_new(2) .le. real(0.d0)) .or. (x_new(2) .ge. real(5.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(0.d0)) .or. (x_new(3) .ge. real(5.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(0.d0)) .or. (x_new(4) .ge. real(5.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(0.d0)) .or. (x_new(5) .ge. real(5.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(0.d0)) .or. (x_new(6) .ge. real(5.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(0.d0)) .or. (x_new(7) .ge. real(5.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(0.d0)) .or. (x_new(8) .ge. real(5.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(0.d0)) .or. (x_new(9) .ge. real(50.d0))
!!$
!!$                      plausibility(10) =  (x_new(10) .le. real(-20.d0)) .or. (x_new(10) .ge. real(0.d0)) 
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'
!!$
!!$                stop
!!$
!!$             End If ! OF ANCHOR
!!$
!!$          End If ! OF CEPHEIDS
!!$    
!!$       Else
!!$
!!$          plausibility(1) = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))
!!$          plausibility(2) = (x_new(2) .le. real(-2.d1)) .or. (x_new(2) .ge. real(0.d0))
!!$          !plausibility(3) =  (x_new(3) .gt. real(0.d0)) .or. (x_new(3) .lt. real(-10.d0))    ! limit log10(sigma_int)
!!$
!!$       End If
!!$
!!$       If (hyperparameters_as_mcmc) then
!!$          ! CHECKING PLAUSIBILITY FOR HYPER-PARAMETERS
!!$          Do n=number_model_parameters+1,number_of_parameters
!!$                   
!!$             plausibility(n) =  (x_new(n) .gt. real(1.d0)) .or. (x_new(n) .lt. real(0.d0))    ! limit alpha_j
!!$
!!$          End Do
!!$            
!!$       End If
!!$
!!$       Do n=1,number_of_parameters
!!$
!!$          If (plausibility(n)) then
!!$
!!$             non_plausible_parameters = .true.
!!$
!!$             exit
!!$
!!$          Else 
!!$
!!$             non_plausible_parameters = .false.
!!$
!!$          End If
!!$
!!$       End Do
!!$
!!$       ! NEW POINT GENERATED 
!!$
!!$       Do n=1,number_of_parameters
!!$
!!$          If (n .gt. number_model_parameters) then
!!$
!!$             If (using_jeffreys_prior) then
!!$
!!$                current_point(n) = 10**(dble(x_new(n))) ! CONVERTING LOG10(alpha_j) to alpha_j 
!!$
!!$             Else
!!$
!!$                current_point(n) = dble(x_new(n))
!!$              
!!$             End If
!!$
!!$          Else
!!$
!!$             current_point(n) = dble(x_new(n))
!!$
!!$          End If
!!$
!!$       End Do
!!$
!!$       ! EVALUATE LOG_LIKELIHOOD FOR CURRENT POINT IN PARAMETER SPACE
!!$       If (testing_Gaussian_likelihood) then
!!$
!!$          current_loglikelihood = log_Gaussian_likelihood(current_point)
!!$
!!$       Else
!!$
!!$          If (non_plausible_parameters) then
!!$
!!$             current_loglikelihood = -1.d10
!!$
!!$          Else
!!$
!!$             If (using_hyperparameters) then    
!!$
!!$                If (doing_R11_analysis) then
!!$
!!$                   If (include_only_cepheids) then
!!$
!!$                      If (all_R11_hosts) then
!!$                        
!!$                         current_loglikelihood = log_R11_likelihood_W_cepheids(current_point(1:number_model_parameters-3),&
!!$                              current_point(number_model_parameters-2),current_point(number_model_parameters-1),&
!!$                              current_point(number_model_parameters),prior_sigma_int)
!!$
!!$                      Else
!!$
!!$                         current_loglikelihood = log_likelihood_only_cepheids(galaxy,current_point(1),&
!!$                              current_point(2),current_point(3),prior_sigma_int)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then !!!HERE
!!$    
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING LMC+MW+NGC4258 AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_LMC_MW_NGC4258(current_point(1:&
!!$                                    number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING LMC AND NGC4258 AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_LMC_NGC4258(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING MW AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_MW_NGC4258(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_MW)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING LMC AND MW AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_LMC_MW(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING MW AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_MW(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_MW)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING LMC AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_LMC(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               current_loglikelihood = log_R11_likelihood_H_NGC4258(current_point(1:number_model_parameters-5),&
!!$                                    current_point(number_model_parameters-4),current_point(number_model_parameters-3),&
!!$                                    current_point(number_model_parameters-2),current_point(number_model_parameters-1),&
!!$                                    current_point(number_model_parameters),prior_sigma_int)
!!$
!!$                            Else
!!$
!!$                               current_loglikelihood = log_R11_likelihood_W(current_point(1:number_model_parameters-5),&
!!$                                    current_point(number_model_parameters-4),current_point(number_model_parameters-3),&
!!$                                    current_point(number_model_parameters-2),current_point(number_model_parameters-1),&
!!$                                    current_point(number_model_parameters),prior_sigma_int)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$
!!$                               current_loglikelihood = log_R11_likelihood_H(current_point(1:number_model_parameters-1),&
!!$                                    current_point(number_model_parameters),prior_sigma_int)
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      End If ! OF ANCHORS 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$!                   current_loglikelihood = log_Efstathiou_likelihood_hyperparameters(current_point(1),&
!!$ !                       current_point(2),prior_sigma_int)
!!$
!!$                   current_loglikelihood = log_Efstathiou_likelihood(current_point(1),current_point(2),prior_sigma_int_LMC)
!!$              
!!$                End If
!!$
!!$             Else
!!$
!!$                current_loglikelihood = log_Efstathiou_likelihood(current_point(1),current_point(2),prior_sigma_int_LMC)
!!$
!!$             End If
!!$
!!$          End If
!!$
!!$       End If
!!$       ! LOG_LIKELIHOOD FOR CURRENT POINT COMPUTED
!!$
!!$       !MAKE DECISION ABOUT CURRENT POINT : ACCEPT OR REJECT IT
!!$       If (current_loglikelihood .ge. old_loglikelihood) then ! ACCEPT CURRENT POINT
!!$
!!$          If (m .gt. steps_taken_before_definite_run) then
!!$
!!$             number_accepted_points = number_accepted_points + 1         
!!$
!!$          End If
!!$
!!$          ! COMPUTING ACCEPTANCE PROBABILITY FOR CURRENT POINT
!!$          acceptance_probability(m) = min(1.d0,exp(current_loglikelihood - old_loglikelihood))    
!!$        
!!$          If (m .le. steps_taken_before_definite_run) then ! WRITE OUT INFORMATION IN TEMPORARY FILE
!!$               
!!$             write(UNIT_MCMC_FILE,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)
!!$
!!$          Else ! WRITE OUT INFORMATION IN DEFINITE FILE
!!$
!!$             write(UNIT_MCMC_FINAL_FILE,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)
!!$
!!$          End If
!!$       
!!$          weight = 1    
!!$
!!$          old_loglikelihood = current_loglikelihood
!!$        
!!$          Do i=1,number_of_parameters 
!!$
!!$             old_point(i) = current_point(i)
!!$
!!$             If (i .gt. number_model_parameters) then
!!$
!!$                If (using_jeffreys_prior) then
!!$
!!$                   x_old(i) = log10(real(old_point(i))) ! CONVERTING alpha_j TO  log10(alpha_j)
!!$
!!$                Else
!!$
!!$                   x_old(i) = real(old_point(i))
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                x_old(i) = real(old_point(i))
!!$
!!$             End If
!!$
!!$          End Do
!!$   
!!$       Else ! ACCEPT OR REJECT THE CURRENT POINT ACCORDING TO :
!!$
!!$          random_uniform = dble(genunf(real(0.),real(1.)))
!!$
!!$          If ( random_uniform .le. exp(current_loglikelihood-old_loglikelihood)) then ! ACCEPT CURRENT POINT
!!$
!!$             If (m .gt. steps_taken_before_definite_run) then
!!$
!!$                number_accepted_points = number_accepted_points + 1         
!!$
!!$             End If
!!$                
!!$             acceptance_probability(m) = min(1.d0,dexp(current_loglikelihood - old_loglikelihood))    
!!$
!!$             If (m .le. steps_taken_before_definite_run) then ! WRITE OUT INFORMATION TO TEMPORARY FILE
!!$
!!$                write(UNIT_MCMC_FILE,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)
!!$
!!$             Else ! WRITE OUT INFORMATION TO DEFINITE FILE
!!$                   
!!$                write(UNIT_MCMC_FINAL_FILE,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)
!!$
!!$             End If
!!$
!!$             weight = 1
!!$
!!$             old_loglikelihood = current_loglikelihood
!!$
!!$             Do i=1,number_of_parameters 
!!$
!!$                old_point(i) = current_point(i)
!!$
!!$                If (i .gt. number_model_parameters) then
!!$                        
!!$                   If (using_jeffreys_prior) then
!!$
!!$                      x_old(i) = real(log10(old_point(i))) ! CONVERTING alpha_j TO log10(alpha_j)
!!$
!!$                   Else
!!$
!!$                      x_old(i) = real(old_point(i))
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   x_old(i) = real(old_point(i))
!!$
!!$                End If
!!$
!!$             End Do
!!$
!!$          Else   ! REJECT CURRENT POINT 
!!$
!!$             If (m .gt. steps_taken_before_definite_run) then
!!$
!!$                number_rejected_points = number_rejected_points + 1            
!!$
!!$             End If
!!$
!!$             acceptance_probability(m) = min(1.d0,exp(current_loglikelihood - old_loglikelihood))    
!!$
!!$             weight = weight + 1
!!$
!!$             Do i=1,number_of_parameters 
!!$
!!$                If (i .gt. number_model_parameters) then
!!$
!!$                   If (using_jeffreys_prior) then
!!$
!!$                      x_old(i) = real(log10(old_point(i))) ! CONVERT alpha_j TO log10(alpha_j)
!!$
!!$                   Else
!!$
!!$                      x_old(i) = real(old_point(i))
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   x_old(i) = real(old_point(i))
!!$
!!$                End If
!!$
!!$             End Do
!!$
!!$          End If
!!$
!!$       End If
!!$       ! DECISION ABOUT CURRENT POINT MADE
!!$
!!$       ! COMPUTING AVERAGE ACCEPTANCE PROBABILITY AND UPDATING BOTH COVARIANCE MATRIX AND JUMPING FACTOR (IF NEEDED)
!!$       If ((mod(m,jumping_factor_update) .eq. 0) .and. (m .le. steps_taken_before_definite_run) ) then
!!$
!!$          average_acceptance_probability = sum(acceptance_probability(m-jumping_factor_update+1:m))&
!!$               /real(jumping_factor_update)
!!$
!!$          !            write(15,*) 'CURRENT AVERAGE ACCEPTANCE PROBABILITY = ',average_acceptance_probability
!!$            
!!$          ! UPDATE JUMPING FACTOR IF NEEDED        
!!$          If (average_acceptance_probability .lt. 0.1) then 
!!$               
!!$             jumping_factor = (1.d0 - step_size_changes)    !    Decreasing step size
!!$
!!$             If (testing_Gaussian_likelihood) then
!!$
!!$                Covgauss = jumping_factor*Covgauss
!!$
!!$             Else
!!$
!!$                Covguess = jumping_factor*Covguess
!!$
!!$             End If
!!$
!!$          Else if (average_acceptance_probability .gt. 0.4) then
!!$
!!$             jumping_factor = (1.d0 + step_size_changes)    !    Increasing step size 
!!$
!!$             If (testing_Gaussian_likelihood) then
!!$
!!$                Covgauss = jumping_factor*Covgauss
!!$
!!$             Else
!!$
!!$                Covguess = jumping_factor*Covguess
!!$
!!$             End If
!!$
!!$          End If
!!$          ! JUMPING FACTOR UPDATED (IF IT WAS NEEDED)
!!$             
!!$          not_good_aap = (average_acceptance_probability .lt. 0.1) .or. (average_acceptance_probability .gt. 0.4)
!!$
!!$          If ( (mod(m,covariance_matrix_update) .eq. 0) .and. not_good_aap) then
!!$                
!!$             call stat('./output/mcmc_output.txt',buff,status1)
!!$
!!$             If ((status1 .eq. 0) .and. (buff(8) .gt. 0)) then
!!$              
!!$                If (testing_Gaussian_likelihood) then
!!$
!!$                   call system('cd output; python compute_covariance_matrix_Gaussian.py')
!!$
!!$                   call read_covariance_matrix_mcmc(Covgauss)
!!$
!!$                   close(UNIT_MCMC_FILE)
!!$
!!$                   call system('rm ./output/mcmc_output.txt')
!!$
!!$                   open(UNIT_MCMC_FILE,file='./output/mcmc_output.txt')
!!$
!!$                Else
!!$
!!$                   If (using_hyperparameters) then
!!$
!!$                      If (hyperparameters_as_mcmc) then
!!$                                
!!$                         call system('cd output; python compute_covariance_matrix_HP.py')
!!$                                
!!$                      Else
!!$
!!$                         call system('cd output; python compute_covariance_matrix.py')
!!$                                
!!$                      End If
!!$
!!$                   Else
!!$                           
!!$                      call system('cd output; python compute_covariance_matrix.py')
!!$
!!$                   End If
!!$
!!$                   call read_covariance_matrix_mcmc(Covguess)
!!$
!!$                   close(UNIT_MCMC_FILE)
!!$
!!$                   call system('rm ./output/mcmc_output.txt')
!!$
!!$                   open(UNIT_MCMC_FILE,file='./output/mcmc_output.txt')
!!$
!!$                End If
!!$
!!$             End If
!!$
!!$          End If
!!$
!!$       End If
!!$
!!$    End Do
!!$    ! LOOP TO EXPLORE PARAMETER SPACE ENDED
!!$
!!$    write(UNIT_EXE_FILE,*) 'NUMBER OF REJECTED POINTS = ', number_rejected_points
!!$
!!$    write(UNIT_EXE_FILE,*) 'ACCEPTANCE RATIO = ', dble(number_iterations - steps_taken_before_definite_run&
!!$    - number_rejected_points)/dble(number_iterations - steps_taken_before_definite_run)
!!$ 
!!$    ! CLOSE FILE STORING CHAIN
!!$    close(UNIT_MCMC_FINAL_FILE)
!!$    ! CLOSE TEMPORARY FILE FOR CHAINS
!!$    close(UNIT_MCMC_FILE)
!!$
!!$    !ANALYZE SAMPLES, MAKE FIGURES, COMPUTE BESTFIT AND HYPER-PARAMETERS (IF NEEDED)
!!$    If (testing_Gaussian_likelihood) then
!!$
!!$        call system('cd analyzer; python analyze.py')
!!$
!!$    Else
!!$
!!$        If (using_hyperparameters) then
!!$
!!$            If (hyperparameters_as_mcmc) then
!!$
!!$                call system('cd analyzer; python analyze_HP_as_MCMC.py')
!!$
!!$            Else
!!$               
!!$               If (doing_R11_analysis) then  !MUST IMPLEMENT OTHER OPTIONS LATER!!!!!!!!!!!!!!!!!
!!$
!!$                  If (include_only_cepheids) then
!!$
!!$                     call system('cd analyzer; python analyze_HP.py')
!!$
!!$                  Else
!!$
!!$                     If (use_H_band) then
!!$
!!$                        call system('cd analyzer; python analyze_HP_R11_H.py')
!!$
!!$                     Else
!!$
!!$                        call system('cd analyzer; python analyze_HP_R11_W.py')
!!$
!!$                     End If
!!$
!!$                  End If
!!$
!!$               Else
!!$
!!$                  call system('cd analyzer; python analyze_HP.py')
!!$
!!$               End If
!!$
!!$            End If
!!$
!!$        Else
!!$
!!$            call system('cd analyzer; python analyze.py')
!!$
!!$        End If    
!!$
!!$    End If
!!$    
!!$    call read_bestfit_mcmc(bestfit)
!!$
!!$    call read_means_mcmc(means)
!!$
!!$    If (doing_R11_analysis) then
!!$
!!$       If (include_only_cepheids) then 
!!$
!!$          write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$          Do m=1,number_model_parameters
!!$
!!$             write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$          End Do
!!$
!!$       Else
!!$
!!$          If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$             write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$             Do m=1,number_model_parameters
!!$
!!$                write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$             End Do
!!$
!!$             chi2SNIabestfit = 0.d0
!!$
!!$             Do n=1,number_of_hosts_galaxies-1 ! SN Ia
!!$
!!$                If (use_HP_in_SNIa) then
!!$                   
!!$                   write(UNIT_EXE_FILE,*) 'NEED TO IMPLEMENT CHI2 VALUE WHEN HPs IN SNIa'
!!$
!!$                   stop
!!$!                   chi2SNIabestfit = log(new_chi2(chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n))) + &
!!$ !                       log(N_tilde_R11_SNIa(n)) + chi2SNIabestfit
!!$
!!$                Else
!!$
!!$                   chi2SNIabestfit = chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n) + chi2SNIabestfit
!!$
!!$                End If
!!$              
!!$             End Do
!!$
!!$             write(UNIT_EXE_FILE,*) 'CHI2 CONTRIBUTION OF SNIa IS', chi2SNIabestfit
!!$
!!$             print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$             Do m=1,number_model_parameters
!!$
!!$                write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$             End Do
!!$
!!$             chi2SNIabestfit = 0.d0
!!$
!!$             Do n=1,number_of_hosts_galaxies-1 ! SN Ia
!!$
!!$                If (use_HP_in_SNIa) then
!!$                   
!!$                   write(UNIT_EXE_FILE,*) 'NEED TO IMPLEMENT CHI2 VALUE WHEN HPs IN SNIa'
!!$
!!$                   stop
!!$!                   chi2SNIabestfit = log(new_chi2(chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n))) + &
!!$ !                       log(N_tilde_R11_SNIa(n)) + chi2SNIabestfit
!!$
!!$                Else
!!$
!!$                   chi2SNIabestfit = chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n) + chi2SNIabestfit
!!$
!!$                End If
!!$              
!!$             End Do
!!$
!!$             write(UNIT_EXE_FILE,*) 'CHI2 CONTRIBUTION OF SNIa IS', chi2SNIabestfit
!!$
!!$             print *,'NGC4258+LMC NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'NGC4258+MW NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$          
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'MW+LMC NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'MW NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$
!!$             write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$             Do m=1,number_model_parameters
!!$
!!$                write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$             End Do
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             If (use_metallicity) then 
!!$
!!$                If (use_H_band) then
!!$
!!$                   write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$                   End Do
!!$
!!$                Else
!!$
!!$                   write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$                   End Do
!!$
!!$                End If
!!$
!!$                chi2SNIabestfit = 0.d0
!!$
!!$                Do n=1,number_of_hosts_galaxies-1 ! SN Ia
!!$
!!$                   If (use_HP_in_SNIa) then
!!$                   
!!$                      write(UNIT_EXE_FILE,*) 'NEED TO IMPLEMENT CHI2 VALUE WHEN HPs IN SNIa'
!!$
!!$                      continue
!!$                      !                   chi2SNIabestfit = log(new_chi2(chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n))) + &
!!$                      !                       log(N_tilde_R11_SNIa(n)) + chi2SNIabestfit
!!$
!!$                   Else
!!$
!!$                      chi2SNIabestfit = chi2R11_SNIa(bestfit(n),bestfit(12),bestfit(14),n) + chi2SNIabestfit
!!$
!!$                      print *, chi2R11_SNIa(bestfit(n),bestfit(12),bestfit(14),n)
!!$
!!$                   End If
!!$              
!!$                End Do
!!$
!!$                write(UNIT_EXE_FILE,*) 'CHI2 CONTRIBUTION OF SNIa IS', chi2SNIabestfit
!!$
!!$             Else
!!$
!!$                If (use_H_band) then
!!$
!!$                   write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$                   write(UNIT_EXE_FILE,*) 'bH = ', bestfit(number_of_parameters)
!!$
!!$                Else
!!$
!!$                   print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                   stop
!!$                        
!!$                End If
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$             print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'
!!$
!!$             stop
!!$
!!$          End If
!!$    
!!$       End If
!!$    
!!$    Else
!!$
!!$       write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$       Do m=1,number_model_parameters
!!$
!!$          write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$       End Do
!!$
!!$       !write(15,*) 'sigma_int = ', bestfit(3)
!!$
!!$    End If
!!$
!!$    If (hyperparameters_as_mcmc .and. using_hyperparameters) then
!!$    ! WRITING BESTFIT FOR HYPER-PARAMETERS
!!$        Do m=number_model_parameters+1,number_of_parameters
!!$
!!$            write(string,'(i2)') m-number_model_parameters
!!$
!!$            write(UNIT_EXE_FILE,*) 'alpha_'//trim(string)//' = ', bestfit(m)
!!$
!!$        End Do
!!$            
!!$    End If
!!$
!!$    If (doing_R11_analysis) then
!!$
!!$       If (include_only_cepheids) then
!!$
!!$          write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$          Do m=1,number_model_parameters
!!$
!!$             write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$          End Do
!!$
!!$       Else
!!$
!!$          If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$             print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             print *,'NGC4258+LMC NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'NGC4258+MW NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$          
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'MW+LMC NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'MW NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             If (use_metallicity) then 
!!$
!!$                If (use_H_band) then
!!$
!!$                   print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                   stop
!!$                
!!$                Else
!!$
!!$                   write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$                   End Do
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                If (use_H_band) then
!!$
!!$                   print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                   stop
!!$
!!$                Else
!!$
!!$                   print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                   stop
!!$                        
!!$                End If
!!$
!!$             End If
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             If (use_metallicity) then 
!!$
!!$                If (use_H_band) then
!!$
!!$                   write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$                   End Do
!!$
!!$                Else
!!$
!!$                   write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$                   End Do
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                If (use_H_band) then
!!$
!!$                   write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$                   write(UNIT_EXE_FILE,*) 'bH = ', means(number_of_parameters)
!!$
!!$                Else
!!$
!!$                   print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                   stop
!!$                        
!!$                End If
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$             print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'
!!$
!!$             stop
!!$
!!$          End If
!!$
!!$       End If
!!$    
!!$    Else
!!$
!!$       write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$       Do m=1,number_model_parameters
!!$
!!$          write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$       End Do
!!$
!!$       !write(15,*) 'sigma_int = ', means(3)
!!$
!!$    End If
!!$
!!$    If (hyperparameters_as_mcmc .and. using_hyperparameters) then
!!$    ! WRITING SAMPLE MEANS FOR HYPER-PARAMETERS
!!$        Do m=number_model_parameters+1,number_of_parameters
!!$
!!$            write(string,'(i2)') m-number_model_parameters
!!$
!!$            write(UNIT_EXE_FILE,*) 'alpha_'//trim(string)//' = ', means(m)
!!$
!!$        End Do
!!$            
!!$    End If
!!$
!!$    If (using_hyperparameters .and. .not.hyperparameters_as_mcmc) then
!!$
!!$        write(UNIT_EXE_FILE,*) 'COMPUTING EFFECTIVE HYPER-PARAMETERS'
!!$
!!$        If (doing_R11_analysis) then
!!$
!!$           If (include_only_cepheids) then
!!$
!!$              print *, 'MUST IMPLEMENT EFFECTIVE HYPER-PARAMETERS WHEN INCLUDING ONLY CEPHEIDS'
!!$
!!$           Else
!!$
!!$              If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                 print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'
!!$
!!$                 stop
!!$
!!$              Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                 print *,'NGC4258+LMC NOT IMPLEMENTED YET'
!!$
!!$                 stop
!!$
!!$              Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                 print *,'NGC4258+MW NOT IMPLEMENTED YET'
!!$
!!$                 stop
!!$          
!!$              Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                 print *,'MW+LMC NOT IMPLEMENTED YET'
!!$
!!$                 stop
!!$
!!$              Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                 print *,'MW NOT IMPLEMENTED YET'
!!$
!!$                 stop
!!$
!!$              Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                 If (use_HP_in_SNIa) then
!!$
!!$                    open(UNIT_HP_FILE,file='./output/chains/effective_hyperparameters_SNIa.txt')
!!$
!!$                    write(UNIT_HP_FILE,*) 'n4258', bestfit(9)-bestfit(9),&
!!$                         bestfit(9)+5.d0*log10(bestfit(13))-25.d0,0.d0,&
!!$                         bestfit(9)+5.d0*log10(bestfit(13))-25.d0,1.d0
!!$
!!$                    Do m=1,number_of_hosts_galaxies-1
!!$
!!$                       If ( chi2R11_SNIa(bestfit(m),bestfit(13),bestfit(15),m) .le. 1.d0) then
!!$
!!$                          write(UNIT_HP_FILE,*) Fieldmvi(m), bestfit(m)-bestfit(9),mvi5av(m), Sigma_mvi5av(m), &
!!$                               bestfit(m)+5.d0*log10(bestfit(13))-25.d0,1.d0
!!$
!!$                       Else
!!$
!!$                          write(UNIT_HP_FILE,*) Fieldmvi(m), bestfit(m)-bestfit(9),mvi5av(m), Sigma_mvi5av(m), &
!!$                               bestfit(m)+5.d0*log10(bestfit(13))-25.d0,&
!!$                               1.d0/chi2R11_SNIa(bestfit(m),bestfit(13),bestfit(15),m)
!!$
!!$                       End If
!!$
!!$                    End Do
!!$
!!$                    close(UNIT_HP_FILE)
!!$
!!$                 Else
!!$
!!$                    continue
!!$
!!$                 End If
!!$
!!$                 If (use_HP_in_anchor) then
!!$
!!$                    If ( chi2R11_anchor_LMC(bestfit(10)) .le. 1.d0) then
!!$
!!$                       write(UNIT_EXE_FILE,*) 'HP FOR ANCHOR LMC IS: ',1.d0
!!$
!!$                    Else
!!$
!!$                       write(UNIT_EXE_FILE,*) 'HP FOR ANCHOR LMC IS: ',1.d0/chi2R11_anchor_LMC(bestfit(10))
!!$
!!$                    End If
!!$
!!$                 Else
!!$
!!$                    continue
!!$
!!$                 End If
!!$
!!$              Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                 If (use_metallicity) then 
!!$
!!$                    If (use_H_band) then
!!$
!!$                       If (use_HP_per_host) then
!!$
!!$                          print *, 'HPs PER HOST NOT IMPLEMENTED WHEN USING H BAND'
!!$                          
!!$                          stop
!!$
!!$!                          open(20,file='./output/effective_hyperparameters_hosts.txt')
!!$
!!$!                          Do n=1,number_of_hosts_galaxies
!!$ 
!!$!                             write(20,*) n, 1.d0/chi2R11_W_host(bestfit(n),bestfit(9),bestfit(10),bestfit(11),&
!!$!                                        bestfit(13),prior_sigma_int,n), host(n)
!!$
!!$!                          End Do
!!$
!!$!                          close(20)
!!$
!!$!                          call system('cd analyzer; python plot_HP_hosts.py')
!!$ 
!!$                       Else If (use_HP_per_cepheid) then
!!$
!!$                          open(UNIT_HP_FILE,file='./output/chains/effective_hyperparameters_cepheids.txt')
!!$
!!$                          Do m=1,size(Field)
!!$
!!$                             Do n=1,number_of_hosts_galaxies
!!$ 
!!$                                If (host(n) .eq. Field(m)) then
!!$                                   
!!$                                   If (PeriodR11(m) .lt. cepheid_Period_limit) then
!!$
!!$                                      If ( chi2R11_H_2(bestfit(n),bestfit(9),bestfit(10),bestfit(11),bestfit(13),prior_sigma_int,m)&
!!$                                        .le. 1.d0 ) then
!!$
!!$                                         write(UNIT_HP_FILE,*) PeriodR11(m), F160WR11(m) - &
!!$                                              P_L_relation_passband_H_2(bestfit(n),bestfit(9),bestfit(10),bestfit(11),bestfit(13),&
!!$                                              OHR11(m),PeriodR11(m)),eF160WR11(m), 1.d0, Field(m)
!!$
!!$                                      Else
!!$
!!$                                         write(UNIT_HP_FILE,*) PeriodR11(m), F160WR11(m) - &
!!$                                              P_L_relation_passband_H_2(bestfit(n),bestfit(9),bestfit(10),bestfit(11),bestfit(13),&
!!$                                              OHR11(m),PeriodR11(m)), eF160WR11(m), 1.d0/chi2R11_W(bestfit(n),bestfit(9),&
!!$                                              bestfit(10),bestfit(11),bestfit(13),prior_sigma_int,m), Field(m)
!!$
!!$                                      End If
!!$
!!$                                   End If
!!$    
!!$                                End If
!!$
!!$                             End Do
!!$
!!$                          End Do
!!$
!!$                          close(UNIT_HP_FILE)
!!$
!!$                          call system('cd analyzer; python plot_HP.py')
!!$
!!$                       End If
!!$
!!$                    Else
!!$
!!$                       If (use_HP_per_host) then
!!$
!!$                          open(UNIT_HP_FILE,file='./output/chains/effective_hyperparameters_hosts.txt')
!!$
!!$                          Do n=1,number_of_hosts_galaxies
!!$ 
!!$                             write(UNIT_HP_FILE,*) n, 1.d0/chi2R11_W_host(bestfit(n),bestfit(9),bestfit(10),bestfit(11),&
!!$                                        bestfit(13),prior_sigma_int,n), host(n)
!!$
!!$                          End Do
!!$
!!$                          close(UNIT_HP_FILE)
!!$
!!$                          call system('cd analyzer; python plot_HP_hosts.py')
!!$ 
!!$                       Else If (use_HP_per_cepheid) then
!!$
!!$                          open(UNIT_HP_FILE,file='./output/chains/effective_hyperparameters_cepheids.txt')
!!$
!!$                          Do m=1,size(Field)
!!$
!!$                             Do n=1,number_of_hosts_galaxies
!!$ 
!!$                                If (host(n) .eq. Field(m)) then
!!$                                   
!!$                                   If (PeriodR11(m) .lt. cepheid_Period_limit) then
!!$
!!$                                      If ( chi2R11_W(bestfit(n),bestfit(9),bestfit(10),bestfit(11),bestfit(13),prior_sigma_int,m)&
!!$                                        .le. 1.d0 ) then
!!$
!!$                                         write(UNIT_HP_FILE,*) PeriodR11(m), observed_m_W(F160WR11(m),VIR11(m)) - &
!!$                                              P_L_relation_passband_W(bestfit(n),bestfit(9),bestfit(10),bestfit(11),bestfit(13),&
!!$                                              OHR11(m),PeriodR11(m)),eF160WR11(m), 1.d0, Field(m)
!!$
!!$                                      Else
!!$
!!$                                         write(UNIT_HP_FILE,*) PeriodR11(m), observed_m_W(F160WR11(m),VIR11(m)) - &
!!$                                              P_L_relation_passband_W(bestfit(n),bestfit(9),bestfit(10),bestfit(11),bestfit(13),&
!!$                                              OHR11(m),PeriodR11(m)), eF160WR11(m), 1.d0/chi2R11_W(bestfit(n),bestfit(9),&
!!$                                              bestfit(10),bestfit(11),bestfit(13),prior_sigma_int,m), Field(m)
!!$
!!$                                      End If
!!$
!!$                                   End If
!!$    
!!$                                End If
!!$
!!$                             End Do
!!$
!!$                          End Do
!!$
!!$                          close(UNIT_HP_FILE)
!!$
!!$                          If (use_HP_in_SNIa) then
!!$
!!$                             open(UNIT_HP_FILE,file='./output/chains/effective_hyperparameters_SNIa.txt')
!!$                          
!!$                             write(UNIT_HP_FILE,*) 'n4258', bestfit(9)-bestfit(9),&
!!$                                  bestfit(9)+5.d0*log10(bestfit(12))-25.d0,0.d0,&
!!$                                  bestfit(9)+5.d0*log10(bestfit(12))-25.d0,1.d0
!!$
!!$                             Do m=1,number_of_hosts_galaxies-1
!!$
!!$                                If ( chi2R11_SNIa(bestfit(m),bestfit(12),bestfit(14),m) .le. 1.d0) then
!!$
!!$                                   write(UNIT_HP_FILE,*) Fieldmvi(m), bestfit(m)-bestfit(9),mvi5av(m), Sigma_mvi5av(m), &
!!$                                        bestfit(m)+5.d0*log10(bestfit(12))-25.d0,1.d0
!!$
!!$                                Else
!!$
!!$                                   write(UNIT_HP_FILE,*) Fieldmvi(m), bestfit(m)-bestfit(9),mvi5av(m), Sigma_mvi5av(m), &
!!$                                        bestfit(m)+5.d0*log10(bestfit(12))-25.d0,&
!!$                                        1.d0/chi2R11_SNIa(bestfit(m),bestfit(12),bestfit(14),m)
!!$
!!$                                End If
!!$
!!$                             End Do
!!$
!!$                             close(UNIT_HP_FILE)
!!$
!!$                          Else
!!$
!!$                             continue
!!$
!!$                          End If
!!$
!!$                          call system('cd analyzer; python plot_HP_SNIa.py')
!!$
!!$                          If (use_HP_in_anchor) then
!!$
!!$                             If ( chi2R11_anchor_NGC4258(bestfit(9)) .le. 1.d0) then
!!$
!!$                                write(UNIT_EXE_FILE,*) 'HP FOR ANCHOR NGC4258 IS: ',1.d0
!!$
!!$                             Else
!!$
!!$                                write(UNIT_EXE_FILE,*) 'HP FOR ANCHOR NGC4258 IS: ',1.d0/chi2R11_anchor_NGC4258(bestfit(9))
!!$
!!$                             End If
!!$
!!$                          Else
!!$                                
!!$                             continue
!!$   
!!$                          End If
!!$                          
!!$                          call system('cd analyzer; python plot_HP_NGC4258_W.py')
!!$
!!$                       End If
!!$
!!$                    End If
!!$
!!$                 Else
!!$
!!$                    If (use_H_band) then
!!$
!!$                       print *, 'EFFECTIVE HYPER-PARAMETERS FOR  H BAND NOT IMPLEMENTED'
!!$
!!$                       stop
!!$
!!$                    Else
!!$
!!$                       print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                       stop
!!$                        
!!$                    End If
!!$
!!$                 End If
!!$
!!$              Else
!!$
!!$                 print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'
!!$
!!$                 stop
!!$
!!$              End If
!!$
!!$           End If
!!$
!!$        Else
!!$
!!$           If (separate_dataA .and. include_dataA) then
!!$
!!$              Do m=1,size(NameA)
!!$
!!$                 If (using_jeffreys_prior) then
!!$
!!$                    write(UNIT_EXE_FILE,*) 'Point ', m,' in data set A = ', 1.d0/chi2A_i(bestfit(1),bestfit(2),prior_sigma_int,m)
!!$                    
!!$                 Else
!!$
!!$                    If (chi2A_i(bestfit(1),bestfit(2),prior_sigma_int,m) .le. 1.d0 ) then
!!$
!!$                       write(UNIT_EXE_FILE,*) 'Point ', m,' in data set A = ', 1.d0
!!$
!!$                    Else
!!$
!!$                       write(UNIT_EXE_FILE,*) 'Point ', m,' in data set A = ', 1.d0/chi2A_i(bestfit(1),bestfit(2),prior_sigma_int,m)
!!$
!!$                    End If
!!$
!!$                 End If
!!$
!!$              End Do
!!$
!!$           Else if (include_dataA .and. .not.separate_dataA) then
!!$
!!$              !    write(15,*) 'For data set A = ', dble(size(NameA))/chi2A(bestfit(1),bestfit(2),bestfit(3))
!!$              write(UNIT_EXE_FILE,*) 'For data set A = ', dble(size(NameA))/chi2A(bestfit(1),bestfit(2),prior_sigma_int)
!!$
!!$           End If
!!$
!!$           If (separate_dataB .and. include_dataB) then
!!$
!!$              Do m=1,size(NameB)
!!$
!!$                 If (using_jeffreys_prior) then
!!$
!!$                    write(UNIT_EXE_FILE,*) 'Point ', m,' in data set B = ', 1.d0/chi2B_i(bestfit(1),bestfit(2),prior_sigma_int,m)
!!$
!!$                 Else
!!$
!!$                    If (chi2B_i(bestfit(1),bestfit(2),prior_sigma_int,m) .le. 1.d0 ) then
!!$
!!$                       write(UNIT_EXE_FILE,*) 'Point ', m,' in data set B = ', 1.d0
!!$
!!$                    Else
!!$
!!$                       write(UNIT_EXE_FILE,*) 'Point ', m,' in data set B = ', 1.d0/chi2B_i(bestfit(1),bestfit(2),prior_sigma_int,m)
!!$
!!$                    End If
!!$
!!$                 End If
!!$
!!$              End Do
!!$
!!$           Else if (include_dataB .and. .not.separate_dataB) then
!!$
!!$              !    write(15,*) 'For data set B = ', dble(size(NameB))/chi2B(bestfit(1),bestfit(2),bestfit(3))
!!$              write(UNIT_EXE_FILE,*) 'For data set B = ', dble(size(NameB))/chi2B(bestfit(1),bestfit(2),prior_sigma_int)
!!$
!!$           End If
!!$
!!$           If (separate_dataC .and. include_dataC) then
!!$
!!$              Do m=1,size(NameC)
!!$
!!$                 If (using_jeffreys_prior) then
!!$
!!$                    write(UNIT_EXE_FILE,*) 'Point ', m,' in data set C = ', 1.d0/chi2C_i(bestfit(1),bestfit(2),prior_sigma_int,m)
!!$
!!$                 Else
!!$
!!$                    If (chi2C_i(bestfit(1),bestfit(2),prior_sigma_int,m) .le. 1.d0) then
!!$
!!$                       write(UNIT_EXE_FILE,*) 'Point ', m,' in data set C = ', 1.d0
!!$
!!$                    Else
!!$
!!$                       write(UNIT_EXE_FILE,*) 'Point ', m,' in data set C = ', 1.d0/chi2C_i(bestfit(1),bestfit(2),prior_sigma_int,m)
!!$
!!$                    End If
!!$
!!$                 End If
!!$
!!$              End Do
!!$
!!$           Else if (include_dataC .and. .not.separate_dataC) then
!!$
!!$              !write(15,*) 'For data set C = ', dble(size(NameC))/chi2C(bestfit(1),bestfit(2),bestfit(3))
!!$              write(UNIT_EXE_FILE,*) 'For data set C = ', dble(size(NameC))/chi2C(bestfit(1),bestfit(2),prior_sigma_int)
!!$
!!$           End If
!!$
!!$           If (using_hyperparameters .and. use_HP_per_cepheid) then
!!$
!!$              open(UNIT_HP_FILE,file='./output/chains/effective_hyperparameters_cepheids.txt')
!!$
!!$              Do m=1,size(Name)
!!$                 
!!$                 If ( (Period(m) .gt. cepheid_lower_Period_limit) .and. (Period(m) .lt. cepheid_Period_limit)) then
!!$
!!$                    If (using_jeffreys_prior) then
!!$
!!$                       write(UNIT_EXE_FILE,*) 'Point ', Name(m),' in data set = ', &
!!$                            1.d0/chi2_i(bestfit(1),bestfit(2),prior_sigma_int,m)
!!$                    
!!$                    Else
!!$
!!$                       If (chi2_i(bestfit(1),bestfit(2),prior_sigma_int_LMC,m) .le. 1.d0 ) then
!!$
!!$                          write(UNIT_EXE_FILE,*) 'Point ', Name(m),' in data set = ', 1.d0
!!$
!!$                          write(UNIT_HP_FILE,*) Period(m), observed_wesenheit_magnitude(H(m),V(m),II(m)),&
!!$                               observed_wesenheit_magnitude(H(m),V(m),II(m)) - &
!!$                               wesenheit_magnitude(bestfit(1),bestfit(2),Period(m)),Sigma_m(m), 1.d0, 'LMC'
!!$
!!$
!!$                       Else
!!$
!!$                          write(UNIT_EXE_FILE,*) 'Point ', Name(m),' in data set = ', &
!!$                               1.d0/chi2_i(bestfit(1),bestfit(2),prior_sigma_int_LMC,m)
!!$
!!$                          write(UNIT_HP_FILE,*) Period(m), observed_wesenheit_magnitude(H(m),V(m),II(m)),&
!!$                               observed_wesenheit_magnitude(H(m),V(m),II(m)) - &
!!$                               wesenheit_magnitude(bestfit(1),bestfit(2),Period(m)), Sigma_m(m), &
!!$                               1.d0/chi2_i(bestfit(1),bestfit(2),prior_sigma_int_LMC,m), 'LMC'
!!$
!!$                       End If
!!$
!!$                    End If
!!$
!!$                 End If
!!$
!!$              End Do
!!$
!!$              close(UNIT_HP_FILE)
!!$
!!$              write(UNIT_EXE_FILE,*) '\ln P(\vec{w},D) at the bestfit is ', &
!!$                   ! log_Efstathiou_likelihood_hyperparameters(bestfit(1),bestfit(2),prior_sigma_int)
!!$                   log_Efstathiou_likelihood(bestfit(1),bestfit(2),prior_sigma_int_LMC)
!!$           
!!$           End If
!!$
!!$        End If
!!$
!!$     End If
!!$
!!$     close(UNIT_EXE_FILE)
!!$
!!$    If ((.not. testing_Gaussian_likelihood) .and. (.not. using_hyperparameters) ) then
!!$
!!$        deallocate (old_point,current_point,acceptance_probability,Name,Period,H,Sigma_m,V,II)
!!$
!!$    End If

End Program vsk




