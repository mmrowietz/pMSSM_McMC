!******************************************************
program HBwithLHClikelihood_SLHA
!
! In this example we evaluate the exclusion likelihoods ATLAS and CMS Higgs searches
! with tautau final states for the mhmod+ scenario. The input is provided from two
! datafiles (8 TeV and 13 Tev) created with SusHi, obtained here:
! https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGMSSMNeutral
!
! We also demonstrate how certain analyses can be deactivated in the standard
! HiggsBounds run. In particular, we deactivate here the latest ATLAS and CMS tau tau
! 95% CL limits from the standard HiggsBounds procedure, as we want to use the
! likelihood instead.
!
! After a successful run the output can be plotted with the python script
! plot_mhmodp_llh.py (needs matplotlib package!)
!
! (TS 24/11/2017)
!******************************************************

    use theory_XS_SM_functions
    use theory_BRfunctions
    use channels, only: HiggsBounds_deactivate_analyses, HiggsBounds_activate_all_analyses
    use output, only: createKey
    implicit none

    integer :: nH, nHplus, i, number_args
    character(len=100)::filename_in_8TeV, filename_in_13TeV, filename_out
    character(len=8) :: istring
    character(len=300) :: inputfilename, outputfilename, stem, temp, tmpstring
    integer :: HBresult, chan, ncombined, HBresult_all, chan_all, ncombined_all
    integer, parameter :: fileid = 78, fileid2 = 80
    double precision :: obsratio, obsratio_all
    integer :: error, status, n, cbin, npoints, ios
    double precision :: M_av, llh_CMS8, llh_CMS13, llh_ATLAS13, llh_ATLAS20, &
        llh_exp_CMS8, llh_exp_CMS13, llh_exp_ATLAS13, llh_exp_ATLAS20
    integer ::  Hindex, nc
! used by set_mass_uncertainties
    double precision, allocatable :: dmhneut(:)
    double precision, allocatable :: dmhch(:)

    nH = 3
    nHplus = 1
!    theory_uncertainty_1s = 1.5D0

    allocate (dmhneut(nH), dmhch(nHplus))

    number_args = IARGC()

    if (number_args .ne. 2) then
        stop "Incorrect number of arguments given to HBwithSLHA"
    endif

    ! Read arguments into text strings.
    i = 1
    temp = ""
    call GETARG(i, temp)
    read (temp, *) npoints

    i = i + 1
    temp = ""
    call GETARG(i, temp)
    stem = ""
    stem = trim(temp)

    call initialize_HiggsBounds(nH, nHplus, 'onlyH')

    filename_in_8TeV = "../example_data/Mh125/mh125_8.tsv"
    filename_in_13TeV = "../example_data/Mh125/mh125_13.tsv"
    filename_out = "Mh125_HBwithLHClikelihood.dat"
    call system('rm -f Mh125_HBwithLHClikelihood.dat')

    open (432, file=trim(adjustl(filename_in_8TeV)), action='read', status='old', iostat=status)
    if (status .ne. 0) then
        write (*, *) 'Bad status', status, 'with the following file:'
        write (*, *) trim(adjustl(filename_in_8TeV))
        stop
    endif
    read(432, *) ! skip header

    open (433, file=trim(adjustl(filename_in_13TeV)), action='read', status='old', iostat=status)
    if (status .ne. 0) then
        write (*, *) 'Bad status', status, 'with the following file:'
        write (*, *) trim(adjustl(filename_in_13TeV))
        stop
    endif
    read(433, *) ! skip header

    open (434, file=trim(adjustl(filename_out)), action='write', status='new')

    do i = 1, npoints

       write (*, *) "number of processed points: ", i

        if (i .gt. 99999999) stop 'need to increase the size of istring in HBwithSLHA'
        write (istring, '(I8)') i

        inputfilename = trim(adjustl(stem))//'.'//trim(adjustl(istring))

        !! Test if input file exists and is non-empty
        open (fileid2, file=inputfilename, form='formatted')
        read (fileid2, '(A)', iostat=ios) tmpstring

        if (ios .eq. 0) then
            close (fileid2)

            call HiggsBounds_input_SLHA(inputfilename)

! Activate all analyses (in case some of them have been deactivated before)
            call HiggsBounds_activate_all_analyses
        
! Run the standard HiggsBounds routine considering all analyses
            call run_HiggsBounds(HBresult_all, chan_all, obsratio_all, ncombined_all)

! Deactivate current CMS and ATLAS searches for non-standard Higgs to tautau
            call HiggsBounds_deactivate_analyses((/14029, 2014049, 20140492, 170907242, 200212223/))
! Standard HiggsBounds run (gives 95% CL limit):
            call run_HiggsBounds(HBresult, chan, obsratio, ncombined)

! Obtain exclusion-likelihood from H->tautau searches:
!
!  The arguments of the following subroutines mean the following:
!        Analysis-ID --> 14029 for CMS 8 TeV results,
!                   17020 for CMS 13 TeV results,
!                   170907242 for ATLAS 13 TeV results.
!                   200212223 for new ATLAS 13 TeV results
!   Hindex  -> Index of the Higgs boson that was selected as most sensitive [int output]
!   M_av    -> mass position where limit is extracted (a signal strength-weighted mass
!              average in case of combined Higgs bosons) [dbl output]
!   nc      -> number of combined Higgs bosons [int output]
!   cbin    -> binary code of the combined Higgs bosons (see manual) [int output]
!   llh     -> -2ln L value [dbl output]
!   obspred -> 'obs' or 'pred' to chose whether the observed or expected likelihood should be
!              extracted. [char input]
!
! Get expected/predicted likelihood
            call HiggsBounds_get_likelihood(14029, Hindex, nc, cbin, M_av, llh_exp_CMS8, 'pred')
! ! Get observed likelihood
            call HiggsBounds_get_likelihood(14029, Hindex, nc, cbin, M_av, llh_CMS8, 'obs')

! Get expected/predicted likelihood
            call HiggsBounds_get_likelihood(17020, Hindex, nc, cbin, M_av, llh_exp_CMS13, 'pred')
! Get observed likelihood
            call HiggsBounds_get_likelihood(17020, Hindex, nc, cbin, M_av, llh_CMS13, 'obs')

! Get expected/predicted likelihood
            call HiggsBounds_get_likelihood(170907242, Hindex, nc, cbin, M_av, llh_exp_ATLAS13, 'pred')
! Get observed likelihood
            call HiggsBounds_get_likelihood(170907242, Hindex, nc, cbin, M_av, llh_ATLAS13, 'obs')

! Get expected/predicted likelihood
            call HiggsBounds_get_likelihood(200212223, Hindex, nc, cbin, M_av, llh_exp_ATLAS20, 'pred')
! Get observed likelihood
            call HiggsBounds_get_likelihood(200212223, Hindex, nc, cbin, M_av, llh_ATLAS20, 'obs')

            write (434, *) HBresult, chan, obsratio, ncombined, &
                 HBresult_all, chan_all, obsratio_all, ncombined_all, &
                 Hindex, M_av, nc, cbin, llh_CMS8, llh_exp_CMS8, llh_CMS13, &
                 llh_exp_CMS13, llh_ATLAS13, llh_exp_ATLAS13, llh_ATLAS20, llh_exp_ATLAS20
            
         else
            close (fileid2)
            call system("rm -f "//inputfilename)
         endif

      enddo

! Write out the key with the used analyses (indicating possibly deactivated analyses)
    call createKey("HB_with_deactivated_analyses_")
    close (432)
    close (433)
    close (434)

    call finish_HiggsBounds
end program HBwithLHClikelihood_SLHA
