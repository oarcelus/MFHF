module dens

  implicit none

  ! Parameters for spinorb = .true. case
  complex*16, save, public, allocatable     :: init_denssoc(:,:,:,:,:), dens_matsoc(:,:,:,:,:)
  real*8    , save, public, allocatable     :: init_magsoc(:,:)
  ! Parameters for spinorb = .false. case
  complex*16, save, public, allocatable     :: init_dens(:,:,:,:), dens_mat(:,:,:,:)
  real*8    , save, public, allocatable     :: init_mag(:)
  real*8    , save, public                  :: efermi, num

  public :: initialize_density_matrix
  public :: set_fermi_level
  public :: set_density_matrix_soc
  public :: set_density_matrix

  contains 

  subroutine initialize_density_matrix()

    use param,                 only : nelec, ntot, nsite, norb, sigma_x, sigma_y, sigma_z, spinorb
    use constants,             only : nspin, elecname, cmplx_0, cmplx_i

    implicit none

    logical                    :: io_find
    integer                    :: i, isite, ispin1, iorb, ispin, neps
    integer     , allocatable  :: nup(:), ndown(:)

    if (.not. spinorb) then

       allocate(init_mag(nsite))
       inquire(file=elecname, exist=io_find)
       if (io_find) then
          open(2, file=elecname, status='old')
             do isite = 1, nsite
                read(2,*) init_mag(isite)
             end do
          close(2)
       end if

    else if (spinorb) then

       allocate(init_magsoc(3,nsite))
       inquire(file=elecname, exist=io_find)
       if (io_find) then
          open(2, file=elecname, status='old')
             do isite = 1, nsite
                read(2,*) init_magsoc(:,isite)
             end do
          close(2)
       end if

    end if

    neps = nelec/nsite

    if (.not.spinorb) then

       allocate(init_dens(norb,norb,nsite,nspin))
       allocate(dens_mat(norb,norb,nsite,nspin))
       ! Initialize collinear spins, this is only provisional
       init_dens = cmplx_0
       do isite = 1, nsite
          do iorb = 1, norb
             do ispin = 1, nspin
                init_dens(iorb,iorb,isite,ispin) = real(neps)/(2.0*real(norb))
             end do
          end do   
          do ispin = 1, nspin
             init_dens(:,:,isite,ispin) = init_dens(:,:,isite,ispin) + &
                                          & + init_mag(isite)*sigma_z(:,:,ispin,ispin)/(2.0*real(norb))  
          end do
       end do


    else if (spinorb) then

       allocate(init_denssoc(norb,norb,nsite,nspin,nspin))
       allocate(dens_matsoc(norb,norb,nsite,nspin,nspin))
       ! Initialize collinear spins, this is only provisional
       init_denssoc = cmplx_0
       do isite = 1, nsite
          do iorb = 1, norb
             do ispin = 1, nspin
                init_denssoc(iorb,iorb,isite,ispin,ispin) = real(neps)/(2.0*real(norb))
             end do
          end do
          do iorb = 1, norb
             init_denssoc(iorb,iorb,isite,1,1) = init_denssoc(iorb,iorb,isite,1,1) + init_magsoc(3,isite)/(2.0*real(norb))
             init_denssoc(iorb,iorb,isite,2,2) = init_denssoc(iorb,iorb,isite,2,2) - init_magsoc(3,isite)/(2.0*real(norb))
             init_denssoc(iorb,iorb,isite,1,2) = init_denssoc(iorb,iorb,isite,1,2) + (init_magsoc(1,isite) + &
                                                                                      cmplx_i*init_magsoc(2,isite))/(2.0*real(norb))
             init_denssoc(iorb,iorb,isite,2,1) = init_denssoc(iorb,iorb,isite,2,1) + (init_magsoc(1,isite) - &
                                                                                      cmplx_i*init_magsoc(2,isite))/(2.0*real(norb))
          end do
       end do
    end if

  end subroutine initialize_density_matrix

  subroutine set_fermi_level(eigenval,beta)
    
    use param,          only : nkpt, ntot, nener, nkpt, nelec
    use utilities,      only : heavyside, fermi_dirac

    real*8  , intent(in)      :: eigenval(ntot,nkpt), beta
    real*8                    :: emax, emin, emax_tmp, emin_tmp, e, etmp, numax, numin
    integer                   :: ikpt, itot, i

    emax = maxval(eigenval(ntot,:))
    emin = minval(eigenval(1,:))
    numax = 0.0
    numin = 0.0
    do itot = 1, ntot
       do ikpt = 1, nkpt
          numax = numax + fermi_dirac((eigenval(itot,ikpt)-emax),beta)/real(nkpt)
          numin = numin + fermi_dirac((eigenval(itot,ikpt)-emin),beta)/real(nkpt)
       end do
    end do

    if ((numax > nelec) .and. (numin < nelec)) then
       do i = 1, 30
          num = 0
          e = (emax-emin)/2.0 + emin
          do itot = 1, ntot
             do ikpt = 1, nkpt
                num = num + fermi_dirac((eigenval(itot,ikpt)-e),beta)/real(nkpt)
             end do
          end do
          if (abs(num-nelec) < 1e-6) then
             efermi = e 
             exit
          else if (num < nelec) then
             emin = e
          else if (num > nelec) then
             emax = e
          end if
       end do
    end if

  end subroutine set_fermi_level

  subroutine set_density_matrix(eigenvec,eigenval,beta)

    use param,          only : nkpt, ntot, norb, nsite
    use constants,      only : nspin, cmplx_0
    use utilities,      only : fermi_dirac

    integer                :: ikpt, itot, iorb, iorb1, ispin, ispin1, isite
    complex*16, intent(in) :: eigenvec(norb*nsite,norb*nsite,nspin,nkpt)
    real*8,     intent(in) :: eigenval(norb*nsite,nspin,nkpt), beta

    dens_mat = cmplx_0
    do itot = 1, norb*nsite
       do ikpt = 1, nkpt
          do isite = 1, nsite
             do iorb = 1, norb
                do iorb1 = 1, norb
                   dens_mat(iorb,iorb1,isite,1) = dens_mat(iorb,iorb1,isite,1) + &
                                                  & + fermi_dirac((eigenval(itot,1,ikpt)-efermi),beta)*   &
                                                      conjg(eigenvec((isite-1)*norb+iorb,itot,1,ikpt))*   &
                                                            eigenvec((isite-1)*norb+iorb1,itot,1,ikpt)/real(nkpt)
                   dens_mat(iorb,iorb1,isite,2) = dens_mat(iorb,iorb1,isite,2) + &
                                                  & + fermi_dirac((eigenval(itot,2,ikpt)-efermi),beta)*   &
                                                      conjg(eigenvec((isite-1)*norb+iorb,itot,2,ikpt))*   &
                                                            eigenvec((isite-1)*norb+iorb1,itot,2,ikpt)/real(nkpt)
                end do
             end do
          end do
       end do
    end do

  end subroutine set_density_matrix

  subroutine set_density_matrix_soc(eigenvec,eigenval,beta)
  
    use param,          only : nkpt, ntot, norb, nsite
    use constants,      only : nspin, cmplx_0
    use utilities,      only : fermi_dirac

    integer                :: ikpt, itot, iorb, iorb1, ispin, ispin1, isite
    complex*16, intent(in) :: eigenvec(ntot,ntot,nkpt)
    real*8,     intent(in) :: eigenval(ntot,nkpt), beta

    dens_matsoc = cmplx_0
    do itot = 1, ntot
       do ikpt = 1, nkpt
          do isite = 1, nsite
             do ispin = 1, nspin
                do ispin1 = 1, nspin
                   do iorb = 1, norb
                      do iorb1 = 1, norb
                         dens_matsoc(iorb,iorb1,isite,ispin,ispin1) = dens_matsoc(iorb,iorb1,isite,ispin,ispin1) + &
                                                                    & + fermi_dirac((eigenval(itot,ikpt)-efermi),beta)                 *   &
                                                                    conjg(eigenvec((isite-1)*norb+iorb+(ispin-1)*norb*nsite,itot,ikpt))*   & 
                                                                    eigenvec((isite-1)*norb+iorb1+(ispin1-1)*norb*nsite,itot,ikpt)/nkpt
                      end do
                   end do
                end do
             end do    
          end do
       end do
    end do

  end subroutine set_density_matrix_soc  

end module dens
 
