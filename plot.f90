module plots

  ! case for spinorb = .true.
  complex*16, save, public, allocatable   :: magsoc(:,:)
  ! case for spinorb = .false.
  complex*16, save, public, allocatable   :: mag(:)

  public :: get_magnetization
  public :: plot_dos
  public :: plot_band

  contains

  subroutine get_magnetization(density)

    use param,         only : ntot, nkpt, norb, nsite, sigma_x, sigma_y, sigma_z
    use constants,     only : nspin, cmplx_0, cmplx_i
    use utilities,     only : trace

    complex*16, intent(in) :: density(norb,norb,nsite,nspin)
    integer                :: isite, iorb, iorb1, ispin, ispin1

    allocate(mag(nsite))

    mag = cmplx_0
    do isite = 1, nsite
       mag(isite) = mag(isite) + trace(density(:,:,isite,1)-density(:,:,isite,2),norb)
    end do

  end subroutine get_magnetization

  subroutine get_magnetization_soc(density)

    use param,         only : ntot, nkpt, norb, nsite, sigma_x, sigma_y, sigma_z
    use constants,     only : nspin, cmplx_0, cmplx_i
    use utilities,     only : trace
   

    complex*16, intent(in) :: density(norb,norb,nsite,nspin,nspin)
    integer                :: isite, iorb, iorb1, ispin, ispin1
 
    allocate(magsoc(3,nsite))

    magsoc = cmplx_0
    do isite = 1, nsite
       magsoc(1,isite) = magsoc(1,isite) +         trace(density(:,:,isite,1,2)+density(:,:,isite,2,1),norb)
       magsoc(2,isite) = magsoc(2,isite) + cmplx_i*trace(density(:,:,isite,2,1)-density(:,:,isite,1,2),norb)
       magsoc(3,isite) = magsoc(3,isite) +         trace(density(:,:,isite,1,1)-density(:,:,isite,2,2),norb)
    end do
  
  end subroutine get_magnetization_soc

  subroutine plot_dos(val,vec,fer)

    use param,         only : ntot, nsite, norb, nkpt, nener, smear
    use utilities,     only : gauss
    use constants,     only : nspin

    real*8                      :: emax, emin, e, num, etmp
    real*8   , intent(in)       :: val(nsite*norb,nspin,nkpt),  fer
    complex*16, intent(in)      :: vec(nsite*norb,nsite*norb,nspin,nkpt)
    integer                     :: itot, ikpt, iener, ispin


    emax = maxval(val(:,:,:))
    emin = minval(val(:,:,:))

    emax = emax+1.0
    emin = emin-1.0

    etmp = emin
    open(unit=11, file='dos.out', status='replace')
       do ispin = 1, nspin
          do iener = 1, nener
             e = (iener*(emax-emin)/nener + emin)
             num = 0.0
             do ikpt = 1, nkpt
                do itot = 1, norb*nsite
                   num = num + dot_product(conjg(vec(:,itot,ispin,ikpt)), &
                                           vec(:,itot,ispin,ikpt))*gauss(e,val(itot,ispin,ikpt),smear)/nkpt
                end do
             end do
             if (ispin == 1) then
                write(11,*) e-fer, num
             else if (ispin == 2) then
                write(11,*) e-fer, -1.0*num
             end if
             etmp = e
          end do
       end do
    close(11)
  
  end subroutine plot_dos

  subroutine plot_band(fer)

    use hamilton,            only : hf_pot, hf_potsoc
    use constants,           only : pathname, cmplx_i, cmplx_0, nspin    
    use param,               only : recipvec, ham_r, ham_rsoc, r_bond, nrpt, norb, nsite, ntot, atom_label, spinorb
    use utilities,           only : utility_frac_to_cart, utility_diag

    implicit none

    logical                  :: io_find
    integer                  :: npth, npnts, totpth, irpt, ipth, itot, itot1,iorb,iorb1, iline, num, io_stat, iat, iat1, ispin, ispin1
    real*8 , intent(in)      :: fer
    real*8                   :: rdotk, rkpt, r(3)
    real*8,  allocatable     :: kpt_pathl(:,:), kpt_pathr(:,:), pathpoints(:,:), band_valsoc(:,:), band_val(:,:,:)
    complex*16, allocatable  :: ham_pth(:,:,:,:,:,:), ham_bands(:,:,:,:), ham_pthsoc(:,:,:,:,:,:,:), ham_bandssoc(:,:,:)
    complex*16               :: fac

    ! Read paths
    inquire(file=pathname, exist=io_find)
    if (io_find) then
       open(40, file=pathname, status='old')
          read(40,*) npth, npnts
          allocate(pathpoints(3,npth))
          totpth = 0
          do ipth = 1, npth
             read(40,*) pathpoints(:,ipth)
          end do
       close(40)
    else
       print*, 'No path file found'
       stop
    end if

    ! Linearly interpolate the end path points to obtain kpoints along high sym
    ! lines
    totpth = (npth-1)*npnts + 1 
    allocate(kpt_pathl(3,totpth))
    allocate(kpt_pathr(3,totpth))

    num = 0
    do ipth = 1, npth - 1
       do iline = 1, npnts
          num = num + 1
          kpt_pathl(:,num) = (iline-1)*(pathpoints(:,ipth+1) - pathpoints(:,ipth))/npnts + pathpoints(:,ipth)
       end do
    end do

    kpt_pathl(:,totpth) = pathpoints(:,npth)

    do ipth = 1, totpth
       call utility_frac_to_cart(kpt_pathl(:,ipth),kpt_pathr(:,ipth),recipvec)
    end do
 
    if (.not.spinorb) then

       allocate(ham_pth(norb,norb,nsite,nsite,totpth,nspin))
       ham_pth = cmplx_0
       do ipth = 1, totpth   
          do irpt = 1, nrpt
             rdotk = dot_product(kpt_pathr(:,ipth),r_bond(:,irpt))
             fac = exp(cmplx_i*rdotk)
             ham_pth(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,1) = ham_pth(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,1) + &
                                                                         fac*ham_r(:,:,irpt,1)
             ham_pth(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,2) = ham_pth(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,2) + &
                                                                         fac*ham_r(:,:,irpt,2)
          end do
       end do

       allocate(ham_bands(norb*nsite,norb*nsite,nspin,totpth))
       !Sum kinetic part to the hamiltonian
       ham_bands = cmplx_0
       do ipth = 1, totpth
          do iat = 1, nsite
             do iat1 = 1, nsite
                ham_bands((iat-1)*norb+1:iat*norb,(iat1-1)*norb+1:iat1*norb,1,ipth) = ham_pth(:,:,iat,iat1,ipth,1)
                ham_bands((iat-1)*norb+1:iat*norb,(iat1-1)*norb+1:iat1*norb,2,ipth) = ham_pth(:,:,iat,iat1,ipth,2)
             end do
          end do
       end do

       !Sum hf_potential
       do ipth = 1, totpth
          do iat = 1, nsite
             ham_bands((iat-1)*norb+1:iat*norb,(iat-1)*norb+1:iat*norb,1,ipth) =   &
             ham_bands((iat-1)*norb+1:iat*norb,(iat-1)*norb+1:iat*norb,1,ipth) +   &
             hf_pot(:,:,iat,1)
             ham_bands((iat-1)*norb+1:iat*norb,(iat-1)*norb+1:iat*norb,2,ipth) =   &
             ham_bands((iat-1)*norb+1:iat*norb,(iat-1)*norb+1:iat*norb,2,ipth) +   &
             hf_pot(:,:,iat,2)
          end do
       end do
 
       allocate(band_val(norb*nsite,nspin,totpth))

       do ipth = 1, totpth
          call utility_diag(ham_bands(:,:,1,ipth),band_val(:,1,ipth),norb*nsite)
          call utility_diag(ham_bands(:,:,2,ipth),band_val(:,2,ipth),norb*nsite)
       end do

       open(unit = 23, file = 'eigpath.out', status = 'replace')
       write(23,*) totpth
       do ipth = 1, totpth
          do ispin = 1, nspin
             write(23,*) kpt_pathr(:,ipth), ispin
             do itot = 1, norb*nsite
                write(23,*) itot, band_val(itot,ispin,ipth)
             end do
             write(23,*) ' '
          end do
       end do
       close(23)

       open(unit = 23, file = 'wavpath.out', status = 'replace')
       do ipth = 1, totpth
          do ispin = 1, nspin
             write(23,*) kpt_pathr(:,ipth), ispin
             do itot = 1, norb*nsite
                do itot1 = 1, norb*nsite
                   write(23,*) itot1, itot, real(ham_bands(itot1,itot,ispin,ipth)), aimag(ham_bands(itot1,itot,ispin,ipth))
                end do
             end do
             write(23,*) ' '
          end do
       end do
       close(23)

       open(unit = 23, file = 'bands.1.out', status = 'replace')
       do itot = 1, norb*nsite
          rkpt = 0.0
          write(23,*) rkpt, band_val(itot,1,1) - fer
          do ipth = 1, totpth - 1
             r(:) = kpt_pathr(:,ipth+1) - kpt_pathr(:,ipth)
             rkpt = rkpt + sqrt(r(1)**2 + r(2)**2 + r(3)**2)
             write(23,*) rkpt, band_val(itot,1,ipth+1) - fer
          end do
          write(23,*) ' '
       end do
       close(23)   

       open(unit = 23, file = 'bands.2.out', status = 'replace')
       do itot = 1, norb*nsite
          rkpt = 0.0
          write(23,*) rkpt, band_val(itot,2,1) - fer
          do ipth = 1, totpth - 1
             r(:) = kpt_pathr(:,ipth+1) - kpt_pathr(:,ipth)
             rkpt = rkpt + sqrt(r(1)**2 + r(2)**2 + r(3)**2)
             write(23,*) rkpt, band_val(itot,2,ipth+1) - fer
          end do
          write(23,*) ' '
       end do
       close(23)


    else if (spinorb) then

       allocate(ham_pthsoc(norb,norb,nsite,nsite,totpth,nspin,nspin))
       ham_pthsoc = cmplx_0
       do ipth = 1, totpth
          do irpt = 1, nrpt
             rdotk = dot_product(kpt_pathr(:,ipth),r_bond(:,irpt))
             fac = exp(cmplx_i*rdotk)
             ham_pthsoc(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,1,1) = ham_pthsoc(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,1,1) + &
                                                                              fac*ham_rsoc(:,:,irpt,1,1)
             ham_pthsoc(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,2,2) = ham_pthsoc(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,2,2) + &
                                                                              fac*ham_rsoc(:,:,irpt,2,2)
             ham_pthsoc(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,1,2) = ham_pthsoc(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,1,2) + &
                                                                              fac*ham_rsoc(:,:,irpt,1,2)
             ham_pthsoc(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,2,1) = ham_pthsoc(:,:,atom_label(1,irpt),atom_label(2,irpt),ipth,2,1) + &
                                                                              fac*ham_rsoc(:,:,irpt,2,1)
          end do
       end do

       allocate(ham_bandssoc(ntot,ntot,totpth))
       !Sum kinetic part to the hamiltonian
       ham_bandssoc = cmplx_0
       do ipth = 1, totpth
          do ispin = 1, nspin
             do ispin1 = 1, nspin
                do iat = 1, nsite
                   do iat1 = 1, nsite
                      ham_bandssoc((iat-1)*norb+1+(ispin-1)*norb*nsite:iat*norb+(ispin-1)*norb*nsite, &
                               (iat1-1)*norb+1+(ispin1-1)*norb*nsite:iat1*norb+(ispin1-1)*norb*nsite,ipth) = &
                               ham_pthsoc(:,:,iat,iat1,ipth,ispin,ispin1)
                   end do
                end do
             end do
          end do
       end do

       !Sum hf_potential
       do ipth = 1, totpth
          do ispin = 1, nspin
             do ispin1 = 1, nspin
                do iat = 1, nsite
                   ham_bandssoc((iat-1)*norb+1+(ispin-1)*norb*nsite:iat*norb+(ispin-1)*norb*nsite, &
                            (iat-1)*norb+1+(ispin1-1)*norb*nsite:iat*norb+(ispin1-1)*norb*nsite,ipth) =   &
                   ham_bandssoc((iat-1)*norb+1+(ispin-1)*norb*nsite:iat*norb+(ispin-1)*norb*nsite, &
                            (iat-1)*norb+1+(ispin1-1)*norb*nsite:iat*norb+(ispin1-1)*norb*nsite,ipth) +   &
                   hf_potsoc(:,:,iat,ispin,ispin1)
                end do
             end do
          end do
       end do

       allocate(band_valsoc(ntot,totpth))

       do ipth = 1, totpth
          call utility_diag(ham_bandssoc(:,:,ipth),band_valsoc(:,ipth),ntot)
       end do

       open(unit = 23, file = 'eigpath.out', status = 'replace')
       write(23,*) totpth
       do ipth = 1, totpth
          write(23,*) kpt_pathr(:,ipth)
          do itot = 1, ntot
             write(23,*) itot, band_valsoc(itot,ipth)
          end do
          write(23,*) ' '
       end do
       close(23)

       open(unit = 23, file = 'wavpath.out', status = 'replace')
       do ipth = 1, totpth
          write(23,*) kpt_pathr(:,ipth)
          do itot = 1, ntot
             do itot1 = 1, ntot
                write(23,*) itot1, itot, real(ham_bandssoc(itot1,itot,ipth)), aimag(ham_bandssoc(itot1,itot,ipth))
             end do
          end do
          write(23,*) ' '
       end do
       close(23)

       open(unit = 23, file = 'bands.out', status = 'replace')
       do itot = 1, ntot
          rkpt = 0.0
          write(23,*) rkpt, band_valsoc(itot,1) - fer
          do ipth = 1, totpth - 1
             r(:) = kpt_pathr(:,ipth+1) - kpt_pathr(:,ipth)
             rkpt = rkpt + sqrt(r(1)**2 + r(2)**2 + r(3)**2)
             write(23,*) rkpt, band_valsoc(itot,ipth+1) - fer
          end do
          write(23,*) ' '
       end do
       close(23)

    end if

  end subroutine plot_band

end module plots
