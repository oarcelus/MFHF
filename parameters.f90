module param 
  implicit none

  ! File names
  character(len=50), save, public                :: atomname, soctype, bondcoord

  ! Calculation parameters
  complex*16     , save, public, allocatable     :: u_mat(:,:,:,:,:), sigma_x(:,:,:,:), sigma_y(:,:,:,:), sigma_z(:,:,:,:)
  integer        , save, public, allocatable     :: atom_label(:,:)
  integer        , save, public                  :: nsite, nelec, norb, nener, nkpt, ndim, ntot, nrpt, nscf
  real*8         , save, public                  :: lattvec(3,3), recipvec(3,3)
  real*8         , save, public, allocatable     :: kpt_latt(:,:), kpt_cart(:,:), atom_dir(:,:), atom_cart(:,:), r_bond(:,:)
  real*8         , save, public                  :: volume, ediff, mix, etemp, smear
  logical        , save, public                  :: spinorb, dosplot, bandplot, lsph
  ! Parameters for spinorb = .true. case
  complex*16     , save, public, allocatable     :: ham_ksoc(:,:,:,:,:,:,:), ham_rsoc(:,:,:,:,:)
  ! Parameters for spinorb = .false. case
  complex*16     , save, public, allocatable     :: ham_k(:,:,:,:,:,:), ham_r(:,:,:,:)


  ! Subroutines
  public    :: get_names
  public    :: get_latt_parameters
  public    :: get_kpt_parameters
  public    :: get_pauli
  public    :: get_kinetic_parameters
  public    :: get_coulomb_matrix
  public    :: spherical_param

  contains

  ! Read filenames
  subroutine get_names()

    open(unit=1, file= 'hf.in', status='old')
       read(1,*) ! Space for separating file names and parameters
       read(1,*) spinorb ! This is true if a DFT + SO calculation is performed
       read(1,*) soctype ! 'f': If read directly from LDA+SO, 's': If read SO potential read separatedly and sum to LDA.
                         ! 'n': no SO potential added, just for testing purposes
       read(1,*) bondcoord ! 'c' if in hf_ham.in bonds are given in cartesian coords, 'f' if fractional coords.
       read(1,*) nelec ! Total number of electrons
       read(1,*) nscf  ! Total number of SCF cycles
       read(1,*) nener ! Total number of points in between emax and emin
       read(1,*) etemp  ! Electron temperature for smearin in K
       read(1,*) ediff ! Convergence criteria
       read(1,*) mix ! Mixing of density
       read(1,*) dosplot ! True if DOS is going to be plot
       read(1,*) smear ! gaussian broadening for DOS plot
       read(1,*) bandplot ! True if BANDS is going to be plot
       read(1,*) lsph ! True if spherical parameterization of coulomb potential
    close(1) 
  end subroutine get_names

  subroutine get_latt_parameters()

    use utilities,                only : utility_recip_lattice, utility_frac_to_cart
    use constants,                only : lattname

    implicit none

    logical                       :: io_find
    integer                       :: ilatt, isite

    inquire(file=lattname, exist=io_find)
    if (io_find) then
       open(2, file=lattname, status='old')
          read(2,*)
          do ilatt = 1, 3
             read(2,*) lattvec(:,ilatt)
          end do
          read(2,*) atomname, nsite
          allocate(atom_dir(3,nsite))
          allocate(atom_cart(3,nsite))
          do isite = 1, nsite
             read(2,*) atom_dir(:,isite)
          end do
       close(2)
    else
       print*, 'No lattvec file found'
       stop
    end if

    call utility_recip_lattice(lattvec, recipvec, volume)
    do isite = 1, nsite
       call utility_frac_to_cart(atom_dir(:,isite),atom_cart(:,isite),lattvec)
    end do

  end subroutine get_latt_parameters

  subroutine get_kpt_parameters()

    use utilities,                only : utility_frac_to_cart
    use constants,                only : kptname

    implicit none

    integer                       :: ikpt
    logical                       :: io_find

    inquire(file=kptname, exist=io_find)
    if (io_find) then
       open(2, file=kptname, status='old')
          read(2,*)
          read(2,*) nkpt
          allocate(kpt_latt(3,nkpt))
          allocate(kpt_cart(3,nkpt))
          do ikpt = 1, nkpt
             read(2,*) kpt_latt(:,ikpt)
             call utility_frac_to_cart(kpt_latt(:,ikpt),kpt_cart(:,ikpt),recipvec)
          end do
       close(2)
    else
       print*, 'No k-points file found'
       stop
    end if
  
  end subroutine get_kpt_parameters

  ! Get parameters

  subroutine get_kinetic_parameters()

    use constants,                 only : cmplx_i, cmplx_0, twopi, hamname, nspin, socname
    use utilities,                 only : utility_frac_to_cart

    implicit none

    logical                        :: io_find
    integer                        :: io_stat, irpt, ikpt, iat, iat1, iorb, iorb1, norbso, isite,isite1
    real*8                         :: rdotk, bond(3)
    complex*16                     :: fac
    real*8         , allocatable   :: rtempnn(:,:), itempnn(:,:)
    complex*16     , allocatable   :: ham(:,:,:), tempnn(:,:), hamsoc(:,:,:,:), vso(:,:,:,:,:)

    ! Case if DFT + SO is NOT used
    if (.not.spinorb) then
       io_stat = 0
       inquire(file=hamname, exist=io_find)
       if (io_find) then
          open(2, file=hamname, status='old')
             read(2,*,iostat=io_stat) iat, iat1, norb, norb
          close(2)
          ntot = norb*nspin*nsite
          allocate(ham(norb,norb,nspin))
          allocate(ham_k(norb,norb,nsite,nsite,nkpt,nspin))
          ham_k = cmplx_0
          do ikpt = 1, nkpt
             open(2, file=hamname, status='old')
             ham = cmplx_0
             nrpt = 0
             do while (io_stat == 0)
                read(2,*,iostat=io_stat) iat, iat1, norb, norb
                if (io_stat < 0) cycle
                read(2,*,iostat=io_stat) bond(:)
                do iorb = 1, norb
                   read(2,*,iostat=io_stat) ham(iorb,:,1)
                end do

                if (bondcoord == 'f') then
                   call utility_frac_to_cart(bond(:),bond(:),lattvec)
                end if

                nrpt = nrpt + 1
                !Expand to spindown subspace
                ham(:,:,2) = ham(:,:,1)
                rdotk = dot_product(kpt_cart(:,ikpt),bond(:))
                fac = exp(cmplx_i*rdotk)
                ham_k(:,:,iat,iat1,ikpt,1) = ham_k(:,:,iat,iat1,ikpt,1) + fac*ham(:,:,1)
                ham_k(:,:,iat,iat1,ikpt,2) = ham_k(:,:,iat,iat1,ikpt,2) + fac*ham(:,:,2)
             end do
             io_stat = 0
             close(2)
          end do

          allocate(ham_r(norb,norb,nrpt,nspin))
          allocate(r_bond(3,nrpt))
          allocate(atom_label(2,nrpt))        

          open(2, file=hamname, status='old')
          ham_r = cmplx_0
          do irpt = 1, nrpt
             read(2,*) atom_label(1,irpt), atom_label(2,irpt)
             read(2,*) r_bond(:,irpt)
           
             if (bondcoord == 'f') then
                call utility_frac_to_cart(r_bond(:,irpt),r_bond(:,irpt),lattvec)
             end if

             do iorb = 1, norb
                read(2,*) ham_r(iorb,:,irpt,1)
             end do
             ham_r(:,:,irpt,2) = ham_r(:,:,irpt,1)
          end do
          close(2)
       else
          print*, 'No kinetic hamiltonian file found'
          stop
       end if
    ! Case where DFT+SO is used
    else if (spinorb.and.(soctype=='f')) then
       io_stat = 0
       inquire(file=hamname, exist=io_find)
       if (io_find) then
          open(2, file=hamname, status='old')
             read(2,*,iostat=io_stat) iat, iat1, norbso, norbso
          close(2)
          ntot = norbso*nsite
          norb = norbso/nspin
          allocate(rtempnn(norbso,norbso))
          allocate(itempnn(norbso,norbso))
          allocate(tempnn(norbso,norbso))
          allocate(hamsoc(norb,norb,nspin,nspin))
          allocate(ham_ksoc(norb,norb,nsite,nsite,nkpt,nspin,nspin))
          ham_ksoc = cmplx_0
          do ikpt = 1, nkpt
             hamsoc = cmplx_0
             nrpt = 0
             open(2, file=hamname, status='old')
                do while (io_stat == 0)
                   read(2,*,iostat=io_stat) iat, iat1, norbso, norbso
                   if (io_stat < 0) cycle
                   read(2,*,iostat=io_stat) bond(:)

                   if (bondcoord == 'f') then
                      call utility_frac_to_cart(bond(:),bond(:),lattvec)
                   end if
 
                   do iorb = 1, norbso
                      read(2,*,iostat=io_stat) rtempnn(iorb,:)
                   end do
                   read(2,*,iostat=io_stat)
                   do iorb = 1, norbso
                      read(2,*,iostat=io_stat) itempnn(iorb,:)
                   end do
                   nrpt = nrpt + 1
                   ! Get from 2n*n to n*n complex matrix
                   do iorb = 1, norbso
                      do iorb1 = 1, norbso
                         tempnn(iorb,iorb1) = cmplx(rtempnn(iorb,iorb1),itempnn(iorb,iorb1))
                      end do
                   end do
                   ! Format the hamiltonian to the ham_r variable
                   hamsoc(:,:,1,1) = tempnn(1:norb,1:norb)
                   hamsoc(:,:,2,2) = tempnn(norb+1:nspin*norb,norb+1:nspin*norb)
                   hamsoc(:,:,1,2) = tempnn(1:norb,norb+1:nspin*norb)
                   hamsoc(:,:,2,1) = tempnn(norb+1:nspin*norb,1:norb)

                   rdotk = dot_product(kpt_cart(:,ikpt),bond(:))
                   fac = exp(cmplx_i*rdotk)
                   ham_ksoc(:,:,iat,iat1,ikpt,1,1) = ham_ksoc(:,:,iat,iat1,ikpt,1,1) + fac*hamsoc(:,:,1,1)
                   ham_ksoc(:,:,iat,iat1,ikpt,2,2) = ham_ksoc(:,:,iat,iat1,ikpt,2,2) + fac*hamsoc(:,:,2,2)
                   ham_ksoc(:,:,iat,iat1,ikpt,1,2) = ham_ksoc(:,:,iat,iat1,ikpt,1,2) + fac*hamsoc(:,:,1,2)
                   ham_ksoc(:,:,iat,iat1,ikpt,2,1) = ham_ksoc(:,:,iat,iat1,ikpt,2,1) + fac*hamsoc(:,:,2,1)
                end do
                io_stat = 0
             close(2)
          end do

          allocate(ham_rsoc(norb,norb,nrpt,nspin,nspin))
          allocate(r_bond(3,nrpt))
          allocate(atom_label(2,nrpt))

          open(2, file=hamname, status='old')
          ham_rsoc = cmplx_0
          do irpt = 1, nrpt
             read(2,*,iostat=io_stat) atom_label(1,irpt), atom_label(2,irpt)
             read(2,*,iostat=io_stat) r_bond(:,irpt)

             if (bondcoord == 'f') then
                call utility_frac_to_cart(r_bond(:,irpt),r_bond(:,irpt),lattvec)
             end if

             do iorb = 1, norbso
                read(2,*,iostat=io_stat) rtempnn(iorb,:)
             end do
             read(2,*,iostat=io_stat)
             do iorb = 1, norbso
                read(2,*,iostat=io_stat) itempnn(iorb,:)
             end do

             do iorb = 1, norbso
                do iorb1 = 1, norbso
                   tempnn(iorb,iorb1) = cmplx(rtempnn(iorb,iorb1),itempnn(iorb,iorb1))
                end do
             end do
             ham_rsoc(:,:,irpt,1,1) = tempnn(1:norb,1:norb)
             ham_rsoc(:,:,irpt,2,2) = tempnn(norb+1:nspin*norb,norb+1:nspin*norb)
             ham_rsoc(:,:,irpt,1,2) = tempnn(1:norb,norb+1:nspin*norb)
             ham_rsoc(:,:,irpt,2,1) = tempnn(norb+1:nspin*norb,1:norb)
          end do
          close(2)
          deallocate(rtempnn)
          deallocate(itempnn)
          deallocate(tempnn)
       else
          print*, 'No kinetic hamiltonian file found'
          stop
       end if

    else if (spinorb.and..not.(soctype=='f')) then

       io_stat = 0
       inquire(file=socname, exist=io_find)
       if (io_find) then
          open(2, file=socname, status='old')
             read(2,*,iostat=io_stat) iat, iat1, norbso, norbso
          close(2)
          norb = norbso/nspin
          allocate(rtempnn(norbso,norbso))
          allocate(itempnn(norbso,norbso))
          allocate(tempnn(norbso,norbso))
          allocate(vso(norb,norb,nsite,nspin,nspin))
          vso = cmplx_0
          tempnn = cmplx_0

          open(2, file=socname, status='old')
          do isite = 1, nsite
             read(2,*,iostat=io_stat)
             read(2,*,iostat=io_stat)
             do iorb = 1, norbso
                read(2,*,iostat=io_stat) rtempnn(iorb,:)
             end do
             read(2,*,iostat=io_stat)
             do iorb = 1, norbso
                read(2,*,iostat=io_stat) itempnn(iorb,:)
             end do

             do iorb = 1, norbso
                do iorb1 = 1, norbso
                   tempnn(iorb,iorb1) = cmplx(rtempnn(iorb,iorb1),itempnn(iorb,iorb1))
                end do
             end do
             vso(:,:,isite,1,1) = tempnn(1:norb,1:norb)
             vso(:,:,isite,2,2) = tempnn(norb+1:nspin*norb,norb+1:nspin*norb)
             vso(:,:,isite,1,2) = tempnn(1:norb,norb+1:nspin*norb)
             vso(:,:,isite,2,1) = tempnn(norb+1:nspin*norb,1:norb)
          end do
          close(2)
       else
          print*, 'No hf_soc.in file found'
          stop
       end if

       if (soctype=='n') vso = cmplx_0

       io_stat = 0
       inquire(file=hamname, exist=io_find)
       if (io_find) then
          open(2, file=hamname, status='old')
             read(2,*,iostat=io_stat) iat, iat1, norb, norb
          close(2)
          ntot = norb*nspin*nsite
          allocate(hamsoc(norb,norb,nspin,nspin))
          allocate(ham_ksoc(norb,norb,nsite,nsite,nkpt,nspin,nspin))
          ham_ksoc = cmplx_0
          do ikpt = 1, nkpt
             open(2, file=hamname, status='old')
             hamsoc = cmplx_0
             nrpt = 0
             do while (io_stat == 0)
                read(2,*,iostat=io_stat) iat, iat1, norb, norb
                if (io_stat < 0) cycle
                read(2,*,iostat=io_stat) bond(:)

                if (bondcoord == 'f') then
                   call utility_frac_to_cart(bond(:),bond(:),lattvec)
                end if

                do iorb = 1, norb
                   read(2,*,iostat=io_stat) hamsoc(iorb,:,1,1)
                end do

                nrpt = nrpt + 1
                !Expand to spindown subspace
                hamsoc(:,:,2,2) = hamsoc(:,:,1,1)
                if (nrpt <= nsite) then
                   hamsoc(:,:,1,1) = hamsoc(:,:,1,1) + vso(:,:,nrpt,1,1)
                   hamsoc(:,:,1,2) = hamsoc(:,:,1,2) + vso(:,:,nrpt,1,2)
                   hamsoc(:,:,2,1) = hamsoc(:,:,2,1) + vso(:,:,nrpt,2,1)
                   hamsoc(:,:,2,2) = hamsoc(:,:,2,2) + vso(:,:,nrpt,2,2)
                end if
                rdotk = dot_product(kpt_cart(:,ikpt),bond(:))
                fac = exp(cmplx_i*rdotk)
                ham_ksoc(:,:,iat,iat1,ikpt,1,1) = ham_ksoc(:,:,iat,iat1,ikpt,1,1) + fac*hamsoc(:,:,1,1)
                ham_ksoc(:,:,iat,iat1,ikpt,2,2) = ham_ksoc(:,:,iat,iat1,ikpt,2,2) + fac*hamsoc(:,:,2,2)
                ham_ksoc(:,:,iat,iat1,ikpt,1,2) = ham_ksoc(:,:,iat,iat1,ikpt,1,2) + fac*hamsoc(:,:,1,2)
                ham_ksoc(:,:,iat,iat1,ikpt,2,1) = ham_ksoc(:,:,iat,iat1,ikpt,2,1) + fac*hamsoc(:,:,2,1)
             end do
             io_stat = 0
             close(2)
          end do

          allocate(ham_rsoc(norb,norb,nrpt,nspin,nspin))
          allocate(r_bond(3,nrpt))
          allocate(atom_label(2,nrpt))

          open(2, file=hamname, status='old')
          ham_rsoc = cmplx_0
          do irpt = 1, nrpt
             read(2,*) atom_label(1,irpt), atom_label(2,irpt)
             read(2,*) r_bond(:,irpt) 

             if (bondcoord == 'f') then
                call utility_frac_to_cart(r_bond(:,irpt),r_bond(:,irpt),lattvec)
             end if

             do iorb = 1, norb
                read(2,*) ham_rsoc(iorb,:,irpt,1,1)
             end do
             ham_rsoc(:,:,irpt,2,2) = ham_rsoc(:,:,irpt,1,1)
             if (irpt <= nsite) then
                ham_rsoc(:,:,irpt,1,1) = ham_rsoc(:,:,irpt,1,1) + vso(:,:,irpt,1,1)
                ham_rsoc(:,:,irpt,1,2) = ham_rsoc(:,:,irpt,1,2) + vso(:,:,irpt,1,2)
                ham_rsoc(:,:,irpt,2,1) = ham_rsoc(:,:,irpt,2,1) + vso(:,:,irpt,2,1)
                ham_rsoc(:,:,irpt,2,2) = ham_rsoc(:,:,irpt,2,2) + vso(:,:,irpt,2,2)
             end if
          end do
          close(2)
       else
          print*, 'No kinetic hamiltonian file found'
          stop
       end if

    end if

  end subroutine get_kinetic_parameters

  subroutine get_pauli()

    use constants,           only : cmplx_0, cmplx_i, nspin

    integer   :: iorb, ispin, ispin1

    allocate(sigma_x(norb,norb,nspin,nspin))
    allocate(sigma_y(norb,norb,nspin,nspin))
    allocate(sigma_z(norb,norb,nspin,nspin))

    do iorb = 1, norb
       do ispin = 1, nspin
          do ispin1 = 1, nspin
             if (ispin == 1 .and. ispin1 == 1) then
                sigma_x(iorb,iorb,ispin,ispin1) = cmplx_0
                sigma_y(iorb,iorb,ispin,ispin1) = cmplx_0
                sigma_z(iorb,iorb,ispin,ispin1) = 1.0
             else if (ispin == 2 .and. ispin1 == 2) then
                sigma_x(iorb,iorb,ispin,ispin1) = cmplx_0
                sigma_y(iorb,iorb,ispin,ispin1) = cmplx_0
                sigma_z(iorb,iorb,ispin,ispin1) = -1.0
             else if (ispin == 1 .and. ispin1 == 2) then
                sigma_x(iorb,iorb,ispin,ispin1) = 1.0
                sigma_y(iorb,iorb,ispin,ispin1) = -1.0*cmplx_i
                sigma_z(iorb,iorb,ispin,ispin1) = cmplx_0
             else if (ispin == 2 .and. ispin1 == 1) then
                sigma_x(iorb,iorb,ispin,ispin1) = 1.0
                sigma_y(iorb,iorb,ispin,ispin1) = cmplx_i
                sigma_z(iorb,iorb,ispin,ispin1) = cmplx_0
             end if
          end do
       end do
    end do

  end subroutine get_pauli  

  subroutine get_coulomb_matrix()

  use constants,            only : nspin, coulombname, cmplx_0

  integer                   :: isite, ispin, i, j, k, l, itot
  complex*16  , allocatable :: u_tmp(:,:,:,:)
  real*8                    :: u_real, u_imag
  logical                   :: io_find

  allocate(u_tmp(ntot/2,ntot/2,ntot/2,ntot/2))
  allocate(u_mat(norb,norb,norb,norb,nsite))

  ! Read screened.dat provisionally
  u_mat = cmplx_0
  u_tmp = cmplx_0

  inquire(file=coulombname, exist=io_find)
  if (io_find) then
     open(unit=1, file=coulombname, status='old')
        do itot = 1, ntot*ntot*ntot*ntot/(2*2*2*2)
           read(1,*) i, j, k, l, u_real, u_imag
           u_tmp(i,j,k,l) = cmplx(u_real,u_imag)
        end do
     close(1)
  end if

  do isite = 1, nsite
     u_mat(:,:,:,:,isite) =  u_tmp((isite-1)*norb+1:isite*norb, &
                                   (isite-1)*norb+1:isite*norb, &
                                   (isite-1)*norb+1:isite*norb, &
                                   (isite-1)*norb+1:isite*norb)
  end do  

  deallocate(u_tmp)

  end subroutine get_coulomb_matrix

  subroutine spherical_param()

  use constants,            only : cmplx_0

  complex*16, allocatable   :: usph(:), jsph(:)
  integer                   :: isite, iorb, iorb1
 
  allocate(u_mat(norb,norb,norb,norb,nsite)) 
  allocate(usph(nsite))
  allocate(jsph(nsite))

  u_mat = cmplx_0
  do isite = 1, nsite
     do iorb = 1, norb
       do iorb1 = 1, norb
           u_mat(iorb,iorb,iorb1,iorb1,isite) = 2.9 !usph(isite)
        end do
     end do
  end do

  do isite = 1, nsite
     do iorb = 1, norb
       do iorb1 = 1, norb
           if (iorb /= iorb1) then
              u_mat(iorb,iorb1,iorb1,iorb,isite) = 0.9 !jsph(isite)
           end if
        end do
     end do
  end do

!  u_mat = cmplx_0
!  u_mat(1,1,1,1,1) = 0.707
!  u_mat(1,1,2,2,1) = 0.509
!  u_mat(1,1,3,3,1) = 0.509
!  u_mat(2,2,1,1,1) = 0.509
!  u_mat(2,2,2,2,1) = 0.676
!  u_mat(2,2,3,3,1) = 0.505
!  u_mat(3,3,1,1,1) = 0.509
!  u_mat(3,3,2,2,1) = 0.505
!  u_mat(3,3,3,3,1) = 0.676

!  u_mat(1,1,1,1,1) = 0.707
!  u_mat(1,2,2,1,1) = 0.086
!  u_mat(1,3,3,1,1) = 0.086
!  u_mat(2,1,1,2,1) = 0.086
!  u_mat(2,2,2,2,1) = 0.676
!  u_mat(2,3,3,2,1) = 0.086
!  u_mat(3,1,1,3,1) = 0.086
!  u_mat(3,2,2,3,1) = 0.086
!  u_mat(3,3,3,3,1) = 0.676

!  u_mat(1,1,1,1,1) = 0.707
!  u_mat(1,2,1,2,1) = 0.086
!  u_mat(1,3,1,3,1) = 0.086
!  u_mat(2,1,2,1,1) = 0.086
!  u_mat(2,2,2,2,1) = 0.676
!  u_mat(2,3,2,3,1) = 0.086
!  u_mat(3,1,3,1,1) = 0.086
!  u_mat(3,2,3,2,1) = 0.086
!  u_mat(3,3,3,3,1) = 0.676


  end subroutine spherical_param

end module param
