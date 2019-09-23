program main 
  use utilities
  use param
  use constants
  use dens
  use hamilton
  use energies
  use plots
 
  implicit none

  real*8                         :: etot_tmp, occ, smearing, tmp
  real*8           , allocatable :: hf_valsoc(:,:), hf_val(:,:,:), val_efer(:,:)
  complex*16       , allocatable :: hf_vecsoc(:,:,:), hf_vec(:,:,:,:), temp(:,:), tempsoc(:,:), temp_denssoc(:,:,:,:,:),  temp_dens(:,:,:,:)
  integer                        :: iscf, itot,j, itot1, itot2, itot3, ikpt, coun, iorb, iorb1, iorb2, iorb3, irpt, ispin, ispin1, isite, i 
 
  open(33, file = 'hf.out', status = 'replace')
     write(33,'(a)') '!-----------------------------------------------------------------------------!'
     write(33,'(a)') '!                                                                             !'                   
     write(33,'(a)') '!            MM     MM   FFFFFFFFF          H       H   FFFFFFFFF             !'
     write(33,'(a)') '!            M M   M M   F                  H       H   F                     !'
     write(33,'(a)') '!            M  M M  M   F                  H       H   F                     !'
     write(33,'(a)') '!            M   M   M   FFFFFF     ---     HHHHHHHHH   FFFFFF                !'
     write(33,'(a)') '!            M       M   F                  H       H   F                     !'
     write(33,'(a)') '!            M       M   F                  H       H   F                     !'
     write(33,'(a)') '!            M       M   F                  H       H   F                     !'
     write(33,'(a)') '!                                                                             !'                        
     write(33,'(a)') '!           Mean-Field Hartree-Fock solver by Oier Arcelus (v.0.0)            !'
     write(33,'(a)') '!                                                                             !' 
     write(33,'(a)') '!----------------------------- CIC Energigune --------------------------------!'
     write(33,'(a)') ' '


     write(33,'(a)') 'Reading input file...'
     write(33,'(a)') ' '
     ! Read input file
     call get_names()
     write(33,'(a,L)') '  spinorb = ', spinorb
     write(33,'(a,i3)') '  nelec = ', nelec
     write(33,'(a,i4)') '  nscf = ', nscf
     write(33,'(a,i4)') '  nener = ', nener
     write(33,'(a,F7.3)') '  etemp = ', etemp
     write(33,'(a,ES12.3,a)') '  ediff = ', ediff, ' eV'
     write(33,'(a,F12.6)') '  mix = ', mix
     write(33,'(a,L)') '  dosplot = ', dosplot
     write(33,'(a,F7.3)') '  smear = ', smear
     write(33,'(a,L)') '  bandplot = ', bandplot
     write(33,'(a,L)') '  lsph = ', lsph
     write(33,'(a)') ' '
     write(33,'(a)') 'Input file read successfuly'


     write(33,'(a)') 'Reading lattice parameters...'
     write(33,'(a)') ' '
     ! Read lattice parameters
     call get_latt_parameters()
     write(33,'(i4,a)') nsite, ' sites found'
     write(33,'(a)') 'Lattice vectors a, b and c in units of Amstrong'
     do i = 1, 3
        write(33,'(a,3F12.6)') '  ',lattvec(:,i)
     end do
     write(33,'(a)') ' '
     write(33,'(a)') 'Site positions in direct coordinates'
     do isite = 1, nsite  
        write(33,'(a,3F12.6)') '  ',atom_dir(:,isite)
     end do
     write(33,'(a)') 'Site positions in cartesian coordinates (Amstrongs)'
     do isite = 1, nsite
        write(33,'(a,3F12.6)') '  ',atom_cart(:,isite)
     end do     
     write(33,'(a)') ' '
     write(33,'(a)') 'Lattice parameters read successfuly'


     write(33,'(a)') 'Reading k-points...'
     write(33,'(a)') ' '
     ! Read Kpoints
     call get_kpt_parameters()
     write(33,'(i4,a)') nkpt, ' k-points found'
     write(33,'(a)') 'k-points in direct and cartesian coordinates'
     write(33,'(a)') '                 direct                                cartesian '
     write(33,'(a)') '    --------------------------------        --------------------------------'
     do ikpt = 1, nkpt
        write(33,'(3F12.6,4x,3F12.6)') kpt_latt(:,ikpt), kpt_cart(:,ikpt)
     end do
     write(33,'(a)') ' '
     write(33,'(a)') 'k-points read successfuly'


     write(33,'(a)') 'Reading hamiltonian in wannier basis from Wannier90...'
     write(33,'(a)') ' '
     ! Read kinetic parameters
     call get_kinetic_parameters()
     write(33,'(a,i3)') '  ntot = ', ntot
     write(33,'(a,i3)') '  nspin = ', nspin
     write(33,'(a,i3)') '  norb = ', norb
     write(33,'(a)') ' '
     write(33,'(a)') 'Hamiltonian read successfuly'
     write(33,'(a)') 'Reading Coulomb parameters...'
     write(33,'(a)') ' '
     call get_pauli() ! Get coulomb matrix
     if (.not.lsph) call get_coulomb_matrix()
     if (lsph) then
        write(33,'(a)') 'Spherical parameterization of coulomb matrices in eV'
        write(33,'(a)') ' '
        call spherical_param()
     end if         

     write(33,'(a)') 'Reduced Coulomb matrices in eV'
     do isite = 1, nsite
        write(33,'(a,i2)') '  site: ', isite
        write(33,'(a)') '  U_IIJJ'
        do iorb = 1, norb
           write(33,'(<norb>F10.3)') (real(u_mat(iorb1,iorb1,iorb,iorb,isite)), iorb1 = 1, norb)
        end do 
        write(33,'(a)') '  U_IJJI'
        do iorb = 1, norb
           write(33,'(<norb>F10.3)') (real(u_mat(iorb1,iorb,iorb,iorb1,isite)), iorb1 = 1, norb)
        end do
        write(33,'(a)') '  U_IJIJ'
        do iorb = 1, norb
           write(33,'(<norb>F10.3)') (real(u_mat(iorb1,iorb,iorb1,iorb,isite)), iorb1 = 1, norb)
        end do
     end do
     write(33,'(a)') ' '
     write(33,'(a)') 'Coulomb parameters read successfuly'
     write(33,'(a)') ' '
     write(33,'(a)') 'Initializing density matrix...'
     write(33,'(a)') ' '
     
     call initialize_density_matrix()

     if (.not.spinorb) then

     do ispin = 1, nspin
        write(33,*) 'spin = ', ispin
        do isite = 1, nsite
           write(33,*) 'site = ', isite
           write(33,'(a)') ' '
           do iorb = 1, norb
              write(33,'(<norb>F7.3)') (real(init_dens(iorb,iorb1,isite,ispin)), iorb1 = 1, norb)
           end do
           write(33,'(a)') ' '
        end do
     end do

     else if (spinorb) then

     do ispin = 1, nspin
        do ispin1 = 1, nspin
           write(33,*) 'spin1 = ', ispin, ' spin2 = ', ispin1
           do isite = 1, nsite
              write(33,*) 'site = ', isite
              write(33,'(a)') ' '
              do iorb = 1, norb
                 write(33,'(<norb>F7.3)') (real(init_denssoc(iorb,iorb1,isite,ispin,ispin1)), iorb1 = 1, norb)
              end do
              write(33,'(a)') ' '
           end do
        end do
     end do

     end if

     write(33,'(a)') ' '
     write(33,'(a)') 'Done'

     if (.not.spinorb) then
        allocate(hf_val(norb*nsite,nspin,nkpt))
        allocate(hf_vec(norb*nsite,norb*nsite,nspin,nkpt))
        allocate(val_efer(ntot,nkpt))
        allocate(temp(norb*nsite,norb*nsite))
        allocate(temp_dens(norb,norb,nsite,nspin))
        temp_dens = init_dens
     else if (spinorb) then
        allocate(hf_valsoc(ntot,nkpt))
        allocate(hf_vecsoc(ntot,ntot,nkpt))
        allocate(tempsoc(ntot,ntot))
        allocate(temp_denssoc(norb,norb,nsite,nspin,nspin))
        temp_denssoc = init_denssoc
     end if

     write(33,'(a)') 'Initializing HF potential from starting density ...'
     write(33,'(a)') ' '

     call init_hf_pot()

     if (.not.spinorb) then

     do ispin = 1, nspin
        write(33,*) 'spin = ', ispin
        do isite = 1, nsite
           write(33,*) 'site = ', isite
           write(33,'(a)') ' '
           do iorb = 1, norb
              write(33,'(<norb>F7.3)') (real(hf_pot(iorb,iorb1,isite,ispin)), iorb1 = 1, norb)
           end do
           write(33,'(a)') ' '
        end do
     end do

     else if (spinorb) then

     do ispin = 1, nspin
        do ispin1 = 1, nspin
           write(33,*) 'spin1 = ', ispin, 'spin2 = ', ispin1
           do isite = 1, nsite
              write(33,*) 'site = ', isite
              write(33,'(a)') ' '
              do iorb = 1, norb
                 write(33,'(<norb>F7.3)') (real(hf_potsoc(iorb,iorb1,isite,ispin,ispin1)), iorb1 = 1, norb)
              end do
              write(33,'(a)') ' '
           end do
        end do
     end do

     end if

     call set_full_ham()

     write(33,'(a)') 'HF potential successfuly initialized'
     write(33,'(a)') 'Starting self-consistent cycle...'
     write(33,'(a)') ' '

     smearing = etemp

     if (.not.spinorb) then

     etot_tmp = 0.0
     ! Start self-consistent cycle ####################################################################################
     do iscf = 1, nscf
        write(33,'(a,i3,a)') ' --------------------- Iteration: ', iscf, ' -------------------------'
        write(33,'(a)') ' '

        ! Diagonalize hamiltonian
        do ikpt = 1, nkpt 
           do ispin = 1, nspin
              temp(:,:) = ham_full(:,:,ispin,ikpt)
              call utility_diag(temp(:,:),hf_val(:,ispin,ikpt),norb*nsite)
              hf_vec(:,:,ispin,ikpt) = temp(:,:)
           end do
        end do

        ! Order hf_val in ascending order, save in val_efer.
        do ikpt = 1, nkpt
           coun = 0 
           do ispin = 1, nspin
              do itot = 1, norb*nsite
                 coun = coun + 1
                 val_efer(coun,ikpt) = hf_val(itot,ispin,ikpt)
              end do
           end do
        end do
        do ikpt = 1, nkpt
           do i = 1, nsite*norb*nspin - 1
              do j = i + 1, nsite*norb*nspin
                 if (val_efer(i,ikpt) > val_efer(j,ikpt)) then
                    ! Swap distance info
                    tmp = val_efer(i,ikpt)
                    val_efer(i,ikpt) = val_efer(j,ikpt)
                    val_efer(j,ikpt) = tmp
                 end if
              end do
           end do
        end do

        call set_fermi_level(val_efer,smearing)
        write(33,'(a,F12.6,a,a,F12.6)') '    Fermi energy : ', efermi,' eV', ' Number elec. :', num
        write(33,'(a)') ' '

        call set_density_matrix(hf_vec,hf_val,smearing)
        write(33,'(a)') 'Mixing old and new density matrix'
        dens_mat = (1.0 - mix)*temp_dens + mix*dens_mat
        temp_dens = dens_mat
        write(33,'(a)') 'On-site density matrix'
        do ispin = 1, nspin
           write(33,*) 'spin = ', ispin
           do isite = 1, nsite
              write(33,*) 'site = ', isite
              write(33,'(a)') ' '
              do iorb = 1, norb
                 write(33,'(<norb>F7.3)') (real(dens_mat(iorb,iorb1,isite,ispin)), iorb1 = 1, norb)
              end do
              write(33,'(a)') ' '
           end do
        end do

        call set_hf_pot(dens_mat) 

        call set_full_ham()

        call evaluate_energy(hf_val,hf_pot,dens_mat)

        write(33,'(a)') ' Calculating HF energy...'
        write(33,'(a)') ' '
        write(33,'(a)') '    E (eV)      dE (eV)'
        write(33,'(ES12.4,E12.4)') etot, etot - etot_tmp
        write(33,'(a)') ' '
       ! Convergence condition

        if (abs(etot - etot_tmp) < ediff) then
           write(33,'(a)') '---------------- Reached required accuracy   ----------------'
           write(33,'(a)') '---------------- Ending self-consistent loop ----------------'
           write(33,'(a)') ' '
           write(33,'(a,F15.8,a)') '!   Converged HF energy : ', etot, ' eV'
           exit
        else
           etot_tmp = etot
        end if
     end do 

     if (dosplot) then

        call plot_dos(hf_val,hf_vec,efermi)    
  
     end if

     if (bandplot) then
 
        call plot_band(efermi)
 
     end if

     call get_magnetization(dens_mat)
     write(33,'(a)') ' '
     write(33,'(a)') 'Occupations'
     write(33,'(a)') '------------------------------'
     do isite = 1, nsite
        occ = real(trace(dens_mat(:,:,isite,1),norb)) + real(trace(dens_mat(:,:,isite,2),norb))
        write(33,'(a,i2,3F10.4)') 'Atom: ', isite,  occ
     end do

     write(33,'(a)') ' '
     write(33,'(a)') 'Magnetization (bohr magnetons)'
     write(33,'(a)') '------------------------------'
     do isite = 1, nsite
        write(33,'(a,i2,F10.4)') 'Atom: ', isite, real(mag(isite))
     end do 

     open(unit=14, file='wav.out', status='replace')
        write(14,*) nkpt
        write(14,*) norb*nsite
        do ikpt = 1, nkpt
           write(14,*) kpt_cart(:,ikpt)
           do ispin = 1, nspin
              write(14,*) ispin
              do itot = 1, norb*nsite
                 do itot1 = 1, norb*nsite
                    write(14,*) itot1, itot, real(hf_vec(itot1,itot,ispin,ikpt)), aimag(hf_vec(itot1,itot,ispin,ikpt))
                 end do
              end do
              write(14,*) ' '
           end do
        end do
     close(14)
     open(unit=15, file='eig.out', status='replace')
        write(15,*) efermi
        do ikpt = 1, nkpt
           write(15,*) kpt_cart(:,ikpt)
           do ispin = 1, nspin
              write(15,*) ispin
              do itot = 1, norb*nsite
                 write(15,*) itot, hf_val(itot,ispin,ikpt)
              end do
              write(15,*) ' '
           end do
        end do
     close(15)
     open(unit=16, file='bond.out', status='replace')
        write(16,*) nrpt
        do irpt = 1, nrpt
           write(16,*) atom_label(:,irpt)
           write(16,*) r_bond(:,irpt)
        end do
     close(16)
     open(unit=17, file='pot.out', status='replace')
        do ispin = 1, nspin
           do isite = 1, nsite
              write(17,*) ispin, isite
              do iorb = 1, norb
                 do iorb1 = 1, norb
                    write(17,*) iorb, iorb1, real(hf_pot(iorb,iorb1,isite,ispin)), &
                                aimag(hf_pot(iorb,iorb1,isite,ispin))
                 end do
              end do
              write(17,'(a)') ' '
           end do
        end do
     close(17)
     open(unit=17, file='dens.out', status='replace')
        do ispin = 1, nspin
           do isite = 1, nsite
              write(17,*) ispin, isite
              do iorb = 1, norb
                 do iorb1 = 1, norb
                    write(17,*) iorb, iorb1, real(dens_mat(iorb,iorb1,isite,ispin)), &
                                aimag(dens_mat(iorb,iorb1,isite,ispin))
                 end do
              end do
              write(17,'(a)') ' '
           end do
        end do
     close(17)
   close(33)

   else if (spinorb) then

     etot_tmp = 0.0
     ! Start self-consistent cycle
     ! ####################################################################################
     do iscf = 1, nscf
        write(33,'(a,i3,a)') ' --------------------- Iteration: ', iscf, ' -------------------------'
        write(33,'(a)') ' '

        ! Diagonalize hamiltonian
        do ikpt = 1, nkpt
           tempsoc(:,:) = ham_fullsoc(:,:,ikpt)
           call utility_diag(tempsoc(:,:),hf_valsoc(:,ikpt),ntot)
           hf_vecsoc(:,:,ikpt) = tempsoc(:,:)
        end do

        call set_fermi_level(hf_valsoc,smearing)
        write(33,'(a,F12.6,a,a,F12.6)') '    Fermi energy : ', efermi,' eV', ' Number elec. :', num
        write(33,'(a)') ' '

        call set_density_matrix_soc(hf_vecsoc,hf_valsoc,smearing)
        write(33,'(a)') 'Mixing old and new density matrix'
        dens_matsoc = (1.0 - mix)*temp_denssoc + mix*dens_matsoc
        temp_denssoc = dens_matsoc
        write(33,'(a)') 'On-site density matrix'
        do ispin = 1, nspin
           do ispin1 = 1, nspin
              write(33,*) 'spin1 = ', ispin, ' spin2 = ', ispin1
              do isite = 1, nsite
                 write(33,*) 'site = ', isite
                 write(33,'(a)') ' '
                 do iorb = 1, norb
                    write(33,'(<2*norb>F7.3)') (dens_matsoc(iorb,iorb1,isite,ispin,ispin1), iorb1 = 1, norb)
                 end do
                 write(33,'(a)') ' '
              end do
           end do
        end do
        call set_hf_pot_soc(dens_matsoc)
        call set_full_ham()
        call evaluate_energy_soc(hf_valsoc,hf_potsoc,dens_matsoc)

        write(33,'(a)') ' Calculating HF energy...'
        write(33,'(a)') ' '
        write(33,'(a)') '    E (eV)      dE (eV)'
        write(33,'(F20.10,F20.10)') etot, etot - etot_tmp
        write(33,'(a)') ' '
       ! Convergence condition

        if (abs(etot - etot_tmp) < ediff) then
           write(33,'(a)') '---------------- Reached required accuracy ----------------'
           write(33,'(a)') '---------------- Ending self-consistent loop ----------------'
           write(33,'(a)') ' '
           write(33,'(a,F20.10,a)') '!   Converged HF energy : ', etot, ' eV'
           exit
        else
           etot_tmp = etot
        end if
     end do

     if (dosplot) then

        write(33,'(a)') ' '
        write(33,'(a)') ' dosplot set to .true., not implemented for spinorb = .true. '

     end if

     if (bandplot) then

        call plot_band(efermi)

     end if

     call get_magnetization_soc(dens_matsoc)
     write(33,'(a)') ' '
     write(33,'(a)') 'Occupations'
     write(33,'(a)') '------------------------------'
     do isite = 1, nsite
        occ = real(trace(dens_matsoc(:,:,isite,1,1),norb)) + real(trace(dens_matsoc(:,:,isite,2,2),norb))
        write(33,'(a,i2,3F10.4)') 'Atom: ', isite,  occ
     end do

     write(33,'(a)') ' '
     write(33,'(a)') 'Magnetization (bohr magnetons)'
     write(33,'(a)') '------------------------------'
     do isite = 1, nsite
        write(33,'(a,i2,3F10.4)') 'Atom: ', isite, real(magsoc(:,isite))
     end do

     open(unit=14, file='wav.out', status='replace')
        write(14,*) nkpt
        write(14,*) ntot
        do ikpt = 1, nkpt
           write(14,*) kpt_cart(:,ikpt)
           do itot = 1, ntot
              do itot1 = 1, ntot
                 write(14,*) itot1, itot, real(hf_vecsoc(itot1,itot,ikpt)), aimag(hf_vecsoc(itot1,itot,ikpt))
              end do
           end do
           write(14,*) ' '
        end do
     close(14)
     open(unit=15, file='eig.out', status='replace')
        write(15,*) efermi
        do ikpt = 1, nkpt
           write(15,*) kpt_cart(:,ikpt)
           do itot = 1, ntot
              write(15,*) itot, hf_valsoc(itot,ikpt)
           end do
           write(15,*) ' '
        end do
     close(15)
     open(unit=16, file='bond.out', status='replace')
        write(16,*) nrpt
        do irpt = 1, nrpt
           write(16,*) atom_label(:,irpt)
           write(16,*) r_bond(:,irpt)
        end do
     close(16)
     open(unit=17, file='pot.out', status='replace')
        do ispin = 1, nspin
           do ispin1 = 1, nspin
              do isite = 1, nsite
                 write(17,*) ispin, ispin1, isite
                 do iorb = 1, norb
                    do iorb1 = 1, norb
                       write(17,*) iorb, iorb1, real(hf_potsoc(iorb,iorb1,isite,ispin,ispin1)), &
                                   aimag(hf_potsoc(iorb,iorb1,isite,ispin,ispin1))
                    end do
                 end do
                 write(17,'(a)') ' '
              end do
           end do
        end do
     close(17)
     open(unit=17, file='dens.out', status='replace')
        do ispin = 1, nspin
           do ispin1 = 1, nspin
              do isite = 1, nsite
                 write(17,*) ispin, ispin1, isite
                 do iorb = 1, norb
                    do iorb1 = 1, norb
                       write(17,*) iorb, iorb1, real(dens_matsoc(iorb,iorb1,isite,ispin,ispin1)), &
                                   aimag(dens_matsoc(iorb,iorb1,isite,ispin,ispin1))
                    end do
                 end do
                 write(17,'(a)') ' '
              end do
           end do
        end do
     close(17)

   close(33)

   end if

end program main
