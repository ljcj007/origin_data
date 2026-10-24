module sse_square_octagon_model
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: nobs = 8
  real(dp), parameter :: coupling = 1.0_dp

  integer :: lx, ly, ncell, nsite, nbond
  integer :: nb_square, nb_link
  integer :: norder, max_order
  real(dp) :: beta

  integer, allocatable :: spin(:), eta(:)
  integer, allocatable :: bonds(:,:), bond_type(:)
  integer, allocatable :: opstring(:)
  integer, allocatable :: first_op(:), last_op(:), vertex(:)

contains

  subroutine read_input(nbins, msteps, isteps, seed, input_file)
    integer, intent(out) :: nbins, msteps, isteps, seed
    character(len=*), intent(in), optional :: input_file
    integer :: unit_id, ios
    character(len=256) :: path

    path = 'sse_square_octagon.in'
    if (present(input_file)) path = input_file

    open(newunit=unit_id, file=trim(path), status='old', &
         action='read', iostat=ios)
    if (ios /= 0) error stop '无法打开 sse_square_octagon.in'

    read(unit_id, *, iostat=ios) lx, ly, beta
    if (ios /= 0) error stop '输入文件第一行格式错误'
    read(unit_id, *, iostat=ios) nbins, msteps, isteps
    if (ios /= 0) error stop '输入文件第二行格式错误'
    read(unit_id, *, iostat=ios) seed
    if (ios /= 0) error stop '输入文件第三行格式错误'
    close(unit_id)

    if (lx < 2 .or. ly < 2) error stop 'Lx 和 Ly 必须至少为 2'
    if (mod(lx, 2) /= 0 .or. mod(ly, 2) /= 0) &
      error stop '周期边界下 Lx 和 Ly 必须为偶数'
    if (beta <= 0.0_dp) error stop 'beta 必须为正数'
    if (nbins < 2 .or. msteps < 1 .or. isteps < 1) &
      error stop 'bin 数至少为 2，步数必须为正数'
    if (seed <= 0) error stop '随机种子必须为正整数'
  end subroutine read_input

  subroutine init_rng(seed)
    integer, intent(in) :: seed
    integer :: i, nseed
    integer, allocatable :: put(:)

    call random_seed(size=nseed)
    allocate(put(nseed))
    do i = 1, nseed
      put(i) = modulo(seed + 104729*i + 37*i*i, huge(1) - 1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
    deallocate(put)
  end subroutine init_rng

  real(dp) function rand01()
    call random_number(rand01)
  end function rand01

  integer function cell_index(x, y)
    integer, intent(in) :: x, y
    cell_index = 1 + modulo(x, lx) + modulo(y, ly)*lx
  end function cell_index

  integer function site_index(cell, local_site)
    integer, intent(in) :: cell, local_site
    site_index = 4*(cell - 1) + local_site
  end function site_index

  subroutine add_bond(b, site1, site2, kind)
    integer, intent(inout) :: b
    integer, intent(in) :: site1, site2, kind

    b = b + 1
    bonds(1, b) = site1
    bonds(2, b) = site2
    bond_type(b) = kind
  end subroutine add_bond

  subroutine make_lattice()
    integer :: x, y, c, cr, cu, a, b

    ncell = lx*ly
    nsite = 4*ncell
    nb_square = 4*ncell
    nb_link = 2*ncell
    nbond = nb_square + nb_link

    allocate(bonds(2, nbond), bond_type(nbond), eta(nsite))
    b = 0

    do y = 0, ly - 1
      do x = 0, lx - 1
        c = cell_index(x, y)
        cr = cell_index(x + 1, y)
        cu = cell_index(x, y + 1)

        call add_bond(b, site_index(c, 1), site_index(c, 2), 1)
        call add_bond(b, site_index(c, 2), site_index(c, 3), 1)
        call add_bond(b, site_index(c, 3), site_index(c, 4), 1)
        call add_bond(b, site_index(c, 4), site_index(c, 1), 1)
        call add_bond(b, site_index(c, 2), site_index(cr, 4), 2)
        call add_bond(b, site_index(c, 3), site_index(cu, 1), 2)

        do a = 1, 4
          if (mod(a + x + y, 2) == 0) then
            eta(site_index(c, a)) = 1
          else
            eta(site_index(c, a)) = -1
          end if
        end do
      end do
    end do

    if (b /= nbond) error stop '建格时键数不一致'
  end subroutine make_lattice

  subroutine validate_lattice()
    integer :: b1, b2, i
    integer, allocatable :: degree(:)
    integer :: a1, a2, c1, c2

    allocate(degree(nsite))
    degree = 0

    do b1 = 1, nbond
      a1 = bonds(1, b1)
      a2 = bonds(2, b1)
      if (a1 < 1 .or. a1 > nsite .or. a2 < 1 .or. a2 > nsite) &
        error stop '存在越界格点编号'
      if (a1 == a2) error stop '存在自连接键'
      if (eta(a1)*eta(a2) /= -1) error stop '晶格不是二分图'
      degree(a1) = degree(a1) + 1
      degree(a2) = degree(a2) + 1

      do b2 = b1 + 1, nbond
        c1 = bonds(1, b2)
        c2 = bonds(2, b2)
        if (min(a1, a2) == min(c1, c2) .and. &
            max(a1, a2) == max(c1, c2)) error stop '存在重复键'
      end do
    end do

    do i = 1, nsite
      if (degree(i) /= 3) error stop '存在配位数不等于 3 的格点'
    end do
    deallocate(degree)

    write(*, '(a,3(i0,a))') '晶格检查通过：N=', nsite, ', Nb=', nbond, &
      ', 晶胞=', ncell, ''
  end subroutine validate_lattice

  subroutine init_configuration()
    allocate(spin(nsite), first_op(nsite), last_op(nsite))
    spin = eta

    max_order = max(64, nsite)
    allocate(opstring(0:max_order - 1))
    allocate(vertex(0:4*max_order - 1))
    opstring = 0
    norder = 0
  end subroutine init_configuration

  subroutine diagonal_update()
    integer :: p, op, b, s1, s2
    real(dp) :: probability, weight

    weight = 0.5_dp*beta*coupling*real(nbond, dp)

    do p = 0, max_order - 1
      op = opstring(p)
      if (op == 0) then
        if (norder == max_order) cycle
        b = 1 + int(rand01()*real(nbond, dp))
        s1 = bonds(1, b)
        s2 = bonds(2, b)
        if (spin(s1) /= spin(s2)) then
          probability = weight/real(max_order - norder, dp)
          if (rand01() <= min(1.0_dp, probability)) then
            opstring(p) = 2*b
            norder = norder + 1
          end if
        end if
      else if (mod(op, 2) == 0) then
        probability = real(max_order - norder + 1, dp)/weight
        if (rand01() <= min(1.0_dp, probability)) then
          opstring(p) = 0
          norder = norder - 1
        end if
      else
        b = op/2
        s1 = bonds(1, b)
        s2 = bonds(2, b)
        spin(s1) = -spin(s1)
        spin(s2) = -spin(s2)
      end if
    end do
  end subroutine diagonal_update

  subroutine loop_update()
    integer :: i, op, b, s1, s2
    integer :: v0, v1, v2

    first_op = -1
    last_op = -1

    do v0 = 0, 4*max_order - 1, 4
      op = opstring(v0/4)
      if (op == 0) then
        vertex(v0:v0 + 3) = -1
        cycle
      end if

      b = op/2
      s1 = bonds(1, b)
      s2 = bonds(2, b)
      v1 = last_op(s1)
      v2 = last_op(s2)

      if (v1 == -1) then
        first_op(s1) = v0
      else
        vertex(v1) = v0
        vertex(v0) = v1
      end if

      if (v2 == -1) then
        first_op(s2) = v0 + 1
      else
        vertex(v2) = v0 + 1
        vertex(v0 + 1) = v2
      end if

      last_op(s1) = v0 + 2
      last_op(s2) = v0 + 3
    end do

    do i = 1, nsite
      if (first_op(i) /= -1) then
        vertex(last_op(i)) = first_op(i)
        vertex(first_op(i)) = last_op(i)
      end if
    end do

    do v0 = 0, 4*max_order - 1, 2
      if (vertex(v0) < 0) cycle
      v1 = v0

      if (rand01() < 0.5_dp) then
        do
          opstring(v1/4) = ieor(opstring(v1/4), 1)
          vertex(v1) = -2
          v2 = ieor(v1, 1)
          v1 = vertex(v2)
          vertex(v2) = -2
          if (v1 == v0) exit
        end do
      else
        do
          vertex(v1) = -1
          v2 = ieor(v1, 1)
          v1 = vertex(v2)
          vertex(v2) = -1
          if (v1 == v0) exit
        end do
      end if
    end do

    do i = 1, nsite
      if (first_op(i) == -1) then
        if (rand01() < 0.5_dp) spin(i) = -spin(i)
      else if (vertex(first_op(i)) == -2) then
        spin(i) = -spin(i)
      end if
    end do
  end subroutine loop_update

  subroutine adjust_cutoff()
    integer :: new_order
    integer, allocatable :: old_string(:)

    new_order = norder + norder/3 + 16
    if (new_order <= max_order) return
    new_order = max(new_order, max_order + 16)

    allocate(old_string(0:max_order - 1))
    old_string = opstring
    deallocate(opstring, vertex)

    allocate(opstring(0:new_order - 1))
    allocate(vertex(0:4*new_order - 1))
    opstring = 0
    opstring(0:max_order - 1) = old_string
    deallocate(old_string)
    max_order = new_order
  end subroutine adjust_cutoff

  subroutine verify_configuration()
    integer :: p, op, b, s1, s2
    integer, allocatable :: propagated(:)

    allocate(propagated(nsite))
    propagated = spin

    do p = 0, max_order - 1
      op = opstring(p)
      if (op == 0) cycle
      b = op/2
      s1 = bonds(1, b)
      s2 = bonds(2, b)
      if (propagated(s1) == propagated(s2)) &
        error stop '算符串作用在平行自旋上'
      if (mod(op, 2) == 1) then
        propagated(s1) = -propagated(s1)
        propagated(s2) = -propagated(s2)
      end if
    end do

    if (any(propagated /= spin)) error stop '算符串不满足虚时间周期边界'
    deallocate(propagated)
  end subroutine verify_configuration

  subroutine reset_accumulators(acc)
    real(dp), intent(out) :: acc(6)
    acc = 0.0_dp
  end subroutine reset_accumulators

  subroutine measure(acc)
    real(dp), intent(inout) :: acc(6)
    integer :: i, p, op, b, ms, mz, nsquare
    real(dp) :: rn

    ms = 0
    mz = 0
    do i = 1, nsite
      ms = ms + eta(i)*spin(i)
      mz = mz + spin(i)
    end do
    ms = ms/2
    mz = mz/2

    nsquare = 0
    do p = 0, max_order - 1
      op = opstring(p)
      if (op == 0) cycle
      b = op/2
      if (bond_type(b) == 1) nsquare = nsquare + 1
    end do

    rn = real(norder, dp)
    acc(1) = acc(1) + rn
    acc(2) = acc(2) + rn*rn
    acc(3) = acc(3) + real(ms, dp)**2
    acc(4) = acc(4) + abs(real(ms, dp))
    acc(5) = acc(5) + real(mz, dp)**2
    acc(6) = acc(6) + real(nsquare, dp)
  end subroutine measure

  subroutine finish_bin(acc, msteps, values)
    real(dp), intent(in) :: acc(6)
    integer, intent(in) :: msteps
    real(dp), intent(out) :: values(nobs)
    real(dp) :: norm, mean_n, mean_n2, mean_nsquare, mean_nlink

    norm = real(msteps, dp)
    mean_n = acc(1)/norm
    mean_n2 = acc(2)/norm
    mean_nsquare = acc(6)/norm
    mean_nlink = mean_n - mean_nsquare

    values(1) = -mean_n/(beta*real(nsite, dp)) + &
      coupling*real(nbond, dp)/(4.0_dp*real(nsite, dp))
    values(2) = (mean_n2 - mean_n*mean_n - mean_n)/real(nsite, dp)
    values(3) = 3.0_dp*(acc(3)/norm)/real(nsite, dp)**2
    values(4) = beta*(acc(5)/norm)/real(nsite, dp)
    values(5) = -mean_nsquare/(beta*real(nsite, dp)) + &
      coupling*real(nb_square, dp)/(4.0_dp*real(nsite, dp))
    values(6) = -mean_nlink/(beta*real(nsite, dp)) + &
      coupling*real(nb_link, dp)/(4.0_dp*real(nsite, dp))
    values(7) = (acc(4)/norm)/real(nsite, dp)
    values(8) = mean_n
  end subroutine finish_bin

  subroutine write_summary(bin_values, nbins)
    integer, intent(in) :: nbins
    real(dp), intent(in) :: bin_values(nobs, nbins)
    character(len=32), parameter :: names(nobs) = [character(len=32) :: &
      'energy_per_site', 'specific_heat_per_site', 'staggered_m2', &
      'uniform_susceptibility', 'square_bond_energy_per_site', &
      'link_bond_energy_per_site', 'abs_staggered_mz', &
      'mean_expansion_order']
    real(dp) :: means(nobs), errors(nobs)
    integer :: i, unit_id

    do i = 1, nobs
      means(i) = sum(bin_values(i, :))/real(nbins, dp)
      errors(i) = sqrt(sum((bin_values(i, :) - means(i))**2)/ &
        real(nbins*(nbins - 1), dp))
    end do

    open(newunit=unit_id, file='results_summary.txt', status='replace', &
         action='write')
    write(unit_id, '(a)') 'Square-octagon S=1/2 Heisenberg SSE'
    write(unit_id, '(a,f8.3)') 'J_square = J_link = ', coupling
    write(unit_id, '(a,2(i0,1x),a,i0,a,f12.6)') 'Lx Ly = ', lx, ly, &
      'N = ', nsite, '  beta = ', beta
    write(unit_id, '(a)') 'observable                              mean                  stderr'
    do i = 1, nobs
      write(unit_id, '(a32,2(2x,es22.14))') names(i), means(i), errors(i)
    end do
    close(unit_id)

    write(*, '(a)') '主要结果：'
    do i = 1, 4
      write(*, '(a32,2x,es14.6,a,es12.4)') names(i), means(i), ' +/- ', errors(i)
    end do
  end subroutine write_summary

  subroutine deallocate_all()
    deallocate(spin, eta, bonds, bond_type, opstring)
    deallocate(first_op, last_op, vertex)
  end subroutine deallocate_all

end module sse_square_octagon_model

program sse_square_octagon
  use sse_square_octagon_model
  implicit none

  integer :: nbins, msteps, isteps, seed
  integer :: step, bin, unit_id
  character(len=256) :: input_file
  real(dp) :: accumulators(6)
  real(dp), allocatable :: bin_values(:,:)

  input_file = 'sse_square_octagon.in'
  if (command_argument_count() >= 1) call get_command_argument(1, input_file)
  call read_input(nbins, msteps, isteps, seed, trim(input_file))
  call init_rng(seed)
  call make_lattice()
  call validate_lattice()
  call init_configuration()

  write(*, '(a,i0,a)') '开始热化：', isteps, ' sweeps'
  do step = 1, isteps
    call diagonal_update()
    call loop_update()
    call adjust_cutoff()
  end do
  call verify_configuration()
  write(*, '(a,i0,a,i0)') '热化完成：展开阶数 n=', norder, ', 截断 M=', max_order

  allocate(bin_values(nobs, nbins))
  open(newunit=unit_id, file='results_bins.csv', status='replace', action='write')
  write(unit_id, '(a)') 'bin,energy_per_site,specific_heat_per_site,'// &
    'staggered_m2,uniform_susceptibility,square_bond_energy_per_site,'// &
    'link_bond_energy_per_site,abs_staggered_mz,mean_expansion_order'

  do bin = 1, nbins
    call reset_accumulators(accumulators)
    do step = 1, msteps
      call diagonal_update()
      call loop_update()
      call adjust_cutoff()
      call measure(accumulators)
    end do
    call verify_configuration()
    call finish_bin(accumulators, msteps, bin_values(:, bin))
    write(unit_id, '(i0,8(",",es24.16))') bin, bin_values(:, bin)
    flush(unit_id)
    write(*, '(a,i0,a,i0)') '完成 bin ', bin, '/', nbins
  end do
  close(unit_id)

  call write_summary(bin_values, nbins)
  call deallocate_all()
  deallocate(bin_values)
end program sse_square_octagon
