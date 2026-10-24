module sse_graph_model
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: nobs = 16

  integer :: lx, ly, nbasis, nedge_cell, nsite, nbond, ntypes, nplus, nminus
  integer :: norder, max_order
  real(dp) :: beta, coupling_sum
  integer, allocatable :: spin(:), eta(:)
  integer, allocatable :: bonds(:,:), bond_type(:), bond_shift(:,:)
  integer, allocatable :: opstring(:), first_op(:), last_op(:), vertex(:)
  real(dp), allocatable :: couplings(:), bond_coupling(:)

contains

  subroutine read_input_header(input_file, graph_file, requested_types)
    character(len=*), intent(in) :: input_file
    character(len=*), intent(out) :: graph_file
    integer, intent(out) :: requested_types
    integer :: unit_id, ios

    open(newunit=unit_id, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop '无法打开 SSE 输入文件'
    read(unit_id, '(a)', iostat=ios) graph_file
    if (ios /= 0) error stop '无法读取格点图路径'
    read(unit_id, *, iostat=ios) lx, ly
    if (ios /= 0) error stop '无法读取 Lx 和 Ly'
    read(unit_id, *, iostat=ios) beta, requested_types
    if (ios /= 0) error stop '无法读取 beta 和键类型数'
    close(unit_id)

    if (lx < 2 .or. ly < 2) error stop 'Lx 和 Ly 必须至少为 2'
    if (beta <= 0.0_dp) error stop 'beta 必须为正数'
    if (requested_types < 1) error stop '键类型数必须为正整数'
  end subroutine read_input_header


  subroutine read_sampling_input(input_file, requested_types, nbins, equil_sweeps, bin_sweeps, seed)
    character(len=*), intent(in) :: input_file
    integer, intent(in) :: requested_types
    integer, intent(out) :: nbins, equil_sweeps, bin_sweeps, seed
    integer :: unit_id, ios, dummy_lx, dummy_ly, dummy_types
    real(dp) :: dummy_beta
    character(len=512) :: dummy_path

    allocate(couplings(requested_types))
    open(newunit=unit_id, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop '无法重新打开 SSE 输入文件'
    read(unit_id, '(a)') dummy_path
    read(unit_id, *) dummy_lx, dummy_ly
    read(unit_id, *) dummy_beta, dummy_types
    read(unit_id, *, iostat=ios) couplings
    if (ios /= 0) error stop '无法读取各类交换常数'
    read(unit_id, *, iostat=ios) nbins, equil_sweeps, bin_sweeps, seed
    if (ios /= 0) error stop '无法读取采样参数'
    close(unit_id)

    if (dummy_types /= requested_types) error stop '输入文件键类型数不一致'
    if (any(couplings <= 0.0_dp)) error stop '无符号问题 SSE 要求所有 J 为正数'
    if (nbins < 2 .or. equil_sweeps < 1 .or. bin_sweeps < 1) &
      error stop '采样参数不合法'
    if (seed <= 0) error stop '随机种子必须为正整数'
  end subroutine read_sampling_input


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


  integer function site_index(x, y, basis)
    integer, intent(in) :: x, y, basis
    site_index = basis + nbasis*(modulo(x, lx) + modulo(y, ly)*lx)
  end function site_index


  subroutine read_graph(graph_file, requested_types)
    character(len=*), intent(in) :: graph_file
    integer, intent(in) :: requested_types
    integer :: unit_id, ios, dimension, edge, cell_x, cell_y, bond
    integer :: i, j, dx, dy, dz, kind, site1, site2
    integer :: phase(3)
    integer, allocatable :: basis1(:), basis2(:), shifts(:,:), types(:)
    integer, allocatable :: eta_basis(:), degree(:)

    open(newunit=unit_id, file=trim(graph_file), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop '无法打开格点图文件'
    read(unit_id, *, iostat=ios) nbasis, nedge_cell, dimension
    if (ios /= 0 .or. nbasis < 1 .or. nedge_cell < 1) error stop '格点图文件头错误'
    if (dimension /= 3) error stop '当前仅支持三分量平移格式'

    allocate(basis1(nedge_cell), basis2(nedge_cell), shifts(3, nedge_cell))
    allocate(types(nedge_cell), eta_basis(nbasis))
    do edge = 1, nedge_cell
      read(unit_id, *, iostat=ios) i, j, dx, dy, dz, kind
      if (ios /= 0) error stop '格点图连边格式错误'
      if (i < 1 .or. i > nbasis .or. j < 1 .or. j > nbasis) &
        error stop '格点图含越界编号'
      if (dz /= 0) error stop '二维 SSE 程序不接受 z 方向连边'
      basis1(edge) = i
      basis2(edge) = j
      shifts(:, edge) = [dx, dy, dz]
      types(edge) = kind
    end do
    read(unit_id, *, iostat=ios) eta_basis
    if (ios /= 0) error stop '无法读取二分子格符号'
    read(unit_id, *, iostat=ios) phase
    if (ios /= 0) error stop '无法读取子格平移相位'
    close(unit_id)

    ntypes = maxval(types)
    if (ntypes /= requested_types .or. minval(types) < 1) &
      error stop '输入交换常数与格点图键类型不匹配'
    if (any(abs(eta_basis) /= 1)) error stop '反铁磁 SSE 需要二分图 eta=±1'
    if (any(phase < 0) .or. any(phase > 1)) error stop '子格平移相位必须为 0 或 1'

    nsite = nbasis*lx*ly
    nbond = nedge_cell*lx*ly
    allocate(bonds(2, nbond), bond_type(nbond), bond_shift(2, nbond))
    allocate(bond_coupling(nbond), eta(nsite), spin(nsite), degree(nsite))

    bond = 0
    do cell_y = 0, ly - 1
      do cell_x = 0, lx - 1
        do i = 1, nbasis
          site1 = site_index(cell_x, cell_y, i)
          eta(site1) = eta_basis(i)*(-1)**(phase(1)*cell_x + phase(2)*cell_y)
        end do
        do edge = 1, nedge_cell
          bond = bond + 1
          bonds(1, bond) = site_index(cell_x, cell_y, basis1(edge))
          bonds(2, bond) = site_index(cell_x + shifts(1, edge), &
            cell_y + shifts(2, edge), basis2(edge))
          bond_type(bond) = types(edge)
          bond_shift(:, bond) = shifts(1:2, edge)
          bond_coupling(bond) = couplings(types(edge))
        end do
      end do
    end do
    if (bond /= nbond) error stop '超胞键数错误'
    nplus = count(eta == 1)
    nminus = count(eta == -1)
    if (nplus + nminus /= nsite .or. min(nplus, nminus) < 1) &
      error stop '二分子格计数错误'

    degree = 0
    do bond = 1, nbond
      site1 = bonds(1, bond)
      site2 = bonds(2, bond)
      if (site1 == site2) error stop '有限超胞产生自连接，请增大尺寸'
      if (eta(site1)*eta(site2) /= -1) &
        error stop '有限超胞不是二分图，请使用兼容尺寸'
      degree(site1) = degree(site1) + 1
      degree(site2) = degree(site2) + 1
    end do
    coupling_sum = sum(bond_coupling)

    write(*, '(a,i0,a,i0,a,i0,a,i0)') '格点检查通过：N=', nsite, ', Nb=', nbond, &
      ', z_min=', minval(degree), ', z_max=', maxval(degree)

    deallocate(basis1, basis2, shifts, types, eta_basis, degree)
  end subroutine read_graph


  subroutine init_configuration()
    allocate(first_op(nsite), last_op(nsite))
    spin = eta
    max_order = max(64, nsite)
    allocate(opstring(0:max_order - 1), vertex(0:4*max_order - 1))
    opstring = 0
    norder = 0
  end subroutine init_configuration


  subroutine diagonal_update()
    integer :: p, op, bond, site1, site2
    real(dp) :: probability, weight

    do p = 0, max_order - 1
      op = opstring(p)
      if (op == 0) then
        if (norder == max_order) cycle
        bond = 1 + int(rand01()*real(nbond, dp))
        site1 = bonds(1, bond)
        site2 = bonds(2, bond)
        if (spin(site1) /= spin(site2)) then
          weight = 0.5_dp*beta*real(nbond, dp)*bond_coupling(bond)
          probability = weight/real(max_order - norder, dp)
          if (rand01() <= min(1.0_dp, probability)) then
            opstring(p) = 2*bond
            norder = norder + 1
          end if
        end if
      else if (mod(op, 2) == 0) then
        bond = op/2
        weight = 0.5_dp*beta*real(nbond, dp)*bond_coupling(bond)
        probability = real(max_order - norder + 1, dp)/weight
        if (rand01() <= min(1.0_dp, probability)) then
          opstring(p) = 0
          norder = norder - 1
        end if
      else
        bond = op/2
        site1 = bonds(1, bond)
        site2 = bonds(2, bond)
        spin(site1) = -spin(site1)
        spin(site2) = -spin(site2)
      end if
    end do
  end subroutine diagonal_update


  subroutine loop_update()
    integer :: site, op, bond, site1, site2
    integer :: vertex0, vertex1, vertex2

    first_op = -1
    last_op = -1
    do vertex0 = 0, 4*max_order - 1, 4
      op = opstring(vertex0/4)
      if (op == 0) then
        vertex(vertex0:vertex0 + 3) = -1
        cycle
      end if
      bond = op/2
      site1 = bonds(1, bond)
      site2 = bonds(2, bond)
      vertex1 = last_op(site1)
      vertex2 = last_op(site2)
      if (vertex1 == -1) then
        first_op(site1) = vertex0
      else
        vertex(vertex1) = vertex0
        vertex(vertex0) = vertex1
      end if
      if (vertex2 == -1) then
        first_op(site2) = vertex0 + 1
      else
        vertex(vertex2) = vertex0 + 1
        vertex(vertex0 + 1) = vertex2
      end if
      last_op(site1) = vertex0 + 2
      last_op(site2) = vertex0 + 3
    end do

    do site = 1, nsite
      if (first_op(site) /= -1) then
        vertex(last_op(site)) = first_op(site)
        vertex(first_op(site)) = last_op(site)
      end if
    end do

    do vertex0 = 0, 4*max_order - 1, 2
      if (vertex(vertex0) < 0) cycle
      vertex1 = vertex0
      if (rand01() < 0.5_dp) then
        do
          opstring(vertex1/4) = ieor(opstring(vertex1/4), 1)
          vertex(vertex1) = -2
          vertex2 = ieor(vertex1, 1)
          vertex1 = vertex(vertex2)
          vertex(vertex2) = -2
          if (vertex1 == vertex0) exit
        end do
      else
        do
          vertex(vertex1) = -1
          vertex2 = ieor(vertex1, 1)
          vertex1 = vertex(vertex2)
          vertex(vertex2) = -1
          if (vertex1 == vertex0) exit
        end do
      end if
    end do

    do site = 1, nsite
      if (first_op(site) == -1) then
        if (rand01() < 0.5_dp) spin(site) = -spin(site)
      else if (vertex(first_op(site)) == -2) then
        spin(site) = -spin(site)
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
    allocate(opstring(0:new_order - 1), vertex(0:4*new_order - 1))
    opstring = 0
    opstring(0:max_order - 1) = old_string
    deallocate(old_string)
    max_order = new_order
  end subroutine adjust_cutoff


  subroutine verify_configuration()
    integer :: p, op, bond, site1, site2
    integer, allocatable :: propagated(:)

    allocate(propagated(nsite))
    propagated = spin
    do p = 0, max_order - 1
      op = opstring(p)
      if (op == 0) cycle
      bond = op/2
      site1 = bonds(1, bond)
      site2 = bonds(2, bond)
      if (propagated(site1) == propagated(site2)) &
        error stop '算符串作用在平行自旋上'
      if (mod(op, 2) == 1) then
        propagated(site1) = -propagated(site1)
        propagated(site2) = -propagated(site2)
      end if
    end do
    if (any(propagated /= spin)) error stop '算符串不满足虚时间周期边界'
    deallocate(propagated)
  end subroutine verify_configuration


  subroutine reset_accumulators(acc)
    real(dp), intent(out) :: acc(11)
    acc = 0.0_dp
  end subroutine reset_accumulators


  subroutine measure(acc)
    real(dp), intent(inout) :: acc(11)
    integer :: p, op, bond, site, site1, site2
    integer :: ms, mz, ma, mb, current_x, current_y
    integer, allocatable :: propagated(:)
    real(dp) :: rn

    ma = 0
    mb = 0
    do site = 1, nsite
      if (eta(site) == 1) then
        ma = ma + spin(site)
      else
        mb = mb + spin(site)
      end if
    end do
    ma = ma/2
    mb = mb/2
    ms = ma - mb
    mz = ma + mb
    current_x = 0
    current_y = 0
    allocate(propagated(nsite))
    propagated = spin
    do p = 0, max_order - 1
      op = opstring(p)
      if (op == 0 .or. mod(op, 2) == 0) cycle
      bond = op/2
      site1 = bonds(1, bond)
      site2 = bonds(2, bond)
      current_x = current_x + propagated(site1)*bond_shift(1, bond)
      current_y = current_y + propagated(site1)*bond_shift(2, bond)
      propagated(site1) = -propagated(site1)
      propagated(site2) = -propagated(site2)
    end do
    if (any(propagated /= spin)) error stop '绕数测量传播失败'
    deallocate(propagated)

    rn = real(norder, dp)
    acc(1) = acc(1) + rn
    acc(2) = acc(2) + rn*rn
    acc(3) = acc(3) + real(ms, dp)**2
    acc(4) = acc(4) + real(ms, dp)**4
    acc(5) = acc(5) + real(mz, dp)**2
    acc(6) = acc(6) + real(current_x, dp)**2
    acc(7) = acc(7) + real(current_y, dp)**2
    acc(8) = acc(8) + abs(real(ms, dp))
    acc(9) = acc(9) + real(ma, dp)**2
    acc(10) = acc(10) + real(mb, dp)**2
    acc(11) = acc(11) + real(ma, dp)*real(mb, dp)
  end subroutine measure


  subroutine finish_bin(acc, bin_sweeps, values)
    real(dp), intent(in) :: acc(11)
    integer, intent(in) :: bin_sweeps
    real(dp), intent(out) :: values(nobs)
    real(dp) :: norm, mean_n, mean_n2, mean_ms2, mean_ms4, mean_mz2
    real(dp) :: mean_ma2, mean_mb2, mean_mamb
    real(dp) :: current_x2, current_y2, winding2, stiffness_x, stiffness_y, stiffness

    norm = real(bin_sweeps, dp)
    mean_n = acc(1)/norm
    mean_n2 = acc(2)/norm
    mean_ms2 = acc(3)/norm
    mean_ms4 = acc(4)/norm
    mean_mz2 = acc(5)/norm
    mean_ma2 = acc(9)/norm
    mean_mb2 = acc(10)/norm
    mean_mamb = acc(11)/norm
    current_x2 = acc(6)/norm
    current_y2 = acc(7)/norm
    winding2 = current_x2/real(lx, dp)**2 + current_y2/real(ly, dp)**2
    stiffness_x = 1.5_dp*current_x2/(beta*real(lx, dp)**2)
    stiffness_y = 1.5_dp*current_y2/(beta*real(ly, dp)**2)
    stiffness = 0.5_dp*(stiffness_x + stiffness_y)

    values(1) = -mean_n/(beta*real(nsite, dp)) + coupling_sum/(4.0_dp*real(nsite, dp))
    values(2) = (mean_n2 - mean_n*mean_n - mean_n)/real(nsite, dp)
    values(3) = 3.0_dp*mean_ms2/real(nsite, dp)**2
    values(4) = mean_ms4/real(nsite, dp)**4
    if (mean_ms2 > 0.0_dp) then
      values(5) = 1.0_dp - mean_ms4/(3.0_dp*mean_ms2**2)
    else
      values(5) = 0.0_dp
    end if
    values(6) = 3.0_dp*mean_mz2/real(nsite, dp)**2
    values(7) = beta*mean_mz2/real(nsite, dp)
    values(8) = 3.0_dp*mean_ma2/real(nplus, dp)**2
    values(9) = 3.0_dp*mean_mb2/real(nminus, dp)**2
    values(10) = 3.0_dp*mean_mamb/(real(nplus, dp)*real(nminus, dp))
    values(11) = stiffness_x
    values(12) = stiffness_y
    values(13) = stiffness
    values(14) = stiffness*sqrt(real(lx*ly, dp))
    values(15) = winding2
    values(16) = mean_n
  end subroutine finish_bin


  subroutine write_outputs(bin_values, nbins)
    integer, intent(in) :: nbins
    real(dp), intent(in) :: bin_values(nobs, nbins)
    character(len=32), parameter :: names(nobs) = [character(len=32) :: &
      'energy_per_site', 'specific_heat_per_site', 'staggered_m2', &
      'staggered_m4', 'binder_z_staggered', 'ferrimagnetic_m2', &
      'uniform_susceptibility', 'sublattice_plus_m2', 'sublattice_minus_m2', &
      'sublattice_cross', 'spin_stiffness_cell_x', 'spin_stiffness_cell_y', &
      'spin_stiffness_cell_average', 'spin_stiffness_times_L', &
      'winding_squared', 'mean_expansion_order']
    real(dp) :: means(nobs), errors(nobs)
    integer :: observable, bin, unit_id

    do observable = 1, nobs
      means(observable) = sum(bin_values(observable, :))/real(nbins, dp)
      errors(observable) = sqrt(sum((bin_values(observable, :) - means(observable))**2)/ &
        real(nbins*(nbins - 1), dp))
    end do

    open(newunit=unit_id, file='sse_graph_bins.csv', status='replace', action='write')
    write(unit_id, '(a)') 'bin,energy_per_site,specific_heat_per_site,staggered_m2,'// &
      'staggered_m4,binder_z_staggered,ferrimagnetic_m2,uniform_susceptibility,'// &
      'sublattice_plus_m2,sublattice_minus_m2,sublattice_cross,'// &
      'spin_stiffness_cell_x,spin_stiffness_cell_y,spin_stiffness_cell_average,'// &
      'spin_stiffness_times_L,winding_squared,mean_expansion_order'
    do bin = 1, nbins
      write(unit_id, '(i0,16(",",es24.16))') bin, bin_values(:, bin)
    end do
    close(unit_id)

    open(newunit=unit_id, file='sse_graph_summary.txt', status='replace', action='write')
    write(unit_id, '(a)') 'S=1/2 antiferromagnetic Heisenberg SSE on a periodic graph'
    write(unit_id, '(a,2(i0,1x),a,i0,a,es18.10)') 'Lx Ly = ', lx, ly, 'N = ', nsite, ' beta = ', beta
    write(unit_id, '(a,2(i0,1x))') 'Nplus Nminus = ', nplus, nminus
    write(unit_id, '(a,*(es18.10,1x))') 'J by bond type = ', couplings
    write(unit_id, '(a)') 'observable                              mean                  stderr'
    do observable = 1, nobs
      write(unit_id, '(a32,2(2x,es22.14))') names(observable), means(observable), errors(observable)
    end do
    close(unit_id)

    write(*, '(a)') '主要结果：'
    do observable = 1, 14
      write(*, '(a32,2x,es14.6,a,es12.4)') names(observable), means(observable), &
        ' +/- ', errors(observable)
    end do
  end subroutine write_outputs


  subroutine deallocate_all()
    deallocate(spin, eta, bonds, bond_type, bond_shift, bond_coupling)
    deallocate(opstring, first_op, last_op, vertex, couplings)
  end subroutine deallocate_all

end module sse_graph_model


program sse_graph
  use sse_graph_model
  implicit none

  character(len=512) :: input_file, graph_file
  integer :: requested_types, nbins, equil_sweeps, bin_sweeps, seed
  integer :: sweep, bin
  real(dp) :: accumulators(11)
  real(dp), allocatable :: bin_values(:,:)

  input_file = 'sse_graph.in'
  if (command_argument_count() >= 1) call get_command_argument(1, input_file)
  call read_input_header(trim(input_file), graph_file, requested_types)
  call read_sampling_input(trim(input_file), requested_types, nbins, equil_sweeps, bin_sweeps, seed)
  call init_rng(seed)
  call read_graph(trim(graph_file), requested_types)
  call init_configuration()

  write(*, '(a,i0,a)') '开始热化：', equil_sweeps, ' sweeps'
  do sweep = 1, equil_sweeps
    call diagonal_update()
    call loop_update()
    call adjust_cutoff()
  end do
  call verify_configuration()
  write(*, '(a,i0,a,i0)') '热化完成：展开阶数 n=', norder, ', 截断 M=', max_order

  allocate(bin_values(nobs, nbins))
  do bin = 1, nbins
    call reset_accumulators(accumulators)
    do sweep = 1, bin_sweeps
      call diagonal_update()
      call loop_update()
      call adjust_cutoff()
      call measure(accumulators)
    end do
    call verify_configuration()
    call finish_bin(accumulators, bin_sweeps, bin_values(:, bin))
    write(*, '(a,i0,a,i0)') '完成 bin ', bin, '/', nbins
  end do

  call write_outputs(bin_values, nbins)
  call deallocate_all()
  deallocate(bin_values)
end program sse_graph
