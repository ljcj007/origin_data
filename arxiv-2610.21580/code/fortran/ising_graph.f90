module ising_graph_model
  use iso_fortran_env, only: int64, real64
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: nobs = 10

  integer :: lx, ly, nbasis, nedge_cell, nsite, nbond
  integer :: max_degree
  real(dp) :: temperature, field, coupling
  character(len=16) :: update_method
  integer, allocatable :: spin(:), eta(:)
  integer, allocatable :: edge1(:), edge2(:), degree(:), neighbors(:,:)
  integer, allocatable :: cluster_stack(:)
  logical, allocatable :: in_cluster(:)

contains

  subroutine read_input(input_file, graph_file, nbins, equil_sweeps, bin_sweeps, seed)
    character(len=*), intent(in) :: input_file
    character(len=*), intent(out) :: graph_file
    integer, intent(out) :: nbins, equil_sweeps, bin_sweeps, seed
    integer :: unit_id, ios

    open(newunit=unit_id, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop '无法打开 Ising 输入文件'
    read(unit_id, '(a)', iostat=ios) graph_file
    if (ios /= 0) error stop '无法读取格点图路径'
    read(unit_id, *, iostat=ios) lx, ly
    if (ios /= 0) error stop '无法读取 Lx 和 Ly'
    read(unit_id, *, iostat=ios) temperature, field, coupling
    if (ios /= 0) error stop '无法读取 T、h 和 J'
    read(unit_id, *, iostat=ios) nbins, equil_sweeps, bin_sweeps, seed
    if (ios /= 0) error stop '无法读取采样参数'
    update_method = 'metropolis'
    read(unit_id, '(a)', iostat=ios) update_method
    if (ios /= 0) update_method = 'metropolis'
    update_method = adjustl(update_method)
    close(unit_id)

    if (lx < 2 .or. ly < 2) error stop 'Lx 和 Ly 必须至少为 2'
    if (temperature <= 0.0_dp) error stop '温度必须为正数'
    if (coupling <= 0.0_dp) error stop '本程序约定反铁磁 J 必须为正数'
    if (nbins < 2 .or. equil_sweeps < 1 .or. bin_sweeps < 1) &
      error stop '采样参数不合法'
    if (seed <= 0) error stop '随机种子必须为正整数'
    if (trim(update_method) /= 'metropolis' .and. trim(update_method) /= 'wolff') &
      error stop '更新算法必须为 metropolis 或 wolff'
    if (trim(update_method) == 'wolff' .and. abs(field) > tiny(field)) &
      error stop 'Wolff 更新当前仅支持零外场'
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


  integer function site_index(x, y, basis)
    integer, intent(in) :: x, y, basis
    site_index = basis + nbasis*(modulo(x, lx) + modulo(y, ly)*lx)
  end function site_index


  subroutine read_graph(graph_file)
    character(len=*), intent(in) :: graph_file
    integer :: unit_id, ios, dimension, edge, cell_x, cell_y, bond
    integer :: i, j, dx, dy, dz, bond_type, site1, site2
    integer :: phase(3)
    integer, allocatable :: basis1(:), basis2(:), shifts(:,:), types(:), eta_basis(:)

    open(newunit=unit_id, file=trim(graph_file), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop '无法打开格点图文件'
    read(unit_id, *, iostat=ios) nbasis, nedge_cell, dimension
    if (ios /= 0 .or. nbasis < 1 .or. nedge_cell < 1) error stop '格点图文件头错误'
    if (dimension /= 3) error stop '当前仅支持三分量平移格式'

    allocate(basis1(nedge_cell), basis2(nedge_cell), shifts(3, nedge_cell))
    allocate(types(nedge_cell), eta_basis(nbasis))
    do edge = 1, nedge_cell
      read(unit_id, *, iostat=ios) i, j, dx, dy, dz, bond_type
      if (ios /= 0) error stop '格点图连边格式错误'
      if (i < 1 .or. i > nbasis .or. j < 1 .or. j > nbasis) &
        error stop '格点图含越界编号'
      if (dz /= 0) error stop '二维 Ising 程序不接受 z 方向连边'
      basis1(edge) = i
      basis2(edge) = j
      shifts(:, edge) = [dx, dy, dz]
      types(edge) = bond_type
    end do
    read(unit_id, *, iostat=ios) eta_basis
    if (ios /= 0) error stop '无法读取二分子格符号'
    read(unit_id, *, iostat=ios) phase
    if (ios /= 0) error stop '无法读取子格平移相位'
    close(unit_id)

    if (any(abs(eta_basis) /= 1)) error stop '反铁磁模拟需要二分图 eta=±1'
    if (any(phase < 0) .or. any(phase > 1)) error stop '子格平移相位必须为 0 或 1'
    if (any(types < 1)) error stop '键类型必须为正整数'

    nsite = nbasis*lx*ly
    nbond = nedge_cell*lx*ly
    allocate(edge1(nbond), edge2(nbond), eta(nsite), spin(nsite))
    allocate(cluster_stack(nsite), in_cluster(nsite))
    in_cluster = .false.

    bond = 0
    do cell_y = 0, ly - 1
      do cell_x = 0, lx - 1
        do i = 1, nbasis
          site1 = site_index(cell_x, cell_y, i)
          eta(site1) = eta_basis(i)*(-1)**(phase(1)*cell_x + phase(2)*cell_y)
        end do
        do edge = 1, nedge_cell
          bond = bond + 1
          edge1(bond) = site_index(cell_x, cell_y, basis1(edge))
          edge2(bond) = site_index(cell_x + shifts(1, edge), &
            cell_y + shifts(2, edge), basis2(edge))
        end do
      end do
    end do
    if (bond /= nbond) error stop '超胞键数错误'

    do bond = 1, nbond
      if (edge1(bond) == edge2(bond)) error stop '有限超胞产生自连接，请增大尺寸'
      if (eta(edge1(bond))*eta(edge2(bond)) /= -1) &
        error stop '有限超胞不是二分图，请使用兼容尺寸'
    end do

    allocate(degree(nsite))
    degree = 0
    do bond = 1, nbond
      degree(edge1(bond)) = degree(edge1(bond)) + 1
      degree(edge2(bond)) = degree(edge2(bond)) + 1
    end do
    max_degree = maxval(degree)
    allocate(neighbors(max_degree, nsite))
    neighbors = 0
    degree = 0
    do bond = 1, nbond
      site1 = edge1(bond)
      site2 = edge2(bond)
      degree(site1) = degree(site1) + 1
      degree(site2) = degree(site2) + 1
      neighbors(degree(site1), site1) = site2
      neighbors(degree(site2), site2) = site1
    end do

    write(*, '(a,i0,a,i0,a,i0,a,i0)') '格点检查通过：N=', nsite, ', Nb=', nbond, &
      ', z_min=', minval(degree), ', z_max=', maxval(degree)

    deallocate(basis1, basis2, shifts, types, eta_basis)
  end subroutine read_graph


  subroutine init_spins()
    integer :: site
    real(dp) :: random_value

    do site = 1, nsite
      call random_number(random_value)
      spin(site) = merge(1, -1, random_value < 0.5_dp)
    end do
  end subroutine init_spins


  subroutine metropolis_sweep(flipped)
    integer(int64), intent(inout) :: flipped
    integer :: attempt, site, neighbor_index, neighbor_sum
    real(dp) :: random_value, delta_energy

    do attempt = 1, nsite
      call random_number(random_value)
      site = 1 + int(random_value*real(nsite, dp))
      neighbor_sum = 0
      do neighbor_index = 1, degree(site)
        neighbor_sum = neighbor_sum + spin(neighbors(neighbor_index, site))
      end do
      delta_energy = -2.0_dp*real(spin(site), dp)* &
        (coupling*real(neighbor_sum, dp) - field)
      call random_number(random_value)
      if (delta_energy <= 0.0_dp .or. random_value < exp(-delta_energy/temperature)) then
        spin(site) = -spin(site)
        flipped = flipped + 1_int64
      end if
    end do
  end subroutine metropolis_sweep


  subroutine wolff_sweep(flipped, clusters)
    integer(int64), intent(inout) :: flipped, clusters
    integer :: seed_site, target_tau, head, cluster_size, site, neighbor_index, neighbor
    real(dp) :: random_value, add_probability

    add_probability = 1.0_dp - exp(-2.0_dp*coupling/temperature)
    call random_number(random_value)
    seed_site = 1 + int(random_value*real(nsite, dp))
    target_tau = eta(seed_site)*spin(seed_site)
    head = 1
    cluster_size = 1
    cluster_stack(1) = seed_site
    in_cluster(seed_site) = .true.

    do while (head <= cluster_size)
      site = cluster_stack(head)
      head = head + 1
      do neighbor_index = 1, degree(site)
        neighbor = neighbors(neighbor_index, site)
        if (in_cluster(neighbor)) cycle
        if (eta(neighbor)*spin(neighbor) /= target_tau) cycle
        call random_number(random_value)
        if (random_value < add_probability) then
          cluster_size = cluster_size + 1
          cluster_stack(cluster_size) = neighbor
          in_cluster(neighbor) = .true.
        end if
      end do
    end do

    do head = 1, cluster_size
      site = cluster_stack(head)
      spin(site) = -spin(site)
      in_cluster(site) = .false.
    end do
    flipped = flipped + int(cluster_size, int64)
    clusters = clusters + 1_int64
  end subroutine wolff_sweep


  subroutine mc_sweep(flipped, clusters)
    integer(int64), intent(inout) :: flipped, clusters

    if (trim(update_method) == 'wolff') then
      call wolff_sweep(flipped, clusters)
    else
      call metropolis_sweep(flipped)
    end if
  end subroutine mc_sweep


  subroutine reset_accumulators(acc)
    real(dp), intent(out) :: acc(7)
    acc = 0.0_dp
  end subroutine reset_accumulators


  subroutine measure(acc)
    real(dp), intent(inout) :: acc(7)
    integer :: bond
    real(dp) :: energy, magnetization, staggered

    energy = -field*real(sum(spin), dp)
    do bond = 1, nbond
      energy = energy + coupling*real(spin(edge1(bond))*spin(edge2(bond)), dp)
    end do
    magnetization = real(sum(spin), dp)
    staggered = real(sum(eta*spin), dp)

    acc(1) = acc(1) + energy
    acc(2) = acc(2) + energy*energy
    acc(3) = acc(3) + magnetization
    acc(4) = acc(4) + magnetization*magnetization
    acc(5) = acc(5) + abs(staggered)
    acc(6) = acc(6) + staggered*staggered
    acc(7) = acc(7) + staggered**4
  end subroutine measure


  subroutine finish_bin(acc, bin_sweeps, flipped, clusters, values)
    real(dp), intent(in) :: acc(7)
    integer, intent(in) :: bin_sweeps
    integer(int64), intent(in) :: flipped, clusters
    real(dp), intent(out) :: values(nobs)
    real(dp) :: norm, mean_energy, mean_energy2, mean_m, mean_m2
    real(dp) :: mean_abs_ms, mean_ms2, mean_ms4

    norm = real(bin_sweeps, dp)
    mean_energy = acc(1)/norm
    mean_energy2 = acc(2)/norm
    mean_m = acc(3)/norm
    mean_m2 = acc(4)/norm
    mean_abs_ms = acc(5)/norm
    mean_ms2 = acc(6)/norm
    mean_ms4 = acc(7)/norm

    values(1) = mean_energy/real(nsite, dp)
    values(2) = (mean_energy2 - mean_energy**2)/(real(nsite, dp)*temperature**2)
    values(3) = mean_abs_ms/real(nsite, dp)
    values(4) = mean_ms2/real(nsite, dp)**2
    values(5) = mean_ms4/real(nsite, dp)**4
    if (mean_ms2 > 0.0_dp) then
      values(6) = 1.0_dp - mean_ms4/(3.0_dp*mean_ms2**2)
    else
      values(6) = 0.0_dp
    end if
    values(7) = (mean_m2 - mean_m**2)/(real(nsite, dp)*temperature)
    values(8) = (mean_ms2 - mean_abs_ms**2)/(real(nsite, dp)*temperature)
    values(9) = real(flipped, dp)/(real(bin_sweeps, dp)*real(nsite, dp))
    values(10) = real(clusters, dp)/real(bin_sweeps, dp)
  end subroutine finish_bin


  subroutine write_outputs(bin_values, nbins)
    integer, intent(in) :: nbins
    real(dp), intent(in) :: bin_values(nobs, nbins)
    character(len=40), parameter :: names(nobs) = [character(len=40) :: &
      'energy_per_site', 'specific_heat_per_site', 'abs_staggered_m', &
      'staggered_m2', 'staggered_m4', 'binder_staggered', 'uniform_susceptibility', &
      'connected_staggered_susceptibility', 'spin_flip_fraction', 'clusters_per_sweep']
    real(dp) :: means(nobs), errors(nobs)
    integer :: observable, bin, unit_id

    do observable = 1, nobs
      means(observable) = sum(bin_values(observable, :))/real(nbins, dp)
      errors(observable) = sqrt(sum((bin_values(observable, :) - means(observable))**2)/ &
        real(nbins*(nbins - 1), dp))
    end do

    open(newunit=unit_id, file='ising_bins.csv', status='replace', action='write')
    write(unit_id, '(a)') 'bin,energy_per_site,specific_heat_per_site,abs_staggered_m,'// &
      'staggered_m2,staggered_m4,binder_staggered,uniform_susceptibility,'// &
      'connected_staggered_susceptibility,spin_flip_fraction,clusters_per_sweep'
    do bin = 1, nbins
      write(unit_id, '(i0,10(",",es24.16))') bin, bin_values(:, bin)
    end do
    close(unit_id)

    open(newunit=unit_id, file='ising_summary.txt', status='replace', action='write')
    write(unit_id, '(a)') 'Classical antiferromagnetic Ising model on a periodic graph'
    write(unit_id, '(a,2(i0,1x),a,i0)') 'Lx Ly = ', lx, ly, 'N = ', nsite
    write(unit_id, '(a,3(es18.10,1x))') 'T h J = ', temperature, field, coupling
    write(unit_id, '(a,a)') 'update method = ', trim(update_method)
    write(unit_id, '(a)') 'observable                              mean                  stderr'
    do observable = 1, nobs
      write(unit_id, '(a32,2(2x,es22.14))') names(observable), means(observable), errors(observable)
    end do
    close(unit_id)

    write(*, '(a)') '主要结果：'
    do observable = 1, 7
      write(*, '(a32,2x,es14.6,a,es12.4)') names(observable), means(observable), &
        ' +/- ', errors(observable)
    end do
  end subroutine write_outputs


  subroutine deallocate_all()
    deallocate(spin, eta, edge1, edge2, degree, neighbors, cluster_stack, in_cluster)
  end subroutine deallocate_all

end module ising_graph_model


program ising_graph
  use iso_fortran_env, only: int64
  use ising_graph_model
  implicit none

  character(len=512) :: input_file, graph_file
  integer :: nbins, equil_sweeps, bin_sweeps, seed
  integer :: sweep, bin
  integer(int64) :: flipped, clusters
  real(dp) :: accumulators(7)
  real(dp), allocatable :: bin_values(:,:)

  input_file = 'ising_graph.in'
  if (command_argument_count() >= 1) call get_command_argument(1, input_file)
  call read_input(trim(input_file), graph_file, nbins, equil_sweeps, bin_sweeps, seed)
  call init_rng(seed)
  call read_graph(trim(graph_file))
  call init_spins()

  flipped = 0_int64
  clusters = 0_int64
  write(*, '(a,i0,a)') '开始热化：', equil_sweeps, ' sweeps'
  do sweep = 1, equil_sweeps
    call mc_sweep(flipped, clusters)
  end do

  allocate(bin_values(nobs, nbins))
  do bin = 1, nbins
    call reset_accumulators(accumulators)
    flipped = 0_int64
    clusters = 0_int64
    do sweep = 1, bin_sweeps
      call mc_sweep(flipped, clusters)
      call measure(accumulators)
    end do
    call finish_bin(accumulators, bin_sweeps, flipped, clusters, bin_values(:, bin))
    write(*, '(a,i0,a,i0)') '完成 bin ', bin, '/', nbins
  end do

  call write_outputs(bin_values, nbins)
  call deallocate_all()
  deallocate(bin_values)
end program ising_graph
