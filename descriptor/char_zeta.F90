module char_zeta

implicit none

integer :: num_atomic_types = 5

integer, dimension(5) :: atomic_numbers = (/3, 7, 15, 16, 32/)
double precision, dimension(5) :: &
atomic_radii = (/1.67d0, 0.56d0, 0.98d0, 0.88d0, 1.25d0/)
double precision, dimension(5) :: &
ionic_radii = (/1.45d0, 0.65d0, 1.00d0, 1.00d0, 1.25d0/)
double precision, dimension(5) :: &
covalent_radii = (/1.34d0, 0.75d0, 1.06d0, 1.02d0, 1.22d0/)
double precision, dimension(5) :: &
vdw_radii = (/1.82d0, 1.55d0, 1.80d0, 1.80d0, 1.86d0/)
double precision, dimension(5) :: &
crystal_radii = (/0.90d0, 0.30d0, 0.31d0, 0.43d0, 0.53d0/)

contains

function get_radii(r_type, atomic_num)

  implicit none

  integer, intent(in) :: atomic_num
  integer, intent(in) :: r_type
  double precision :: get_radii

  integer :: local_index

  do local_index = 1, num_atomic_types
    if (atomic_numbers(local_index) .eq. atomic_num) then
      if (r_type .eq. 1) then
        get_radii = atomic_radii(local_index)
      elseif (r_type .eq. 2) then
        get_radii = ionic_radii(local_index)
      elseif (r_type .eq. 3) then
        get_radii = covalent_radii(local_index)
      elseif (r_type .eq. 4) then
        get_radii = vdw_radii(local_index)
      elseif (r_type .eq. 5) then
        get_radii = crystal_radii(local_index)
      end if
    end if
  end do

end function get_radii

subroutine scalar(& 
  r_scalar, &
  el_i, el_j, el_k, el_l, &
  interaction, &
  r_type, &
  factor)
  
  implicit none
  
  double precision :: r_scalar
  integer :: interaction
  integer, intent(in) :: r_type
  integer :: el_i, el_j; integer, optional :: el_k, el_l
  double precision :: r_i, r_j, r_k, r_l
  double precision :: denominator
  double precision, dimension(6), intent(out) :: factor

  r_i = get_radii(r_type, el_i)
  r_j = get_radii(r_type, el_j)
  if (present(el_k)) then
    r_k = get_radii(r_type, el_k)
  end if
  if (present(el_l)) then
    r_l = get_radii(r_type, el_l)
  end if

  denominator = 0.0d0

  if (interaction .eq. 2) then
    factor(1) = r_scalar * 0.5d0 / (r_i ** 2 + r_j ** 2)
  elseif (interaction .eq. 3) then
    denominator = denominator + &
      (r_i * r_j) ** 2 + (r_j * r_k) ** 2 + (r_i * r_k) ** 2
    factor(1) = r_scalar * 0.5d0 * r_k ** 2 / denominator ! R_ij
    factor(2) = r_scalar * 0.5d0 * r_j ** 2 / denominator ! R_ik
    factor(3) = r_scalar * 0.5d0 * r_i ** 2 / denominator ! R_jk
  elseif (interaction .eq. 4) then
    denominator = denominator + &
      (r_j * r_k * r_l) ** 2 + (r_i * r_k * r_l) ** 2 + &
      (r_i * r_j * r_l) **2 + (r_i * r_j * r_k) ** 2
    factor(1) = r_scalar * 0.5d0 * (r_k * r_l) ** 2 / denominator ! R_ij
    factor(2) = r_scalar * 0.5d0 * (r_j * r_l) ** 2 / denominator ! R_ik
    factor(3) = r_scalar * 0.5d0 * (r_j * r_k) ** 2 / denominator ! R_il
    factor(4) = r_scalar * 0.5d0 * (r_i * r_l) ** 2 / denominator ! R_jk
    factor(5) = r_scalar * 0.5d0 * (r_i * r_k) ** 2 / denominator ! R_jl
    factor(6) = r_scalar * 0.5d0 * (r_i * r_j) ** 2 / denominator ! R_kl
  end if

end subroutine scalar

function neighbor_one(elin, elexp)

  implicit none
  
  integer, intent(in) :: elin, elexp
  integer :: neighbor_one
!f2py intent(in) :: elin, elexp
  
  neighbor_one = 0
  if (elin .eq. elexp) then
    neighbor_one = 1
  end if
  
end function neighbor_one

function neighbor_two(& 
  elin1, elin2, &
  elexp1, elexp2)

  implicit none
    
  integer, intent(in) :: elin1, elin2
  integer, intent(in) :: elexp1, elexp2
  integer :: nelin1, nelin2
  integer :: nelexp1, nelexp2
  integer :: neighbor_two
!f2py intent(in) :: elin1, elin2, elexp1, elexp2  

  neighbor_two = 0
  nelin1 = min(elin1, elin2)
  nelin2 = max(elin1, elin2)
  
  nelexp1 = min(elexp1, elexp2)
  nelexp2 = max(elexp1, elexp2)
    
  if ((nelin1 .eq. nelexp1) .and. & 
        (nelin2 .eq. nelexp2)) then
    neighbor_two = 1
  end if
    
end function neighbor_two     

function neighbor_three(& 
  elin1, elin2, elin3, &
  elexp1, elexp2, elexp3)

  implicit none
    
  integer, intent(in) :: elin1, elin2, elin3
  integer, intent(in) :: elexp1, elexp2, elexp3
  integer :: nelin1, nelin2, nelin3
  integer :: nelexp1, nelexp2, nelexp3
  integer :: neighbor_three 
!f2py intent(in) :: elin1, elin2, elin3
!f2py intent(in) :: elexp1, elexp2, elexp3
  
  
  neighbor_three = 0
  nelin1 = min(elin1, elin2, elin3)
  nelin3 = max(elin1, elin2, elin3)
  nelin2 = elin1 + elin2 + elin3 - nelin1 - nelin3
  
  nelexp1 = min(elexp1, elexp2, elexp3)
  nelexp3 = max(elexp1, elexp2, elexp3)
  nelexp2 = elexp1 + elexp2 + elexp3 - nelexp1 - nelexp3
    
  if ((nelin1 .eq. nelexp1) .and. & 
        (nelin2 .eq. nelexp2) .and. & 
        (nelin3 .eq. nelexp3)) then
    neighbor_three = 1
  end if
    
end function neighbor_three

function cutoff_fxn(r, rc, cutofffn_code, p_gamma)

  implicit none
  
  double precision :: r, rc, pi, cutoff_fxn
  integer :: p_gamma
  integer :: cutofffn_code

!f2py intent(in) :: r
!f2py intent(in) :: rc
!f2py intent(in) :: cutofffn_code
!f2py intent(in) :: p_gamma
    
  cutoff_fxn = 0.0d0
        
  if (r < rc) then
    if (cutofffn_code == 1) then
      pi = 4.0d0 * datan(1.0d0)
      cutoff_fxn = 0.5d0 * (cos(pi*r/rc) + 1.0d0)
    elseif (cutofffn_code == 2) then
      cutoff_fxn = 1. + p_gamma &
        * (r / rc) ** (p_gamma + 1) &
        - (p_gamma + 1) * (r / rc) ** p_gamma
    end if
  end if
  
end function cutoff_fxn

function cutoff_fxn_prime(r, rc, cutofffn_code, p_gamma)

  implicit none
  
  double precision :: r, rc, cutoff_fxn_prime, pi
  integer :: p_gamma
  integer :: cutofffn_code
  
!f2py intent(in) :: r
!f2py intent(in) :: rc
!f2py intent(in) :: cutofffn_code
!f2py intent(in) :: p_gamma
    
  cutoff_fxn_prime = 0.0d0
    
  if (r < rc) then
    if (cutofffn_code == 1) then
      pi = 4.0d0 * datan(1.0d0)
      cutoff_fxn_prime = -0.5d0 * pi * sin(pi*r/rc) / rc
    elseif (cutofffn_code == 2) then
      cutoff_fxn_prime = (p_gamma * (p_gamma + 1) / rc) &
        * ((r / rc) ** p_gamma - (r / rc) ** (p_gamma - 1))
    end if
  end if
  
end function cutoff_fxn_prime 

function Kronecker(i, j)

  implicit none
  
  integer :: i, j
  integer :: Kronecker
  
  if (i == j) then
    Kronecker = 1
  else
    Kronecker = 0
  end if
  
end function Kronecker

function dRij_dRmd(i, j, &
  Ri, Rj, &
  m, d)
  
  implicit none
  
  integer :: xyz
  integer :: i, j, m, d
  double precision, dimension(3) :: Ri, Rj, Rij_vector
  double precision :: dRij_dRmd, Rij
    
  ! m is the index of the center atom, d is the force direction
  do xyz = 1, 3
    Rij_vector(xyz) = Rj(xyz) - Ri(xyz)
  end do
  Rij = sqrt(dot_product(Rij_vector, Rij_vector))
  if ((m == i) .AND. (i /= j)) then
    ! indexing starts from 0 in python, but from 1 in fortran
    ! so we need d + 1 in the following expression
    ! R_ij = (\sum_{d = 1}^3(R_{id} - R_{jd})**2)^{frac{1}{2}}
    dRij_dRmd = - (Rj(d + 1) - Ri(d + 1)) / Rij
  else if ((m == j) .AND. (i /= j)) then
    dRij_dRmd = (Rj(d + 1) - Ri(d + 1)) / Rij
  else
    dRij_dRmd = 0.0d0
  end if
  
end function dRij_dRmd

function dRij_dRml_vector(i, j, m, l)

  implicit none

  integer:: i, j, m, l, c1
  integer, dimension(3):: dRij_dRml_vector

  if ((m /= i) .AND. (m /= j)) then
    dRij_dRml_vector(1) = 0
    dRij_dRml_vector(2) = 0
    dRij_dRml_vector(3) = 0
  else
    c1 = Kronecker(m, j) - Kronecker(m, i)
    dRij_dRml_vector(1) = c1 * Kronecker(0, l)
    dRij_dRml_vector(2) = c1 * Kronecker(1, l)
    dRij_dRml_vector(3) = c1 * Kronecker(2, l)
  end if

end function dRij_dRml_vector

function dCos_ijk_dR_ml(i, j, k, &
  Ri, Rj, Rk, &
  m, l)

  implicit none

  integer:: xyz
  integer:: i, j, k, m, l
  double precision:: dCos_ijk_dR_ml
  double precision, dimension(3):: Ri, Rj, Rk
  integer, dimension(3):: dRijdRml, dRikdRml
  double precision:: scalar_dRijdRml, scalar_dRikdRml
  double precision, dimension(3):: Rij_vector, Rik_vector
  double precision:: Rij, Rik

  do xyz = 1, 3
    Rij_vector(xyz) = Rj(xyz) - Ri(xyz)
    Rik_vector(xyz) = Rk(xyz) - Ri(xyz)
  end do
  Rij = sqrt(dot_product(Rij_vector, Rij_vector))
  Rik = sqrt(dot_product(Rik_vector, Rik_vector))

  ! Costheta = dot_product(Rij_vector, Rik_vector) / (Rij * Rik)
  ! the dot_product part
  dRijdRml = dRij_dRml_vector(i, j, m, l)
  dRikdRml = dRij_dRml_vector(i, k, m, l)
  dCos_ijk_dR_ml = 0.0d0
  if ((dRijdRml(1) /= 0.0d0) .OR. (dRijdRml(2) /= 0.0d0) .OR. (dRijdRml(3) /= 0.0d0)) then
    dCos_ijk_dR_ml = dCos_ijk_dR_ml + 1.0d0 / (Rij * Rik) * dot_product(dRijdRml, Rik_vector)
  end if
  if ((dRikdRml(1) /= 0.0d0) .OR. (dRikdRml(2) /= 0.0d0) .OR. (dRikdRml(3) /= 0.0d0)) then
    dCos_ijk_dR_ml = dCos_ijk_dR_ml + 1.0d0 / (Rij * Rik) * dot_product(dRikdRml, Rij_vector)
  end if
  ! the fraction part
  scalar_dRijdRml = dRij_dRmd(i, j, Ri, Rj, m, l)
  scalar_dRikdRml = dRij_dRmd(i, k, Ri, Rk, m, l)
  if (scalar_dRijdRml /= 0.0d0) then
    dCos_ijk_dR_ml = dCos_ijk_dR_ml - (dot_product(Rij_vector, Rik_vector) * &
      1.0d0 / Rik) * 1.0d0 / (Rij * Rij) * scalar_dRijdRml
  end if 
  if (scalar_dRikdRml /= 0.0d0) then
    dCos_ijk_dR_ml = dCos_ijk_dR_ml - (dot_product(Rij_vector, Rik_vector) * &
      1.0d0 / Rij) * 1.0d0 / (Rik * Rik) * scalar_dRikdRml
  end if

end function dCos_ijk_dR_ml
  
function calculate_g2(& 
  num_rscalars, r_scalars, &
  num_n, &
  atomic_num, ri, &
  atomic_num_n, pos_n, &
  atomicn_first, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
    
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer :: num_n
  integer :: atomic_num
  double precision, dimension(3) :: ri
  integer, dimension(num_n) :: atomic_num_n
  double precision, dimension(num_n, 3) :: pos_n
  integer :: atomicn_first
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision ::  rc; integer :: cutofffn_code, p_gamma
  
  integer :: j, xyz, r_count, s_count, t_count
  double precision, dimension(3) :: Rij_vector
  double precision :: Rij
  
  double precision, dimension(6) :: factor
  integer :: match
  
  double precision, dimension(num_rscalars * num_rtypes) :: ridge
    
!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n
!f2py intent(in) :: atomic_num, ri
!f2py intent(in) :: atomic_num_n, pos_n
!f2py intent(in) :: atomicn_first
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, cutofffn_code, p_gamma
!f2py intent(out) :: ridge

  ridge = 0.0d0

  do j = 1, num_n

    if (atomicn_first .ne. -1) then
      match = neighbor_one(atomic_num_n(j), atomicn_first)
    else
      match = 1
    end if

    if (match .eq. 1) then
      do xyz = 1, 3
        Rij_vector(xyz) = pos_n(j, xyz) - ri(xyz)
      end do
      Rij = sqrt(dot_product(Rij_vector, Rij_vector))

      t_count = 0
      do r_count = 1, num_rtypes
        do s_count = 1, num_rscalars
        
          t_count = t_count + 1
          call scalar(&
            r_scalar = r_scalars(s_count), &
            el_i = atomic_num, el_j = atomic_num_n(j), &
            interaction = 2, &
            r_type = r_types(r_count), &
            factor = factor)

          ridge(t_count) = ridge(t_count) + &
            exp(- factor(1) * (Rij  / rc ) ** 2.0d0) * &
            cutoff_fxn(Rij, rc, cutofffn_code, p_gamma)
            
        end do 
      end do
    end if 
    
  end do
  
end function calculate_g2

function calculate_g3(& 
  num_rscalars, r_scalars, &
  num_n, &
  atomic_num, ri, &
  atomic_num_n, pos_n, &
  atomicn_first, atomicn_second, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
    
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer :: num_n
  integer :: atomic_num
  double precision, dimension(3) :: ri
  integer, dimension(num_n) :: atomic_num_n
  double precision, dimension(num_n, 3) :: pos_n
  integer :: atomicn_first, atomicn_second
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision ::  rc; integer :: cutofffn_code, p_gamma
    
  integer :: j, k, xyz, r_count, s_count, t_count
  double precision, dimension(3) :: Rij_vector, Rik_vector, Rjk_vector
  double precision :: Rij, Rik, Rjk
  
  double precision, dimension(6) :: factor
  
  integer :: match 
    
  double precision, dimension(num_rscalars * num_rtypes) :: ridge
    
!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n
!f2py intent(in) :: atomic_num, ri
!f2py intent(in) :: atomic_num_n, pos_n
!f2py intent(in) :: atomicn_first, atomicn_second
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, cutofffn_code, p_gamma
!f2py intent(out) :: ridge

  ridge = 0.0d0

  do j = 1, num_n
    do k = (j + 1), num_n

      if (atomicn_first .ne. -1) then
        match = neighbor_two(&
          atomic_num_n(j), atomic_num_n(k), &
          atomicn_first, atomicn_second)
      else
        match = 1
      end if
        
      if (match .eq. 1) then 
        do xyz = 1, 3
          Rij_vector(xyz) = pos_n(j, xyz) - ri(xyz)
          Rik_vector(xyz) = pos_n(k, xyz) - ri(xyz)
          Rjk_vector(xyz) = pos_n(k, xyz) - pos_n(j, xyz)
        end do
        Rij = sqrt(dot_product(Rij_vector, Rij_vector))
        Rik = sqrt(dot_product(Rik_vector, Rik_vector))
        Rjk = sqrt(dot_product(Rjk_vector, Rjk_vector))

        t_count = 0
        do r_count = 1, num_rtypes
          do s_count = 1, num_rscalars
            
            t_count = t_count + 1
            call scalar(&
              r_scalar = r_scalars(s_count), &
              el_i = atomic_num, el_j = atomic_num_n(j), el_k = atomic_num_n(k), &
              interaction = 3, &
              r_type = r_types(r_count), &
              factor = factor)

            ridge(t_count) = ridge(t_count) + &
              exp(- factor(1) * (Rij /  rc)  ** 2.0d0 &
                - factor(2) * (Rik /  rc)  ** 2.0d0 &
                - factor(3) * (Rjk /  rc)  ** 2.0d0) * &
              cutoff_fxn(Rij, rc, cutofffn_code, p_gamma) * &
              cutoff_fxn(Rik, rc, cutofffn_code, p_gamma)
                
          end do 
        end do
      end if 
      
    end do
  end do
  
end function calculate_g3

function calculate_g3_angle(& 
  num_rscalars, r_scalars, &
  num_n, &
  atomic_num, ri, &
  atomic_num_n, pos_n, &
  atomicn_first, atomicn_second, &
  g_gamma, g_zeta, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
  
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer:: num_n
  integer:: atomic_num
  double precision, dimension(3):: ri
  integer, dimension(num_n):: atomic_num_n
  double precision, dimension(num_n, 3):: pos_n
  integer :: atomicn_first, atomicn_second
  double precision:: g_gamma, g_zeta
  
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision:: rc; integer:: cutofffn_code, p_gamma
    
  integer:: j, k, xyz, r_count, s_count, t_count
  double precision, dimension(3):: Rij_vector, Rik_vector, Rjk_vector
  double precision:: Rij, Rik, Rjk, costheta
  
  double precision, dimension(6):: factor
      
  integer::  match
  
  double precision, dimension(num_rscalars * num_rtypes) :: ridge

!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n
!f2py intent(in) :: atomic_num, ri
!f2py intent(in) :: atomic_num_n, pos_n
!f2py intent(in) :: atomicn_first, atomicn_second
!f2py intent(in) :: g_gamma, g_zeta
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, cutofffn_code, p_gamma
!f2py intent(out) :: ridge

  ridge = 0.0d0
  
  do j = 1, num_n
    do k = (j + 1), num_n ! for a pair of atoms 
    
      if (atomicn_first .ne. -1) then
        match = neighbor_two(&
          atomic_num_n(j), atomic_num_n(k), &
          atomicn_first, atomicn_second)
      else
        match = 1
      end if
        
      if (match .eq. 1) then
        do xyz = 1, 3
          Rij_vector(xyz) = pos_n(j, xyz) - ri(xyz)
          Rik_vector(xyz) = pos_n(k, xyz) - ri(xyz)
          Rjk_vector(xyz) = pos_n(k, xyz) - pos_n(j, xyz)
        end do
        Rij = sqrt(dot_product(Rij_vector, Rij_vector))
        Rik = sqrt(dot_product(Rik_vector, Rik_vector))
        Rjk = sqrt(dot_product(Rjk_vector, Rjk_vector))
        costheta = dot_product(Rij_vector, Rik_vector) / (Rij * Rik)
        
        t_count = 0
        do r_count = 1, num_rtypes
          do s_count = 1, num_rscalars
          
            t_count = t_count + 1
            call scalar(&
              r_scalar = r_scalars(s_count), &
              el_i = atomic_num, el_j = atomic_num_n(j), el_k = atomic_num_n(k), &
              interaction = 3, &
              r_type = r_types(r_count), &
              factor = factor)
        
            ridge(t_count) = ridge(t_count) + &
              2.0d0 ** (1.0d0 - g_zeta) * (1.0d0 + g_gamma * costheta) ** g_zeta * &
              exp(- factor(1) * (Rij /  rc)  ** 2.0d0 &
                - factor(2) * (Rik /  rc)  ** 2.0d0 &
                - factor(3) * (Rjk /  rc)  ** 2.0d0) * &
              cutoff_fxn(Rij, rc, cutofffn_code, p_gamma) * &
              cutoff_fxn(Rik, rc, cutofffn_code, p_gamma) * &
              cutoff_fxn(Rjk, rc, cutofffn_code, p_gamma)
              
          end do 
        end do
        
      end if
    end do
  end do

end function calculate_g3_angle

function calculate_g3_angle_partial(& 
  num_rscalars, r_scalars, &
  num_n, &
  atomic_num, ri, &
  atomic_num_n, pos_n, &
  atomicn_first, atomicn_second, &
  g_gamma, g_zeta, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
  
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer:: num_n
  integer:: atomic_num
  double precision, dimension(3):: ri
  integer, dimension(num_n):: atomic_num_n
  double precision, dimension(num_n, 3):: pos_n
  integer :: atomicn_first, atomicn_second
  double precision:: g_gamma, g_zeta
  
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision:: rc; integer:: cutofffn_code, p_gamma
    
  integer:: j, k, xyz, r_count, s_count, t_count
  double precision, dimension(3):: Rij_vector, Rik_vector
  double precision:: Rij, Rik, costheta
  
  double precision, dimension(6):: factor_ij
  double precision, dimension(6):: factor_ik
      
  integer::  match
  
  double precision, dimension(num_rscalars * num_rtypes) :: ridge

!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n
!f2py intent(in) :: atomic_num, ri
!f2py intent(in) :: atomic_num_n, pos_n
!f2py intent(in) :: atomicn_first, atomicn_second
!f2py intent(in) :: g_gamma, g_zeta
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, cutofffn_code, p_gamma
!f2py intent(out) :: ridge

  ridge = 0.0d0
  
  do j = 1, num_n
    do k = (j + 1), num_n ! for a pair of atoms 
    
      if (atomicn_first .ne. -1) then
        match = neighbor_two(&
          atomic_num_n(j), atomic_num_n(k), &
          atomicn_first, atomicn_second)
      else
        match = 1
      end if
        
      if (match .eq. 1) then
        do xyz = 1, 3
          Rij_vector(xyz) = pos_n(j, xyz) - ri(xyz)
          Rik_vector(xyz) = pos_n(k, xyz) - ri(xyz)
        end do
        Rij = sqrt(dot_product(Rij_vector, Rij_vector))
        Rik = sqrt(dot_product(Rik_vector, Rik_vector))
        costheta = dot_product(Rij_vector, Rik_vector) / (Rij * Rik)
        
        t_count = 0
        do r_count = 1, num_rtypes
          do s_count = 1, num_rscalars
          
            t_count = t_count + 1
            call scalar(&
              r_scalar = r_scalars(s_count), & 
              el_i = atomic_num, el_j = atomic_num_n(j), &
              interaction = 2, &
              r_type = r_types(r_count), &
              factor = factor_ij)
            call scalar(&
              r_scalar = r_scalars(s_count), &
              el_i = atomic_num, el_j = atomic_num_n(k), &
              interaction = 2, &
              r_type = r_types(r_count), &
              factor = factor_ik)
        
            ridge(t_count) = ridge(t_count) + &
              2.0d0 ** (1.0d0 - g_zeta) * (1.0d0 + g_gamma * costheta) ** g_zeta * &
              exp(- factor_ij(1) * (Rij /  rc)  ** 2.0d0 &
                - factor_ik(1) * (Rik /  rc)  ** 2.0d0) * &
              cutoff_fxn(Rij, rc, cutofffn_code, p_gamma) * &
              cutoff_fxn(Rik, rc, cutofffn_code, p_gamma)
              
          end do 
        end do
        
      end if
    end do
  end do

end function calculate_g3_angle_partial

function calculate_g4(& 
  num_rscalars, r_scalars, &
  num_n, &
  atomic_num, ri, &
  atomic_num_n, pos_n, &
  atomicn_first, atomicn_second, atomicn_third, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
    
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer :: num_n
  integer :: atomic_num
  double precision, dimension(3) :: ri
  integer, dimension(num_n) :: atomic_num_n
  double precision, dimension(num_n, 3) :: pos_n
  integer :: atomicn_first, atomicn_second, atomicn_third
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision ::  rc; integer :: cutofffn_code, p_gamma
    
  integer :: j, k, l, xyz, r_count, s_count, t_count
  double precision, dimension(3) :: Rij_vector, Rik_vector, Ril_vector
  double precision, dimension(3) :: Rjk_vector, Rjl_vector
  double precision, dimension(3) :: Rkl_vector
  
  double precision :: Rij, Rik, Ril
  double precision :: Rjk, Rjl
  double precision :: Rkl
  
  double precision, dimension(6) :: factor

  integer :: match
    
  double precision, dimension(num_rscalars * num_rtypes) :: ridge

!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n
!f2py intent(in) :: atomic_num, ri
!f2py intent(in) :: atomic_num_n, pos_n
!f2py intent(in) :: atomicn_first, atomicn_second, atomicn_third
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, cutofffn_code, p_gamma
!f2py intent(out) :: ridge

  ridge = 0.0d0

  do j = 1, num_n
    do k = (j + 1), num_n
      do l = (k + 1), num_n

        if (atomicn_first .ne. -1) then
          match = neighbor_three(&
            atomic_num_n(j), atomic_num_n(k), atomic_num_n(l), &
            atomicn_first, atomicn_second, atomicn_third)
        else
          match = 1
        end if
          
        if (match .eq. 1) then 
          do xyz = 1, 3
            Rij_vector(xyz) = pos_n(j, xyz) - ri(xyz)
            Rik_vector(xyz) = pos_n(k, xyz) - ri(xyz)
            Ril_vector(xyz) = pos_n(l, xyz) - ri(xyz)
            Rjk_vector(xyz) = pos_n(k, xyz) - pos_n(j, xyz)
            Rjl_vector(xyz) = pos_n(l, xyz) - pos_n(j, xyz)
            Rkl_vector(xyz) = pos_n(l, xyz) - pos_n(k, xyz)
          end do
          Rij = sqrt(dot_product(Rij_vector, Rij_vector))
          Rik = sqrt(dot_product(Rik_vector, Rik_vector))
          Ril = sqrt(dot_product(Ril_vector, Ril_vector))
          Rjk = sqrt(dot_product(Rjk_vector, Rjk_vector))
          Rjl = sqrt(dot_product(Rjl_vector, Rjl_vector))
          Rkl = sqrt(dot_product(Rkl_vector, Rkl_vector))

          t_count = 0
          do r_count = 1, num_rtypes
            do s_count = 1, num_rscalars
            
              t_count = t_count + 1
              call scalar(&
                r_scalar = r_scalars(s_count), &
                el_i = atomic_num, &
                el_j = atomic_num_n(j), &
                el_k = atomic_num_n(k), &
                el_l = atomic_num_n(l), &
                interaction = 4, &
                r_type = r_types(r_count), &
                factor = factor)

              ridge(t_count) = ridge(t_count) + &
                exp(- factor(1) * (Rij / rc) ** 2.0d0 &
                  - factor(2) * (Rik / rc) ** 2.0d0 &
                  - factor(3) * (Ril / rc) ** 2.0d0 &
                  - factor(4) * (Rjk / rc) ** 2.0d0 &
                  - factor(5) * (Rjl / rc) ** 2.0d0 &
                  - factor(6) * (Rkl / rc) ** 2.0d0) * &
                cutoff_fxn(Rij, rc, cutofffn_code, p_gamma) * &
                cutoff_fxn(Rik, rc, cutofffn_code, p_gamma) * &
                cutoff_fxn(Ril, rc, cutofffn_code, p_gamma)
              
            end do 
          end do
        end if 
        
      end do 
    end do
  end do

end function calculate_g4

function calculate_g2_prime(& 
  num_rscalars, r_scalars, &
  num_n_n, &
  ind_n_n, atomic_num_n_n, pos_n_n, &
  ind_n, atomic_num_n, pos_n, &
  ind_center, d, & 
  atomicn_first, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
  
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer :: num_n_n
  integer, dimension(num_n_n) :: ind_n_n
  integer, dimension(num_n_n) :: atomic_num_n_n
  double precision, dimension(num_n_n, 3) :: pos_n_n
  integer :: ind_n
  integer :: atomic_num_n
  double precision, dimension(3) :: pos_n
  integer :: ind_center, d
  integer :: atomicn_first
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision ::  rc; integer :: cutofffn_code, p_gamma
    
  integer :: j, xyz, r_count, s_count, t_count
  double precision, dimension(3) :: Rj
  double precision, dimension(3) :: Rij_vector
  double precision :: Rij, dRijdRmd
  
  double precision, dimension(6) :: factor

  integer :: match 
  
  double precision :: derivates
  double precision, dimension(num_rscalars * num_rtypes) :: ridge
    
!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n_n
!f2py intent(in) :: ind_n_n, atomic_num_n_n, pos_n_n
!f2py intent(in) :: ind_n, atomic_num_n, pos_n
!f2py intent(in) :: ind_center, d
!f2py intent(in) :: atomicn_first
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, p_gamma, cutofffn_code
!f2py intent(out) :: ridge

  ridge = 0.0d0

  do j = 1, num_n_n

    if (atomicn_first .ne. -1) then
      match = neighbor_one(atomic_num_n_n(j), atomicn_first)
    else
      match = 1
    end if

    if (match .eq. 1) then 
      do xyz = 1, 3
        Rj(xyz) = pos_n_n(j, xyz)
        Rij_vector(xyz) = Rj(xyz) - pos_n(xyz)
      end do
      dRijdRmd = dRij_dRmd(ind_n, ind_n_n(j), pos_n, Rj, ind_center, d)
      Rij = sqrt(dot_product(Rij_vector, Rij_vector))

      t_count = 0
      do r_count = 1, num_rtypes
        do s_count = 1, num_rscalars
        
          t_count = t_count + 1
          call scalar(&
            r_scalar = r_scalars(s_count), &
            el_i = atomic_num_n, el_j = atomic_num_n_n(j), &
            interaction = 2, &
            r_type = r_types(r_count), &
            factor = factor)

          derivates = 0.0d0
          if (dRijdRmd /= 0.0d0) then
            derivates = derivates + &
              cutoff_fxn_prime(Rij, rc, cutofffn_code, p_gamma) * dRijdRmd - &
              cutoff_fxn(Rij, rc, cutofffn_code, p_gamma) * &
              2.0d0 * factor(1) / (rc ** 2.0d0) * Rij * dRijdRmd
          end if
          
          ridge(t_count) = ridge(t_count) + derivates * &
            exp(- factor(1)  * (Rij / rc) ** 2.0d0)
        
        end do
      end do
      
    end if 
  end do
  
end function calculate_g2_prime

function calculate_g3_prime(& 
  num_rscalars, r_scalars, &
  num_n_n, &
  ind_n_n, atomic_num_n_n, pos_n_n, &
  ind_n, atomic_num_n, pos_n, &
  ind_center, d, &
  atomicn_first, atomicn_second, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none

  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer :: num_n_n
  integer, dimension(num_n_n) :: ind_n_n
  integer, dimension(num_n_n) :: atomic_num_n_n
  double precision, dimension(num_n_n, 3) :: pos_n_n
  integer :: ind_n
  integer :: atomic_num_n
  double precision, dimension(3) :: pos_n
  integer :: ind_center, d
  integer :: atomicn_first, atomicn_second
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision ::  rc; integer :: cutofffn_code, p_gamma
    
  integer :: j, k, xyz, r_count, s_count, t_count
  double precision, dimension(3) :: Rj, Rk 
  double precision, dimension(3) :: Rij_vector, Rik_vector, Rjk_vector
  double precision :: Rij, Rik, Rjk
  
  double precision :: fcRij, fcRik, fcRijfcRik
  double precision :: dRijdRmd, dRikdRmd, dRjkdRmd
  
  double precision, dimension(6) :: factor

  integer :: match
    
  double precision :: derivates
  double precision, dimension(num_rscalars * num_rtypes) :: ridge

!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n_n
!f2py intent(in) :: ind_n_n, atomic_num_n_n, pos_n_n
!f2py intent(in) :: ind_n, atomic_num_n, pos_n
!f2py intent(in) :: ind_center, d
!f2py intent(in) :: atomicn_first, atomicn_second
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, p_gamma, cutofffn_code
!f2py intent(out) :: ridge
  
  ridge = 0.0d0

  do j = 1, num_n_n
    do k = (j + 1), num_n_n

      if (atomicn_first .ne. -1) then
        match = neighbor_two(&
          atomic_num_n_n(j), atomic_num_n_n(k), &
          atomicn_first, atomicn_second)
      else
        match = 1
      end if
        
      if (match .eq. 1) then 
        do xyz = 1, 3
          Rj(xyz) = pos_n_n(j, xyz) 
          Rk(xyz) = pos_n_n(k, xyz) 
          Rij_vector(xyz) = pos_n_n(j, xyz) - pos_n(xyz)
          Rik_vector(xyz) = pos_n_n(k, xyz) - pos_n(xyz)
          Rjk_vector(xyz) = pos_n_n(k, xyz) - pos_n_n(j, xyz)
        end do
        Rij = sqrt(dot_product(Rij_vector, Rij_vector))
        Rik = sqrt(dot_product(Rik_vector, Rik_vector))
        Rjk = sqrt(dot_product(Rjk_vector, Rjk_vector))
        dRijdRmd = dRij_dRmd(ind_n, ind_n_n(j), pos_n, Rj, ind_center, d)
        dRikdRmd = dRij_dRmd(ind_n, ind_n_n(k), pos_n, Rk, ind_center, d)
        dRjkdRmd = dRij_dRmd(ind_n_n(j), ind_n_n(k), Rj, Rk, ind_center, d)
        
        fcRij = cutoff_fxn(Rij, rc, cutofffn_code, p_gamma)
        fcRik = cutoff_fxn(Rik, rc, cutofffn_code, p_gamma)
        fcRijfcRik = fcRij * fcRik
        
        t_count = 0
        do r_count = 1, num_rtypes
          do s_count = 1, num_rscalars
          
            t_count = t_count + 1
            call scalar(&
              r_scalar = r_scalars(s_count), &
              el_i = atomic_num_n, el_j = atomic_num_n_n(j), el_k = atomic_num_n_n(k), &
              interaction = 3, &
              r_type = r_types(r_count), &
              factor = factor)
            ! fcRij * fcRik * exp(ij) * exp(ik) * exp(jk)
            derivates = 0.0d0
            if (dRijdRmd /= 0.0d0) then
              derivates = derivates + &
                cutoff_fxn_prime(Rij, rc, cutofffn_code, p_gamma) * dRijdRmd * fcRik - &
                fcRijfcRik * 2.0d0 * factor(1) / (rc ** 2.0d0) * Rij * dRijdRmd
            end if
            if (dRikdRmd /= 0.0d0) then
              derivates = derivates + &
                fcRij * cutoff_fxn_prime(Rik, rc, cutofffn_code, p_gamma) * dRikdRmd - &
                fcRijfcRik * 2.0d0 * factor(2) / (rc ** 2.0d0) * Rik * dRikdRmd
            end if
            if (dRjkdRmd /= 0.0d0) then
              derivates = derivates - &
                fcRijfcRik * 2.0d0 * factor(3) / (rc ** 2.0d0) * Rjk * dRjkdRmd
            end if
          
            ridge(t_count) = ridge(t_count) + derivates * &
              exp(- factor(1) * (Rij / rc) ** 2.0d0 &
                - factor(2) * (Rik / rc) ** 2.0d0 &
                - factor(3) * (Rjk / rc) ** 2.0d0)
              
          end do 
        end do
        
      end if 
    end do
  end do

end function calculate_g3_prime

function calculate_g3_angle_prime(&
  num_rscalars, r_scalars, &
  num_n_n, &
  ind_n_n, atomic_num_n_n, pos_n_n, &
  ind_n, atomic_num_n, pos_n, &
  ind_center, d, &
  atomicn_first, atomicn_second, &
  g_gamma, g_zeta, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
  
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer:: num_n_n
  integer, dimension(num_n_n):: ind_n_n
  integer, dimension(num_n_n):: atomic_num_n_n
  double precision, dimension(num_n_n, 3):: pos_n_n
  integer:: ind_n
  integer:: atomic_num_n
  double precision, dimension(3):: pos_n
  integer:: ind_center, d
  integer :: atomicn_first, atomicn_second
  double precision:: g_gamma, g_zeta
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision:: rc; integer:: cutofffn_code, p_gamma
    
  integer:: j, k, xyz, r_count, s_count, t_count
  double precision, dimension(3):: Rj, Rk
  double precision, dimension(3):: Rij_vector, Rik_vector, Rjk_vector
  double precision:: Rij, Rik, Rjk, costheta

  double precision:: fcRij, fcRik, fcRjk, fcRijfcRikfcRjk
  double precision:: dCosthetadRml
  double precision:: dRijdRml, dRikdRml, dRjkdRml

  double precision, dimension(6):: factor

  integer:: match

  double precision:: angular, radial
  double precision, dimension(num_rscalars * num_rtypes) :: ridge

!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n_n
!f2py intent(in) :: ind_n_n, atomic_num_n_n, pos_n_n
!f2py intent(in) :: ind_n, atomic_num_n, pos_n
!f2py intent(in) :: ind_center, d
!f2py intent(in) :: atomicn_first, atomicn_second
!f2py intent(in) :: g_gamma, g_zeta
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, p_gamma, cutofffn_code
!f2py intent(out) :: ridge
  
  ridge = 0.0d0
  
  do j = 1, num_n_n
    do k = (j + 1), num_n_n
    
      if (atomicn_first .ne. -1) then
        match = neighbor_two(&
          atomic_num_n_n(j), atomic_num_n_n(k), &
          atomicn_first, atomicn_second)
      else
        match = 1
      end if
        
      if (match .eq. 1) then
        do xyz = 1, 3
          Rj(xyz) = pos_n_n(j, xyz)
          Rk(xyz) = pos_n_n(k, xyz)
          Rij_vector(xyz) = Rj(xyz) - pos_n(xyz)
          Rik_vector(xyz) = Rk(xyz) - pos_n(xyz)
          Rjk_vector(xyz) = Rk(xyz) - Rj(xyz)
        end do
        Rij = sqrt(dot_product(Rij_vector, Rij_vector))
        Rik = sqrt(dot_product(Rik_vector, Rik_vector))
        Rjk = sqrt(dot_product(Rjk_vector, Rjk_vector))
        costheta = dot_product(Rij_vector, Rik_vector) / Rij / Rik
        dRijdRml = dRij_dRmd(ind_n, ind_n_n(j), pos_n, Rj, ind_center, d)
        dRikdRml = dRij_dRmd(ind_n, ind_n_n(k), pos_n, Rk, ind_center, d)
        dRjkdRml = dRij_dRmd(ind_n_n(j), ind_n_n(k), Rj, Rk, ind_center, d)
        dCosthetadRml = dCos_ijk_dR_ml(ind_n, &
          ind_n_n(j), ind_n_n(k), &
          pos_n, Rj, Rk, &
          ind_center, d)
          
        fcRij = cutoff_fxn(Rij, rc, cutofffn_code, p_gamma)
        fcRik = cutoff_fxn(Rik, rc, cutofffn_code, p_gamma)
        fcRjk = cutoff_fxn(Rjk, rc, cutofffn_code, p_gamma)
        fcRijfcRikfcRjk = fcRij * fcRik * fcRjk
        
        t_count = 0
        do r_count = 1, num_rtypes
          do s_count = 1, num_rscalars
          
            t_count = t_count + 1
            call scalar(&
              r_scalar = r_scalars(s_count), &
              el_i = atomic_num_n, el_j = atomic_num_n_n(j), el_k = atomic_num_n_n(k), &
              interaction = 3, &
              r_type = r_types(r_count), &
              factor = factor)

            ! 2 ** (1 - zeta) * (1 + gamma cos) ** zeta * & 
            ! fcRij * fcRik * fcRjk * exp(ij) * exp(ik) * exp(jk)
            angular = 0.0d0; radial = 0.0d0
            if (g_zeta == 1.0d0) then
              angular = g_zeta * g_gamma * dCosthetadRml
            else
              angular = g_zeta * g_gamma * dCosthetadRml * &
                (1.0d0 + g_gamma * costheta) ** (g_zeta - 1.0d0)
            end if
            angular = angular * fcRijfcRikfcRjk

            if (dRijdRml /= 0.0d0) then
              radial = radial + &
                cutoff_fxn_prime(Rij, rc, cutofffn_code, p_gamma) * dRijdRml * fcRik * fcRjk - &
                fcRijfcRikfcRjk * 2.0d0 * factor(1) / (rc ** 2.0d0) * Rij * dRijdRml
            end if
            if (dRikdRml /= 0.0d0) then
              radial = radial + &
                fcRij * cutoff_fxn_prime(Rik, rc, cutofffn_code, p_gamma) * dRikdRml * fcRjk - &
                fcRijfcRikfcRjk * 2.0d0 * factor(2) / (rc ** 2.0d0) * Rik * dRikdRml
            end if
            if (dRjkdRml /= 0.0d0) then
              radial = radial + &
                fcRij * fcRik * cutoff_fxn_prime(Rjk, rc, cutofffn_code, p_gamma) * dRjkdRml - &
                fcRijfcRikfcRjk * 2.0d0 * factor(3) / (rc ** 2.0d0) * Rjk * dRjkdRml
            end if
            radial = radial * (1.0d0 + g_gamma * costheta) ** g_zeta

            ridge(t_count) = ridge(t_count) + (angular + radial) * &
              2.0d0 ** (1.0d0 - g_zeta) * &
              exp(- factor(1) * (Rij / rc) ** 2.0d0 &
                - factor(2) * (Rik / rc) ** 2.0d0 &
                - factor(3) * (Rjk / rc) ** 2.0d0)
                
          end do 
        end do
        
      end if
    end do
  end do
  
end function calculate_g3_angle_prime

function calculate_g3_angle_partial_prime(&
  num_rscalars, r_scalars, &
  num_n_n, &
  ind_n_n, atomic_num_n_n, pos_n_n, &
  ind_n, atomic_num_n, pos_n, &
  ind_center, d, &
  atomicn_first, atomicn_second, &
  g_gamma, g_zeta, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
  
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer:: num_n_n
  integer, dimension(num_n_n):: ind_n_n
  integer, dimension(num_n_n):: atomic_num_n_n
  double precision, dimension(num_n_n, 3):: pos_n_n
  integer:: ind_n
  integer:: atomic_num_n
  double precision, dimension(3):: pos_n
  integer:: ind_center, d
  integer :: atomicn_first, atomicn_second
  double precision:: g_gamma, g_zeta
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision:: rc; integer:: cutofffn_code, p_gamma
    
  integer:: j, k, xyz, r_count, s_count, t_count
  double precision, dimension(3):: Rj, Rk
  double precision, dimension(3):: Rij_vector, Rik_vector
  double precision:: Rij, Rik, costheta

  double precision:: fcRij, fcRik, fcRijfcRik
  double precision:: dCosthetadRml
  double precision:: dRijdRml, dRikdRml

  double precision, dimension(6):: factor_ij
  double precision, dimension(6):: factor_ik

  integer:: match

  double precision:: angular, radial
  double precision, dimension(num_rscalars * num_rtypes) :: ridge

!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n_n
!f2py intent(in) :: ind_n_n, atomic_num_n_n, pos_n_n
!f2py intent(in) :: ind_n, atomic_num_n, pos_n
!f2py intent(in) :: ind_center, d
!f2py intent(in) :: atomicn_first, atomicn_second
!f2py intent(in) :: g_gamma, g_zeta
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, p_gamma, cutofffn_code
!f2py intent(out) :: ridge
  
  ridge = 0.0d0
  
  do j = 1, num_n_n
    do k = (j + 1), num_n_n
    
      if (atomicn_first .ne. -1) then
        match = neighbor_two(&
          atomic_num_n_n(j), atomic_num_n_n(k), &
          atomicn_first, atomicn_second)
      else
        match = 1
      end if
        
      if (match .eq. 1) then
        do xyz = 1, 3
          Rj(xyz) = pos_n_n(j, xyz)
          Rk(xyz) = pos_n_n(k, xyz)
          Rij_vector(xyz) = Rj(xyz) - pos_n(xyz)
          Rik_vector(xyz) = Rk(xyz) - pos_n(xyz)
        end do
        Rij = sqrt(dot_product(Rij_vector, Rij_vector))
        Rik = sqrt(dot_product(Rik_vector, Rik_vector))
        costheta = dot_product(Rij_vector, Rik_vector) / Rij / Rik
        dRijdRml = dRij_dRmd(ind_n, ind_n_n(j), pos_n, Rj, ind_center, d)
        dRikdRml = dRij_dRmd(ind_n, ind_n_n(k), pos_n, Rk, ind_center, d)
        dCosthetadRml = dCos_ijk_dR_ml(ind_n, &
          ind_n_n(j), ind_n_n(k), &
          pos_n, Rj, Rk, &
          ind_center, d)
          
        fcRij = cutoff_fxn(Rij, rc, cutofffn_code, p_gamma)
        fcRik = cutoff_fxn(Rik, rc, cutofffn_code, p_gamma)
        fcRijfcRik = fcRij * fcRik
        
        t_count = 0
        do r_count = 1, num_rtypes
          do s_count = 1, num_rscalars
          
            t_count = t_count + 1 
            call scalar(&
              r_scalar = r_scalars(s_count), &
              el_i = atomic_num_n, el_j = atomic_num_n_n(j), &
              interaction = 2, &
              r_type = r_types(r_count), &
              factor = factor_ij)
            call scalar(&
              r_scalar = r_scalars(s_count), &
              el_i = atomic_num_n, el_j = atomic_num_n_n(k), &
              interaction = 2, &
              r_type = r_types(r_count), &
              factor = factor_ik)

            ! 2 ** (1 - zeta) * (1 + gamma cos) ** zeta * & 
            ! fcRij * fcRik * exp(ij) * exp(ik) 
            angular = 0.0d0; radial = 0.0d0
            if (g_zeta == 1.0d0) then
              angular = g_zeta * g_gamma * dCosthetadRml
            else
              angular = g_zeta * g_gamma * dCosthetadRml * &
                (1.0d0 + g_gamma * costheta) ** (g_zeta - 1.0d0)
            end if
            angular = angular * fcRijfcRik

            if (dRijdRml /= 0.0d0) then
              radial = radial + &
                cutoff_fxn_prime(Rij, rc, cutofffn_code, p_gamma) * dRijdRml * fcRik - &
                fcRijfcRik * 2.0d0 * factor_ij(1) / (rc ** 2.0d0) * Rij * dRijdRml
            end if
            if (dRikdRml /= 0.0d0) then
              radial = radial + &
                fcRij * cutoff_fxn_prime(Rik, rc, cutofffn_code, p_gamma) * dRikdRml - &
                fcRijfcRik * 2.0d0 * factor_ik(1) / (rc ** 2.0d0) * Rik * dRikdRml
            end if
            radial = radial * (1.0d0 + g_gamma * costheta) ** g_zeta

            ridge(t_count) = ridge(t_count) + (angular + radial) * &
              2.0d0 ** (1.0d0 - g_zeta) * &
              exp(- factor_ij(1) * (Rij / rc) ** 2.0d0 &
                - factor_ik(1) * (Rik / rc) ** 2.0d0)
              
          end do 
        end do

      end if
    end do
  end do
  
end function calculate_g3_angle_partial_prime

function calculate_g4_prime(& 
  num_rscalars, r_scalars, &
  num_n_n, &
  ind_n_n, atomic_num_n_n, pos_n_n, &
  ind_n, atomic_num_n, pos_n, &
  ind_center, d, &
  atomicn_first, atomicn_second, atomicn_third, &
  num_rtypes, r_types, &
  rc, cutofffn_code, p_gamma) result(ridge)
  
  implicit none
  
  integer :: num_rscalars
  double precision, dimension(num_rscalars) :: r_scalars
  
  integer :: num_n_n
  integer, dimension(num_n_n) :: ind_n_n
  integer, dimension(num_n_n) :: atomic_num_n_n
  double precision, dimension(num_n_n, 3) :: pos_n_n
  integer :: ind_n
  integer :: atomic_num_n
  double precision, dimension(3) :: pos_n
  integer :: ind_center, d
  integer :: atomicn_first, atomicn_second, atomicn_third
  integer :: num_rtypes
  integer, dimension(num_rtypes) :: r_types
  double precision ::  rc; integer :: cutofffn_code, p_gamma
    
  integer :: j, k, l, xyz, r_count, s_count, t_count
  double precision, dimension(3) :: Rj, Rk, Rl
  double precision, dimension(3) :: Rij_vector, Rik_vector, Ril_vector
  double precision, dimension(3) :: Rjk_vector, Rjl_vector
  double precision, dimension(3) :: Rkl_vector
  
  double precision :: Rij, Rik, Ril
  double precision :: Rjk, Rjl
  double precision :: Rkl
  
  double precision:: fcRij, fcRik, fcRil, fcRijfcRikfcRil
  double precision :: dRijdRmd, dRikdRmd, dRildRmd
  double precision :: dRjkdRmd, dRjldRmd
  double precision :: dRkldRmd
    
  double precision, dimension(6) :: factor
  
  integer :: match
  
  double precision :: derivates
  double precision, dimension(num_rscalars * num_rtypes) :: ridge

!f2py intent(in) :: num_rscalars, r_scalars
!f2py intent(in) :: num_n_n
!f2py intent(in) :: ind_n_n, atomic_num_n_n, pos_n_n
!f2py intent(in) :: ind_n, atomic_num_n, pos_n
!f2py intent(in) :: ind_center, d
!f2py intent(in) :: atomicn_first, atomicn_second. atomicn_third
!f2py intent(in) :: num_rtypes, r_types
!f2py intent(in) :: rc, p_gamma, cutofffn_code
!f2py intent(out) :: ridge
  
  ridge = 0.0d0
  
  do j = 1, num_n_n
    do k = (j + 1), num_n_n
      do l = (k + 1), num_n_n

        if (atomicn_first .ne. -1) then
          match = neighbor_three(&
            atomic_num_n_n(j), atomic_num_n_n(k), atomic_num_n_n(l), &
            atomicn_first, atomicn_second, atomicn_third)
        else
          match = 1
        end if
          
        if (match .eq. 1) then 
          do xyz = 1, 3
            Rj(xyz) = pos_n_n(j, xyz) 
            Rk(xyz) = pos_n_n(k, xyz) 
            Rl(xyz) = pos_n_n(l, xyz) 
            Rij_vector(xyz) = pos_n_n(j, xyz) - pos_n(xyz)
            Rik_vector(xyz) = pos_n_n(k, xyz) - pos_n(xyz)
            Ril_vector(xyz) = pos_n_n(l, xyz) - pos_n(xyz)
            Rjk_vector(xyz) = pos_n_n(k, xyz) - pos_n_n(j, xyz)
            Rjl_vector(xyz) = pos_n_n(l, xyz) - pos_n_n(j, xyz)
            Rkl_vector(xyz) = pos_n_n(l, xyz) - pos_n_n(k, xyz)
          end do
          Rij = sqrt(dot_product(Rij_vector, Rij_vector))
          Rik = sqrt(dot_product(Rik_vector, Rik_vector))
          Ril = sqrt(dot_product(Ril_vector, Ril_vector))
          Rjk = sqrt(dot_product(Rjk_vector, Rjk_vector))
          Rjl = sqrt(dot_product(Rjl_vector, Rjl_vector))
          Rkl = sqrt(dot_product(Rkl_vector, Rkl_vector))
          dRijdRmd = dRij_dRmd(ind_n, ind_n_n(j), pos_n, Rj, ind_center, d)
          dRikdRmd = dRij_dRmd(ind_n, ind_n_n(k), pos_n, Rk, ind_center, d)
          dRildRmd = dRij_dRmd(ind_n, ind_n_n(l), pos_n, Rl, ind_center, d)
          dRjkdRmd = dRij_dRmd(ind_n_n(j), ind_n_n(k), Rj, Rk, ind_center, d)
          dRjldRmd = dRij_dRmd(ind_n_n(j), ind_n_n(l), Rj, Rl, ind_center, d)
          dRkldRmd = dRij_dRmd(ind_n_n(k), ind_n_n(l), Rk, Rl, ind_center, d)
            
          fcRij = cutoff_fxn(Rij, rc, cutofffn_code, p_gamma)
          fcRik = cutoff_fxn(Rik, rc, cutofffn_code, p_gamma)
          fcRil = cutoff_fxn(Ril, rc, cutofffn_code, p_gamma)
          fcRijfcRikfcRil = fcRij * fcRik * fcRil
          
          t_count = 0
          do r_count = 1, num_rtypes
            do s_count = 1, num_rscalars
            
              t_count = t_count + 1
              call scalar(&
                r_scalar = r_scalars(s_count), &
                el_i = atomic_num_n, &
                el_j = atomic_num_n_n(j), &
                el_k = atomic_num_n_n(k), &
                el_l = atomic_num_n_n(l), &
                interaction = 4, &
                r_type = r_types(r_count), &
                factor = factor)
              ! fcRij * fcRik * fcRil  * &
              ! exp( - (factor(1) * (Rij / rc) ** 2)
              ! exp( - (factor(2) * (Rik / rc) ** 2)
              ! exp( - (factor(3) * (Ril / rc) ** 2)
              ! exp( - (factor(4) * (Rjk / rc) ** 2)
              ! exp( - (factor(5) * (Rjl / rc) ** 2)
              ! exp( - (factor(6) * (Rkl / rc) ** 2)
              derivates = 0.0d0
              if (dRijdRmd /= 0.0d0) then
                derivates = derivates + &
                  cutoff_fxn_prime(Rij, rc, cutofffn_code, p_gamma) * dRijdRmd * fcRik * fcRil - &
                  fcRijfcRikfcRil * 2.0d0 * factor(1) / (rc ** 2.0d0) * Rij * dRijdRmd
              end if
              if (dRikdRmd /= 0.0d0) then
                derivates = derivates + &
                  fcRij * cutoff_fxn_prime(Rik, rc, cutofffn_code, p_gamma) * dRikdRmd * fcRil - &
                  fcRijfcRikfcRil * 2.0d0 * factor(2) / (rc ** 2.0d0) * Rik * dRikdRmd
              end if
              if (dRildRmd /= 0.0d0) then
                derivates = derivates + &
                  fcRij * fcRik * cutoff_fxn_prime(Ril, rc, cutofffn_code, p_gamma) * dRildRmd - &
                  fcRijfcRikfcRil * 2.0d0 * factor(3) / (rc ** 2.0d0) * Ril * dRildRmd
              end if
              if (dRjkdRmd /= 0.0d0) then
                derivates = derivates - &
                  fcRijfcRikfcRil * 2.0d0 * factor(4) / (rc ** 2.0d0) * Rjk * dRjkdRmd
              end if
              if (dRjldRmd /= 0.0d0) then
                derivates = derivates - &
                  fcRijfcRikfcRil * 2.0d0 * factor(5) / (rc ** 2.0d0) * Rjl * dRjldRmd
              end if
              if (dRkldRmd /= 0.0d0) then
                derivates = derivates - &
                  fcRijfcRikfcRil * 2.0d0 * factor(6) / (rc ** 2.0d0) * Rkl * dRkldRmd
              end if
            
              ridge(t_count) = ridge(t_count) + derivates * &
                exp(- factor(1) * (Rij / rc) ** 2.0d0 &
                  - factor(2) * (Rik / rc) ** 2.0d0 &
                  - factor(3) * (Ril / rc) ** 2.0d0 &
                  - factor(4) * (Rjk / rc) ** 2.0d0 &
                  - factor(5) * (Rjl / rc) ** 2.0d0 &
                  - factor(6) * (Rkl / rc) ** 2.0d0)
                
            end do 
          end do
          
        end if  
      end do
    end do
  end do
 
end function calculate_g4_prime

end module char_zeta

