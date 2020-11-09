module neuralnetwork

implicit none

type:: two_d_placeholder
  sequence
  double precision, allocatable:: twodarray(:,:)
end type two_d_placeholder

type:: one_d_placeholder
  sequence
  double precision, allocatable:: onedarray(:)
end type one_d_placeholder

type:: structured_parameter
  sequence
  type(two_d_placeholder), allocatable:: weights(:)
end type structured_parameter

contains

function calculate_atomic_energy(num_gs, &
  fingerprint, &
  min_fingerprints, &
  max_fingerprints, &
  num_parameters, &
  parameters, &
  no_layers, &
  no_nodes)
  
  implicit none
  
  integer:: num_parameters, num_gs
  double precision:: fingerprint(num_gs)
  double precision:: min_fingerprints(num_gs)
  double precision:: max_fingerprints(num_gs)
  double precision:: parameters(num_parameters)
  integer:: no_layers
  integer:: no_nodes(no_layers)
  
  double precision:: calculate_atomic_energy
  
  integer:: parameter_ind
  integer:: afp_ind, layer_ind, node_ind
  integer:: num_rows, num_cols
  
  integer:: q, p
  integer:: aux_ind
  character(len=10):: aux_char
  
  ! net is the value before activation
  double precision, allocatable:: net(:)
  type(one_d_placeholder), allocatable:: o(:), ohat(:)
  type(structured_parameter):: unraveled_parameters
  
  ! used for storing the scaled fingerprint
  double precision:: fingerprint_(num_gs)
    
  do afp_ind = 1, num_gs
    if ((max_fingerprints(afp_ind) - min_fingerprints(afp_ind)) .GT. &
      (10.0d0 ** (-8.0d0))) then
      fingerprint_(afp_ind) = &
        (fingerprint(afp_ind) - min_fingerprints(afp_ind)) / &
        (max_fingerprints(afp_ind) - min_fingerprints(afp_ind))
    else
      fingerprint_(afp_ind) = fingerprint(afp_ind)
    endif
  end do
  
  ! structure the model parameter
  allocate(unraveled_parameters%weights(no_layers - 1))
  parameter_ind = 0
    
  do layer_ind = 1, no_layers - 1
    ! because there is also a bias node
    num_rows = no_nodes(layer_ind) + 1
    num_cols = no_nodes(layer_ind + 1)
    allocate(unraveled_parameters%weights(layer_ind)%twodarray(num_rows, num_cols))
    
    do p = 1, num_rows
      do q = 1, num_cols
        unraveled_parameters%weights(layer_ind)%twodarray(p, q) = &
        parameters(parameter_ind + (p - 1) * num_cols + q)
      end do
    end do
    parameter_ind = parameter_ind + num_rows * num_cols
  end do
  
  allocate(o(no_layers))
  allocate(ohat(no_layers))
  
  allocate(o(1)%onedarray(num_gs))
  allocate(ohat(1)%onedarray(num_gs + 1))
  do afp_ind = 1, num_gs
    o(1)%onedarray(afp_ind) = fingerprint_(afp_ind)
  end do
  
  do layer_ind = 1, no_layers - 1
    do node_ind = 1, size(unraveled_parameters%weights(layer_ind)%twodarray, dim=1) - 1
      ohat(layer_ind)%onedarray(node_ind) = o(layer_ind)%onedarray(node_ind)
    end do
    ohat(layer_ind)%onedarray(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=1)) = 1.0d0
    
    allocate(net(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2)))
    allocate(o(layer_ind + 1)%onedarray(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2)))
    allocate(ohat(layer_ind + 1)%onedarray(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2) + 1))
    
    ! hardcode matrix multiplication
    ! input(row vec) @ W^(l, l + 1) ==> the output vec (row vec)
    ! so we iterate through each entry of the output row vec
    do q = 1, size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2)
      net(q) = 0.0d0
      do p = 1, size(unraveled_parameters%weights(layer_ind)%twodarray, dim=1)
        net(q) =  net(q) + &
          ohat(layer_ind)%onedarray(p) * unraveled_parameters%weights(layer_ind)%twodarray(p, q)
      end do
      
      if (layer_ind .lt. no_layers - 1) then 
        o(layer_ind + 1)%onedarray(q) = tanh(net(q))
        ohat(layer_ind + 1)%onedarray(q) = o(layer_ind + 1)%onedarray(q)
      elseif (layer_ind .eq. no_layers - 1) then
        o(layer_ind + 1)%onedarray(q) = net(q)
        ohat(layer_ind + 1)%onedarray(q) = o(layer_ind + 1)%onedarray(q)
      end if 
      
    end do
    ohat(layer_ind + 1)%onedarray(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2) + 1) =  1.0d0
    deallocate(net)
  end do
  
  calculate_atomic_energy = o(layer_ind)%onedarray(1)
  
  do p = 1, size(o)
    deallocate(o(p)%onedarray)
  end do
  deallocate(o)
  
  do p = 1, size(ohat)
    deallocate(ohat(p)%onedarray)
  end do
  deallocate(ohat)
  
  do layer_ind = 1, no_layers - 1
    deallocate(unraveled_parameters%weights(layer_ind)%twodarray)
  end do
  deallocate(unraveled_parameters%weights)
end function calculate_atomic_energy

function calculate_force(num_gs, &
  fingerprint, &
  fingerprintprime, &
  min_fingerprints, &
  max_fingerprints, &
  num_parameters, &
  parameters, &
  no_layers, &
  no_nodes)
  
  implicit none
  
  integer:: parameter_ind
  integer:: num_gs
  integer:: num_parameters
  integer:: afp_ind, layer_ind, node_ind
  integer:: num_rows, num_cols
  integer:: p, q
  
  double precision:: fingerprint(num_gs)
  double precision:: fingerprintprime(num_gs)
  double precision:: min_fingerprints(num_gs)
  double precision:: max_fingerprints(num_gs)
  double precision:: parameters(num_parameters)
  integer:: no_layers
  integer:: no_nodes(no_layers)
  
  double precision:: calculate_force

  double precision, allocatable:: net(:)
  type(one_d_placeholder), allocatable:: o(:), ohat(:)
  type(one_d_placeholder), allocatable:: doutputs_dinputs(:)
  type(structured_parameter):: unraveled_parameters
  
  ! hold the scaled fingerprint and fingerprint primes
  double precision:: fingerprint_(num_gs)
  double precision:: fingerprintprime_(num_gs)
  
  do afp_ind = 1, num_gs
    if ((max_fingerprints(afp_ind) - min_fingerprints(afp_ind)) .GT. &
      (10.0d0 ** (-8.0d0))) then
      fingerprint_(afp_ind) = &
        (fingerprint(afp_ind) - min_fingerprints(afp_ind)) / &
        (max_fingerprints(afp_ind) - min_fingerprints(afp_ind))
    else
      fingerprint_(afp_ind) = fingerprint(afp_ind)
    endif
  end do
  
  ! scaling fingerprintprimes
  do afp_ind = 1, num_gs
    if ((max_fingerprints(afp_ind) - min_fingerprints(afp_ind)) .GT. &
      (10.0d0 ** (-8.0d0))) then
      fingerprintprime_(afp_ind) = &
        fingerprintprime(afp_ind) / (max_fingerprints(afp_ind) - min_fingerprints(afp_ind))
    else
      fingerprintprime_(afp_ind) = fingerprintprime(afp_ind)
    endif
  end do
  
  ! structure the model parameter
  allocate(unraveled_parameters%weights(no_layers - 1))
  parameter_ind = 0
  do layer_ind = 1, no_layers - 1
    ! because there is also a bias node
    num_rows = no_nodes(layer_ind) + 1
    num_cols = no_nodes(layer_ind + 1)
    allocate(unraveled_parameters%weights(layer_ind)%twodarray(num_rows, num_cols))
    do p = 1, num_rows
      do q = 1, num_cols
        unraveled_parameters%weights(layer_ind)%twodarray(p, q) = &
        parameters(parameter_ind + (p - 1) * num_cols + q)
      end do
    end do
    parameter_ind = parameter_ind + num_rows * num_cols
  end do
  
  allocate(o(no_layers))
  allocate(ohat(no_layers))
  
  allocate(o(1)%onedarray(num_gs))
  allocate(ohat(1)%onedarray(num_gs + 1))
  do afp_ind = 1, num_gs
    o(1)%onedarray(afp_ind) = fingerprint_(afp_ind)
  end do
  
  do layer_ind = 1, no_layers - 1
    do node_ind = 1, size(unraveled_parameters%weights(layer_ind)%twodarray, dim=1) - 1
      ohat(layer_ind)%onedarray(node_ind) = o(layer_ind)%onedarray(node_ind)
    end do
    ohat(layer_ind)%onedarray(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=1)) = 1.0d0
    
    allocate(net(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2)))
    allocate(o(layer_ind + 1)%onedarray(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2)))
    allocate(ohat(layer_ind + 1)%onedarray(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2) + 1))
    
    ! hardcode matrix multiplication
    ! input(row vec) @ W^(l, l + 1) ==> the output vec (row vec)
    ! so we iterate through each entry of the output row vec
    do q = 1, size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2)
      net(q) = 0.0d0
      do p = 1, size(unraveled_parameters%weights(layer_ind)%twodarray, dim=1)
        net(q) =  net(q) + &
          ohat(layer_ind)%onedarray(p) * unraveled_parameters%weights(layer_ind)%twodarray(p, q)
      end do
      
      if (layer_ind .lt. no_layers - 1) then 
        o(layer_ind + 1)%onedarray(q) = tanh(net(q))
        ohat(layer_ind + 1)%onedarray(q) = o(layer_ind + 1)%onedarray(q)
      elseif (layer_ind .eq. no_layers - 1) then
        o(layer_ind + 1)%onedarray(q) = net(q)
        ohat(layer_ind + 1)%onedarray(q) = o(layer_ind + 1)%onedarray(q)
      end if 
    end do
    ohat(layer_ind + 1)%onedarray(size(unraveled_parameters%weights(layer_ind)%twodarray, dim=2) + 1) =  1.0d0
    deallocate(net)
  end do
  
  ! then it comes something interesting
  allocate(doutputs_dinputs(no_layers))
  allocate(doutputs_dinputs(1)%onedarray(num_gs))
  
  do afp_ind = 1, num_gs
    doutputs_dinputs(1)%onedarray(afp_ind) = fingerprintprime_(afp_ind)
  end do
  
  do layer_ind = 1, no_layers - 1
    allocate(doutputs_dinputs(layer_ind + 1)%onedarray(size(o(layer_ind + 1)%onedarray)))
    do p = 1, size(o(layer_ind + 1)%onedarray)
      doutputs_dinputs(layer_ind + 1)%onedarray(p) = 0.0d0
      ! perform forward pass
      do q = 1, size(o(layer_ind)%onedarray)
        doutputs_dinputs(layer_ind + 1)%onedarray(p) = & 
          doutputs_dinputs(layer_ind + 1)%onedarray(p) + &
          doutputs_dinputs(layer_ind)%onedarray(q) * &
          unraveled_parameters%weights(layer_ind)%twodarray(q, p)
      end do
      ! the derivate through the layer_ind + 1
      if (layer_ind .lt. no_layers - 1) then 
        doutputs_dinputs(layer_ind + 1)%onedarray(p) = & 
          doutputs_dinputs(layer_ind + 1)%onedarray(p) * & 
          (1.0d0 - o(layer_ind + 1)%onedarray(p) * o(layer_ind + 1)%onedarray(p))
      end if  
    end do
  end do
  
  calculate_force = - 6.9d0 * doutputs_dinputs(layer_ind)%onedarray(1)
  
  
  do p = 1, size(o)
    deallocate(o(p)%onedarray)
  end do
  deallocate(o)
  do p = 1, size(ohat)
    deallocate(ohat(p)%onedarray)
  end do
  deallocate(ohat)
  
  do p = 1, size(doutputs_dinputs)
    deallocate(doutputs_dinputs(p)%onedarray)
  end do
  deallocate(doutputs_dinputs)
  
  do layer_ind = 1, no_layers - 1
    deallocate(unraveled_parameters%weights(layer_ind)%twodarray)
  end do
  deallocate(unraveled_parameters%weights)
end function calculate_force

end module neuralnetwork
