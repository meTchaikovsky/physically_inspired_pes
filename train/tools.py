
import numpy as np


def weight_bias_initializer(hidden_layers, input_dim):
    
    weight, bias = {}, {}
    nn_structure = [input_dim] + [layer for layer in hidden_layers] + [1]
    for count in range(len(nn_structure) - 1):
        epsilon = np.sqrt(6. / (nn_structure[count] + nn_structure[count + 1]))
        normalized_arg_range = 2. * epsilon

        weight[count + 1] = \
            np.random.random(
                    (nn_structure[count], nn_structure[count + 1])
                    ).astype('float32') * \
            normalized_arg_range - normalized_arg_range / 2.
        bias[count + 1] = \
            np.random.random(nn_structure[count + 1]).astype('float32') * \
            normalized_arg_range - normalized_arg_range / 2.
            
    return weight, bias

def check_convergence(iteration, convergence, func_vals, file_name):
        
    converged = [
            convergence['energy rmse'] > func_vals[1],
            convergence['energy mae'] > func_vals[2],
            convergence['energy maxresid'] > func_vals[3],
             convergence['force rmse'] > func_vals[4],
            convergence['force mae'] > func_vals[5], 
            convergence['force maxresid'] > func_vals[6]
            ]
    
    verbose_converge = [' C' if item else ' -' for item in converged]    
    verbose_converge.insert(0, ' ')
    output = []
    for val, result in zip(func_vals, verbose_converge):
        output.append('%.6f %s' % (val, result))
    output.insert(0, iteration)
    
    recorder = open(file_name, 'a')
    recorder.write('{:^10}{:^15}{:^15}{:^15}{:^15}{:^15}{:^15}{:^15}\n'.format(*output))
    recorder.flush(); recorder.close()
    
    return np.all(converged)

class Convergence_Occurred(Exception):
    pass

class Early_stopping(Exception):
    pass