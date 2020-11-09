import pickle
from kyzhang.train.model import neural_network

db_path = '/THL7/home/lcyin/KYZhang/bp_nn_pes/playground/create_db/database'
species_partitions = [['Li', 'N']]
hidden_layers = [[100, 80, 60]]
activation = 'tanh'

model_path = './model'
test_size = 0.1
random_state = 1231324
num_iterations = 2000
earlystopping_patience = 15
xavier_initializer = True
retrain = False

train_params = {
        'rate': 0, 
        'force_coefficient': 1, 
        'energy_coefficient': 1, 
        'regularization_strength': 0
        }
convergence = {
        'energy rmse': 0.03, 'energy maxresid': 0.03, 'energy mae': 0.02, 
        'force rmse': 0.03, 'force maxresid': 0.03 , 'force mae': 0.8
        }

##############################################################
with open('%s/des_summary' % db_path, 'rb') as source:
    des_dict = pickle.load(source)
fp_dim = des_dict['fp_dim']
nn_model = neural_network(
        species_partitions = species_partitions,
        fp_dim = fp_dim, hidden_layers = hidden_layers,
        activation = activation,
        )
nn_model.lbfgs(
        model_path=model_path, 
        db_path = db_path, portion = -1, test_size = test_size, random_state = random_state,
        num_iterations=num_iterations, earlystopping_patience = earlystopping_patience,
        train_params = train_params,
        convergence = convergence,
        xavier_initializer = xavier_initializer,
        retrain = retrain
        )
