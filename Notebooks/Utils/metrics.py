import sklearn.metrics as metrics
import scipy.stats
import numpy as np

def compute_metrics(y_predicted, y_test, y_train, fc_scale=10):
    '''
    Function to compute regression evaluation metrics.
    https://stats.stackexchange.com/questions/142873/how-to-determine-the-accuracy-of-regression-which-measure-should-be-used
    https://stats.stackexchange.com/questions/131267/how-to-interpret-error-measures
    https://homepages.inf.ed.ac.uk/ckiw/postscript/Chalupka2011diss.pdf
    
    Inputs:
    -y_predicted = predited y values from model, np.array(n_samples,1)
    -y_test = the ground truth labels for the prediction, np.array(n_samples,1)
    -y_train = the truth values from the test set, used in smse, np.array(n_samples,1)
    -fc_scale = is the log base used to transform protein measurments for fold change, default is 10
    
    Return:
    -evaluation (dict of metrics)
    
    Notes:
    -naive/normalizing method would be to guess the avg training target for every test case
    -For explained variance and r2 differences
     https://stats.stackexchange.com/questions/241329/difference-between-coefficient-of-determination-and-explained-variance
    
    Metrics comptued:
    -MAE
    -MSE
    -sMSE
    -R2
    -Explained variance score
    -Spearman rho
    -Pearson r
    '''
    assert isinstance(y_predicted,np.ndarray) and len(y_predicted.shape)==2
    assert isinstance(y_test,np.ndarray) and len(y_test.shape)==2
    assert isinstance(y_train,np.ndarray) and len(y_train.shape)==2
    assert y_predicted.shape==y_test.shape
    
    evaluation = {}
    evaluation['mae'] = metrics.mean_absolute_error(y_test,y_predicted)
    evaluation['mse'] = metrics.mean_squared_error(y_test,y_predicted)
    #Denominator is just normalized mse (using avg train output as guess)
    evaluation['smse'] = metrics.mean_squared_error(y_test,y_predicted)/(y_test.var() + (y_train.mean()-y_test.mean())**2)
    evaluation['r2'] = metrics.r2_score(y_test,y_predicted)
    evaluation['evs'] = metrics.explained_variance_score(y_test,y_predicted)
    
    rho, rho_p = scipy.stats.spearmanr(y_test.squeeze(),y_predicted.squeeze())
    evaluation['spearmanrho'] = rho
    evaluation['spearmanrho_p'] = rho_p
    
    r, r_p = scipy.stats.pearsonr(y_test.squeeze(),y_predicted.squeeze())
    evaluation['pearsonr'] = r
    evaluation['pearsonr_p'] = r_p
    
    abs_difference = np.abs(y_test - y_predicted)
    evaluation['median_abs_fc'] = fc_scale**np.median(abs_difference)
    evaluation['mean_abs_fc'] = fc_scale**np.mean(abs_difference)
    
    return evaluation