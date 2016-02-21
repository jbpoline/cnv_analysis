from __future__ import print_function
from __future__ import division

import numpy as np
# import scipy.stats as sst
# import matplotlib.pyplot as plt
# import os
# import os.path as osp
# import six


def str2floats(value, sep=',', comma2point=False):
    """
    takes a string with one or several p-values 
    returns:
    combined score dangerosite
    """
    if comma2point:
        value = value.replace(',', '.')
    lst = value.split(sep)
    lst = [l.rstrip().lstrip() for l in lst if (l != '')]
    try:
        return [float(f) for f in lst]
    except:
        print("value: ", value, "value.split(sep): ", lst)
        raise ValueError


def pH1_with_apriori(alpha, pi=0.5, beta=.8):
    """
    Returns the probability of H1 if we know power and a priori
    """
    return (1 - beta) * pi / ((1 - beta) * pi + alpha * (1 - pi))

def pH1_simple(alpha, pi=None, beta=None):
    """
    pi and beta unsused:
    """
    return (1. - alpha)


#-------------------
defaultarg = {'pi':0.5, 'beta':0.8}
#-------------------

def danger_score(pvalues, pval2score=pH1_simple, argvals=defaultarg):
    """
    take a series of pvalue (one line), return a dangerosity score
    """
    
    pvalues = np.asarray(str2floats(pvalues))
    
    assert np.all( np.logical_and(pvalues <= 1., pvalues >=0)), \
            [pvalues, np.logical_and((pvalues <= 1.), (pvalues >=0))]
    
    """ how to combine the p-values?  
    I want a score that scales with the number of gene affected. Per gene, if the p is zero, 
    I want my score to be 1. and if the p value is close to 1, the score should be 0.
    To start with, I'm setting the score to be 1-p:
        
        res = (1. - pvalues).sum()

    The problem with this approach: I would like something closer to the probability of 
    CNV to be dangerous, which requires an apriori on the probablility that the CNV is dangerous, 
    and the power of the test (leading to the p value that we are using).
    
        p(CNV dangerous) =  (1 - beta) * pi / ((1 - beta) * pi + alpha * (1 - pi))
    
    """
    res = np.asarray([pval2score(p, **argvals) for p in pvalues]).sum()
    return res


def test_danger_score(tests, verbose=True):
    """
    tests: the list of strings that look like "danger"
    results: expected results
    """
    for idx, test in enumerate(tests):
        pstring = test[0]
        func = test[1]
        arg = test[2]
        result = test[3]
        status = True
        
        try:
            dger = danger_score(pstring, func, arg)
            if verbose: print('dger:',dger)
            if abs(dger - result) < np.finfo('float').eps:
                pass
            else:
                if verbose: print(test, " is failing", dger, results[idx])
                status = False
        except:
            if result == "Failed":
                pass
            else:
                if verbose: print(test, " is failing", results[idx])
                status = False
                
    return(status)


def _test_danger_score_1():

    assert abs(danger_score('1., 0., .5', pH1_with_apriori) - 1.45238095238) < 1.e-6
    assert abs(danger_score('1.', pH1_with_apriori) - 0.166666666667) < 1.e-6 
    assert abs(danger_score('1.', pH1_with_apriori, {'pi':0.5, 'beta':0.8}) - 0.166666666667) < 1.e-6 
    return True

def _test_danger_score():

    one_beta_pi = (1. - defaultarg['beta'])*defaultarg['pi']
    pH1_res_1 = one_beta_pi / (one_beta_pi + 1. - defaultarg['pi']) 
    assert abs(pH1_res_1 - pH1_with_apriori(1.)) < np.finfo('float').eps

    tests = [('1., 0., .5', pH1_simple, defaultarg,          1.5),
             ('.2, -.01',    pH1_simple, defaultarg,         'Failed'),
             ('0.1, 0.1, .8',pH1_simple, defaultarg,         2.), 
             ('.0',          pH1_simple, defaultarg,         1.),
             ('.0',          pH1_with_apriori, defaultarg,   1.),
             ('1.',          pH1_with_apriori, defaultarg,   pH1_res_1)]

    assert test_danger_score(tests, verbose=False)
    return True
