from __future__ import print_function
from __future__ import division

import numpy as np
import six
# from datetime import datetime
# import scipy.stats as sst
# import matplotlib.pyplot as plt
# import os
# import os.path as osp


def str2floats(value, sep=',', comma2point=False, val_len=None):
    """
    takes a string (in this use case, the string contains one or several
    values (eg, p-values): 

    parameters: 
    ------------
    value: string
        the string 
    sep: cut the string with this separator 
    returns:
    -------
        list of float

    >>> str2floats('1,3 4,2, 5', sep=' ', comma2point=True) 
    [1.3, 4.2, 5.0]
    >>> str2floats('1.3, 4.2, 5')
    [1.3, 4.2, 5.0]
    """

    # somme strings have , for the . 
    if comma2point and sep==',':
        raise ValueError, "comma2point {} and sep :' {} ' not compatible".format(
                                                            comma2point, sep)
    if comma2point:
        value = value.replace(',', '.')

    lst = value.split(sep)
    lst = [l.rstrip().lstrip() for l in lst if (l != '')]
    if val_len is not None:
        assert len(lst) == val_len, "val_len not ok: {}".format(lst)

    try:
        return [float(f) for f in lst]
    except:
        print("cannot make a float of : ", value, "value.split(sep): ", lst)
        raise ValueError


def get_one_col_val(i_col, line, comma2point=True, sep=' ', exclude=None, val_len=None):
    """
    extract the column index in line 
    i_col: index of column
    line: string
    comma2point: bool
        replace commas by points 
    sep: string
        within this column (string) split string with sep
    """

    # assume col_split is \t
    col_split = '\t'
    # first line is the name of the columns
    str_val = line.split(col_split)[i_col]
    
    if isinstance(exclude, six.string_types):
        # something to exclude
        if str_val == exclude:
            return False

    val = str2floats(str_val, sep=sep, comma2point=comma2point, val_len=val_len)[0] 
    
    return val


def get_one_col_str(i_col, line):
    """
    extract the column index in line 
    i_col: index of column
    line: string
    """
    col_split = '\t'
    return line.split(col_split)[i_col]
    
def get_col_vals(col_name, array, comma2point=True, sep=' ', exclude=None, val_len=None):
    """
    extract the column with name col_name in array
    
    parameters
    ----------
    col_name: string
        the first elt of array should have this col_name separated by \t
    array: numpy array
        results of reading with np.loadtxt, array.dtype is 'S????' (strings)
    comma2point: bool
    exclude: string
        the string to exclude
    val_len: int or None
        if int, indicates max number of values can be read in this columns
        eg: 1 will constrain to 1 value

    returns
    -------
    array of values for this column 

    """
    # assume col_split is \t
    col_split = '\t'
    # first line is the name of the columns
    line0 = array[0].split(col_split)
    i_col = line0.index(col_name)
    str_vals = np.asarray([line.split(col_split)[i_col] for line in array[1:]])
    
    if isinstance(exclude, six.string_types):
        # something to exclude
        to_keep = str_vals != exclude
    else:
        to_keep = np.full_like(str_vals, True, dtype=bool)
    
    print('array read shape:', to_keep.shape, 'number to keep: ', to_keep.sum()) 
    
    #  
    vals = [str2floats(val, sep=sep, comma2point=comma2point, val_len=val_len)[0] \
            for val in str_vals[to_keep]]
    
    assert len(vals) > 0, "nothing in this column : {}, empty result".format(vals)
    
    return np.asarray(vals)



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


#-------------------------------------
defaultarg = {'pi':0.5, 'beta':0.8}
#-------------------------------------

def danger_score(pvalues, pval2score=pH1_simple, sep=' ', comma2point=True, 
                                                            argvals=defaultarg):
    """
    take a series of pvalue (one line), return a dangerosity score

    How to combine the p-values?  
    I want a score that scales with the number of gene affected. Per gene, if the p is zero, 
    I want my score to be 1. and if the p value is close to 1, the score should be 0.
    To start with, I'm setting the score to be 1-p:
        
        res = (1. - pvalues).sum()

    The problem with this approach: I would like something closer to the probability of 
    CNV to be dangerous, which requires an apriori on the probablility that the CNV is dangerous, 
    and the power of the test: 
    
        p(CNV dangerous) =  (1 - beta) * pi / ((1 - beta) * pi + alpha * (1 - pi))
        (see power notebook)
        where:
        beta    == power of the test,
        pi      == apriori probability that the cnv is dangerous
        alpha   == risk of false positive (the p-value given by the test)


    parameters:
    -----------
    pvalue: string 
        list of p values, one per gene, separated by space and with 
        comma instead of point.
    pval2score: callable
        one of (pH1_simple, pH1_with_apriori), default to pH1_simple
        the function that takes a p value as input and spits out a score
    argvals: dict
        default arguments for pval2score function
        default: defaultarg = {'pi':0.5, 'beta':0.8}
    """
   
    # 
    pvalues = np.asarray(str2floats(pvalues, sep=sep, comma2point=comma2point))
   
    # check these values are in [0..1]
    assert np.all( np.logical_and(pvalues <= 1., pvalues >=0)), \
            [pvalues, np.logical_and((pvalues <= 1.), (pvalues >=0))]
    
    # sum transform each p value to score and sum
    res = np.asarray([pval2score(p, **argvals) for p in pvalues]).sum()
    return res


def danger_score_fl(pvalues, pval2score=pH1_simple, argvals=defaultarg):
    """
    pvalues: one float

    """
    raise
    

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
                if verbose: print(test, " is failing", dger, result[idx])
                status = False
        except:
            if result == "Failed":
                pass
            else:
                if verbose: print(test, " is failing", result[idx])
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





#===============================================================================
# functions to go from a cnv score to a probability of this cnv to be real
#===============================================================================

from scipy.interpolate import UnivariateSpline as spl
from collections import OrderedDict
import numbers

def _build_dict_prob_cnv():
    """
    construct the dictionary used to build the 
    function that will take a cnv score and spits out 
    a proba of being true cnv

    
    """
    # Values come from Guillaumes H. 

    prob_cnv_vrai_str = \
                    "15	20	80	25\n" + \
                    "17,5	24	76	25\n" + \
                    "20	36	64	25\n" + \
                    "22,5	44	56	25\n" + \
                    "25	83	17	60\n" + \
                    "27,5	77	23	44\n" + \
                    "30	95	5	60"

    p_cnv = OrderedDict()
    for line in prob_cnv_vrai_str.splitlines():
        line = line.split('\t')
        #print(line)
        # only consider the index 0 and 1 values of line:
        p_cnv[float(line[0].replace(',','.'))] = float(line[1])/100.
    
    return p_cnv


def aff(bi, bs, vbi, vbs, x0, y0):
    """
    returns the affine function
    """
    a = (vbs - vbi)/(bs - bi)
    b = y0  - a*x0
    
    def this_aff(x):
        return a*x + b
    
    return this_aff

def create_score2prob_lin_piecewise(p_cnv):

    # create a dict that has the affine function
    pf_cnv = OrderedDict({})
    kcnv = p_cnv.keys()

    for idx,k in enumerate(kcnv[:-1]):
        k1 = kcnv[idx+1]
        pf_cnv[k] = aff(k, k1, p_cnv[k], p_cnv[k1], k, p_cnv[k])
    
    # testing just this 
    def test_pf_cnv():
        assert pf_cnv[15](15) == .20
        assert abs(pf_cnv[15](17.5) - pf_cnv[17.5](17.5)) < np.finfo(float).eps
        assert abs(pf_cnv[17.5](20) - pf_cnv[20.](20.)) < np.finfo(float).eps
        assert abs(pf_cnv[20.](22.5) - pf_cnv[22.5](22.5)) < np.finfo(float).eps
        #assert pf_cnv[30](45) == 1.
        return True
    test_pf_cnv()
    
    # function linear piecewise:
    def sc2prob(sc):
        """
        returns the proba that cnv with score sc is a cnv; 
        """
        # 
        assert sc >= 15.
        kcnv = pf_cnv.keys()
        
        # if greater than the last key, returns 1.0
        if sc >= kcnv[-1]:
            return 1.0
        
        idx = 0
        while sc > kcnv[idx+1]: idx += 1
        
        return pf_cnv[kcnv[idx]](sc)
    
    return sc2prob
    



def create_score2prob_lin(p_cnv):
    """
    create the function that will transform a cnv score in a p-value 
    of being a cnv

    Just a linear approximation
    """


    kpcnv = p_cnv.keys()
    lin_funct = spl(kpcnv, p_cnv.values(), k=1)
    #x = np.arange(0,50,1)
    #plt.plot(x, lin_funct(x), '-', p_cnv.keys(), p_cnv.values(), '+')

    (x1, x2) = (kpcnv[0], kpcnv[-1])
    a = (lin_funct(x2) - lin_funct(x1))/(x2 - x1)
    b = lin_funct(x2) - a*x2
    x_where_y_is_1 = (1. - b)/a
    x_where_y_is_0 = (-b)/a

    
    def sc2prob(x):
        # x should not be a string
        # assert not isinstance(x, six.string_types), "{} is str".format(x)
        # but should be a number
        assert isinstance(x, numbers.Real), \
                                "{} is not a number: type {}".format(x,type(x))

        if x < x_where_y_is_0:
            return 0.
        elif x < x_where_y_is_1:
            return float(lin_funct(x))
        else:
            return 1.
    
    return sc2prob


def process_scores(scores):
    """
    scores: a list of strings 

    return the cleaned scores : floats
    """

    # 1. comma to point, ' ' is separator
    scores = np.asarray([str2floats(s, comma2point=True, sep=' ')[0] 
                                                             for s in scores])
    # 2. if the score is zero, put it to maximum value : 
    # means the CNV has a maximum score
    clean_score = scores
    clean_score[scores==0] = scores.max()

    return clean_score

def process_one_score(score, max_score):

    assert isinstance(score, six.string_types),\
                        "{} is of type {}".format(score, type(score))

    score = str2floats(score, comma2point=True, sep=' ')
    assert len(score) == 1
    score = score[0]
    if score == 0:
        score = max_score

    return score
        

def make_uiid(tsv_arr, columns_names, first_line=None):
    """
    columns_names: list of column names with which the uiid are constructed
    returns a list if given the full array, 
    returns a string if given only one line (and the first line)
    """
    
    chr2rm = ''.join([',', '.', ' ']) # others: ['!', '?', ...] ?
    
    if isinstance(tsv_arr, six.string_types) and (first_line != None):
        # assume we are given one line of the tsv array
        indexes = [first_line.split('\t').index(colname) for colname in columns_names]
        uiid = '_'.join([(tsv_arr.split('\t')[ind]).translate(None,chr2rm) for ind in indexes])
        
    else: 
        uiid = []
        indexes = [tsv_arr[0].split('\t').index(colname) for colname in columns_names]
        for line in tsv_arr[1:]:
            ll = line.split('\t')
            uiid.append('_'.join([ll[ind] for ind in indexes]))

    return uiid














