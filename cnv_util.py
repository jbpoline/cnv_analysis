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
