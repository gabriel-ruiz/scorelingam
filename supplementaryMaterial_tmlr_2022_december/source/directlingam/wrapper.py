#!/usr/bin/python3
import numpy as np
import pandas as pd
import graphviz
import lingam
# from lingam import utils
from lingam.utils import make_dot

def lingam_fit(X,kernel=False):
  # model = DirectLiNGAM()
  model = lingam.DirectLiNGAM()
  #
  if kernel == True:
    model = lingam.DirectLiNGAM(measure='kernel')
  #
  order_est = model.fit(X)
  return order_est

# with prior knowledge of causal paths (see paper)
def lingam_fit_prior(X,prior_knowledge):
  model = lingam.DirectLiNGAM(prior_knowledge=prior_knowledge)
  order_est = model.fit(X)
  return order_est
