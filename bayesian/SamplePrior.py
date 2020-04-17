#!/usr/bin/env python
# coding: utf-8

import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
import pandas as pd
import pickle
from optparse import OptionParser
import os
import errno

parser = OptionParser(description="Samples the prior distribution of the parameters in the basic model with E | B uniform, and a single growth rate")
parser.add_option("-a", "--alpha",
                  type="float", dest="alpha", default=1,
                  help="Concentration parameter of tilted Dirichlet prior on w")
parser.add_option("-b", "--num_burnin_steps",
                  type="int", dest="num_burnin_steps", default=1,
                  help="Number of burn-in steps")
parser.add_option("-n", "--num_samples",
                  type="int", dest="num_samples", default=1,
                  help="Number of samples")
parser.add_option("-t", "--thinning",
                  type="int", dest="thinning", default=1,
                  help="Number of steps between samples")

(options, args) = parser.parse_args()

# Global params and data, as tensors
maxS = 30
gam = tfp.distributions.Gamma(6*1.5,1.5)
mean_dist = np.array([(gam.cdf(i)-gam.cdf(i-1)).numpy() for i in range(1,maxS+1)])
mean_dist = mean_dist/mean_dist.sum()
alpha = options.alpha
num_burnin_steps = options.num_burnin_steps
num_samples = options.num_samples
thinning = options.thinning

def logit(p):
    return tf.math.log(p/(1-p))

def inv_logit(x):
    return tf.math.exp(x)/(1+tf.math.exp(x))
            
def log_prior_fn(logitPi, logitLamw, logitLamv, logr ,g):
    alpha_r = 1
    beta_r = 1
    gamma = 1
    lp =  tf.math.log(inv_logit(logitPi)) + tf.math.log(1-inv_logit(logitPi))       # log Jacobian for logit
    lp += tf.math.log(inv_logit(logitLamw)) + tf.math.log(1-inv_logit(logitLamw))   # log Jacobian for logit
    lp += tf.math.log(inv_logit(logitLamv)) + tf.math.log(1-inv_logit(logitLamv))   # log Jacobian for logit 
    lp += alpha_r * logr - tf.math.exp(logr)/beta_r                                 # log-density of a log Gamma 
    lp += tf.reduce_sum(alpha * mean_dist * g - tf.math.exp(g))                     # log-density of a log Gamma  
    # Add regularisation for log concavity
    logw = tf.math.log(tf.math.exp(g)/tf.reduce_sum(tf.math.exp(g)))
    summands = tf.minimum(2*logw[1:-1] - logw[0:-2] - logw[2:], 0)
    lp += gamma * tf.reduce_sum(summands)
    return lp
    
rwmh = tfp.mcmc.RandomWalkMetropolis(
    log_prior_fn, new_state_fn = tfp.mcmc.random_walk_normal_fn(scale=0.1, name=None)
)

initial_state = [tf.zeros(1,dtype=tf.float32),
                 tf.zeros(1,dtype=tf.float32),
                 tf.zeros(1,dtype=tf.float32),
                 tf.zeros(1,dtype=tf.float32), 
                 tf.zeros(maxS,dtype=tf.float32)]

@tf.function
def run_chain_fn():
    return tfp.mcmc.sample_chain(
        num_results=num_samples,
        num_steps_between_results=thinning,
        num_burnin_steps=num_burnin_steps,
        current_state= initial_state,
        kernel=rwmh,
)

samples, _ = run_chain_fn()

samples = [x.numpy() for x in samples]

# Save samples and summary

r = np.exp(samples[3])
w = np.exp(samples[4])
w = np.transpose(np.transpose(w)/w.sum(axis=1))
w_means = (w*np.arange(maxS)).sum(1)
p14 = w[:,14:].sum(1)

output_nm = "prior_alpha_%.3f"%(alpha)
    
f = open("summary_"+output_nm+".txt",'w')
f.write("                Est     95% C.I.\n")
f.write("r               %.4f (%.4f,%.4f) \n"%(r.mean(),np.quantile(r,q=0.025),np.quantile(r,0.975)))
f.write("mean of w       %.4f (%.4f,%.4f) \n"%(w_means.mean(),np.quantile(w_means,q=0.025),np.quantile(w_means,q=0.975)))
f.write("p(S-T > 13)     %.4f (%.4f,%.4f) \n"%(p14.mean(),np.quantile(p14,q=0.025),np.quantile(p14,q=0.975)))
for i in range(30):
    f.write("w%i              %.4f (%.4f,%.4f) \n"%(i,w[:,i].mean(),np.quantile(w[:,i],q=0.025),np.quantile(w[:,i],q=0.975)))
f.close()

pickle.dump(samples, open("out_"+output_nm+".pkl",'wb'))
