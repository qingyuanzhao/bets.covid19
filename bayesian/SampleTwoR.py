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

parser = OptionParser(description="Samples the posterior distribution in the model with E | B uniform, and two different growth rates, before and after January 21")
parser.add_option("-l", "--location",
                  type="str", dest="location", default=None,
                  help="Include samples only from this location")
parser.add_option("-e", "--exclude",
                  type="float", dest="exclude", default=None,
                  help="Exclude samples from this location")
parser.add_option("-w", "--wuhan",
                  action="store_true", dest="wuhan",
                  help="Uses only people in Wuhan on day 0")
parser.add_option("-v", "--traveler",
                  action="store_true", dest="traveler",
                  help="Uses only people arriving in Wuhan after day 0")
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
parser.add_option("-i", "--initLogR",
                  type="float", dest="initLogR", default=0,
                  help="Initial value of log(r)")
parser.add_option("-c", "--copy",
                  type="int", dest="copy", default=0,
                  help="Index appended to simulation output")


(options, args) = parser.parse_args()


# Read data
wuhan_data = pd.read_csv("wuhan_exported.csv")

if (options.location!=None)+(options.exclude!=None)+(options.wuhan==True)+(options.traveler==True)>1:
    raise ValueError

# Subset samples
if options.location:
    wuhan_data = wuhan_data.loc[wuhan_data.Location==options.location,]
if options.exclude:
    wuhan_data = wuhan_data.loc[wuhan_data.Location!=options.exclude,]
if options.wuhan:
    wuhan_data = wuhan_data.loc[wuhan_data.B==0,]
if options.traveler:
    wuhan_data = wuhan_data.loc[wuhan_data.B>0,]

# Global params and data, as tensors
B = tf.constant(wuhan_data.B.values,dtype=tf.float32)
E = tf.constant(wuhan_data.E.values,dtype=tf.float32)
S = tf.constant(wuhan_data.S.values,dtype=tf.float32)
n = len(B)
L = 54
maxS = 30
alpha = options.alpha
gam = tfp.distributions.Gamma(6*1.5,1.5)
mean_dist = np.array([(gam.cdf(i)-gam.cdf(i-1)).numpy() for i in range(1,maxS+1)])
mean_dist = mean_dist/mean_dist.sum()
num_burnin_steps = options.num_burnin_steps
num_samples = options.num_samples
thinning = options.thinning

def logit(p):
    return tf.math.log(p/(1-p))

def inv_logit(x):
    return tf.math.exp(x)/(1+tf.math.exp(x))

def log_prior_fn(logitPi, logitLamw, logitLamv, logr1, r2 ,g):
    alpha_r = 1
    beta_r = 1
    alpha = 1
    gamma = 1
    lp =  tf.math.log(inv_logit(logitPi)) + tf.math.log(1-inv_logit(logitPi))       # log Jacobian for logit
    lp += tf.math.log(inv_logit(logitLamw)) + tf.math.log(1-inv_logit(logitLamw))   # log Jacobian for logit
    lp += tf.math.log(inv_logit(logitLamv)) + tf.math.log(1-inv_logit(logitLamv))   # log Jacobian for logit 
    lp += alpha_r * logr1 - tf.math.exp(logr1)/beta_r                                  # log-density of a log Gamma 
    lp += -r2**2/(2*4)                                                                 # normal prior 
    lp += tf.reduce_sum(alpha * mean_dist * g - tf.math.exp(g))                       # log-density of a log Gamma  
    # Add regularisation for log concavity
    logw = tf.math.log(tf.math.exp(g)/tf.reduce_sum(tf.math.exp(g)))
    summands = tf.minimum(2*logw[1:-1] - logw[0:-2] - logw[2:], 0)
    lp += gamma * tf.reduce_sum(summands)
    return lp
    
def log_posterior_fn(logitPi, logitLamw, logitLamv, logr1, r2, g):
    pi = inv_logit(logitPi)
    lamw = inv_logit(logitLamw)
    lamv = inv_logit(logitLamv)
    r1 = tf.math.exp(logr1)
    w = tf.math.exp(g)/tf.reduce_sum(tf.math.exp(g))
    w.padded = tf.concat([tf.zeros(1),w,tf.zeros(100)],axis=0)
    f = tf.zeros(n)
    for t in range(52):
        a = ( (1-pi)*lamw/(L+1)*tf.cast(tf.equal(B,0.0),dtype=tf.float32)+   \
                pi/L*lamv/L*tf.cast(tf.greater(B,0.0),dtype=tf.float32) )*   \
                tf.math.exp(r1*t)*tf.cast(tf.greater_equal(E, float(t)),dtype=tf.float32)*    \
                tf.cast(tf.greater_equal(float(t),B),dtype=tf.float32)* \
                tf.gather(w.padded,tf.maximum(tf.cast(S-t+1,dtype=tf.int32),0)) 
        f += a
    for t in range(52,L+1):
        a = ( (1-pi)*lamw/(L+1)*tf.cast(tf.equal(B,0.0),dtype=tf.float32)+    \
                pi/L*lamv/L*tf.cast(tf.greater(B,0.0),dtype=tf.float32) )*    \
                tf.math.exp(r1*51+r2*(t-51))*tf.cast(tf.greater_equal(E, float(t)),dtype=tf.float32)*  \
                tf.cast(tf.greater_equal(float(t),B),dtype=tf.float32)*   \
                tf.gather(w.padded,tf.maximum(tf.cast(S-t+1,dtype=tf.int32),0)) 
        f += a
    f = tf.reduce_sum(tf.math.log(f))
    # Sample Selection Conditioning
    p = 0
    for b in range(L+1):
        for t in range(b,52):
            p += ((1-pi)*lamw/(L+1)*(b==0) +   \
                    pi/L*lamv/L*(b>0)) *   \
                    tf.math.exp(r1*t)*(L+1-t)
        for t in range(max(52,b),L+1):
            p += ((1-pi)*lamw/(L+1)*(b==0) +  \
                    pi/L*lamv/L*(b>0)) *  \
                    tf.math.exp(r1*51+r2*(t-51))*(L+1-t)
    p = tf.math.log(p)
    log_post = f - n*p + log_prior_fn(logitPi,logitLamw,logitLamv,logr1, r2,g)
    return log_post


rwmh = tfp.mcmc.RandomWalkMetropolis(
    log_posterior_fn, new_state_fn = tfp.mcmc.random_walk_normal_fn(scale=0.1, name=None)
)

initial_state = [tf.zeros(1,dtype=tf.float32),
                 tf.zeros(1,dtype=tf.float32),
                 tf.zeros(1,dtype=tf.float32),
                 tf.zeros(1,dtype=tf.float32)*options.initLogR, 
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

r1 = np.exp(samples[3])
r2 = samples[4]
w = np.exp(samples[5])
w = np.transpose(np.transpose(w)/w.sum(axis=1))
w_means = (w*np.arange(maxS)).sum(1)
p14 = w[:,14:].sum(1)

output_nm = "twoR_alpha_%.3f_copy_%i"%(alpha,options.copy)
if options.location:
    output_nm = "%s_%s"%(output_nm,options.location)
elif options.exclude:
    output_nm = "%s_no_%s"%(output_nm,options.exclude)
elif options.wuhan:
    output_nm = "%s_from_wuhan"%(output_nm)
elif options.traveler:
    output_nm = "%s_traveler"%(output_nm)
else:
    output_nm = "%s_all"%(output_nm)

f = open("summary_"+output_nm+".txt",'w')
f.write("                Est     95% C.I.\n")
f.write("r1              %.4f (%.4f,%.4f) \n"%(r1.mean(),np.quantile(r1,q=0.025),np.quantile(r1,0.975)))
f.write("r2              %.4f (%.4f,%.4f) \n"%(r2.mean(),np.quantile(r2,q=0.025),np.quantile(r2,0.975)))
f.write("mean of w       %.4f (%.4f,%.4f) \n"%(w_means.mean(),np.quantile(w_means,q=0.025),np.quantile(w_means,q=0.975)))
f.write("p(S-T > 13)     %.4f (%.4f,%.4f) \n"%(p14.mean(),np.quantile(p14,q=0.025),np.quantile(p14,q=0.975)))
for i in range(30):
    f.write("w%i              %.4f (%.4f,%.4f) \n"%(i,w[:,i].mean(),np.quantile(w[:,i],q=0.025),np.quantile(w[:,i],q=0.975)))
f.close()

pickle.dump(samples, open("out_"+output_nm+".pkl",'wb'))


