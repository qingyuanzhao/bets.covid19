## MCMC for sampling the posterior distribution in the BETS model for Wuhan-exported cases of COVID-19

These scripts implement random walk Metropolis-Hastings sampling for the prior and posterior distribution in the models described in Section 5 of:

Qingyuan Zhao, Nianqiao Ju, and Sergio Bacallado. "BETS: The dangers
of selection bias in early analyses of the coronavirus
epidemic", 2020. arXiv:2004.07743.

For script usage see, for example:

```
$ python SamplePrior.py -h
```

### Requirements

The code was tested on:

* Python 3.5
* Tensorflow 2.1.0
* Tensorflow-probability 0.9.0
