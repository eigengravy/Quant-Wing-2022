#!/usr/bin/env python3
from scipy.stats import norm
from math import log, sqrt, pi, exp

d1 = lambda S, K, sigma, r, T, delta: (log(S/K)+(r-delta+sigma**2/2.)*T)/(sigma*sqrt(T))
d2 = lambda S, K, sigma, r, T, delta: d1(S, K, sigma, r, T, delta)-sigma*sqrt(T)

def call(S, K, sigma, r, T, delta):
    S*exp(-delta*T)*norm.cdf(d1(S, K, sigma, r, T, delta))-K*exp(-r*T)*norm.cdf(d2(S, K, sigma, r, T, delta))

def put(S, K, sigma, r, T, delta):
    K*exp(-r*T)*norm.cdf(-d2(S, K, sigma, r, T, delta)) - S*exp(-delta*T)*norm.cdf(-d1(S, K, sigma, r, T, delta))