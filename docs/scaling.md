---
title: ""
author: "Masatoshi Katabuchi"
date: "`r format(Sys.time(), '%B %d, %Y')`"
fontsize: 12pt
format:
  html:
    theme: coderpro
    toc: true
    toc-depth: 2
    number-sections: true
    smooth-scroll: true
    standalone: true
    embed-resources: true
---



# Weibull

$$
y = a \left\{1 - \exp(-b D^{k})\right\}
$$

Log-scale

$$
\log y = \log a + \log(1 - \exp(-b D^{k}))
$$

Log-scale with scaled D

$$
\log y = \log a\prime + \log(1 - \exp(-b\prime D^{k\prime}/D_{max}^{k\prime} ))
$$

Relationships between the parameters

$$
\log a = \log a\prime
$$

$$
k = k\prime
$$

$$
bD^k = -b\prime D^{k\prime}/D_{max}^{k\prime}
$$

$$
b = b\prime D_{max}^{-k}
$$

# gMM (generalized Michaelis–Menten)

$$
y = \frac{aD^k}{k + D^b}
$$

Log-scale

$$
\log y = \log a + b\log D - \log (k + D^b)
$$


Log-scale with scaled D

$$
\log y = \log a\prime + b\prime(\log D - \log D_{max}) - \log (k\prime + D^{b\prime}/D_{max}^{b\prime})
$$

Relationships between the parameters

$$
\log a = \log a\prime - b\prime \log D_{max}
$$

$$
b = b\prime
$$

$$
k + D^b = k\prime + D^b/D_{max}^b
$$

$$
k = k\prime + D^b (D_{max}^b - 1)
$$
