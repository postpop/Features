# Features
For estimating LN-type models from spiking data.

## Model types include:
- standard STA and STC analysis
- STA and STC with a parametric information-theoretic objective 
- wrapper for MID
- wrapper for sparse GLMs
- ML estimation of linear and quadratic models with regularization

## Demo

```Matlab
load('demo.mat')% load sample spikeTimes and stim
resolution = 1;%ms
filterDuration = 64;%samples
sr = SpikeResp(spikeTimes, resolution);
sr.setStim(stim, resolution);
sf = FeaturesSTC(sr, filterDuration);
sf.getFeat();
plot(sf.STA);
```

## Code base used in:
Jan Clemens, Sandra Wohlgemuth, Bernhard Ronacher  
Nonlinear computations underlying temporal and population sparseness in the auditory system of the grasshopper  
[_2012_, Journal of Neuroscience, 32(29):10053-10062](http://www.jneurosci.org/content/32/29/10053.abstract) | [pdf](http://www.princeton.edu/~janc/pdf/clemens_2012_nonlinear.pdf)

Jan Clemens, Olaf Kutzki, Bernhard Ronacher, Susanne Schreiber, Sandra Wohlgemuth   
Efficient transformation of an auditory population code  
[_2011_, PNAS, 108:13812-13817](http://www.pnas.org/cgi/doi/10.1073/pnas.1104506108) | [pdf](http://www.princeton.edu/~janc/pdf/clemens_2011_efficient.pdf)


