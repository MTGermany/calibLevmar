# calibTraj 

The repo is now complete and working but not yet documented.
Now just stubs.

## Compiling

Assume that you are on a linux system and gnu g++ is
implemented, just enter 

```
make calibTraj
```
When using another compiler, change the makefile
accordingly.

## Running a calibration project


## Files of a calibration poroject


### Top-level input files

The following two files are the only ones that must be present in
every `calibTraj` project.

* `.run`



### Model specification

This is the *Unique selling point* of the FC-model calibration  `calibTraj`: The
generality in specifying several different models in one and the same 
simulation. Even mixing of discrete CA-like and time-continuous models
is allowed.

* `.<MOD><n>` Parameter specification file of the chosen car-following
  model, e.g., `.IDM1`, `.OVM3`, or `.BIDM1` (Brownian IDM).
  `<MOD>` is a model abbreviation consiting of two to four
  uppercase letters associated to the model numbers as follows:
  0=IDM, 1=VW, 2=HDM, 3=FVDM, 4=FPE, 5=HSM, 6=VDT, 7=FVDM, 9=CDDA,
  10=GIP, 11=KCA, 12=PT, 13=ASGM, 14=ADAS, 16=CACC, 17=PCF, 18=LCM,
  19=BIDM.  

 
  *Notice*: Some models come in different variants such as the
  OVM/FVDM and the variant is selected by a `choice_variant` entry in
  the parameter file. For example, in the `.FVDM` files,
  `choice_variant=0` denotes the original OV function based on tanh
  functions while `choice_variant=1` denotes an OV function making up
  a triangular fundamental diagram. 


### Output control and output files



## Numerical Integration inside each evaluation of the objective function

By default, `calibTraj` uses the _ballistic scheme_ with following pseudo-code:

```
speed(t+dt)=speed(t)+acc(t)*dt;
pos(t+dt)=pos(t)+speed(t)*dt+1/2*acc(t)*dt^2;
```

where `acc(t)` is the acceleration calculated by the car-following model
at the (old) time t.

Although technically of first consistency order, it turned out to be
the best overall choice in terms of robustness and precision for a
given computation load. Specifically, it is significantly more
performant than simple Euler (where the `1/2*acc(t)*dt^2` is omitted
in the positional update). 

## References

[1] M. Treiber and A. Kesting. [_Traffic Flow Dynamics, Data, Models and Simulation_](http://www.traffic-flow-dynamics.org). Springer 2013. [Link](http://www.springer.com/physics/complexity/book/978-3-642-32459-8)
    
[2] A. Kesting, M. Treiber, and D. Helbing. _Enhanced intelligent driver model to access the impact of driving strategies on traffic capacity_. Philosophical Transactions of the Royal Society A, 4585-4605 (2010). [Preprint](http://arxiv.org/abs/0912.3613)
    
[3] Calibration refs and levmar
