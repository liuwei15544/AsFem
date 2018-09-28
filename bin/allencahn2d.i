# the bechmark test for 2d allen-cahn equation
[mesh]
  type=asfem
  dim=2
  xmin=-2.0
  xmax= 2.0
  ymin=-2.0
  ymax= 2.0
  nx=250
  ny=250
  meshtype=quad9
[]

[dofs]
name = c
[]

[kernel]
  type=allencahn
  params= 0.0625
[]

[ics]
  [ringic]
    type=ringic
    dof=c
    params=0.4 1.0 0.0625
  [end]
[]


[run]
  type=transient
  proj=true
  debug=true
  dt=1.0e-4
  step=5000
  interval=100
[]
