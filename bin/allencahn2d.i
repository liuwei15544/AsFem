# the bechmark test for 2d allen-cahn equation
[mesh]
  type=asfem
  dim=2
  xmin=-2.0
  xmax= 2.0
  ymin=-2.0
  ymax= 2.0
  nx=80
  ny=80
  meshtype=quad9
[]

[dofs]
name = c
[]

[kernel]
  type=allencahn
  params= 1.0
[]



[run]
  type=transient
  proj=true
  debug=true
  dt=1.0e-3
  step=150
  interval=10
[]
