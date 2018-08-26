# test input file for 2-dofs coupled example
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=20
  ny=20
  meshtype=quad9
[]

[dofs]
name = c mu
[]

[kernel]
  type=cahnhilliard
  params=1.0 2.0e-2
[]

[ics]
  [random]
    type=randomnoise
    dof=c
    params=0.63 0.02
  [end]
[]

[run]
  type=transient
  proj=false
  dt=5.0e-6
  step=10
  output=10
[]


