# test input file for 2-dofs coupled example
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=80
  ny=80
  meshtype=quad9
[]

[dofs]
name = c mu
[]

[kernel]
  type=cahnhilliard
  params=1.0 2.0e-2
  mate=freeenergy
[]

[ics]
  [random]
    type=randomnoiseic
    dof=c
    params=0.1 0.02
  [end]
[]

[bcs]
  [flux]
    type=neumann
    dof=c
    value=0.01
    side=right
  [end]
[]

[run]
  type=transient
  proj=false
  debug=true
  dt=5.0e-6
  step=5000
  output=10
  nr=linesearch
[]


