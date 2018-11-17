# test input file for cahnhilliard equation in 3D case
[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  zmin=0.0
  zmax=10.0
  nx=4
  ny=4
  nz=80
  meshtype=hex20
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
    params=0.63 0.02
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


