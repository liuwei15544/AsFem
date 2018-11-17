# test input file
[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  zmin=0.0
  zmax=5.0
  nx=2
  ny=2
  nz=20
  meshtype=hex8
[]

[dofs]
name = c
[]

[kernel]
  type=poisson
  params= 1 2 4.0
[]

[bcs]
  [left]
    type=dirichlet
    dof=c
    value=-2.0
    side=back
  [end]

  [right]
   type=dirichlet
   dof=c
   value=2.0
   side=front
  [end]
[]

[run]
  type=static
  proj=false
  debug=true
[]
