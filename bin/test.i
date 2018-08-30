# test input file
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=500
  ny=500
  meshtype=quad9
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
    value=1
    side=left
  [end]

  [right]
   type=dirichlet
   dof=c
   value=2
   side=right
  [end]
[]

[run]
  type=static
  proj=false
  debug=true
[]
