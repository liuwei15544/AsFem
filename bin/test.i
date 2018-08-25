# test input file
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=10
  ny=10
  meshtype=quad4
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

  [right1]
   type=dirichlet
   dof=c
   value=5
   side=top
  [end]
[]
