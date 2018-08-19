# test input file
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=10
  ny=20
  meshtype=quad4
[]

[dofs]
name = ux   uy   uz
[]

[kernel]

  type=poisson

  mate=   user1

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
   type=neumann
   dof=u
   value=2
   side=right
  [end]
[]