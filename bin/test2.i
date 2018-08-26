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
name = u v
[]

[kernel]
  type=uel1
  params= 1 2 4.0
[]

[bcs]
  [leftu]
    type=dirichlet
    dof=u
    value=1
    side=left
  [end]
  [rightu]
   type=dirichlet
   dof=u
   value=0
   side=right
  [end]
  [leftv]
    type=dirichlet
    dof=v
    value=2
    side=left
  [end]
  [rightv]
   type=dirichlet
   dof=v
   value=0
   side=right
  [end]
[]

[run]
  type=static
  proj=true
  output=10
[]


