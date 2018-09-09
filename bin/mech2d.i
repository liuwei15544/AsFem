# test input file
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=2.0
  nx=4
  ny=1
  meshtype=quad4
[]

[dofs]
name = ux uy
[]

[kernel]
  type=mechanics
  strain=small
  params= 10.0e5 0.3
  mate=linearelastic
[]

[bcs]
  [left_ux]
    type=dirichlet
    dof=ux
    value=0.0
    side=left
  [end]

  [left_uy]
   type=dirichlet
   dof=uy
   value=0.0
   side=left
  [end]
  
  [right_ux]
    type=dirichlet
    dof=ux
    value=0.1
    side=right
  [end]
[]

[run]
  type=static
  proj=true
  debug=true
  maxiter=2
[]
