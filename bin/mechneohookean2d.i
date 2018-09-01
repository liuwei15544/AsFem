# test input file
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=2.0
  nx=20
  ny=100
  meshtype=quad9
[]

[dofs]
name = ux uy
[]

[kernel]
  type=mechanics
  strain=finite
  params= 10.0e5 0.3
  mate=neohookean
[]

[bcs]
  [left_ux]
    type=dirichlet
    dof=ux
    value=0.0
    side=left
  [end]

  [bottom_uy]
   type=dirichlet
   dof=uy
   value=0.0
   side=bottom
  [end]
  
  [right_ux]
    type=tdirichlet
    dof=ux
    value=0.01
    side=right
  [end]
  [right_uy]
    type=tdirichlet
    dof=uy
    value=0.04
    side=right
  [end]
[]

[run]
  type=transient
  proj=true
  debug=true
  dt=1.0
  step=10
[]
