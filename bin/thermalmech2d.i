# the bechmark test for 2d linear elastic problem
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=2.0
  nx=100
  ny=20
  meshtype=quad9
[]

[dofs]
name = ux uy c
[]

[kernel]
  type=thermalmechanics
  strain=small
  params= 1.0e5 0.3 1.0 0.08
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
    dof=c
    value=0.01
    side=bottom
  [end]
  [right_uy]
    type=dirichlet
    dof=c
    value=0.0
    side=top
  [end]
[]

[run]
  type=transient
  proj=true
  debug=true
  dt=1.0e-3
  step=50
  output=10
[]
