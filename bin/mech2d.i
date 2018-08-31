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
  params= 10.0e5 0.3
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
    type=dirichlet
    dof=ux
    value=0.1
    side=top
  [end]
[]

[run]
  type=static
  proj=true
  debug=true
[]
