# test input file
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=2.0
  nx=100
  ny=20
  meshtype=quad4
[]

[dofs]
name = disp_x disp_y
[]

[kernel]
  type=mechanics
  strain=finite
  params= 10.0e5 0.3
  mate=deformgrad
[]

[bcs]
  [left_ux]
    type=dirichlet
    dof=disp_x
    value=0.0
    side=left
  [end]

  [bottom_uy]
   type=dirichlet
   dof=disp_y
   value=0.0
   side=left
  [end]
  
  [right_ux]
    type=dirichlet
    dof=disp_x
    value=-0.1
    side=right
  [end]
  [right_uy]
    type=dirichlet
    dof=disp_y
    value=-0.1
    side=right
  [end]
[]

[run]
  type=static
  proj=true
  debug=true
[]
