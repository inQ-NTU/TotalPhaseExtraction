function z_axis = get_symmetric_z_axis(no_gridpoints,grid_spacing)

grid_extend = (no_gridpoints-1)*grid_spacing;
z_axis = linspace(-grid_extend/2,grid_extend/2,no_gridpoints);

end