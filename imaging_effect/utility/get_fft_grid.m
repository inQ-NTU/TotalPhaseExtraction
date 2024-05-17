function grid = get_fft_grid(spacing,no_points) 

extend = no_points*spacing;

if mod(no_points,2)
    grid = ((1:no_points)-1/2)*spacing - extend/2;
else
    grid = ((1:no_points)-1)*spacing - extend/2;
end

end