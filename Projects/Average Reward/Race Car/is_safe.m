function safe = is_safe(x, v)
    left_wall  = (x >= -2 & x <= -1) & (v >= 0 & v <= 2);
    right_wall = (x >= 1 & x <= 2) & (v >= 0 & v <= 2);
    bottom_bar = (x >= -2 & x <= 2) & (v >= -1 & v <= 0);
    safe = left_wall | right_wall | bottom_bar;
end