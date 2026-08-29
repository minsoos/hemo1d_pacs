import common as c
print( c.cfl_max_dt(c.element_size(128),2))
c.check_cfl(c.SPATIAL_DT, c.element_size(128), 2, 'x')