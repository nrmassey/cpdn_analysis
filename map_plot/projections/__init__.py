from matplotlib.projections import register_projection
from north_polar_proj import *
from south_polar_proj import *
from cylin_proj import *
from robin_proj import *
from winkel_iii_proj import *

register_projection(NorthPolarAxes)
register_projection(SouthPolarAxes)
register_projection(CylinAxes)
register_projection(RobinAxes)
register_projection(WinkelIIIAxes)
