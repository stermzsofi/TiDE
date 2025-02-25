"""

    Python wrapper to TiDE C++ code

    If use this package, please cite Kovács-Stermeczky & Vinkó, 2023PASP..135c4102K

    Calculate the monochromatic light curve of Tidal Disruption Events.

    Example usage: Plot a TDE light curve around an M=3e6 M_sun Blackhole, with
    standard options, from t= tmin days to t= 1000 days.

    ```python
    >>> import tidepy
    >>> p=tidepy.Parameters()
        
    >>> p.bh_M6=3
    >>> p.init()
    >>> p.param_init()
    >>> lc=tidepy.Light_curve_of_tde(p)
    >>> p.tend=1000            
    >>> res=lc.light_curve()
    >>> from matplotlib import pyplot as plt
    >>> plt.plot(res[0],res[1]*p.nu)
    [<matplotlib.lines.Line2D object at 0x7f2a85607d90>]
    >>> plt.ylim(bottom=1e36)
    (1e+36, 7.556604716127834e+41)
    >>> plt.yscale('log')
    >>> plt.show()
    ```
"""

from .tide import Parameters,Wind_part_of_lc,Disk_part_of_lc,Light_curve_of_tde