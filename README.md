# TiDE

version 1.1.1

TiDE (TIdal Disruption Event) is a C++ code that computes the light curves or spectrum of tidal disruption events.

# Installation of TiDE

- Download the files from github
- Create a directory for the installed program (eg.: mkdir ∼/TiDE)
- In the download directory use autoreconf ( autoreconf --install) (if needed)
- Run the created configure file: (./configure --prefix ∼/TiDE)
- make
- make install

# Using TiDE

For detailed information to using TiDE, please read TiDE_documentation.pdf file.
Now it is possible to use TiDE through python. For python usage, install the python/tidepy package and set up TIDE_LIB
environment parameter to the installation directory/lib. Details are in the GettingStarted jupyter-notebook.

# Citing

If you are use TiDE through command line or through the TiDEpy package, please cite the original instrument paper of
[Kovács-Stermeczky & Vinkó (2023) ](https://ui.adsabs.harvard.edu/abs/2023PASP..135c4102K/abstract). The bibtex entry is
this:
```bibtex
@ARTICLE{2023PASP..135c4102K,
       author = {{Kov{\'a}cs-Stermeczky}, Zs{\'o}fia V. and {Vink{\'o}}, J{\'o}zsef},
        title = "{Comparison of Different Tidal Disruption Event Light Curve Models with TiDE, a New Modular Open Source Code}",
      journal = {\pasp},
     keywords = {Supermassive black holes, Accretion, Tidal disruption, 1663, 14, 1696, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2023,
        month = mar,
       volume = {135},
       number = {1045},
          eid = {034102},
        pages = {034102},
          doi = {10.1088/1538-3873/acb9bb},
archivePrefix = {arXiv},
       eprint = {2302.08441},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023PASP..135c4102K},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

