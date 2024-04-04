# LAPS

LAPS (UC**LA**-**P**seudo-**S**pectral) is a pseudo-spectral 3D Hall-MHD code writen in FORTRAN, with corotating-expanding-box-model implemented. Brief notes on the numerical method are included as a PDF file here. A more detailed description can also be found in the paper: https://iopscience.iop.org/article/10.3847/1538-4357/ab5fce

To compile the code, MPI and FFTW are needed.

Samples of Python scripts to read the output are included.

<hr>
TODO: Instructions on how to set up initial conditions.

<hr>
Copyright 2024 Chen Shi

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

If you use the code to generate any publications, you should cite the following articles
@article{shi2020propagation,
  title={Propagation of Alfv{\'e}n waves in the expanding solar wind with the fast--slow stream interaction},
  author={Shi, Chen and Velli, Marco and Tenerani, Anna and Rappazzo, Franco and R{\'e}ville, Victor},
  journal={The Astrophysical Journal},
  volume={888},
  number={2},
  pages={68},
  year={2020},
  publisher={IOP Publishing}
}
