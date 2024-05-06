# LAPS

LAPS (UC**LA**-**P**seudo-**S**pectral) is a 3D MPI-parallelized Fourier-transform-based pseudo-spectral Hall-MHD code writen in FORTRAN, with corotating-expanding-box-model implemented. 

Brief notes on the numerical method are included as a PDF file: /src_compressible/notes_MHD3D_EBM.pdf. Some more detailed description can be found in the paper: https://iopscience.iop.org/article/10.3847/1538-4357/ab5fce. A thorough description of all the numerical methods along with fundamental tests will be published as a standalone paper soon.

## Installation & Compilation

To compile the code, MPI and FFTW are necessary. 

### OpenMPI
You can download the latest OpenMPI from https://www.open-mpi.org/software/ompi/v5.0/. The instruction to install OpenMPI can be found here: https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/index.html.

### FFTW
You can download FFTW from https://www.fftw.org/download.html and install it follow the instruction here https://www.fftw.org/fftw3_doc/Installation-and-Customization.html to install.

### Compile LAPS
After installing MPI and FFTW, LAPS is ready to be compiled. From terminal, go to the directory of the desired version of the code, i.e. /src_compressible for 3D compressible version, /src_compressible/2D for 2D compressible version, /src_incompressible for 3D incompressible version, and /src_incompressible/2D for 2D incompressible version. Then run 
> make
The only thing that you may need to modify in "**makefile**" is "**fftwpath**", which by default is "/usr/local". You need to change it to the location your FFTW is installed. In addition, depending on the specific FORTRAN compiler you use, you may need to change the "**-r8**" flag in "**OPTIONS**" to "**-fdefault-real-8**". If successfully compiled, a "**mhd.exe**" will be generated.


## Run LAPS

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
