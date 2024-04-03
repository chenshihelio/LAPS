# LAPS

LAPS (UC**LA**-**P**seudo-**S**pectral) is a pseudo-spectral 3D Hall-MHD code writen in FORTRAN, with corotating-expanding-box-model implemented. Brief notes on the numerical method are included as a PDF file here. A more detailed description can also be found in the paper: https://iopscience.iop.org/article/10.3847/1538-4357/ab5fce

To compile the code, MPI and FFTW are needed.

Samples of Python scripts to read the output are included.

The code is open-source and can be used freely. Please cite the following articles if you use the code to generate any publications:

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


<hr>
TODO: Instructions on how to set up initial conditions.

------
Copyright 2024 Chen Shi

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
