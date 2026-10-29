# AMReX
DIM = 2
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Fernando-Clem: enable the CLEM subgrid model (forces AMReX particles)
CLEM_MODEL = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
PELE_CVODE_FORCE_YCORDER = FALSE
PELE_USE_MAGMA = FALSE
PELE_COMPILE_AJACOBIAN = TRUE
Eos_Model := Fuego
Chemistry_Model := LiDryer
Transport_Model := Simple

# GNU Make
Bpack := ./Make.package
Blocs := .
PELE_HOME := /home/fernando/FPCE/ParticleCLEM/PeleC
include $(PELE_HOME)/Exec/Make.PeleC
