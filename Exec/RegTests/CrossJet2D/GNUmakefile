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

# CLEM particle model (Source/Clem). Both builds live side by side:
#   make -j16                    -> PeleC2d.gnu.MPI.CLEM.ex   (CLEM)
#   make -j16 CLEM_MODEL=FALSE   -> PeleC2d.gnu.MPI.ex        (pure PeleC)
CLEM_MODEL ?= TRUE
ifeq ($(CLEM_MODEL), TRUE)
  USERSuffix = .CLEM
endif

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
Chemistry_Model := BurkeDryer
Transport_Model := Simple

# GNU Make
Bpack := ../CrossJet/Make.package
Blocs := ../CrossJet
PELE_HOME := ../../..
include $(PELE_HOME)/Exec/Make.PeleC
