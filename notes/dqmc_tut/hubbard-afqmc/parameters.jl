# model parameters
Larr = [8]
ham_t = 1.0
ham_U = 4.0
trial_delta = 0.1

# QMC parameters
Theta = 10.0
beta = 0.0
dtau = 0.05
nstab = 10
Naf = 4
Nbin = 20
Nsweep = 200
Nana = Nbin
rerun = false

# BLAS / MKL threads per rank
nblas = 1

# Hubbard-Stratonovich decoupling channel
hs_channel = :charge

# Time-displaced Green's-function measurement
Mtau = false
Nmes = -1
