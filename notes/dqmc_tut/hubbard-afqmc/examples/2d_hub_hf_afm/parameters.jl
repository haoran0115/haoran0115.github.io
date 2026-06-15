const DTYPE = ComplexF64

# model parameters
Larr = [6, 6]
ham_t = 1.0
ham_U = 4.0
ham_mu = 0.0
trial_delta = 0.1

# Hubbard-Stratonovich decoupling channel (:spin or :SU2)
hs_channel = :SU2

# QMC parameters
Theta = 10.0
beta = 0.0
dtau = 0.05
nstab = 10
Naf = 4
Nbin = 20
Nsweep = 100
Nana = Nbin
rerun = false
Mtau = false

# BLAS / MKL threads per rank
nblas = 1
