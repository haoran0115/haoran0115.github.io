include("hubbard-afqmc.jl")
include("parameters.jl")

const hpar = ham_par(;
    Larr = Larr,
    ham_t = ham_t,
    ham_U = ham_U,
    trial_delta = trial_delta,
    hs_channel = hs_channel,
    Theta = Theta,
    beta = beta,
    dtau = dtau,
)

const qpar = qmc_par(;
    dim = length(Larr),
    Theta = Theta,
    beta = beta,
    dtau = dtau,
    nstab = nstab,
    Naf = Naf,
    Nbin = Nbin,
    Nsweep = Nsweep,
    Nana = Nana,
    rerun = rerun,
    nblas = nblas,
    Mtau = Mtau,
)

analysis(hpar, qpar)
