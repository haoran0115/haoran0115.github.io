include("hubbard-afqmc.jl")
include("seeds.jl")

function run_qmc(hpar::ham_par, qpar::qmc_par, seeds::Vector{Int})

    # MPI initialization
    # each MPI process corresponds to a Markov chain
    MPI.Init()
    comm = MPI.COMM_WORLD
    irank = MPI.Comm_rank(comm)
    nrank = MPI.Comm_size(comm)

    if nrank > length(seeds)
        error("Need at least $nrank seeds")
    end

    # multi-threading within a single Markov chain
    # useful for accelerating matrix multiplications
    BLAS.set_num_threads(qpar.nblas)

    # lattice object, useful for lattice Fourier transformations
    latt = hubbard_latt(hpar)

    # nsblock = Ntau / nstab, number of stabliation blocs
    nsblock = div(qpar.Ntau, qpar.nstab)

    # kinetic matrices
    # Tmat = T_{ij}
    Tmat = zeros(DTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    # Ttrial is the one-body matrix for 
    Ttrial = zeros(DTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    # e^{+dtau * Tmat}
    expT = zeros(DTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    # e^{-dtau * Tmat}
    expmT = zeros(DTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)

    # HS-decoupled interaction matrices
    # Vmat = Vi
    Vmat = zeros(DTYPE, qpar.NVdim, qpar.NVdim, hpar.Nf, qpar.Nv)
    # expV  = exp(+sqrt(-λi * Δτ) * Vi)
    expV = zeros(DTYPE, qpar.NVdim, qpar.NVdim, hpar.Nf, qpar.Nv, qpar.Naf)
    # expmV = exp(-sqrt(-λi * Δτ) * Vi)
    expmV = zeros(DTYPE, qpar.NVdim, qpar.NVdim, hpar.Nf, qpar.Nv, qpar.Naf)
    # expVidx 
    expVidx = zeros(Int, qpar.NVdim, qpar.Nv, hpar.Ns)
    # auxiliary field configuration C
    afconf = zeros(Int, qpar.Nv, hpar.Ns, qpar.Ntau)

    # left and right trial wavefunction
    Pl = zeros(DTYPE, hpar.Npart, hpar.Ndim, hpar.Nf)
    Pr = zeros(DTYPE, hpar.Ndim, hpar.Npart, hpar.Nf)

    # left and right propagators
    Bl = zeros(DTYPE, hpar.Npart, hpar.Ndim, hpar.Nf)
    Br = zeros(DTYPE, hpar.Ndim, hpar.Npart, hpar.Nf)

    # what's this? check
    B12 = zeros(DTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)

    # stablized left and right propagators
    Bl_stab = zeros(DTYPE, hpar.Npart, hpar.Ndim, hpar.Nf, nsblock + 1)
    Br_stab = zeros(DTYPE, hpar.Ndim, hpar.Npart, hpar.Nf, nsblock + 1)
    
    # G_{ij} = ⟨c_i c_j^†⟩ is the time-equal Green's function
    G = zeros(DTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    G_ = zeros(DTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    Δ = zeros(DTYPE, qpar.NVdim, qpar.NVdim, hpar.Nf)

    # time-equal observables
    obeq = obs_eq(hpar, qpar)
    # time-displaced observables
    obtau = qpar.Mtau ? obs_tau(hpar, qpar) : nothing

    # execution time statistics
    t_propose_bin = zeros(DTYPE, qpar.Nbin)
    t_mulTmat_bin = zeros(DTYPE, qpar.Nbin)
    t_measure_bin = zeros(DTYPE, qpar.Nbin)
    t_stablize_bin = zeros(DTYPE, qpar.Nbin)
    
    # acceptance rate
    acc_loc_bin = zeros(DTYPE, qpar.Nbin)

    # scalar observables
    mean_sign_bin = zeros(DTYPE, qpar.Nbin)
    mean_sign_ke_bin = zeros(DTYPE, qpar.Nbin)
    mean_sign_pe_bin = zeros(DTYPE, qpar.Nbin)

    # initialize fermion sign
    sign = one(DTYPE)

    # initialize random seed
    Random.seed!(seeds[irank + 1])

    # initialize HS field
    # af_gam: γ, af_eta: η
    af_gam, af_eta = discrete_hs_parameters(qpar)

    # build hopping matrices
    build_hopping!(Tmat, Ttrial, hpar)

    # check degeneracies of trial wavefunction
    degen = initialize_trial_state!(Pl, Pr, Ttrial, hpar)

    # compute e^(-dtau * T) and e^(+dtau * T)
    for ifl = 1:hpar.Nf
        expT[:, :, ifl] .= exp(-qpar.dtau .* Tmat[:, :, ifl])
        expmT[:, :, ifl] .= exp(+qpar.dtau .* Tmat[:, :, ifl])
    end

    # build interaction matrices
    build_interaction!(Vmat, expV, expmV, expVidx, hpar, qpar, af_eta)
    afconf .= rand(1:qpar.Naf, size(afconf))

    # irank==0: only let master process print the information
    if irank == 0
        @printf("Trial WF degen = %f\n", degen)
        mkpath("data/chk")
        obs_h5_init("data/data.h5", obeq, qpar.rerun, latt)
        if qpar.Mtau
            obs_h5_init("data/data_tau.h5", obtau, qpar.rerun, latt)
        end
    end
    MPI.Barrier(comm)

    # stablization error check
    G_err_max = 0.0
    G_err_mean = 0.0
    N_stab_checks = 0

    # construct stablized projected operators
    construct_Bl_Br!(Bl, Br, Bl_stab, Br_stab, qpar.nstab, Pl, Pr, expT, expV, expVidx, afconf, qpar.Ntau, hpar.Ns, hpar.Nf, qpar.Nv, qpar.NVdim)
    # calculate equal-time Green's function at itau = 0
    recalc_G_stable!(G, Bl, Br, Bl_stab, Br_stab, 0, qpar.nstab, Pl, Pr, expT, expV, expVidx, afconf, qpar.Ntau, hpar.Ns, hpar.Nf, qpar.Nv, qpar.NVdim)

    # for each bin
    for ibin = 1:qpar.Nbin

        # reset observables to zero
        obs_reset!(obeq)
        if qpar.Mtau
            obs_reset!(obtau)
        end

        # re-calculate equal-time Green's function at itau = 0
        recalc_G_stable!(G, Bl, Br, Bl_stab, Br_stab, 0, qpar.nstab, Pl, Pr, expT, expV, expVidx, afconf, qpar.Ntau, hpar.Ns, hpar.Nf, qpar.Nv, qpar.NVdim)

        # for each sweep
        time = @elapsed for _ = 1:qpar.Nsweep

            # forward sweep
            # for each imaginary time slice
            for itau = 1:qpar.Ntau
                # propagate Green's function, kinetic part
                t_mulTmat_bin[ibin] += @elapsed @views @inbounds for ifl = 1:hpar.Nf
                    mul!(B12[:, :, ifl], expT[:, :, ifl], G[:, :, ifl])
                    mul!(G[:, :, ifl], B12[:, :, ifl], expmT[:, :, ifl])
                end

                # for each spatial index
                t_propose_bin[ibin] += @elapsed for is = 1:hpar.Ns
                    for iv = 1:qpar.Nv
                        # propagate Green's function
                        af_curr = afconf[iv, is, itau]
                        fidx = expVidx[1, iv, is]
                        @views @inbounds for ifl = 1:hpar.Nf
                            v = expV[1, 1, ifl, iv, af_curr]
                            mv = expmV[1, 1, ifl, iv, af_curr]
                            prop_left_1x1!(G, v, fidx, ifl)
                            prop_right_1x1!(G, mv, fidx, ifl)
                        end

                        # propose new auxiliary field
                        af_old = afconf[iv, is, itau]
                        af_new = mod1(af_old + rand(1:qpar.Naf - 1), qpar.Naf)
                        r, fidx = propose_r_1x1(af_new, afconf, af_gam, af_eta, itau, is, iv, G, expV, expmV, Vmat, expVidx, hpar.Nf, hpar.NSUN, hpar.hs_channel, Δ)
                        r_real = real(r)
                        r_weight = abs(r_real)
                        # Metropolis-Hastings accept/reject
                        # if accept, update sign, field configs, Green's function
                        # the code is implemented in the general case, where r might be not non-negative
                        if r_weight > 0 && r_weight > rand()
                            acc_loc_bin[ibin] += 1.0
                            sign *= r_real / r_weight
                            afconf[iv, is, itau] = af_new
                            update_G_1x1!(G, Δ, fidx, hpar.Nf)
                        end
                    end
                end

                # measure time-equal Green's function
                # ⟨e^{-(Θ+β/2)H} c_i c_j^† e^{-(Θ+β/2)H}⟩ / ⟨e^{-(2Θ+β)H}⟩
                if itau == div(qpar.Ntau, 2)
                    t_measure_bin[ibin] += @elapsed begin
                        measure_obs_eq!(obeq, hpar, qpar, latt, Tmat, sign, G)
                    end
                end

                # stablization
                if mod(itau, qpar.nstab) == 0
                    t_stablize_bin[ibin] += @elapsed begin
                        copyto!(G_, G)
                        stable_G!(G, 1, Bl_stab, Br_stab, Pl, Pr, itau, qpar.nstab, expT, expV, expVidx, afconf, qpar.Ntau, hpar.Ns, hpar.Nf, qpar.Nv, qpar.NVdim)
                        dG = abs.(G .- G_)
                        max_err = maximum(dG)
                        G_err_max = max(G_err_max, max_err)
                        G_err_mean += sum(dG) / length(dG)
                        N_stab_checks += 1
                        if max_err > G_STAB_ERR
                            error("Loss of precision: Green's function error ($max_err) exceeds threshold")
                        end
                    end
                end
            end

            # downward sweep, almost same as the upward sweep
            for itau = qpar.Ntau:-1:1
                # measure time-equal observables
                # e.g. time-equal Green's function ⟨e^{-(Θ+β/2)H} c_i c_j^† e^{-(Θ+β/2)H}⟩ / ⟨e^{-(2Θ+β)H}⟩
                if itau == div(qpar.Ntau, 2)
                    t_measure_bin[ibin] += @elapsed begin
                        measure_obs_eq!(obeq, hpar, qpar, latt, Tmat, sign, G)
                    end
                end

                # measure time-dispalced Green's function
                # e.g. time-dispalced Green's function ⟨e^{-(Θ+β/2)H} c_i(τ) c_j^† e^{-(Θ+β/2)H}⟩ / ⟨e^{-(2Θ+β)H}⟩
                if qpar.Mtau && itau == hpar.Nmes0
                    t_measure_bin[ibin] += @elapsed begin
                        measure_obs_tau!(obtau, hpar, qpar, latt, sign, G, afconf, Bl_stab, Br_stab, expT, expmT, expV, expmV, expVidx)
                    end
                end

                t_propose_bin[ibin] += @elapsed for is = hpar.Ns:-1:1
                    for iv = qpar.Nv:-1:1
                        # propose new auxiliary field
                        af_old = afconf[iv, is, itau]
                        af_new = mod1(af_old + rand(1:qpar.Naf - 1), qpar.Naf)
                        r, fidx = propose_r_1x1(af_new, afconf, af_gam, af_eta, itau, is, iv, G, expV, expmV, Vmat, expVidx, hpar.Nf, hpar.NSUN, hpar.hs_channel, Δ)
                        r_real = real(r)
                        r_weight = abs(r_real)
                        # Metropolis-Hastings accept/reject
                        # if accept, update sign, field configs, Green's function
                        if r_weight > 0 && r_weight > rand()
                            acc_loc_bin[ibin] += 1.0
                            sign *= r_real / r_weight
                            afconf[iv, is, itau] = af_new
                            update_G_1x1!(G, Δ, fidx, hpar.Nf)
                        end

                        # propagate Green's function
                        af_curr = afconf[iv, is, itau]
                        fidx = expVidx[1, iv, is]
                        @views @inbounds for ifl = 1:hpar.Nf
                            mv = expmV[1, 1, ifl, iv, af_curr]
                            v = expV[1, 1, ifl, iv, af_curr]
                            prop_left_1x1!(G, mv, fidx, ifl)
                            prop_right_1x1!(G, v, fidx, ifl)
                        end
                    end
                end

                # propagate Green's function, kinetic part
                t_mulTmat_bin[ibin] += @elapsed @views @inbounds for ifl = 1:hpar.Nf
                    mul!(B12[:, :, ifl], expmT[:, :, ifl], G[:, :, ifl])
                    mul!(G[:, :, ifl], B12[:, :, ifl], expT[:, :, ifl])
                end

                # stablization
                if mod(itau - 1, qpar.nstab) == 0
                    t_stablize_bin[ibin] += @elapsed begin
                        copyto!(G_, G)
                        stable_G!(G, 0, Bl_stab, Br_stab, Pl, Pr, itau - 1, qpar.nstab, expT, expV, expVidx, afconf, qpar.Ntau, hpar.Ns, hpar.Nf, qpar.Nv, qpar.NVdim)
                        dG = abs.(G .- G_)
                        max_err = maximum(dG)
                        G_err_max = max(G_err_max, max_err)
                        G_err_mean += sum(dG) / length(dG)
                        N_stab_checks += 1
                        if max_err > G_STAB_ERR
                            error("Loss of precision: Green's function error ($max_err) exceeds threshold")
                        end
                    end
                end
            end
        end

        # average observables over a single bin
        obs_avg!(obeq)
        if qpar.Mtau
            obs_avg!(obtau)
        end

        # average observables over all Markov chains
        obeq_avg = obs_mpi_avg(obeq, hpar, qpar)
        obtau_avg = qpar.Mtau ? obs_mpi_avg(obtau, hpar, qpar) : nothing

        # print information and save data
        if irank == 0
            t_io = @elapsed obs_h5_append("data/data.h5", obeq_avg, latt)
            if qpar.Mtau
                t_io += @elapsed obs_h5_append("data/data_tau.h5", obtau_avg, latt)
            end
            mean_sign = obeq_avg.sign[1]
            mean_ke = abs(obeq_avg.sign[1]) > 0 ? obeq_avg.KE[1] / obeq_avg.sign[1] : zero(DTYPE)
            mean_pe = abs(obeq_avg.sign[1]) > 0 ? obeq_avg.PE[1] / obeq_avg.sign[1] : zero(DTYPE)
            mean_sign_bin[ibin] = mean_sign
            mean_sign_ke_bin[ibin] = obeq_avg.KE[1]
            mean_sign_pe_bin[ibin] = obeq_avg.PE[1]
            @printf(
                "Bin = %3d, Time = %2.3fs, sign = %s, KE = %s, PE = %s\n",
                ibin,
                time,
                fmtc(mean_sign),
                fmtc(mean_ke),
                fmtc(mean_pe),
            )
            @printf(
                "propose_r = %7.3fs, mulT = %7.3fs, measure = %7.3fs, stablize = %7.3fs, io = %7.3fs\n",
                t_propose_bin[ibin],
                t_mulTmat_bin[ibin],
                t_measure_bin[ibin],
                t_stablize_bin[ibin],
                t_io,
            )
            @printf("acc_loc = %.6f\n", acc_loc_bin[ibin] / (2 * qpar.Nsweep * hpar.Ns * qpar.Ntau))
            @printf("G_err: max = %e, mean = %e\n\n", G_err_max, N_stab_checks > 0 ? G_err_mean / N_stab_checks : 0.0)

            if mod(ibin, qpar.Nana) == 0
                analysis(hpar, qpar)
            end
        end

        MPI.Barrier(comm)
        G_err_max = 0.0
        G_err_mean = 0.0
        N_stab_checks = 0

        rand_state = copy(Random.default_rng())
        @save "data/chk/rand_state_$irank.jld2" rand_state
        @save "data/chk/afconf_$irank.jld2" afconf
    end

    # print some scalar measurement at the end of simulation
    if irank == 0
        sign_jk = jackknife_mean(mean_sign_bin)
        ke_jk = jackknife_ratio(mean_sign_ke_bin, mean_sign_bin)
        pe_jk = jackknife_ratio(mean_sign_pe_bin, mean_sign_bin)
        @printf("⟨sign⟩ = %s ± %s\n", fmtc(sign_jk[1]), fmtc(sign_jk[2]))
        @printf("⟨KE⟩   = %s ± %s\n", fmtc(ke_jk[1]), fmtc(ke_jk[2]))
        @printf("⟨PE⟩   = %s ± %s\n", fmtc(pe_jk[1]), fmtc(pe_jk[2]))
    end

    MPI.Barrier(comm)
    MPI.Finalize()
end

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

if abspath(PROGRAM_FILE) == @__FILE__
    run_qmc(hpar, qpar, seeds)
end
