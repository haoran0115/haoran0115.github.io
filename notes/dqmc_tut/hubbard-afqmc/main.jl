include("hubbard-afqmc.jl")
include("seeds.jl")

function run_qmc(hpar::ham_par, qpar::qmc_par, seeds::Vector{Int})
    MPI.Init()
    comm = MPI.COMM_WORLD
    irank = MPI.Comm_rank(comm)
    nrank = MPI.Comm_size(comm)

    if nrank > length(seeds)
        error("Need at least $nrank seeds")
    end

    BLAS.set_num_threads(qpar.nblas)

    latt = hubbard_latt(hpar)
    nsblock = div(qpar.Ntau, qpar.nstab)
    FTYPE = DTYPE

    Ttrial = zeros(FTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    Tmat = zeros(FTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    expT = zeros(FTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    expmT = zeros(FTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    Vmat = zeros(FTYPE, qpar.NVdim, qpar.NVdim, hpar.Nf, qpar.Nv)
    expV = zeros(FTYPE, qpar.NVdim, qpar.NVdim, hpar.Nf, qpar.Nv, qpar.Naf)
    expmV = zeros(FTYPE, qpar.NVdim, qpar.NVdim, hpar.Nf, qpar.Nv, qpar.Naf)
    expVidx = zeros(Int, qpar.NVdim, qpar.Nv, hpar.Ns)
    afconf = zeros(Int, qpar.Nv, hpar.Ns, qpar.Ntau)
    Pl = zeros(FTYPE, hpar.Npart, hpar.Ndim, hpar.Nf)
    Pr = zeros(FTYPE, hpar.Ndim, hpar.Npart, hpar.Nf)
    Bl = zeros(FTYPE, hpar.Npart, hpar.Ndim, hpar.Nf)
    Br = zeros(FTYPE, hpar.Ndim, hpar.Npart, hpar.Nf)
    Bl_stab = zeros(FTYPE, hpar.Npart, hpar.Ndim, hpar.Nf, nsblock + 1)
    Br_stab = zeros(FTYPE, hpar.Ndim, hpar.Npart, hpar.Nf, nsblock + 1)
    G = zeros(FTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    G_ = zeros(FTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)
    Δ = zeros(FTYPE, qpar.NVdim, qpar.NVdim, hpar.Nf)
    B12 = zeros(FTYPE, hpar.Ndim, hpar.Ndim, hpar.Nf)

    obeq = obs_eq(hpar, qpar)
    obtau = qpar.Mtau ? obs_tau(hpar, qpar) : nothing
    t_propose_bin = zeros(DTYPE, qpar.Nbin)
    t_mulTmat_bin = zeros(DTYPE, qpar.Nbin)
    t_measure_bin = zeros(DTYPE, qpar.Nbin)
    t_stablize_bin = zeros(DTYPE, qpar.Nbin)
    acc_loc_bin = zeros(DTYPE, qpar.Nbin)
    mean_sign_bin = zeros(DTYPE, qpar.Nbin)
    mean_sign_ke_bin = zeros(DTYPE, qpar.Nbin)
    mean_sign_pe_bin = zeros(DTYPE, qpar.Nbin)
    sign = one(FTYPE)

    Random.seed!(seeds[irank + 1])

    af_gam, af_eta = discrete_hs_parameters(qpar)
    build_hopping!(Tmat, Ttrial, hpar)
    degen = initialize_trial_state!(Pl, Pr, Ttrial, hpar)

    for spin = 1:hpar.Nf
        expT[:, :, spin] .= exp(-qpar.dtau .* Tmat[:, :, spin])
        expmT[:, :, spin] .= exp(+qpar.dtau .* Tmat[:, :, spin])
    end
    build_interaction!(Vmat, expV, expmV, expVidx, hpar, qpar, af_eta)
    afconf .= rand(1:qpar.Naf, size(afconf))

    if irank == 0
        @printf("Degen = %f\n", degen)
        mkpath("data/chk")
        obs_h5_init("data/data.h5", obeq, qpar.rerun, latt)
        if qpar.Mtau
            obs_h5_init("data/data_tau.h5", obtau, qpar.rerun, latt)
        end
    end
    MPI.Barrier(comm)

    G_err_max = 0.0
    G_err_mean = 0.0
    N_stab_checks = 0

    construct_Bl_Br!(Bl, Br, Bl_stab, Br_stab, qpar.nstab, Pl, Pr, expT, expV, expVidx, afconf, qpar.Ntau, hpar.Ns, hpar.Nf, qpar.Nv, qpar.NVdim)
    recalc_G_stable!(G, Bl, Br, Bl_stab, Br_stab, 0, qpar.nstab, Pl, Pr, expT, expV, expVidx, afconf, qpar.Ntau, hpar.Ns, hpar.Nf, qpar.Nv, qpar.NVdim)

    for ibin = 1:qpar.Nbin
        obs_reset!(obeq)
        if qpar.Mtau
            obs_reset!(obtau)
        end

        recalc_G_stable!(G, Bl, Br, Bl_stab, Br_stab, 0, qpar.nstab, Pl, Pr, expT, expV, expVidx, afconf, qpar.Ntau, hpar.Ns, hpar.Nf, qpar.Nv, qpar.NVdim)
        time = @elapsed for _ = 1:qpar.Nsweep
            for itau = 1:qpar.Ntau
                t_mulTmat_bin[ibin] += @elapsed @views @inbounds for spin = 1:hpar.Nf
                    mul!(B12[:, :, spin], expT[:, :, spin], G[:, :, spin])
                    mul!(G[:, :, spin], B12[:, :, spin], expmT[:, :, spin])
                end

                t_propose_bin[ibin] += @elapsed for is = 1:hpar.Ns
                    for iv = 1:qpar.Nv
                        af_curr = afconf[iv, is, itau]
                        fidx = expVidx[1, iv, is]

                        @views @inbounds for spin = 1:hpar.Nf
                            v = expV[1, 1, spin, iv, af_curr]
                            mv = expmV[1, 1, spin, iv, af_curr]
                            prop_left_1x1!(G, v, fidx, spin)
                            prop_right_1x1!(G, mv, fidx, spin)
                        end

                        af_old = afconf[iv, is, itau]
                        af_new = mod1(af_old + rand(1:qpar.Naf - 1), qpar.Naf)
                        r, fidx = propose_r_1x1(af_new, afconf, af_gam, itau, is, iv, G, expV, expmV, expVidx, hpar.Nf, hpar.NSUN, Δ)

                        if abs(r) > rand()
                            acc_loc_bin[ibin] += 1.0
                            sign *= r / abs(r)
                            afconf[iv, is, itau] = af_new
                            update_G_1x1!(G, Δ, fidx, hpar.Nf)
                        end
                    end
                end

                if itau == div(qpar.Ntau, 2)
                    t_measure_bin[ibin] += @elapsed begin
                        measure_obs_eq!(obeq, hpar, qpar, latt, Tmat, sign, G)
                        if qpar.Mtau && itau == qpar.Nmes0
                            measure_obs_tau!(obtau, hpar, qpar, latt, sign, G, afconf, Bl_stab, Br_stab, expT, expmT, expV, expmV, expVidx)
                        end
                    end
                end
                if qpar.Mtau && itau == qpar.Nmes0 && itau != div(qpar.Ntau, 2)
                    t_measure_bin[ibin] += @elapsed begin
                        measure_obs_tau!(obtau, hpar, qpar, latt, sign, G, afconf, Bl_stab, Br_stab, expT, expmT, expV, expmV, expVidx)
                    end
                end

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

            for itau = qpar.Ntau:-1:1
                if itau == div(qpar.Ntau, 2)
                    t_measure_bin[ibin] += @elapsed begin
                        measure_obs_eq!(obeq, hpar, qpar, latt, Tmat, sign, G)
                    end
                end

                t_propose_bin[ibin] += @elapsed for is = hpar.Ns:-1:1
                    for iv = qpar.Nv:-1:1
                        af_old = afconf[iv, is, itau]
                        af_new = mod1(af_old + rand(1:qpar.Naf - 1), qpar.Naf)
                        r, fidx = propose_r_1x1(af_new, afconf, af_gam, itau, is, iv, G, expV, expmV, expVidx, hpar.Nf, hpar.NSUN, Δ)

                        if abs(r) > rand()
                            acc_loc_bin[ibin] += 1.0
                            sign *= r / abs(r)
                            afconf[iv, is, itau] = af_new
                            update_G_1x1!(G, Δ, fidx, hpar.Nf)
                        end

                        af_curr = afconf[iv, is, itau]
                        fidx = expVidx[1, iv, is]
                        @views @inbounds for spin = 1:hpar.Nf
                            mv = expmV[1, 1, spin, iv, af_curr]
                            v = expV[1, 1, spin, iv, af_curr]
                            prop_left_1x1!(G, mv, fidx, spin)
                            prop_right_1x1!(G, v, fidx, spin)
                        end
                    end
                end

                t_mulTmat_bin[ibin] += @elapsed @views @inbounds for spin = 1:hpar.Nf
                    mul!(B12[:, :, spin], expmT[:, :, spin], G[:, :, spin])
                    mul!(G[:, :, spin], B12[:, :, spin], expT[:, :, spin])
                end

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

        obs_avg!(obeq)
        obeq_avg = obs_mpi_avg(obeq, hpar, qpar)
        obtau_avg = qpar.Mtau ? obs_mpi_avg(obtau, hpar, qpar) : nothing

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
    hs_channel = hs_channel,
    Mtau = Mtau,
)

if abspath(PROGRAM_FILE) == @__FILE__
    run_qmc(hpar, qpar, seeds)
end
