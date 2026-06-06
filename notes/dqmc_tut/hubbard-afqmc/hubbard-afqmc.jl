using LinearAlgebra
using Statistics
using Random
using Base.Threads
using Printf
using JLD2
using HDF5
using MPI

const DTYPE = Float64
const G_STAB_ERR = 1e-3

function fmtc(z)
    if abs(imag(z)) < 1e-10
        return @sprintf("%+10.6f", real(z))
    end
    return @sprintf("%+10.6f%+10.6fim", real(z), imag(z))
end

struct ham_par
    dim::Int
    Larr::Vector{Int}
    Ns::Int
    Ndim::Int
    Npart::Int
    Nf::Int
    NSUN::Int
    ham_t::DTYPE
    ham_U::DTYPE
    trial_delta::DTYPE
end

function ham_par(; Larr = [8], ham_t = 1.0, ham_U = 4.0, trial_delta = 0.1, hs_channel = :spin)
    dim = length(Larr)
    if dim < 1 || dim > 2
        error("Only 1D and 2D Hubbard lattices are supported")
    end
    if any(L -> L < 1, Larr)
        error("All lattice sizes should be positive")
    end
    Ns = prod(Larr)
    if hs_channel == :charge
        Nf = 1
        NSUN = 2
    elseif hs_channel == :spin
        Nf = 2
        NSUN = 1
    else
        error("hs_channel must be :charge or :spin")
    end
    return ham_par(
        dim,
        collect(Larr),
        Ns,
        Ns,
        div(Ns, 2),
        Nf,
        NSUN,
        DTYPE(ham_t),
        DTYPE(ham_U),
        DTYPE(trial_delta),
    )
end

struct qmc_par
    Theta::DTYPE
    beta::DTYPE
    dtau::DTYPE
    Ntau::Int
    Mtau::Bool
    Nmes::Int
    Nmes0::Int
    Nv::Int
    NVdim::Int
    nstab::Int
    Naf::Int
    Nbin::Int
    Nsweep::Int
    Nana::Int
    rerun::Bool
    nblas::Int
    hs_channel::Symbol
end

function qmc_par(;
    dim = 1,
    Theta = 5.0,
    beta = 0.0,
    dtau = 0.05,
    Mtau = false,
    nstab = 10,
    Naf = 4,
    Nbin = 20,
    Nsweep = 20,
    Nana = Nbin,
    rerun = false,
    nblas = 1,
    hs_channel = :spin,
)
    Ntau = Int(round((Theta * 2 + beta) / dtau))
    Nmes0 = Int(round(Theta / dtau))
    Nmes = Int(round(beta / dtau))

    return qmc_par(
        DTYPE(Theta),
        DTYPE(beta),
        DTYPE(dtau),
        Ntau,
        Mtau,
        Nmes,
        Nmes0,
        1,
        1,
        nstab,
        Naf,
        Nbin,
        Nsweep,
        Nana,
        rerun,
        nblas,
        hs_channel,
    )
end

Base.@kwdef struct obs_eq
    sign::Vector{DTYPE}
    KE::Vector{DTYPE}
    PE::Vector{DTYPE}
    G0::Array{DTYPE, 4}
    SzSz::Array{DTYPE, 3}
    count::Vector{Int}
end

function obs_eq(hpar::ham_par, qpar::qmc_par)
    return obs_eq(
        sign = zeros(DTYPE, 1),
        KE = zeros(DTYPE, 1),
        PE = zeros(DTYPE, 1),
        G0 = zeros(DTYPE, 1, 1, hpar.Ns, hpar.Nf),
        SzSz = zeros(DTYPE, 1, 1, hpar.Ns),
        count = zeros(Int, 1),
    )
end

Base.@kwdef struct obs_tau
    sign::Vector{DTYPE}
    G00::Array{DTYPE, 4}
    GT0::Array{DTYPE, 5}
    G0T::Array{DTYPE, 5}
    count::Vector{Int}
end

function obs_tau(hpar::ham_par, qpar::qmc_par)
    return obs_tau(
        sign = zeros(DTYPE, 1),
        G00 = zeros(DTYPE, 1, 1, hpar.Ns, hpar.Nf),
        GT0 = zeros(DTYPE, 1, 1, hpar.Ns, hpar.Nf, qpar.Nmes + 1),
        G0T = zeros(DTYPE, 1, 1, hpar.Ns, hpar.Nf, qpar.Nmes + 1),
        count = zeros(Int, 1),
    )
end

struct hubbard_latt
    dim::Int
    Larr::Vector{Int}
    Ns::Int
    A::Matrix{DTYPE}
    B::Matrix{DTYPE}
    xpts::Matrix{DTYPE}
    kpts::Matrix{DTYPE}
    imj::Matrix{Int}
end

function i2vec(i::Int, Larr::Vector{Int}, dim::Int)
    if dim == 1
        return [i]
    elseif dim == 2
        return [div(i - 1, Larr[2]) + 1, mod1(i, Larr[2])]
    else
        error("Unsupported dimension")
    end
end

function vec2i(vec::Vector{Int}, Larr::Vector{Int}, dim::Int)
    if dim == 1
        return mod1(vec[1], Larr[1])
    elseif dim == 2
        return (mod1(vec[1], Larr[1]) - 1) * Larr[2] + mod1(vec[2], Larr[2])
    else
        error("Unsupported dimension")
    end
end

function hubbard_latt(hpar::ham_par)
    A = Matrix{DTYPE}(I, hpar.dim, hpar.dim)
    B = (2pi) .* inv(A')
    xpts = zeros(DTYPE, hpar.dim, hpar.Ns)
    kpts = zeros(DTYPE, hpar.dim, hpar.Ns)
    imj = zeros(Int, hpar.Ns, hpar.Ns)

    for i = 1:hpar.Ns
        vec = i2vec(i, hpar.Larr, hpar.dim)
        xpts[:, i] = A * (vec .- 1)
        kpts[:, i] = B * ((vec .- 1) ./ hpar.Larr)
    end

    for i = 1:hpar.Ns
        vi = i2vec(i, hpar.Larr, hpar.dim)
        for j = 1:hpar.Ns
            vj = i2vec(j, hpar.Larr, hpar.dim)
            imj[i, j] = vec2i(vi - vj .+ 1, hpar.Larr, hpar.dim)
        end
    end

    return hubbard_latt(hpar.dim, hpar.Larr, hpar.Ns, A, B, xpts, kpts, imj)
end

function sinv(A)
    return inv(lu(A))
end

@views @inbounds function prop_left_1x1!(G, v, fidx, i)
    for col = 1:size(G, 2)
        G[fidx, col, i] *= v
    end
end

@views @inbounds function prop_right_1x1!(G, mv, fidx, i)
    for row = 1:size(G, 1)
        G[row, fidx, i] *= mv
    end
end

@views @inbounds function GReq!(GR, Bl, Br, Nf)
    for i = 1:Nf
        GR[:, :, i] = Br[:, :, i] * sinv(Bl[:, :, i] * Br[:, :, i]) * Bl[:, :, i]
    end
end

@views @inbounds function Geq!(G, Bl, Br, Nf)
    for i = 1:Nf
        G[:, :, i] = I - Br[:, :, i] * sinv(Bl[:, :, i] * Br[:, :, i]) * Bl[:, :, i]
    end
end

@views @inbounds function recalc_G_old!(
    G,
    Bl,
    Br,
    itau,
    Pl,
    Pr,
    expT,
    expV,
    expVidx,
    afconf,
    Ntau,
    Ns,
    Nf,
    Nv,
    NVdim,
    sbatch,
)
    Bl[:] = Pl[:]
    Br[:] = Pr[:]
    tmp_Br = similar(Br[:, :, 1])
    tmp_Bl = similar(Bl[:, :, 1])

    for i = 1:Nf
        for itau_r = 1:itau
            mul!(tmp_Br, expT[:, :, i], Br[:, :, i])
            Br[:, :, i] .= tmp_Br
            for is = 1:Ns
                for iv = 1:Nv
                    Br[expVidx[:, iv, is], :, i] =
                        expV[:, :, i, iv, afconf[iv, is, itau_r]] * Br[expVidx[:, iv, is], :, i]
                end
            end

            if mod(itau_r, sbatch) == 0
                Br[:, :, i] = Matrix(qr(Br[:, :, i], ColumnNorm()).Q)
            end
        end
    end

    for i = 1:Nf
        for itau_l = Ntau:-1:itau + 1
            for is = Ns:-1:1
                for iv = Nv:-1:1
                    Bl[:, expVidx[:, iv, is], i] =
                        Bl[:, expVidx[:, iv, is], i] * expV[:, :, i, iv, afconf[iv, is, itau_l]]
                end
            end
            mul!(tmp_Bl, Bl[:, :, i], expT[:, :, i])
            Bl[:, :, i] .= tmp_Bl

            if mod(itau_l, sbatch) == 0
                Bl[:, :, i] = Matrix(qr(Bl[:, :, i]', ColumnNorm()).Q)'
            end
        end
    end

    Geq!(G, Bl, Br, Nf)
end

@views @inbounds function construct_Bl_Br!(
    Bl,
    Br,
    Bl_stab,
    Br_stab,
    nstab,
    Pl,
    Pr,
    expT,
    expV,
    expVidx,
    afconf,
    Ntau,
    Ns,
    Nf,
    Nv,
    NVdim,
)
    Bl[:] = Pl[:]
    Br[:] = Pr[:]
    tmp_Br = similar(Br[:, :, 1])
    tmp_Bl = similar(Bl[:, :, 1])

    Br_stab[:, :, :, 1] = Pr[:, :, :]
    for i = 1:Nf
        for itau_r = 1:Ntau
            mul!(tmp_Br, expT[:, :, i], Br[:, :, i])
            Br[:, :, i] .= tmp_Br
            for is = 1:Ns
                for iv = 1:Nv
                    Br[expVidx[:, iv, is], :, i] =
                        expV[:, :, i, iv, afconf[iv, is, itau_r]] * Br[expVidx[:, iv, is], :, i]
                end
            end

            if mod(itau_r, nstab) == 0
                Br[:, :, i] = Matrix(qr(Br[:, :, i], ColumnNorm()).Q)
                Br_stab[:, :, i, div(itau_r, nstab) + 1] = Br[:, :, i]
            end
        end
    end

    Bl_stab[:, :, :, end] = Pl[:, :, :]
    for i = 1:Nf
        for itau_l = Ntau:-1:1
            for is = Ns:-1:1
                for iv = Nv:-1:1
                    Bl[:, expVidx[:, iv, is], i] =
                        Bl[:, expVidx[:, iv, is], i] * expV[:, :, i, iv, afconf[iv, is, itau_l]]
                end
            end
            mul!(tmp_Bl, Bl[:, :, i], expT[:, :, i])
            Bl[:, :, i] .= tmp_Bl

            if mod(itau_l - 1, nstab) == 0
                Bl[:, :, i] = Matrix(qr(Bl[:, :, i]', ColumnNorm()).Q)'
                Bl_stab[:, :, i, div(itau_l - 1, nstab) + 1] = Bl[:, :, i]
            end
        end
    end
end

@views @inbounds function stable_G!(
    G,
    dir,
    Bl_stab,
    Br_stab,
    Pl,
    Pr,
    itau,
    nstab,
    expT,
    expV,
    expVidx,
    afconf,
    Ntau,
    Ns,
    Nf,
    Nv,
    NVdim,
)
    istab = div(itau, nstab)
    if mod(itau, nstab) != 0
        error("itau should be a multiple of nstab")
    end

    tmp_Br = similar(Br_stab[:, :, 1, 1])
    tmp_Bl = similar(Bl_stab[:, :, 1, 1])

    if dir == 1
        if istab > 0
            Br_stab[:, :, :, istab + 1] = Br_stab[:, :, :, istab]
            for i = 1:Nf
                for itau_r = (istab - 1) * nstab + 1:istab * nstab
                    mul!(tmp_Br, expT[:, :, i], Br_stab[:, :, i, istab + 1])
                    Br_stab[:, :, i, istab + 1] .= tmp_Br
                    for is = 1:Ns
                        for iv = 1:Nv
                            Br_stab[expVidx[:, iv, is], :, i, istab + 1] =
                                expV[:, :, i, iv, afconf[iv, is, itau_r]] *
                                Br_stab[expVidx[:, iv, is], :, i, istab + 1]
                        end
                    end
                end
                Br_stab[:, :, i, istab + 1] = Matrix(qr(Br_stab[:, :, i, istab + 1], ColumnNorm()).Q)
            end
        end
    elseif dir == 0
        if istab < div(Ntau, nstab)
            Bl_stab[:, :, :, istab + 1] = Bl_stab[:, :, :, istab + 2]
            for i = 1:Nf
                for itau_l = (istab + 1) * nstab:-1:istab * nstab + 1
                    for is = Ns:-1:1
                        for iv = Nv:-1:1
                            Bl_stab[:, expVidx[:, iv, is], i, istab + 1] =
                                Bl_stab[:, expVidx[:, iv, is], i, istab + 1] *
                                expV[:, :, i, iv, afconf[iv, is, itau_l]]
                        end
                    end
                    mul!(tmp_Bl, Bl_stab[:, :, i, istab + 1], expT[:, :, i])
                    Bl_stab[:, :, i, istab + 1] .= tmp_Bl
                end
                Bl_stab[:, :, i, istab + 1] = Matrix(qr(Bl_stab[:, :, i, istab + 1]', ColumnNorm()).Q)'
            end
        end
    else
        error("dir should be 0 or 1")
    end

    Geq!(G, Bl_stab[:, :, :, istab + 1], Br_stab[:, :, :, istab + 1], Nf)
end

@views @inbounds function recalc_G_stable!(
    G,
    Bl,
    Br,
    Bl_stab,
    Br_stab,
    itau,
    nstab,
    Pl,
    Pr,
    expT,
    expV,
    expVidx,
    afconf,
    Ntau,
    Ns,
    Nf,
    Nv,
    NVdim,
)
    istab = div(itau, nstab)
    jstab = mod(itau, nstab)
    tmp_Br = similar(Br_stab[:, :, 1, 1])
    tmp_Bl = similar(Bl_stab[:, :, 1, 1])

    Br[:] = Br_stab[:, :, :, istab + 1]
    if jstab == 0
        Bl[:] = Bl_stab[:, :, :, istab + 1]
    else
        Bl[:] = Bl_stab[:, :, :, istab + 2]
        for i = 1:Nf
            for itau_r = istab * nstab + 1:itau
                mul!(tmp_Br, expT[:, :, i], Br[:, :, i])
                Br[:, :, i] .= tmp_Br
                for is = 1:Ns
                    for iv = 1:Nv
                        Br[expVidx[:, iv, is], :, i] =
                            expV[:, :, i, iv, afconf[iv, is, itau_r]] * Br[expVidx[:, iv, is], :, i]
                    end
                end
            end
        end

        for i = 1:Nf
            for itau_l = (istab + 1) * nstab:-1:itau + 1
                for is = Ns:-1:1
                    for iv = Nv:-1:1
                        Bl[:, expVidx[:, iv, is], i] =
                            Bl[:, expVidx[:, iv, is], i] * expV[:, :, i, iv, afconf[iv, is, itau_l]]
                    end
                end
                mul!(tmp_Bl, Bl[:, :, i], expT[:, :, i])
                Bl[:, :, i] .= tmp_Bl
            end
        end
    end

    Geq!(G, Bl, Br, Nf)
end

@views @inbounds function propose_r_1x1(
    af_new::Int,
    afconf::Array{Int, 3},
    af_gam,
    itau::Int,
    is::Int,
    iv::Int,
    G,
    expV,
    expmV,
    expVidx::Array{Int, 3},
    Nf::Int,
    NSUN::Int,
    Δ,
)
    af_old = afconf[iv, is, itau]
    fidx = expVidx[1, iv, is]
    r = one(eltype(Δ))
    for i = 1:Nf
        delta = expV[1, 1, i, iv, af_new] * expmV[1, 1, i, iv, af_old] - 1.0
        Δ[1, 1, i] = delta
        r *= 1.0 + delta * (1.0 - G[fidx, fidx, i])
    end
    r = r^NSUN
    r *= af_gam[af_new] / af_gam[af_old]
    return r, fidx
end

@views @inbounds function update_G_1x1!(G, Δ, fidx, Nf)
    Ndim = size(G, 1)
    for i = 1:Nf
        delta = Δ[1, 1, i]
        g_loc = G[fidx, fidx, i]
        d = delta / (1.0 + delta * (1.0 - g_loc))

        for y = 1:Ndim
            if y == fidx
                continue
            end
            t = d * G[fidx, y, i]
            for x = 1:Ndim
                G[x, y, i] += G[x, fidx, i] * t
            end
        end

        scale = 1.0 - d * (1.0 - g_loc)
        for x = 1:Ndim
            G[x, fidx, i] *= scale
        end
    end
end

@views @inbounds function kinetic_energy(G0, Tmat, sign, NSUN::Int)
    ret = 0.0
    for k = 1:size(G0, 3)
        for i = 1:size(G0, 1)
            for j = 1:size(G0, 2)
                ret += sign * (((i == j) ? 1.0 : 0.0) - G0[i, j, k]) * Tmat[j, i, k]
            end
        end
    end
    return NSUN * ret
end

@views @inbounds function potential_energy(G0, hpar::ham_par, qpar::qmc_par, sign)
    ret = 0.0
    if qpar.hs_channel == :charge
        for i = 1:hpar.Ns
            nup = 1.0 - G0[i, i, 1]
            ndn = 1.0 - G0[i, i, 1]
            ret += sign * hpar.ham_U * nup * ndn
        end
    else
        for i = 1:hpar.Ns
            ret += sign * hpar.ham_U * (1.0 - G0[i, i, 1]) * (1.0 - G0[i, i, 2])
        end
    end
    return ret
end

@views @inbounds function szsz_corr(G0, P0, i, j)
    pup_i = P0[i, i, 1]
    pdn_i = P0[i, i, 2]
    pup_j = P0[j, j, 1]
    pdn_j = P0[j, j, 2]

    if i == j
        return 0.25 * (pup_i + pdn_i - 2.0 * pup_i * pdn_i)
    end

    nn_uu = pup_i * pup_j - P0[j, i, 1] * G0[i, j, 1]
    nn_dd = pdn_i * pdn_j - P0[j, i, 2] * G0[i, j, 2]
    nn_ud = pup_i * pdn_j
    nn_du = pdn_i * pup_j
    return 0.25 * (nn_uu + nn_dd - nn_ud - nn_du)
end

@views @inbounds function measure_obs_eq!(
    obeq::obs_eq,
    hpar::ham_par,
    qpar::qmc_par,
    latt::hubbard_latt,
    Tmat,
    sign,
    G0,
)
    P0 = similar(G0)
    for i = 1:hpar.Nf
        P0[:, :, i] .= I - G0[:, :, i]
    end

    obeq.sign .+= sign
    obeq.KE .+= kinetic_energy(G0, Tmat, sign, hpar.NSUN)
    obeq.PE .+= potential_energy(G0, hpar, qpar, sign)
    for j = 1:hpar.Ns
        for i = 1:hpar.Ns
            idx = latt.imj[i, j]
            obeq.G0[1, 1, idx, :] .+= sign .* G0[i, j, :]
            if qpar.hs_channel == :charge
                nup_i = P0[i, i, 1]
                nup_j = P0[j, j, 1]
                ndn_i = P0[i, i, 1]
                ndn_j = P0[j, j, 1]

                if i == j
                    obeq.SzSz[1, 1, idx] += sign * 0.25 * (nup_i + ndn_i - 2.0 * nup_i * ndn_i)
                else
                    nn_uu = nup_i * nup_j - P0[j, i, 1] * G0[i, j, 1]
                    nn_dd = ndn_i * ndn_j - P0[j, i, 1] * G0[i, j, 1]
                    nn_ud = nup_i * ndn_j
                    nn_du = ndn_i * nup_j
                    obeq.SzSz[1, 1, idx] += sign * 0.25 * (nn_uu + nn_dd - nn_ud - nn_du)
                end
            elseif qpar.hs_channel == :spin
                obeq.SzSz[1, 1, idx] += sign * szsz_corr(G0, P0, i, j)
            end
        end
    end
    obeq.count .+= 1
end

@views @inbounds function measure_obs_tau!(
    obtau::obs_tau,
    hpar::ham_par,
    qpar::qmc_par,
    latt::hubbard_latt,
    sign,
    G0,
    afconf,
    Bl_stab,
    Br_stab,
    expT,
    expmT,
    expV,
    expmV,
    expVidx,
)
    GT0 = copy(G0)
    G0T = similar(G0)
    GTT = copy(G0)
    B12 = similar(G0)
    for spin = 1:hpar.Nf
        G0T[:, :, spin] .= -(I - G0[:, :, spin])
    end

    for it = qpar.Nmes0:(qpar.Nmes0 + qpar.Nmes)
        it_idx = it - qpar.Nmes0 + 1

        if it_idx == 1
            for j = 1:hpar.Ns
                for i = 1:hpar.Ns
                    idx = latt.imj[i, j]
                    obtau.G00[1, 1, idx, :] .+= sign .* G0[i, j, :]
                end
            end
        else
            for spin = 1:hpar.Nf
                mul!(B12[:, :, spin], expT[:, :, spin], GT0[:, :, spin])
                GT0[:, :, spin] .= B12[:, :, spin]
                mul!(B12[:, :, spin], G0T[:, :, spin], expmT[:, :, spin])
                G0T[:, :, spin] .= B12[:, :, spin]
                mul!(B12[:, :, spin], expT[:, :, spin], GTT[:, :, spin])
                mul!(GTT[:, :, spin], B12[:, :, spin], expmT[:, :, spin])
            end

            for is = 1:hpar.Ns
                for iv = 1:qpar.Nv
                    af_curr = afconf[iv, is, it]
                    fidx = expVidx[1, iv, is]
                    for spin = 1:hpar.Nf
                        v = expV[1, 1, spin, iv, af_curr]
                        mv = expmV[1, 1, spin, iv, af_curr]
                        prop_left_1x1!(GT0, v, fidx, spin)
                        prop_left_1x1!(GTT, v, fidx, spin)
                        prop_right_1x1!(G0T, mv, fidx, spin)
                        prop_right_1x1!(GTT, mv, fidx, spin)
                    end
                end
            end
        end

        if mod(it, qpar.nstab) == 0
            istab = div(it, qpar.nstab)
            Geq!(GTT, Bl_stab[:, :, :, istab + 1], Br_stab[:, :, :, istab + 1], hpar.Nf)
            for spin = 1:hpar.Nf
                GT0[:, :, spin] .= GTT[:, :, spin] * GT0[:, :, spin]
                G0T[:, :, spin] .= G0T[:, :, spin] * (I - GTT[:, :, spin])
            end
        end

        for j = 1:hpar.Ns
            for i = 1:hpar.Ns
                idx = latt.imj[i, j]
                obtau.GT0[1, 1, idx, :, it_idx] .+= sign .* GT0[i, j, :]
                obtau.G0T[1, 1, idx, :, it_idx] .+= sign .* G0T[i, j, :]
            end
        end
    end

    obtau.sign .+= sign
    obtau.count .+= 1
end

function obs_reset!(obs)
    for fname in fieldnames(typeof(obs))
        field = getfield(obs, fname)
        fill!(field, zero(eltype(field)))
    end
end

function obs_avg!(obs)
    n = obs.count[1]
    if n == 0
        return
    end
    for fname in fieldnames(typeof(obs))
        if fname != :count
            getfield(obs, fname) ./= n
        end
    end
end

function obs_mpi_avg(obs, hpar::ham_par, qpar::qmc_par)
    irank = MPI.Comm_rank(MPI.COMM_WORLD)
    nrank = MPI.Comm_size(MPI.COMM_WORLD)
    obs_avg = typeof(obs)(hpar, qpar)

    for name in fieldnames(typeof(obs))
        arr = getfield(obs, name)
        arr_avg = getfield(obs_avg, name)
        MPI.Reduce!(arr, arr_avg, MPI.SUM, 0, MPI.COMM_WORLD)
        if irank == 0
            arr_avg ./= nrank
        end
    end

    return obs_avg
end

function obs_realspace_to_kspace(obs::AbstractArray{<:Number, 3}, latt::hubbard_latt)
    out = zeros(ComplexF64, size(obs))
    for k = 1:latt.Ns
        for i = 1:latt.Ns
            phase = exp(-1im * dot(latt.kpts[:, k], latt.xpts[:, i])) / latt.Ns
            out[1, 1, k] += phase * obs[1, 1, i]
        end
    end
    return out
end

function obs_realspace_to_kspace(obs::AbstractArray{<:Number, 4}, latt::hubbard_latt)
    out = zeros(ComplexF64, size(obs))
    for spin = 1:size(obs, 4)
        for k = 1:latt.Ns
            for i = 1:latt.Ns
                phase = exp(-1im * dot(latt.kpts[:, k], latt.xpts[:, i])) / latt.Ns
                out[1, 1, k, spin] += phase * obs[1, 1, i, spin]
            end
        end
    end
    return out
end

function obs_realspace_to_kspace(obs::AbstractArray{<:Number, 5}, latt::hubbard_latt)
    out = zeros(ComplexF64, size(obs))
    for it = 1:size(obs, 5)
        for spin = 1:size(obs, 4)
            for k = 1:latt.Ns
                for i = 1:latt.Ns
                    phase = exp(-1im * dot(latt.kpts[:, k], latt.xpts[:, i])) / latt.Ns
                    out[1, 1, k, spin, it] += phase * obs[1, 1, i, spin, it]
                end
            end
        end
    end
    return out
end

function obs_h5_init(fname::String, obs, rerun::Bool, latt::hubbard_latt)
    if isfile(fname) && rerun
        return
    end

    h5open(fname, "w") do file
        for name in fieldnames(typeof(obs))
            val = getfield(obs, name)
            dims = size(val)
            init_dims = (dims..., 0)
            max_dims = (dims..., -1)
            chunk_size = (dims..., 1)
            dspace = dataspace(init_dims, max_dims = max_dims)
            create_dataset(file, string(name), eltype(val), dspace, chunk = chunk_size)

            if ndims(val) >= 3 && size(val, 3) == latt.Ns
                val_k = obs_realspace_to_kspace(val, latt)
                dims_k = size(val_k)
                init_dims_k = (dims_k..., 0)
                max_dims_k = (dims_k..., -1)
                chunk_size_k = (dims_k..., 1)
                dspace_k = dataspace(init_dims_k, max_dims = max_dims_k)
                create_dataset(file, string(name) * "k", eltype(val_k), dspace_k, chunk = chunk_size_k)
            end
        end
    end
end

function obs_h5_append(fname::String, obs, latt::hubbard_latt)
    h5open(fname, "r+") do file
        for name in fieldnames(typeof(obs))
            val = getfield(obs, name)
            dset = file[string(name)]
            new_bins = size(dset, ndims(dset)) + 1
            HDF5.set_extent_dims(dset, (size(val)..., new_bins))
            dset[ntuple(_ -> Colon(), ndims(val))..., new_bins] = val

            if ndims(val) >= 3 && size(val, 3) == latt.Ns
                val_k = obs_realspace_to_kspace(val, latt)
                dset_k = file[string(name) * "k"]
                new_bins_k = size(dset_k, ndims(dset_k)) + 1
                HDF5.set_extent_dims(dset_k, (size(val_k)..., new_bins_k))
                dset_k[ntuple(_ -> Colon(), ndims(val_k))..., new_bins_k] = val_k
            end
        end
    end
end

function jackknife_mean(data)
    dim = ndims(data)
    nbins = size(data, dim)
    mean_data = mean(data, dims = dim)
    err_data = similar(mean_data)

    if nbins < 2
        err_data .= NaN
        return cat(mean_data, err_data, dims = dim)
    end

    total = sum(data, dims = dim)
    leave_one_out = (total .- data) ./ (nbins - 1)
    leave_one_out_mean = mean(leave_one_out, dims = dim)
    err_data .= sqrt.((nbins - 1) / nbins .* sum(abs2.(leave_one_out .- leave_one_out_mean), dims = dim))
    return cat(mean_data, err_data, dims = dim)
end

function jackknife_ratio(num, den)
    dim = ndims(num)
    nbins = size(num, dim)
    den_bins = reshape(den, ntuple(i -> i == dim ? size(den, ndims(den)) : 1, dim))
    ratio_mean = sum(num, dims = dim) ./ sum(den_bins, dims = dim)
    ratio_err = similar(ratio_mean)

    if nbins < 2
        ratio_err .= NaN
        return cat(ratio_mean, ratio_err, dims = dim)
    end

    num_total = sum(num, dims = dim)
    den_total = sum(den_bins, dims = dim)
    leave_one_out = (num_total .- num) ./ (den_total .- den_bins)
    leave_one_out_mean = mean(leave_one_out, dims = dim)
    ratio_err .= sqrt.((nbins - 1) / nbins .* sum(abs2.(leave_one_out .- leave_one_out_mean), dims = dim))
    return cat(ratio_mean, ratio_err, dims = dim)
end

function analysis(hpar::ham_par, qpar::qmc_par)
    mkpath("analysis")
    h5open("data/data.h5", "r") do in_f
        h5open("analysis/data.h5", "w") do out_f
            sign_data = read(in_f["sign"])
            for dname in keys(in_f)
                data = read(in_f[dname])
                if dname == "sign" || dname == "count"
                    out_f[dname] = jackknife_mean(data)
                else
                    out_f[dname] = jackknife_ratio(data, sign_data)
                end
            end
        end
    end

    if qpar.Mtau && isfile("data/data_tau.h5")
        h5open("data/data_tau.h5", "r") do in_f
            h5open("analysis/data_tau.h5", "w") do out_f
                sign_data = read(in_f["sign"])
                for dname in keys(in_f)
                    data = read(in_f[dname])
                    if dname == "sign" || dname == "count"
                        out_f[dname] = jackknife_mean(data)
                    else
                        out_f[dname] = jackknife_ratio(data, sign_data)
                    end
                end
            end
        end
    end
end

function discrete_hs_parameters(qpar::qmc_par)
    if qpar.Naf == 2
        af_gam = DTYPE[1.0, 1.0]
        af_eta = DTYPE[1.0, -1.0]
    elseif qpar.Naf == 4
        af_gam = DTYPE[
            0.25 * (1 - sqrt(6) / 3),
            0.25 * (1 - sqrt(6) / 3),
            0.25 * (1 + sqrt(6) / 3),
            0.25 * (1 + sqrt(6) / 3),
        ]
        af_eta = DTYPE[
            sqrt(2 * (3 + sqrt(6))),
            -sqrt(2 * (3 + sqrt(6))),
            sqrt(2 * (3 - sqrt(6))),
            -sqrt(2 * (3 - sqrt(6))),
        ]
    else
        error("Naf should equals to 2 or 4")
    end
    return af_gam, af_eta
end

function build_hopping!(Tmat, Ttrial, hpar::ham_par)
    fill!(Tmat, 0.0)
    fill!(Ttrial, 0.0)

    if hpar.dim == 1
        L = hpar.Larr[1]
        for spin = 1:hpar.Nf
            for x = 1:L
                x1 = mod1(x + 1, L)
                Tmat[x, x1, spin] = -hpar.ham_t
                Tmat[x1, x, spin] = -hpar.ham_t

                δ = hpar.trial_delta * (-1)^x
                Ttrial[x, x1, spin] = -hpar.ham_t + δ
                Ttrial[x1, x, spin] = -hpar.ham_t + δ
            end
        end
    elseif hpar.dim == 2
        Lx, Ly = hpar.Larr
        for spin = 1:hpar.Nf
            for x = 1:Lx
                for y = 1:Ly
                    i = (x - 1) * Ly + y
                    ix = (mod1(x + 1, Lx) - 1) * Ly + y
                    iy = (x - 1) * Ly + mod1(y + 1, Ly)

                    Tmat[i, ix, spin] = -hpar.ham_t
                    Tmat[ix, i, spin] = -hpar.ham_t
                    Tmat[i, iy, spin] = -hpar.ham_t
                    Tmat[iy, i, spin] = -hpar.ham_t

                    δx = hpar.trial_delta * cos(pi * (x + y))
                    δy = hpar.trial_delta
                    Ttrial[i, ix, spin] = -hpar.ham_t * (1.0 + δx)
                    Ttrial[ix, i, spin] = -hpar.ham_t * (1.0 + δx)
                    Ttrial[i, iy, spin] = -hpar.ham_t * (1.0 - δy)
                    Ttrial[iy, i, spin] = -hpar.ham_t * (1.0 - δy)
                end
            end
        end
    else
        error("Larr should be 1 or 2 dimension")
    end
end

function build_interaction!(Vmat, expV, expmV, expVidx, hpar::ham_par, qpar::qmc_par, af_eta)
    fill!(Vmat, 0.0)
    if qpar.hs_channel == :charge
        if qpar.Naf == 2
            λ = acosh(exp(-qpar.dtau * hpar.ham_U / 2.0))
        elseif qpar.Naf == 4
            λ = sqrt(-qpar.dtau * hpar.ham_U / 2.0)
        else
            error("Naf should equals to 2 or 4")
        end
        for spin = 1:hpar.Nf
            Vmat[1, 1, spin, 1] = λ
        end
    elseif qpar.hs_channel == :spin
        if qpar.Naf == 2
            λ = acosh(exp(qpar.dtau * hpar.ham_U / 2.0))
        elseif qpar.Naf == 4
            λ = sqrt(qpar.dtau * hpar.ham_U / 2.0)
        else
            error("Naf should equals to 2 or 4")
        end
        Vmat[1, 1, 1, 1] = λ
        for spin = 2:hpar.Nf
            Vmat[1, 1, spin, 1] = -λ
        end
    else
        error("hs_channel should equals to :spin or :charge")
    end

    for spin = 1:hpar.Nf
        for iv = 1:qpar.Nv
            for iaf = 1:qpar.Naf
                expV[:, :, spin, iv, iaf] .= exp.(af_eta[iaf] .* Vmat[:, :, spin, iv])
                expmV[:, :, spin, iv, iaf] .= exp.(-af_eta[iaf] .* Vmat[:, :, spin, iv])
            end
        end
    end

    for is = 1:hpar.Ns
        for iv = 1:qpar.Nv
            expVidx[:, iv, is] .= is
        end
    end
end

function initialize_trial_state!(Pl, Pr, Ttrial, hpar::ham_par)
    eig = eigen(Hermitian(Ttrial[:, :, 1]))
    for spin = 1:hpar.Nf
        Pr[:, :, spin] = eig.vectors[:, 1:hpar.Npart]
        Pl[:, :, spin] = Pr[:, :, spin]'
    end
    return eig.values[hpar.Npart + 1]
end
