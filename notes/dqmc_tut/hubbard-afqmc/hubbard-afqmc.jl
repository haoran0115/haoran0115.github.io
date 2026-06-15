using LinearAlgebra
using Statistics
using Random
using Base.Threads
using Printf
using JLD2
using HDF5
using MPI

const DTYPE = ComplexF64
const G_STAB_ERR = 1e-3

# formatting helper
# string format function for real and complex number
function fmtc(z)
    if abs(imag(z)) < 1e-10
        return @sprintf("%+10.6f", real(z))
    end
    return @sprintf("%+10.6f%+10.6fim", real(z), imag(z))
end

# data structures and constructors
# ham_par data struct: parameters related to Hamiltonian and system size
struct ham_par
    # system geometries
    dim::Int           # system dimension
    Larr::Vector{Int}  # system size, Larr = [L1, ..., Ldim], dim = 1, 2
    Ns::Int            # number of unit cells, Ns = L1 * L2 * ... * Ldim
    Ndim::Int          # dimension of fermion matrices, in the Hubbard model Ndim = Ns
    Npart::Int         # number of particles (per spin index), the code assumes N↑ = N↓

    # flavor symmetries
    Nf::Int            # number of fermion flavors (spin), Nf = 2 (↑ and ↓) if the HS decoupling breaks the SU(2) symmetry
    NSUN::Int          # N in SUN symmetry, for SU(2)-symmetric HS decoupling, NSUN = 2

    # Hamiltonian parameters
    ham_t::DTYPE       # hopping t in Hubbard model
    ham_U::DTYPE       # interaction U in Hubbard model
    ham_mu::DTYPE      # chemical potential μ in Hubbard model

    # trial wavefunction parameters
    trial_delta::DTYPE # trial wavefunction parameter, which opens a gap in the trial kinetic matrix to avoid degeneracies
    hs_channel::Symbol
    Nmes::Int
    Nmes0::Int
end

# ham_par initialization function
function ham_par(; Larr = [8], ham_t = 1.0, ham_U = 4.0, ham_mu = 0.0, trial_delta = 0.1, hs_channel = :spin, Theta = 5.0, beta = 0.0, dtau = 0.05)
    dim = length(Larr)
    if dim < 1 || dim > 2
        error("Only 1D and 2D Hubbard lattices are supported")
    end
    if any(L -> L < 1, Larr)
        error("All lattice sizes should be positive")
    end
    Ns = prod(Larr)
    if hs_channel == :SU2
        Nf = 1
        NSUN = 2
    elseif hs_channel == :spin
        Nf = 2
        NSUN = 1
    else
        error("hs_channel must be :SU2 or :spin")
    end
    Nmes0 = Int(round(Theta / dtau))
    Nmes = Int(round(beta / dtau))
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
        DTYPE(ham_mu),
        DTYPE(trial_delta),
        hs_channel,
        Nmes,
        Nmes0,
    )
end

# qmc_par data struct: quantum monte carlo simulation parameters
struct qmc_par
    # imaginary time projection parameters
    # ⟨O⟩ = ⟨Ψ| e^{-(Θ+β/2)H} O e^{-(Θ+β/2)H} |Ψ⟩ / ⟨Ψ| e^{-(2Θ+β)H} |Ψ⟩
    Theta::DTYPE  # Θ
    beta::DTYPE   # β
    dtau::DTYPE   # imaginary time discritization Δτ
    Ntau::Int     # Nτ = (2Θ+β)/Δτ

    # measurement variables
    Mtau::Bool    # if measure time-dependent variables or not
    Nv::Int       # number of interaction vertices per unit cell
    NVdim::Int    # dimension of interaction vertices
    nstab::Int    # perform stablization every nstab*Δτ time
    Naf::Int      # number of auxiliary field components, 2 for Ising-type and 4 for Gauss-Hermite type
    Nbin::Int     # number of statistical bins
    Nsweep::Int   # number of sweeps (forward+backward) within one bin
    Nana::Int     # postprocessing and analysis will be performed every Nana bins
    rerun::Bool   # if resume from the previous run or not
    nblas::Int    # number of BLAS/MKL threads per process
end

# qmc_par initialization function
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
)
    Ntau = Int(round((Theta * 2 + beta) / dtau))

    return qmc_par(
        DTYPE(Theta),
        DTYPE(beta),
        DTYPE(dtau),
        Ntau,
        Mtau,
        1,
        1,
        nstab,
        Naf,
        Nbin,
        Nsweep,
        Nana,
        rerun,
        nblas,
    )
end

# equal-time observables struct
Base.@kwdef struct obs_eq
    sign::Vector{DTYPE}    # averaged sign of fermion weight Re W[C]
    KE::Vector{DTYPE}      # kinetic energy
    PE::Vector{DTYPE}      # potential energy
    G0::Array{DTYPE, 4}    # time-equal Green's function
    SzSz::Array{DTYPE, 3}  # SzSz[i-j] = ⟨Sz_i Sz_j⟩, where Sz_i = ni↑ - ni↓
    count::Vector{Int}     # number of measurements
end

# obs_eq initialization function
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

# time-displaced observables struct
Base.@kwdef struct obs_tau
    sign::Vector{DTYPE}  # averaged sign of fermion weight Re W[C]
    G00::Array{DTYPE, 4} #  ⟨ci(0) cj†(0)⟩
    GT0::Array{DTYPE, 5} #  ⟨ci(t) cj†(0)⟩
    G0T::Array{DTYPE, 5} # -⟨cj†(t) ci(0)⟩
    count::Vector{Int}
end

# obs_tau initialization function
function obs_tau(hpar::ham_par, qpar::qmc_par)
    return obs_tau(
        sign = zeros(DTYPE, 1),
        G00 = zeros(DTYPE, 1, 1, hpar.Ns, hpar.Nf),
        GT0 = zeros(DTYPE, 1, 1, hpar.Ns, hpar.Nf, hpar.Nmes + 1),
        G0T = zeros(DTYPE, 1, 1, hpar.Ns, hpar.Nf, hpar.Nmes + 1),
        count = zeros(Int, 1),
    )
end

# lattice object, useful for lattice observable collections
struct hubbard_latt
    dim::Int            # dimesion of lattice
    Larr::Vector{Int}   # [L1, L2, ..., Ldim]
    Ns::Int             # number of unit cells
    A::Matrix{DTYPE}    # A = [a1, a2, ..., adim] the real-space basis vectors
    B::Matrix{DTYPE}    # B = [b1, b2, ..., bdim] = 2πA^{-1} the reciprocal space basis vectors
    xpts::Matrix{DTYPE} # real space x points: x1, x2, ..., xNs
    kpts::Matrix{DTYPE} # reciprocal space k points
    imj::Matrix{Int}    # index of xpts[i] - xpts[j] in xpts
                        # e.g. xpts[imj[i, j]] = xpts[i] - xpts[j]
end

# maps x index to x vector
function i2vec(i::Int, Larr::Vector{Int}, dim::Int)
    if dim == 1
        return [i]
    elseif dim == 2
        return [div(i - 1, Larr[2]) + 1, mod1(i, Larr[2])]
    else
        error("Unsupported dimension")
    end
end

# maps x vector to x index
function vec2i(vec::Vector{Int}, Larr::Vector{Int}, dim::Int)
    if dim == 1
        return mod1(vec[1], Larr[1])
    elseif dim == 2
        return (mod1(vec[1], Larr[1]) - 1) * Larr[2] + mod1(vec[2], Larr[2])
    else
        error("Unsupported dimension")
    end
end

# hubbard_latt initialization
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

# model construction helpers
# returns the discrete HS quadrature weights γ and nodes η
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

# builds the hopping matrix and the trial Hamiltonian, including the chemical potential
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

                Tmat[x, x, spin] += -hpar.ham_mu
                Ttrial[x, x, spin] += -hpar.ham_mu
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

                    Tmat[i, i, spin] += -hpar.ham_mu
                    Ttrial[i, i, spin] += -hpar.ham_mu
                end
            end
        end
    else
        error("Larr should be 1 or 2 dimension")
    end
end

# builds all local interaction matrices used in the HS update
function build_interaction!(Vmat, Vconst, expV, expmV, expVidx, hpar::ham_par, qpar::qmc_par, af_eta)
    fill!(Vmat, 0.0)
    fill!(Vconst, 0.0)
    if hpar.hs_channel == :SU2
        if qpar.Naf == 2
            λ = acosh(exp(-qpar.dtau * hpar.ham_U / hpar.NSUN))
        elseif qpar.Naf == 4
            λ = sqrt(-qpar.dtau * hpar.ham_U / hpar.NSUN)
        else
            error("Naf should equals to 2 or 4")
        end
        for spin = 1:hpar.Nf
            Vmat[1, 1, spin, 1] = λ
            Vconst[spin, 1] = -0.5 * λ
        end
    elseif hpar.hs_channel == :spin
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
        error("hs_channel should equals to :spin or :SU2")
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

# initializes the left and right trial Slater determinants
function initialize_trial_state!(Pl, Pr, Ttrial, hpar::ham_par)
    eig = eigen(Hermitian(Ttrial[:, :, 1]))
    for spin = 1:hpar.Nf
        Pr[:, :, spin] = eig.vectors[:, 1:hpar.Npart]
        Pl[:, :, spin] = Pr[:, :, spin]'
    end
    return eig.values[hpar.Npart + 1]
end

# low-level propagation helpers

# numiercal-stable matrix inversion
function sinv(A)
    return inv(lu(A))
end

# propagate equal-time Green's function on the left
# G = e^v * G
@views @inbounds function prop_left_1x1!(G, v, fidx, i)
    for col in axes(G, 2)
        G[fidx, col, i] *= v
    end
end

# propagate equal-time Green's function on the right
# G = G * e^(-v)
@views @inbounds function prop_right_1x1!(G, mv, fidx, i)
    for row in axes(G, 1)
        G[row, fidx, i] *= mv
    end
end

# calculate equal-time Green's function from projectors
@views @inbounds function Geq!(G, Bl, Br, Nf)
    for i = 1:Nf
        G[:, :, i] = I - Br[:, :, i] * sinv(Bl[:, :, i] * Br[:, :, i]) * Bl[:, :, i]
    end
end

# legacy function to calculate the equal-time Green's function at an arbitrary time slice
# ⟨e^{-dtau*(Ntau-itau)} ci cj† e^{-dtau*itau}⟩
@views @inbounds function recalc_G_old!(
    G,         # equal-time Green's function at time slice itau
    Bl,        # left projected wavefunction B⟨
    Br,        # right projected wavefunction B⟩
    itau,      # target time slice
    Pl,        # ⟨ΨL|
    Pr,        # |ΨR⟩
    expT,      # e^(-Δτ T)
    expV,      # local interaction propagators
    expVidx,   # site index of each local interaction vertex
    afconf,    # auxiliary field configuration
    Ntau,      # number of time slices
    Ns,        # number of unit cells
    Nf,        # number of stored fermion flavors
    Nv,        # number of interaction vertices per unit cell
    NVdim,     # dimension of a local interaction vertex
    sbatch,    # QR re-orthonormalization period in time slices
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

# stabilized projector propagation helpers
# constructs the cached projected wavefunctions B⟨ and B⟩ at each stabilization slice
@views @inbounds function construct_Bl_Br!(
    Bl,         # Bl = B⟨
    Br,         # Br = B⟩
    Bl_stab,    # stablized Bl at each stablization time slices
    Br_stab,    # stablized Br at each stablization time slices
    nstab,      # perform stablization every nstab*Δτ time
    Pl,         # ⟨ΨL|
    Pr,         # |Ψ_R⟩
    expT,       # e^(-Δτ T)
    expV,       # local interaction propagators exp(+ηl * Vi)
    expVidx,    # index of Vi
    afconf,     # auxiliary field configuration
    Ntau,       # (2Θ+β)/Δτ
    Ns,         # number of spatial sites
    Nf,         # number of explicitly stored fermion flavor blocks
    Nv,         # number of interaction operators attached to each site
    NVdim,      # local matrix dimension of each interaction operator
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

# updates the cached projectors at one stabilization slice and recomputes G
@views @inbounds function stable_G!(
    G,         # equal-time Green's function at the stabilized slice
    dir,       # sweep direction, 1 for forward and 0 for backward
    Bl_stab,   # cached left projected wavefunctions
    Br_stab,   # cached right projected wavefunctions
    Pl,        # ⟨ΨL|
    Pr,        # |ΨR⟩
    itau,      # stabilization time slice
    nstab,     # perform stabilization every nstab*Δτ time
    expT,      # e^(-Δτ T)
    expV,      # local interaction propagators
    expVidx,   # site index of each local interaction vertex
    afconf,    # auxiliary field configuration
    Ntau,      # number of time slices
    Ns,        # number of unit cells
    Nf,        # number of stored fermion flavors
    Nv,        # number of interaction vertices per unit cell
    NVdim,     # dimension of a local interaction vertex
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

# reconstructs the equal-time Green's function from the nearest cached stabilization slice
@views @inbounds function recalc_G_stable!(
    G,         # equal-time Green's function at time slice itau
    Bl,        # left projected wavefunction B⟨
    Br,        # right projected wavefunction B⟩
    Bl_stab,   # cached left projected wavefunctions
    Br_stab,   # cached right projected wavefunctions
    itau,      # target time slice
    nstab,     # perform stabilization every nstab*Δτ time
    Pl,        # ⟨ΨL|
    Pr,        # |ΨR⟩
    expT,      # e^(-Δτ T)
    expV,      # local interaction propagators
    expVidx,   # site index of each local interaction vertex
    afconf,    # auxiliary field configuration
    Ntau,      # number of time slices
    Ns,        # number of unit cells
    Nf,        # number of stored fermion flavors
    Nv,        # number of interaction vertices per unit cell
    NVdim,     # dimension of a local interaction vertex
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

# local Metropolis update helpers
# computes the local Metropolis ratio for a proposed HS-field update
@views @inbounds function propose_r_1x1(
    af_new::Int,            # proposed new HS-field value
    afconf::Array{Int, 3},  # current auxiliary field configuration
    af_gam,                 # quadrature weights γ
    af_eta,                 # quadrature nodes η
    itau::Int,              # imaginary-time index
    is::Int,                # spatial index
    iv::Int,                # interaction-vertex index
    G,                      # equal-time Green's function before the update
    expV,                   # local interaction propagators
    expmV,                  # inverse local interaction propagators
    Vmat,                   # local HS-coupling matrices
    Vconst,                 # scalar HS shifts vi in V̂i = c† Vi c + vi
    expVidx::Array{Int, 3}, # site index of each local interaction vertex
    Nf::Int,                # number of stored fermion flavors
    NSUN::Int,              # SU(N) multiplicity represented by each stored block
    Δ,                      # workspace storing the local rank-one update
)
    af_old = afconf[iv, is, itau]
    fidx = expVidx[1, iv, is]
    det_ratio = one(eltype(Δ))
    const_ratio = one(eltype(Δ))
    for i = 1:Nf
        delta = expV[1, 1, i, iv, af_new] * expmV[1, 1, i, iv, af_old] - 1.0
        Δ[1, 1, i] = delta
        det_ratio *= 1.0 + delta * (1.0 - G[fidx, fidx, i])
        dη = af_eta[af_new] - af_eta[af_old]
        const_ratio *= exp(Vconst[i, iv] * dη)
    end
    r = (const_ratio * det_ratio)^NSUN
    r *= af_gam[af_new] / af_gam[af_old]
    return r, fidx
end

# updates G after an accepted local rank-one HS-field update
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

# equal-time measurement helpers
# calculates the kinetic energy contribution sign * tr[(I - G) * T]
@views @inbounds function kinetic_energy(G0, Tmat, sign, NSUN::Int)
    ret = zero(eltype(G0))
    for k in axes(G0, 3)
        Gk = G0[:, :, k]
        Tk = Tmat[:, :, k]
        ret += sign * tr((I - Gk) * Tk)
    end
    return NSUN * ret
end

# calculates the potential energy contribution for the chosen HS channel
@views @inbounds function potential_energy(G0, hpar::ham_par, qpar::qmc_par, sign)
    ret = 0.0
    if hpar.hs_channel == :SU2
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

# evaluates the equal-time spin correlation ⟨Sz_i Sz_j⟩ for the spin channel
@views @inbounds function szsz_corr(G0, P0, i, j)
    pup_i = P0[i, i, 1]
    pdn_i = P0[i, i, 2]
    pup_j = P0[j, j, 1]
    pdn_j = P0[j, j, 2]

    if i == j
        return 0.25 * (pup_i + pdn_i - 2.0 * pup_i * pdn_i)
    end

    nn_uu = pup_i * pup_j + P0[j, i, 1] * G0[i, j, 1]
    nn_dd = pdn_i * pdn_j + P0[j, i, 2] * G0[i, j, 2]
    nn_ud = pup_i * pdn_j
    nn_du = pdn_i * pup_j
    return 0.25 * (nn_uu + nn_dd - nn_ud - nn_du)
end

# accumulates one equal-time measurement into the observable container
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
            if hpar.hs_channel == :SU2
                nup_i = P0[i, i, 1]
                nup_j = P0[j, j, 1]
                ndn_i = P0[i, i, 1]
                ndn_j = P0[j, j, 1]

                if i == j
                    obeq.SzSz[1, 1, idx] += sign * 0.25 * (nup_i + ndn_i - 2.0 * nup_i * ndn_i)
                else
                    nn_uu = nup_i * nup_j + P0[j, i, 1] * G0[i, j, 1]
                    nn_dd = ndn_i * ndn_j + P0[j, i, 1] * G0[i, j, 1]
                    nn_ud = nup_i * ndn_j
                    nn_du = ndn_i * nup_j
                    obeq.SzSz[1, 1, idx] += sign * 0.25 * (nn_uu + nn_dd - nn_ud - nn_du)
                end
            elseif hpar.hs_channel == :spin
                obeq.SzSz[1, 1, idx] += sign * szsz_corr(G0, P0, i, j)
            end
        end
    end
    obeq.count .+= 1
end

# accumulates one time-displaced Green's function measurement
@views @inbounds function measure_obs_tau!(
    obtau::obs_tau,      # time-displaced observable container
    hpar::ham_par,       # Hamiltonian and system-size parameters
    qpar::qmc_par,       # QMC simulation parameters
    latt::hubbard_latt,  # lattice geometry and Fourier-transform metadata
    sign,                # current Monte Carlo sign / phase estimator
    G0,                  # equal-time Green's function at τ = Nmes0
    afconf,              # auxiliary field configuration
    Bl_stab,             # cached left projected wavefunctions
    Br_stab,             # cached right projected wavefunctions
    expT,                # e^(-Δτ T)
    expmT,               # e^(+Δτ T)
    expV,                # local interaction propagators
    expmV,               # inverse local interaction propagators
    expVidx,             # site index of each local interaction vertex
)
    GT0 = copy(G0)
    G0T = similar(G0)
    GTT = copy(G0)
    B12 = similar(G0)
    for spin = 1:hpar.Nf
        G0T[:, :, spin] .= -(I - G0[:, :, spin])
    end

    for it = hpar.Nmes0:(hpar.Nmes0 + hpar.Nmes)
        it_idx = it - hpar.Nmes0 + 1

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

# observable accumulation and I/O helpers
# resets all accumulated observable arrays to zero before a new bin starts
function obs_reset!(obs)
    for fname in fieldnames(typeof(obs))
        field = getfield(obs, fname)
        fill!(field, zero(eltype(field)))
    end
end

# divides accumulated observables by the number of measurements in the bin
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

# averages one observable container over all MPI ranks / Markov chains
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

# Fourier transforms a real-space two-point observable with no flavor or time indices
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

# Fourier transforms a real-space two-point observable with a flavor index
function obs_realspace_to_kspace(obs::AbstractArray{<:Number, 4}, latt::hubbard_latt)
    out = zeros(ComplexF64, size(obs))
    for spin in axes(obs, 4)
        for k = 1:latt.Ns
            for i = 1:latt.Ns
                phase = exp(-1im * dot(latt.kpts[:, k], latt.xpts[:, i])) / latt.Ns
                out[1, 1, k, spin] += phase * obs[1, 1, i, spin]
            end
        end
    end
    return out
end

# Fourier transforms a real-space two-point observable with flavor and time indices
function obs_realspace_to_kspace(obs::AbstractArray{<:Number, 5}, latt::hubbard_latt)
    out = zeros(ComplexF64, size(obs))
    for it in axes(obs, 5)
        for spin in axes(obs, 4)
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

# initializes the HDF5 datasets used to store raw bin-by-bin observables
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

# appends one bin of observables to the existing HDF5 datasets
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

# error analysis helpers
# computes the jackknife mean and error bar along the last dimension
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

# computes the jackknife estimate for a ratio observable num / den
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

# reads raw measurements and writes jackknife-processed observables to analysis/
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
