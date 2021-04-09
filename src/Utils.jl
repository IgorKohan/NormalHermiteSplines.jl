function _estimate_accuracy(spline::NormalSpline{T, RK}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    n = size(spline._nodes, 1)
    m = size(spline._nodes, 2)
    nodes = similar(spline._nodes)
    @inbounds for i = 1:n
        for j = 1:m
            nodes[i,j] = spline._min_bound[i] + spline._compression * spline._nodes[i,j]
        end
    end
    σ = evaluate(spline, nodes)
    # calculating a value of the Relative Maximum Absolute Error (RMAE) of interpolation
    # at the function value interpolation nodes.
    fun_max = maximum(abs.(spline._values))
    if fun_max > 0.0
        rmae = maximum(abs.(spline._values .- σ)) / fun_max
    else
        rmae = maximum(abs.(spline._values .- σ))
    end
    rmae = rmae > eps(T) ? rmae : eps(T)
    res = -floor(log10(rmae)) - 1
    if res <= 0
        res = 0
    end
    return trunc(Int, res)
end

function _estimate_ε(nodes::Matrix{T}) where T <: AbstractFloat
    n = size(nodes, 1)
    n_1 = size(nodes, 2)
    ε = T(0.0)
    @inbounds for i = 1:n_1
        for j = i:n_1
            ε += norm(nodes[:,i] .- nodes[:,j])
        end
    end
    if ε > T(0.0)
        ε *= T(n)^(1.0/n) / T(n_1)^(5.0/3.0)
    else
        ε = T(1.0)
    end
    return ε
end

function _estimate_ε(nodes::Matrix{T},
                     d_nodes::Matrix{T}
                    ) where T <: AbstractFloat
    return _estimate_ε([nodes 0.1 .* d_nodes])
end

function _estimate_epsilon(nodes::Matrix{T},
                           kernel::RK) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    n = size(nodes, 1)
    n_1 = size(nodes, 2)
    min_bound = Vector{T}(undef, n)
    compression::T = 0
    @inbounds for i = 1:n
        min_bound[i] = nodes[i,1]
        maxx::T = nodes[i,1]
        for j = 2:n_1
            min_bound[i] = min(min_bound[i], nodes[i,j])
            maxx = max(maxx, nodes[i,j])
        end
        compression = max(compression, maxx - min_bound[i])
    end

    if compression <= eps(T(1.0))
        error("Cannot estimate_epsilon: `nodes` data are not correct.")
    end

    t_nodes = similar(nodes)
    @inbounds for i = 1:n
        for j = 1:n_1
            t_nodes[i,j] = (nodes[i,j] - min_bound[i]) / compression
        end
    end
    ε = _estimate_ε(t_nodes)
    if isa(kernel, RK_H1)
        ε *= T(1.5)
    elseif isa(kernel, RK_H2)
        ε *= T(2.0)
    end
    return ε
end

function _estimate_epsilon(nodes::Matrix{T},
                           d_nodes::Matrix{T},
                           kernel::RK) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
    n = size(nodes, 1)
    n_1 = size(nodes, 2)
    n_2 = size(d_nodes, 2)

    min_bound = Vector{T}(undef, n)
    compression::T = 0
    @inbounds for i = 1:n
        min_bound[i] = nodes[i,1]
        maxx::T = nodes[i,1]
        for j = 2:n_1
            min_bound[i] = min(min_bound[i], nodes[i,j])
            maxx = max(maxx, nodes[i,j])
        end
        for j = 1:n_2
            min_bound[i] = min(min_bound[i], d_nodes[i,j])
            maxx = max(maxx, d_nodes[i,j])
        end
        compression = max(compression, maxx - min_bound[i])
    end

    if compression <= eps(T(1.0))
        error("Cannot estimate_epsilon: `nodes`, `d_nodes` data are not correct.")
    end

    t_nodes = similar(nodes)
    t_d_nodes = similar(d_nodes)
    @inbounds for i = 1:n
        for j = 1:n_1
            t_nodes[i,j] = (nodes[i,j] - min_bound[i]) / compression
        end
        for j = 1:n_2
            t_d_nodes[i,j] = (d_nodes[i,j] - min_bound[i]) / compression
        end
    end

    ε = _estimate_ε(t_nodes, t_d_nodes)
    if isa(kernel, RK_H1)
        ε *= T(2.0)
    elseif isa(kernel, RK_H2)
        ε *= T(2.5)
    end

    return ε
end

function _get_gram(nodes::Matrix{T},
                   kernel::RK
                  ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    n = size(nodes, 1)
    n_1 = size(nodes, 2)
    min_bound = Vector{T}(undef, n)
    compression::T = 0
    @inbounds for i = 1:n
        min_bound[i] = nodes[i,1]
        maxx::T = nodes[i,1]
        for j = 2:n_1
            min_bound[i] = min(min_bound[i], nodes[i,j])
            maxx = max(maxx, nodes[i,j])
        end
        compression = max(compression, maxx - min_bound[i])
    end

    if compression <= eps(T(1.0))
        error("Cannot prepare the spline: `nodes` data are not correct.")
    end

    t_nodes = similar(nodes)
    @inbounds for i = 1:n
        for j = 1:n_1
            t_nodes[i,j] = (nodes[i,j] - min_bound[i]) / compression
        end
    end

    if T(kernel.ε) == T(0.0)
        ε = _estimate_ε(t_nodes)
        if isa(kernel, RK_H0)
            kernel = RK_H0(ε)
        elseif isa(kernel, RK_H1)
            ε *= T(1.5)
            kernel = RK_H1(ε)
        elseif isa(kernel, RK_H2)
            ε *= T(2.0)
            kernel = RK_H2(ε)
        else
            error("incorrect `kernel` type.")
        end
    end

    return _gram(t_nodes, kernel)
end

function _get_gram(nodes::Matrix{T},
                   d_nodes::Matrix{T},
                   es::Matrix{T},
                   kernel::RK
                  ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
    n = size(nodes, 1)
    n_1 = size(nodes, 2)
    n_2 = size(d_nodes, 2)

    if(size(es, 2) != n_2)
        error("Number of derivative directions does not correspond to the number of derivative nodes.")
    end

    t_es = similar(es)
    try
        @inbounds for i = 1:n_2
            t_es[:,i] = es[:,i] ./ norm(es[:,i])
        end
    catch
        error("Cannot normalize derivative direction: zero direction vector.")
    end

    min_bound = Vector{T}(undef, n)
    compression::T = 0
    @inbounds for i = 1:n
        min_bound[i] = nodes[i,1]
        maxx::T = nodes[i,1]
        for j = 2:n_1
            min_bound[i] = min(min_bound[i], nodes[i,j])
            maxx = max(maxx, nodes[i,j])
        end
        for j = 1:n_2
            min_bound[i] = min(min_bound[i], d_nodes[i,j])
            maxx = max(maxx, d_nodes[i,j])
        end
        compression = max(compression, maxx - min_bound[i])
    end

    if compression <= eps(T(1.0))
        error("Cannot prepare the spline: `nodes`, `d_nodes` data are not correct.")
    end

    t_nodes = similar(nodes)
    t_d_nodes = similar(d_nodes)
    @inbounds for i = 1:n
        for j = 1:n_1
            t_nodes[i,j] = (nodes[i,j] - min_bound[i]) / compression
        end
        for j = 1:n_2
            t_d_nodes[i,j] = (d_nodes[i,j] - min_bound[i]) / compression
        end
    end

    if T(kernel.ε) == T(0.0)
        ε = _estimate_ε(t_nodes, t_d_nodes)
        if isa(kernel, RK_H1)
            ε *= T(2.0)
            kernel = RK_H1(ε)
        elseif isa(kernel, RK_H2)
            ε *= T(2.5)
            kernel = RK_H2(ε)
        else
            error("incorrect `kernel` type.")
        end
    end

    return _gram(t_nodes, t_d_nodes, t_es, kernel)
end

function _get_cond(nodes::Matrix{T},
                   kernel::RK
                  ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
        cond = T(0.0)
        gram = _get_gram(nodes, kernel)
        try
            evs = svdvals!(gram)
            maxevs = maximum(evs)
            minevs = minimum(evs)
            if minevs > T(0.0)
               cond = maxevs / minevs
               cond = 10.0^floor(log10(cond))
            end
        catch
        end
        return cond
end

function _get_cond(nodes::Matrix{T},
                   d_nodes::Matrix{T},
                   es::Matrix{T},
                   kernel::RK
                  ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
        cond = T(0.0)
        gram = _get_gram(nodes, kernel)
        try
            evs = svdvals!(gram)
            maxevs = maximum(evs)
            minevs = minimum(evs)
            if minevs > T(0.0)
               cond = maxevs / minevs
               cond = 10.0^floor(log10(cond))
            end
        catch
        end
        return cond
end

# ```
# Get estimation of the Gram matrix condition number
# Brás, C.P., Hager, W.W. & Júdice, J.J. An investigation of feasible descent algorithms for estimating the condition number of a matrix. TOP 20, 791–809 (2012).
# https://link.springer.com/article/10.1007/s11750-010-0161-9
# ```
function _estimate_cond(gram::Matrix{T},
                        chol::LinearAlgebra.Cholesky{T,Array{T,2}},
                        nit = 3
                       ) where T <: AbstractFloat
    if isnothing(gram)
        throw(DomainError(gram, "Parameter `gram` is `nothing'."))
    end
    if isnothing(chol)
        throw(DomainError(chol, "Parameter `chol` is `nothing'."))
    end
    mat_norm = norm(gram, 1)
    n = size(gram, 1)
    x = Vector{T}(undef, n)
    @. x = T(1.0) / T(n)
    z = Vector{T}(undef, n)
    gamma = T(0.0)
    for it = 1:nit
        z = ldiv!(z, chol, x)
        gamma = T(0.0);
        for i = 1:n
            gamma += abs(z[i])
            z[i] = sign(z[i])
        end
        z = ldiv!(z, chol, copy(z))
        zx = z ⋅ x
        idx = 1
        for i = 1:n
            z[i] = abs(z[i])
            if z[i] > z[idx]
                idx = i
            end
        end
        if z[idx] <= zx
            break
        end
        @. x = T(0.0)
        x[idx] = T(1.0)
    end
    cond = T(10.0)^floor(log10(mat_norm * gamma))
    return cond
end
