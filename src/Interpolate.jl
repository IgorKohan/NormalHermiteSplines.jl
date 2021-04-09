function _prepare(nodes::Matrix{T},
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
    @inbounds for j = 1:n_1
        for i = 1:n
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

    gram = _gram(t_nodes, kernel)
    chol = nothing
    try
        chol = cholesky(gram)
    catch
        error("Cannot prepare the spline: Gram matrix is degenerate.")
    end

    cond = _estimate_cond(gram, chol)

    spline = NormalSpline(kernel,
                          compression,
                          t_nodes,
                          nothing,
                          nothing,
                          nothing,
                          nothing,
                          min_bound,
                          gram,
                          chol,
                          nothing,
                          cond
                         )
    return spline
end

function _construct(spline::NormalSpline{T, RK},
                    values::Vector{T},
                    cleanup::Bool = false
                   ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    if(length(values) != size(spline._nodes, 2))
        error("Number of data values does not correspond to the number of nodes.")
    end
    if isnothing(spline._chol)
        error("Gram matrix was not factorized.")
    end

    mu = Vector{T}(undef, size(spline._gram, 1))
    ldiv!(mu, spline._chol, values)

    spline = NormalSpline(spline._kernel,
                          spline._compression,
                          spline._nodes,
                          values,
                          nothing,
                          nothing,
                          nothing,
                          spline._min_bound,
                          cleanup ? nothing : spline._gram,
                          cleanup ? nothing : spline._chol,
                          mu,
                          spline._cond
                         )
    return spline
end

###################

function _prepare(nodes::Matrix{T},
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
        error("Cannot prepare the spline: `nodes` data are not correct.")
    end

    t_nodes = similar(nodes)
    t_d_nodes = similar(d_nodes)

    @inbounds for j = 1:n_1
        for i = 1:n
            t_nodes[i,j] = (nodes[i,j] - min_bound[i]) / compression
        end
    end

    @inbounds for j = 1:n_2
        for i = 1:n
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

    gram = _gram(t_nodes, t_d_nodes, t_es, kernel)
    chol = nothing
    try
        chol = cholesky(gram)
    catch
        error("Cannot prepare the spline: Gram matrix is degenerate.")
    end

    cond = _estimate_cond(gram, chol)

    spline = NormalSpline(kernel,
                          compression,
                          t_nodes,
                          nothing,
                          t_d_nodes,
                          t_es,
                          nothing,
                          min_bound,
                          gram,
                          chol,
                          nothing,
                          cond
                         )
    return spline
end

function _construct(spline::NormalSpline{T, RK},
                    values::Vector{T},
                    d_values::Vector{T},
                    cleanup::Bool = false
                   ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    if(length(values) != size(spline._nodes, 2))
        error("Number of data values does not correspond to the number of nodes.")
    end
    if(length(d_values) != size(spline._d_nodes, 2))
        error("Number of derivative values does not correspond to the number of derivative nodes.")
    end
    if isnothing(spline._chol)
        error("Gram matrix was not factorized.")
    end

    mu = Vector{T}(undef, size(spline._gram, 1))
    ldiv!(mu, spline._chol, [values; spline._compression .* d_values])

    spline = NormalSpline(spline._kernel,
                          spline._compression,
                          spline._nodes,
                          values,
                          spline._d_nodes,
                          spline._es,
                          d_values,
                          spline._min_bound,
                          cleanup ? nothing : spline._gram,
                          cleanup ? nothing : spline._chol,
                          mu,
                          spline._cond
                         )
    return spline
end

function _evaluate(spline::NormalSpline{T, RK},
                   points::Matrix{T},
                   do_parallel::Bool = false
                  ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    if isnothing(spline)
        error("Spline was not prepared.")
    end
    if isnothing(spline._mu)
        error("Spline coefficients were not calculated.")
    end

    if size(points, 1) != size(spline._nodes, 1)
        if size(points, 1) == 1 && size(points, 2) > 1
            error("Incorrect first dimension of the `points` parameter (use 'evaluate_one' function for evaluating the spline at one point).")
        else
            error("Incorrect first dimension of the `points` parameter (the spline was built in the space of different dimension).")
        end
    end

    n = size(spline._nodes, 1)
    n_1 = size(spline._nodes, 2)
    m = size(points, 2)

    pts = similar(points)
    @inbounds for j = 1:m
        for i = 1:n
            pts[i,j] = (points[i,j] - spline._min_bound[i]) / spline._compression
        end
    end

    spline_values = Vector{T}(undef, m)
    iend1 = iend2 = iend3 = 0
    istart2 = istart3 = istart4 = 0
    mu = spline._mu[1:n_1]
    d_mu = spline._mu[(n_1 + 1):end]
    if do_parallel && m >= 1000 && Threads.nthreads() >= 4
        step = m ÷ 4
        iend1 = 1 + step
        istart2 = iend1 + 1
        iend2 = istart2 + step
        istart3 = iend2 + 1
        iend3 = (istart3 + step) < m ? (istart3 + step) : m
        istart4 = iend3 + 1
        @inbounds Threads.@threads for t = 1:4
            if t == 1
                _do_work(1, iend1, pts, spline._nodes, mu, spline._kernel, spline_values)
                if !isnothing(spline._d_nodes)
                    _do_work_d(1, iend1, pts, spline._d_nodes, spline._es, d_mu, spline._kernel, spline_values)
                end
            elseif t == 2
                _do_work(istart2, iend2, pts, spline._nodes, mu, spline._kernel, spline_values)
                if !isnothing(spline._d_nodes)
                    _do_work_d(istart2, iend2, pts, spline._d_nodes, spline._es, d_mu, spline._kernel, spline_values)
                end
            elseif t == 3
                _do_work(istart3, iend3, pts, spline._nodes, mu, spline._kernel, spline_values)
                if !isnothing(spline._d_nodes)
                    _do_work_d(istart3, iend3, pts, spline._d_nodes, spline._es, d_mu, spline._kernel, spline_values)
                end
            elseif t == 4
                if istart4 <= m
                    _do_work(istart4, m, pts, spline._nodes, mu, spline._kernel, spline_values)
                    if !isnothing(spline._d_nodes)
                        if istart4 <= m
                            _do_work_d(istart4, m, pts, spline._d_nodes, spline._es, d_mu, spline._kernel, spline_values)
                        end
                    end
                end
            end
        end
    else
        _do_work(1, m, pts, spline._nodes, mu, spline._kernel, spline_values)
        if !isnothing(spline._d_nodes)
            _do_work_d(1, m, pts, spline._d_nodes, spline._es, d_mu, spline._kernel, spline_values)
        end
    end

    return spline_values
end

function _do_work(istart::Int,
                  iend::Int,
                  points::Matrix{T},
                  nodes::Matrix{T},
                  mu::Vector{T},
                  kernel::RK,
                  spline_values::Vector{T}
                ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    n_1 = size(nodes, 2)
    h_values = Vector{T}(undef, n_1)
    @inbounds for p = istart:iend
        for i = 1:n_1
            h_values[i] = _rk(kernel, points[:,p], nodes[:,i])
        end
        spline_values[p] = sum(mu .* h_values)
    end
end

function _do_work_d(istart::Int,
                    iend::Int,
                    pts::Matrix{T},
                    d_nodes::Matrix{T},
                    es::Matrix{T},
                    d_mu::Vector{T},
                    kernel::RK,
                    spline_values::Vector{T}
                   ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
   n_2 = size(d_nodes, 2)
   d_h_values = Vector{T}(undef, n_2)
   @inbounds for p = istart:iend
       for i = 1:n_2
           d_h_values[i] = _∂rk_∂e(kernel, pts[:,p], d_nodes[:,i], es[:,i])
       end
       spline_values[p] += sum(d_mu .* d_h_values)
    end
end

function _evaluate_gradient(spline::NormalSpline{T, RK},
                            point::Vector{T}
                           ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    if isnothing(spline._mu)
        error("Spline coefficients were not calculated.")
    end

    n = size(spline._nodes, 1)
    n_1 = size(spline._nodes, 2)

    pt = Vector{T}(point)
    @inbounds for i = 1:n
        pt[i] = (point[i] - spline._min_bound[i]) / spline._compression
    end

    d_h_values = Vector{T}(undef, n_1)
    grad = Vector{T}(undef, n)
    mu = spline._mu[1:n_1]
    @inbounds for k = 1:n
        for i = 1:n_1
            d_h_values[i] = _∂rk_∂η_k(spline._kernel, pt, spline._nodes[:,i], k)
        end
        grad[k] = sum(mu .* d_h_values)
    end

    if !isnothing(spline._d_nodes)
        n_2 = size(spline._d_nodes, 2)
        d_h_values = Vector{T}(undef, n_2)
        d_mu = spline._mu[n_1+1:end]
        @inbounds for k = 1:n
            for i = 1:n_2
                d_h_values[i] = T(0.0)
                for l = 1:n
                    d_h_values[i] += (_∂²rk_∂η_r_∂ξ_k(spline._kernel, pt, spline._d_nodes[:,i], k, l) * spline._es[l,i])
                end
            end
            grad[k] += sum(d_mu .* d_h_values)
        end
    end
    @inbounds for k = 1:n
        grad[k] /= spline._compression
    end
    return grad
end
