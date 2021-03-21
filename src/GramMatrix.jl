function _gram(nodes::Matrix{T},
               kernel::RK
              ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    n_1 = size(nodes, 2)
    mat = Matrix{T}(undef, n_1, n_1)
    @inbounds for l = 1:n_1
        for i = 1:l
            mat[i,l] = _rk(kernel, nodes[:,i], nodes[:,l])
            mat[l,i] = mat[i,l]
        end
    end
    return mat
end

function _gram(nodes::Matrix{T},
               d_nodes::Matrix{T},
               es::Matrix{T},
               kernel::RK
              ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
    if !(isa(kernel, RK_H1) || isa(kernel, RK_H2))
        error("incorrect `kernel` type.")
    end

    n = size(nodes, 1)
    m_1 = size(nodes, 2)
    m_2 = size(d_nodes, 2)
    m = m_1 + m_2
    mat = Matrix{T}(undef, m, m)
    @inbounds for j = 1:m_1
              for i = 1:j
                  mat[i,j] = _rk(kernel, nodes[:,i], nodes[:,j])
                  mat[j,i] = mat[i,j]
              end
    end
    m_1_p1 = m_1 + 1
    @inbounds for j = m_1_p1:m
                  j1 = j - m_1
                  for i = 1:m_1
                      mat[i,j] = _∂rk_∂e(kernel, nodes[:,i], d_nodes[:,j1], es[:,j1])
                      mat[j,i] = mat[i,j]
                  end
    end
    ε2 = kernel.ε^2
    @inbounds for j = m_1_p1:m
        j1 = j - m_1
        for i = j:m
            if i == j
                mat[j,j] = ε2
                continue
            end
            i1 = i - m_1
            s::T = T(0.0)
            for r = 1:n
                for k = 1:n
                    s += _∂²rk_∂η_r_∂ξ_k(kernel, d_nodes[:,j1], d_nodes[:,i1], r, k) * es[:,j1][k] * es[:,i1][r]
                end
            end
            mat[j,i] = s
            mat[i,j] = s
        end
    end
    return mat
end
