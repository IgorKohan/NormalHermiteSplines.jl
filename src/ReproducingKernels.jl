@doc raw"
`struct RK_H0{T <: AbstractFloat} <: ReproducingKernel_0`

Defines a type of reproducing kernel of Bessel Potential space ``H^{n/2 + 1/2}_ε (R^n)`` ('Basic Matérn kernel'):
```math
V(\eta , \xi, \varepsilon) = \exp (-\varepsilon |\xi - \eta|) \, .
```
# Fields
- `ε::T`: 'scaling parameter' from the Bessel Potential space definition,
           it may be omitted in the struct constructor otherwise it must be greater than zero
"
struct RK_H0{T <: AbstractFloat} <: ReproducingKernel_0
     ε::T
     RK_H0() = new{Float64}(0.0)
     function RK_H0(ε::T) where T <: AbstractFloat
        if ε <= 0
          throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{T}(ε)
     end
     function RK_H0(ε::Integer)
        if ε <= 0
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{Float64}(convert(Float64, ε))
     end
end

@doc raw"
`struct RK_H1{T <: AbstractFloat} <: ReproducingKernel_1`

Defines a type of reproducing kernel of Bessel Potential space ``H^{n/2 + 3/2}_ε (R^n)`` ('Linear Matérn kernel'):
```math
V(\eta , \xi, \varepsilon) = \exp (-\varepsilon |\xi - \eta|)
             (1 + \varepsilon |\xi  - \eta|) \, .
```
# Fields
- `ε::T`: 'scaling parameter' from the Bessel Potential space definition,
           it may be omitted in the struct constructor otherwise it must be greater than zero
"
struct RK_H1{T <: AbstractFloat} <: ReproducingKernel_1
     ε::T
     RK_H1() = new{Float64}(0.0)
     function RK_H1(ε::T) where T <: AbstractFloat
        if ε <= T(0.0)
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{T}(ε)
     end
     function RK_H1(ε::Integer)
        if ε <= 0
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{Float64}(convert(Float64, ε))
     end
end

@doc raw"
`struct RK_H2{T <: AbstractFloat} <: ReproducingKernel_2`

Defines a type of reproducing kernel of Bessel Potential space ``H^{n/2 + 5/2}_ε (R^n)`` ('Quadratic Matérn kernel'):
```math
V(\eta , \xi, \varepsilon) = \exp (-\varepsilon |\xi - \eta|)
             (3 + 3\varepsilon |\xi  - \eta| + \varepsilon ^2 |\xi - \eta| ^2 ) \, .
```
# Fields
- `ε::T`: 'scaling parameter' from the Bessel Potential space definition,
           it may be omitted in the struct constructor otherwise it must be greater than zero
"
struct RK_H2{T <: AbstractFloat} <: ReproducingKernel_2
     ε::T
     RK_H2() = new{Float64}(0.0)
     function RK_H2(ε::T) where T <: AbstractFloat
        if ε <= T(0.0)
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{T}(ε)
     end
     function RK_H2(ε::Integer)
        if ε <= 0
           throw(DomainError(ε, "Parameter ε must be a positive number."))
        end
        new{Float64}(convert(Float64, ε))
     end
end

@inline function _rk(kernel::RK,
                     η::Vector{T},
                     ξ::Vector{T}
                    ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
   defined::Bool  = false
   x::T = kernel.ε * norm(ξ .- η)
   if isa(kernel, RK_H2)
      defined = true
      value = (T(3.0) + x * (T(3.0) + x)) * exp(-x)
   end
   if isa(kernel, RK_H1)
      defined = true
      value = (T(1.0) + x) * exp(-x)
   end
   if isa(kernel, RK_H0)
      defined = true
      value = exp(-x)
   end
   if !defined
      Throw(ArgumentError("`kernel`: Incorrect parameter type."))
   end
   return value
end

@inline function _∂rk_∂e(kernel::RK,
                         η::Vector{T},
                         ξ::Vector{T},
                         e::Vector{T},
                        ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
   value::T = T(0.0)
   defined::Bool  = false
   t = η .- ξ
   x = kernel.ε * norm(t)
   s = sum(t .* e)
   if isa(kernel, RK_H2)
      defined = true
      value = kernel.ε^2 * exp(-x) * (T(1.0) + x) * s
   end
   if isa(kernel, RK_H1)
      defined = true
      value = kernel.ε^2 * exp(-x) * s
   end
   if !defined
      Throw(ArgumentError("`kernel`: Incorrect parameter type."))
   end
   return value
end

@inline function _∂rk_∂η_k(kernel::RK,
                           η::Vector{T},
                           ξ::Vector{T},
                           k::Int
                          ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
#  Note: Derivative of spline built with reproducing kernel RK_H0 does not exist at the spline nodes.
   value::T = T(0.0)
   defined::Bool  = false
   normt = norm(η .- ξ)
   x = kernel.ε * normt
   if isa(kernel, RK_H2)
      defined = true
      value = kernel.ε^2 * exp(-x) * (T(1.0) + x) * (ξ[k] - η[k])
   end
   if isa(kernel, RK_H1)
      defined = true
      value = kernel.ε^2 * exp(-x) * (ξ[k] - η[k])
   end
   if isa(kernel, RK_H0)
      if normt < sqrt(eps(T))
         Throw(ArgumentError("Derivative does not exist at this point."))
      end
      defined = true
      value = kernel.ε * exp(-x) * (ξ[k] - η[k]) / normt
   end
   if !defined
      Throw(ArgumentError("`kernel`: Incorrect parameter type."))
   end
   return value
end

@inline function _∂²rk_∂η_r_∂ξ_k(kernel::RK,
                                η::Vector{T},
                                ξ::Vector{T},
                                r::Int,
                                k::Int
                               ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
   value::T = T(0.0)
   defined::Bool  = false
   if isa(kernel, RK_H2)
      defined = true
      x = kernel.ε * norm(η .- ξ)
      if r == k
         if x > T(0.0)
            value = kernel.ε^2 * exp(-x) * (T(1.0) + x - (kernel.ε * (ξ[r] - η[r]))^2)
         else
            value = kernel.ε^2
         end
      else
         if x > T(0.0)
            value = -kernel.ε^4 * exp(-x) * (ξ[r] - η[r]) * (ξ[k] - η[k])
         else
            value = T(0.0)
         end
      end
   end
   if isa(kernel, RK_H1)
      defined = true
      t = norm(η .- ξ)
      x = kernel.ε * t
      if r == k
         if t > T(0.0)
            value = kernel.ε^2 * exp(-x) * (T(1.0) - kernel.ε * (ξ[r] - η[r])^2 / t)
         else
            value = kernel.ε^2
         end
      else
         if t > T(0.0)
            value = -kernel.ε^3 * exp(-x) * (ξ[r] - η[r]) * (ξ[k] - η[k]) / t
         else
            value = T(0.0)
         end
      end
   end
   if !defined
      Throw(ArgumentError("`kernel`: Incorrect parameter type."))
   end
   return value
end
