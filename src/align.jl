export Align

"""
    Align(inner_layout :: AbstractLayout{2}, angle :: Real = zero(Float64))
    
Align the vertex positions of `inner_layout` so that the principal axis of the resulting
layout makes an `angle` with the **x**-axis.

Only supports two-dimensional inner layouts.
"""
struct Align{Ptype, L <: AbstractLayout{2, Ptype}} <: AbstractLayout{2, Ptype}
    inner_layout :: L
    angle :: Ptype
    function Align(inner_layout::L, angle::Real) where {L <: AbstractLayout{2, Ptype}} where Ptype
        new{Ptype, L}(inner_layout, convert(Ptype, angle))
    end
end
Align(inner_layout::AbstractLayout{2, Ptype}) where Ptype = Align(inner_layout, zero(Ptype))

function layout(algo::Align{Ptype, <:AbstractLayout{2, Ptype}}, adj_matrix::AbstractMatrix) where {Ptype}
    # compute "inner" layout
    rs = layout(algo.inner_layout, adj_matrix)

    # align the "inner" layout to have its principal axis make `algo.angle` with x-axis
    # step 1: compute covariance matrix for PCA analysis: 
    #       C = ∑ᵢ (rᵢ - ⟨r⟩) (rᵢ - ⟨r⟩)ᵀ
    # for vertex positions rᵢ, i = 1, …, N, and center of mass ⟨r⟩ = N⁻¹ ∑ᵢ rᵢ.
    centerofmass = sum(rs) / length(rs)
    C = zeros(SMatrix{2, 2, Ptype})
    for r in rs
        C += (r - centerofmass) * (r - centerofmass)'
    end
    vs = eigen(C).vectors

    # step 2: pick principal axis (largest eigenvalue → last eigenvalue/vector)
    axis = vs[:, end]
    axis_angle = atan(axis[2], axis[1])

    # step 3: rotate positions `rs` so that new axis is aligned with `algo.angle`
    s, c = sincos(-axis_angle + algo.angle)
    R = @SMatrix [c -s; s c] # [cos(θ) -sin(θ); sin(θ) cos(θ)]
    for (i, r) in enumerate(rs)
        rs[i] = Point2{Ptype}(R * r) :: Point2{Ptype}
    end

    return rs
end