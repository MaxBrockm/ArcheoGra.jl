
"""
    distance(u::AbstractArray{Float64},v::AbstractArray{Float64})

Function to determine Haversine distance between two pairs of Coordinates
"""
function distance(u::AbstractArray{Float64},v::AbstractArray{Float64})
    d2g(x)=deg2rad(x)
    return 2*earth_radius*asin(sqrt(0.5*(1-cos(d2g(v[2])-d2g(u[2]))+
            cos(d2g(u[2]))*cos(d2g(v[2]))*(1-cos(d2g(v[1])-d2g(u[1]))))))
end

######################################################################
######################################################################
