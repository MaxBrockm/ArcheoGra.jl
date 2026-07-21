"""
    Tobler_Hiking(S::Float64)

Calculates the hiking velocity using Tobler's Hiking function.

The function estimates walking speed based on terrain slope. It is generally 
defined for slopes between -60 and 50 degrees. The formula used is 
`v = 6 * exp(-3.5 * |S + 0.05|)`, where `S` is the slope.

# Arguments
- `S::Float64`: The terrain slope (rise over run), where `S = tan(angle)`.

# Output
- `Float64`: The estimated hiking velocity in km/h (bounded within `(0, 6]`).
"""
function Tobler_Hiking(S::Float64)
    return 6*exp(-3.5*abs(S+0.05))
end

"""
    KS_Hiking(S::Float64)

Calculates the hiking velocity based on the Kondo-Seino model.

This piecewise exponential model is derived from GPS-aided walking experiments 
on the historical Nakasendō-Kisoji road in Japan (Kondo & Seino, 2010). It 
differentiates between uphill/gentle downhill and steeper downhill slopes.

# Arguments
- `S::Float64`: The terrain slope (rise over run).

# Output
- `Float64`: The estimated hiking velocity in km/h.
"""
function KS_Hiking(S::Float64)
    if S >= -0.07
        return 5.1*exp(-2.25*abs(S+0.07))
    else
        return 5.1*exp(-1.5*abs(S+0.07))
    end
end

"""
    adjusted_Hiking(S::Float64)

Computes an adjustment factor for hiking speed as a function of slope `S`.

For slopes `S >= -0.07`, the function uses a steeper logistic decay;
for steeper downhill slopes (`S < -0.07`), a gentler decay is applied.
This provides a relative cost factor describing the hiking speed, required for least-cost path calculations.

# Arguments
- `S::Float64`: Terrain slope (rise over run), positive for uphill, negative for downhill.

# Output
- `Float64`: Adjustment factor (typically between 0 and 1) describing the relative hiking speed.
"""
function adjusted_Hiking(S::Float64)
    if S >= -0.07
        return 1/(1+5.1*exp(-2.25*abs(S+0.07)))
    else
        return 1/(1+5.1*exp(-1.5*abs(S+0.07)))
    end
end

"""
    f_LloberaSluckin(S; cS=15)

Calculates a symmetric travel cost using the Llobera-Sluckin function.

# Arguments
- `S::Number`: The terrain slope (rise over run).

# Keywords
- `cS::Number`: Critical slope parameter (default: `15`).

# Output
- `Float64`: The calculated symmetric travel cost.
"""
function f_LloberaSluckin(S;cS=15) 
    return 1+(S/cS)^2
end

const LloberaSluckin = ConnectionFunction(
    f_LloberaSluckin , ##Least-cost Paths – Some Methodological Issues: 5-1-4; Irmela Herzog
    true
)

const Tobler_Hiking_adjusted = ConnectionFunction(
    adjusted_Hiking,
    false
)
