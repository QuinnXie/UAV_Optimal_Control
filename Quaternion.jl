function ToQuaternion(ϕ, θ, ψ)
    # Abbreviations for the various angular functions
    cϕ = cos(ϕ * 0.5) 
    sϕ = sin(ϕ * 0.5) 
    cθ = cos(θ * 0.5) 
    sθ = sin(θ * 0.5) 
    cψ = cos(ψ * 0.5) 
    sψ = sin(ψ * 0.5)

    w = cψ * cθ * cϕ + sψ * sθ * sϕ 
    x = cψ * cθ * sϕ - sψ * sθ * cϕ 
    y = sψ * cθ * sϕ + cψ * sθ * cϕ 
    z = sψ * cθ * cϕ - cψ * sθ * sϕ 

    q = [w, x, y, z]    
    return q
end