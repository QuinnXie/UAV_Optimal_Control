# Trapezoidal interpolation: linear spline for dynamics and control
function interpu(fLow, fUpp, tLow, tUpp)
    c0 = fLow - (fUpp-fLow)/(tUpp-tLow) * tLow
    c1 = (fUpp-fLow)/(tUpp-tLow)
    # Interp = t -> (fUpp-fLow)/(tUpp-tLow)*(t-tLow)+fLow
    Interp = Polynomial([c0, c1])
    return Interp
end

# Hermite-Simpson interpolation: quadratic spline for control
function interpu(uLow, uMid, uUpp, tLow, tUpp)
    tk = tUpp-tLow;

    a0 = uLow;
    a1 = (-3*uLow + 4*uMid - uUpp) / tk;
    a2 = (2*uLow - 4*uMid + 2*uUpp) / (tk^2);

    c0 = a2 * tLow^2 - a1 * tLow + a0;
    c1 = -2 * a2 * tLow + a1;
    c2 = a2;

    # interp = t -> c0 + c1*t + c2*t^2
    interp = Polynomial([c0, c1, c2])
    return interp
end

# Trapezoidal interpolation: quadratic spline for state
function interpx(xLow, fLow, fUpp, tLow, tUpp)
    # c3 = 0.5 * (value.(δq₁)[i+1] - value.(δq₁)[i]) / (ts[i+1] - ts[i])
    # c1 = value.(q₁)[i] - value.(δq₁)[i] * ts[i] + c3 * ts[i]^2
    # c2 = value.(δq₁)[i] - 2 * c3 * ts[i]
    # f_position[i] = Polynomial([c1, c2, c3])

    c3 = 0.5 * (fUpp - fLow) / (tUpp - tLow)
    c1 = xLow - fLow * tLow + c3 * tLow^2
    c2 = fLow - 2 * c3 * tLow

    # interp = t -> c0 + c1*t + c2*t^2
    Interp = Polynomial([c1, c2, c3])
    return Interp
end

# Hermite-Simpson interpolation: cubic spline for state
function interpx(xLow, fLow, fMid, fUpp, tLow, tUpp)
    # tk = ts[2*i+1]-ts[2*i-1];

    # a0 = value.(q₁)[2*i-1];
    # a1 = value.(δq₁)[2*i-1];
    # a2 = 0.5 * (-3*value.(δq₁)[2*i-1] + 4*value.(δq₁)[2*i] - value.(δq₁)[2*i+1]) / tk;
    # a3 = (2*value.(δq₁)[2*i-1] - 4*value.(δq₁)[2*i] + 2*value.(δq₁)[2*i+1]) / (3*tk^2);

    # c0 = -a3*ts[2*i-1]^3 + a2*ts[2*i-1]^2 - a1*ts[2*i-1] + a0;
    # c1 = 3*a3*ts[2*i-1]^2 - 2*a2*ts[2*i-1] + a1;
    # c2 = -3*a3*ts[2*i-1] + a2;
    # c3 = a3;

    # f_position[i] = Polynomial([c0, c1, c2, c3])

    tk = tUpp-tLow;

    a0 = xLow;
    a1 = fLow;
    a2 = 0.5 * (-3*fLow + 4*fMid - fUpp) / tk;
    a3 = (2*fLow - 4*fMid + 2*fUpp) / (3*tk^2);

    c0 = -a3*tLow^3 + a2*tLow^2 - a1*tLow + a0;
    c1 = 3*a3*tLow^2 - 2*a2*tLow + a1;
    c2 = -3*a3*tLow + a2;
    c3 = a3;

    # interp = t -> c0 + c1*t + c2*t^2 + c3*t^3
    Interp = Polynomial([c0, c1, c2, c3])
    return Interp
end