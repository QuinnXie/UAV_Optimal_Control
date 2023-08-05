using JuMP
import Interpolations
import Ipopt
using Polynomials
using Plots

include("interp.jl")
include("Quaternion.jl")

## Global variables

## Physical parameters for the cart-pole example (Standard units)
const L = 18.0        # length of map
const W = 4.0         # width of map
const H = 10.0        # height of map

const l = 1.0         # length of table
const w = 1.6         # width of table
const T₁ = 0.0        # table right below corner x - position
const T₂ = 0.0        # table right below corner y - position

# UAV parameters
const m = 1.540820046
const Ixx = 0.0053423608
const Iyy = 0.0039935015
const Izz = 0.0452023222
const r₁ = 0.3       # radius of UAV

const g = 9.8

# # initial position
# const A₁ = -8.5
# const A₂ = 1.0
# const A₃ = 1.5

# # final position
# const B₁ = 8.0
# const B₂ = 1.0
# const B₃ = 1.5

# initial position
const A₁ = -8.5
const A₂ = 1.0
const A₃ = 0.0

# final position
const B₁ = -8.5
const B₂ = 1.0
const B₃ = 1.5

const γ = 1.0
## Initial conditions
const q₁_s = A₁       # x - position of UAV
const q₂_s = A₂       # y - position of UAV
const q₃_s = A₃       # z - position of UAV
const q₄_s = 0.0      # ϕ x - orientation of UAV
const q₅_s = 0.0      # θ y - orientation of UAV
const q₆_s = 0.0      # ψ z - orientation of UAV

const δq₁_s = 0.0     # v - x - velocity of UAV
const δq₂_s = 0.0     # v - y - velocity of UAV
const δq₃_s = 0.0     # v - z - velocity of UAV

u₁_s = 0.0
u₂_s = 0.0
u₃_s = 0.0
u₄_s = 0.0

## State Limit
const x₁_low = -L/2 + r₁
const x₁_upp = L/2 - r₁
const y₁_low = -W/2 + r₁
const y₁_upp = W/2 - r₁
const z₁_low = 0.0
const z₁_upp = H

const v_max = 20.0           # maximum velocity
const w_max = 5.0            # maximum angular velocity
const T_max = 19.35          # maximum acceleration
const M_max = 10.0           # maximum angular acceleration

## Final conditions, the so-called Terminal Area Energy Management (TAEM)
const q₁_t = B₁       # x - position of UAV
const q₂_t = B₂       # y - position of UAV
const q₃_t = B₃       # z - position of UAV
const q₄_t = 0.0      # ϕ x - orientation of UAV
const q₅_t = 0.0      # θ y - orientation of UAV
const q₆_t = 0.0      # ψ z - orientation of UAV

const δq₁_t = 0.0     # v - x - velocity of UAV
const δq₂_t = 0.0     # v - y - velocity of UAV
const δq₃_t = 0.0     # v - z - velocity of UAV

u₁_t = 0.0
u₂_t = 0.0
u₃_t = 0.0
u₄_t = 0.0
## Number of mesh points (knots) to be used
N = 200
n = N+1

Time = 10.0
t_s = Time/N

## Integration scheme to be used for the dynamics
# const integration_rule = "rectangular" 
const integration_rule = "trapezoidal"

# !!! tip "Choose a good linear solver"
#     Picking a good linear solver is **extremely important**
#     to maximize the performance of nonlinear solvers.
#     For the best results, it is advised to experiment different linear solvers.
#
#     For example, the linear solver `MA27` is outdated and can be quite slow.
#     `MA57` is a much better alternative, especially for highly-sparse problems
#     (such as trajectory optimization problems).

## Uncomment the lines below to pass user options to the solver
user_options = (
## "mu_strategy" => "monotone",
## "linear_solver" => "ma27",
)

## Create JuMP model, using Ipopt as the solver
model = Model(optimizer_with_attributes(Ipopt.Optimizer, user_options...))

@variables(model, begin
    x₁_low ≤ q₁[1:n] ≤ x₁_upp
    y₁_low ≤ q₂[1:n] ≤ y₁_upp
    z₁_low ≤ q₃[1:n] ≤ z₁_upp
    q₄[1:n]
    q₅[1:n]
    q₆[1:n]
    -v_max ≤ δq₁[1:n] ≤ v_max
    -v_max ≤ δq₂[1:n] ≤ v_max
    -v_max ≤ δq₃[1:n] ≤ v_max
    0 ≤ u₁[1:n] ≤ T_max
    -w_max ≤ u₂[1:n] ≤ w_max
    -w_max ≤ u₃[1:n] ≤ w_max
    -w_max ≤ u₄[1:n] ≤ w_max
    # Δt[1:n] ≥ 0
    ##        3.5 ≤       Δt[1:n] ≤ 4.5     # time step (sec)
    Δt[1:n] == t_s                          # time step (sec)
end)

## Fix initial conditions
fix(q₁[1], q₁_s; force = true)
fix(q₂[1], q₂_s; force = true)
fix(q₃[1], q₃_s; force = true)
fix(q₄[1], q₄_s; force = true)
fix(q₅[1], q₅_s; force = true)
fix(q₆[1], q₆_s; force = true)

fix(δq₁[1], δq₁_s; force = true)
fix(δq₂[1], δq₂_s; force = true)
fix(δq₃[1], δq₃_s; force = true)

## Fix final conditions
fix(q₁[n], q₁_t; force = true)
fix(q₂[n], q₂_t; force = true)
fix(q₃[n], q₃_t; force = true)
fix(q₄[n], q₄_t; force = true)
fix(q₅[n], q₅_t; force = true)
fix(q₆[n], q₆_t; force = true)

fix(δq₁[n], δq₁_t; force = true)
fix(δq₂[n], δq₂_t; force = true)
fix(δq₃[n], δq₃_t; force = true)

# ## Initial guess: linear interpolation
x_s = [q₁_s, q₂_s, q₃_s, q₄_s, q₅_s, q₆_s, δq₁_s, δq₂_s, δq₃_s, u₁_s, u₂_s, u₃_s, u₄_s, t_s]
x_t = [q₁_t, q₂_t, q₃_t, q₄_t, q₅_t, q₆_t, δq₁_t, δq₂_t, δq₃_t, u₁_t, u₂_t, u₃_t, u₄_t, t_s]
interp_linear = Interpolations.LinearInterpolation([1, n], [x_s, x_t])
initial_guess = mapreduce(transpose, vcat, interp_linear.(1:n))
# initial_guess[:,2] = [LinRange(1.0, -1.3, 500); LinRange(-1.3, 1.0, 501)]
set_start_value.(all_variables(model), vec(initial_guess))

## System dynamics
# state = [x, y, z, ϕ, θ, ψ, v_x, v_y, v_z]
# input = [T, p, q, r]

# ̇x = v_x
# ̇y = v_y
# ̇z = v_z
# ̇ϕ = p + (q⋅sinϕ + r⋅cosϕ)⋅tanθ
# ̇θ = q⋅conϕ - r⋅sinϕ
# ̇ψ = (q⋅sinϕ + r⋅cosϕ)/cosθ
# ̇v_x = (sinψ⋅sinϕ + cosψ⋅sinθ⋅cosϕ)⋅T/m
# ̇v_y = (-cosψ⋅sinϕ + sinψ⋅sinθ⋅cosϕ)⋅T/m
# ̇v_z = -g + cosθ⋅cosϕ⋅T/m
# ̇p = 1/Ixx⋅L + (Iyy - Izz)/Ixx⋅r⋅q
# ̇q = 1/Iyy⋅M + (Izz - Ixx)/Iyy⋅p⋅r
# ̇r = 1/Izz⋅N + (Ixx - Iyy)/Izz⋅p⋅q

## Motion of the UAV as a differential-algebraic system of equations (DAEs)
@NLexpression(model, δ²q₁[j = 1:n], (sin(q₆[j])*sin(q₄[j]) + cos(q₆[j])*sin(q₅[j])*cos(q₄[j]))*u₁[j]/m) # ̇dynamics for δ²q₁
@NLexpression(model, δ²q₂[j = 1:n], (-cos(q₆[j])*sin(q₄[j]) + sin(q₆[j])*sin(q₅[j])*cos(q₄[j]))*u₁[j]/m) # ̇dynamics for δ²q₂
@NLexpression(model, δ²q₃[j = 1:n], - g + cos(q₅[j])*cos(q₄[j])*u₁[j]/m) # ̇dynamics for δ²q₃
@NLexpression(model, δq₄[j = 1:n], u₂[j] + (u₃[j]*sin(q₄[j]) + u₄[j]*cos(q₄[j]))*tan(q₅[j])) # ̇dynamics for δq₄
@NLexpression(model, δq₅[j = 1:n], u₃[j]*cos(q₄[j]) - u₄[j]*sin(q₄[j])) # ̇dynamics for δq₅
@NLexpression(model, δq₆[j = 1:n], (u₃[j]*sin(q₄[j]) + u₄[j]*cos(q₄[j]))/cos(q₅[j])) # ̇dynamics for δq₆

# time in equal step
# @NLconstraint(model, [j=1:N], Δt[j] == Δt[j+1])


for j in 2:n
    i = j - 1  # index of previous knot
    if integration_rule == "trapezoidal"
        ## Trapezoidal integration
        @NLconstraint(model, q₁[j] == q₁[i] + 0.5 * Δt[i] * (δq₁[j] + δq₁[i]))
        @NLconstraint(model, q₂[j] == q₂[i] + 0.5 * Δt[i] * (δq₂[j] + δq₂[i]))
        @NLconstraint(model, q₃[j] == q₃[i] + 0.5 * Δt[i] * (δq₃[j] + δq₃[i]))
        @NLconstraint(model, q₄[j] == q₄[i] + 0.5 * Δt[i] * (δq₄[j] + δq₄[i]))
        @NLconstraint(model, q₅[j] == q₅[i] + 0.5 * Δt[i] * (δq₅[j] + δq₅[i]))
        @NLconstraint(model, q₆[j] == q₆[i] + 0.5 * Δt[i] * (δq₆[j] + δq₆[i]))

        @NLconstraint(model, δq₁[j] == δq₁[i] + 0.5 * Δt[i] * (δ²q₁[j] + δ²q₁[i]))
        @NLconstraint(model, δq₂[j] == δq₂[i] + 0.5 * Δt[i] * (δ²q₂[j] + δ²q₂[i]))
        @NLconstraint(model, δq₃[j] == δq₃[i] + 0.5 * Δt[i] * (δ²q₃[j] + δ²q₃[i]))
    else
        @error "Unexpected integration rule '$(integration_rule)'"
    end
end


@NLconstraint(model, [j = 1:n], ((q₁[j]-l/2)^2 / ((l/2+r₁)^2*2) + (q₂[j]-w/2)^2 / ((w/2+r₁)^2*2)) ≥ 1.0)

## Objective: minimize control effort
# @objective(model, Min, sum(Δt[j] for j in 1:N))
# @objective(model, Min, sum(Δt[j] for j in 1:N) + 0.05*sum((u₁[j]^2 + u₂[j]^2 + u₃[j]^2 + u₄[j]^2 + u₅[j]^2 + u₆[j]^2) for j in 1:n))
@NLobjective(model, Min, sum((u₁[j]^2 + u₂[j]^2 + u₃[j]^2 + u₄[j]^2) for j in 1:n))

# set_silent(model)  # Hide solver's verbose output
optimize!(model)  # Solve for the control and state
@assert termination_status(model) == LOCALLY_SOLVED

## Show final cross-range of the solution
println(
    "Final total control effort = ",
    round(objective_value(model); digits = 2),
)

Δx = q₁_t .- value.(q₁)
Δy = q₂_t .- value.(q₂)
Δz = q₃_t .- value.(q₃)

dist = sqrt.(Δx.^2 + Δy.^2 + Δz.^2)

q = zeros(n,4)
for i in 1:n 
    q[i,:] = ToQuaternion(value.(q₄)[i], value.(q₅)[i], value.(q₆)[i])
end

ts = cumsum([0;value.(Δt)])[1:end-1]

# position and angle
plt_position = plot(
    value.(q₁), value.(q₂);
    # seriestype=:scatter, markersize=2, 
    xlims=(-L/2, L/2),
    ylims=(-W/2, W/2), 
    label="UAV",
    title = "UAV trajectory",
    ylabel="y (m)", xlabel="x (m)",
    )

plt_position = plot!(plt_position,
    [0.0, l, l, 0.0, 0.0], [0.0, 0.0, w, w, 0.0];
    # seriestype=:scatter, markersize=2, 
    label="obstacle",
    )


# velocity and agular velocity

velocity_robot = sqrt.(value.(δq₁).^2 + value.(δq₂).^2 + value.(δq₃).^2)

plt_velocity = plot(
    ts, value.(δq₁), 
    title="UAV velocity",
    label="v_x", 
    ylabel="v (m/s)", 
    xlabel="time (s)",
    )

plt_velocity = plot!(plt_velocity,
    ts, value.(δq₂), 
    label="v_y", 
)

plt_velocity = plot!(plt_velocity,
    ts, value.(δq₃),  
    label="v_z", 
)

plt_angle = plot(
    ts, value.(q₄), 
    title="UAV angle",
    label="ϕ", 
    ylabel="angle (rad)", 
    xlabel="time (s)",
    )

plt_angle = plot!(plt_angle,
    ts, value.(q₅), 
    label="θ", 
)

plt_angle = plot!(plt_angle,
    ts, value.(q₆),  
    label="ψ", 
)

plt_angular = plot(
    ts, value.(u₂), 
    title="UAV angle",
    label="p", 
    ylabel="angular velocity (rad/s)", 
    xlabel="time (s)",
    )

plt_angular = plot!(plt_angular,
    ts, value.(u₃), 
    label="q", 
)

plt_angular = plot!(plt_angular,
    ts, value.(u₄),  
    label="r", 
)

plt_Thrust = plot(
    ts[1:n], value.(u₁)[1:n], 
    title="UAV Thrust",
    label="T", 
    ylabel="Thrust (N)", 
    xlabel="time (s)",
    )


display(plt_position)
display(plt_velocity)
display(plt_angle)
display(plt_angular)
display(plt_Thrust)

savefig(plt_position, "plt//plt_position_UAV")
savefig(plt_velocity, "plt//plt_velocity_UAV")
savefig(plt_angle, "plt//plt_angle_UAV")
savefig(plt_angular, "plt//plt_angular_UAV")
savefig(plt_Thrust, "plt//plt_Thrust_UAV")

using CSV
using DataFrames

Thrust_scale = value.(u₁) ./ T_max
df = DataFrame(time = ts, 
               x = value.(q₁),
               y = value.(q₂),
               z = value.(q₃),
               phi = value.(q₄),
               theta = value.(q₅),
               psi = value.(q₆),
               v_x = value.(δq₁),
               v_y = value.(δq₂),
               v_z = value.(δq₃),
               phi_dot = value.(δq₄),
               theta_dot = value.(δq₅),
               psi_dot = value.(δq₆),
               a_x = value.(δ²q₁),
               a_y = value.(δ²q₂),
               a_z = value.(δ²q₃),
               Thrust = value.(u₁),
               Thrust_Scale = Thrust_scale,
               p = value.(u₂),
               q = value.(u₃),
               r = value.(u₄),
               dist_x = Δx,
               dist_y = Δy,
               dist_z = Δz,
               distance = dist,
               q_w = q[:,1],
               q_x = q[:,2],
               q_y = q[:,3],
               q_z = q[:,4]
               )

# CSV.write("export_uav_data.csv", df, newline='\n') 
CSV.write("export_uav_data1.csv", df, newline='\n') 