### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 786426e2-9713-499b-b3a6-f1abdc1a8c2a
using Pkg

# ╔═╡ 35443bfa-2604-4eaa-954b-b57d08090383
# ╠═╡ show_logs = false
Pkg.add("DifferentialEquations")

# ╔═╡ bce1b73f-20b5-429c-b282-91c3287ced1a
# ╠═╡ show_logs = false
Pkg.add("Plots")

# ╔═╡ 6961e14d-bdc5-431e-a8c8-3bece1babb84
using DifferentialEquations

# ╔═╡ 424e42d9-c6c6-471d-ad8c-e75709b68e29
using Plots

# ╔═╡ 80da5392-1db7-486a-aa55-ed6b152d09f5
g = 9.81

# ╔═╡ e15fd4dc-5cbb-48bf-9e96-5e3f884d3b9d
h1 = 0.22

# ╔═╡ 34c00ea3-b6dd-4444-b29b-8ab5636b6a66
w1 = 0.1

# ╔═╡ 227d8898-48e7-4647-8f84-cb0eb9467e27
L = 0.15

# ╔═╡ e4f6111d-a490-4ffd-bfa4-08a6cf26d24d
m = 0.1

# ╔═╡ 588daec4-11b1-43a4-b23f-3a610cf7cfde
function spinning_pendulum!(du, u, p, t)
    Ω = p
    θ, θdot = u
    du[1] = θdot
    du[2] = (Ω^2 / L) * (w1 + L*sin(θ)) * cos(θ) - (g / L) * sin(θ)
end

# ╔═╡ c7da7154-5724-4a0f-aae4-4996a5fd9f5c
function position(θ, t, Ω)
    r = w1 + L*sin(θ)
    x = r*cos(Ω*t)
    y = r*sin(Ω*t)
    z = h1 - L*cos(θ)

    return x, y, z
end

# ╔═╡ a47bbb6a-d9c4-48cc-8ac0-039b4f838b4d
θ0 = 0.05

# ╔═╡ a7451fa8-7aa2-4ac4-8e89-685c7abad1ae
θdot0 = 0.0

# ╔═╡ 0275c055-35d7-4755-b601-2f6fd97d119d
u0 = [θ0, θdot0]

# ╔═╡ d271d2af-46cb-4f4a-9e4e-9ff884757bf2
tspan = (0.0, 8.0)

# ╔═╡ 009e2b4e-9f69-438b-a3be-aa966f1a3ea9
Ω_slow = 2.0

# ╔═╡ ef90a3e6-0e7e-45e0-9675-eb1f31b06930
Ω_fast = 12.0

# ╔═╡ 39d73d33-62ee-4186-95bb-24ab6319f87a
prob_slow = ODEProblem(spinning_pendulum!, u0, tspan, Ω_slow)

# ╔═╡ b7649c24-8cb2-4f4d-9bc3-aea576b8f38a
prob_fast = ODEProblem(spinning_pendulum!, u0, tspan, Ω_fast)

# ╔═╡ 8c38f861-3f37-4740-ab64-4d7586737d1c
sol_slow = solve(prob_slow, Tsit5(), reltol=1e-9, abstol=1e-9)

# ╔═╡ 6b0985bf-48f0-4308-8ef0-e5f55fe96dfd
sol_fast = solve(prob_fast, Tsit5(), reltol=1e-9, abstol=1e-9)

# ╔═╡ 9dcbd25e-693a-4d1f-8bf2-24192f4ae942
p1 = plot(sol_slow.t, getindex.(sol_slow.u,1),
    xlabel="t (s)", ylabel="θ (rad)",
    label="slow Ω = $Ω_slow rad/s", lw=2)

# ╔═╡ a254bc18-7a5d-417f-8d33-f2ccf7da4cd1
p2 = plot(sol_fast.t, getindex.(sol_fast.u,1),
    xlabel="t (s)", ylabel="θ (rad)",
    label="fast Ω = $Ω_fast rad/s", lw=2)

# ╔═╡ 8ee6c502-fb67-4681-aa7d-f10714995221
plot(p1, p2, layout=(2,1), size=(800,600))

# ╔═╡ 022d3257-ae12-45fd-89d0-e6e0ec788405
xs = Float64[]

# ╔═╡ 5a1f228b-34c7-4aac-8af2-d7f435054627
ys = Float64[]

# ╔═╡ 47f1e26e-ac49-4aef-a44e-f4ba8745ffa0
zs = Float64[]

# ╔═╡ 8ad69100-f0c9-441d-8033-940702fd2fbd
for (ti, ui) in zip(sol_slow.t, sol_slow.u)
    x, y, z = position(ui[1], ti, Ω_slow)
    push!(xs, x)
    push!(ys, y)
    push!(zs, z)
end

# ╔═╡ 2d8c9c6a-2cde-41a5-bfeb-84849a85ca33
plot(xs, ys, zs,
    lw=2,
    xlabel="x",
    ylabel="y",
    zlabel="z",
    label="slow case",
    title="Mass trajectory")

# ╔═╡ 17f01baf-f74a-4256-b7b1-24928701b560
xslow = Float64[]

# ╔═╡ f7d80cfd-2a79-4462-98d6-35bb69ee0815
yslow = Float64[]

# ╔═╡ 5486d2c2-709c-4d0a-ab65-1534f60246ca
zslow = Float64[]

# ╔═╡ 9452ef53-441f-4e18-9820-261c183624cd
for (ti, ui) in zip(sol_slow.t, sol_slow.u)
    x, y, z = position(ui[1], ti, Ω_slow)
    push!(xslow, x)
    push!(yslow, y)
    push!(zslow, z)
end

# ╔═╡ a18d645e-bf60-4099-8853-8246f30d1b03
anim = @animate for i in eachindex(sol_slow.t)
    plot(
        [0, xslow[i]], [0, yslow[i]], [h1, zslow[i]],
        lw=3,
        xlabel="x", ylabel="y", zlabel="z",
        xlims=(-0.3, 0.3), ylims=(-0.3, 0.3), zlims=(-0.1, 0.3),
        label="pendulum",
        title="Slow rotation"
    )
    scatter!([xslow[i]], [yslow[i]], [zslow[i]], label="mass")
end

# ╔═╡ 20fa462e-12e6-4ee1-904e-a2c8c719e1bf
gif(anim, "slow_pendulum.gif", fps=30)

# ╔═╡ Cell order:
# ╠═786426e2-9713-499b-b3a6-f1abdc1a8c2a
# ╠═35443bfa-2604-4eaa-954b-b57d08090383
# ╠═bce1b73f-20b5-429c-b282-91c3287ced1a
# ╠═6961e14d-bdc5-431e-a8c8-3bece1babb84
# ╠═424e42d9-c6c6-471d-ad8c-e75709b68e29
# ╠═80da5392-1db7-486a-aa55-ed6b152d09f5
# ╠═e15fd4dc-5cbb-48bf-9e96-5e3f884d3b9d
# ╠═34c00ea3-b6dd-4444-b29b-8ab5636b6a66
# ╠═227d8898-48e7-4647-8f84-cb0eb9467e27
# ╠═e4f6111d-a490-4ffd-bfa4-08a6cf26d24d
# ╠═588daec4-11b1-43a4-b23f-3a610cf7cfde
# ╠═c7da7154-5724-4a0f-aae4-4996a5fd9f5c
# ╠═a47bbb6a-d9c4-48cc-8ac0-039b4f838b4d
# ╠═a7451fa8-7aa2-4ac4-8e89-685c7abad1ae
# ╠═0275c055-35d7-4755-b601-2f6fd97d119d
# ╠═d271d2af-46cb-4f4a-9e4e-9ff884757bf2
# ╠═009e2b4e-9f69-438b-a3be-aa966f1a3ea9
# ╠═ef90a3e6-0e7e-45e0-9675-eb1f31b06930
# ╠═39d73d33-62ee-4186-95bb-24ab6319f87a
# ╠═b7649c24-8cb2-4f4d-9bc3-aea576b8f38a
# ╠═8c38f861-3f37-4740-ab64-4d7586737d1c
# ╠═6b0985bf-48f0-4308-8ef0-e5f55fe96dfd
# ╠═9dcbd25e-693a-4d1f-8bf2-24192f4ae942
# ╠═a254bc18-7a5d-417f-8d33-f2ccf7da4cd1
# ╠═8ee6c502-fb67-4681-aa7d-f10714995221
# ╠═022d3257-ae12-45fd-89d0-e6e0ec788405
# ╠═5a1f228b-34c7-4aac-8af2-d7f435054627
# ╠═47f1e26e-ac49-4aef-a44e-f4ba8745ffa0
# ╠═8ad69100-f0c9-441d-8033-940702fd2fbd
# ╠═2d8c9c6a-2cde-41a5-bfeb-84849a85ca33
# ╠═17f01baf-f74a-4256-b7b1-24928701b560
# ╠═f7d80cfd-2a79-4462-98d6-35bb69ee0815
# ╠═5486d2c2-709c-4d0a-ab65-1534f60246ca
# ╠═9452ef53-441f-4e18-9820-261c183624cd
# ╠═a18d645e-bf60-4099-8853-8246f30d1b03
# ╠═20fa462e-12e6-4ee1-904e-a2c8c719e1bf
