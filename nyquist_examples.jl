### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

macro bind(def, element)
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 63f6d141-e3e3-4b2f-b197-201cfd481111
begin
    using CairoMakie
    using PlutoUI
    using LinearAlgebra
end

# ╔═╡ 088fcb9e-b1ec-4ec8-b4bd-b5d9e3ec2222
md"""
# Nyquist Plots from Bode Data (without `ControlSystems.jl`)

This notebook computes frequency responses directly from polynomial coefficients and renders linked
**Bode + Nyquist (polar)** plots using **Makie**.

## Why this refactor?
1. `ControlSystems.jl` was removed to reduce environment size/startup cost.
2. `GR()` polar plots were replaced with `CairoMakie.PolarAxis` for better viewport control.

Use the sliders in each section to move a marker and to zoom the polar view (radius + angle window).
"""

# ╔═╡ 82b5ea57-f7d6-4609-8ca1-b7e1ef8b3333
begin
    struct RationalTF
        num::Vector{Float64}  # descending powers
        den::Vector{Float64}  # descending powers
    end

    polyval_desc(c::Vector{Float64}, s) = foldl((acc, a) -> acc*s + a, c[2:end]; init=c[1])
    eval_tf(G::RationalTF, s) = polyval_desc(G.num, s) / polyval_desc(G.den, s)

    const ω_range = 10 .^ range(log10(0.005), log10(500), length=1400)

    function freq_data(G::RationalTF, ω_arr)
        h = eval_tf.(Ref(G), im .* ω_arr)
        mag = abs.(h)
        phase_rad = angle.(h)
        phase_deg = rad2deg.(phase_rad)

        p_cross_ω = Float64[]
        p_cross_mag = Float64[]
        for i in 2:length(ω_arr)
            re1, re2 = real(h[i-1]), real(h[i])
            im1, im2 = imag(h[i-1]), imag(h[i])
            if re1 < 0 && re2 < 0 && im1 * im2 <= 0
                denom = abs(im1) + abs(im2)
                α = denom > 0 ? abs(im1) / denom : 0.5
                ω_c = exp(log(ω_arr[i-1]) * (1-α) + log(ω_arr[i]) * α)
                m_c = mag[i-1] * (1-α) + mag[i] * α
                push!(p_cross_ω, ω_c)
                push!(p_cross_mag, m_c)
            end
        end

        g_cross_ω = Float64[]
        g_cross_ph = Float64[]
        for i in 2:length(ω_arr)
            m1, m2 = mag[i-1], mag[i]
            if (m1 - 1) * (m2 - 1) <= 0
                denom = abs(m1 - 1) + abs(m2 - 1)
                α = denom > 0 ? abs(m1 - 1) / denom : 0.5
                ω_c = exp(log(ω_arr[i-1]) * (1-α) + log(ω_arr[i]) * α)
                p_c = phase_deg[i-1] * (1-α) + phase_deg[i] * α
                push!(g_cross_ω, ω_c)
                push!(g_cross_ph, p_c)
            end
        end

        return h, mag, phase_deg, phase_rad, p_cross_ω, p_cross_mag, g_cross_ω, g_cross_ph
    end

    function section_figure(G::RationalTF, label, ω_idx; neg_real_pts=[], rmax=1.0, θmin_deg=-210, θmax_deg=210)
        h, mag, phase_deg, phase_rad, p_cross_ω, p_cross_mag, g_cross_ω, g_cross_ph = freq_data(G, ω_range)

        ω_c = ω_range[ω_idx]
        mag_c = mag[ω_idx]
        ph_c = phase_deg[ω_idx]
        ph_rad_c = phase_rad[ω_idx]
        mag_dB = 20 .* log10.(max.(mag, 1e-12))
        mag_dB_c = 20log10(max(mag_c, 1e-12))

        fig = Figure(size=(1200, 680))
        Label(
            fig[0, 1:2],
            "Readout @ ω=$(round(ω_c, sigdigits=4)) rad/s   |G(jω)|=$(round(mag_c, sigdigits=5)) ($(round(mag_dB_c, digits=2)) dB)   ∠G(jω)=$(round(ph_c, digits=2))°",
            fontsize=18
        )

        ax_mag = Axis(
            fig[1, 1],
            xscale=log10,
            xlabel="ω (rad/s)",
            ylabel="Magnitude (dB)",
            title="Bode Magnitude — $label  (marker: $(round(mag_dB_c, digits=2)) dB)"
        )
        lines!(ax_mag, ω_range, mag_dB, color=:steelblue, linewidth=2)
        hlines!(ax_mag, [0.0], color=:black, linestyle=:dash)
        scatter!(ax_mag, [ω_c], [mag_dB_c], color=:green, markersize=12)
        for ωpc in p_cross_ω
            vlines!(ax_mag, [ωpc], color=:crimson, linestyle=:dot)
        end
        for ωgc in g_cross_ω
            vlines!(ax_mag, [ωgc], color=:darkorange, linestyle=:dot)
        end

        ax_phase = Axis(
            fig[2, 1],
            xscale=log10,
            xlabel="ω (rad/s)",
            ylabel="Phase (deg)",
            title="Bode Phase — $label  (marker: $(round(ph_c, digits=2))°)"
        )
        lines!(ax_phase, ω_range, phase_deg, color=:steelblue, linewidth=2)
        hlines!(ax_phase, [-180.0], color=:crimson, linestyle=:dash)
        hlines!(ax_phase, [0.0], color=:black, linestyle=:dot)
        scatter!(ax_phase, [ω_c], [ph_c], color=:green, markersize=12)
        for (ωgc, phgc) in zip(g_cross_ω, g_cross_ph)
            scatter!(ax_phase, [ωgc], [phgc], color=:darkorange, marker=:star5, markersize=13)
        end

        pax = PolarAxis(fig[:, 2], title="Nyquist Polar — $label",
            thetalimits=(deg2rad(θmin_deg), deg2rad(θmax_deg)),
            rlimits=(0.0, rmax))
        lines!(pax, phase_rad, mag, color=:steelblue, linewidth=2)
        lines!(pax, -phase_rad, mag, color=:steelblue, linewidth=1.5, linestyle=:dash)

        θ_unit = range(0, 2π, length=240)
        lines!(pax, θ_unit, ones(length(θ_unit)), color=:gray, linestyle=:dot)

        scatter!(pax, [π], [1.0], color=:black, marker=:diamond, markersize=10)
        for (_, r_nb, _) in neg_real_pts
            scatter!(pax, [π], [r_nb], color=:purple, marker=:star5, markersize=13)
        end
        for mpc in p_cross_mag
            scatter!(pax, [π], [mpc], color=:crimson, marker=:xcross, markersize=12)
        end
        scatter!(pax, [ph_rad_c], [mag_c], color=:green, markersize=12)
        scatter!(pax, [-ph_rad_c], [mag_c], color=:limegreen, markersize=9)

        fig
    end
end

# ╔═╡ 56f5984f-c0d6-4dd2-83c8-a77f877d4444
md"""
## Section 1 — \(G_1(s)=\frac{1}{(s+1)(s+2)}\), stable for \(K>-2\)
"""

# ╔═╡ 71f042f8-f79e-4707-9275-01f1ba7f5555
G1 = RationalTF([1.0], [1.0, 3.0, 2.0])

# ╔═╡ b913b8b6-f74a-4d8c-b6a6-547fc3e56666
@bind ω_idx1 Slider(1:length(ω_range), default=280, show_value=true)

# ╔═╡ d16fbe37-efb8-4c26-9b65-378f59487777
@bind rmax1 Slider(0.2:0.01:1.2, default=0.6, show_value=true)

# ╔═╡ e65be8ba-5723-4efb-8ec8-fca719f28888
@bind θzoom1 Slider(40:5:220, default=140, show_value=true)

# ╔═╡ 98abf2c4-17ef-428e-b4f6-53347e349999
section_figure(G1, "G₁(s)", ω_idx1; neg_real_pts=[], rmax=rmax1, θmin_deg=-θzoom1, θmax_deg=θzoom1)

# ╔═╡ e2be7a3b-a73d-4d8c-9d6a-8d4f5770aaaa
md"""
## Section 2 — \(G_2(s)=\frac{1}{(s-1)(s^2+2s+3)}\), stable for \(3<K<4\)
"""

# ╔═╡ c9f7b88d-2b63-4f9f-8664-56374790bbbb
G2 = RationalTF([1.0], [1.0, 1.0, 1.0, -3.0])

# ╔═╡ 2fb9d965-9c0a-4e60-8847-d6fbb593cccc
@bind ω_idx2 Slider(1:length(ω_range), default=390, show_value=true)

# ╔═╡ 6f56197a-3b8e-4d1f-8cdd-fef3b8a2dddd
@bind rmax2 Slider(0.2:0.01:0.8, default=0.45, show_value=true)

# ╔═╡ 875ce4c8-9f57-4f2f-8f8d-e47a99c8eeee
@bind θzoom2 Slider(40:5:220, default=120, show_value=true)

# ╔═╡ 9ffd6c95-d6a7-4f6a-9f32-88e58d39ffff
section_figure(G2, "G₂(s)", ω_idx2;
    neg_real_pts=[(0.0, 1/3, "K=3"), (1.0, 1/4, "K=4")],
    rmax=rmax2, θmin_deg=-θzoom2, θmax_deg=θzoom2)

# ╔═╡ 14dbde3f-3fc9-40e7-8d75-f66cbf620001
md"""
## Section 3 — \(G_3(s)=\frac{s-1}{(s+2)(s^2-s+1)}\), stable for \(\frac{3}{2}<K<2\)
"""

# ╔═╡ a3b8b48e-cb73-4f9c-b267-c2ba8bf60002
G3 = RationalTF([1.0, -1.0], [1.0, 1.0, -1.0, 2.0])

# ╔═╡ 2f20678d-e7f4-4731-a9a1-da58df280003
@bind ω_idx3 Slider(1:length(ω_range), default=300, show_value=true)

# ╔═╡ 1a768a4e-f2c7-4ec6-8a39-37e86b610004
@bind rmax3 Slider(0.2:0.01:1.0, default=0.72, show_value=true)

# ╔═╡ ef48c0f6-49b9-46cb-ac53-12fcfd520005
@bind θzoom3 Slider(40:5:220, default=130, show_value=true)

# ╔═╡ 0c0805e3-62a5-4858-bab3-69a248090006
section_figure(G3, "G₃(s)", ω_idx3;
    neg_real_pts=[(1/sqrt(2), 2/3, "K=3/2"), (0.0, 1/2, "K=2")],
    rmax=rmax3, θmin_deg=-θzoom3, θmax_deg=θzoom3)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
CairoMakie = "~0.15.7"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# Intentionally omitted: let Pluto resolve a lean manifest from the Project TOML.
"""

# ╔═╡ Cell order:
# ╠═63f6d141-e3e3-4b2f-b197-201cfd481111
# ╟─088fcb9e-b1ec-4ec8-b4bd-b5d9e3ec2222
# ╠═82b5ea57-f7d6-4609-8ca1-b7e1ef8b3333
# ╟─56f5984f-c0d6-4dd2-83c8-a77f877d4444
# ╠═71f042f8-f79e-4707-9275-01f1ba7f5555
# ╠═b913b8b6-f74a-4d8c-b6a6-547fc3e56666
# ╠═d16fbe37-efb8-4c26-9b65-378f59487777
# ╠═e65be8ba-5723-4efb-8ec8-fca719f28888
# ╠═98abf2c4-17ef-428e-b4f6-53347e349999
# ╟─e2be7a3b-a73d-4d8c-9d6a-8d4f5770aaaa
# ╠═c9f7b88d-2b63-4f9f-8664-56374790bbbb
# ╠═2fb9d965-9c0a-4e60-8847-d6fbb593cccc
# ╠═6f56197a-3b8e-4d1f-8cdd-fef3b8a2dddd
# ╠═875ce4c8-9f57-4f2f-8f8d-e47a99c8eeee
# ╠═9ffd6c95-d6a7-4f6a-9f32-88e58d39ffff
# ╟─14dbde3f-3fc9-40e7-8d75-f66cbf620001
# ╠═a3b8b48e-cb73-4f9c-b267-c2ba8bf60002
# ╠═2f20678d-e7f4-4731-a9a1-da58df280003
# ╠═1a768a4e-f2c7-4ec6-8a39-37e86b610004
# ╠═ef48c0f6-49b9-46cb-ac53-12fcfd520005
# ╠═0c0805e3-62a5-4858-bab3-69a248090006
# ╠═00000000-0000-0000-0000-000000000001
# ╠═00000000-0000-0000-0000-000000000002
