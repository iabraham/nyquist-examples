### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a1b2c3d4-1111-1111-1111-111111111111
begin
    using ControlSystems
    using Plots
    using PlutoUI
    gr()
end

# ╔═╡ a1b2c3d4-4444-4444-4444-444444444444
md"""
# Nyquist Plots from Bode Plots: Stability Analysis

This notebook illustrates how to construct Nyquist plots directly from Bode magnitude and phase data, and how to apply the Nyquist stability criterion.

**The three open-loop transfer functions studied:**

| System | Transfer Function | Open-loop RHP poles (P) | Stable range of K |
|--------|-------------------|--------------------------|-------------------|
| G₁(s)  | 1 / ((s+1)(s+2))  | 0  | K > −2 |
| G₂(s)  | 1 / ((s−1)(s²+2s+3)) | 1 | 3 < K < 4 |
| G₃(s)  | (s−1) / ((s+2)(s²−s+1)) | 2 | 3/2 < K < 2 |

**Unity-feedback loop**: characteristic equation is **1 + K·Gᵢ(s) = 0**.

**How to read the linked plots**: Use the **frequency slider** in each section to move a green
marker simultaneously on:
- **Bode magnitude** (top-left): reads off |G(jω)| in dB → this is the **radius r** on the polar plot
- **Bode phase** (bottom-left): reads off ∠G(jω) in degrees → this is the **angle θ** on the polar plot
- **Nyquist polar** (right): the green dot sits at (r = |G(jω)|, θ = ∠G(jω))

The polar representation makes the Bode-to-Nyquist mapping direct and visual.
"""

# ╔═╡ a1b2c3d4-2222-2222-2222-222222222222
# Shared frequency range: ω ∈ [0.005, 500] rad/s, log-spaced
const ω_range = 10 .^ range(log10(0.005), log10(500), length=1000)

# ╔═╡ a1b2c3d4-3333-3333-3333-333333333333
begin
    """
    Compute Bode data and Nyquist markers for a transfer function G over ω_arr.
    Returns: mag (linear), phase_deg, phase_rad, p_cross_ω, p_cross_mag, g_cross_ω, g_cross_ph

    Phase crossovers are found as points where Im(G(jω)) changes sign while Re(G(jω)) < 0
    (i.e., the Nyquist curve crosses the negative real axis for ω > 0).  This avoids the
    phase-wrapping artifacts that break a naïve ±180° sign-change test for systems whose
    Nyquist curve starts on the negative real axis (e.g. G₂ and G₃).
    """
    function freq_data(G, ω_arr)
        h = [G(im * w) for w in ω_arr]
        mag = abs.(h)
        phase_rad = angle.(h)
        phase_deg = rad2deg.(phase_rad)

        # Phase crossover: Im(G) changes sign AND Re(G) < 0
        p_cross_ω  = Float64[]
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

        # Gain crossover: where |G| = 1
        g_cross_ω  = Float64[]
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

        return mag, phase_deg, phase_rad, p_cross_ω, p_cross_mag, g_cross_ω, g_cross_ph
    end

    """
    Build the three-panel combined plot for one transfer function.

    Arguments:
    - G         : ControlSystems transfer function
    - G_label   : String label like "G₁(s)"
    - ω_arr     : frequency array
    - ω_idx     : current slider index
    - neg_real_pts : list of (ω_val, r_val, K_val) for negative real axis crossings
                     (shown as stability boundary markers on the polar plot)
    - show_margins : if true, annotate GM / PM
    - P_open    : number of open-loop RHP poles (for Nyquist criterion caption)
    """
    function section_plot(G, G_label, ω_arr, ω_idx, neg_real_pts; show_margins=true, P_open=0)
        mag, phase_deg, phase_rad, p_cross_ω, p_cross_mag, g_cross_ω, g_cross_ph =
            freq_data(G, ω_arr)

        ω_c     = ω_arr[ω_idx]
        mag_c   = mag[ω_idx]
        ph_c    = phase_deg[ω_idx]
        phr_c   = phase_rad[ω_idx]
        mag_dB  = 20 .* log10.(max.(mag, 1e-10))
        mag_dB_c = 20 * log10(max(mag_c, 1e-10))

        # ── Bode Magnitude ──────────────────────────────────────────────
        p_mag = plot(
            ω_arr, mag_dB,
            xscale  = :log10,
            xlabel  = "ω (rad/s)",
            ylabel  = "Magnitude (dB)",
            title   = "Bode Magnitude — $(G_label)",
            label   = "",
            linewidth = 2,
            color   = :steelblue,
            grid    = true,
            legend  = :topright,
            minorgrid = true
        )
        hline!(p_mag, [0], linestyle=:dash, color=:black, linewidth=1, label="0 dB")

        # Phase crossover vertical lines → gain margin annotation
        if show_margins && !isempty(p_cross_ω)
            for (ωpc, mpc) in zip(p_cross_ω, p_cross_mag)
                gm_dB = -20 * log10(max(mpc, 1e-10))
                vline!(p_mag, [ωpc], linestyle=:dot, color=:crimson, linewidth=1.5,
                    label="ωpc = $(round(ωpc, sigdigits=3)) (GM = $(round(gm_dB, digits=1)) dB)")
            end
        end
        # Gain crossover vertical lines → phase margin annotation
        if show_margins && !isempty(g_cross_ω)
            for ωgc in g_cross_ω
                vline!(p_mag, [ωgc], linestyle=:dot, color=:darkorange, linewidth=1.5,
                    label="ωgc = $(round(ωgc, sigdigits=3))")
            end
        end

        # Moving frequency marker
        scatter!(p_mag, [ω_c], [mag_dB_c],
            marker = :circle, markersize = 9, color = :green,
            label  = "ω=$(round(ω_c, sigdigits=3)), |G|=$(round(mag_dB_c, digits=1)) dB")

        # ── Bode Phase ───────────────────────────────────────────────────
        p_phase = plot(
            ω_arr, phase_deg,
            xscale  = :log10,
            xlabel  = "ω (rad/s)",
            ylabel  = "Phase (deg)",
            title   = "Bode Phase — $(G_label)",
            label   = "",
            linewidth = 2,
            color   = :steelblue,
            grid    = true,
            legend  = :topright,
            minorgrid = true
        )
        hline!(p_phase, [-180], linestyle=:dash, color=:crimson, linewidth=1.5, label="−180°")
        hline!(p_phase, [0],    linestyle=:dot,  color=:black,  linewidth=1,   label="0°")

        # Gain crossover markers on phase plot → phase margin
        if show_margins && !isempty(g_cross_ω)
            for (ωgc, phgc) in zip(g_cross_ω, g_cross_ph)
                pm = 180 + phgc
                scatter!(p_phase, [ωgc], [phgc],
                    marker=:star5, markersize=10, color=:darkorange,
                    label="PM = $(round(pm, digits=1))° at ω=$(round(ωgc, sigdigits=3))")
            end
        end
        # Phase crossover markers on phase plot
        if show_margins && !isempty(p_cross_ω)
            for ωpc in p_cross_ω
                vline!(p_phase, [ωpc], linestyle=:dot, color=:crimson, linewidth=1.5, label="")
            end
        end

        # Moving frequency marker
        scatter!(p_phase, [ω_c], [ph_c],
            marker = :circle, markersize = 9, color = :green,
            label  = "ω=$(round(ω_c, sigdigits=3)), ∠G=$(round(ph_c, digits=1))°")

        # ── Nyquist Polar Plot ───────────────────────────────────────────
        # Positive ω (ω > 0): solid curve
        # Negative ω (ω < 0): dashed mirror (same |G|, negated phase)
        # Unit circle for PM reference
        θ_unit  = range(0, 2π, length=200)
        r_unit  = ones(200)

        p_nyq = plot(
            phase_rad, mag,
            proj      = :polar,
            linewidth = 2,
            color     = :steelblue,
            label     = "ω > 0",
            title     = "Nyquist Polar — $(G_label)\n(r = |G(jω)|, θ = ∠G(jω))",
            legend    = :outertop
        )
        # Mirror image for ω < 0
        plot!(p_nyq,
            -phase_rad, mag,
            proj      = :polar,
            linewidth = 1.5,
            color     = :steelblue,
            linestyle = :dash,
            label     = "ω < 0 (mirror)"
        )
        # Unit circle
        plot!(p_nyq, θ_unit, r_unit,
            proj      = :polar,
            color     = :gray,
            linestyle = :dot,
            linewidth = 1,
            label     = "Unit circle (|G|=1)"
        )

        # Mark the −1 point (r=1, θ=π) — critical for K=1
        scatter!(p_nyq, [π], [1.0],
            proj       = :polar,
            marker     = :diamond,
            markersize = 9,
            color      = :black,
            label      = "−1 (K=1)"
        )

        # Mark stability boundary points on negative real axis (θ = π)
        colors_nb = [:crimson, :darkorange, :purple, :brown]
        for (k_idx, (ω_nb, r_nb, K_nb)) in enumerate(neg_real_pts)
            clr = colors_nb[min(k_idx, length(colors_nb))]
            scatter!(p_nyq, [π], [r_nb],
                proj       = :polar,
                marker     = :star5,
                markersize = 11,
                color      = clr,
                label      = "−1/K (K=$(K_nb)), r=$(round(r_nb, digits=3))"
            )
        end

        # Phase crossover points on polar (where curve crosses negative real axis, ω > 0)
        for (ωpc, mpc) in zip(p_cross_ω, p_cross_mag)
            gm_lin = 1/mpc
            scatter!(p_nyq, [π], [mpc],
                proj       = :polar,
                marker     = :xcross,
                markersize = 10,
                color      = :crimson,
                label      = "Phase cross. ω=$(round(ωpc, sigdigits=3)), GM=$(round(gm_lin, digits=2))×"
            )
        end

        # Moving marker on polar plot (positive and negative ω)
        scatter!(p_nyq, [phr_c], [mag_c],
            proj = :polar, marker = :circle, markersize = 10, color = :green, label="")
        scatter!(p_nyq, [-phr_c], [mag_c],
            proj = :polar, marker = :circle, markersize = 7, color = :limegreen,
            label="Marker (ω=$(round(ω_c, sigdigits=3)) rad/s)")

        # ── Combine into [[mag;phase], polar] layout ─────────────────────
        l = @layout [grid(2,1){0.5w} a{0.5w}]
        plot(p_mag, p_phase, p_nyq,
            layout = l,
            size   = (1100, 620),
            margin = 4Plots.mm
        )
    end
end

# ════════════════════════════════════════════════════════════════════════════
# SECTION 1 — G₁(s) = 1 / ((s+1)(s+2))
# ════════════════════════════════════════════════════════════════════════════

# ╔═╡ b1b2c3d4-1111-1111-1111-111111111111
md"""
---
## Section 1 — G₁(s) = 1 / ((s+1)(s+2))

---
### Routh–Hurwitz Analysis

Closed-loop characteristic equation with gain K:
```
1 + K · G₁(s) = 0   ⟹   s² + 3s + (2 + K) = 0
```

**Routh array:**

| Row | Col 1 | Col 2 |
|:---:|:-----:|:-----:|
| s²  |  1    | 2+K   |
| s¹  |  3    |   0   |
| s⁰  | 2+K   |       |

**Stability conditions** (all first-column entries > 0):
- 3 > 0  ✓  (always)
- 2 + K > 0   ⟹   **K > −2**

---
### Nyquist Analysis

Open-loop poles: s = −1, −2  → both in LHP → **P = 0** (no RHP open-loop poles).

Nyquist stability criterion:  **Z = N + P**,  so for stability (Z = 0): **N = 0** (zero net encirclements of the critical point −1/K).

**Key observations from the polar plot:**

- **K > 0**: the critical point −1/K lies on the **negative real axis** (r = 1/K, θ = π).
  G₁(jω) starts at G₁(0) = 0.5 (r = 0.5, θ = 0) and the phase travels from 0° to −180°
  asymptotically as ω → ∞, never crossing the negative real axis at finite ω.
  Therefore **N = 0** for all K > 0  →  always **stable**.

- **K < 0**: the critical point −1/K lies on the **positive real axis** (r = 1/|K|, θ = 0).
  The Nyquist curve passes through G₁(0) = 0.5 on the positive real axis.
  - K = −2: critical point at r = 0.5 = G₁(0)  → **marginal stability**.
  - K < −2: the critical point (r > 0.5) is encircled once CW (N = +1, Z = 1)  → **unstable**.

**Result**: K > −2  ✓  (matches Routh–Hurwitz)

**Gain/Phase Margins** (for nominal K = 1):
Phase never reaches −180° at finite ω → **GM = ∞**.
|G₁(jω)| ≤ G₁(0) = 0.5 < 1 for all ω ≥ 0 → no gain crossover → **PM = ∞**.
"""

# ╔═╡ b1b2c3d4-2222-2222-2222-222222222222
G1 = tf([1], [1, 3, 2])   # G₁(s) = 1 / (s² + 3s + 2)

# ╔═╡ b1b2c3d4-3333-3333-3333-333333333333
@bind ω_idx1 PlutoUI.Slider(1:length(ω_range), default=300, show_value=false)

# ╔═╡ b1b2c3d4-4444-4444-4444-444444444444
md"""
**Frequency slider** → ω = **$(round(ω_range[ω_idx1], sigdigits=4)) rad/s**  |  |G₁(jω)| = **$(round(abs(G1(im*ω_range[ω_idx1])), sigdigits=4))**  |  ∠G₁(jω) = **$(round(rad2deg(angle(G1(im*ω_range[ω_idx1]))), digits=2))°**
"""

# ╔═╡ b1b2c3d4-5555-5555-5555-555555555555
begin
    # For G₁ the only stability boundary from below is K = −2.
    # G₁(0) = 0.5 → on the *positive* real axis; there is no finite negative real axis crossing
    # for ω > 0. We mark K = −2 via the DC point on the positive real axis separately.
    # No negative-real-axis crossings in [0.005, 500] → empty list.
    g1_neg_real = []   # no negative real-axis crossings for ω > 0

    section_plot(G1, "G₁(s) = 1/((s+1)(s+2))", ω_range, ω_idx1,
                 g1_neg_real; show_margins=true, P_open=0)
end

# ════════════════════════════════════════════════════════════════════════════
# SECTION 2 — G₂(s) = 1 / ((s−1)(s²+2s+3))
# ════════════════════════════════════════════════════════════════════════════

# ╔═╡ c1b2c3d4-1111-1111-1111-111111111111
md"""
---
## Section 2 — G₂(s) = 1 / ((s−1)(s²+2s+3))

---
### Routh–Hurwitz Analysis

Expand the denominator: (s−1)(s²+2s+3) = s³ + s² + s − 3.

Closed-loop characteristic equation:
```
1 + K · G₂(s) = 0   ⟹   s³ + s² + s + (K − 3) = 0
```

**Routh array:**

| Row | Col 1 | Col 2 |
|:---:|:-----:|:-----:|
| s³  |  1    |   1   |
| s²  |  1    |  K−3  |
| s¹  |  4−K  |   0   |
| s⁰  |  K−3  |       |

**Stability conditions:**
- 1 > 0  ✓
- 1 > 0  ✓
- 4 − K > 0  ⟹  K < 4
- K − 3 > 0  ⟹  K > 3

⟹  **3 < K < 4**

---
### Nyquist Analysis

Open-loop pole at s = +1 (RHP)  → **P = 1**.
For stability (Z = 0):  **N = −1** (one CCW encirclement of −1/K required).

**Key observations from the polar plot:**

G₂(0) = 1/((−1)(3)) = **−1/3**  → the Nyquist curve *starts* on the **negative real axis**
at r = 1/3 (θ = π) when ω = 0.

There is a second negative real-axis crossing at **ω = 1 rad/s** where G₂(j·1) = **−1/4**
(r = 1/4, θ = π).  This is the **phase crossover frequency**.

- **K = 3**: critical point −1/K = −1/3 = G₂(0)  → system enters stable region from below.
- **K = 4**: critical point −1/K = −1/4 = G₂(j·1)  → system leaves stable region (upper boundary).
- **3 < K < 4**: critical point lies *between* the two crossings  → the full Nyquist contour
  makes exactly **one CCW encirclement** (N = −1 → Z = 0 → **stable**).
- **K > 4 or K < 3**: N ≠ −1 → Z ≠ 0 → **unstable**.

**Gain Margin** (from phase crossover at ω = 1 rad/s):
|G₂(j·1)| = 1/4, so **GM = 4 (= 12 dB)** — the maximum multiplicative factor by which K can be
increased from the lower stability boundary before becoming unstable.

**Phase Margin**: |G₂(jω)| < 1 for all ω  → no gain crossover  → **PM = ∞** in the conventional sense.
"""

# ╔═╡ c1b2c3d4-2222-2222-2222-222222222222
G2 = tf([1], [1, 1, 1, -3])   # G₂(s) = 1 / (s³ + s² + s − 3)

# ╔═╡ c1b2c3d4-3333-3333-3333-333333333333
@bind ω_idx2 PlutoUI.Slider(1:length(ω_range), default=400, show_value=false)

# ╔═╡ c1b2c3d4-4444-4444-4444-444444444444
md"""
**Frequency slider** → ω = **$(round(ω_range[ω_idx2], sigdigits=4)) rad/s**  |  |G₂(jω)| = **$(round(abs(G2(im*ω_range[ω_idx2])), sigdigits=4))**  |  ∠G₂(jω) = **$(round(rad2deg(angle(G2(im*ω_range[ω_idx2]))), digits=2))°**
"""

# ╔═╡ c1b2c3d4-5555-5555-5555-555555555555
begin
    # G₂: negative real axis crossings at ω≈0 (r=1/3, K=3) and ω=1 (r=1/4, K=4)
    # ω=0 is represented by the DC limit (not in our log range), so we add it manually.
    # The phase crossover crossing at ω=1 is found automatically by freq_data.
    # We annotate the two boundary K values:
    g2_neg_real = [
        (0.0,  1/3, "3 (lower boundary)"),
        (1.0,  1/4, "4 (upper boundary / GM)")
    ]

    section_plot(G2, "G₂(s) = 1/((s−1)(s²+2s+3))", ω_range, ω_idx2,
                 g2_neg_real; show_margins=true, P_open=1)
end

# ════════════════════════════════════════════════════════════════════════════
# SECTION 3 — G₃(s) = (s−1) / ((s+2)(s²−s+1))
# ════════════════════════════════════════════════════════════════════════════

# ╔═╡ d1b2c3d4-1111-1111-1111-111111111111
md"""
---
## Section 3 — G₃(s) = (s−1) / ((s+2)(s²−s+1))

---
### Routh–Hurwitz Analysis

Expand denominator: (s+2)(s²−s+1) = s³ + s² − s + 2.

Closed-loop characteristic equation:
```
1 + K · G₃(s) = 0
  ⟹  s³ + s² + (K−1)s + (2−K) = 0
```

**Routh array:**

| Row | Col 1 | Col 2 |
|:---:|:-----:|:-----:|
| s³  |   1   |  K−1  |
| s²  |   1   |  2−K  |
| s¹  |  2K−3 |   0   |
| s⁰  |  2−K  |       |

**Stability conditions:**
- 1 > 0  ✓
- 1 > 0  ✓
- 2K − 3 > 0  ⟹  K > 3/2
- 2 − K > 0   ⟹  K < 2

⟹  **3/2 < K < 2**

---
### Nyquist Analysis

Poles of G₃: s = −2 (LHP) and roots of s²−s+1 = 0, i.e. s = (1 ± j√3)/2
→ Re(s) = +1/2 > 0  →  **two RHP poles → P = 2**.

Also a RHP zero at s = +1 (non-minimum phase).

For stability (Z = 0):  **N = −2** (two CCW encirclements required).

**Key observations from the polar plot:**

The Nyquist curve has **two** crossings of the negative real axis:

1. **ω = 0**:  G₃(0) = (−1)/((2)(1)) = **−1/2** → r = 1/2, θ = π → K = 2 boundary.
2. **ω = 1/√2 ≈ 0.707 rad/s**:  G₃(j/√2) = **−2/3** → r = 2/3, θ = π → K = 3/2 boundary.

(These can be verified algebraically: setting Im(G₃(jω)) = 0 gives ω = 0 or ω = 1/√2.)

- **3/2 < K < 2**: critical point −1/K ∈ (−2/3, −1/2) → the full Nyquist contour makes
  exactly **two CCW encirclements** (N = −2 → Z = 0 → **stable**).
- **K > 2 or K < 3/2**: N ≠ −2 → Z ≠ 0 → **unstable**.

**Gain Margins** (from the two phase crossovers):
- At ω = 0: r = 1/2  → GM₁ = 1/(1/2) = **2** (= 6 dB)  → upper bound (K < 2).
- At ω ≈ 0.707 rad/s: r = 2/3 → GM₂ = 1/(2/3) = **3/2** (≈ 3.5 dB) → lower bound (K > 3/2).

**Phase Margin**: |G₃(jω)| < 1 for all ω  → no gain crossover  → **PM = ∞** in the conventional sense.
"""

# ╔═╡ d1b2c3d4-2222-2222-2222-222222222222
G3 = tf([1, -1], [1, 1, -1, 2])   # G₃(s) = (s−1) / (s³ + s² − s + 2)

# ╔═╡ d1b2c3d4-3333-3333-3333-333333333333
@bind ω_idx3 PlutoUI.Slider(1:length(ω_range), default=300, show_value=false)

# ╔═╡ d1b2c3d4-4444-4444-4444-444444444444
md"""
**Frequency slider** → ω = **$(round(ω_range[ω_idx3], sigdigits=4)) rad/s**  |  |G₃(jω)| = **$(round(abs(G3(im*ω_range[ω_idx3])), sigdigits=4))**  |  ∠G₃(jω) = **$(round(rad2deg(angle(G3(im*ω_range[ω_idx3]))), digits=2))°**
"""

# ╔═╡ d1b2c3d4-5555-5555-5555-555555555555
begin
    # G₃: negative real axis crossings at ω=1/√2 (r=2/3, K=3/2) and ω→0 (r=1/2, K=2)
    g3_neg_real = [
        (1/sqrt(2),  2/3, "3/2 (lower boundary)"),
        (0.0,        1/2, "2 (upper boundary)")
    ]

    section_plot(G3, "G₃(s) = (s−1)/((s+2)(s²−s+1))", ω_range, ω_idx3,
                 g3_neg_real; show_margins=true, P_open=2)
end

# ╔═╡ Cell order:
# ╟─a1b2c3d4-1111-1111-1111-111111111111
# ╟─a1b2c3d4-4444-4444-4444-444444444444
# ╟─a1b2c3d4-2222-2222-2222-222222222222
# ╟─a1b2c3d4-3333-3333-3333-333333333333
# ╟─b1b2c3d4-1111-1111-1111-111111111111
# ╠═b1b2c3d4-2222-2222-2222-222222222222
# ╠═b1b2c3d4-3333-3333-3333-333333333333
# ╟─b1b2c3d4-4444-4444-4444-444444444444
# ╟─b1b2c3d4-5555-5555-5555-555555555555
# ╟─c1b2c3d4-1111-1111-1111-111111111111
# ╠═c1b2c3d4-2222-2222-2222-222222222222
# ╠═c1b2c3d4-3333-3333-3333-333333333333
# ╟─c1b2c3d4-4444-4444-4444-444444444444
# ╟─c1b2c3d4-5555-5555-5555-555555555555
# ╟─d1b2c3d4-1111-1111-1111-111111111111
# ╠═d1b2c3d4-2222-2222-2222-222222222222
# ╠═d1b2c3d4-3333-3333-3333-333333333333
# ╟─d1b2c3d4-4444-4444-4444-444444444444
# ╟─d1b2c3d4-5555-5555-5555-555555555555
# ╟─00000000-0000-0000-0000-000000000001

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ControlSystems = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
ControlSystems = "1"
Plots = "1"
PlutoUI = "0.7"
julia = "1.9"
"""
