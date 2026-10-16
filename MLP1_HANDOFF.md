# MLP1 Handoff

Context for continuing work on `MLP1_Simulations.ipynb`. Written at the end of the session on
branch `claude/blissful-edison-793yk7` (6 commits ahead of `main`, working tree clean).

---

## 1. What the model is

A 2 → 2 → 1 molecular MLP.

- **Hidden layer:** two STAR perceptrons, nodes `j` and `k`.
- **Output layer:** one trans-splicing ribozyme perceptron, which uses the **General Model**
  kinetics (identical ODEs to `GeneralModel_Simulations.ipynb`).

### Wiring between the layers

The released transcripts of layer 1 **are** the input species of layer 2 — they are not
re-transcribed from a plasmid. Specifically:

- `Y⁽¹⁾_j` plays the role of the ribozyme's **input RNA**.
- `Y⁽¹⁾_k` plays the role of the ribozyme's **sequestering RNA**.

Written this way, the downstream loading (retroactivity) terms appear automatically in
`dy_j` / `dy_k` with no double counting. **Do not** model this by composing two isolated
activation functions — see §5.

### State vector (16 species, `OUT = 14`)

```
[m_j, z_j, p_j, c_j, y_j, q_j,   m_k, z_k, p_k, c_k, y_k, q_k,   p_2, c_2, y_2, q_2]
```

### Why `X0` seeds `p_j`/`p_k` but not `p_2`

Layer 1 (STAR) targets are **catalytic**: `dp_n = d*c - a*p*m + v*c`. No production term, no decay
term, so `p_n + c_n` is a conserved quantity fixed entirely by the initial condition. Setting
`p_j(0) = 0` kills node j permanently (verified: `y_j = 0.00000`, `y_2 = 0.0000` everywhere).
`pTot1 = 0.5` is a real model parameter that happens to live in `X0` because that is the only
place a conserved quantity can live.

Layer 2's target is produced at `theta` and decays at `phi`, so it self-fills and its initial
value is irrelevant (verified: `p_2(0)` of 0.0 / 1.0 / 5.0 all converge to `p_2 = 1.7003`,
`y_2 = 1.1816`, identical to 4 d.p.).

---

## 2. Notation

Decided during the session, applied throughout:

| Concept | Form | Note |
|---|---|---|
| Layer | superscript `⁽¹⁾` / `⁽²⁾` | square brackets `^[l]` recommended for the writeup |
| Node / edge | subscript | |
| Weight | `w^[l]_ji` | **destination-first**: into node `j`, from input `i` |
| Physical realization | `α⁺`, `α⁻` | promoter strengths; `w ∝ α⁺ − α⁻` |

**Transpose trap.** The reaction-scheme figure labels edges source-first (`α⁺_{1,j}` = from input 1
into node j). The code is destination-first `[node, input]`, so the matrix is the **transpose** of the
figure's subscript ordering. `alphaPos[0,1]` is `α⁺_{2,j}` — input 2 into node j. With the current
symmetric values this reads identically either way; it will not once the symmetry is broken.

There are **8 alphas**: `alphaPos` holds 4 (2 nodes × 2 inputs) and `alphaNeg` holds 4.

---

## 3. Notebook structure (24 cells)

Mirrors `GeneralModel_Simulations.ipynb` on `main`, including the alpha-through-`ratesDefault`
pattern and all the features ported from it.

| Cell | Contents |
|---|---|
| 0–2 | Title, imports, model markdown (ODEs + species table) |
| 3–6 | Analytic composition section (`starNodeOutput`, `ribozymeOutput`, `networkActivationFunction`) + boundary and activation plots |
| 7–9 | `contourf` helper, `mlpODEs` |
| 10–11 | `K_C`/`K_Q` markdown, grid + `ratesDefault` + `getVariantRates` / `getVariantLabel` |
| 12–14 | Variant loop (heatmaps + activation fns), results display, **sanity checks** |
| 15–16 | `plotActivationFunctionVariants` + example call |
| 17–18 | `plotDecisionBoundaryVariants` + example call |
| 19 | Alpha comparison (5 variants, black/reds/blues) |
| 20 | `plotHeatmap` + `heatmapForIllustrator` |
| 21 | 3D surface |
| 22–23 | Separability analysis |

### Conventions carried over from the General Model

- `mlpODEs(X, t, x1, x2, rates)` — `alphaPos`/`alphaNeg` are unpacked from `rates`, not passed
  separately. In `ratesDefault` they are stacked per variant as `(nVariants, 2, 2)` arrays, so
  `getVariantRates` indexes them like any other rate and returns a `[node, input]` matrix.
- `w1` / `w2` extra input weights; call sites pass `(w1 * X1[j, i], w2 * X2[j, i], rates)`.
- `getVariantLabel` supports `"alpha"`, `"xi2"`, `"K_C2"`, `"K_Q2"`, `"theta"`, `"None"`.
  The `"alpha"` label shows the **net drive** `αPos − αNeg`, because printing both matrices came to
  50 characters and was unreadable in a legend.
- SVG export uses the bare two-line form (`p.output_backend = "svg"` then `export_svgs(...)`),
  matching the other notebooks. **This raises `RuntimeError` when no webdriver is on PATH.**
  A robust fallback helper existed briefly and was deliberately removed for consistency
  (commit `79976f5`). Fix on the user's machine: `conda install -c conda-forge firefox geckodriver`,
  and launch Jupyter from a terminal so the kernel inherits the shell PATH — a kernel started from
  Anaconda Navigator does not, which is why the PATH appears to "break again" intermittently.

---

## 4. Established results (all measured, not assumed)

### No bias anywhere — this is deliberate

Commit `d20010a` removed all `α_{n,0}` terms so the model matches the drawn reaction scheme exactly.
Consequence: every production term is proportional to an input, so `net_n(0,0) = 0` and **every
hidden decision boundary is pinned to the origin.**

Measured cost (cell 23, random search over the alpha matrices, margin = `min(ON) − max(OFF)`
against an output range of ~2.4):

| target | trials | best margin | verdict |
|---|---|---|---|
| XOR | 120 | +0.0086 | **not reachable** |
| OR | 50 | +2.6969 | reachable |
| AND | 50 | +1.7508 | reachable |

With a bias cassette added, XOR reached **33.5× contrast** (`(0,1)` and `(1,0)` → 2.386,
`(1,1)` → 0.071). That configuration is **not** in the notebook; it needed `gain_k = 6.0`,
`biasNeg = [0.5, 9.0]`, `b_2 = 2000`, `a_2 = 1.0`.

**A bias cassette must sit on DNA that is not an input** — its own plasmid, or the backbone already
carrying the target templates — held at fixed concentration. Placed on X1 or X2 it is *not* a bias:
it just re-weights that edge and the boundary stays pinned to the origin. Verified: a constitutive
STAR cassette on its own plasmid gives `y_2(0,0)` of 0.858 / 2.384 / 2.692 at strengths
0.1 / 0.3 / 0.6, while the same cassette on X1 gives 0.0000 at every strength.

`theta` is **not** a bias — it produces the ribozyme *target* `P⁽²⁾`, so it is a gain. Tripling it
moves max output 33% but the relative threshold only 11%.

### Two rules for shaping decision boundaries

1. **Equal weights give exactly zero curvature.** Every `[[1,1],[1,1]]` config scored `curv = 0.000`
   regardless of kinetics. Both nodes see identical drive, output depends only on `x1 + x2`, and the
   level set is a straight anti-diagonal. The weight matrix rows must differ or no kinetics will
   bend the heatmap.
2. **Larger α makes boundaries straighter, not curvier** (the opposite of the intuitive guess):
   scaling a split-weight config ×1 / ×2 / ×4 / ×8 / ×16 gave `curv` = 0.131 / 0.116 / 0.061 /
   0.026 / 0.009. Once both nodes saturate hard the output becomes a step in `sign(x1 − x2)`.
   Curvature needs *partial, asymmetric* saturation: keep α near 1 and create the asymmetry via
   `a_j` vs `a_k` and `v_j` vs `v_k`. At ×4 and above the integrator also emits
   `ODEintWarning: Excess work done on this call` — large α makes the system stiff.

### Not reachable

A hot region near the origin, or a diagonal stripe — both saturate the entire square
(`hotfrac = 0.99`, no contour). With no bias the hot region always runs to a corner or an edge.
A closed/ring-shaped region requires the bias cassette.

### Low decision-boundary thresholds produce blank plots

At `threshold = 0.01` exactly **one** grid cell (the origin) falls below the level; one grid step
away (0.0204) the output is already 0.053. `find_contours` returns a 2-point degenerate contour,
invisible at plot scale. `plotDecisionBoundaryVariants` only warns when `len(contours) == 0`, so
this fails silently. Use `relative=True` with 0.3–0.5.

---

## 5. Known problem: the analytic section disagrees with the ODEs

The closed-form composition in cells 3–6 does **not** match the coupled ODE results:

| x1 | x2 | ODE `y_2` | analytic `y_2` |
|---|---|---|---|
| 0.2 | 0.2 | 0.6080 | 0.0071 |
| 0.6 | 0.2 | 1.6258 | 0.0199 |
| 1.0 | 1.0 | 0.9656 | 0.0121 |

**The root cause is upstream, not in this notebook.** `GeneralModel_Simulations.ipynb`'s own
single-node polynomial does not reproduce that notebook's own ODE simulation:

| x1 | ODE y | analytic y | analytic × (1/tau) |
|---|---|---|---|
| 0.2 | 0.8004 | 0.0108 | 0.1082 |
| 1.0 | 0.9069 | 0.0816 | 0.8158 |
| 2.0 | 0.9081 | 0.0980 | 0.9797 |

Rescaling by `1/tau` fixes saturation but not the low-input shape (0.108 vs 0.800), so it is not
a normalization slip — the cubic's input dependence looks wrong. **The ODE sections are the
reliable ones.** Re-deriving the single-node polynomial is an open task that would fix the
General Model, STAR, TMSE and MLP1 notebooks together.

---

## 6. Open items

1. **Re-derive the single-node activation polynomial** (§5). Highest value — it unblocks the
   analytic half of four notebooks.
2. **`MLP1_Reaction.png` is referenced by cell 2 but not in the repo.** No PNGs are tracked.
   If the bias is ever added, the figure needs two more arrows per hidden node
   (`∅ --α⁺_{n,0}--> M⁽¹⁾_n` and `∅ --α⁻_{n,0}--> Z⁽¹⁾_n`).
3. **`pTot1` is a bare module-level scalar**, shared by both hidden nodes and absent from
   `ratesDefault`, so the variant sweep cannot vary it and nodes j/k cannot have different target
   supplies. It materially scales output (0.5 → 1.0 took `y_2` from 1.18 to 2.45).
4. **The degenerate-contour warning** in `plotDecisionBoundaryVariants` (§4).
5. **Runtime**: the variant loop is 50×50×3 solves of a 16-species system (~25 ms each, a few
   minutes), and `plotDecisionBoundaryVariants` repeats it. Drop the grid to 30 if iterating.

---

## 7. Where the session left off

The user was choosing three variants for a figure showing the range of **heatmaps** (not decision
boundary plots) the MLP can produce. They had locally redesigned `ratesDefault` to take three
distinct matrix sets — **that redesign is in the chat only, not committed.** The committed
`ratesDefault` still sweeps `theta` with identical weights across variants.

### Recommended trio (measured at 40×40, base kinetics otherwise)

| | weights | kinetics overrides | curv | ON frac | max y₂ |
|---|---|---|---|---|---|
| **V1** rounded front | `αP=[[1.0,0.6],[0.6,1.0]]`, `αN=[[0.2,0.4],[0.4,0.2]]` | `b_2=20, a_k=20, d_k=2, v_k=2` | **0.255** | 0.40 | 0.250 |
| **V2** dome | `αP=[[0.6,1.0],[1.0,0.6]]`, `αN=[[0.4,0.2],[0.2,0.4]]` | `a_j=40, v_j=4` | 0.155 | 0.54 | 0.271 |
| **V3** tight band | `αP=[[1.2,0.5],[0.5,1.2]]`, `αN=[[0.25,0.5],[0.5,0.25]]` | `a_k=40, v_k=4` | 0.133 | 0.21 | 0.267 |

`curv` = mean perpendicular deviation of the 50% contour from its best-fit line, normalized by the
contour's extent. ON fraction matters for heatmaps: above ~0.9 or below ~0.05 the panel renders as a
flat block.

Base kinetics for all three (identical to the committed `ratesDefault` defaults):
`phi=0.1, theta=0.3, a_j=10, d_j=1, v_j=1, b_j=200, u_j=0.2, a_k=10, d_k=1, v_k=1, b_k=200,
u_k=0.2, a_2=10, d_2=1, v_2=1, delta_2=1, b_2=200, u_2=0.2`

### Caveats on that trio

- Four of the six candidates explored differed **only in kinetics**, not weights. If the figure
  claims "range of boundaries this architecture can form", decide whether the claim is about weight
  tuning or kinetic tuning. The trio above gives each variant its own matrices.
- `[[1.0,0.6],[0.6,1.0]]` and its transpose `[[0.6,1.0],[1.0,0.6]]` produce the **same surface
  reflected** (identical stats to 3 d.p.) — the model is symmetric under swapping
  (node j ↔ node k) with (X1 ↔ X2). Do not use both as evidence of range.
- The user's own draft `alphaPos2`/`alphaNeg2` gave net `W = [[0.05, 0.1], [-0.1, -0.05]]` — node k's
  drive is negative on both inputs, so node k never activates and the output node loses its
  inhibitor, degenerating toward a single perceptron (`curv = 0.014`). Flagged, not yet resolved.

---

## 8. Verification habits used here

Every claim in this document was measured. Useful checks to keep:

- **Cell 14 sanity checks**: layer-1 target conservation (`p + c = pTot1` — a real invariant that
  catches an accidental decay or production term), steady-state residual (`max |dX/dt| ≈ 1.7e-9`),
  and the analytic-vs-ODE comparison.
- **The notebook is validated by executing every code cell** with a reduced grid, stubbing
  `biocircuits`/`show`/`export_svgs`. Two bugs were caught only this way: a stub returning `None`
  where bokeh returns a list (`for fn in export_svgs(...)`), and a "passing" SVG export that wrote
  no file because the stub swallowed it. **If a test passes suspiciously easily, check the stub.**
- **Scalar metrics undersell shape.** `curv` ranked an L-shaped corner below a near-straight
  diagonal. Render the heatmaps and look before choosing.
- **Ratio metrics mislead near zero**: an early XOR score of "2.2× contrast" was 0.048 vs 0.015.
  Absolute margin is the honest metric.
