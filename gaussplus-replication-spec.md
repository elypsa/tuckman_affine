# Gauss+ Replication Spec

Implementation spec for coding up Tuckman & Serrat, *Fixed Income Securities*, Chapter 9 (Gauss+) and Appendix A9.2, in R. Goal: get the thing running end to end so the Annex can be read against working code.

**Read this instead of the book.** Every equation needed is restated here in implementable form. Several equations as printed in the Annex are wrong — coding from the printed Annex gives wrong loadings. Corrections are flagged inline and listed in §3.

---

## 1. Model definition

### 1.1 Cascade form (risk-neutral dynamics)

Three factors: `r` (short rate), `m` (medium), `l` (long). Two Brownian motions.

```
dr_t = α_r (m_t − r_t) dt
dm_t = α_m (l_t − m_t) dt + σ_m (ρ dW¹_t + sqrt(1−ρ²) dW²_t)
dl_t = α_l (μ − l_t) dt + σ_l dW¹_t
E[dW¹ dW²] = 0
```

**Correction.** The book prints (9.9)–(9.11) as `−α(target − current)`, which is explosive. Use the form above. Half-life = ln(2)/α.

`r` has **no diffusion term** — that's the "+" in Gauss+: three state variables, two risk sources.

Parameters `P = (α, σ, μ)` with `α = (α_r, α_m, α_l)`, `σ = (σ_m, σ_l, ρ)`, `μ` scalar. The Annex writes `α_s` for `α_r` — same thing.

### 1.2 Matrix form

With `x_t = (r_t, m_t, l_t)'` and `1` a column of ones:

```
d(x_t − μ1) = −K (x_t − μ1) dt + Ω dW_t
```

```
        ⎡ α_r   −α_r     0  ⎤              ⎡   0            0        0 ⎤
  K  =  ⎢  0     α_m   −α_m ⎥      Ω(σ) =  ⎢ ρσ_m   sqrt(1−ρ²)σ_m   0 ⎥
        ⎣  0      0     α_l ⎦              ⎣ σ_l            0        0 ⎦
```

`Ω` is (A9.8), correct as printed. Row 1 is zero because `r` has no shock. Column 3 is zero because there are only two Brownian motions; carry `dW_t` as `(dW¹, dW², 0)'` to keep everything 3×3.

### 1.3 Reduced form

`K` is upper triangular with distinct eigenvalues `α_r, α_m, α_l`:

```
K = A(α) diag(α) A(α)⁻¹
X_t = A(α)⁻¹ (x_t − μ1)      ⟺      x_t = A(α) X_t + μ1       [A9.5]
dX_t = −diag(α) X_t dt + A(α)⁻¹ Ω dW_t                        [A9.7]
```

Columns of `A(α)` are eigenvectors of `K` normalised so each top entry is 1, which makes the first row `(1,1,1)` and hence `r_t = μ + 1'X_t`.

---

## 2. Closed-form objects

### 2.1 A(α) — **corrected**

```
          ⎡ 1        1                      1                   ⎤
  A(α) =  ⎢ 0   (α_r−α_m)/α_r        (α_r−α_l)/α_r              ⎥
          ⎣ 0        0         (α_r−α_l)(α_m−α_l)/(α_r·α_m)     ⎦
```

**Correction.** (A9.6) prints the (3,3) entry with denominator `α_r`, omitting `α_m`.

### 2.2 A(α)⁻¹ — **corrected**

```
            ⎡ 1   −α_r/(α_r−α_m)      α_r α_m /((α_r−α_m)(α_r−α_l))  ⎤
  A(α)⁻¹ =  ⎢ 0    α_r/(α_r−α_m)     −α_r α_m /((α_r−α_m)(α_m−α_l))  ⎥
            ⎣ 0         0             α_r α_m /((α_r−α_l)(α_m−α_l))  ⎦
```

**Correction.** (A9.12) omits the minus sign on the (2,3) entry.

Simplest implementation: build `A` from the formula and use `solve()`. The closed form above is here for reference.

### 2.3 Loading functions

Two versions of `B` are needed; the book uses one symbol for both, which is a trap.

```
B_i(τ)      = (1 − exp(−α_i τ)) / (α_i τ)        # yield loading, "averaged"
Btilde_i(τ) = (1 − exp(−α_i τ)) / α_i            # = τ · B_i(τ), "unnormalised"
```

For small `α_i τ` use `-expm1(-a*tau)/(a*tau)`.

### 2.4 Yields

```
y_t(τ) = μ − C(τ,α,σ) + B(τ,α) X_t                           [A9.10]
       = μ (1 − Υ(τ,α) 1) − C(τ,α,σ) + Υ(τ,α) x_t            [A9.13]

Υ(τ,α) = B(τ,α) A(α)⁻¹      (1×3 row)                        [A9.14]
```

`Υ` depends on `α` only — not `σ`, `ρ` or `μ`. That's what makes the staged estimation work.

### 2.5 Convexity — **corrected sign**

```
                3   3    σ_ij        ⎛                        1 − exp(−(α_i+α_j)τ) ⎞
C(τ,α,σ)  =    ∑   ∑   ───────   ×   ⎜ 1 − B_i(τ) − B_j(τ) + ───────────────────── ⎟
               i=1 j=1  2 α_i α_j    ⎝                            (α_i+α_j) τ      ⎠

σ_ij = (i,j) entry of   S = A(α)⁻¹ Ω Ω' (A(α)⁻¹)'
```

**Corrections.** (a) (A9.11) prints a **minus** before the final fraction; it must be a **plus** — otherwise the bracket tends to `−2` rather than `0` as `τ → 0`. (b) The book writes `A⁻¹ΩΩ'A⁻¹`; the last factor must be transposed.

`C` depends on maturity and parameters only, never the state — so it disappears from anything written in differences.

### 2.6 Forward rates

Forward starting at `τ` with tenor `τ'`:

```
f_t(τ) = μ (1 − Υ'(τ,α,τ') 1) + Υ'(τ,α,τ') x_t − C'(τ,α,σ,τ')    [A9.15]

Υ'(τ,α,τ') = ( Btilde(τ+τ') − Btilde(τ) ) A(α)⁻¹ / τ'            [A9.16, corrected]
C'(τ,α,σ,τ') = ( (τ+τ') C(τ+τ') − τ C(τ) ) / τ'                  [A9.17, corrected]
```

**Correction.** (A9.16) and (A9.17) as printed omit the `τ` weights. They only work read as `Btilde` and `τ·C`. Implementing them literally gives wrong loadings everywhere.

**Instantaneous forwards** (`τ' → 0`):

```
Υ'_inst(τ,α) = [ exp(−α_r τ), exp(−α_m τ), exp(−α_l τ) ] A(α)⁻¹
```

Figure 9.7 plots the **instantaneous** version — it starts at `(1,0,0)` at `τ=0`, which the one-year-tenor version does not. Use `Υ'_inst` for that figure, `Υ'` with `τ'=1` everywhere forwards are used as data.

---

## 3. Errata summary

| Equation | As printed | Corrected |
| --- | --- | --- |
| (9.9)–(9.11) | `−α(target − current)` | `α(target − current)` |
| (A9.6) entry (3,3) | `(α_r−α_l)(α_m−α_l)/α_r` | `… /(α_r α_m)` |
| (A9.11) last bracket term | subtracted | added |
| (A9.11) `σ_ij` | `A⁻¹ΩΩ'A⁻¹` | `A⁻¹ΩΩ'(A⁻¹)'` |
| (A9.12) entry (2,3) | positive | negative |
| (A9.16) | uses `B` | uses `Btilde = τB` |
| (A9.17) | `(C(τ+τ') − C(τ,·,τ'))/τ'` | `((τ+τ')C(τ+τ') − τC(τ))/τ'` |
| (A9.29) | `f − λ·RP` | `f − λ·RP/Δτ` |

---

## 4. Reference parameters

Table 9.1, for optimiser starting values and for comparison once things run:

```r
par_ref <- list(
  a_r = 1.0547, a_m = 0.6358, a_l = 0.0165,
  s_m = 0.01092,   # 109.2 bp
  s_l = 0.00964,   #  96.4 bp
  rho = 0.212,
  mu  = 0.10555    # 10.555%
)
```

Half-lives: 0.66, 1.09, 42.0 years.

A few numbers the model should produce at these parameters, useful if something looks off: the 10-year yield loads about 0.70 on the long factor; the 10-year yield volatility comes out near 73 bp; the instantaneous forward loading on the long factor peaks around 0.90 near 7 years and is about 0.64 at 30 years.

---

## 5. Data

### 5.1 Yield curve — GSW / Federal Reserve Board

File: `feds200628.csv` from federalreserve.gov. (The book says New York Fed; it's a Board staff dataset.) Header rows precede the data — skip to the row starting `Date`. **Values are in percent, so divide by 100.**

Compounding differs by series:

| Series | Mnemonic | Compounding |
| --- | --- | --- |
| Zero-coupon yield | `SVENYxx` | **continuous** |
| Par yield | `SVENPYxx` | coupon-equivalent |
| Instantaneous forward | `SVENFxx` | **continuous** |
| One-year forward | `SVEN1Fxx` | coupon-equivalent |

**Don't use `SVEN1Fxx`** — it's coupon-equivalent and the model is continuous, and converting means guessing whether that's annual or semiannual. Build one-year forwards from the continuously compounded zeros instead:

```
f_cc(τ, τ'=1) = (τ+1)·SVENY{τ+1}/100 − τ·SVENY{τ}/100
```

Exact and consistent with the model's algebra. Maturities available are `SVENY01`–`SVENY30`. You need `SVENY02, SVENY03, SVENY10, SVENY11` for the benchmarks and `SVENY14, SVENY15, SVENY16` for λ.

### 5.2 Short rate — fed funds target

The sample is entirely in the target-*range* era, so use the midpoint:

```
r_t = (DFEDTARU + DFEDTARL) / 200          # FRED, daily, percent
```

Optionally convert to continuous with `log(1+r)` — the difference is under a basis point at these levels.

### 5.3 Sample and weights

```
2014-01-05  to  2022-01-21
```

Business days; inner-join GSW and target on date, drop dates missing either.

Exponential decay weights, 0.8 per year, normalised at the last observation:

```r
age_years <- as.numeric(max(dates) - dates) / 365.25
w <- 0.8 ^ age_years
```

Apply in all three estimation stages.

Factors are also extracted back to January 2007 with parameters held fixed, for Figure 9.8 — no weights, no re-estimation.

---

## 6. Estimation

Three stages, each conditioning on the one before.

### 6.1 Shapes

| Object | Shape | Meaning |
| --- | --- | --- |
| `Y` | T × N | zero yields, T days × N maturities |
| `dY` | (T−1) × N | daily changes |
| `Yb` | T × 2 | benchmarks (2y, 10y) |
| `Ups` | N × 2 | netted loadings, columns (m, l) |
| `Ups_b` | 2 × 2 | the two benchmark rows |

### 6.2 Stage 0 — netting

The short-rate loading is `Υ_s(τ,α) = B_1(τ)`, because the first column of `A⁻¹` is `(1,0,0)'`. It depends on `α_r` alone.

```
yhat_t(τ) = y_t(τ) − Υ_s(τ, α_r) · r_t
```

After this the factor vector is `(m_t, l_t)` and loading rows are `(Υ_m, Υ_l)` — the short column is dropped everywhere. Recompute the netting inside the stage-one optimiser, since it depends on `α_r`.

### 6.3 Stage 1 — α from regression slopes

Work in daily changes; the constants (`μ(1−Υ1)` and `C`) difference away, so this stage needs neither `σ` nor `μ`.

```r
# empirical: weighted OLS of all netted yield changes on the two benchmarks
beta_hat <- solve(t(dYb) %*% (w * dYb)) %*% t(dYb) %*% (w * dY)   # 2 × N

# model-implied
slopes <- Ups(alpha) %*% solve(Ups_b(alpha))                      # N × 2

# objective
min over alpha of   norm(t(slopes) - beta_hat, "F")
```

Constrain `α_r > α_m > α_l > 0`. Easiest parameterisation:

```r
a_l <- exp(p3); a_m <- a_l + exp(p2); a_r <- a_m + exp(p1)
```

Start from `par_ref`. The objective is flat in `α_l`, so if it wanders, fix `α_l` on a small grid and optimise the other two.

Gives Figure 9.5.

### 6.4 Stage 2 — σ from the volatility term structure

With `α̂` fixed:

```r
Sigma_model <- Ups_full %*% Omega(sigma) %*% t(Omega(sigma)) %*% t(Ups_full)   # N × N
v_emp       <- colSums(w * dY^2) / sum(w) * 252                               # length N

min over sigma of   norm(diag(Sigma_model) - v_emp, "2")
```

`Ups_full` is N × 3 (unnetted — `r` carries no shock so its column contributes nothing, but keeping it 3-wide avoids special cases).

**Match diagonals only.** The book writes a full matrix minus a diagonal, which would also push model yield correlations to zero; that isn't the intent.

Constraints `σ_m, σ_l > 0` and `−1 < ρ < 1`; parameterise `ρ = tanh(z)`.

Gives Figure 9.6.

### 6.5 Stage 3 — μ from the level

Only `μ` left. For each candidate: extract factors daily (§7), compute model yields at all maturities, score

```
min over mu of   Σ_t w_t · || Y_t − y_t(mu) ||²
```

Scalar optimisation — `optimize()` over `(0, 0.30)`. The benchmark maturities contribute zero error by construction, so `μ` is pinned by the other maturities.

Expect something near 10.5%. It's a risk-neutral long-run target, not a rate forecast — the long factor reverts to it with a 42-year half-life.

---

## 7. Factor extraction

Daily, with `P` fixed. `r_t` is observed, so only `(m_t, l_t)` are extracted, by fitting the **two- and ten-year one-year-tenor forwards exactly**:

```r
targets <- c(f_mkt_2y, f_mkt_10y)
const_i <- mu * (1 - sum(Ups_fwd(tau_i))) - C_fwd(tau_i) + Ups_fwd(tau_i)[1] * r_t
Lmat    <- rbind(Ups_fwd(2)[2:3], Ups_fwd(10)[2:3])     # 2 × 2
ml      <- solve(Lmat, targets - const)
```

Exact fit, no residual, one 2×2 solve per day. Cache everything that doesn't depend on the date.

(The book uses yields in the stage-one regressions and forwards for extraction. That's its own inconsistency; forwards here matches the text and Figure 9.8.)

For plotting, transform the long factor to its ten-year-ahead conditional mean:

```
L(l_t) = μ (1 − exp(−10 α_l)) + l_t · exp(−10 α_l)       [A9.24]
```

Just the OU conditional mean at a ten-year horizon; with `α_l = 0.0165` the weight on `l_t` is about 0.848.

---

## 8. Risk premia

### 8.1 Notation change

**From here `τ` is a fixed future date, not an interval.** Time to maturity is `τ − t`. Only the long factor is priced; `λ_t` is a scalar assumed highly persistent.

### 8.2 The two expressions for the same expected return

Strategy: buy the `τ+Δτ` bond at `t`, sell the `τ` bond, unwind at `τ`.

```
market:  E_t[R] = ( f_t(τ) − E_t[r_τ] ) · Δτ        (definition of forwards, no model)
model:   E_t[R] = ∫ of g over the long-end window    (integrate the drift of log P)
```

with, for remaining maturity `u`,

```
g(u) = λ_t · u · Υ₃(u,α) · σ_l  −  ½ · u² · Υ(u,α) Ω Ω' Υ(u,α)'
```

First term is the risk premium from (A9.25); second is the Itô correction from `dP/P` → `d(log P)`, i.e. convexity. The short rate cancels between the legs.

Collapsing to the endpoint and dividing by `Δτ`:

```
f_t(τ) − E_t[r_τ]  ≈  λ_t · RP_unit(τ)  −  C'(τ,α,σ,Δτ)
RP_unit(τ) = τ · Υ₃(τ,α) · σ_l
```

### 8.3 Extracting λ

Impose flat expectations beyond some long maturity, then difference two long maturities:

```
λ_t = ( f_t(τ') − f_t(τ) ) / ( RP_unit(τ') − RP_unit(τ) )      [A9.28]
```

The chapter uses `τ = 14`, `τ' = 15` with one-year tenor. Numerator from market data, denominator from the model.

### 8.4 Recovering the expectation

```
E_t[r_τ] = f_t(τ) − λ_t · RP_unit(τ) + C'(τ,α,σ,Δτ)            [A9.29, corrected]
```

In words: **expected future rate = forward − risk premium + convexity.**

Apply at `τ = 9`, `Δτ = 1` for the expected one-year rate nine years forward (Figure 9.11). Sense check from the chapter: a 3% forward, ~63 bp premium, ~24 bp convexity → 2.61%.

### 8.5 Comparison series (Figure 9.11)

The external benchmark is a Cleveland Fed ten-year real rate forecast plus long-run inflation (average of the Cleveland Fed ten-year inflation forecast, the Philadelphia Fed ATSIX ten-year forecast, and a trend-inflation EWMA with decay 0.987). This is the least reproducible part of the chapter — treat as optional.

---

## 9. Outputs

| Output | Content | From |
| --- | --- | --- |
| Table 9.1 | `α, σ, ρ, μ` + half-lives `ln2/α` | §6 |
| Figure 9.5 | regression coefficients on 2y and 10y, empirical vs model | §6.3 |
| Figure 9.6 | yield volatility in bp, empirical vs model | §6.4 |
| Figure 9.7 | instantaneous forward loadings per factor vs term | §2.6 |
| Figure 9.8 | extracted factors + 2y forward, 2007–2022, with `L(l_t)` | §7 |
| Figure 9.9 | 9y minus 5y fitting-error signal (5d minus 40d MA of demeaned errors) | §6.5 residuals |
| Figure 9.10 | risk premium on the 10y forward vs its level | §8 |
| Figure 9.11 | model-implied long-run expected short rate vs external estimate | §8.4–8.5 |

---

## 10. Package layout

```
R/
  01_data.R          load_gsw(), load_target(), build_panel()
  02_model.R         A_mat(), A_inv(), Omega(), B(), Btilde(),
                     Ups_yield(), Ups_fwd(), Ups_fwd_inst(), C_conv(), C_fwd()
  03_estimate.R      stage1_alpha(), stage2_sigma(), stage3_mu(), estimate_all()
  04_factors.R       extract_factors(), L_transform()
  05_premia.R        rp_unit(), solve_lambda(), expected_rate()
  06_figures.R       one function per figure
data-raw/            feds200628.csv, fed_funds_target.csv   (committed)
data/                panel.rds                              (built)
output/              figures, parameter table
run_all.R
renv.lock
```

Signatures:

```r
A_mat(alpha); A_inv(alpha); Omega(sigma)              # 3×3
B(tau, alpha); Btilde(tau, alpha)                     # length(tau) × 3
Ups_yield(tau, alpha)                                 # length(tau) × 3
Ups_fwd(tau, alpha, tenor = 1)                        # length(tau) × 3
Ups_fwd_inst(tau, alpha)                              # length(tau) × 3
C_conv(tau, alpha, sigma)                             # length(tau)
C_fwd(tau, alpha, sigma, tenor = 1)                   # length(tau)
```

Vectorise over `tau` — all of these get called at ~15 maturities inside optimiser loops.

renv dependencies, kept small:

```r
install.packages(c("data.table", "nloptr", "ggplot2"))
renv::snapshot()
```

`data.table` for the panel, `nloptr` for the bounded stage-one and stage-two optimisations (or `stats::optim` with `L-BFGS-B` — the problems are small), `ggplot2` for figures. `stats::optimize` covers stage three.

---

## 11. Pitfalls

**Divide GSW by 100.** A silent factor of 100 makes `σ` absurd and sends `μ` to the bracket edge.

**`B` versus `Btilde`.** The most likely source of a wrong answer. Forwards use `Btilde`; yields use `B`.

**Instantaneous versus one-year-tenor forwards.** Figure 9.7 uses instantaneous; data and factor extraction use one-year tenor.

**Recompute the netting inside the optimiser** — `Υ_s` depends on `α_r`.

**Annualisation.** 252 trading days, `τ` in years, yields in decimals.

**`expm1` for small arguments** — `α_l·τ` is about 0.016 at one year, where `(1-exp(-x))/x` loses precision.

**Naming.** The Annex uses `λ` for the price of risk, and GSW's Svensson parameters use `τ₁, τ₂`. Use explicit names (`price_of_risk`, `alpha`, `svensson_tau`) to avoid collisions.
