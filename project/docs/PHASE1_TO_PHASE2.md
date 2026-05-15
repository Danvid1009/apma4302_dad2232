# How Phase 1 ties into Phase 2 (for your report)

## One-sentence bridge

**Phase 1** builds and validates the **parallel Monte Carlo engine** (correctness, scaling, variance behavior) under **controlled inputs**; **Phase 2** keeps that **same engine** and replaces synthetic assumptions with **data-informed inputs** (mainly volatility and spot from history), so you can discuss **model risk** and **statistical vs computational cost** on realistic-looking settings.

## What stays identical (the “spine”)

- **Path generation:** GBM on a discrete time grid; same `simulate_gbm_paths` machinery.
- **Parallelism:** path batches per MPI rank, **global reductions** for sums / counts; same HPC story (strong/weak scaling, efficiency).
- **Payoffs and discounts:** European / Asian / barrier logic unchanged; discount factor \(e^{-rT}\) unchanged unless you extend the model.
- **Uncertainty quantification:** standard errors and CIs from batch means—Phase 2 does not remove sampling error; it changes **which \(\sigma\)** (and often **which \(S_0\)**) you plug in.

So Phase 2 is **not** a different numerical method; it is a **different calibration layer** on top of the Phase 1 **solver**.

## What Phase 1 is responsible for (evidence you must have first)

1. **Correctness:** European vs Black–Scholes, CI contains truth, convergence slope near \(N^{-1/2}\).
2. **HPC behavior:** strong/weak scaling curves, interpretation of overhead vs rank count.
3. **Algorithm knobs:** e.g. antithetic pairs, (optional) control variate—these are **variance–compute** tradeoffs that exist **before** you touch market data.

Without Phase 1, Phase 2 is just “numbers from Yahoo”—you could not separate **implementation bugs**, **Monte Carlo error**, and **model mismatch**.

## What Phase 2 adds (interpretation layer)

1. **Empirical dynamics:** use **adjusted closes** to estimate **annualized realized volatility** (log-return standard deviation, scaled by \(\sqrt{252}\) for daily data).
2. **Spot alignment:** often take \(S_0\) as the **last quoted price** in the window (or a date you label in the report).
3. **Side-by-side comparison:** run the **same option** with a **synthetic baseline** \(\sigma\) (Phase 1 style) vs **calibrated** \(\sigma\) (Phase 2). The **spread in price and CI width** is the object of study.
4. **Optional “market face”:** compare to **observed option mid quotes** (Tier B)—then Phase 2 also speaks to **smile / jump risk** not captured by GBM.

## How to phrase the “nice tie” in prose

- Paragraph 1: Phase 1 proves the **pipeline is correct and scalable** under idealized parameters.  
- Paragraph 2: Phase 2 **reuses that pipeline** and asks how sensitive outputs are when \(\sigma\) and \(S_0\) come from **the same era of data** you show in a price/vol figure.  
- Paragraph 3: Remaining gap (quotes, dividends, stochastic vol) is **modeling depth**, not **MPI depth**—keeps the HPC narrative honest.

## Repo pointers

- Phase 1 drivers: `run_pricing.py`, `validate_european.py`, `run_convergence.py`, `run_scaling.py`, registry + gallery plots under `visuals/`.
- Phase 2 drivers: `phase2_calibrate_and_price.py`, `fetch_yahoo_ohlcv.py`, `plot_spy_volatility_regime.py`, optional `compare_mc_to_quotes.py`, `run_phase2_roll_forward.py`.

---

## Deeper tie: what is “the same” vs what is “new uncertainty”

### Same object: a discounted payoff expectation

Both phases estimate \(\mathbb{E}[e^{-rT}\,\Phi(S_\cdot)]\) under a **GBM law** on a discrete grid. Phase 1 emphasizes that the **discretization + RNG + MPI reduction** reproduces known benchmarks (Black–Scholes for Europeans). Phase 2 does not change that expectation operator; it changes the **parameters** \((S_0,\sigma)\) (and optionally \(r,T,K\)) fed into the same simulator.

So any Phase 2 headline (“calibrated price moved by X bps”) decomposes cleanly in the report:

1. **Monte Carlo noise** — standard error from finite \(N\); shrink with paths, antithetic, or control variates (still Phase-1-class knobs).
2. **Parameter movement** — \(\Delta\) price from moving \(\sigma\) or \(S_0\); this is **sensitivity analysis**, not parallelism.
3. **Model gap** — if you add **market quotes**, residual vs mid is **smile / jump / dividend / microstructure** not captured by GBM; that belongs in Phase 2 **interpretation**, not in “MPI efficiency.”

That decomposition is the “nice tie”: Phase 1 earns the right to trust (1); Phase 2 spends most words on (2)–(3).

### Statistical depth: variance reduction stays in the Phase 1 toolbox

- **Antithetic pairs** (engine flag): negative correlation across paired paths often cuts variance of **symmetric** payoffs at similar cost per step; `run_convergence.py --antithetic` logs the same \(N^{-1/2}\) story with a different constant.
- **Control variate (Asian vs European)** — `petsc_asian_call_cv_european_control` uses the analytic European value on the **same terminal stock** as a control for the Asian payoff, with a pooled estimate of the optimal coefficient \(c^\*\). Same MPI path batches; **less variance per path** ⇒ fewer paths for a target CI — an HPC-relevant “effective \(N\)” story without changing the market-data layer.

### Phase 2 depth: rolling calibration and “many runs along time”

`run_phase2_roll_forward.py` steps through history: for each snapshot \(t_i\), estimate trailing realized \(\sigma_i\) from the last `vol_window` closes, set \(S_{0,i}\) to the spot at \(t_i\), and re-price. That is **Phase 1 MC repeated** with **time-varying inputs** — good for discussing **non-stationarity** (vol clusters; ATM strike drifts with spot) while keeping the **same parallel code path** each time. Plot output with `scripts/plot_phase2_roll_forward.py`.

### One paragraph you can paste into Discussion

We first verified the implementation under synthetic parameters (European vs Black–Scholes, convergence in \(N\), and PETSc scaling). We then held the estimator fixed and replaced \((S_0,\sigma)\) with values implied from adjusted closes, optionally comparing to quoted option mids. Discrepancies after calibration isolate **model error** from **numerical error**, while variance-reduction experiments report whether **statistical efficiency** or **raw rank scaling** is the binding constraint for a target confidence interval.
