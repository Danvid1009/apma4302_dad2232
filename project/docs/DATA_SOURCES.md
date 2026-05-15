# Data sources (optional, mostly Phase 2)

## What the assistant / automation needs

Nothing proprietary is required to **develop or run** the Monte Carlo + PETSc code in this repo. There is no separate “training dataset” for the assistant to operate on.

## Public market data (for realism / calibration)

For **Phase 2** (historical \(\sigma\), optional comparison to quotes), you can use:

- **Yahoo-style OHLCV** — often accessed in Python via **`yfinance`** (unofficial API; fine for coursework experiments, not for production compliance). Good for **log-return volatility** estimates on an equity index or single name.
- **Option chains** — `yfinance` can return limited chain snapshots; quality and fields vary. For a serious smile study you may prefer a vendor or exchange documented feed.

## Tier A: step-by-step (get OHLCV + run calibration)

Do this **once** on your machine inside the project repo (with your PETSc Python or a venv that has `yfinance` for the fetch step only).

### Step 1 — Pick the underlying

Choose one liquid ticker for a clean story, e.g. **SPY** (S&P 500 ETF) or **QQQ**. One file = one ticker.

### Step 2 — Pick the history window

- **Long sample (full-sample vol):** e.g. `--start 2018-01-01 --end 2024-12-31` (many trading days).  
- **Recent vol only:** keep a long download but pass **`--window 252`** to `phase2_calibrate_and_price.py` (last ~1 year of *rows* in the CSV after sort).

### Step 3 — Install fetch dependencies (not in core `requirements.txt`)

From the project root:

```bash
pip install -r requirements-phase2.txt
```

(Uses `yfinance` + `pandas` only for downloading; core Monte Carlo still uses `requirements.txt`.)

### Step 4 — Download daily OHLCV to `data/raw/`

Example:

```bash
python scripts/fetch_yahoo_ohlcv.py \
  --ticker SPY \
  --start 2018-01-01 \
  --end 2024-12-31 \
  --out data/raw/SPY.csv
```

**Check:** open `data/raw/SPY.csv` and confirm you see columns including **Date** and **Adj Close** (or **Close**).

### Step 5 — Pick your option “scenario” (write these down for the report)

Decide and record:

| Input | Meaning | Example |
| --- | --- | --- |
| \(S_0\) | Spot for pricing | Often **last row’s Adj Close** in the CSV (script default), or a specific print → use `--s0` |
| \(K\) | Strike | e.g. **450** for an SPY call near the money (adjust to your date’s level) |
| \(T\) | Years to expiry | e.g. **0.25** (~3 months), **1.0** for one year |
| \(r\) | Annual risk-free, continuous-ish | e.g. **0.04**–**0.05** with a one-line justification |

You will pass **`--strike`**, **`--t`**, **`--r`**, and optionally **`--s0`** into `phase2_calibrate_and_price.py`.

### Step 6 — Run Phase 2 (baseline σ vs calibrated σ from your CSV)

**Full-history vol** (no window):

```bash
PYTHONPATH=src python scripts/phase2_calibrate_and_price.py \
  --csv data/raw/SPY.csv \
  --strike 450 \
  --t 0.25 \
  --r 0.045 \
  --paths 200000 \
  --sigma-baseline 0.20 \
  --out-json output/phase2_SPY_summary.json
```

**Last 252 rows only** for “recent” annualized vol:

```bash
PYTHONPATH=src python scripts/phase2_calibrate_and_price.py \
  --csv data/raw/SPY.csv \
  --window 252 \
  --strike 450 \
  --t 0.25 \
  --r 0.045 \
  --paths 200000 \
  --out-json output/phase2_SPY_recent252.json
```

Omit `--s0` to use the **last close** in the file as spot; add `--s0 480` if you want a fixed print.

### Step 7 — Log the run (optional but good for your registry)

```bash
python scripts/log_experiment_run.py \
  --key phase2_spy_vol \
  --output data/raw/SPY.csv \
  --output output/phase2_SPY_summary.json \
  --notes "Tier A OHLCV + phase2_calibrate_and_price"
```

### Step 8 — Report one sentence of provenance

State: **ticker, source (Yahoo via yfinance), date range, whether vol is full-sample or `--window 252`, and chosen \(S_0,K,T,r\)**.

---

### Phase 2 pipeline (short reference)

Same commands as **Tier A Step 6** above. Offline demo:

```bash
PYTHONPATH=src python scripts/phase2_calibrate_and_price.py --csv data/raw/example_synthetic_ohlcv.csv --paths 8000
```

Respect Yahoo’s terms of service and rate limits; cache CSVs under `data/` and cite the source in the report.

### If you want to hand files to the project

- Drop **CSV OHLCV** (columns like `Date,Open,High,Low,Close,Adj Close,Volume`) under `data/raw/` and document the ticker and date range in `output/registry/experiment_runs.csv` or your report appendix.
- For **option quotes**, a small CSV (strike, expiry, bid, ask, mid, underlying) is enough for a toy comparison to MC. Copy **`data/raw/quotes_template.csv`**, fill numeric rows, then run `PYTHONPATH=src python scripts/compare_mc_to_quotes.py --help`.

No upload to the model is required unless you paste a snippet in chat for debugging.
