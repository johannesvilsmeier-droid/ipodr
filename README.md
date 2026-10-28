# ipodr

Option-implied distress measures for single stocks, together with the ranges
of values that the option quotes actually determine.

The package implements:

* a penalised maximum-entropy version of the option-implied probability of
  default (iPoD; Capuano 2008, Vilsmeier 2016), with barrier profiles and
  barrier-selection rules;
* barrier-free tail probabilities \(Q(S_T/F < a)\) and the anchored distress
  probability \(Q(S_T < \theta \cdot \text{peak})\);
* sharp model-free bounds on tail probabilities, call prices and at-the-money
  implied volatility from bid and ask quotes (Dillschneider 2026), and a
  linear-programming version for general payoffs;
* benchmarks (Black implied volatility, lognormal tails, drawdown, distance
  to threshold) and evaluation tools (event windows, bank-cluster bootstrap
  AUC comparisons, AUC bounds, first-alarm dates).

The probability of zero equity is never identified by option quotes: its
lower bound is zero whenever the lowest quoted strike is positive. Report
`tail_bounds(chain, 0)` alongside any PoD.

## Example

```r
library(ipodr)
q  <- sim_quotes(sigma = 0.4, pod = 0.01)   # long data frame: one row per quote
ch <- ipod_chain(q)                          # parity forward, OTM calls and puts
fit <- ipod_fit(ch, d = 0.1)
tail_prob(fit, c(0.5, 0.8))
tail_bounds(ch, c(0, 0.5, 0.8))
atm_iv(ch)
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("<user>/ipodr")
```
