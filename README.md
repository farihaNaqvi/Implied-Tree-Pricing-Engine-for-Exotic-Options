# Implied Tree Pricing Engine for Exotic Options

A Java implementation of an implied trinomial tree calibration and pricing engine for exotic options. The project builds an arbitrage-free trinomial lattice from market-observed European option prices, infers state prices and transition probabilities, and uses the calibrated tree to price an American up-and-out put.

## Overview

This project is a Java port of a Python workflow that:

- calibrates an implied trinomial tree using European call/put prices,
- recovers Arrow-Debreu state prices Q[j, i],
- derives forward transition probabilities pu, pm, and pd,
- prices a path-dependent American barrier option on the resulting lattice.

The implementation is self-contained and can run in two modes:

1. Market-data mode via CSV input files for call and put prices.
2. Black-Scholes fallback mode when no CSV files are supplied.

## Key Features

- Implied trinomial tree construction from option prices
- Calibration of state prices and transition probabilities
- Recombining logarithmic tree structure
- American up-and-out barrier put pricing by backward induction
- Black-Scholes validation fallback for sanity checks
- CSV-based market data support for practical calibration workflows

## Mathematical Model

The engine uses a trinomial lattice with log-spaced asset levels and a one-step gross risk-free factor:

- S0 = 100
- r = 5%
- T = 1 year
- N = 4 time steps
- dx = 0.2524

The code computes:

- terminal asset prices on the lattice,
- state prices Q[j, i] from European option prices,
- transition probabilities pu, pm, pd,
- an American up-and-out put value via dynamic programming.

## Repository Contents

- `ImpliedTrinomialTree.java` - main implementation
- `ImpliedTrinomialTree.class` - compiled Java class output
- `.gitignore` - standard repository ignore rules

## Prerequisites

- Java JDK 8 or later

## Build and Run

Compile:

```bash
javac ImpliedTrinomialTree.java
```

Run with Black-Scholes fallback:

```bash
java ImpliedTrinomialTree
```

Run with market CSV data:

```bash
java ImpliedTrinomialTree --call calls.csv --put puts.csv
```

## CSV Input Format

The project expects CSV files containing a grid of European call or put prices in the same indexing convention used by the original spreadsheet.

The code reads each CSV as a matrix of size:

- rows = 2N + 1
- cols = N + 1

For the current configuration, this is:

- 9 rows x 5 columns

Each file should contain prices aligned to the implied tree nodes for each step.

## Example Usage

### Default mode

```bash
javac ImpliedTrinomialTree.java
java ImpliedTrinomialTree
```

This generates synthetic Black-Scholes European option prices and prices an American up-and-out put on the implied tree.

### Market calibration mode

```bash
java ImpliedTrinomialTree --call data/calls.csv --put data/puts.csv
```

This reads actual market prices from CSV files and calibrates the tree to those inputs.

## Output

The program prints:

- the terminal asset-price matrix,
- the state-price matrix,
- the implied up/middle/down transition probability matrices,
- the American up-and-out put price.

## Notes

- The code is educational and research-oriented.
- The fallback Black-Scholes mode is intended as a sanity check rather than a production pricing model.
- CSV mode is useful when calibrating to observed option market data.

## License

This repository does not currently include an explicit license file. If you plan to share or use this code publicly, consider adding an open-source license such as MIT or Apache-2.0.

## Author

This repository appears to be maintained by `farihaNaqvi`.

## Contact / Repository

- GitHub: https://github.com/farihaNaqvi/Implied-Tree-Pricing-Engine-for-Exotic-Options
