# sex-limited-imperfect-Müllerian-mimicry
Codes and raw data accompanying the manuscript "Sexual difference in defense can drive the evolution of imperfect Müllerian mimicry in the less defended sex"

**Files:**

1. male-limited mimicry.R: Simulation R codes

2. df_FLR.txt.gz: Compressed raw data for simulations with a monandrous mating system

- **tf:** female trait value, ranging from 0 (identical to male value) to 1 (maximally different from male value)
- **M:** proportion of offpring that are male
- **Pm:** predator atatck probability towards males
- **Pf:** predator attack probability towards females
- **Nm:** number of males at equilibrium
- **Nf:** number of females at equilibrium
- **N:** population size at equilibrium (Nm + Mf)
- **lag:** mismatch between male preference and female trait (NA in this dataset)
- **MPmPf:** indexing variable for calculating relative fitness and graphing purposes
- **w:** relative fitness

3. df_notFLR.txt.gz: Compressed raw data for simulations with a polyandrous mating system. Columns are the same as df_FLR.txt

4. df_lag.txt.gz: Compressed raw data for simulations with mismatching male preference and female trait under a monandrous mating system and even offspring sex ratio. Columns are the same as df_FLR.txt, with the following two exceptions:

- **lag:** lag values in this dataset ranges from 0 to 0.1
- **lagPmPf:** indexing variable for calculating relative fitness and graphing purposes


