## Code refactor

### Loglikelihood

Current implementation uses double for loop to compute. I have so far implemented a single loop. I need to test and do speed comparisons.
	- The R implementaion slightly between the package

### Optimization algorithm 

Current implementation uses proximal gradient descent which inferior compared to cyclic coordinate descent. My main challenge with this approach is computing the residuals after every iteration along the coordinates for survival data. Still trying to move this forward.

2023 Aug 09 (Wed)
=================

Reading https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4824408/ (a key glmnet paper)
* Need to agree on what is meant by “warm start”
* What are their “risk-set update” tricks?
* How well do we understand the duality between penalties and constraints?
	* The set of penalized and constrained solutions is probably the same set
	* Here they seem to switch just so they can switch back!
* What is z and what is our logical framework for the NR step?
* LESS important: do we understand what is “convenient” about 2/n?

