## Code refactor

### Loglikelihood

Current implementation uses double for loop to compute. I have so far implemented a single loop. I need to test and do speed comparisons.
	- The R implementaion slightly between the package


### Optimization algorithm 

Current implementation uses proximal gradient descent which inferior compared to cyclic coordinate descent. My main challenge with this approach is computing the residuals after every iteration along the coordinates for survival data. Still trying to move this forward.
