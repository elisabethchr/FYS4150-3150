# Project 5
FYS4150/FYS3150: Project 5 by Eirik Bratli and Elisabeth Christensen

#### Project description:
This project aims to construct an algorithm for simulating the wealth distribution in a closed economy, when allowing pairs of financial agents to trade money through $10^7$ transactions as Monte Carlo cycles. Here, we take into account different financial restrictions, such as including a saving parameter $\lambda$, parameter $\alpha$ for probability of interaction with neighbouring agents, and the parameter $\gamma$ for the number of times an agent has made transactions. The tails of the wealth distribution can be fitted by the Pareto distribution $w_m \propto m^{-1-\alpha-\lambda-\gamma}$.

#### Files included:

  - main.cpp
  - stockmarked.h/cpp
