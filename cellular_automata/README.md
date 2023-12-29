## Game of Life ##
The evolution of the Game of Life on a lattice is fully deterministic
and obeys the following rules:
- Any live cell with fewer than 2 live neighbours dies
- Any live cell with 2 or 3 live neighbours lives on to the next step
- Any live cell with more than 3 live neighbours dies
- Any dead cell with exactly 3 live neighbours comes to life

The game of life is simulated by the script `game_of_life.py`.
When run, the script will prompt the user to define the lattice size.
The results contained in this folder are for a 50 x 50 lattice.

## SIRS model for epidemic spreading ##
The SIRS model is a stochastic model consisting of agents 
arranged on a 2d lattice. 

These agents can be in 1 of the 3 following states:
- S: susceptible to infection
- I: infected
- R: recovered from infection (and therefore immune)

Agents pass through the available states cyclically: S --> I when a 
susceptible agent comes into contact with an infected one, I --> R
through recovery, and, in the absence of contact with an infected agent,
the immune response will wear off and R --> S.

`sirs.py` simulates this model through a random sequential updating scheme.
Only one site is updated in each timestep, and the updating site is chosen
at random. The SIRS model obeys the following rules:
- S --> I with probability $p_1$ if minimum 1 neighbour of the 
updating site is I; otherwise, it is unchanged.
- I --> R with probability $p_2$.
- R --> S with probability $p_3$.