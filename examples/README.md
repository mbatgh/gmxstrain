
# Example for calculation of elastic constants gmxstrain

1. Generate super-cell

2. Generate topology

3. Run MD simulations

./gmxstrain.awk -v id=dps -v nmax=3 -v dd=0.004 -v t0=1000 1 > er-dps 2>&1 &

4. Analyze

