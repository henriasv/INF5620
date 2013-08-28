from ParachuteProblem import ParachuteProblem
import matplotlib.pyplot as plt
# Take command line args

# Do calculations
T = 10
dt = 0.01

parachuter = ParachuteProblem(80.0, 1.5)
parachuter.set_initial_condition(100)
t, u = parachuter.solve(T, dt)

# Plot results
plt.plot(t, u)
plt.show()



