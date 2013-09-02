from ParachuteProblem import ParachuteProblem
import matplotlib.pyplot as plt
import sys
# Take command line args
def read_command_line():
	if (len(sys.argv) < 6):
		print 'Usage %s I b m T dt' % sys.argv[0]; sys.exit(1)
	I 	= 	float(sys.argv[1])
	b 	= 	float(sys.argv[2])
	m 	= 	float(sys.argv[3])
	T 	= 	float(sys.argv[4])
	dt 	= 	float(sys.argv[5])
	return I, b, m, T, dt

if __name__=="__main__":
	read_command_line()
	# Do calculations
	I, b, m, T, dt = read_command_line()

	parachuter = ParachuteProblem(m, b)
	parachuter.set_initial_condition(I)
	parachuter.plot(T, dt);