from scipy import weave
from scipy.weave import converters
from numpy import random
from time import time;

def mySum(a):
	n = int(len(a))

	# Start scipy.weave
	code = """
	double counter = 0;
	for (int i = 0; i<n; i++ ) {
		counter = counter + a(i);
	}
	return_val = counter;
	"""
	ret = weave.inline(code, ['a', 'n'], type_converters = converters.blitz, compiler='gcc')
	# Stop scipy.weave

	return ret;

numbers = random.random(10000000);
t_0 = time();
sum_numbers = mySum(numbers)
t_1 = time()
sum_numbers2 = sum(numbers)
t_2 = time()

t_spent_1 = t_1-t_0;
t_spent_2 = t_2-t_1;

print "Speedup factor %f " % (t_spent_2/t_spent_1)

print sum_numbers
print sum_numbers2