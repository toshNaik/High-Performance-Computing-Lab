# Python file for plotting

import matplotlib.pyplot as plt
import subprocess
import numpy as np
import re

np.random.seed(10)

def generate_input(no_rows):
	file = open('input.txt', 'w')
	file.write(f'{no_rows}\n')
	file.close()
	v = np.random.randint(1, 100, size=(no_rows, no_rows))
	v = np.where(v > 80, -1, v)
	np.fill_diagonal(v, 0)
	with open('input.txt', 'a') as f:
		np.savetxt(f, v, fmt='%d', delimiter=' ')

def empty_output_file():
	with open('output.txt', 'w') as f:
		pass

if __name__ == '__main__':
	# Empty contents of output file
	empty_output_file()
	
	# Generate inputs and run the c program for all process sizes
	for i in range(4, 11):
		generate_input(i)
		for j in ['1','2','3','4']:
			subprocess.run(["mpirun", "-use-hwthread-cpus", "-n", j, "expt2", "0"])
	
	# Get output from file
	with open("output.txt", 'r') as f:
		data = f.readlines()
	
	pattern = re.compile(r'\d*\.?\d*')
	X, Y = [], []
	for line in data:
		_, x, y = list(filter(None, pattern.findall(line)))
		X.append(float(x))
		Y.append(float(y))
	
	# Plotting the results
	plt.plot(X[::4], Y[::4], label='1')
	plt.plot(X[1::4], Y[1::4], label='2')
	plt.plot(X[2::4], Y[2::4], label='3')
	plt.plot(X[3::4], Y[3::4], label='4')
	plt.legend()
	plt.xlabel('Size of Matrix')
	plt.ylabel('Time of execution')
	plt.title('Time of execution for number of processes')
	plt.savefig('4.png')