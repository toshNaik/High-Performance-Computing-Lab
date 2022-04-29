import matplotlib.pyplot as plt
import subprocess
import numpy as np
import re

np.random.seed(10)

def generate_input(file_name, no_rows, no_columns):
	'''
	Generates input for cpp program
	'''
	file = open(file_name, 'w')
	file.write(f'{no_rows} {no_columns}\n')
	file.close()
	v = np.random.randint(1, 100, size=(no_rows, no_columns))
	with open(file_name, 'a') as f:
		np.savetxt(f, v, fmt='%d', delimiter=' ')

def empty_output_file():
	'''
	Empties the output file
	'''
	with open('output.txt', 'w') as f:
		pass

if __name__ == '__main__':
	# Empty contents of output file
	empty_output_file()
	
	# Generate inputs and run the c program for all process sizes
	for i in range(100, 2000, 300):
		generate_input("matrix1.txt", i, i)
		generate_input("matrix2.txt", i, i)	
		for j in ['1','2','3','4']:
			subprocess.run(["mpirun", "-use-hwthread-cpus", "-n", j, "expt5", "0"])
	
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
	plt.savefig('3.png')
