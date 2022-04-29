import matplotlib.pyplot as plt
import subprocess
import numpy as np
import re

np.random.seed(10)

def generate_input(file_name, n):
	'''
	Generates input for cpp program
	'''
	file = open(file_name, 'w')
	file.write(f'{n}\n')
	file.close()
	v = np.random.randint(1, 100, size=(1,n))
	with open(file_name, 'a') as f:
		np.savetxt(f, v, fmt='%d', delimiter=' ')

def empty_output_file(file_name):
	'''
	Empties the output file
	'''
	with open(file_name, 'w') as f:
		pass

if __name__ == '__main__':
	# Empty contents of output file
	empty_output_file('expt8_output.txt')
	
	# Generate inputs and run the c program for all process sizes
	for i in range(16, 10000000, 100000):
		generate_input("input.txt", i)
		for j in ['2','4','8']:
			subprocess.run(["mpirun", "--oversubscribe", "-n", j, "expt8"])
	
	# Get output from file
	with open("expt8_output.txt", 'r') as f:
		data = f.readlines()
	
	pattern = re.compile(r'\d*\.?\d*')
	X, Y = [], []
	for line in data:
		_, x, y = list(filter(None, pattern.findall(line)))
		X.append(float(x))
		Y.append(float(y))
	
	# Plotting the results
	plt.plot(X[::3], Y[::3], label='2')
	plt.plot(X[1::3], Y[1::3], label='4')
	plt.plot(X[2::3], Y[2::3], label='8')
	plt.legend()
	plt.xlabel('Size of Matrix')
	plt.ylabel('Time of execution')
	plt.title('Time of execution for number of processes')
	plt.savefig('1.png')

	# empty_output_file("expt6_output.txt")
	# n = int(input("Enter matrix size: "))
	# generate_input("matrix.txt", n);
	# for j in range(1, 11):
		# subprocess.run(["mpirun", "--oversubscribe", "-n", str(j), "expt6"])
	
	# with open("expt6_output.txt", 'r') as f:
		# data = f.readlines()
	
	# pattern = re.compile(r'\d*\.?\d*')
	# X, Y = [], []
	# for line in data:
		# p, x, y = list(filter(None, pattern.findall(line)))
		# X.append(float(p))
		# Y.append(float(y))
	# plt.plot(X, Y)
	# plt.title(f'No of processes vs. time for matrix size {n}')
	# plt.savefig('2.png')