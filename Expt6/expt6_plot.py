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
	v = np.random.randint(1, 100, size=(n,n))
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
	empty_output_file('expt6_output.txt')
	
	choice = int(input("Enter choice:\n1. Custom input\n2. Plot graph"))
	if choice == 1:
		array = np.loadtxt('matrix.txt', dtype='f', skiprows=1)
		if(np.linalg.det(array[:, :-1]) == 0):
			print('Error determinant is zero')
			exit()
		j = input('Enter number of processes: ')
		subprocess.run(["mpirun", "--oversubscribe", "-n", j, "expt6", "0"])
	if choice == 2:
		# # Generate inputs and run the c program for all process sizes
		# for i in range(100, 2000, 300):
			# generate_input("matrix.txt", i)
			# for j in ['1','2','3','4']:
				# subprocess.run(["mpirun", "--oversubscribe", "-n", j, "expt6", "0"])
		
		# # Get output from file
		# with open("expt6_output.txt", 'r') as f:
			# data = f.readlines()
		
		# pattern = re.compile(r'\d*\.?\d*')
		# X, Y = [], []
		# for line in data:
			# _, x, y = list(filter(None, pattern.findall(line)))
			# X.append(float(x))
			# Y.append(float(y))
		
		# # Plotting the results
		# plt.plot(X[::4], Y[::4], label='1')
		# plt.plot(X[1::4], Y[1::4], label='2')
		# plt.plot(X[2::4], Y[2::4], label='3')
		# plt.plot(X[3::4], Y[3::4], label='4')
		# plt.legend()
		# plt.xlabel('Size of Matrix')
		# plt.ylabel('Time of execution')
		# plt.title('Time of execution for number of processes')
		# plt.savefig('1.png')

		empty_output_file("expt6_output.txt")
		n = int(input("Enter matrix size: "))
		generate_input("matrix.txt", n);
		for j in range(1, 11):
			subprocess.run(["mpirun", "--oversubscribe", "-n", str(j), "expt6"])
		
		with open("expt6_output.txt", 'r') as f:
			data = f.readlines()
		
		pattern = re.compile(r'\d*\.?\d*')
		X, Y = [], []
		for line in data:
			p, x, y = list(filter(None, pattern.findall(line)))
			X.append(float(p))
			Y.append(float(y))
		plt.plot(X, Y)
		plt.title(f'No of processes vs. time for matrix size {n}')
		plt.savefig('2.png')