import argparse
import numpy as np
import copy

"""
Algorithm:
ToDo
"""

def kadane_1d(arr):
    in_arr = [-1]+arr
    position = [0]*len(in_arr)
    sum_arr = [-1]*len(in_arr)
    max_sum = 0
    max_sum_ind = 0
    for i in range(1, len(in_arr)):
        if in_arr[i] > sum_arr[i-1]+in_arr[i]:
            #Начинаем новую цепочку
            position[i] = i-1
            sum_arr[i] = in_arr[i]
            if sum_arr[i] > max_sum:
                max_sum = sum_arr[i]
                max_sum_ind = i-1
        else:
            #Продолжаем текущую цепочку
            position[i] = position[i-1]
            sum_arr[i] = in_arr[i] + sum_arr[i-1]
            if sum_arr[i] > max_sum:
                max_sum = sum_arr[i]
                max_sum_ind = i-1
    return max_sum, [position[max_sum_ind], max_sum_ind]

print(kadane_1d([1,-3,2,1,-1]))

def prefix_trick(matrix):
    current_sum = 0
    max_sum = 0
    max_left = 0
    max_right = 0
    max_up = 0
    max_down = 0

    for l in range(matrix.shape[0]):
        buff = np.zeros(matrix.shape[1])
        for r in range(l, matrix.shape[1]):
            buff = buff + matrix[r,]
            current_sum, c_c = kedane_1d(buff)
            if current_sum > max_sum:
                max_sum = current_sum
                max_left, max_right, max_up, max_down = l, r, c_coords[0], c_coords[1]

    return max_sum, [(max_left, max_up),(max_right, max_down)]


def prepro(A):
    p_array = np.zeros((4,A.shape[0],A.shape[1]))
    p_array[0,0,:] = A[0,:]
    p_array[1,:,A.shape[1]-1] = A[:,A.shape[1]-1]
    p_array[2,A.shape[0]-1,:] = A[A.shape[0]-1,:]
    p_array[3,:,A.shape[1]-1] = A[:,A.shape[1]-1]
    for i in range(1,A.shape[0]):
            p_array[0,i,:] = p_array[0,i-1,:]+A[i,:]
    for j in range(1,A.shape[1]):
            p_array[0,:,j] = p_array[0,:,j-1]+p_array[0,:,j]
    for j in range(A.shape[1]-1,0,-1):
            p_array[1,:,j-1] = p_array[1,:,j]+A[:,j-1]
    for i in range(1,A.shape[0]):
            p_array[1,i,:] = p_array[1,i-1,:]+p_array[1,i,:]
    for i in range(A.shape[0]-1,0,-1):
        p_array[2,i-1,:] = p_array[2,i,:]+A[i-1,:]
    for j in range(1,A.shape[1]):
            p_array[2,:,j] = p_array[2,:,j-1]+p_array[2,:,j]
    for j in range(A.shape[1]-1,0,-1):
            p_array[3,:,j-1] = p_array[3,:,j]+A[:,j-1]
    for i in range(A.shape[0]-1,0,-1):
            p_array[3,i-1,:] = p_array[3,i-1,:]+p_array[3,i,:]
    return p_array

class Prepro():
    def __init__(self,A):
        self.p_array = prepro(A)
        self.sum = np.sum(A)

    def get_area(self,t1,t2,a,b):
        area = 0.0
        x1,y1 = t1
        x2,y2 = t2
        if x1>0 and y2<b:
            area += self.p_array[0,x1-1,y2]
        if x2<a and y2<b:
            area += self.p_array[1,x2,y2]
        if x1>0 and y1>0:
            area += self.p_array[2,x1-1,y1-1]
        if x2<a and y1>0:
            area += self.p_array[3,x2,y1-1]
        area = self.sum - area
        return area

a = Prepro(A)
a.get_area((0,0),(2,2),3,3)

def get_sum():
    pass

def dummy_solution():
    pass



def print_array(arr, max_ind):
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if ((i >= max_ind[0][0]) and (i < max_ind[1][0]) and \
                (j >= max_ind[0][1]) and (j < max_ind[1][1])):
                print('\033[92m'+str(int(arr[i][j]))+' ', end='')
            else:
                print('\033[0m'+str(int(arr[i][j]))+' ', end='')
        print('\n')

def read_array(file_path):
    with open(file_path) as f:
        content = f.readlines()
        arr = np.zeros((len(content), len(content[0].split(' '))))
        for i,line in enumerate(content):
            arr[i,] = [int(el) for el in line.split(' ')]
    return arr

parser = argparse.ArgumentParser(description='Finding subarray with maximum sum element.')
parser.add_argument('--file', type=str, default='./example_data/array_problem2.txt', help='path to .txt array file')
args = parser.parse_args()

arr = np.array(read_array(args.file))
A, result, max_ind = max_area(arr)
print_array(arr, max_ind)
print("Maximum area: ", result)
