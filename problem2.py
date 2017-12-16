"""
Задача 2
Дана матрица А, размера MxN, заполненная произвольными действительными числами.
Необходимо написать программу, которая находит в этой матрице квадратную
подматрицу с максимальной суммой элементов.

Решение:
Динамическое программирование, алгоритм приведен ниже.
"""

import argparse
import numpy as np
import copy

"""
Algorithm:
Step 1: Find maximum area for row[0]
Step 2:
    for each row in 1 to N - 1
        for each column in that row
            if A[row][column] == 1
              update A[row][column] with
                A[row][column] += A[row - 1][column]
    find area for that row
    and update maximum area so far
"""

def histogram_max_area(height):
    increasing, area, i = [], 0, 0
    coords = []
    while i <= len(height):
        if not increasing or (i < len(height) and height[i] > height[increasing[-1]]):
            increasing.append(i)
            i += 1
        else:
            last = increasing.pop()
            if not increasing:
                p_area = height[last] * i
                if p_area > area:
                    area = p_area
                    coords = [0, i, height[last]]
            else:
                p_area = height[last] * (i - increasing[-1] - 1 )
                if p_area > area:
                    area = p_area
                    coords = [increasing[-1]+1, i, height[last]]
    return area, coords


def max_area(arr):
    A = copy.deepcopy(arr)
    result,max_ind = histogram_max_area(A[0])
    argmax = 0
    for i in range(1, A.shape[0]):
        for j in range(0, A.shape[1]):
            if A[i][j] == 1:
                A[i][j] += A[i-1][j]

            area, coords = histogram_max_area(A[i])
            if area > result:
                argmax = i
                result = area
                max_ind = [(i - coords[2]+1, coords[0]),(i+1, coords[1])]
    return A, result, max_ind


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


parser = argparse.ArgumentParser(description='Finding maximum subarray with all ones.')
parser.add_argument('--file', type=str, default='./example_data/array_problem2.txt', help='path to .txt array file')
args = parser.parse_args()

arr = np.array(read_array(args.file))
A, result, max_ind = max_area(arr)
print_array(arr, max_ind)
print("Maximum area: ", result)
