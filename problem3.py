import argparse
import numpy as np
import copy

"""
Задача 3

Дана матрица А, размера MxN, заполненная произвольными действительными числами.
Необходимо написать программу, которая находит в этой матрице квадратную
подматрицу с максимальной суммой элементов.

Алгоритм: сложность N^3
Сначала считаем за N^2 сумму во всех полосках k*1.
Далее для любого размера квадрата максимальный квадрат можно найти за N^2.
И внешний цикл по размерам квадратов.
"""

def get_k_strips(arr, k):
    strips = np.zeros(arr.shape)
    for j in range(arr.shape[0]):
        sum_strip=0
        for i in range(k):
            sum_strip += arr[i,j]
            strips[0,j] = sum_strip
        for i in range(1, arr.shape[1]-k+1):
            sum_strip += (arr[i+k-1,j] - arr[i-1,j]);
            strips[i, j] = sum_strip;
    return strips


def max_sub_square(arr, k):
    strip_sum = get_k_strips(arr, k)
    max_sum = float('-inf')
    n = arr.shape[0]
    for i in range(n-k+1):
        temp_sum = 0
        for j in range(k):
            temp_sum += strip_sum[i,j]

        if temp_sum > max_sum:
            max_sum = temp_sum
            pos = [i,0]

        for j in range(1, n-k+1):
            temp_sum += (strip_sum[i, j+k-1] - strip_sum[i, j-1])

            if temp_sum > max_sum:
                max_sum = temp_sum
                pos = [i,j]

    i, j = pos
    return max_sum, [[i,j],[i+k,j+k]],


def solution(arr):
    max_sum = float('-inf')
    max_arr = None
    for n in range(1, arr.shape[0]+1):
        temp_sum, temp_arr = max_sub_square(arr, n)
        if temp_sum > max_sum:
            max_sum = temp_sum
            max_arr = temp_arr

    return max_sum, max_arr


def print_array(arr, max_ind):
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if ((i >= max_ind[0][0]) and (i < max_ind[1][0]) and \
                (j >= max_ind[0][1]) and (j < max_ind[1][1])):
                print('\033[92m'+str(int(arr[i][j])).zfill(2)+' ', end='')
            else:
                print('\033[0m'+str(int(arr[i][j])).zfill(2)+' ', end='')
        print('\n')


def read_array(file_path):
    with open(file_path) as f:
        content = f.readlines()
        arr = np.zeros((len(content), len(content[0].split(' '))))
        for i,line in enumerate(content):
            arr[i,] = [int(el) for el in line.split(' ')]
    return arr

parser = argparse.ArgumentParser(description='Finding subarray with maximum sum element.')
parser.add_argument('--file', type=str, default='./example_data/array_problem3.txt', help='path to .txt array file')
args = parser.parse_args()

arr = np.array(read_array(args.file))
solution = solution(arr)
print('Maximum area: ', solution[0])
print('\n')
print_array(arr, solution[1])
