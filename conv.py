from NTT_2D import ntt_2d, intt_2d


# Define two 4x4 matrices to convolve
matrix_a = [
    [1, 2, 0, 0],
    [3, 4, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0]
]
matrix_b = [
    [0, 0, 0, 0],
    [0, 5, 6, 0],
    [0, 7, 8, 0],
    [0, 0, 0, 0]
]

def convolve_2d(matrix_a, matrix_b, p):
    """
    Perform 2D convolution between two matrices of the same dimensions.
    :param matrix_a: First input matrix.
    :param matrix_b: Second input matrix.
    :param p: prime number for modulo operation to fit in NTT domain
    :return: Resultant matrix after convolution.
    """
    n = len(matrix_a)  # Assuming square matrices of the same size
    result = [[0] * n for _ in range(n)]

    # Flip matrix_b for convolution
    matrix_b_flipped = [row[::-1] for row in matrix_b[::-1]]

    for i in range(n):
        for j in range(n):
            sum = 0  # Accumulate result in sum
            for k in range(n):
                for l in range(n):
                    # Wrap around the indices for convolution
                    i_k = (i - k) % n
                    j_l = (j - l) % n
                    sum += matrix_a[i_k][j_l] * matrix_b_flipped[k][l]
                    sum %= p  # Modulo p for NTT compatibility
            result[i][j] = sum

    return result


# Convolution in the spatial domain
convolution_spatial = convolve_2d(matrix_a, matrix_b)
print("Convolution (Spatial Domain):")
for row in convolution_spatial:
    print(row)

# NTT of both matrices
ntt_matrix_a = ntt_2d(matrix_a, p, w)
ntt_matrix_b = ntt_2d(matrix_b, p, w)

# Element-wise multiplication in frequency domain
ntt_product = [[0]*4 for _ in range(4)]
for i in range(4):
    for j in range(4):
        ntt_product[i][j] = (ntt_matrix_a[i][j] * ntt_matrix_b[i][j]) % p

# Inverse NTT to get the convolution result back in spatial domain
convolution_frequency = intt_2d(ntt_product, p, w_inv)

print("\nConvolution (Frequency Domain):")
for row in convolution_frequency:
    print(row)
