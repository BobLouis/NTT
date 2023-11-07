def ntt(a, p, w):
    """
    Perform the Number Theoretic Transform (NTT).
    :param a: list of integers (polynomial coefficients)
    :param p: prime modulus of the form p = c*2^k + 1
    :param w: primitive root of the unity modulo p
    :return: list of transformed coefficients
    """
    n = len(a)
    if n == 1:
        return a
    even = ntt(a[0::2], p, pow(w, 2, p))
    odd = ntt(a[1::2], p, pow(w, 2, p))
    t = 1
    y = [0] * n
    for i in range(n // 2):
        y[i] = (even[i] + t * odd[i]) % p
        y[i + n // 2] = (even[i] - t * odd[i]) % p
        t = t * w % p
    return y

def intt(y, p, w_inv):
    """
    Perform the Inverse Number Theoretic Transform (INTT).
    :param y: list of NTT-transformed coefficients
    :param p: prime modulus of the form p = c*2^k + 1
    :param w_inv: inverse of the primitive root of the unity modulo p
    :return: list of original polynomial coefficients
    """
    n = len(y)
    y_inv = ntt(y, p, w_inv)
    n_inv = pow(n, p - 2, p)  # Modular multiplicative inverse of n
    return [(val * n_inv) % p for val in y_inv]




def ntt_2d(matrix, p, w):
    """
    Perform the 2D Number Theoretic Transform (NTT) on a matrix.
    :param matrix: 2D list of integers (matrix of polynomial coefficients)
    :param p: prime modulus of the form p = c*2^k + 1
    :param w: primitive root of the unity modulo p
    :return: 2D list of transformed coefficients
    """
    # Apply NTT to each row
    transformed_rows = [ntt(row, p, w) for row in matrix]

    # Apply NTT to each column by transposing, transforming, and transposing back
    transformed_columns = list(map(list, zip(*transformed_rows)))  # Transpose
    transformed_columns = [ntt(col, p, w) for col in transformed_columns]
    transformed_matrix = list(map(list, zip(*transformed_columns)))  # Transpose back

    return transformed_matrix

def intt_2d(transformed_matrix, p, w_inv):
    """
    Perform the 2D Inverse Number Theoretic Transform (INTT) on a matrix.
    :param transformed_matrix: 2D list of NTT-transformed coefficients
    :param p: prime modulus of the form p = c*2^k + 1
    :param w_inv: inverse of the primitive root of the unity modulo p
    :return: 2D list of original polynomial coefficients
    """
    # Apply INTT to each column by transposing, transforming, and transposing back
    inv_transformed_columns = list(map(list, zip(*transformed_matrix)))  # Transpose
    inv_transformed_columns = [intt(col, p, w_inv) for col in inv_transformed_columns]
    inv_transformed_rows = list(map(list, zip(*inv_transformed_columns)))  # Transpose back

    # Apply INTT to each row
    original_matrix = [intt(row, p, w_inv) for row in inv_transformed_rows]

    return original_matrix

# Example usage with a 4x4 matrix
p = 17  # Example prime modulus, p = 1 + 4*2^2
w = 3   # Example primitive root for NTT modulo 17
w_inv = pow(w, p - 2, p)  # Modular multiplicative inverse of w modulo p

# 4x4 matrix of polynomial coefficients
matrix = [
    [1, 2, 3, 4],
    [5, 6, 7, 8],
    [9, 10, 11, 12],
    [13, 14, 15, 16]
]

# Compute 2D NTT
ntt_2d_result = ntt_2d(matrix, p, w)
print("2D NTT result:")
for row in ntt_2d_result:
    print(row)

# Compute 2D INTT
intt_2d_result = intt_2d(ntt_2d_result, p, w_inv)
print("\n2D INTT (original matrix):")
for row in intt_2d_result:
    print(row)
