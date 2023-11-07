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

# Example usage:
p = 17  # Example prime modulus, p = 1 + 4*2^2
w = 3   # Example primitive root for NTT modulo 17

# Polynomial coefficients
a = [1, 2, 3, 4]  # Example polynomial: 1 + 2x + 3x^2 + 4x^3

# Compute NTT
ntt_result = ntt(a, p, w)
print("NTT result:", ntt_result)

# Compute INTT
intt_result = intt(ntt_result, p, pow(w, p - 2, p))
print("INTT (original coefficients):", intt_result)
