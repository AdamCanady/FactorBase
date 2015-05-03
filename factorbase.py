def b_smooth_factors(i, B):
    ''' Returns the factorization of i in terms of the factors in B,
        or False if no such factorization exists.
    '''

def factor_base(n):
    ''' Returns a factor of n. '''
    def pick_factor_list(k):
        p = Primes()
        B = [-1]
        for i in xrange(k):
            B.append(p.next(B[-1]))
        return B

    def pick_candidates(k):
        return xrange(n - k, n + k)

    b_factor_list = pick_factor_list(20)
    b_candidates = pick_candidates(1000)

    # b_smooth is a list of tuples (b, factorization of b)
    b_smooth = []
    for b in b_candidates:
        b_squared = Mod(b, n)^2.centerlift()
        b_factors = b_smooth_factors(b_squared, b_factor_list)
        if b_factors:
            b_smooth.append(b, b_factors)

    A = Matrix(GF(2), [[(1 if factor in b_factors else 0) for factor in B] for (b, b_factors) in b_smooth])

    A_transpose = A.transpose()

    kernel_basis = A.right_kernel().basis() # set of vectors that span the space of ker(A)
    for kernel_vector in kernel_basis:
        rows = [b_tuple for i, b_tuple in enumerate(b_smooth) if kernel_vector[i]]
        x = reduce(prod, map(lambda x: x[0], rows))
        y = reduce(prod, map(lambda x: x[1], rows)) # Todo: take product of factor objects
        if x != y:
            # We win!
            return gcd(n, x - y)
    return 1 # but really add to B and keep going.
