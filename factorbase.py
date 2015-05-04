def b_smooth_factors(i, list_of_primes):
    ''' Returns the factorization of i in terms of the factors in list_of_primes,
        or False if no such factorization exists.
    '''
    list_of_factors = [] # [2, 2, 3, ...]

    cur_prime_i = 0
    while cur_prime_i < len(list_of_primes):
        cur_prime = list_of_primes[cur_prime_i]
        if cur_prime == -1:
            if i < 0:
                list_of_factors.append(cur_prime)
                i = -i
            cur_prime_i += 1
            continue
        if i % cur_prime == 0:
            list_of_factors.append(cur_prime)
            i /= cur_prime
        else:
            cur_prime_i += 1

    if i == 1:
        return list_of_factors
    else:
        return False

def factor_base(n, quad_sieve = False):
    ''' Returns a factor of n. '''
    def pick_prime_list(k):
        p = Primes()
        B = [-1]
        for i in xrange(k):
            B.append(p.next(B[-1]))
        return B

    def pick_candidates(k):
        lower = round(n^(0.5)) + 1
        upper = min(round(2 * n^(0.5)), lower + k)

        i = lower
        while i < upper:
            yield i
            i += 1
            
    power_limit = 10
    number_of_primes = 20
    candidate_size = 2000

    while True:
        prime_list = pick_prime_list(number_of_primes)
        if quad_sieve:
            prime_list = filter(lambda x: x in [-1, 2] or legendre_symbol(n, x) != -1, prime_list)
        print "primes list", prime_list
        b_candidates = pick_candidates(2000)
        
        #print "prime_list:", prime_list
        #print "b_candidates:", b_candidates
    
        # b_smooth is a list of tuples (b, factorization of b)
        b_smooth = []
        if quad_sieve: 
            b_sieve = [(b, b^2 - n, []) for b in b_candidates]
            for a_prime in prime_list:
                if a_prime < 0: continue
                
                prime_power = 1
                while prime_power < power_limit:
                    prime_mod = a_prime ^ prime_power
                    prime_power += 1
                    roots = Mod(n, prime_mod).nth_root(2, all=True)
                    if len(roots) < 1: break
                    for root in roots:
                        current_index = int((root - b_sieve[0][0]) % prime_mod)
                        while current_index < len(b_sieve):
                            (b_orig, b_remaining, b_factors) = b_sieve[current_index]
                            assert b_remaining % a_prime == 0
                            if b_remaining % a_prime == 0:
                                b_sieve[current_index] = (b_orig, b_remaining / a_prime, b_factors + [a_prime])
                            current_index += int(prime_mod)
            # print b_sieve
            print "len b_sieve", len(b_sieve)
            print "b_sieve", b_sieve

            b_smooth = [(b_orig, b_factors)
                        for (b_orig, b_remaining, b_factors) in b_sieve if b_remaining == 1]
        else:
            for b in b_candidates:
                #print b
                b_squared = (Mod(b, n)^2).centerlift()
                b_factors = b_smooth_factors(b_squared, prime_list)
                if b_factors:
                    b_smooth.append((b, b_factors))
                
        print "b_smooth:", b_smooth
        print "len bsmooth", len(b_smooth)
    
        A = Matrix(GF(2), [[(1 if factor in b_factors else 0) for factor in prime_list] for (b, b_factors) in b_smooth])
    
        A_transpose = A.transpose()
    
        kernel_basis = A_transpose.right_kernel().basis() # set of vectors that span the space of ker(A)
        print "kernel basis length", len(kernel_basis)
        for kernel_vector in kernel_basis:
            rows = [b_tuple for i, b_tuple in enumerate(b_smooth) if kernel_vector[i]]
    
            x = prod(map(lambda x: x[0], rows))
            y = prod(map(lambda x: prod(x[1]), rows)) # Todo: take product of factor objects
            if x != y:
                # We win!
                gcd_result = gcd(n, x - y)
                if gcd_result != 1:
                    return gcd_result
        candidate_size *= 2
        number_of_primes *= 2
        # power_limit *= 2
