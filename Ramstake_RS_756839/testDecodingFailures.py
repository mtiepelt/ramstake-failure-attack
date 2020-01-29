import sys

n = 756839
p = 2**n - 1


def readEstimations(fname):
    """
        Target file:
        g, a, b, c, d (first line is text, rest is numbers)
    """
    
    with open(fname, 'r') as f:
        # Skip g, a, b, c, d, h
        line = f.readline()
        

        for line in f:
            str_g, str_a, str_b, str_c, str_d, str_H = line.split(' ')

            G = int(str_g)
            
            a = int(str_a)
            b = int(str_b)
            c = int(str_c)
            d = int(str_d)
            
            H = int(str_H)

            if ((a * G + b) % p) != H:
                print "Cant not reconstruct pk!"
            else:
                print "PK successful!" 

            sa = ((a * G + b) * c) % p
            sb = ((c * G + d) * a) % p

            hex_sa = list(reversed(hex(sa)))
            hex_sb = list(reversed(hex(sb)))

            num_errors = 0

            for pos in range(255):
                if (hex_sa[2*pos] != hex_sb[2*pos]) or hex_sa[2*pos+1] != hex_sb[2*pos+1]:
                    num_errors += 1

            if num_errors < 110:
                print "Not enough errors: " + str(num_errors)
            else:
                print "All good!"

readEstimations(sys.argv[1])