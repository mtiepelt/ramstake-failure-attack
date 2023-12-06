import math

n = 756839
w = 128

# Bytes per codeword
num_bytes = 255 

prec = 3


# Single bit is zero
p_bit_zero = 1 - w/n

# Byte is zero, if all bits are zero 
p_byte_zero = p_bit_zero**8

# Byte is non-zero, if not all bits zero
p_byte_non_zero = 1 - p_byte_zero


def prob_any_x_bytes_non_zero(p_byte_non_zero, num_bytes_non_zero, num_bytes):
	p = 1
	for i in range(num_bytes_non_zero):
		p *= (num_bytes - i) * p_byte_non_zero 

	p *= (1-p_byte_non_zero)**(num_bytes - num_bytes_non_zero)
	return p 


##
## PROBABILITY OF COLLISION IN 2^i b_i c
##
print(f"Pr[bit is zero] = {p_bit_zero} --> {round(p_bit_zero * 100, prec)}%")
print(f"Pr[byte is zero] = {p_byte_zero} --> {round(p_byte_zero * 100, prec)}%")
print(f"Pr[byte non zero] = {1 - p_bit_zero} --> {round((1 - p_bit_zero) * 100, prec)}%")
print()

for i in range(8):
	p_i = prob_any_x_bytes_non_zero(p_byte_non_zero, i, num_bytes)
	print(f"Pr[HW = {i}] = Pr[all but {i} bytes zero] = {p_i}--> {round(p_i * 100, prec)}%")


##
## PROBABILITY OF COLLISION IN snotp TERM
##
avg_byte_HW = 80.68 / 255 

print(f"Average byte HW = {avg_byte_HW} --> {round(avg_byte_HW * 100, prec)}%")

##
## PROBABILITY OF COLLISION OF 2^i b_i c and snotp TERM
## --> Requires at least one 1 in 2^i b_i c
##
prob_collision = 0
for i in range(1,num_bytes):
	p = avg_byte_HW * prob_any_x_bytes_non_zero(p_byte_non_zero, i, num_bytes)
	prob_collision += p
	#print(f"prob_collision {i}: {p}" )

print(f"Probability collision = {prob_collision} --> {round(prob_collision * 100, prec)}%")
