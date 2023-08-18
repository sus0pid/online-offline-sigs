
from falcon import *
from fft import fft, ifft
from random import *
from numpy import linalg

sk = SecretKey(32)
pk = PublicKey(sk)

a = sk.B0_fft[0][0]
b = sk.B0_fft[0][1]
c = sk.B0_fft[1][0]
d = sk.B0_fft[1][1]

def v_add(a, b):
    return list(map(lambda x,y:x+y,a,b))

def v_sub(a, b):
    return list(map(lambda x,y:x-y, a, b))

def v_mul(a, b):
    return list(map(lambda x,y:x*y, a, b))

#def v_inv(a):
#    return list(map(lambda x: complex(*operator.inv(x)), vector(CC, a)))

def v_round(a):
    a = ifft(a)
    a = [round(a_) for a_ in a]
    a = fft(a)
    return a

def mat_mul(A, x):
    y = [[0]*32 for _ in range(2)]
    for i in range(2):
        for j in range(2):
            y[i] = v_add(y[i], v_mul(A[i][j], x[j]))
    return y

a_inv = a
bc = v_mul(b, c)
bc_a = v_mul(v_mul(b, c), a_inv)
aa = v_mul(a, a)
bc_a_d = v_sub(bc_a, d)

a_ = v_sub(a_inv, v_mul(bc, (v_mul(aa, bc_a_d))))
b_ = v_mul(b, (v_mul(a, bc_a_d)))
c_ = v_mul(c, (v_mul(a, bc_a_d)))
d_ = v_sub([0]*32, (bc_a_d))

B0_inv_fft = [[a_, b_], [c_, d_]]

# construct some input x
x = [fft([randint(-1000,1000) for _ in range(32)]), fft([randint(-1000, 1000) for _ in range(32)])]
# inverse
y = mat_mul(B0_inv_fft, x)
# round
y = [v_round(y[0]), v_round(y[1])]
# forward
x_ = mat_mul(sk.B0_fft, y)
print(x_)
print(len(x_))
print(len(x_[0]))
# check distance
dist = linalg.norm(v_sub(ifft(x_[0]), ifft(x[0]))), linalg.norm(ifft(x[0]))
print(dist)