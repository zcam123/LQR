using LinearAlgebra
using MatrixEquations

#define system
A = [1 -6.66e-13 -2.03e-9 -4.14e-6;
     9.83e-4 1 -4.09e-8 -8.32e-5;
     4.83e-7 9.83e-4 1 -5.34e-4;
     1.58e-10 4.83e-7 9.83e-4 .9994;
]

B = [9.83e-4, 4.83e-7, 1.58e-10, 3.89e-14]

C = [-.0096 .0135 .005 -.0095]

D = [0.0]

#for getting firing rate Z
function y_to_z(y)
    return ℯ^(61.4*y - 5.468);
end


z2y(z) = (log(z) + 5.468)/61.4
y2z(y) = exp(61.4*y - 5.468)

function sim(steps, x0, u=0)
    z = zeros(0);
    x = x0;
    for t in 1:steps
        x = A*x + B .* u
        y = C*x
        z1 = y_to_z(y[1]);
        append!(z, z1)
    end
    return z
end

#run simulation and plot results
zs = sim(50000, [-1.548; 2.18; .806; -1.53]);
x = [i for i in 1:50000]
plot(x, zs, xlabel = "time steps", ylabel = "Z")


#make LQR controller using MatrixEquations.jl
Q = C'*C
R = (1e-3)*I
P = ared(A, B, R, Q)[1] #index of 1 to retrieve only the real solution
K = inv(R + B'*P*B) * (B'*P*A)

print(P)
print("\n",K)

function control(T, zD, x0, A, B, K)
    # yD = z2y(zD)
    # uD = inv(C*inv(I - A)*B) * yD
    # xD = inv(I - A) * B * uD
    xD = [75.03983654, 481.45636789, 540.76706284, 886.34429676]

    zs = []
    us = []
    append!(zs, C*x0)
    x = x0
    for t in 1:T
        u = -K*(x - xD) .+ uD
        append!(us, u)
        x = A*x + B.*u
        y = C*x
        append!(zs, y2z(y[1]))
    end
    return zs, us
end

zs, us = control(10^5, 0.2, zeros(4), A, B, K)

zDs = [0.2 for _ in 1:10^5]
plot([zs, zDs])