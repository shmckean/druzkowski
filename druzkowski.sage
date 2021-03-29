# This is a modification of a1-degree.sage (written jointly with Thomas Brazelton and Sabrina Pauli).
# Input a matrix A as a list of row vectors.
# druzkowski.sage will compute the Jacobian and Bezoutian of the Druzkowski morphism with matrix A.

k = QQbar;

A = [[0,1,0],[-1,0,1],[0,1,0]]; # INPUT: matrix of Druzkowski morphism
q = [1,2,3]; # INPUT: point at which to check injectivity

n = len(A);
R = PolynomialRing(k,n,'x');
x = R.gens();

f = []; # (note: computes Druzkowski morphism f_A-q)
for i in range(n):
    f.append(-q[i]+x[i]+(sum([A[i][j]*x[j] for j in range(n)]))^3)

print('Jacobian :',jacobian(f,x).det())
    
Rxy = PolynomialRing(k,n,var_array=['X','Y']); # (note: computes Bezoutian)
L = Rxy.gens();
X = [L[2*i] for i in range(n)];
Y = [L[2*i+1] for i in range(n)];
W = X.copy();
Z = [X];
for i in range(n):
    W = W.copy();
    W[i] = Y[i];
    Z.append(W);
D = [];
for i in range(n):
    for j in range(n):
        D.append((f[i](Z[j])-f[i](Z[j+1]))/(X[j]-Y[j]));
Bez = matrix(Rxy,n,D).det();

F = []; # (note: compute ideal in k[X,Y])
for i in range(n):
    F.append(f[i](X));
    F.append(f[i](Y));
J = F*Rxy;

Q = Rxy.quo(J); # (note: computes reduced Bezoutian)
Bez = Q(Bez).lift();
print('Bezoutian:',Bez)
