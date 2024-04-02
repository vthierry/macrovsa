#
# Checks VTB algebra operations
#

with(LinearAlgebra):
with(Statistics):

ok := r -> evalb(0 = Norm(normal(r))):

# Defines the dimension as a square

d_:= 3: d := d_^2:

# Defines literal vectors

v := s -> Vector(d, (i) -> cat(s, '_', i)):

x := v('x'): y := v('y'): z := v('z'):

# Defines the binding matrix

B := v -> Matrix(d, d, (j, i) -> if 1 <= i - d_ * iquo(j-1, d_) and i - d_ * iquo(j-1, d_) <= d_ then sqrt(d_) * v[1 + irem(i-1, d_) + d_ * irem(j-1, d_)] else 0 fi):

# Dual forms of the binding matrix

B_ := v -> Matrix(d, d, (j, i) -> if 1 <= j - d_ * iquo(i-1, d_) and j - d_ * iquo(i-1, d_) <= d_ then sqrt(d_) * v[1 + irem(i-1, d_) + d_ * irem(j-1, d_)] else 0 fi):

ok_binding_dual := ok(B(y) - B_(y));

# Numerical comparison of complexity

for d0 in [10, 100, 500, 1000]
do [d0, evalf(log10(d0 * log(d0))), evalf(log10(d0 * sqrt(d0))), evalf((d0 * sqrt(d0)/(d0 * log(d0))))]
od;

# Defines the binding explicit formula

alpha := i -> d_ * iquo(i-1, d_):
beta := i -> d_ * irem(i-1, d_):

b := proc(y, x) local k: Vector(d, (i) -> sqrt(d_) * add(y[k + beta(i)] * x[k + alpha(i)], k = 1 .. d_)) end:

# Checks that B(y) x = b(y, x)

ok_binding := ok(Multiply(B(y), x) - b(y, x));

# Defines the explicit transpose formula

sigma := (i) -> 1 + d_ * irem(i-1, d_) + iquo(i-1, d_):

t := (x) -> Vector(d, (i) -> x[sigma(i)]):

ok_idempotent := ok(Vector(d, (i) -> sigma(sigma(i)) - i));
ok_transpose := ok(B(t(y))- Transpose(B(y)));

# Defines the identity vector

u := Vector(d, (i) -> if i = sigma(i) then 1/sqrt(d_) else 0 fi):

ok_identity := ok(B(u) - IdentityMatrix(d, d));

# Defines the left mirror matrix

M := Matrix(d, d, (i, j) -> if j = sigma(i) then 1 else 0 fi):

ok_mirror            := ok(Multiply(M, Multiply(B(x), y)) - Multiply(B(y), x));
ok_mirror_idempotent := ok(Multiply(M, Multiply(M, B(y))) - B(y));
ok_mirror_explicit   := ok(Multiply(M, x) - t(x));

# Verifies the c(y, x) = B_y~ B_<-> x~ formula

c := proc(y, x) local k: Vector(d, (i) -> sqrt(d_) * add(y[sigma(k + beta(i))] * x[sigma(k + alpha(i))], k = 1 .. d_)) end:

ok_inverse := ok(Multiply(B(t(y)), Multiply(M, x)) - c(y, x));

# Verifying that mirror matrix cannot be a binding matrix

ok_non_solvable := evalb(1 in convert(M - B(x),set));

ok_non_solvable := evalb(1 in convert(Multiply(B(u), M) - B(x), set));

ok_non_solvable := evalb(1 in convert(Multiply(M, B(u)) - B(x), set));

# Defines the explicit composition formula

o := proc(y, x) local k: Vector(d, (i) -> sqrt(d_) * add(y[k + d_ * iquo(i-1, d_)] * x[1 + d_ * (k - 1) + irem(i-1, d_)], k = 1 .. d_)) end:

ok_composition := ok(B(o(y, x)) - Multiply(B(y), B(x)));
ok_composition := ok(t(o(y, x)) - o(t(x), t(y)));
ok_composition := ok(o(z, o(y, x)) - o(o(z, y), x));
ok_composition := ok(o(l * x + y, z) - (l * o(x, z) + o(y, z)));

# Checks the normalization factor

d_ := 16: d := d_^2:

draw := proc() convert(Sample(RandomVariable(Normal(0, 1 / evalf(sqrt(d)))), d), Vector) end:

DataSummary(map(proc() Norm(b(draw(), draw()), Euclidean) end, [$100]));

