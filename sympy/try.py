import sympy as sy

x, y, z = sy.symbols("x y z")

f_xyz = x**2 + y - z

print(f_xyz.diff(x))

print(f_xyz.subs(x, 2.0555))

print((f_xyz.subs(x, 2.1).subs(y, 8.0)))

print(f_xyz.subs([(x, 2), (y, 4), (z, 0)]))