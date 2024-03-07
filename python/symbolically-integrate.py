# Symbolically derive definite integral of AA over a range of log(g_sw) values
import sympy as sp

# Define symbols
x, b0a, b1a, b2a, b0h, b1h, b2h = sp.symbols('x b0a b1a b2a b0h b1h b2h')

# Define the expression
expr = sp.ln((b0a + b1a * x + b2a * x ** 2) / (b0h + b1h * x + b2h * x ** 2))

# Integrate with respect to x
integral_expr = sp.integrate(expr, x)

# Print the result
print(integral_expr.simplify)
