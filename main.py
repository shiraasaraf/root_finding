import sympy

def bisection_method(polynomial, start_point, end_point, epsilon=0.0001):
    x = sympy.symbols('x')
    f = sympy.utilities.lambdify(x, polynomial)
    iterations = 0

    if f(start_point) * f(end_point) >= 0:
        # Cannot guarantee a root in this interval
        return None

    while (end_point - start_point) / 2 > epsilon:
        iterations += 1
        midpoint = (start_point + end_point) / 2
        f_mid = f(midpoint)
        if abs(f_mid) < epsilon:
            return midpoint, iterations
        elif f(start_point) * f_mid < 0:
            end_point = midpoint
        else:
            start_point = midpoint

    root = (start_point + end_point) / 2
    return root, iterations


def newton_raphson(polynomial, start_point, end_point, epsilon=0.0001, max_iter=1000):
    x = sympy.symbols('x')
    f = sympy.utilities.lambdify(x, polynomial)
    f_prime = sympy.utilities.lambdify(x, sympy.diff(polynomial, x))
    guess = (start_point + end_point) / 2
    iterations = 0

    while abs(f(guess)) > epsilon and iterations < max_iter:
        iterations += 1
        derivative = f_prime(guess)
        if derivative == 0:
            # Failed due to zero derivative
            return None
        guess = guess - f(guess) / derivative

    if iterations == max_iter:
        return None

    return guess, iterations


def secant_method(polynomial, start_point, end_point, epsilon=0.0001, max_iter=1000):
    x = sympy.symbols('x')
    f = sympy.utilities.lambdify(x, polynomial)
    x0, x1 = start_point, end_point
    iterations = 0

    while abs(f(x1)) > epsilon and iterations < max_iter:
        iterations += 1
        denominator = f(x1) - f(x0)
        if denominator == 0:
            return None
        x_temp = x1 - f(x1) * (x1 - x0) / denominator
        x0, x1 = x1, x_temp

    if iterations == max_iter:
        return None

    return x1, iterations


def main():
    x = sympy.symbols('x')
    polynomial = x**3 - 6*x**2 + 11*x - 6  # example polynomial
    f = sympy.utilities.lambdify(x, polynomial)

    start_interval = 0
    end_interval = 4
    sub_interval_size = 0.5
    epsilon = 0.0001

    print("Choose a method to find roots:")
    print("1. Bisection Method")
    print("2. Newton-Raphson Method")
    print("3. Secant Method")
    choice = int(input("Enter your choice (1/2/3): "))

    found_roots = []  # list of found roots to avoid duplicates
    current_start = start_interval
    while current_start < end_interval:
        current_end = current_start + sub_interval_size
        print(f"Checking interval: [{current_start}, {current_end}]")

        f_start = f(current_start)
        f_end = f(current_end)

        # If function value at interval edge is close to zero, there is an exact root - add if not duplicate
        if abs(f_start) < epsilon:
            root = current_start
            # Check if root was already found
            if not any(abs(root - r) < epsilon for r in found_roots):
                found_roots.append(root)
                print(f"Root found at x = {root} (exact root).")

        if abs(f_end) < epsilon:
            root = current_end
            if not any(abs(root - r) < epsilon for r in found_roots):
                found_roots.append(root)
                print(f"Root found at x = {root} (exact root).")

        # If there is a sign change, search for root in the interval
        if f_start * f_end < 0:
            if choice == 1:
                result = bisection_method(polynomial, current_start, current_end, epsilon)
            elif choice == 2:
                result = newton_raphson(polynomial, current_start, current_end, epsilon)
            elif choice == 3:
                result = secant_method(polynomial, current_start, current_end, epsilon)
            else:
                print("Invalid choice.")
                return

            if result is not None:
                root, iterations = result
                # If root is close to any found root, don't print duplicate
                if not any(abs(root - r) < epsilon for r in found_roots):
                    found_roots.append(root)
                    print(f"Root found in interval [{current_start}, {current_end}]: {root} in {iterations} iterations.")

        current_start = current_end


if __name__ == '__main__':
    print("Starting root finding program...")
    main()
