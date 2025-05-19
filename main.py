import sympy


def bisection_method(polynomial, start_point, end_point, epsilon=0.0001):
    x = sympy.symbols('x')
    f = sympy.lambdify(x, polynomial)
    iterations = 0

    if f(start_point) * f(end_point) >= 0:
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
    f = sympy.lambdify(x, polynomial)
    f_prime = sympy.lambdify(x, sympy.diff(polynomial, x))
    guess = (start_point + end_point) / 2
    iterations = 0

    while abs(f(guess)) > epsilon and iterations < max_iter:
        iterations += 1
        derivative = f_prime(guess)
        if derivative == 0:
            return None
        guess = guess - f(guess) / derivative
        if not (start_point <= guess <= end_point):
            return None  # guess out of bounds

    if iterations == max_iter:
        return None

    return guess, iterations


def secant_method(polynomial, start_point, end_point, epsilon=0.0001, max_iter=1000):
    x = sympy.symbols('x')
    f = sympy.lambdify(x, polynomial)
    x0, x1 = start_point, end_point
    iterations = 0

    while abs(f(x1)) > epsilon and iterations < max_iter:
        iterations += 1
        denominator = f(x1) - f(x0)
        if denominator == 0:
            return None
        x_temp = x1 - f(x1) * (x1 - x0) / denominator
        x0, x1 = x1, x_temp
        if not (start_point <= x1 <= end_point):
            return None  # out of bounds

    if iterations == max_iter:
        return None

    return x1, iterations


def main():
    x = sympy.symbols('x')
    polynomial = polynomial = x**3 - x - 2


    f = sympy.lambdify(x, polynomial)
    f_prime_expr = sympy.diff(polynomial, x)
    f_prime = sympy.lambdify(x, f_prime_expr)

    start_interval = 0
    end_interval = 4
    sub_interval_size = 0.1
    epsilon = 0.0001

    print("Choose a method to find roots:")
    print("1. Bisection Method")
    print("2. Newton-Raphson Method")
    print("3. Secant Method")
    try:
        choice = int(input("Enter your choice (1/2/3): "))
        if choice not in [1, 2, 3]:
            print("Invalid choice.")
            return
    except ValueError:
        print("Invalid input.")
        return

    found_roots = []
    current_start = start_interval

    while current_start < end_interval:
        current_end = current_start + sub_interval_size
        print(f"Checking interval: [{current_start:.2f}, {current_end:.2f}]")

        f_start = f(current_start)
        f_end = f(current_end)
        fp_start = f_prime(current_start)
        fp_end = f_prime(current_end)

        def root_already_found(r):
            return any(abs(r - existing) < epsilon for existing in found_roots)

        for val, label in [(f_start, current_start), (f_end, current_end)]:
            if abs(val) < epsilon and not root_already_found(label):
                found_roots.append(label)
                print(f"Root found at x = {label} (exact root).")

        should_search = False
        if f_start * f_end < 0 or fp_start * fp_end < 0:
            should_search = True

        if should_search:
            if choice == 1:
                result = bisection_method(polynomial, current_start, current_end, epsilon)
            elif choice == 2:
                result = newton_raphson(polynomial, current_start, current_end, epsilon)
            elif choice == 3:
                result = secant_method(polynomial, current_start, current_end, epsilon)

            if result is not None:
                root, iterations = result
                if not root_already_found(root):
                    found_roots.append(root)
                    print(f"Root found in interval [{current_start}, {current_end}]: {root} in {iterations} iterations.")
            else:
                print(f"Method did not converge in interval [{current_start}, {current_end}].")

        current_start = current_end

    print("\nSummary of found roots:")
    for i, r in enumerate(found_roots, 1):
        print(f"Root {i}: x = {r}")


if __name__ == '__main__':
    print("Starting root finding program...")
    main()
