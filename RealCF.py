import numpy as np

# ==============================
#           GENERATORS
# ==============================

T1 = np.array([[1, 1], [0, 1]])   # x -> x + 1
S = np.array([[0, -1], [1,  0]])  # x -> -1/x

epsilon = 1e-10

# ==============================
#         MOBIUS ACTION
# ==============================

def mobius(M, x):
    a, b, c, d = M.ravel()
    if abs(c*x + d) < epsilon:
        return np.inf
    else:
        return (a*x + b) / (c*x + d)
    
# ==========================================
#               MAIN ALGORITHM
# ==========================================

def cont_frac_exp_real(x, max_steps = 10):
    M = np.eye(2)
    digits = []
    steps = []

    def record_step():
        steps.append({"M" : M.copy()})

    for _ in range(max_steps):
        if abs(x - round(x)) < epsilon:
            x = round(x)
        if abs(x) < epsilon:
            break

        if x > 0:
            a = int(np.floor(x + epsilon))
            digits.append(a)

            T_inv_a = np.array([[1, -a], [0,  1]])
            x = mobius(T_inv_a, x)
            M = T_inv_a @ M

            if abs(x) < epsilon:
                record_step()
                break

            x = mobius(S, x)
            M = S @ M
            record_step()

        else:
            a = int(np.floor(-x + epsilon))
            digits.append(a)

            T_a = np.array([[1, a], [0,  1]])
            x = mobius(T_a, x)
            M = T_a @ M

            if abs(x) < epsilon:
                record_step()
                break

            x = mobius(S, x)
            M = S @ M
            record_step()

    return digits, steps

# ==========================================
#           ALLOWS NUMERICAL INPUT
# ==========================================

if __name__ == "__main__":
    user_input = input("Enter a real number (allows for 'np' functions): ")
    steps_input = input("Enter maximum iterations of the algorithm (press Enter for default = 10): ")

    try:
        x = eval(user_input, {"__builtins__": None}, {"np": np})
    except Exception as e:
        print("Invalid input:", e)
        exit()

    try:
        max_steps = int(steps_input) if steps_input.strip() != "" else 10
    except:
        print("Invalid step count, using default = 10")
        max_steps = 10
    
    if max_steps <= 0:
        print("Steps must be positive, using default = 10")
        max_steps = 10

    digits, steps = cont_frac_exp_real(x, max_steps = max_steps)
    print("x =", x)

    if digits:
      output = "[" + str(digits[0])
      if len(digits) > 1:
          output += "; " + ", ".join(str(d) for d in digits[1:])
      output += "]"
    else:
      output = "[]"

    print("Continued fraction digits:", output)

    array = digits.copy()
    text=''

    def make_frac(t, i, a):
        t=r''+t[:-i]+str(a)+'+\\frac{1}{}'+t[-i:]
        return t

    i = 0
    for a in array:
        text = make_frac(text, i, a)
        i+=1

    text = r'\['+text[:-10-i]+text[-i+1:]+r'\]'
    print(text)