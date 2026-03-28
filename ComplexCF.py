import numpy as np

# Allow mathematical inputs to find continued fractions of
def safe_eval(input):
    input = input.replace("i", "j")    # allow i instead of j for complex
    return complex(eval(input, {"__builtins__": None}, {"np": np}))

# ==============================
#           GENERATORS
# ==============================

T1 = np.array([[1, 1], [0, 1]])
Ti = np.array([[1, 1j], [0, 1]])
phi = np.array([[0, 1], [1,  0]])
kappa = np.array([[1j, 0], [0, 1]])

T1_inv = np.array([[1, -1], [0, 1]])
Ti_inv = np.array([[1, -1j], [0, 1]])

epsilon = 1e-10

# ==================================
#       EXTENDED MOBIUS ACTION
# ==================================

def mobius(M, z):
    a, b, c, d = M.ravel()
    if abs(c*z + d) < epsilon:
        return complex('inf')
    return (a*z + b) / (c*z + d)

# =======================================
#       DEFINING TRANSLATION DIGITS
# =======================================

translation_map = {
    "T1": 1,
    "T1_inv": -1,
    "Ti": 1j,
    "Ti_inv": -1j,
}

# =============================================
#    REDUCING THE WORD TO TRANSLATIONS & PHI
# =============================================

# phi^2=id so can cancel consecutive appearances
def cancel_phi_pairs(word):
    new_word = []
    i = 0
    while i < len(word):
        if i+1 < len(word) and word[i] == "phi" and word[i+1] == "phi":
            i += 2
        else:
            new_word.append(word[i])
            i += 1
    return new_word

# After cancelling, any consecutive translations can be combined into a single one
def combine_translations(word):
    new_word = []
    i = 0
    rotation = 0

    def rotate(val, r):
        return val * (1j ** r)

    while i < len(word):

        g = word[i]

        # Accounting for rotations before a translation by conjugating
        if g == "kappa":
            rotation = (rotation + 1) % 4
            new_word.append(g)
            i += 1
            continue

        if (i+1 < len(word) and
            word[i] in translation_map and
            word[i+1] in translation_map):

            v1 = rotate(translation_map[word[i]], rotation)
            v2 = rotate(translation_map[word[i+1]], rotation)

            val = v1 + v2
            val = val * (1j ** (-rotation))

            a = int(round(val.real))
            b = int(round(val.imag))

            if a > 0:
                new_word += ["T1"] * a
            elif a < 0:
                new_word += ["T1_inv"] * (-a)

            if b > 0:
                new_word += ["Ti"] * b
            elif b < 0:
                new_word += ["Ti_inv"] * (-b)

            i += 2
            continue

        new_word.append(g)
        i += 1

    return new_word

# Use Lemma 6.6 to shift the kappa rotations to the end of the word
def push_kappa_to_end(word):
    new_word = []

    current = 0
    block_index = 0
    kappa_power = 0

    def flush_block():
        nonlocal current, block_index

        if current != 0:
            factor = (1j ** kappa_power) * ((-1) ** (block_index * kappa_power))
            val = current * factor

            a = int(round(val.real))
            b = int(round(val.imag))

            if a > 0:
                new_word.extend(["T1"] * a)
            elif a < 0:
                new_word.extend(["T1_inv"] * (-a))

            if b > 0:
                new_word.extend(["Ti"] * b)
            elif b < 0:
                new_word.extend(["Ti_inv"] * (-b))

        current = 0

    for g in word:
        if g in translation_map:
            current += translation_map[g]

        elif g == "phi":
            flush_block()
            new_word.append("phi")
            block_index += 1

        elif g == "kappa":
            flush_block()
            kappa_power = (kappa_power + 1) % 4

    new_word.extend(["kappa"] * kappa_power)
    return new_word

# Truncate the word to remove unnecessary rotations and translations. Theses are invariant at infinity
def truncate_at_last_phi(word):
    for i in range(len(word)-1, -1, -1):
        if word[i] == "phi":
            return word[:i+1]
    return word

# =================================================================
#       OUTPUTTING THE STAGES OF THE `WORD TO DIGITS` PROCESS
# =================================================================

# Neaten the presentation fo the output 
def present(word):
    result = []
    current = 0 + 0j

    def flush():
        nonlocal current
        if abs(current) > 0:
            a = int(round(current.real))
            b = int(round(current.imag))

            val = complex(a, b)
            result.append(f"T({format_complex(val)})")

            current = 0 + 0j

    for g in word:
        if g in translation_map:
            current += translation_map[g]
        else:
            flush()
            result.append(g)

    flush()
    return " ".join(result)

# ============================================
#       WORD → DIGITS
# ============================================

def word_to_digits(word):
    digits = []
    current = 0

    if len(word) > 0 and word[0] == "phi":
        digits.append(0)

    for g in word:
        if g == "phi":
            if current != 0:
                digits.append(current)
                current = 0
        elif g in translation_map:
            current += translation_map[g]

    if current != 0:
        digits.append(current)

    return digits

def format_complex(z):
    a = int(round(z.real))
    b = int(round(z.imag))

    parts = []

    # Real part
    if a != 0:
        parts.append(str(a))

    # Imaginary part
    if b != 0:
        if b == 1:
            imag = "i"
        elif b == -1:
            imag = "-i"
        else:
            imag = f"{b}i"

        if a != 0 and b > 0:
            parts.append("+" + imag)
        else:
            parts.append(imag)

    if not parts:
        return "0"

    return "".join(parts)

# ==========================
#       SUBREGION MAPS
# ==========================

f1 = phi @ T1 @ Ti_inv @ phi @ kappa**3
f2 = phi @ Ti_inv @ phi
f3 = phi @ T1 @ Ti_inv @ phi @ T1_inv @ Ti_inv @ kappa
f4 = T1 @ Ti @ phi @ T1_inv @ Ti

f1_inv = kappa @ phi @ Ti @ T1_inv @ phi
f2_inv = phi @ Ti @ phi
f3_inv = kappa**3 @ Ti @ T1 @ phi @ Ti @ T1_inv @ phi
f4_inv = Ti_inv @ T1 @ phi @ Ti_inv @ T1_inv

f1_word = ["phi", "T1", "Ti_inv", "phi", "kappa", "kappa", "kappa"]
f2_word = ["phi", "Ti_inv", "phi"]
f3_word = ["phi", "T1", "Ti_inv", "phi", "T1_inv", "Ti_inv", "kappa"]
f4_word = ["T1", "Ti", "phi", "T1_inv", "Ti"]

# ===============================
#       SUBREGION SELECTION
# ===============================

def choose_map(z):
    x, y = z.real, z.imag
    eps = 1e-12

    if abs(z) < epsilon:
        return None, [], None
    if (y <= x + eps) and (y < 1 - x - eps):
        return f1_inv, f1_word, "f1"
    elif (y > x + eps) and (y <= 1 - x + eps):
        return f2_inv, f2_word, "f2"
    elif (y >= x - eps) and (y >= 1 - x - eps):
        return f3_inv, f3_word, "f3"
    else:
        return f4_inv, f4_word, "f4"

# ==========================================
#               MAIN ALGORITHM
# ==========================================

def cont_frac_exp_complex(z, max_steps=10):

    steps = []
    final_word = []
    symbolic_word = []

    for _ in range(max_steps):

        if not np.isfinite(z.real) or not np.isfinite(z.imag):
            break
        if abs(z) < epsilon:
            break

        a = int(np.floor(z.real + 1e-12))
        b = int(np.floor(z.imag + 1e-12))

        if a != 0 or b != 0:
            symbolic_word.append(f"T({format_complex(complex(a, b))})")

        if a > 0:
            final_word += ["T1"] * a
        elif a < 0:
            final_word += ["T1_inv"] * (-a)

        if b > 0:
            final_word += ["Ti"] * b
        elif b < 0:
            final_word += ["Ti_inv"] * (-b)

        T_inv = np.array([[1, -(a + b*1j)], [0, 1]])
        z = mobius(T_inv, z)

        if abs(z) < epsilon:
            break

        f_inv, f_word, f_label = choose_map(z)

        if f_inv is None:
            break

        symbolic_word.append(f_label)

        z = mobius(f_inv, z)
        final_word += f_word

        if abs(z) < epsilon:
            break

        steps.append(z)

    if abs(z) < epsilon:
        final_word.append("phi")
        symbolic_word.append("phi")

    reduced_word = final_word
    while True:
        new_word = cancel_phi_pairs(reduced_word)
        new_word = combine_translations(new_word)

        if new_word == reduced_word:
            break

        reduced_word = new_word

    kappa_reduced = push_kappa_to_end(reduced_word)
    truncated = truncate_at_last_phi(kappa_reduced)

    print("\n--- SYMBOLIC WORD ---")
    print(" ".join(symbolic_word))

    print("\n--- DECOMPOSED WORD ---")
    print(present(final_word))

    print("\n--- REDUCED WORD ---")
    print(present(kappa_reduced))

    print("\n--- TRUNCATED WORD ---")
    print(present(truncated))

    digits = word_to_digits(truncated)
    return digits, steps

# ==========================================
#           ALLOWS NUMERICAL INPUT
# ==========================================
if __name__ == "__main__":
    user_input = input("Enter a complex number (allows for 'np' functions): ")
    
    try:
        z = safe_eval(user_input)
    except Exception as e:
        print("Invalid input:", e)
        exit()

    digits, steps = cont_frac_exp_complex(z)
    print("\nz =", str(z).replace('j', 'i').replace('(', '').replace(')', ''))

    formatted_digits = [format_complex(d) for d in digits]
    if formatted_digits:
        output = "[" + formatted_digits[0]
        if len(formatted_digits) > 1:
            output += "; " + ", ".join(formatted_digits[1:])
        output += "]"
    else:
        output = "[]"

    print("Continued Fraction Digits:", output)

    # Write the continued fraction
    array = formatted_digits.copy()
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