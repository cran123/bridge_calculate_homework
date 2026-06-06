from pathlib import Path

p = Path.home() / "bridge3_stappp_sparse/data/data-3/Bridge-3.rcm.generated.out"

keys = (
    "TOTAL SPARSE",
    "NUMBER OF EQUATIONS",
    "NUMBER OF SPARSE",
    "ESTIMATED SPARSE",
    "TIME FOR INPUT PHASE",
    "TIME FOR CALCULATION OF STIFFNESS MATRIX",
    "TIME FOR FACTORIZATION",
    "TIME FOR LOAD CASE SOLUTIONS",
    "TIME FOR RESULT OUTPUT",
    "T O T A L   S O L U T I O N",
)

print("summary_lines")
for line in p.open(errors="ignore"):
    if any(k in line for k in keys):
        print(line.rstrip())

max_umag = 0.0
max_node = None
max_vals = None
in_disp = False
for line in p.open(errors="ignore"):
    if "D I S P L A C E M E N T S" in line:
        in_disp = True
        continue
    if in_disp:
        s = line.split()
        if len(s) >= 7 and s[0].isdigit():
            node = int(s[0])
            vals = [float(x) for x in s[1:7]]
            umag = (vals[0] ** 2 + vals[1] ** 2 + vals[2] ** 2) ** 0.5
            if umag > max_umag:
                max_umag = umag
                max_node = node
                max_vals = vals
        elif "S T R E S S" in line:
            break

print("max_displacement")
print(max_node)
print(max_umag)
print(max_vals)
