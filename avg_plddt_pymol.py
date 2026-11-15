from pymol import cmd

def avg_plddt(selection="all"):
    global stored
    stored = []  # Define 'stored' in global scope so PyMOL can access it
    cmd.iterate(selection, "stored.append(b)")
    total = sum(stored)
    count = len(stored)
    avg = total / count if count else 0
    print(f"Average pLDDT for {selection}: {avg:.2f}")

cmd.extend("avg_plddt", avg_plddt)
