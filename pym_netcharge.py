from pymol import cmd
from collections import Counter

charges = {'K': 1, 'R': 1, 'D': -1, 'E': -1, 'H': 1}

def net_charge(selection='all'):
    seq = cmd.get_fastastr(selection).split('\n', 1)[1].replace('\n', '')
    residue_counts = Counter(seq)
    net_chg = sum(charges.get(res, 0) * residue_counts.get(res, 0) for res in charges)

    print(f"Residue counts in '{selection}':")
    for res in ['K', 'R', 'H', 'D', 'E']:
        print(f"  {res}: {residue_counts.get(res, 0)}")
    
    print(f"\nNet charge for '{selection}': {net_chg:.2f}")

cmd.extend('net_charge', net_charge)
