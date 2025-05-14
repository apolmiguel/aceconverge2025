import ase.io
import pyjulip
import sys


potential = sys.argv[1]
structure = sys.argv[2]
# calc = pyjulip.ACE1("acejulia/test_pot.json")
# calc = pyjulip.ACE1("acejulia/Tr124_dim_kappa1_pot.json")
calc = pyjulip.ACE1(potential)
print(type(calc))

ats = ase.io.read(structure)
print(type(ats))
ats.calc = calc
print(type(ats.calc))

try:
    print("Try ats.get_potential_energy()")
    print(ats.get_potential_energy())
except Exception as e:
    print(f"Error computing potential energy: {e}")
    try:
        print("Try ats.calc.get_potential_energy()")
        print(ats.calc.get_potential_energy())
    except Exception as e2:
            print(f"Error in ats.calc.get_forces(): {e2}")



# if __name__ == "__main__":

