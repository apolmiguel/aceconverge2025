import ase.io
import pyjulip
import sys


potential = sys.argv[1]
# calc = pyjulip.ACE1("acejulia/test_pot.json")
# calc = pyjulip.ACE1("acejulia/Tr124_dim_kappa1_pot.json")
calc = pyjulip.ACE1(potential)
ats = ase.io.read("/leonardo_work/")


# ### version is adapted from `pyjulip.py`, modified to suit minimal needs ### 
# import numpy as np
# from ase.calculators.calculator import Calculator
# from ase.optimize.optimize import Optimizer 

# from julia.api import Julia
# jl = Julia(compiled_modules=False)

# from julia import Main
# Main.eval("ENV[\"JULIA_PKG_USE_CLI_GIT\"] = \"true\"")
# Main.eval("using ASE, JuLIP, ACE1")

# from julia.JuLIP import energy, forces

# ASEAtoms = Main.eval("ASEAtoms(a) = ASE.ASEAtoms(a)")
# ASECalculator = Main.eval("ASECalculator(c) = ASE.ASECalculator(c)")
# convert = Main.eval("julip_at(a) = JuLIP.Atoms(a)")

# class JulipCalculator_bare(Calculator):
#     """ 
#     ASE-compatible Calculator using bare pyjulip.py. 
#     """
#     implemented_properties = ['forces', 'energy']
#     default_parameters = {}
#     name = 'JulipCalculator_bare'

#     def __init__(self, julip_calculator):
#         Calculator.__init__(self)
#         self.julip_calculator = Main.eval(julip_calculator)

#     def calculate(self, atoms, properties, system_changes):
#         Calculator.calculate(self, atoms, properties, system_changes)
#         julia_atoms = ASEAtoms(atoms)
#         julia_atoms = convert(julia_atoms)
#         self.results = {}
#         if 'energy' in properties:
#             E = energy(self.julip_calculator, julia_atoms)
#             self.results['energy'] = E
#         if 'forces' in properties:
#             self.results['forces'] = np.array(forces(self.julip_calculator, julia_atoms))


# if __name__ == "__main__":
#     # load the potential file 
#     potname = '/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/acejulia/Tr124_dim_kappa1_pot.json'
#     print('Loading potential file:', potname,'\nAlso loading IP from the potential file.')
#     # Main.eval("D = load_dict(\"" + potname + "\"); IP = read_dict(D, \"IP\")")
#     !julia -e "D = load_dict(\"/leonardo_work/Sis25_degironc_0/apol/aceconverge2025/acejulia/Tr124_dim_kappa1_pot.json\"); IP = read_dict(D, \"IP\")"
#     # Main.eval("IP = read_dict(D, \"IP\")")
#     print('Loading IP into ASE.')
#     ASE_IP = JulipCalculator_bare("IP")
